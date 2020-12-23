# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""

.. include:: ../include/parset_usage.rst

.. todo::
    - Add range and length parameters allowing one to define the range
      allowed for the parameter values and number of elements required
      (if the parameter is an array)
    - Allow for a from_par_file classmethod to initialize the parameter
      set based on a yanny parameter file.
    - Save the defaults and allow for a revert_to_default function.
    - Write an __add__ function that will all you to add multiple
      parameter sets.

----

.. include license and copyright
.. include:: ../include/copy.rst

----

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

import os
import inspect
import warnings
import textwrap

from IPython import embed

import numpy

from configobj import ConfigObj

from .util import recursive_dict_evaluate

class ParSet:
    """
    Generic base class to handle and manipulate a list of operational
    parameters.  A glorified dictionary that constrains and types its
    components.

    Args:
        pars (:obj:`list`):
            A list of keywords for a list of parameter values.
        values (:obj:`list`, optional):
            Initialize the parameters to these values.  If not provided,
            all parameters are initialized to `None` or the provided
            default.
        defaults (:obj:`list`, optional):
            For any parameters not provided in the *values* list, use
            these default values.  If not provided, no defaults are
            assumed.
        options (:obj:`list`, optional):
            Force the parameters to be one of a list of options.  Each
            element in the list can be a list itself.  If not provided,
            all parameters are allowed to take on any value within the
            allowed data type.
        dtypes (:obj:`list`, optional):
            Force the parameter to be one of a list of data types.  Each
            element in the list can be a list itself.  If not provided,
            all parameters are allowed to have any data type.
        can_call (:obj:`list`, optional): Flag that the parameters are
            callable operations.  Default is False.
        descr (:obj:`list`, optional):
            A list of parameter descriptions.  Empty strings by default.
        cfg_section (:obj:`str`, optional): 
            The top-level designation for a configuration section
            written based on the contents of this parameter set.
        cfg_comment (:obj:`str`, optional): 
            Comment to be placed at the top-level of the configuration
            section written based on the contents of this parameter set.

    Raises:
        TypeError:
            Raised if the input parameters are not lists or if the input
            keys are not strings.
        ValueError:
            Raised if any of the optional arguments do not have the same
            length as the input list of parameter keys.

    Attributes:
        npar (:obj:`int`):
            Number of parameters
        data (:obj:`dict`):
            Dictionary with the parameter values
        default (:obj:`dict`):
            Dictionary with the default values
        options (:obj:`dict`):
            Dictionary with the allowed options for the parameter values
        dtype (:obj:`dict`):
            Dictionary with the allowed data types for the parameters
        can_call (:obj:`dict`):
            Dictionary with the callable flags
        descr (:obj:`dict`):
            Dictionary with the description of each parameter.
        cfg_section (:obj:`str`): 
            The top-level designation for a configuration section
            written based on the contents of this parameter set.
        cfg_comment (:obj:`str`): 
            Comment to be placed at the top-level of the configuration
            section written based on the contents of this parameter set.
    """

    prefix = 'PAR'
    """
    Class prefix for header keywords when writing the parset to an
    `astropy.io.fits.Header`_ object.
    """

    def __init__(self, pars, values=None, defaults=None, options=None, dtypes=None, can_call=None,
                 descr=None, cfg_section=None, cfg_comment=None):
        # Check that the list of input parameters is a list of strings
        if not isinstance(pars, list):
            raise TypeError('Input parameter keys must be provided as a list.')
        for key in pars:
            if not isinstance(key, str):
                raise TypeError('Input parameter keys must be strings.')
        
        # Get the length of the parameter list and make sure the list
        # has unique values
        self.npar = len(pars)
        if len(numpy.unique(numpy.array(pars))) != self.npar:
            raise ValueError('All input parameter keys must be unique.')

        # Check that the other lists, if provided, have the correct type
        # and length
        if values is not None and (not isinstance(values, list) or len(values) != self.npar):
            raise ValueError('Values must be a list with the same length as the keys list.')
        if defaults is not None and (not isinstance(defaults, list) or len(defaults) != self.npar):
            raise ValueError('Defaults must be a list with the same length as the keys list.')
        if options is not None and (not isinstance(options, list) or len(options) != self.npar):
            raise ValueError('Options must be a list with the same length as the keys list.')
        if dtypes is not None and (not isinstance(dtypes, list) or len(dtypes) != self.npar):
            raise ValueError('Data types list must have the same length as the keys list.')
        if can_call is not None and (not isinstance(can_call, list) or len(can_call) != self.npar):
            raise ValueError('List of callable flags must have the same length as keys list.')
        if descr is not None and (not isinstance(descr, list) or len(descr) != self.npar):
            raise ValueError('List of parameter descriptions must have the same length as '
                             'keys list.')

        # Set up dummy lists for no input
        _values = [None]*self.npar if values is None else values
        _defaults = [None]*self.npar if defaults is None else defaults
        _options = [None]*self.npar if options is None else options
        _dtypes = [None]*self.npar if dtypes is None else dtypes
        _can_call = [False]*self.npar if can_call is None else can_call
        _descr = ['']*self.npar if descr is None else descr

        # Set the defaults
        self.default = dict([ (p, d) for p, d in zip(pars, _defaults) ])

        # Set the valid options
        self.options = dict([ (p, [o]) if o is not None and not isinstance(o, list) else (p, o) \
                                       for p, o in zip(pars, _options) ])
        # Set the valid types
        self.dtype = dict([ (p, [t]) if not isinstance(t, list) else (p, t) \
                                     for p, t in zip(pars, _dtypes) ])

        # Set the calling flags
        self.can_call = dict([ (p, t) for p, t in zip(pars, _can_call) ])

        # Set the calling flags
        self.descr = dict([ (p, t) for p, t in zip(pars, _descr) ])

        # Set the data dictionary using the overloaded
        # __setitem__function so that value checking is performed
        self.data = dict.fromkeys(pars)
        for p, d, v, t in zip(pars, _defaults, _values, _dtypes):
            # Check if 'None' is an allowed option
            none_allowed = False
            if type(t) is list:
                if type(None) in t:
                    none_allowed = True
            if v is None and not none_allowed:
                self.__setitem__(p, d)
                continue
            self.__setitem__(p, v)

        # Save the configuration file section details
        self.cfg_section = cfg_section
        self.cfg_comment = cfg_comment

    def __getitem__(self, key):
        """
        Return the value of the designated key.

        Args:
            key (:obj:`str`):
                Key for new parameter
        """
        return self.data[key]

    def __setitem__(self, key, value):
        """
        Set the value for a key.

        Args:
            key (:obj:`str`):
                Key for new parameter
            value (object):
                Parameter value, must have the data type specified by
                :attr:`dtype` when instantiating the object, if any were
                provided.

        Raises:
            KeyError:
                Raised if the keyword is not valid for this class.
            ValueError:
                Raised if the parameter value is not among the allowed
                options (:attr:`options`).
            TypeError:
                Raised if the parameter value does not have an allowed
                data type (:attr:`dtype`) or if the provided value is
                not a callable object, but is expected to be by
                :attr:`can_call`.
        """
        if key not in self.keys():
            raise KeyError('{0} is not a valid key for {1}.'.format(key, self.__class__.__name__))

        if value is None:
            self.data[key] = value
            return

        if isinstance(value, list):
            is_parset_or_dict = [ isinstance(v, (ParSet, dict)) for v in value ]
            if numpy.any(is_parset_or_dict) and not numpy.all(is_parset_or_dict):
                warnings.warn('List includes a mix of ParSet and dicts with other types.  '
                              'Displaying and writing the ParSet will not be correct!')

        if self.options[key] is not None and value not in self.options[key]:
            raise ValueError('Input value for {0} invalid: {1}.\nOptions are: {2}'.format(
                                                                    key, value, self.options[key]))

        if self.dtype[key] is not None \
                and not any([ d is None or isinstance(value, d) for d in self.dtype[key]]):
            raise TypeError('Input value for {0} has incorrect type: {1}.'.format(key, value) +
                            '\nValid types are: {0}'.format(self.dtype[key]))

        if self.can_call[key] and not callable(value):
            raise TypeError('{0} is not a callable object.'.format(value))

        self.data[key] = value

    def __len__(self):
        """Return the number of parameters."""
        return self.npar
        
    def __iter__(self):
        """Return an iterable to the parameter values."""
        return iter(self.data.values())

    def __repr__(self):
        """Return a string representation of the parameters."""
        return self._output_string(header=self.cfg_section)

    def _output_string(self, header=None, value_only=False):
        """
        Constructs the short-format table strings for the
        ``__repr__`` method.

        Args:
            header (:obj:`str`, optional):
                String header to provide for the table.  This is
                typically the name of the configuration section.
            value_only (:obj:`bool`, optional):
                By default, the table includes the parameter key, its
                current value, the default value, its data type, and if
                the value can be a callable function.  If
                `value_only=True`, only the parameter key and current
                value are returned.

        Returns:
            str: Single long string with the parameter table for the
            ``__repr__`` method.
        """
        additional_par_strings = []
        ncol = 2 if value_only else 5
        data_table = numpy.empty((self.npar+1, ncol), dtype=object)
        data_table[0,:] = ['Parameter', 'Value'] if value_only \
                            else ['Parameter', 'Value', 'Default', 'Type', 'Callable']
        for i, k in enumerate(self.keys()):
            data_table[i+1,0] = k
            if isinstance(self.data[k], ParSet):
                _header = k if header is None else '{0}:{1}'.format(header, k)
                additional_par_strings += [ self.data[k]._output_string(header=_header,
                                                                        value_only=value_only) ]
                data_table[i+1,1] = 'see below'
                if not value_only:
                    data_table[i+1,2] = 'see below'
            else:
                data_table[i+1,1] = ParSet._data_string(self.data[k])
                if not value_only:
                    data_table[i+1,2] = ParSet._data_string(self.default[k])
            if value_only:
                continue
            data_table[i+1,3] = ', '.join(self._types_list(k))
            data_table[i+1,4] = self.can_call[k].__repr__()

        output = [ParSet._data_table_string(data_table)]
        if header is not None:
            output = [header] + output
        if len(additional_par_strings) > 0:
            output += additional_par_strings
        return '\n'.join(output)

    @staticmethod
    def _data_table_string(data_table, delimiter='print'):
        """
        Provided the array of data, format it with equally spaced
        columns and add a header (first row) and contents delimiter.

        Args:
            data_table (`numpy.ndarray`_):
                Array of string representations of the data to print.

        Returns:
            :obj:`str`: Single long string with the data table.
        """
        nrows, ncols = data_table.shape
        col_width = [ numpy.amax([ len(dij) for dij in dj]) for dj in data_table.T ]
        row_string = ['']*(nrows+1) if delimiter == 'print' else ['']*(nrows+3)
        start = 2 if delimiter == 'print' else 3
        for i in range(start,nrows+start-1):
            row_string[i] = '  '.join([ data_table[1+i-start,j].ljust(col_width[j]) 
                                                                        for j in range(ncols)])
        if delimiter == 'print':
            # Heading row
            row_string[0] = '  '.join([ data_table[0,j].ljust(col_width[j]) for j in range(ncols)])
            # Delimiter
            row_string[1] = '-'*len(row_string[0])
            return '\n'.join(row_string)+'\n'

        # For an rst table
        row_string[0] = '  '.join([ '='*col_width[j] for j in range(ncols)])
        row_string[1] = '  '.join([ data_table[0,j].ljust(col_width[j]) for j in range(ncols)])
        row_string[2] = row_string[0]
        row_string[-1] = row_string[0]
        return '\n'.join(row_string)+'\n'

    @staticmethod
    def _data_string(data, use_repr=True, verbatum=False):
        """
        Convert a single datum into a string

        Simply return strings, recursively convert the elements of
        any objects with a ``__len__`` attribute, or use the object's
        own ``__repr__`` attribute for all other objects.

        Args:
            data (object):
                The object to stringify.
        """
        if isinstance(data, str):
            return data if not verbatum else '``' + data + '``'
        if hasattr(data, '__len__'):
            return '[]' if isinstance(data, list) and len(data) == 0 \
                        else ', '.join([ ParSet._data_string(d, use_repr=use_repr,
                                                             verbatum=verbatum) for d in data ])
        if use_repr:
            return data.__repr__()
        return str(data)

    def _wrap_print(self, head, output, tcols):
        """
        Wrap the contents of an output string for a fixed terminal
        width.  Used for the long-format :func:`info` method.

        Args:
            head (:obj:`str`):
                The inline header for the output.  Can be an empty
                string, but cannot be ``None``.
            output (:obj:`str`):
                The main body of the text to write.
            tcols (:obj:`int`):
                The allowed width for the output.
        """
        tail = ' '*len(head)
        if tcols is not None:
            lines = textwrap.wrap('{0}'.format(output), tcols-len(head))
            if len(lines) == 0:
                print('{0}None'.format(head))
            else:
                _head = [ head ] + [ tail ]*(len(lines)-1)
                print('\n'.join([ h+l for h,l in zip(_head, lines)]))
        else:
            print(head+'{0}'.format(output))

    def _types_list(self, key):
        """Return the string names for the specified data types."""
        return ['Undefined' if t is None else t.__name__ for t in self.dtype[key]]

    @staticmethod
    def config_lines(par, section_name=None, section_comment=None, section_level=0,
                     exclude_defaults=False, include_descr=True):
        """
        Recursively generate the lines of a configuration file based on
        the provided ParSet or dict (par).

        Args:
            section_name (:obj:`str`, optional):
                Name to give to the top-level of the configuration
                output.
            section_comment (:obj:`str`, optional):
                Description to provide for the top-level configuration
                output.
            section_level (:obj:`int`, optional):
                The level for the configuration output.  Sets the
                indentation level and the number of square brackets
                assigned to the section name.
            exclude_defaults (:obj:`bool`, optional):
                Do not include any parameters that are identical to the
                defaults.
            include_descr (:obj:`bool`, optional):
                Include the descriptions of each parameter as comments.

        Returns:
            :obj:`list`: The list of the lines to write to a
            configuration file.
        """
        # Get the list of parameters that are ParSets
        parset_keys = [ k for k in par.keys() if isinstance(par[k], (ParSet, dict)) ]
        n_parsets = len(parset_keys)

        # Set the top-level comment and section name
        section_indent = ' '*4*section_level
        component_indent = section_indent + ' '*4
        lines = [] if section_comment is None \
                            else ParSet._config_comment(section_comment, section_indent)
        lines += [ section_indent + '['*(section_level+1) + section_name
                   + ']'*(section_level+1) ]

        min_lines = len(lines)

        # Add all the parameters that are not ParSets
        for k in par.keys():
            # Skip it if this element is a ParSet
            if n_parsets > 0 and k in parset_keys:
                continue

            # If the value is a list, determine if all the elements of
            # the list are also dictionaries or ParSets
            if isinstance(par[k], list) and len(par[k]) > 0:
                is_parset_or_dict = [ isinstance(v, (ParSet, dict)) for v in par[k] ]
                if numpy.all(is_parset_or_dict):
                    ndig = int(numpy.log10(len(par[k])))+1
                    for i, v in enumerate(par[k]):
                        indx = str(i+1).zfill(ndig)
                        # Try to add the section comment
                        section_comment = None
                        if include_descr:
                            try:
                                section_comment = par.descr[k] + ': ' + indx
                            except:
                                pass
                        lines += ParSet.config_lines(v, section_name=k+indx,
                                                     section_comment=section_comment,
                                                     section_level=section_level+1,
                                                     exclude_defaults=exclude_defaults,
                                                     include_descr=include_descr)
                    continue

            # Working with a single element
            # Try to add the description for this parameter
            try:
                if par.descr[k] is not None and include_descr:
                    lines += ParSet._config_comment(par.descr[k], component_indent)
            except:
                pass
            if not exclude_defaults or par[k] != par.default[k]:
                lines += [ component_indent + k + ' = ' + ParSet._data_string(par[k]) ]

        # Then add the items that are ParSets as subsections
        for k in parset_keys:
            section_comment = None
            if include_descr:
                try:
                    section_comment = par.descr[k]
                except:
                    pass
            lines += ParSet.config_lines(par[k], section_name=k, section_comment=section_comment,
                                         section_level=section_level+1,
                                         exclude_defaults=exclude_defaults,
                                         include_descr=include_descr)
        return lines if len(lines) > min_lines else []

    @staticmethod
    def _config_comment(comment, indent, full_width=72):
        """
        Create the list of lines for the description of a given
        parameter in the configuration file.

        Args:
            comment (:obj:`str`):
                The description of the configuration parameter.
            indent (:obj:`str`):
                The string used to indent the text.
            full_width (:obj:`int`, optional):
                The full width allowed for each output string in the
                returned list.

        Returns:
            :obj:`list`: List of the strings to write to the output
            configuration file.
        """
        head = indent + '# '
        lines = textwrap.wrap('{0}'.format(comment), full_width-len(head))
        return [ head + l for l in lines ]
   
    def info(self, basekey=None):
        """
        A long-form version of __repr__ that includes the parameter descriptions.
        """
        # Try to get the width of the available space to print
        try:
            tr, tcols = numpy.array(os.popen('stty size', 'r').read().split()).astype(int)
            tcols -= int(tcols*0.1)
        except:
            tr = None
            tcols = None

        for k in self.data.keys():
            if isinstance(self.data[k], ParSet):
                self.data[k].info(basekey=k)
                continue
            print('{0}'.format(k) if basekey is None else '{0}:{1}'.format(basekey,k))
            self._wrap_print('        Value: ', self.data[k], tcols)
            self._wrap_print('      Default: ', self.default[k], tcols)
            self._wrap_print('      Options: ', 'None' if self.options[k] is None
                                                else ', '.join(self.options[k]), tcols)
            self._wrap_print('  Valid Types: ', 'None' if self.dtype[k] is None
                                                else ', '.join(self._types_list(k)), tcols)
            self._wrap_print('     Callable: ', self.can_call[k], tcols)
            self._wrap_print('  Description: ', self.descr[k], tcols)
            print(' ')

    def keys(self):
        """Return the list of parameter set keys."""
        return list(self.data.keys())
    
    def add(self, key, value, default=None, options=None, dtype=None, can_call=None, descr=None):
        """
        Add a new parameter.

        Args:
            key (:obj:`str`):
                Key for new parameter
            value (object):
                Parameter value, must have a type in the list provided
                by ``dtype``, if the list is provided
            default (object, optional):
                Define a default value for the parameter, must have a
                type in the list provided by ``dtype``, if the list is
                provided.  No default if not provided.
            options (:obj:`list`, optional):
                List of discrete values that the parameter is allowed to
                have.  Allowed to be anything if not provided.
            dtype (:obj:`list`, optional):
                List of allowed data types that the parameter can have.
                Allowed to be anything if not provided.
            can_call (:obj:`bool`, optional):
                Flag that the parameters are callable operations.
                Default is False.
            descr (:obj:`str`, optional):
                Parameter description.  Default is that no description
                is added.

        Raises:
            ValueError:
                Raised if the keyword alread exists.
        """
        if key in self.data.keys():
            raise ValueError('Keyword {0} already exists and cannot be added!')
        self.npar += 1
        self.default[key] = None if default is None else default
        self.options[key] = [options] if options is not None and not isinstance(options, list) \
                                      else options
        self.dtype[key] = [dtype] if dtype is not None and not isinstance(dtype, list) else dtype
        self.can_call[key] = False if can_call is None else can_call
        self.descr[key] = None if descr is None else descr
        self.data[key] = None
        try:
            self.__setitem__(key, value)
        except:
            # Delete the added components
            del self.default[key]
            del self.options[key]
            del self.dtype[key]
            del self.can_call[key]
            del self.descr[key]
            # Re-raise the exception
            raise

    def to_config(self, cfg_file=None, section_name=None, section_comment=None, section_level=0,
                  append=False, quiet=False, exclude_defaults=False, include_descr=True):
        """
        Write/Append the parameter set to a configuration file.

        Args:
            cfg_file (:obj:`str`, optional):
                The name of the file to write/append to.  If None
                (default), the function will just return the list of
                strings that would have been written to the file.  These
                lines can be used to construct a `configobj.ConfigObj`_
                instance.
            section_name (:obj:`str`, optional):
                The top-level name for the config section.  This must be
                provided if :attr:`cfg_section` is None or any of the
                parameters are not also :class:`ParSet` instances
                themselves.
            section_comment (:obj:`str`, optional):
                The top-level comment for the config section based on
                this :class:`ParSet`.
            section_level (:obj:`int`, optional):
                The top level of this :class:`ParSet`.  Used for
                recursive output of nested :class:`ParSet` objects.
            append (:obj:`bool`, optional):
                Append this configuration output of this :class:`ParSet`
                to the file.  False by default.  If not appending and
                the file exists, the file is automatically overwritten.
            quiet (:obj:`bool`, optional):
                Suppress all standard output from the function.
            exclude_defaults (:obj:`bool`, optional):
                Do not include any parameters that are identical to the
                defaults.
            include_descr (:obj:`bool`, optional):
                Include the descriptions of each parameter as comments.

        Raises:
            ValueError:
                Raised if there are types other than :class:`ParSet` in
                the parameter list, :attr:`cfg_section` is ``None``, and
                no section_name argument was provided.
        """
        if cfg_file is not None and os.path.isfile(cfg_file) and not append and not quiet:
            warnings.warn('Selected configuration file already exists and will be overwritten!')

        config_output = []
        if numpy.all([ isinstance(d, ParSet) or d is None for d in self.data.values() ]):
            # All the elements are ParSets themselves, so just iterate
            # through each one
            for k in self.keys():
                if self.data[k] is None:
                    continue
                section_comment = self.descr[k] if include_descr else None
                config_output += ParSet.config_lines(self.data[k], section_name=k,
                                                     section_comment=section_comment,
                                                     section_level=section_level,
                                                     exclude_defaults=exclude_defaults,
                                                     include_descr=include_descr)
        else:
            # Cannot write the parameters as a configuration file
            # without a top-level configuration section
            if section_name is None and self.cfg_section is None:
                warnings.warn('No top-level section name available; using [default].')
                section_name = 'default'

            _section_name = self.cfg_section if section_name is None else section_name
            _section_comment = self.cfg_comment if section_comment is None else section_comment
            config_output += ParSet.config_lines(self, section_name=_section_name,
                                                 section_comment=_section_comment,
                                                 section_level=section_level,
                                                 exclude_defaults=exclude_defaults,
                                                 include_descr=include_descr)

        if cfg_file is None:
            # Only return the list of lines for the output file.  Useful
            # if you want to use instantly create a new ConfigObj
            # instance without having to write a file
            return config_output

        # Write the file
        with open(cfg_file, 'a' if append else 'w') as f:
            f.write('\n'.join(config_output))

    @classmethod
    def from_config(cls, cfg, section_name='default', evaluate=True):
        """
        Construct the parameter set using a configuration file.

        Args:
            cfg (:obj:`str`, :obj:`list`):
                Either a single string with a file name to read, or a
                list of configuration-file-style strings with the
                parameters.
            section_name (:obj:`str`, optional):
                The configuration file section with the parameters.
            evaluate (:obj:`bool`, optional):
                Evaluate the values in the config object before
                assigning them in the subsequent parameter sets.  The
                parameters in the config file are *always* read as
                strings, so this should almost always be true; however,
                see the warning below.

                .. warning::

                    When ``evaluate`` is true, the function runs
                    ``eval()`` on all the entries in the `ConfigObj`
                    dictionary, done using
                    :func:`mangadap.par.util.recursive_dict_evaluate`.
                    This has the potential to go haywire if the name of
                    a parameter unintentionally happens to be identical
                    to an imported or system-level function.  Of course,
                    this can be useful by allowing one to define the
                    function to use as a parameter, but it also means
                    one has to be careful with the values that the
                    parameters should be allowed to have.  The current
                    way around this is to provide a list of strings that
                    should be ignored during the evaluation, done using
                    :func:`mangadap.par.util._eval_ignore`.
                
        Returns:
            :class:`ParSet`: The instance of the parameter set.
        """
        # Instantiate the ConfigObj instance
        _cfg = ConfigObj(cfg)[section_name]

        # Evaluate the strings, if requested
        if evaluate:
            _cfg = recursive_dict_evaluate(_cfg)
        
        # Instantiate the object based on the configuration dictionary
        return cls.from_dict(_cfg)

    @staticmethod
    def _rst_class_name(p):
        return ':class:`' +  type(p).__module__ + '.' + type(p).__name__ + '`'

    def to_rst_table(self, parsets_listed=[], header=True, class_link=True):
        r"""
        Construct a reStructuredText table with the :class:`ParSet`
        data.

        This method is mostly meant for documentation purposes, as way
        of showing the format and default parameters of a given single
        :class:`ParSet` or nested set of :class:`ParSet` objects.

        Args:
            parsets_listed (:obj:`list`, optional):
                For nested :class:`ParSet` objects, this is used to keep
                a log of :class:`ParSet` objects that have already been
                included in the rst table, forcing the table to only
                appear once.  
            header (:obj:`bool`, optional):
                Include a section header
            class_link (:obj:`bool`, optional):
                Include an rst-style link to the class instantiation
                documentation.

        Returns:
            :obj:`list`: A list of strings containing each line of the
            rst table.  To print the table::

                print('\n'.join(p.to_rst_table()))
            
            where ``p`` is a :class:`ParSet` instance.
        """
        new_parsets = []
        data_table = numpy.empty((self.npar+1, 5), dtype=object)
        data_table[0,:] = ['Key', 'Type', 'Options', 'Default', 'Description']
        for i,k in enumerate(self.keys()):
            data_table[i+1,0] = ParSet._data_string(k, use_repr=False, verbatum=True)
            if isinstance(self.data[k], ParSet):
                if type(self.data[k]).__name__ not in parsets_listed:
                    new_parsets += [k]
                parsets_listed += [ type(self.data[k]).__name__ ]
                data_table[i+1,1] = ParSet._rst_class_name(self.data[k])
                data_table[i+1,3] = '`{0} Keywords`_'.format(type(self.data[k]).__name__)
            else: 
                data_table[i+1,1] = ', '.join(self._types_list(k))
                data_table[i+1,3] = '..' if self.default[k] is None \
                                    else ParSet._data_string(self.default[k], use_repr=False,
                                                             verbatum=True)

            data_table[i+1,2] = '..' if self.options[k] is None \
                                    else ParSet._data_string(self.options[k], use_repr=False,
                                                             verbatum=True)
            data_table[i+1,4] = '..' if self.descr[k] is None \
                                    else ParSet._data_string(self.descr[k])

        output = ['']
        if header:
            output += ['{0} Keywords'.format(type(self).__name__)]
            output += ['-'*len(output[0])]
            output += ['']
        if class_link:
            output += ['Class Instantiation: ' + ParSet._rst_class_name(self)]
            output += ['']
        output += [ParSet._data_table_string(data_table, delimiter='rst')]
        output += ['']
        for k in new_parsets:
            output += ['----']
            output += ['']
            output += self.data[k].to_rst_table(parsets_listed=parsets_listed)

        return output

    def validate_keys(self, required=None, can_be_None=None):
        """

        Validate the keys in the :class:`ParSet`.

        Args:
            required (:obj:`list`, optional):
                A list of required keys.
            can_be_None (:obj:`list`, optional):
                A list of keys with values that are allowed to be None.
                All other keys are expected to have defined values.
        
        Raises:
            ValueError:
                Raised if required keys are not present or if keys have
                ``None`` values and are not expected to be.
        """
        if required is None and can_be_None is None:
            # No validation rules, so implicitly valid
            return

        if required is not None:
            not_defined = numpy.array([ k not in self.keys() for k in required ])
            if numpy.any(not_defined):
                raise ValueError('Required keys were not defined: {0}'.format(
                                    numpy.asarray(required)[not_defined].tolist()))

        if can_be_None is not None:
            should_not_be_None = numpy.array([ self.data[k] is None and k not in can_be_None 
                                                                    for k in self.keys()])
            if numpy.any(should_not_be_None):
                raise ValueError('These keys should not be None: {0}'.format(
                                    numpy.asarray(self.keys())[should_not_be_None].tolist()))

    def to_dict(self):
        """
        Return a dictionary with the parameters.

        .. warning::

            This simply returns a pointer to the internal object
            dictionary, :attr:`data`.

        """
        # TODO: Return a copy?
        return self.data

    @classmethod
    def from_dict(cls, data):
        """
        Create a :class:`ParSet` from a dictionary.
        
        Objects built in this way are nearly identical to a normal
        dictionary, except that one cannot add keys in the same way.
        """
        return cls([*data.keys()], values=[*data.values()])

    def to_header(self, hdr, prefix=None, quiet=False):
        """
        Write the parameters to a fits header.

        Any element that has a value of ``None`` or is a :class:`ParSet`
        itself is *not* written to the header.

        Args:
            hdr (`astropy.io.fits.Header`_):
                Header object for the parameters. **Modified in-place**.
            prefix (:obj:`str`, optional):
                Prefix to use for the header keywords, which overwrites
                the string defined for the class. If ``None``, uses the
                default for the class.
            quiet (:obj:`bool`, optional):
                Suppress print statements.
        """
        if prefix is None:
            prefix = self.prefix
        ndig = int(numpy.log10(self.npar))+1 
        for i, (key, value) in enumerate(self.data.items()):
#            if value is None:
#                # Don't write Nones
#                continue
            if isinstance(value, ParSet):
                if verbose:
                    warnings.warn('ParSets within ParSets are not written to headers!  '
                                  'Skipping {0}.'.format(key))
                continue
            _value = str(value) if value is None or isinstance(value, (list, tuple)) else value
            hdr['{0}{1}'.format(prefix, str(i+1).zfill(ndig))] \
                    = (_value, '{0}: {1}'.format(self.__class__.__name__, key))

    @classmethod
    def from_header(cls, hdr, prefix=None):
        """
        Instantiate the :class:`ParSet` using data parsed from a fits
        header.

        Args:
            hdr (`astropy.io.fits.Header`_):
                Header object with the parameters.
            prefix (:obj:`str`, optional):
                Prefix of the relevant header keywords, which overwrites
                the string defined for the class. If None, uses the
                default for the class.
        """
        if prefix is None:
            prefix = cls.prefix
        return cls.from_dict(recursive_dict_evaluate(ParSet.parse_par_from_hdr(hdr, prefix)))

    @staticmethod
    def parse_par_from_hdr(hdr, prefix):
        """
        Parse the dictionary of parameters written to a header.

        Args:
            hdr (`astropy.io.fits.Header`_):
                Header object to parse.
            prefix (:obj:`str`):
                The prefix used for the header keywords.
        
        Returns:
            :obj:`dict`: A dictionary with the parameter keywords and
            values.
        """
        par = {}
        for k, v in hdr.items():
            # Check if this header keyword starts with the required
            # prefix
            if k[:len(prefix)] == prefix:
                try:
                    # Try to convert the keyword without the prefix into
                    # an integer.
                    i = int(k[len(prefix):])-1
                except ValueError:
                    # Assume the value is some other random keyword that
                    # starts with the prefix but isn't a parameter
                    continue
                # Assume we've found a parameter. Parse the parameter
                # name from the header comment and add it to the
                # dictionary.
                par_key = hdr.comments[k].split(':')[-1].strip()
                par[par_key] = v
        return par


class KeywordParSet(ParSet):
    """
    An abstract class that uses :class:`ParSet` as its base.

    The main purpose of this class is to redefine the
    :func:`ParSet.from_dict` method and disallow adding new parameters.
    """
    @classmethod
    def from_dict(cls, data, ignore_extra=True):
        """
        Construct the object using a dictionary.
        """
        k = numpy.array([*data.keys()])
        parkeys = [*inspect.signature(cls).parameters.keys()]

        if not ignore_extra:
            badkeys = numpy.array([pk not in parkeys for pk in k])
            if numpy.any(badkeys):
                raise ValueError('{0} not recognized key(s) for {1}.'.format(k[badkeys],
                                 cls.__class__.__name__))

        kwargs = {}
        for pk in parkeys:
            kwargs[pk] = data[pk] if pk in k else None
        return cls(**kwargs)

    def add(self, *args, **kwargs):
        """
        Disallow functionality of base class that adds new parameters.
        """
        raise NotImplementedError('Cannot add parameters to a {0} instance.'.format(
                                  self.__class__.__name__))


# TODO: Change this to a DataTable?
class ParDatabase:
    """
    Class used as a list of ParSets in a glorified structured numpy
    array.

    Very similar to yanny when converted to a numpy array.

    Can be initialized using a list of ParSet objects, or an SDSS
    parameter file.

    .. todo::

        - Check that the data types are the same for all ParSet objects
          in the list
        - Better handle the NaN values when converting None to a float
          type
        - Add from_par_file classmethod?
    
    """
    def __init__(self, inp):
        """
        nsets - number of parameter sets
        npar - number of parameters in each parameter set
        data - parameter set data
        options - allowed options for values
        dtype - allowed datatypes
        can_call - parameter is a callable function
        """
        _inp = [inp] if isinstance(inp, ParSet) else inp
        if not isinstance(_inp, list):
            raise TypeError('Input must be a list.')
        for i in _inp:
            if not isinstance(i, ParSet):
                raise TypeError('Input must be a list of ParSet objects.')
        self.npar = _inp[0].npar
        self.nsets = len(_inp)
        keys = _inp[0].keys()
        for i in range(1,self.nsets):
            if _inp[i].npar != self.npar:
                raise ValueError('Not all ParSet objects have the same number of parameters.')
            if _inp[i].keys() != keys:
                raise ValueError('Not all ParSet objects have the same keys.')
            # Other checks?

        record_dtypes = self._set_dtypes(_inp, 0)

        data = []
        for i in range(self.nsets):
            data += [ tuple([_inp[i][k] for k in keys]) ]

        # WARNING: None values are converted to nan if data type is
        # float
        self.data = numpy.array(data, dtype=record_dtypes ).view(numpy.recarray)
        self.options = inp[0].options.copy()
        self.dtype = inp[0].dtype.copy()
        self.can_call = inp[0].can_call.copy()
   
    def __getitem__(self, key):
        """
        Return the value of the designated key.

        Args:
            key (str) : Key for new parameter
        """
        return self.data[key]

    def __len__(self):
        return self.nsets

    @staticmethod
    def _set_dtypes(inp, i):
        keys = inp[i].keys()
        dtypes = []
        for k in keys:
            if inp[i].dtype[k] is None:
                dtypes += [(k,object)]
                continue
            # inp.dtype is always a list
            if any([t in inp[i].dtype[k] for t in [int , float]]) \
                and any([t in inp[i].dtype[k] for t in [list, numpy.ndarray]]):
                warnings.warn('Parameter set has elements that can be either individual ' \
                              'ints/floats or lists/arrays.  Database column {0} will have type ' \
                              '\'object\'.'.format(k))
                dtypes += [(k,object)]
            elif len(list({int, float} - set(inp[i].dtype[k]))) == 0:
                dtypes += [(k,float)]
            elif len(list({list, numpy.ndarray} - set(inp[i].dtype[k]))) == 0 \
                    or inp[i].dtype[k] == numpy.ndarray:
                _inp = numpy.asarray(inp[i][k])
                dtypes += [(k,_inp.dtype,_inp.shape)]
            elif isinstance(inp[i][k], str):
                if any([ _inp[k] is None for _inp in inp]):
                    dtypes += [(k, object)]
                else:
                    dtypes += [(k,'<U{0:d}'.format(max([ len(_inp[k]) for _inp in inp])))]
            else:
                dtypes += [(k,type(inp[i][k]))]
        return dtypes


    def append(self, db):
        if not isinstance(db, ParDatabase):
            raise TypeError('Can only append ParDatabase object.')

        try:
            self.data = numpy.append(self.data, db.data)
        except TypeError as e:
            raise TypeError('Could not append data:: {0}'.format(e))
            



