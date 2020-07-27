# -*- coding: utf-8 -*-
"""
Utility functions for parameter sets.

----

.. include license and copyright
.. include:: ../include/copy.rst

----

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import numpy 


def _eval_ignore():
    """Provides a list of strings that should not be evaluated."""
    return [ 'open', 'file', 'dict', 'list', 'tuple' ]


def recursive_dict_evaluate(d):
    """
    Recursively run `eval`_ on each element of the provided
    dictionary.

    A raw read of a configuration file with `configobj.ConfigObj`_
    results in a dictionary that contains strings or lists of strings.
    However, when assigning the values for the various
    :class:`mangadap.par.parset.ParSet` objects, the
    :func:`mangadap.par.parset.ParSet.from_dict` methods expect the
    dictionary values to have the appropriate type.  E.g., the
    `configobj.ConfigObj`_ will have something like ``d['foo'] = '1'``,
    when the :func:`mangadap.par.parset.ParSet.from_dict` method expects
    the value to be an integer (``d['foo'] = 1``).

    This function tries to evaluate *all* dictionary values, except for
    those listed by the :func:`_eval_ignore` function.  Any value in
    this list or where::

        eval(d[k]) for k in d.keys()

    raises an exception and is returned as the original string.

    Args:
        d (:obj:`dict`):
            Dictionary of values to evaluate

    Returns:
        :obj:`dict`: Identical to input dictionary, but with all string
        values replaced with the result of ``eval(d[k])`` for all ``k``
        in ``d.keys()``.
    """
    ignore = _eval_ignore()
    for k in d.keys():
        if isinstance(d[k], dict):
           d[k] = recursive_dict_evaluate(d[k])
        elif isinstance(d[k], list):
            replacement = []
            for v in d[k]:
                if v in ignore:
                    replacement += [ v ]
                else:
                    try:
                        replacement += [ eval(v) ]
                    except:
                        replacement += [ v ]
            d[k] = replacement
        else:
            try:
                d[k] = eval(d[k]) if d[k] not in ignore else d[k]
            except:
                pass

    return d


def get_parset_list(cfg, pk, parsetclass):
    """
    Create a list of :class:`mangadap.par.parset.ParSet` objects based
    on a root keyword for a set of defined groups in the configuration
    file.
    
    For example, the a parameter group could allow for a list of parset
    objects instead of a single one.  This function parses the provided
    configuration object (``cfg``) to find any sections with a string
    keyword (`pk`) as its root.  The remainder of the section name must
    be able to be converted to an integer and the section itself must be
    able to setup an instance of `parsetclass`.  The sections must be
    number sequentially from 1..N.  

    Args:
        cfg (`configobj.ConfigObj`_, :obj:`dict`):
            The top-level configuration that defines a list of
            sub-ParSets.
        pk (:obj:`str`):
            The root of the keywords used to set a list of sub-ParSets.
        parsetclass (:class:`mangadap.par.parset.ParSet`):
            The class used to construct each element in the list of
            parameter subsets.  The class **must** have a ``from_dict``
            method that instantiates the
            :class:`mangadap.par.parset.ParSet` based on the provided
            subsection/subdict from ``cfg``.

    Returns:
        :obj:`list`: A list of instances of ``parsetclass`` parsed from
        the provided configuration data.

    Raises:
        ValueError:
            Raised if the indices of the subsections are not sequential
            and 1-indexed.
    """
    # Get the full list of keys
    k = cfg.keys()

    # Iterate through the list of keys to find the appropriate sub
    # parameter sets and their order.
    par = []
    order = []
    for _k in k:
        if _k == pk and cfg[_k] is None:
            continue
        if pk in _k:
            try:
                # Get the order for this subgroup (e.g., 2 for
                # 'detector2'
                order += [ int(_k.replace(pk,'')) ]
                # And instantiate the parameter set
                par += [ parsetclass.from_dict(cfg[_k]) ]
            except:
                continue

    if len(par) > 0:
        # Make sure the instances are correctly sorted and sequential
        srt = numpy.argsort(order)
        if numpy.any(numpy.array(order)[srt]-1 != numpy.arange(order[srt[-1]])):
            raise ValueError('Parameter set series must be sequential and 1-indexed.')
        # Return the sorted instances
        return [par[i] for i in srt]

    # No such subsets were defined, so return a null result
    return None


