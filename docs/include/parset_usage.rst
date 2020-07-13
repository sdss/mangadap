
Parameter sets are used to control the behavior of internal methods from
top-level input.  Parameter sets:

    * provide I/O methods for inspection, documentation, and logging
    * force parameters to have specific data types, including if the
      parameter should be callable
    * specify allowed options for method keywords
    * provide parameter descriptions

Beyond this and once instantiated, though, they are similar to a normal
python dictionary.

ParSet base class
-----------------

As an example, we can create a bogus parameter set as follows:

.. code-block:: python

    from mangadap.par.parset import ParSet
    pars = ['test', 'par', 'lst', 'junk']
    defaults = ['this', 0, 0.0, None]
    options = [['this', 'that'], None, None, None]
    dtypes = [str, int, [int, float], None]
    descr = ['test parameter description',
             'par parameter description',
             None,
             'this is a rather long description of the junk parameter.  You can include ' \
                'rst-style references like pointing back to the ' \
                ':class:`~mangadap.par.parset.ParSet` class, for when this description is ' \
                'written to an rst table using :func:`~mangadap.par.parset.ParSet.to_rst_table` ' \
                'and included in an rst doc synthesized into html using sphinx.']
    p = ParSet(pars, defaults=defaults, options=options, dtypes=dtypes, descr=descr)

Once defined, accessing/showing the values is done similarly to a
dictionary, with some convenience printing methods::

    >>> p['test']
    'that'
    >>> p['lst']
    0.0
    >>> p['junk'] is None
    True
    >>> p
    Parameter  Value  Default  Type        Callable
    -----------------------------------------------
    test       that   this     str         False
    par        3      0        int         False
    lst        0.0    0.0      int, float  False
    junk       None   None     Undefined   False
    >>> p.info()
    test
            Value: this
          Default: this
          Options: this, that
      Valid Types: str
         Callable: False
      Description: test parameter description

    par
            Value: 0
          Default: 0
          Options: None
      Valid Types: int
         Callable: False
      Description: par parameter description

    lst
            Value: 0.0
          Default: 0.0
          Options: None
      Valid Types: int, float
         Callable: False
      Description: None

    junk
            Value: None
          Default: None
          Options: None
      Valid Types: Undefined
         Callable: False
      Description: this is a rather long description of the junk parameter.
                   You can include rst-style references like pointing back
                   to the :class:`~mangadap.par.parset.ParSet` class, for
                   when this description is written to an rst table using
                   :func:`~mangadap.par.parset.ParSet.to_rst_table` and
                   included in an rst doc synthesized into html using
                   sphinx.
    
Restrictions are placed on the allowed values and types for the
parameters and new keys cannot be added without doing so explicitly::

    >>> p['test'] = 'the other'
    ...
    ValueError: Input value for test invalid: the other.
    Options are: ['this', 'that']
    >>> p['par'] = 'foo'
    ...
    TypeError: Input value for par has incorrect type: foo.
    Valid types are: [<class 'int'>]
    >>> p['new'] = 8.3
    ...
    KeyError: 'new is not a valid key for ParSet.'
    >>> p.add('new', 8.3, dtype=float)
    >>> p
    Parameter  Value  Default  Type        Callable
    -----------------------------------------------
    test       that   this     str         False
    par        3      0        int         False
    lst        0.0    0.0      int, float  False
    junk       None   None     Undefined   False
    new        8.3    None     float       False
    >>> p['new'] = 8
    ...
    TypeError: Input value for new has incorrect type: 8.
    Valid types are: [<class 'float'>]
    >>> p['new'] = 8.
    >>> p['new']
    8.0

There are also a number of IO methods:

    - To convert to or instantiate from a dictionary::

        >>> p.to_dict()
        {'test': 'that', 'par': 3, 'lst': 0.0, 'junk': None}
        >>> p.data
        {'test': 'that', 'par': 3, 'lst': 0.0, 'junk': None}
        >>> p.to_dict() is p.data
        True
        >>> ParSet.from_dict(p.to_dict())
        Parameter  Value  Default  Type       Callable
        ----------------------------------------------
        test       that   None     Undefined  False
        par        3      None     Undefined  False
        lst        0.0    None     Undefined  False
        junk       None   None     Undefined  False

    - To write to or read from an `astropy.io.fits.Header`_::

        >>> from astropy.io import fits
        >>> hdr = fits.Header()
        >>> p.to_header(hdr)
        >>> hdr
        PAR1    = 'that    '           / ParSet: test
        PAR2    =                    3 / ParSet: par
        PAR3    =                  0.0 / ParSet: lst
        PAR4    = 'None    '           / ParSet: junk
        >>> ParSet.from_header(hdr)
        Parameter  Value  Default  Type       Callable
        ----------------------------------------------
        test       that   None     Undefined  False
        par        3      None     Undefined  False
        lst        0.0    None     Undefined  False
        junk       None   None     Undefined  False

    - To write to or read from a configuration file::

        >>> print('\n'.join(p.to_config()))
        [default]
            # test parameter description
            test = that
            # par parameter description
            par = 3
            lst = 0.0
            # this is a rather long description of the junk parameter.  You can
            # include rst-style references like pointing back to the
            # :class:`~mangadap.par.parset.ParSet` class, for when this
            # description is written to an rst table using
            # :func:`~mangadap.par.parset.ParSet.to_rst_table` and included in an
            # rst doc synthesized into html using sphinx.
            junk = None
        >>> ParSet.from_config(p.to_config())
        Parameter  Value  Default  Type       Callable
        ----------------------------------------------
        test       that   None     Undefined  False
        par        3      None     Undefined  False
        lst        0.0    None     Undefined  False
        junk       None   None     Undefined  False

Note that in all of the IO methods above, the instantiation method loses
essentially all of the differences between the
:class:`~mangadap.par.parset.ParSet` and a normal dictionary.  For this
and other reasons, we've implemented an abstract class called
:class:`~mangadap.par.parset.KeywordParSet`.

KeywordParSet class
-------------------

The :class:`~mangadap.par.parset.KeywordParSet` class is derived from
:class:`~mangadap.par.parset.ParSet` and does two things:

    1. overwrites the :func:`~mangadap.par.parset.ParSet.add` method so
       that no new parameters can be added and

    2. overwrites the :func:`~mangadap.par.parset.ParSet.from_dict`
       method with the expectation that any class derived from
       :class:`~mangadap.par.parset.KeywordParSet` has an ``__init__``
       method that takes a fixed set of keyword arguments.

By overwriting the base class definition,
:func:`~mangadap.par.parset.KeywordParSet.from_dict` takes care of all
of the other "from" methods because they in turn use this "from_dict"
method to instantiate the object.

All of the parameter-set classes defined and used by the DAP use
:class:`~mangadap.par.parset.KeywordParSet` as their base. We can
rewrite the :class:`~mangadap.par.parset.ParSet` example above to use
this new base class and construct a relevant demonstration class:

.. code-block:: python

    from mangadap.par.parset import KeywordParSet

    class DemoPar(KeywordParSet):
        def __init__(self, test=None, par=None, lst=None, junk=None):

            pars = ['test', 'par', 'lst', 'junk']
            values = [test, par, lst, junk]
            defaults = ['this', 0, 0.0, None]
            options = [ ['this', 'that'], None, None, None ]
            dtypes = [ str, int, [int, float], None ]
            descr = ['test parameter description',
                     'par parameter description',
                     None,
                     'this is a rather long description of the junk parameter.  You can include ' \
                        'rst-style references like pointing back to the ' \
                        ':class:`~mangadap.par.parset.ParSet` class, for when this description ' \
                        'is written to an rst table using ' \
                        ':func:`~mangadap.par.parset.ParSet.to_rst_table` ' \
                        'and included in an rst doc synthesized into html using sphinx.']
            super(DemoPar, self).__init__(pars, values=values, defaults=defaults, options=options,
                                          dtypes=dtypes, descr=descr)


The instantiation method for the derived class looks nearly identical
to how we originally defined the :class:`~mangadap.par.parset.ParSet`
instance. However, we can now define the instance using keyword
arguments directly, and the ancillary information is now propagated
to all the IO methods::

    >>> p = DemoPar(par=3, test='that')
    >>> p
    Parameter  Value  Default  Type        Callable
    -----------------------------------------------
    test       that   this     str         False
    par        3      0        int         False
    lst        0.0    0.0      int, float  False
    junk       None   None     Undefined   False
    >>> p['test'] = 'the other'
    ...
    ValueError: Input value for test invalid: the other.
    Options are: ['this', 'that']
    >>> p.add('new', 8.3, dtype=float)
    ...
    NotImplementedError: Cannot add parameters to a DemoPar instance.
    >>> DemoPar.from_dict(p.to_dict())
    Parameter  Value  Default  Type        Callable
    -----------------------------------------------
    test       that   this     str         False
    par        3      0        int         False
    lst        0.0    0.0      int, float  False
    junk       None   None     Undefined   False
    >>> from astropy.io import fits
    >>> hdr = fits.Header()
    >>> p.to_header(hdr)
    >>> DemoPar.from_header(hdr)
    Parameter  Value  Default  Type        Callable
    -----------------------------------------------
    test       that   this     str         False
    par        3      0        int         False
    lst        0.0    0.0      int, float  False
    junk       None   None     Undefined   False
    >>> DemoPar.from_config(p.to_config())
    Parameter  Value  Default  Type        Callable
    -----------------------------------------------
    test       that   this     str         False
    par        3      0        int         False
    lst        0.0    0.0      int, float  False
    junk       None   None     Undefined   False



