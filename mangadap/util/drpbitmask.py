r"""
Defines a :class:`~mangadap.util.bitmask.BitMask` objects for use in the DAP
based on the definition from the MaNGA DRP.


----

.. include license and copyright
.. include:: ../include/copy.rst

----

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
from IPython import embed

from .bitmask import BitMask
from ..config import defaults


class DRPFitsBitMask(BitMask):
    r"""
    Structure with the DRP mask bits.

    The defined bits are listed at :ref:`metadatamodel-drp3pixmask`.

    Mode is not checked
    """
    def __init__(self, sdss_maskbits=None, mode='CUBE'):
        _sdss_maskbits = defaults.sdss_maskbits_file() if sdss_maskbits is None else sdss_maskbits
        tmp = BitMask.from_par_file(str(_sdss_maskbits),
                                    'MANGA_DRP3PIXMASK' if mode == 'CUBE' else 'MANGA_DRP2PIXMASK')
        keys, descr = tmp._init_objs()
        super(DRPFitsBitMask, self).__init__(keys, descr=descr)


class DRPQuality3DBitMask(BitMask):
    r"""
    Structure with the definition of the DRP3QUAL mask bits.

    The defined bits are listed at :ref:`metadatamodel-drp3qual`.
    """
    def __init__(self, sdss_maskbits=None):
        _sdss_maskbits = defaults.sdss_maskbits_file() if sdss_maskbits is None else sdss_maskbits
        tmp = BitMask.from_par_file(str(_sdss_maskbits), 'MANGA_DRP3QUAL')
        keys, descr = tmp._init_objs()
        super(DRPQuality3DBitMask, self).__init__(keys, descr=descr)


