import pytest

from collections import OrderedDict

import warnings
warnings.simplefilter("ignore", UserWarning)
warnings.simplefilter("ignore", RuntimeWarning)

import numpy

from astropy.io import fits

from mangadap.util.bitmask import BitMask
from .util import data_test_file

#-----------------------------------------------------------------------------

class ImageBitMask(BitMask):
    def __init__(self):
        bits = OrderedDict([('BPM', 'Pixel is part of a bad-pixel mask'),
                            ('COSMIC', 'Pixel is contaminated by a cosmic ray'),
                            ('SATURATED', 'Pixel is saturated.')])
        super(ImageBitMask, self).__init__(list(bits.keys()), descr=list(bits.values()))


class ImageBitMaskFromPar(BitMask):
    def __init__(self):
        tmp = BitMask.from_par_file(data_test_file('imagebitmask.par'), 'IMAGEMASK')    
        super(ImageBitMaskFromPar, self).__init__(tmp.keys(), descr=tmp.descr)


def test_new():
    image_bm = ImageBitMask()
    assert list(image_bm.bits.keys()) == ['BPM', 'COSMIC', 'SATURATED']
    assert image_bm.keys() == ['BPM', 'COSMIC', 'SATURATED']
    assert list(image_bm.bits.values()) == [0, 1, 2]


def test_flagging():

    n = 1024
    shape = (n,n)

    image_bm = ImageBitMask()
    mask = numpy.zeros(shape, dtype=image_bm.minimum_dtype())

    cosmics_indx = numpy.zeros(shape, dtype=bool)
    cosmics_indx[numpy.random.randint(0,high=n,size=100),
                 numpy.random.randint(0,high=n,size=100)] = True
    mask[cosmics_indx] = image_bm.turn_on(mask[cosmics_indx], 'COSMIC')

    saturated_indx = numpy.zeros(shape, dtype=bool)
    saturated_indx[numpy.random.randint(0,high=n,size=100),
                   numpy.random.randint(0,high=n,size=100)] = True
    mask[saturated_indx] = image_bm.turn_on(mask[saturated_indx], 'SATURATED')

    assert numpy.sum(image_bm.flagged(mask, flag='BPM')) == 0
    assert numpy.sum(image_bm.flagged(mask, flag='COSMIC')) == numpy.sum(cosmics_indx)
    assert numpy.sum(image_bm.flagged(mask, flag='SATURATED')) == numpy.sum(saturated_indx)

    assert image_bm.flagged_bits(1) == ['BPM']
    assert image_bm.flagged_bits(2) == ['COSMIC']
    assert image_bm.flagged_bits(4) == ['SATURATED']

    unique_flags = numpy.sort(numpy.unique(numpy.concatenate(
                        [image_bm.flagged_bits(b) for b in numpy.unique(mask)]))).tolist()
    assert unique_flags == ['COSMIC', 'SATURATED']

    mask[saturated_indx] = image_bm.turn_off(mask[saturated_indx], 'SATURATED')
    assert numpy.sum(image_bm.flagged(mask, flag='COSMIC')) == numpy.sum(cosmics_indx)
    assert numpy.sum(image_bm.flagged(mask, flag='SATURATED')) == 0

    unique_flags = numpy.sort(numpy.unique(numpy.concatenate(
                        [image_bm.flagged_bits(b) for b in numpy.unique(mask)]))).tolist()
    assert unique_flags == ['COSMIC']

    b_indx, c_indx, s_indx = image_bm.unpack(mask)
    assert numpy.sum(b_indx) == 0
    assert numpy.sum(c_indx) == numpy.sum(cosmics_indx)
    assert numpy.sum(s_indx) == 0


def test_hdr_io():
    
    image_bm = ImageBitMask()
    hdr = fits.Header()
    image_bm.to_header(hdr)

    assert list(hdr.keys()) == ['BIT0', 'BIT1', 'BIT2']
    assert list(hdr.values()) == ['BPM', 'COSMIC', 'SATURATED']

    bm = BitMask.from_header(hdr)
    assert list(bm.bits.keys()) == ['BPM', 'COSMIC', 'SATURATED']
    assert list(bm.bits.values()) == [0, 1, 2]


def test_par_io():
    
    bm = BitMask.from_par_file(data_test_file('imagebitmask.par'), 'IMAGEMASK')
    assert list(bm.keys()) == ['BPM', 'COSMIC', 'SATURATED']
    assert list(bm.bits.values()) == [0, 1, 2]

    bm = ImageBitMaskFromPar()
    assert list(bm.keys()) == ['BPM', 'COSMIC', 'SATURATED']
    assert list(bm.bits.values()) == [0, 1, 2]


def test_ini_io():
    
    bm = BitMask.from_ini_file(data_test_file('imagebitmask.ini'))
    assert list(bm.bits.keys()) == ['BPM', 'COSMIC', 'SATURATED']
    assert list(bm.bits.values()) == [0, 1, 2]


def test_rst():

    bm = ImageBitMask()
    lines = bm.to_rst_table()
    assert lines[0] == 'ImageBitMask Bits'


