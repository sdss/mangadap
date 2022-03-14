#!/usr/bin/env python3

import numpy
import os
import warnings
import time
import argparse

from IPython import embed

from astropy.io import fits
from astropy import constants

from mangadap.util.bitmask import BitMask
from mangadap.config import defaults

from mangadap.scripts import scriptbase


def double_print(ostream, line, **kwargs):
    if ostream is not None:
        print(line, file=ostream, **kwargs)
    print(line, **kwargs)


def find_repeat_observations(dapall_file, ofile=None):

    dapall = fits.open(dapall_file)[1].data

#    # Get the maskbits file
#    sdssMaskbits = defaults.sdss_maskbits_file()
#
#    # Instantiate the BitMask object
#    mngtarg1_bm = BitMask.from_par_file(sdssMaskbits, 'MANGA_TARGET1')
#    mngtarg3_bm = BitMask.from_par_file(sdssMaskbits, 'MANGA_TARGET3')

    cz = dapall['Z'] * constants.c.to('km/s').value
    valid_galaxy = ((dapall['MNGTARG1'] > 0) | (dapall['MNGTARG3'] > 0)) & (cz > -500)
    ngal = numpy.sum(valid_galaxy)
    print('# Number of valid galaxies: {0}'.format(ngal))

    mid = dapall['MANGAID'][valid_galaxy].copy()
    plt = dapall['PLATE'][valid_galaxy].copy()
    ifu = dapall['IFUDESIGN'][valid_galaxy].copy()

    # Find the unique MaNGA IDs
    unid, unid_indx, unid_invs, unid_cnts = numpy.unique(mid, return_index=True,
                                                         return_inverse=True, return_counts=True)

    f = None if ofile is None else open(ofile, 'w')
    double_print(f, '# Number of DAPall entries: {0}'.format(mid.size))
    double_print(f, '# Number of unique MaNGA IDs: {0}'.format(unid.size))
    max_rep = numpy.amax(unid_cnts)
    double_print(f, '# Maximum number of repeated observations: {0}'.format(max_rep))
    double_print(f, '#')
    double_print(f, '#{0:>9s}'.format('MANGAID'), end='')
    for i in range(max_rep):
        double_print(f, ' {1:>4s}{0} {2:>4s}{0}'.format(i+1, 'PLT','IFU'), end='')
    double_print(f, '')
    for i in range(len(unid)):
        if unid_cnts[i] == 1:
            continue
        indx = numpy.where(mid == unid[i])[0]
        nempty = max_rep - len(indx)

        double_print(f, ' {0:>9s}'.format(unid[i]), end='')
        for j in range(len(indx)):
            double_print(f, ' {0:>5} {1:>5}'.format(plt[indx[j]], ifu[indx[j]]), end='')
        if nempty > 0:
            for j in range(nempty):
                double_print(f, ' {0:>5} {1:>5}'.format(-1, -1), end='')
        double_print(f, '')
    if f is not None:
        f.close()


class FindRepeatObservations(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):

        parser = super().get_parser(description='Use the DRPall file to find repeat observations',
                                    width=width)
        parser.add_argument('--drpver', type=str, help='DRP version', default=None)
        parser.add_argument('--dapver', type=str, help='DAP version', default=None)
        parser.add_argument('--analysis_path', type=str, help='main DAP output path', default=None)
        parser.add_argument('--dapall', type=str, help='full path to DAPall file', default=None)
        parser.add_argument('--output_file', type=str, default=None,
                            help='output file name, including path; no file is written by default')
        parser.add_argument('-o', '--overwrite', default=False, action='store_true',
                            help='Overwrite any existing files.')
        return parser

    @staticmethod
    def main(args):

        t = time.perf_counter()

        dapall_file = defaults.dapall_file(drpver=args.drpver, dapver=args.dapver,
                                           analysis_path=args.analysis_path) \
                        if args.dapall is None else args.dapall
        if not os.path.isfile(dapall_file):
            raise FileNotFoundError(f'File does not exist: {dapall_file}')

        if args.output_file is not None and os.path.isfile(args.output_file):
            if not overwrite:
                raise FileExistsError('Output file already exists, use the -o option to overwrite.')
            warnings.warn('Output file already exists and will be overwritten.')

        find_repeat_observations(dapall_file, ofile=args.output_file)

        print('Elapsed time: {0} seconds'.format(time.perf_counter() - t))



