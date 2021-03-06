#!/usr/bin/env python3

import numpy
import os
import warnings
import time

import argparse

from astropy.io import fits

from mangadap.util.bitmask import BitMask
from mangadap.survey.drpcomplete import DRPComplete
from mangadap.config import defaults


def double_print(ostream, line, **kwargs):
    print(line, **kwargs)
    print(line, file=ostream, **kwargs)


def find_repeat_observations(output_file, drpver, redux_path, dapver, analysis_path,
                             directory_path):

    # Get the DRPComplete database
    drpc = DRPComplete(drpver=drpver, redux_path=redux_path, dapver=dapver,
                       analysis_path=analysis_path, directory_path=directory_path, readonly=True)
    try:
        drpc._confirm_access()
    except FileNotFoundError:
        warnings.warn('DRPComplete file must be present to find repeat observations!')
        raise
    except:
        raise

    # Get the maskbits file
    sdssMaskbits = defaults.sdss_maskbits_file()

    # Instantiate the BitMask object
    mngtarg1_bm = BitMask.from_par_file(sdssMaskbits, 'MANGA_TARGET1')
    mngtarg3_bm = BitMask.from_par_file(sdssMaskbits, 'MANGA_TARGET3')

    print('# DRPComplete file: {0}'.format(drpc.file_path()))
    print('# Total number of DRPComplete entries: {0}'.format(drpc.nobs))
    valid_galaxy = ((drpc['MANGA_TARGET1'] > 0) | (drpc['MANGA_TARGET3'] > 0)) & (drpc['VEL'] > 0)
    ngal = numpy.sum(valid_galaxy)
    print('# Number of valid galaxies: {0}'.format(ngal))

    mid = drpc['MANGAID'][valid_galaxy].copy()
    plt = drpc['PLATE'][valid_galaxy].copy()
    ifu = drpc['IFUDESIGN'][valid_galaxy].copy()

    # Find the unique MaNGA IDs
    unid, unid_indx, unid_invs, unid_cnts = numpy.unique(mid, return_index=True,
                                                         return_inverse=True, return_counts=True)

    with open(output_file, 'w') as f:
        double_print(f, '# Number of DRPComplete entries: {0}'.format(mid.size))
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

def parse_args(options=None, return_parser=False):
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--drpver', type=str, help='DRP version', default=None)
    parser.add_argument('--redux_path', type=str, help='main DRP output path', default=None)
    parser.add_argument('--dapver', type=str, help='DAP version', default=None)
    parser.add_argument('--analysis_path', type=str, help='main DAP output path', default=None)
    parser.add_argument('--directory_path', type=str, help='directory with DRP complete file',
                        default=None)
    parser.add_argument('--output_file', type=str, help='output file name, including path',
                        default=None)
    if return_parser:
        return parser
    return parser.parse_args() if options is None else parser.parse_args(options)

def main(args):

    t = time.perf_counter()

    drpver = defaults.drp_version() if args.drpver is None else args.drpver
    dapver = defaults.dap_version() if args.dapver is None else args.dapver

    if args.output_file is None:
        path_root = defaults.dap_common_path(drpver=drpver, dapver=dapver,
                                             analysis_path=args.analysis_path)
        ofile = 'repeat-observations-{0}-{1}.db'.format(drpver, dapver)
        output_file = os.path.join(path_root, ofile)
    else:
        output_file = args.output_file

    find_repeat_observations(output_file, drpver, args.redux_path, dapver, args.analysis_path,
                             args.directory_path)

    print('Elapsed time: {0} seconds'.format(time.perf_counter() - t))



