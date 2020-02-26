import os
import time

from argparse import ArgumentParser

import numpy

from mangadap.util.log import init_DAP_logging, module_logging
from mangadap.par.analysisplan import AnalysisPlanSet
from mangadap.survey.dapall import DAPall

def parse_args(options=None):
    parser = ArgumentParser()

    parser.add_argument('plan', type=str, help='SDSS parameter file with analysis plan')

    parser.add_argument('--drpver', type=str, help='DRP version', default=None)
    parser.add_argument('-r', '--redux_path', type=str,
                        help='Top-level directory with the DRP products; defaults to '
                             '$MANGA_SPECTRO_REDUX/$MANGADRP_VER', default=None)
    parser.add_argument('--dapver', type=str, help='DAP version', default=None)
    parser.add_argument('-a', '--analysis_path', type=str,
                        help='Top-level output directory for the DAP results; defaults to '
                             '$MANGA_SPECTRO_ANALYSIS/$MANGADRP_VER/$MANGADAP_VER', default=None)
    parser.add_argument('-m', '--methods', type=str, nargs='+', default=None,
                        help='Only include output from this DAP method designation in the output')
    parser.add_argument('-v', '--verbose', action='count',
                        help='Set verbosity level; can be omitted and set up to -vv', default=0)
    parser.add_argument('--quiet', action='store_true', default=False,
                        help='suppress all terminal output')
    parser.add_argument('--double', dest='single_precision', action='store_false', default=True,
                        help='Output the floating-point data in double precision '
                             '(default is single precision)')
    return parser.parse_args() if options is None else parse.parse_args(options)

def main(args):
    t = time.perf_counter()
    analysisplan = AnalysisPlanSet.from_par_file(args.plan)

    # Initialize the logging objects and start the log
    init_DAP_logging(None)#, simple_warnings=False)
    loggers = module_logging(__name__, args.verbose)

    DAPall(analysisplan, methods=args.methods, drpver=args.drpver, redux_path=args.redux_path,
           dapver=args.dapver, analysis_path=args.analysis_path, loggers=loggers, quiet=args.quiet,
           single_precision=args.single_precision)

    print('Elapsed time: {0} seconds'.format(time.perf_counter() - t))



