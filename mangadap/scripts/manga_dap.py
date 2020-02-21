import os
import time

from argparse import ArgumentParser

from mangadap.par.obsinput import ObsInputPar
from mangadap.par.analysisplan import AnalysisPlanSet
from mangadap.survey.manga_dap import manga_dap

def parse_args(options=None):
    parser = ArgumentParser()

    parser.add_argument('obs', type=str, help='SDSS parameter file with observational input')
    parser.add_argument('plan', type=str, help='SDSS parameter file with analysis plan')
    parser.add_argument('--dbg', help='Run manga_dap in debug mode', action='store_true',
                        default=False)
    parser.add_argument('--log', type=str, help='File name for runtime log', default=None)
    parser.add_argument('-v', '--verbose', action='count',
                        help='Set verbosity level; can be omitted and set up to -vv', default=0)

    parser.add_argument('--drpver', type=str, help='DRP version', default=None)
    parser.add_argument('-r', '--redux_path', type=str,
                        help='Top-level directory with the DRP products; defaults to '
                             '$MANGA_SPECTRO_REDUX/$MANGADRP_VER', default=None)
    parser.add_argument('-d', '--directory_path', type=str,
                        help='Path directly to directory with DRP file to analyze', default=None)
    parser.add_argument('--dapver', type=str, help='DAP version', default=None)
    parser.add_argument('-a', '--analysis_path', type=str,
                        help='Top-level output directory for the DAP results; defaults to '
                             '$MANGA_SPECTRO_ANALYSIS/$MANGADRP_VER/$MANGADAP_VER', default=None)
    return parser.parse_args() if options is None else parser.parse_args(options)

def main(args):
    t = time.perf_counter()

    obspar = ObsInputPar.from_par_file(args.obs)
    analysisplan = AnalysisPlanSet.from_par_file(args.plan)

    status = manga_dap(obspar, analysisplan, dbg=args.dbg, log=args.log, verbose=args.verbose,
                       drpver=args.drpver, redux_path=args.redux_path,
                       directory_path=args.directory_path, dapver=args.dapver,
                       analysis_path=args.analysis_path)

    print('Elapsed time: {0} seconds'.format(time.perf_counter() - t))



