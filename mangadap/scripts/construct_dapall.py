
from mangadap.scripts import scriptbase

class ConstructDapAll(scriptbase.ScriptBase):

    @classmethod
    def name(cls):
        """
        Return the name of the executable.
        """
        return 'construct_dapall'

    @classmethod
    def get_parser(cls, width=None):

        parser = super().get_parser(description='Compile metadata for the DAPall file',
                                    width=width)

        parser.add_argument('--drpver', type=str, help='DRP version', default=None)
        parser.add_argument('-r', '--redux_path', type=str,
                            help='Top-level directory with the DRP products; defaults to '
                                '$MANGA_SPECTRO_REDUX/$MANGADRP_VER', default=None)
        parser.add_argument('--dapver', type=str, help='DAP version', default=None)
        parser.add_argument('-a', '--analysis_path', type=str, default=None,
                            help='Top-level output directory for the DAP results; defaults to '
                                '$MANGA_SPECTRO_ANALYSIS/$MANGADRP_VER/$MANGADAP_VER')
        parser.add_argument('--plan_file', type=str, help='toml file with the MaNGA DAP '
                            'execution plan to use instead of the default' , default=None)

        parser.add_argument('-m', '--methods', type=str, nargs='*', default=None,
                            help='Only include output for these DAP method designations in the '
                                 'output')
        parser.add_argument('-v', '--verbose', action='count', default=0,
                            help='Set verbosity level; can be omitted and set up to -vv')
        parser.add_argument('--quiet', action='store_true', default=False,
                            help='suppress all terminal output')
        parser.add_argument('--double', dest='single_precision', action='store_false', default=True,
                            help='Output the floating-point data in double precision '
                                '(default is single precision)')
        return parser

    @staticmethod
    def main(args):

        import time

        from mangadap.util.log import init_DAP_logging, module_logging
        from mangadap.config.manga import MaNGAAnalysisPlan
        from mangadap.survey.dapall import DAPall

        t = time.perf_counter()
        plan = manga.MaNGAAnalysisPlan.default() if args.plan_file is None \
                            else manga.MaNGAAnalysisPlan.from_toml(args.plan_file)

        # Initialize the logging objects and start the log
        init_DAP_logging(None)#, simple_warnings=False)
        loggers = module_logging(__name__, args.verbose)

        DAPall(plan, methods=args.methods, drpver=args.drpver, redux_path=args.redux_path,
               dapver=args.dapver, analysis_path=args.analysis_path, loggers=loggers,
               quiet=args.quiet, single_precision=args.single_precision)

        print('Elapsed time: {0} seconds'.format(time.perf_counter() - t))



