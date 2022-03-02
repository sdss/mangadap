

from IPython import embed

from mangadap.scripts import scriptbase
from mangadap.util.pkg import load_object


class MangaDap(scriptbase.ScriptBase):

    @classmethod
    def name(cls):
        """
        Return the name of the executable.
        """
        return 'manga_dap'

    @classmethod
    def get_parser(cls, width=None):

        parser = super().get_parser(description='Perform analysis of integral-field data.',
                                    width=width)

        group = parser.add_mutually_exclusive_group(required=True)
        group.add_argument('-c', '--config', type=str, default=None,
                           help='Configuration file used to instantiate the relevant DataCube '
                                'derived class.')
        group.add_argument('-f', '--cubefile', type=str, default=None,
                           help='Name of the file with the datacube data.  Must be possible to '
                                'instantiate the relevant DataCube derived class directly from '
                                'the file only.')
        parser.add_argument('--cube_module', nargs='*',
                            default='mangadap.datacube.MaNGADataCube',
                            help='The name of the module that contains the DataCube derived class used to read the data.')
#        parser.add_argument('-o', '--cube_object', type=str, default='MaNGADataCube',
#                            help='The name of the DataCube derived class object.')

        parser.add_argument('--plan_module', nargs='*',
                            default='mangadap.config.MaNGAAnalysisPlan',
                            help='The name of the module used to define the analysis plan and the output paths.')

        parser.add_argument('-p', '--plan', type=str,
                            help='SDSS parameter file with analysis plan.  If not provided, a '
                                 'default plan is used.')

        parser.add_argument('--dbg', help='Run manga_dap in debug mode', action='store_true',
                            default=False)
        parser.add_argument('--log', type=str, help='File name for runtime log', default=None)
        parser.add_argument('-v', '--verbose', action='count', default=0,
                            help='Set verbosity level; can be omitted and set up to -vv')

        # TODO: Given the config file and the write_dap_config script,
        # should I continue to allow for these? They won't be meaningful
        # for other datacubes

        parser.add_argument('--drpver', type=str, help='DRP version', default=None)
        parser.add_argument('-r', '--redux_path', type=str,
                            help='Top-level directory with the DRP products; defaults to '
                                '$MANGA_SPECTRO_REDUX/$MANGADRP_VER', default=None)
        parser.add_argument('-d', '--directory_path', type=str, default=None,
                            help='Path directly to directory with DRP file to analyze')
        parser.add_argument('--dapver', type=str, help='DAP version', default=None)
        parser.add_argument('-a', '--analysis_path', type=str, default=None,
                            help='Top-level output directory for the DAP results; defaults to '
                                '$MANGA_SPECTRO_ANALYSIS/$MANGADRP_VER/$MANGADAP_VER')
        return parser

    @staticmethod
    def main(args):

#        import warnings
#        warnings.simplefilter('error', DeprecationWarning)
        from mangadap.par.analysisplan import AnalysisPlanSet
        from mangadap.survey.manga_dap import manga_dap
        from mangadap.datacube import DataCube

        # Instantiate the DataCube
        #   - Try to import the module used to read the datacube
        UserDataCube = load_object(args.cube_module, args.cube_object)
        #   - Check that the class is derived from DataCube
        if not issubclass(UserDataCube, DataCube):
            raise TypeError('Defined cube object must subclass from mangadap.datacube.DataCube.')
        #   - Instantiate using either the datacube file directly or a
        #     configuration file
        cube = UserDataCube(args.cubefile) if args.config is None \
                    else UserDataCube.from_config(args.config)

        # Read the analysis plan
        analysisplan = AnalysisPlanSet.default() if args.plan is None \
                            else AnalysisPlanSet.from_par_file(args.plan)

        # Run the pipeline
        status = manga_dap(cube, analysisplan, dbg=args.dbg, log=args.log, verbose=args.verbose,
                           drpver=args.drpver, redux_path=args.redux_path,
                           directory_path=args.directory_path, dapver=args.dapver,
                           analysis_path=args.analysis_path)

