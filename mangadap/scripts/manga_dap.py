

from IPython import embed

from mangadap.scripts import scriptbase


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
                            help='The name of the module that contains the DataCube derived '
                                 'class used to read the data.')

        parser.add_argument('-p', '--plan', type=str,
                            help='TOML file with analysis plan.  If not provided, a '
                                 'default plan is used.')
        parser.add_argument('--plan_module', nargs='*',
                            default='mangadap.config.manga.MaNGAAnalysisPlan',
                            help='The name of the module used to define the analysis plan and '
                                 'the output paths.')

        parser.add_argument('--dbg', help='Run manga_dap in debug mode', action='store_true',
                            default=False)
        parser.add_argument('--log', type=str, help='File name for runtime log', default=None)
        parser.add_argument('-v', '--verbose', action='count', default=0,
                            help='Set verbosity level; can be omitted and set up to -vv')

        parser.add_argument('-o', '--output_path', type=str, default=None,
                            help='Top-level directory for the DAP output files; default path is '
                                 'set by the provided analysis plan object (see plan_module).')
        return parser

    @staticmethod
    def main(args):

#        import warnings
#        warnings.simplefilter('error', DeprecationWarning)
        from mangadap.config.analysisplan import AnalysisPlan
        from mangadap.survey.manga_dap import manga_dap
        from mangadap.datacube import DataCube
        from mangadap.util.pkg import load_object

        # Execute the MaNGA DAP for a single datacube

        #   - Import the module used to read the datacube
        if isinstance(args.cube_module, list) and len(args.cube_module) > 2:
            raise ValueError('Provided cube module must be one or two strings.')
        if isinstance(args.cube_module, str) or len(args.cube_module) == 1:
            UserDataCube = load_object(args.cube_module if isinstance(args.cube_module, str) 
                                       else args.cube_module[0])
        else:
            UserDataCube = load_object(args.cube_module[0], obj=args.cube_module[1])
        #   - Check that the class is derived from DataCube
        if not issubclass(UserDataCube, DataCube):
            raise TypeError('Defined cube object must subclass from mangadap.datacube.DataCube.')

        #   - Import the module used to set the analysis plan
        if isinstance(args.plan_module, list) and len(args.plan_module) > 2:
            raise ValueError('Provided plan module must be one or two strings.')
        if isinstance(args.plan_module, str) or len(args.plan_module) == 1:
            UserPlan = load_object(args.plan_module if isinstance(args.plan_module, str) 
                                       else args.plan_module[0])
        else:
            UserPlan = load_object(args.plan_module[0], obj=args.plan_module[1])
        #   - Check that the class is derived from AnalysisPlan
        if not issubclass(UserPlan, AnalysisPlan):
            raise TypeError('Defined plan object must subclass from '
                            'mangadap.config.analysisplan.AnalysisPlan')

        #   - Instantiate the datacube object using either the datacube file
        #     directly or a configuration file
        cube = UserDataCube(args.cubefile) if args.config is None \
                    else UserDataCube.from_config(args.config)

        #   - Read the analysis plan
        plan = UserPlan.default(cube=cube, analysis_path=args.output_path) if args.plan is None \
                    else UserPlan.from_toml(args.plan, cube=cube, analysis_path=args.output_path)

        #   - Run the pipeline
        status = manga_dap(cube, plan, dbg=args.dbg, log=args.log, verbose=args.verbose)


