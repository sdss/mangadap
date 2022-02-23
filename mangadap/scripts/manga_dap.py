

from IPython import embed

from mangadap.scripts import scriptbase

def load_cube_object(cube_module, cube_object):
    """
    Load the module and object used to read the datacube to analyze.

    Thanks to: https://stackoverflow.com/questions/67631/how-to-import-a-module-given-the-full-path?rq=1

    Args:
        cube_module (:obj:`str`):
            The name of a global python module, or the root name of a local file
            with the :class:`~mangadap.datacube.datacube.DataCube` subclass.
        cube_object (:obj:`str`):
            The name of the subclass

    Return:
        :obj:`type`: The object used to read the data to analyze.

    Raises:
        ImportError:
            Raised if unable to import ``cube_module``.
    """
    import importlib
    from pathlib import Path

    try:
        CubeModule = importlib.import_module(cube_module)
    except (ModuleNotFoundError, ImportError) as e:
        p = Path(cube_module + '.py').resolve()
        if not p.exists():
            raise ImportError(f'Unable to load module {cube_module}!') from e
        spec = importlib.util.spec_from_file_location(cube_module, str(p))
        CubeModule = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(CubeModule)

    return getattr(CubeModule, cube_object)


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
                                'instantiate the relevant DataCube derived class directly from the '
                                'file only.')
        parser.add_argument('-p', '--plan', type=str,
                            help='SDSS parameter file with analysis plan.  If not provided, a '
                                'default plan is used.')
        parser.add_argument('-m', '--cube_module', type=str, default='mangadap.datacube',
                            help='The name of the module that contains the DataCube derived class.')
        parser.add_argument('-o', '--cube_object', type=str, default='MaNGADataCube',
                            help='The name of the DataCube derived class object.')

        parser.add_argument('--dbg', help='Run manga_dap in debug mode', action='store_true',
                            default=False)
        parser.add_argument('--log', type=str, help='File name for runtime log', default=None)
        parser.add_argument('-v', '--verbose', action='count',
                            help='Set verbosity level; can be omitted and set up to -vv', default=0)

        # TODO: Given the config file and the write_dap_config script,
        # should I continue to allow for these? They won't be meaningful
        # for other datacubes

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
        return parser

    @staticmethod
    def main(args):

        from mangadap.par.analysisplan import AnalysisPlanSet
        from mangadap.survey.manga_dap import manga_dap
        from mangadap.datacube import DataCube

        # Instantiate the DataCube
        #   - Try to import the module
        UserDataCube = load_cube_object(args.cube_module, args.cube_object)
        #   - Check that the class is derived from DataCube
        if not issubclass(UserDataCube, DataCube):
            raise TypeError('Defined cube object must subclass from mangadap.datacube.DataCube.')
        #   - Instantiate using either the datacube file directly or a
        #     configuration file
        cube = UserDataCube(args.cubefile) if args.config is None \
                    else UserDataCube.from_config(args.config, drpver=args.drpver,
                                                redux_path=args.redux_path,
                                                directory_path=args.directory_path)

        # Read the analysis plan
        analysisplan = AnalysisPlanSet.default() if args.plan is None \
                            else AnalysisPlanSet.from_par_file(args.plan)

        # Run the pipeline
        status = manga_dap(cube, analysisplan, dbg=args.dbg, log=args.log, verbose=args.verbose,
                        drpver=args.drpver, redux_path=args.redux_path,
                        directory_path=args.directory_path, dapver=args.dapver,
                        analysis_path=args.analysis_path)

