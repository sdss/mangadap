
from IPython import embed

from mangadap.scripts import scriptbase


class RunDap(scriptbase.ScriptBase):
    """
    Simple wrapper class for the rundap script.
    """
    @classmethod
    def name(cls):
        """
        Return the name of the executable.
        """
        return 'rundap'

    @classmethod
    def get_parser(cls, width=None):

        parser = super().get_parser(description='Perform analysis of integral-field data.',
                                    width=width)

        # Read the optional run-mode arguments
        parser.add_argument("--overwrite",
                            help="if all selected, will run dap for all plates/ifudesigns/modes "
                            " regardless of state", action="store_true", default=False)
        parser.add_argument('-v', '--verbose', action='count',
                            help='Set verbosity level for manga_dap; can be omitted and set up '
                                 'to -vv', default=0)
        parser.add_argument("--quiet", help="suppress screen output", action="store_true",
                            default=False)
        parser.add_argument("--print_version", help="print DAP version and stop",
                            action="store_true", default=False)
       
        # These arguments are used to override default behavior
        parser.add_argument("--drpver", type=str, default=None,
                            help='MaNGA DRP version for analysis; $MANGADRP_VER by default')
        parser.add_argument("--redux_path", type=str, help="main DRP output path", default=None)
        parser.add_argument("--dapver", type=str, default=None,
                            help='optional output version, different from product version.  This '
                                 '*only* affects the output directory structure.  It does *not* '
                                 'select the version of the DAP to use.')
        parser.add_argument("--analysis_path", type=str, help="main DAP output path", default=None)

        parser.add_argument("--plan_file", type=str, help="parameter file with the MaNGA DAP "
                            "execution plan to use instead of the default" , default=None)
        parser.add_argument("--platelist", type=str, help="set list of plates to reduce",
                            default=None)
        parser.add_argument("--ifudesignlist", type=str, help="set list of ifus to reduce",
                            default=None)
        parser.add_argument("--list_file", type=str,
                            help='A file with the list of plates and ifudesigns to analyze',
                            default=None)
        parser.add_argument("--combinatorics", help="force execution of all permutations of the "
                            "provided lists", action="store_true", default=False)

        parser.add_argument('--sres_ext', type=str, default=None,
                            help='Spectral resolution extension to use.  Default set by '
                                 'MaNGADataCube class.')
        parser.add_argument('--sres_fill', type=str, default=None,
                            help='If present, use interpolation to fill any masked pixels in the '
                                 'spectral resolution vectors. Default set by MaNGADataCube '
                                 'class.')
        parser.add_argument('--covar_ext', type=str, default=None,
                            help='Use this extension to define the spatial correlation matrix.  '
                                 'Default set by MaNGADataCube class.')

        parser.add_argument('--on_disk', action='store_true', default=False,
                            help='When using the DRPall file to collate the data for input to '
                                 'the DAP, search for available DRP files on disk instead of '
                                 'using the DRPall file content.')
        parser.add_argument('--can_analyze', action='store_true', default=False,
                            help='Only construct script files for datacubes that can/should be '
                                 'analyzed by the DAP.  See '
                                 ':func:`~mangadap.survey.drpcomplete.DRPComplete.can_analyze`.')

        parser.add_argument("--log", help="Have the main DAP executable produce a log file",
                            action="store_true", default=False)
        parser.add_argument("--no_proc", help="Do NOT perform the main DAP processing steps",
                            action="store_true", default=False)
        parser.add_argument("--no_plots", help="Do NOT create QA plots", action="store_true",
                            default=False)
        parser.add_argument("--post", help="Create/Submit the post-processing scripts",
                            action="store_true", default=False)
        parser.add_argument("--post_plots", action="store_true", default=False,
                            help="Create/Submit the post-processing plotting scripts")

        # Read arguments specific to the cluster submission behavior
        parser.add_argument("--label", type=str, help='label for cluster job', default='mangadap')
        parser.add_argument("--nodes", type=int, help='number of nodes to use in cluster',
                            default=1)
        parser.add_argument("--cpus", type=int,
                            help='number of cpus to use per node.  Default is to use all available'
                                 '; otherwise, set to minimum of provided number and number of '
                                 'processors per node', default=None)
        parser.add_argument("--fast", dest='qos', type=str, help='qos state', default=None)
        parser.add_argument("--umask", type=str, help='umask bit for cluster job', default='0027')
        parser.add_argument("--walltime", type=str, help='walltime for cluster job',
                            default='240:00:00')
        parser.add_argument("--toughness", dest='hard', action='store_false', default=True,
                            help='turn off hard keyword for cluster submission')

        parser.add_argument('--create', action='store_true', default=False,
                            help='use the pbs package to create the cluster scripts')
        parser.add_argument('--submit', action='store_true', default=False,
                            help='submit the scripts to the cluster')
        parser.add_argument('--progress', action='store_true', default=False,
                            help='instead of closing the script, report the progress of the '
                                 'analysis on the cluster; this is required if you want to submit '
                                 'the DAPall script immediately after completing the individual '
                                 'cube analysis')

        parser.add_argument("--queue", dest='queue', type=str, help='set the destination queue',
                            default=None)
        return parser

    @staticmethod
    def main(args):

        if args.print_version:
            from mangadap import __version__
            print(f'DAP Version: {__version__}')
            return

        from mangadap.survey.rundap import rundap

        # Handle special argument cases
        if args.qos is not None and args.nodes > 1:
            warnings.warn('Requesting the fast node requires node=1.  Ignoring input node number.')
            nodes = 1               # Force the number of nodes to be 1
        else:
            nodes = args.nodes

        _rundap = rundap(overwrite=args.overwrite, quiet=args.quiet, drpver=args.drpver,
                         redux_path=args.redux_path, dapver=args.dapver,
                         analysis_path=args.analysis_path, plan_file=args.plan_file,
                         platelist=args.platelist, ifudesignlist=args.ifudesignlist,
                         combinatorics=args.combinatorics, list_file=args.list_file,
                         sres_ext=args.sres_ext, sres_fill=args.sres_fill,
                         covar_ext=args.covar_ext, on_disk=args.on_disk,
                         can_analyze=args.can_analyze, log=args.log, dapproc=not args.no_proc,
                         pltifu_plots=not args.no_plots, post_process=args.post,
                         post_plots=args.post_plots, report_progress=args.progress,
                         verbose=args.verbose, label=args.label, nodes=nodes, cpus=args.cpus,
                         qos=args.qos, umask=args.umask, walltime=args.walltime, hard=args.hard,
                         create=args.create, submit=args.submit, queue=args.queue)



