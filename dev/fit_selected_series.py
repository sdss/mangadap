import os
import time
import warnings
import argparse

import numpy

try:
    import pbs.queue
except:
    warnings.warn('Could not import pbs.queue!  Any cluster submission will fail!', ImportWarning)

def parse_args():

    # Declare the ArgumentParser object
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--output_path', type=str, help='path for output', default=os.getcwd())
    parser.add_argument('--output_root', type=str, help='root name for outputfiles', default=None)

    parser.add_argument('--label', type=str, help='label for cluster job', default='daptest')
    parser.add_argument('--nodes', type=int, help='number of nodes to use in cluster', default=1)
    parser.add_argument('--cpus', type=int,
                        help='number of cpus to use per node.  Default is to use all available'
                             '; otherwise, set to minimum of provided number and number of '
                             'processors per node', default=None)
    parser.add_argument('--fast', dest='qos', type=str, help='qos state', default=None)
    parser.add_argument('--umask', type=str, help='umask bit for cluster job', default='0027')
    parser.add_argument('--walltime', type=str, help='walltime for cluster job',
                        default='240:00:00')
    parser.add_argument('--toughness', dest='hard', action='store_false', default=True,
                        help='turn off hard keyword for cluster submission')

    parser.add_argument('--create', action='store_true', default=False,
                        help='use the pbs package to create the cluster scripts')
    parser.add_argument('--submit', action='store_true', default=False,
                        help='submit the scripts to the cluster')
    
    return parser.parse_args()


#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
def main():

    # Read the arguments
    args = parse_args()
    if args.submit and not args.create:
        args.create = True

    # Hard-wired cases to run
    sc_tpl = ['MILESHC']*16
    sc_vsr = [2] + [4]*15
    sc_deg = [8]*16

    el_tpl = [None, 'MILESHC'] + ['MASTARHC2']*14
    el_vsr = [None, 4] + [1]*14
    el_deg = [None, 8] + (numpy.arange(14)+1).tolist()
    el_band = [None] + ['ELBMPL9']*15
    el_list = [None] + ['ELPMPL9']*15

    # Hard-wired location of script
    fit_script = os.path.abspath('fit_selected_spectra.py')
    # Hard-wired files to use
    benchmark = os.path.abspath('benchmark_spectra_v2.fits')
    flags = os.path.abspath('representative_spectra_flags_v2_fitall.db')

    # Output path for results
    output_path = os.path.abspath(args.output_path)
    if not os.path.isdir(output_path):
        os.makedirs(output_path)

    # Path for script and log files
    log_path = os.path.join(output_path, 'log')
    if not os.path.isdir(log_path):
        os.makedirs(log_path)

    # Build the queue if requested
    if args.create:
        queue = pbs.queue(verbose=True)
        ppn = 16
        cpus = ppn if args.cpus is None else min(args.cpus, ppn)
        queue.create(label=args.label, nodes=args.nodes, qos=args.qos, umask=args.umask,
                     walltime=args.walltime, ppn=ppn, cpus=cpus)

    nfits = len(sc_tpl)

    for i in range(nfits):
        ofile_root = '{0}_{1}_{2}_{3}'.format(args.output_root, sc_tpl[i], sc_vsr[i], sc_deg[i])
        if el_tpl[i] is not None:
            ofile_root += '_{0}_{1}_{2}_{3}_{4}'.format(el_tpl[i], el_vsr[i], el_deg[i],
                                                        el_band[i], el_list[i])

        # Set the names for the script, stdout, and stderr files
        scriptfile = os.path.join(log_path, ofile_root)
        stdoutfile = '{0}.out'.format(scriptfile)
        stderrfile = '{0}.err'.format(scriptfile)

        # Write the file
        file = open(scriptfile, 'w')
        file.write('# Auto-generated batch file\n')
        file.write('# {0}\n'.format(time.strftime("%a %d %b %Y %H:%M:%S",time.localtime())))
        file.write('\n')

        file.write('touch {0}\n'.format('{0}.started'.format(scriptfile)))
        file.write('\n')

        command = 'OMP_NUM_THREADS=1 python3 {0}'.format(fit_script)
        command += ' {0} {1}.fits.gz'.format(benchmark, ofile_root)
        command += ' --output_root {0}'.format(output_path)
        command += ' --spec_flags {0}'.format(flags)
        command += ' --sc_tpl {0}'.format(sc_tpl[i])
        command += ' --sc_vsr {0}'.format(sc_vsr[i])
        command += ' --sc_deg {0}'.format(sc_deg[i])
        if el_tpl[i] is None:
            command += ' --sc_only'
        else:
            command += ' --el_tpl {0}'.format(el_tpl[i])
            command += ' --el_vsr {0}'.format(el_vsr[i])
            command += ' --el_deg {0}'.format(el_deg[i])
            command += ' --el_band {0}'.format(el_band[i])
            command += ' --el_list {0}'.format(el_list[i])
        file.write('{0}\n'.format(command))
        file.write('\n')

        file.write('touch {0}\n'.format('{0}.done'.format(scriptfile)))
        file.write('\n')
        file.close()
        
        if args.create:
            queue.append('source {0}'.format(scriptfile), outfile=stdoutfile, errfile=stderrfile)

    # Submit the queue to the cluster
    if args.create:
        queue.commit(hard=args.hard, submit=args.submit)

if __name__ == '__main__':
    main()