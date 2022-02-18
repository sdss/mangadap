import os
import datetime
import glob
import mmap

from IPython import embed

import numpy

from mangadap.scripts import scriptbase


def has_error(log_file):
    with open(log_file, 'rb', 0) as f, mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ) as s:
        if s.find(b'Traceback') != -1:
            return True
    return False


def sort_errors(err_list):
    tracebacks = []
    for e in err_list:
        with open(e, 'r') as f:
            lines = f.readlines()
        for i,l in enumerate(lines):
            if 'Traceback' in l:
                break
        if i == len(lines):
            continue
        tracebacks += [(lines[0] if i == 1 
                            else (lines[i-2] if 'ERROR' in lines[i-2] else lines[i-1])).strip()]

    unique_traces, inv, counts = numpy.unique(tracebacks, return_inverse=True, return_counts=True)
    return unique_traces, counts, numpy.arange(len(err_list))[inv]
    

def error_report(tracebacks, counts, trace_in_err, err_pltifu, analysis_path, daptypes, ofile):
    with open(ofile, 'w') as f:
        # Print the Tracebacks
        f.write('# Unique Tracebacks\n#\n')
        f.write('# {0:>4} {1:>4} {2}\n'.format('INDX', 'CNT', 'Traceback'))
        for i in range(len(tracebacks)):
            f.write('  {0:>4} {1:>4} {2}\n'.format(i, counts[i], tracebacks[i]))
        f.write('\n\n')

        # Associate each PLATEIFU that failed with its Traceback
        f.write('# Traceback for each PLATEIFU ending in error\n')
        _err_pltifu = numpy.asarray(err_pltifu)
        for i in range(len(tracebacks)):
            f.write('# Traceback {0}: {1}\n'.format(i, tracebacks[i]))
            indx = numpy.asarray(trace_in_err) == i
            for pltifu in _err_pltifu[indx]:
                f.write('{0:>5} {1:>5}\n'.format(*pltifu.split('-')))
            f.write('\n')
        f.write('\n\n')

        # Find if any of the MAPS files were produced for the errored
        # plateifus
        f.write('# MAPS written for errored PLATEIFUs\n')
        s = '# {0:>5} {1:>5}'.format('PLATE', 'IFU')
        for d in daptypes:
            s += ' {0:>30}'.format(d)
        f.write(s + '\n')
        for i in range(len(err_pltifu)):
            plt, ifu = err_pltifu[i].split('-')
            s = '  {0:>5} {1:>5}'.format(plt, ifu)
            for d in daptypes:
                s += ' {0:>30}'.format(os.path.isfile(os.path.join(analysis_path, d, plt, ifu,
                                        'manga-{0}-{1}-MAPS-{2}.fits.gz'.format(plt, ifu, d))))
            f.write(s + '\n')
        f.write('\n\n')


def build_lists(analysis_path, logdir=None):

    if logdir is not None:
        log_path = os.path.join(analysis_path, 'log', logdir)
        print('Searching for status files at:\n    {0}'.format(log_path))

        # Need to ensure the returned lists are sorted and coincident
        ready = glob.glob(os.path.join(log_path, '*', '*', '*.ready'))

        root = ['.'.join(r.split('.')[:-1]) for r in ready]
        started = [r+'.started' if os.path.isfile(r+'.started') else None for r in root]
        done = [r+'.done' if os.path.isfile(r+'.done') else None for r in root]
        err = [r+'.err' if os.path.isfile(r+'.err') else None for r in root]

        return numpy.array(ready, dtype=object), numpy.array(started, dtype=object), \
                numpy.array(done, dtype=object), numpy.array(err, dtype=object)

    # Get and sort the list of directories
    lp = [d for d in os.listdir(os.path.join(analysis_path, 'log')) if 'UTC' in d]
    indx = numpy.argsort(numpy.array([datetime.datetime.strptime(t, '%d%b%YT%H.%M.%SUTC')
                                      for t in lp]))

    # Find all the relevant files
    ndir = len(lp)
    ready = [None]*ndir
    start = [None]*ndir
    done = [None]*ndir
    err = [None]*ndir
    for i,j in enumerate(indx):
        ready[i], start[i], done[i], err[i] = build_lists(analysis_path, logdir=lp[j])

    # Get the list of all possible pltifus
    pltifu = numpy.unique([ '-'.join(os.path.split(f)[1].split('.')[0].split('-')[1:3]) 
                             for f in numpy.concatenate(ready)])
    ncube = len(pltifu)

    # Initialize the lists of the most recent status files
    last_ready = numpy.empty(ncube, dtype=object)
    last_start = numpy.empty(ncube, dtype=object)
    last_done = numpy.empty(ncube, dtype=object)
    last_err = numpy.empty(ncube, dtype=object)

    # Find the status files for the executions that most recently
    # finished
    print('Finding most recently finished executions.')
    for j in range(ncube):
        for i in range(ndir-1,-1,-1):
            indx = numpy.array([e is not None and pltifu[j] in e for e in done[i]])
            if not numpy.any(indx):
                continue
            if numpy.sum(indx) > 1:
                raise ValueError('Degenerate plate ifu')
            k = numpy.where(indx)[0][0]
            last_ready[j] = ready[i][k]
            last_start[j] = start[i][k]
            last_done[j] = done[i][k]
            last_err[j] = err[i][k]
            break

    # For those that did not finish, find the most recent status files
    # for the executions that were only started
    not_finished = numpy.where(last_ready == None)[0]
    print('Finding most recently started executions that did not finish.')
    for j in not_finished:
        for i in range(ndir-1,-1,-1):
            indx = numpy.array([e is not None and pltifu[j] in e for e in start[i]])
            if not numpy.any(indx):
                continue
            if numpy.sum(indx) > 1:
                raise ValueError('Degenerate plate ifu')
            k = numpy.where(indx)[0][0]
            last_ready[j] = ready[i][k]
            last_start[j] = start[i][k]
            last_err[j] = err[i][k]
            break
        
    # For those that were not started, find the most recent status files
    # for the executions that were only ready
    not_started = numpy.where(last_ready == None)[0]
    print('Finding most recent ready executions that did not start.')
    for j in not_started:
        for i in range(ndir-1,-1,-1):
            indx = numpy.array([e is not None and pltifu[j] in e for e in ready[i]])
            if not numpy.any(indx):
                continue
            if numpy.sum(indx) > 1:
                raise ValueError('Degenerate plate ifu')
            k = numpy.where(indx)[0][0]
            last_ready[j] = ready[i][k]
            break

    return last_ready, last_start, last_done, last_err


def dap_status(analysis_path, daptypes, logdir=None):

    ready, started, done, err = build_lists(analysis_path, logdir=logdir)

    # There should be no remaining None values in last_ready
    if numpy.any(ready == None):
        raise ValueError('Should not be any None entries in list of ready files.')

    # Remove any Nones
    done = done[numpy.logical_not(done == None)]
    err = err[numpy.logical_not(err == None)]
    started = started[numpy.logical_not(started == None)]

    log_path = os.path.join(analysis_path, 'log')
    if logdir is not None:
        log_path = os.path.join(log_path, logdir)

    print('Constructing plateifu lists.')
    ready_pltifu = [ '{0}-{1}'.format(*(r.split('/')[-3:-1])) for r in ready ]
    started_pltifu = [ '{0}-{1}'.format(*(r.split('/')[-3:-1])) for r in started ]
    done_pltifu = [ '{0}-{1}'.format(*(r.split('/')[-3:-1])) for r in done ]
    err_pltifu = ['{0}-{1}'.format(*(r.split('/')[-3:-1])) for r in err]

    print('Determine which executions finished in error.')
    errored = [ has_error(e) for e in err ]
    err_pltifu = (numpy.array(err_pltifu)[errored]).tolist() if numpy.any(errored) else []

    print('STATUS:')
    print('      Ready: {0}'.format(len(ready_pltifu)))
    print('    Started: {0}'.format(len(started_pltifu)))
    print('       DONE: {0}'.format(len(done_pltifu)))
    print('    Errored: {0}'.format(len(err_pltifu)))

    tracebacks, counts, trace_in_err = sort_errors(numpy.array(err)[errored])
    print('Errors:')
    for i in range(len(tracebacks)):
        print('{0:>4} {1:>4} {2}'.format(i, counts[i], tracebacks[i]))

    # Print the error report
    ofile = os.path.join(log_path, 'error_report.txt')
    print('Error report written to:\n    {0}'.format(ofile))
    error_report(tracebacks, counts, trace_in_err, err_pltifu, analysis_path, daptypes, ofile)

    print('Comparing lists')
    ready_but_not_started = list(set(ready_pltifu)-set(started_pltifu))
    started_but_not_done = list(set(started_pltifu)-set(done_pltifu))

    print('Ready but not started: {0}'.format(len(ready_but_not_started)))
    print('Started but not DONE: {0}'.format(len(started_but_not_done)))

    if len(ready_but_not_started) > 0:
        ofile = os.path.join(log_path, 'not_started.lst')
        print('PLATEIFUs not started written to:\n    {0}'.format(ofile))
        with open(ofile, 'w') as f:
            for rns in ready_but_not_started:
                f.write('{0:>5} {1:>5}\n'.format(*rns.split('-')))

    if len(started_but_not_done) > 0:
        ofile = os.path.join(log_path, 'not_finished.lst')
        print('PLATEIFUs not completed written to:\n    {0}'.format(ofile))
        with open(ofile, 'w') as f:
            for snd in started_but_not_done:
                f.write('{0:>5} {1:>5}\n'.format(*snd.split('-')))

    if len(err_pltifu) > 0:
        ofile = os.path.join(log_path, 'errored.lst')
        print('PLATEIFUs finished in error written to:\n    {0}'.format(ofile))
        with open(ofile, 'w') as f:
            for err in err_pltifu:
                f.write('{0:>5} {1:>5}\n'.format(*err.split('-')))


class DapStatus(scriptbase.ScriptBase):

    @classmethod
    def name(cls):
        """
        Return the name of the executable.
        """
        return 'dap_status'

    @classmethod
    def get_parser(cls, width=None):

        parser = super().get_parser(description='Check status and log files for DAP progress',
                                    width=width)

        parser.add_argument('plan_file', type=str, help='parameter file with the MaNGA DAP '
                            'execution plan to use instead of the default')

        parser.add_argument('--logdir', type=str, help='Log path (e.g., 31Jan2019T19.16.28UTC)')

        parser.add_argument('--drpver', type=str, help='DRP version', default=None)
        parser.add_argument('--dapver', type=str, help='DAP version', default=None)
        parser.add_argument('--analysis_path', type=str, help='main DAP output path', default=None)

        return parser

    @staticmethod
    def main(args):

        import time
        from mangadap.config import defaults
        from mangadap.par.analysisplan import AnalysisPlanSet
        from mangadap.proc.spatiallybinnedspectra import SpatiallyBinnedSpectra
        from mangadap.proc.stellarcontinuummodel import StellarContinuumModel
        from mangadap.proc.emissionlinemodel import EmissionLineModel

        t = time.perf_counter()

        # Set the the analysis path and make sure it exists
        analysis_path = defaults.dap_analysis_path(drpver=args.drpver, dapver=args.dapver) \
                                if args.analysis_path is None else args.analysis_path

        analysisplan = AnalysisPlanSet.from_par_file(args.plan_file)
        daptypes = []
        for p in analysisplan:
            bin_method = SpatiallyBinnedSpectra.define_method(p['bin_key'])
            sc_method = StellarContinuumModel.define_method(p['continuum_key'])
            el_method = EmissionLineModel.define_method(p['elfit_key'])
            daptypes += [defaults.dap_method(bin_method['key'],
                                            sc_method['fitpar']['template_library_key'],
                                            el_method['continuum_tpl_key'])]

        dap_status(analysis_path, daptypes, logdir=args.logdir)

        print('Elapsed time: {0} seconds'.format(time.perf_counter() - t))

