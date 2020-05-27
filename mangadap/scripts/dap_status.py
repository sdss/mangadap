import os
import time
import datetime
import glob
import mmap
import numpy

import argparse

from mangadap.config import defaults
from mangadap.par.analysisplan import AnalysisPlanSet
from mangadap.proc.spatiallybinnedspectra import SpatiallyBinnedSpectra
from mangadap.proc.stellarcontinuummodel import StellarContinuumModel
from mangadap.proc.emissionlinemodel import EmissionLineModel


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
        return glob.glob(os.path.join(log_path, '*', '*', '*.ready')), \
                glob.glob(os.path.join(log_path, '*', '*', '*.started')), \
                glob.glob(os.path.join(log_path, '*', '*', '*.done')), \
                glob.glob(os.path.join(log_path, '*', '*', '*.err'))

    # Get and sort the list of directories
    lp = [d for d in os.listdir(os.path.join(analysis_path, 'log')) if 'UTC' in d]
    indx = numpy.argsort(numpy.array([datetime.datetime.strptime(t, '%d%b%YT%H.%M.%SUTC')
                                      for t in lp]))

    ndir = len(lp)
    ready = [None]*ndir
    start = [None]*ndir
    done = [None]*ndir
    err = [None]*ndir
    for i,j in enumerate(indx):
        ready[i], start[i], done[i], err[i] = build_lists(analysis_path, logdir=lp[j])

    pltifu = numpy.unique([ '-'.join(os.path.split(f)[1].split('.')[0].split('-')[1:3]) 
                             for f in numpy.concatenate(ready)])

    ncube = len(pltifu)

    last_ready = [None]*ncube
    last_start = [None]*ncube
    last_done = [None]*ncube
    last_err = [None]*ncube
    for j in range(ncube):
        for i in range(ndir-1,-1,-1):
            indx = numpy.array([pltifu[j] in e for e in done[i]])
            if not numpy.any(indx):
                continue
            k = numpy.where(indx)[0][0]
            last_ready[j] = ready[i][k]
            last_start[j] = start[i][k]
            last_done[j] = done[i][k]
            last_err[j] = err[i][k]
            break

    return last_ready, last_start, last_done, last_err


def dap_status(analysis_path, daptypes, logdir=None):

    ready_list, started_list, done_list, err_list = build_lists(analysis_path, logdir=logdir)
    log_path = os.path.join(analysis_path, 'log')
    if logdir is not None:
        log_path = os.path.join(log_path, logdir)

    print('Constructing plateifu lists.')
    ready_pltifu = [ '{0}-{1}'.format(*(r.split('/')[-3:-1])) for r in ready_list ]
    started_pltifu = [ '{0}-{1}'.format(*(r.split('/')[-3:-1])) for r in started_list ]
    done_pltifu = [ '{0}-{1}'.format(*(r.split('/')[-3:-1])) for r in done_list ]
    err_pltifu = ['{0}-{1}'.format(*(r.split('/')[-3:-1])) for r in err_list]

    print('Determine which executions finished in error.')
    errored = [ has_error(e) for e in err_list ]
    err_pltifu = (numpy.array(err_pltifu)[errored]).tolist() if numpy.any(errored) else []

    print('STATUS:')
    print('      Ready: {0}'.format(len(ready_pltifu)))
    print('    Started: {0}'.format(len(started_pltifu)))
    print('       DONE: {0}'.format(len(done_pltifu)))
    print('    Errored: {0}'.format(len(err_pltifu)))

    tracebacks, counts, trace_in_err = sort_errors(numpy.array(err_list)[errored])
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


def parse_args(options=None):
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('plan_file', type=str, help='parameter file with the MaNGA DAP '
                        'execution plan to use instead of the default')

    parser.add_argument('--logdir', type=str, help='Log path (e.g., 31Jan2019T19.16.28UTC)')

    parser.add_argument('--drpver', type=str, help='DRP version', default=None)
    parser.add_argument('--dapver', type=str, help='DAP version', default=None)
    parser.add_argument('--analysis_path', type=str, help='main DAP output path', default=None)
    return parser.parse_args() if options is None else parser.parse_args(options)


def main(args):
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




