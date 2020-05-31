#!/usr/bin/env python3

"""
Dynamically build the rst documentation of the bitmasks.
"""

import os
import time

#-----------------------------------------------------------------------------

def write_help(script_mod, opath, prepend_dap=False):
    exe = script_mod.__name__.split('.')[-1]
    if prepend_dap:
        exe = 'dap_'+exe
    ofile = os.path.join(opath, '{0}.rst'.format(exe))
    lines = ['.. code-block:: console', '']
    lines += ['    $ {0} -h'.format(exe)]
    parser = script_mod.parse_args(return_parser=True)
    parser.prog = exe
    lines += ['    ' + l for l in parser.format_help().split('\n')]
    with open(ofile, 'w') as f:
        f.write('\n'.join(lines))

if __name__ == '__main__':
    t = time.perf_counter()

    import mangadap

    path = os.path.join(os.environ['MANGADAP_DIR'], 'docs', 'help')
    if not os.path.isdir(path):
        os.makedirs(path)

    from mangadap.scripts import calculate_covariance
    write_help(calculate_covariance, path, prepend_dap=True)

    from mangadap.scripts import construct_dapall
    write_help(construct_dapall, path)

    from mangadap.scripts import dap_status
    write_help(dap_status, path)

    from mangadap.scripts import dapall_qa
    write_help(dapall_qa, path)

    from mangadap.scripts import find_repeat_observations
    write_help(find_repeat_observations, path, prepend_dap=True)

    from mangadap.scripts import fit_residuals
    write_help(fit_residuals, path, prepend_dap=True)

    from mangadap.scripts import manga_dap_inspector
    write_help(manga_dap_inspector, path)

    from mangadap.scripts import manga_dap
    write_help(manga_dap, path)

    from mangadap.scripts import plate_fit_qa
    write_help(plate_fit_qa, path, prepend_dap=True)

    from mangadap.scripts import ppxffit_qa
    write_help(ppxffit_qa, path, prepend_dap=True)

    from mangadap.scripts import rundap
    write_help(rundap, path)

    from mangadap.scripts import spotcheck_dap_maps
    write_help(spotcheck_dap_maps, path)

    from mangadap.scripts import template_flux_norm
    write_help(template_flux_norm, path, prepend_dap=True)

    from mangadap.scripts import write_dap_config
    write_help(write_dap_config, path)

    print('Elapsed time: {0} seconds'.format(time.perf_counter() - t))


