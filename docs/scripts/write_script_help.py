#!/usr/bin/env python3

"""
Dynamically build the rst documentation of the bitmasks.
"""

import time
from pathlib import Path

from mangadap.config.defaults import dap_source_dir
from mangadap.scripts import script_classes

#-----------------------------------------------------------------------------

# TODO: rundap isn't updated because it doesn't yet inherit from ScriptBase

def write_help(script_cls, opath, width=80):
    exe = script_cls.name()
    ofile = opath / f'{exe}.rst'
    lines = ['.. code-block:: console', '']
    lines += ['    $ {0} -h'.format(exe)]
    parser = script_cls.get_parser(width=80)
    parser.prog = exe
    lines += ['    ' + l for l in parser.format_help().split('\n')]
    print(f'Writing: {ofile}')
    with open(str(ofile), 'w') as f:
        f.write('\n'.join(lines))

if __name__ == '__main__':
    t = time.perf_counter()

    path = dap_source_dir() / 'docs' / 'help'
    if not path.is_dir():
        path.mkdir()

    # Get the list of script names and script classes
    scr_clss = script_classes()

    for name, script_cls in scr_clss.items():
        write_help(script_cls, path)

    print('Elapsed time: {0} seconds'.format(time.perf_counter() - t))

