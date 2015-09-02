#!/usr/bin/env python3

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

import sys
if sys.version > '3':
    long = int

#-----------------------------------------------------------------------
from os import remove, environ
from time import clock
from mangadap.util.defaults import default_drp_version, default_dap_version
from mangadap.util.defaults import default_analysis_path, default_dap_directory_path
#-----------------------------------------------------------------------

def build_map_files(plt, ifu, mpl, output_path):
    dap_source = default_dap_source()
    file_root = default_dap_file_root(plate, ifudesign, mode)
        
    scr = os.path.join(dap_source, 'bin', 'dap_data_to_fits_dube.py')

    bintype = [ 'NONE', 'STON' ]

    for bt in bintype:
        file_search_str = os.path.join(output_path, file_root+'_BIN-'+bt+'*fits')
        file_list = glob.glob(file_search_str)
        for f in file_list:
            bintype_index = f.find(bt)
            ofile = '{0}_MAP-{1}'.format(f[:f.find('BIN')],f[bintype_index:])
            indx = int(f[bintype_indx+4:bintype_indx+7])
            # TODO: Make this a system call
            print('{0} {1} {2} {3} {4} {5} {6} {7}\n\n'.format(scr, plate, ifudesign, bt, mpl,
                                                               ofile, output_path, indx))


#-----------------------------------------------------------------------
if __name__ == '__main__':
    t = clock()

    if len(sys.argv) != 5:
        print('Usage: build_map_files.py <plate> <ifudesign> <MPL-3 or MPL-4> <directory path>')
        exit(1)

    plt = str(sys.argv[1])
    ifu = str(sys.argv[2])
    mpl = str(sys.argv[3])
    if len(sys.argv) == 5:
        output_path = str(sys.argv[4])
    else :
        drpver = default_drp_version()
        dapver = default_dap_version()
        analysis_path = default_analysis_path(drpver, dapver)
        output_path = default_dap_directory_path(drpver, dapver, analysis_path, plt, ifu)

    build_map_files(plt, ifu, mpl, output_path)

    print('Elapsed time: {0} seconds'.format(clock() - t))


