#!/usr/bin/env python3

import numpy
from os import remove
import os.path
from time import clock

from mangadap.drpfile import drpfile
from mangadap.TplLibrary import TplLibrary
from matplotlib import pyplot

#-----------------------------------------------------------------------------

def run_test():

    drpf = drpfile(7495, 12703, 'CUBE')

    print(drpf.file_path())

    tpl_lib = TplLibrary('M11-MILES', drpf=drpf, directory_path='.', force=True)

#    tpl_lib.hdu.info()
#    print(tpl_lib.hdu['WAVE'].shape)
#    print(tpl_lib.hdu['FLUX'].shape)
#    pyplot.plot(tpl_lib.hdu['WAVE'].data[101,:], tpl_lib.hdu['FLUX'].data[101,:])
#    pyplot.show()


#-----------------------------------------------------------------------------

if __name__ == '__main__':
    t = clock()

    run_test()

    print('Elapsed time: {0} seconds'.format(clock() - t))



