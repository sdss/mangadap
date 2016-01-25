#!/usr/bin/env python3

import numpy
from os import remove
import os.path
from time import clock

from mangadap.util.par import ParList

#-----------------------------------------------------------------------------

def run_test():

    keys = [ 'test', 'par', 'list']

    defaults = [ 'this', 0, 0.0 ]

    options = [ 'this', None, None ]

    dtypes = [ str, int, [int, float] ]

    par = ParList(keys, defaults=defaults, options=options, dtypes=dtypes)

    for p in par:
        print(p)

    print(par.dtype)
    
    print(par['test'])
    print(par['par'])
    print(par['list'])

    par['test'] = 'this'

    par['list'] = 3

    print(par['test'])
    print(par['par'])
    print(par['list'])

    par.add_new('echo', 10, dtype=int)

    for p in par:
        print(p)
    print(par.data)

    par['echo'] = 1.3


#-----------------------------------------------------------------------------

if __name__ == '__main__':
    t = clock()
    run_test()
    print('Elapsed time: {0} seconds'.format(clock() - t))



