#!/usr/bin/env python3

import numpy
from os import remove
import os.path
from time import clock

#-----------------------------------------------------------------------------

class BasePar:
    def __init__(self, inp):
        print('running base init')
        print('')
        self.pars = inp
        self.write()
    def write(self):
        print('running base write')
        print(self.pars)
        print('')

class DerivedPar(BasePar):
    def __init__(self, item1, item2, item3):
        print('running derived init')
        print('')
        inp = [ item1, item2, item3 ]
        BasePar.__init__(self, inp)

def run_test():
    x = DerivedPar('test', 1.0, 0)
    x.write()
    print('')

    x.pars[0] = 'derivedtest'
    x.write()
    print('')

    y = BasePar(['basetest', 1.0, 0])
    y.write()
    print('')

#-----------------------------------------------------------------------------

if __name__ == '__main__':
    t = clock()
    run_test()
    print('Elapsed time: {0} seconds'.format(clock() - t))



