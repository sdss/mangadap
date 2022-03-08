#!/usr/bin/env python3

"""
Dynamically build the rst documentation of the bitmasks.
"""

import os
import time
import numpy

#-----------------------------------------------------------------------------

if __name__ == '__main__':
    t = time.perf_counter()

    from pkg_resources import resource_filename
    root = os.path.dirname(resource_filename('mangadap', ''))
    path = os.path.join(root, 'docs', 'tables')
    if not os.path.isdir(path):
        os.makedirs(path)

    from mangadap.config.analysisplan import AnalysisPlanSet
    from mangadap.par.parset import ParSet

    ap = AnalysisPlanSet.default()

    data = numpy.empty((ap.npar//2+1,2), dtype=object)
    data[0] = ['Key', 'Value']
    j = 0
    for n in ap.keys():
        if n == 'key' or 'clobber' in n:
            continue
        data[j+1,0] = n
        data[j+1,1] = str(ap[0][n])
        j += 1
    lines = ParSet._data_table_string(data, delimiter='rst')

    ofile = os.path.join(path, 'default_analysisplan.rst')
    with open(ofile, 'w') as f:
        f.write(lines)

    print('Elapsed time: {0} seconds'.format(time.perf_counter() - t))



