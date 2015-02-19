# Force Python 3 behavior
from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

# Check long/int definition
import sys
if sys.version > '3':
    long = int

# Imports
import numpy

class mangampl:

        def __init__(self, version=None):

            # Default to MPL-2
            self.mplver = 'MPL-2' if version is None else version

            mpls = self._available_mpls()
            mpli = numpy.where(mpls[:,0] == self.mplver)
            if len(mpli[0]) == 0:
                mpls = self._available_mpls(write=True)
                raise Exception('{0} is not an available MPL!'.format(self.mplver))
                
            self.mplver = str(mpls[mpli].reshape(3)[0])
            self.idlver = str(mpls[mpli].reshape(3)[1])
            self.drpver = str(mpls[mpli].reshape(3)[2])


        def _available_mpls(self, write=False):
            """
            Return a list of the available MPLs to analyze, providing a
            list if requested.
            """

            nmpl = 2
            #                        MPL      IDLUTILS   DRP
            mpl_def = numpy.array([
                                    ['MPL-1', 'v5_5_16', 'v1_0_0'],
                                    ['MPL-2', 'v5_5_17', 'v1_1_2']
                                  ])
    
            if write:
                print('Available MPLs:')
                for x in mpl_def[0:nmpl,:]:
                    print('{0}: IDLUTILS:{1}; DRPVER:{2}'.format(x[0], x[1], x[2]))

            return mpl_def


        def module_file(self):
            """
            Return the name of the module file specific to the MPL to analyze.
            
            TODO: For now this ALWAYS uses the python3 versions.
            """
            if self.mplver == 'MPL-1':
                return 'manga/westfall3_mpl1'
            if self.mplver == 'MPL-2':
                return 'manga/westfall3_mpl2'


        def show(self):
            """
            Print the versions to the screen.
            """

            print('{0}: IDLUTILS:{1}; DRPVER:{2}'.format(self.mplver, self.idlver, self.drpver))
    


