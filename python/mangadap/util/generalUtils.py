#!usr/bin/env python

'''
    NAME:
        generalUtils

    METHODS:
        getMangaVersions - Retrieves the Manga DRP pipeline version number being used for current reductions.
        getSlitmap - Retrieves a slitmap
        traceset2xy - Converts a traceset into a 2d array of y-positions
        
    REVISION HISTORY:
        18-Apr-2014 - Written by Brian Cherinka (Toronto)
        18-Apr-2014 - Added getMangaVersion
        19-Apr-2014 - Added getSlitmap
        19-Apr-2014 - Added traceset2xy
            
'''

from __future__ import division
from __future__ import print_function
import numpy as np, os
import subprocess,glob

from mangadap.util.yanny import yanny, read_yanny

# 13 Feb 2015: K. Westfall, python3 compatibilty
import sys
if sys.version > '3':
    long = int

__author__='Brian Cherinka'

def getMangaVersion(simple=False, drp=False, core=False):
    '''Gets the Manga Pipeline Version for the DRP or mangaCORE.  Default is DRP.'''
    
    if drp: 
        ver=subprocess.check_output('mangadrp_version',shell=True)
    elif core:
        ver=subprocess.check_output('mangacore_version',shell=True)
    else:
        ver=subprocess.check_output('mangadrp_version',shell=True)    
    
    if simple:
        return ver.split(' ')[0].rstrip('\n')
    else:
        return ver.rstrip('\n') 

def getSlitmap(filename=None, plate=None, mjd=None, specid=None):
    '''Gets the slitmap from mangacore'''
    
    # Read slitmap
    if filename != None:
        fullslitmap = read_yanny(filename)
    else:
        if plate == None:
            print('Must provide the Plate')
            return None
        if mjd == None: mjd='*'              
        name = 'slitmap-{0}-{1}-*.par'.format(plate,mjd)
        if long(plate) >= 10000: pltgroup = '0'+plate[:-2]+'XX'
        if long(plate) <= 9999: pltgroup = '00'+plate[:-2]+'XX'
        slitpath = os.path.join(os.getenv('MANGACORE_DIR'),'slitmaps',pltgroup,str(plate))
        filename = max(glob.glob(os.path.join(slitpath,name)))
        fullslitmap = read_yanny(filename)
        
    # Grab slitmap structure    
    slitmap = fullslitmap['SLITMAP']
    
    # Split on specid if possible
    if specid != None:
        endone = slitmap['spectrographid'].index(2)-1
        if int(specid) == 1: subset={key:slitmap[key][0:endone+1] for key in slitmap.keys()}
        if int(specid) == 2: subset={key:slitmap[key][endone+1:] for key in slitmap.keys()}
        slitmap = subset
    
    return slitmap
    
def traceset2xy(tset):
    '''Convert a traceset to x,y positions'''
    
    ndim = len(tset[0][3].shape)
    dims = tset[0][3].shape
    
    if ndim == 1:
        ncoeff = dims[1]
        ntrace = 1
    elif ndim == 2:
        ncoeff = dims[1]
        ntrace = dims[0]    
    
    nx = (tset[0][2]-tset[0][1]+1)
    xmid = 0.5 * (tset[0][1]+tset[0][2])
    xrange = tset[0][2]-tset[0][1]
    
    xpos = np.arange(nx*ntrace,dtype=long).reshape(ntrace,nx) % nx
    ypos = xpos*0.0
    
    for i in range(ntrace):
        xinput = xpos[i]
        xvec = 2.0 * (xinput - xmid)/xrange
        ypos[i] = np.polynomial.legendre.legval(xvec,tset[0][3][i])
        
    arr = np.power(10.,ypos)
        
    return arr 
    
def getPlateGrp(plate):
    '''Plate group - 0073XX, useful for MANGACORE directory structure'''
    
    plate = str(plate)
    if long(plate) >= 10000: pltgroup = '0'+plate[:-2]+'XX'
    if long(plate) <= 9999: pltgroup = '00'+plate[:-2]+'XX'

    return pltgroup
    
def getNullVal(typ):
    ''' Return a standard null value based on type of data '''
    
    if typ == str: nullval = 'NULL'
    if typ == int: nullval = -999
    if typ == long: nullval = long(-999)
    if typ == float: nullval = -999.0
    
    return nullval

def getModules(name='manga'):
    ''' Return the module versions for the specified name'''    
    
    try: modules = os.environ['LOADEDMODULES']
    except: modules = None
    
    if modules:
        versions = [module.split('/')[1] for module in modules.split(':') if module.split('/')[0]==name] if modules else []
        if len(versions) == 1: version = versions[0]
        elif len(versions):
            print('Multiple versions found for module {0}'.format(name))
            version = None 
        else: 
            print('Module {0} is not loaded'.format(name))
            version = None
    else: version = None
    
    return version
    
def readSDSSMaskBits(name=None):
    ''' Read the SDSS Mask Bits from the yanny file and return a dictionary'''
    
    path = os.path.join(os.getenv('IDLUTILS_DIR'),'data','sdss','sdssMaskbits.par')

    print(path)
    
    data = read_yanny(path)['MASKBITS']
    flags = data['flag']
    
    # Select out named subset
    if name:
        data={key:[v for i,v in enumerate(val) if name in flags[i]] for key,val in data.items()}    
# 13 Feb 2015: K. Westfall, python3 compatibilty
#       data={key:[v for i,v in enumerate(val) if name in flags[i]] for key,val in data.iteritems()}
    
    return data
    
def getSDSSFlagName(bits, name=None, data=None):
    ''' Retrieve the flag names from a bit flag '''
   
    if data is None:
        data = readSDSSMaskBits(name=name)
    
    # if bits not a digit, return None
    if not str(bits).isdigit(): return 'NULL'
    else: bits=int(bits)
    
    # Convert the integer value to list of bits
    bitlist=[int(i) for i in '{0:08b}'.format(bits)]
    bitlist.reverse()
    indices = [i for i,bit in enumerate(bitlist) if bit]
    
    # Make new dictionary for mask bits
    bitdict = {key:[v for i,v in enumerate(val) if i in indices] for key,val in data.items()}
# 13 Feb 2015: K. Westfall, python3 compatibilty
#   bitdict = {key:[v for i,v in enumerate(val) if i in indices] for key,val in data.iteritems()}
    
    return bitdict['label']
    
    
    
    
    
    
    
    
    
    
