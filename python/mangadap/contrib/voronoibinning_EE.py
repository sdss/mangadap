#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Version 0.0.6 - August 2015 - Included Penalisation on binning, default M=1
# Version 0.0.5 - August 2015 - Change np.average to np.ma.average to avoid Zero
#                               Error message 
# Version 0.0.4 - August 2015 - Added the use of griddata and masks as an option
# Version 0.0.3 - April 2015 - Fixed small bug in the type of data
# Version 0.0.2 - Sep  2012 -  Fixed small bug in output format
# Version 0.0.1 - Jan 2011 - Going from C (Copin/Emsellem) and idl (Cappellari) to python
#
# Copyright (C) 2011  Eric Emsellem
#
#                          Voronoi binning
#
# Code taking, as an input a 2D grid of data+noise points and rebin them
#   to reach a target Signal to Noise.
#
# Freely inspired by codes written by Yannick Copin and Michele Cappellari
# Following the paper by Cappellari, Copin, 2003, MNRAS 342, 345
# 
# From first version written in C in 2000 - 2003
# (mostly by Yannick Copin in the C-Lyon-lib, with modifications by Eric Emsellem)
# Programmes accret_table.c, bin_makeima.c, pack_makeima.c, pack_tiger.c, pack_table.c
# 
# Expanded version written in idl by Michele Cappellari and distributed on the web
# And later translated in Python by Michele Cappellari (also on his web)
#
# The present version was closely inspired by the idl version of Michele Cappellari
#   rewritten and expanded to follow python rules - thanks to its availability
#   on the web.
#
# -----------------------------------------------------------------------------
import numpy as np
from numpy import sum, sqrt, min, max, any
from numpy import argmax, argmin, mean, abs
from numpy import int32 as Nint
from numpy import float32 as Nfloat
import copy

class Error(Exception):
    """Base class for exceptions in this module."""
    pass

class ModuleError(Error):
    """Exception raised for errors in the array sizes
    Attributes:
       msg  -- explanation of the error
    """
    def __init__(self, msg):
        self.msg = msg

    def __str__(self) :
        return repr(self.msg)

def Penalisation_Factor(Nbins=1, M=1.0) :
    """ Penalisation factor to be used when 
    wishing to include systematics in Datacubes
    A standard value is around 20

    Default is 1, namely -> no change
    """
    return sqrt(M * Nbins / (M + Nbins - 1.0))

def bin2d_roundness(x, y, size) :
    """ Return the roundness as defined in Cappellari \& Copin 2003
    """
    npix = x.size
    eR = sqrt(npix / np.pi) * size
    xb = mean(x)
    yb = mean(y)
    maxD = sqrt(max((x-xb)**2 + (y-yb)**2))
    return maxD / eR - 1.0

def dist2(x1,y1,x2,y2, scale=1.0) :
    return ((x1-x2)**2 + (y1-y2)**2)/scale**2

def derive_pixelsize(x, y, verbose=0) :
    """ Find the pixelsize by looking at the minimum distance between
        pairs of x,y
        x: xaxis coordinates
        y: yaxis coordinates
        Return: pixelsize
    """
    pixelsize = 1.e30
    for i in range(len(x)-1) :
        mindist = np.min(dist2(x[i], y[i], x[i+1:], y[i+1:]))
        pixelsize = np.minimum(mindist, pixelsize)
    pixelsize = np.sqrt(pixelsize)
    if verbose:
        print "Pixelsize will be: ", pixelsize
    return pixelsize

def guess_regular_grid(xnodes, ynodes, pixelsize=None) :
    """
    Return a regular grid guessed on an irregular one (Voronoi)
    xnodes, ynodes: arrays of Voronoi bins

    Return: xunb, yunb = regular grid for x and y (unbinned)
    """
    ## First deriving a pixel size
    xn_rav, yn_rav = xnodes.ravel(), ynodes.ravel()
    if pixelsize is None :
        pixelsize = derive_pixelsize(xnodes, ynodes)
    minxn = np.int(np.min(xn_rav) / pixelsize) * pixelsize
    minyn = np.int(np.min(yn_rav) / pixelsize) * pixelsize
    xunb, yunb = np.meshgrid(np.arange(minxn, np.max(xn_rav)+pixelsize, pixelsize), 
                           np.arange(minyn, np.max(yn_rav)+pixelsize, pixelsize))

    return xunb, yunb

def derive_unbinned_field(xnodes, ynodes, data, xunb=None, yunb=None, mask=None, scipy=True) :
    """
       Provide an array of the same shape as the input xunb, and yunb
       with the values derived from the Voronoi binned data

       xnodes, ynodes: 2 arrays providing the nodes from the binning
       data : values for each node
       xunb, yunb: x and y coordinates of the unbinned data
                 if not provided (default) they will be guessed from the nodes
       mask: array with the same shape than xunb providing mask values for positions to ignore
       scipy: by default using scipy.griddata. If not True, using a loop on bins.

       Return: xunb, yunb, and unbinned_data arrays with the same shape as xunb,
    """
    if xunb is None :
        xunb, yunb = guess_regular_grid(xnodes, ynodes)

    x_rav, y_rav = xunb.ravel(), yunb.ravel()
    xnodes_rav, ynodes_rav = xnodes.ravel(), ynodes.ravel()
    data_rav = data.ravel()
    unbinned_data = np.zeros_like(x_rav, dtype=data_rav.dtype)

    ## Getting the mask array
    if mask == None : 
        mask_ravel = np.zeros_like(x_rav, dtype=bool)
    else :
        mask_ravel = mask.mask.ravel()
    
    ## If scipy is True (default) use griddata
    if scipy :
        import scipy
        unbinned_data[-mask_ravel] = scipy.interpolate.griddata((xnodes_rav, ynodes_rav), data_rav, 
             (x_rav[-mask_ravel], y_rav[-mask_ravel]), method='nearest')
    ## Otherwise loop through the points
    else :
        interm_data = np.zeros_like(x_rav[-mask_ravel])
        sizein_M = len(x_rav[-mask_ravel])
        sizeout = len(xnodes_rave)
        indclosestBin = np.argmin(dist2(x_rav[-mask_ravel].reshape(1,sizein_M),
           y_rav[-mask_ravel].reshape(1, sizein_M), xnodes_rav.reshape(sizeout, 1),
           ynodes_rav.reshape(sizeout, 1)), axis=0)
        unbinned_data[-mask_ravel] = data_rav[indclosestBin]

#         for i in xrange(len(x_rav[-mask_ravel])) :
#             indclosestBin = argmin(dist2(x_rav[-mask_ravel][i], y_rav[-mask_ravel][i], xnodes_rav, ynodes_rav))
#             interm_data[i] = data_rav[indclosestBin]
#         unbinned_data[-mask_ravel] = interm_data
    
    return xunb, yunb, unbinned_data.reshape(xunb.shape)

listmethods = ["voronoi", "quadtree"]
class bin2D :
    """ 
    Class for Voronoi binning of a set of x and y coordinates
    using given data and potential associated noise
    """
    def __init__(self, xin, yin, data, noise=None, targetSN=1.0, pixelsize=None, method="Voronoi", cvt=1, wvt=0, M_Penalise=1.0) :
        self.xin = xin.ravel()
        self.yin = yin.ravel()
        self.data = data.ravel()
        self.noise = noise.ravel()
        self.SN = np.where(self.noise > 0, self.data/self.noise, 0.0)

        ## Sort all pixels by their distance to the maximum SN
        if pixelsize is None :
            self.pixelsize = derive_pixelsize(self.xin, self.yin, verbose=1)
        else :
            self.pixelsize = pixelsize

        self.targetSN = targetSN

        ## Binning method and options
        self.method = str.lower(method)
        self.cvt = cvt
        self.wvt = wvt
        self.scale = 1.0
        self.M_Penalise = M_Penalise

        self._check_input()
        self._check_data()

    def _warning(self, text) :
        """ 
        Warning message for 2D Binning class
        """
        print "WARNING [2D Binning]: %s"%(text)

    def _error(self, text) :
        """ 
        Error message for 2D Binning class
        Exit after message
        """
        print "ERROR [2D Binning]: %s"%(text)
        return

    def _check_input(self) :
        """ 
        Check consistency of input data
        """
        if self.noise is None :
            if any(self.data < 0) :
                self._error("No Noise given, and some pixels have negative data value")
            ## If noise is not there then just use SQRT(data)
            self.noise == sqrt(self.data)

        # Basic checks of the data
        # First about the dimensions of the datasets
        self.npix = len(self.xin)
        if self.npix == 0 :
            raise ModuleError("xin as 0 pixels")

        if (len(self.yin) != self.npix) | (len(self.data) != self.npix) :
            raise ModuleError("x, y and data do not have the same size")
            return

        try :
            if self.targetSN <= 0 :
                self._error("targetSN must be positive")
        except :
            self._error("targetSN must be a number")

        if self.method not in listmethods :
            self._error("Given method [%s] not supported"%(self.method))

    ### Check input data prior to binning ============================================
    def _check_data(self) :
        """ 
        Check data before binning (noise versus signal)
        """
        ## Check if enough S/N in the sum of All pixels
        if sum(self.data) / sqrt(sum(self.noise**2)) < self.targetSN :
            self._warning("Not enough S/N in data for this targetSN")
            return

        ## Check of the any binning is needed
        if min(self.data/self.noise) > self.targetSN :
            self._warning("All pixels have enough S/N - no binning performed")
            self.xb = self.xin
            self.yb = self.yin
            self.datab = self.data
            self.noiseb = self.noise
            self.SNb = np.where(self.noiseb > 0, self.datab/self.noiseb, 0.0)

    ### Signal to Noise ============================================
    def SN_indexArray(self, listbins) :
        """ Return the SN array given an array of data values and noise values
        taking into account a list of bins
        This includes the Penalisation function with the M factor given in self
        """
        return sum(self.data[listbins]) / (Penalisation_Factor(len(listbins), self.M_Penalise) * sqrt(sum(self.noise[listbins])**2))

    ### Do Quadtree binning ============================================
    def bin_quadtree(self, verbose=1) :
        """ Actually do the Quadreee binning if that method is chosen
        Not implemented YET
        """
        pass

    ### Accrete the bins ============================================
    def bin2d_accretion(self, verbose=0) :
        """ Accrete the bins according to their SN
        """
        ## Status holds the bin number when assigned
        self.status = np.zeros(self.npix, dtype=Nint)
        ## Good is 1 if the bin was accepted
        self.good = np.zeros(self.xin.size, dtype=Nint)

        ## Start with the max SN in the data
        currentBin = [argmax(self.SN)]

        for ind in range(1,self.npix+1) :  ## Running over the index of the Voronoi BIN
            ## Only one pixel at this stage
            currentSN = self.SN[currentBin]
            if verbose : print "Bin %d"%(ind)

            self.status[currentBin] = ind   # only one pixel at this stage
            ## Barycentric centroid for 1 pixel...
            xbar, ybar = self.xin[currentBin], self.yin[currentBin]

            ## Indices of remaining unbinned data
            unBinned = np.where(self.status == 0)[0]
            ##+++++++++++++++++++++++++++++++++++++++++++++++++
            ## STOP THE WHILE Loop if all pixels are now binned
            while len(unBinned) > 0 :
            ##+++++++++++++++++++++++++++++++++++++++++++++++++

                ## Coordinates of the Remaining unbinned pixels
                xunBin = self.xin[unBinned]
                yunBin = self.yin[unBinned]

                ## Location of current bin
                xcurrent = self.xin[currentBin]
                ycurrent = self.yin[currentBin]

                ## Closest unbinned pixel to the centroid making sure the highest SN is used
                indclosestBar = argmin(dist2(xbar,ybar, xunBin, yunBin))
                xclosest = xunBin[indclosestBar]
                yclosest = yunBin[indclosestBar]

                ## Distance between this closest pixel and the current pixel
                currentSqrtDist= sqrt(min(dist2(xclosest,yclosest, xcurrent,ycurrent)))

                ## Testing the roundness for this potential new Bin
                possibleBin = currentBin + [unBinned[indclosestBar]]
                roundness = bin2d_roundness(self.xin[possibleBin], self.yin[possibleBin], self.pixelsize)

                ## Transfer new SN to current value
                oldSN = currentSN
                currentSN = self.SN_indexArray(possibleBin)

                ##++++++++++++++++++++++++++++++++++++++++++
                ## Test the new potential Bin
                ## Criterion 1 = Connected?
                ## Criterion 2 = round enough?
                ## Criterion 3 = getting better for SN?
                ##++++++++++++++++++++++++++++++++++++++++++
                ## If the pixel is not ok  for one of these reasons
                ##    The accretion process stops and we test if the new
                ##    bin has the right SN
                ##++++++++++++++++++++++++++++++++++++++++++
                if (currentSqrtDist > 1.2 * self.pixelsize) | (roundness > 0.3) | (abs(currentSN-self.targetSN) > abs(oldSN - self.targetSN)) :
                    if (oldSN > 0.8 * self.targetSN) :
                        self.good[currentBin] = 1
                    break
                ##++++++++++++++++++++++++++++++++++++++++++

                ## If the new Bin is ok we associate the number of the bin to that one
                self.status[unBinned[indclosestBar]] = ind
                ##   ... and we use that now as the current Bin
                currentBin = possibleBin

                ## And update the values
                ## First the centroid
                xbar, ybar = mean(self.xin[currentBin]), mean(self.yin[currentBin])

                ## New set of unbinned pixels
                unBinned = np.where(self.status == 0)[0]
                ## ----- End of While Loop --------------------------------------------

            ## Unbinned pixels
            unBinned = np.where(self.status == 0)[0]
            ## Break if all pixels are binned
            if len(unBinned) == 0 :
                break
            ## When the while loop is finished for this new BIN
            ## Find the centroid of all Binned pixels
            Binned = np.where(self.status != 0)[0]
            xbar, ybar = mean(self.xin[Binned]), mean(self.yin[Binned])

            ## Closest unbinned pixel to the centroid of all Binned pixels
            xunBin = self.xin[unBinned]
            yunBin = self.yin[unBinned]
            indclosestBar = argmin(dist2(xbar,ybar, xunBin, yunBin))
            ## Take now this closest pixel as the next one to look at
            currentBin = [unBinned[indclosestBar]]

        ## Set to zero all bins that did not reach the target SN
        self.status *= self.good

    ### Compute centroid of bins ======================================
    def bin2d_centroid(self, verbose=0):
        ## Compute the area for each node
        self.Areanode, self.bins = np.histogram(self.status, bins=np.arange(np.max(self.status+0.5))+0.5)
        ## Select the ones which have at least one bin
        self.indgoodbins = np.where(self.Areanode > 0)[0]
        self.Areanode = self.Areanode[self.indgoodbins]
        ngoodbins = self.indgoodbins.size
        ## Reset the xnode, ynode, SNnode, and statusnode (which provides the number for the node)
        self.xnode = np.zeros(ngoodbins, Nfloat)
        self.ynode = np.zeros_like(self.xnode)
        self.SNnode = np.zeros_like(self.xnode)
        self.statusnode = np.zeros_like(self.xnode)
        self.listbins = []
        for i in range(ngoodbins) :
            ## indgoodbins[i] provides the bin of the Areanode, so indgoodbins[i] + 1 is the status number
            self.statusnode[i] = self.indgoodbins[i]+1
            ## List of bins which have statusnode as a status
            listbins = np.where(self.status==self.statusnode[i])[0]
            ## Centroid of the node
            self.xnode[i], self.ynode[i] = mean(self.xin[listbins]), mean(self.yin[listbins])
            self.SNnode[i] = self.SN_indexArray(listbins)
            self.listbins.append(listbins)

    ### Compute WEIGHTED centroid of bins ======================================
    def bin2d_weighted_centroid(self, weight=None, verbose=0):
        if weight is not None : self.weight = weight

        self.Areanode, self.bins = np.histogram(self.status, bins=np.arange(np.max(self.status+0.5))+0.5)
        ## Select the ones which have at least one bin
        self.indgoodbins = np.where(self.Areanode > 0)[0]
        self.Areanode = self.Areanode[self.indgoodbins]
        ngoodbins = self.indgoodbins.size
        ## Reset the xnode, ynode, SNnode, and statusnode (which provides the number for the node)
        self.xnode = np.zeros(ngoodbins, Nfloat)
        self.ynode = np.zeros_like(self.xnode)
        self.SNnode = np.zeros_like(self.xnode)
        self.statusnode = np.zeros_like(self.xnode)
        self.listbins = []
        for i in range(ngoodbins) :
            ## indgoodbins[i] provides the bin of the Areanode, so indgoodbins[i] + 1 is the status number
            self.statusnode[i] = self.indgoodbins[i]+1
            ## List of bins which have statusnode as a status
            listbins = np.where(self.status==self.statusnode[i])[0]
            ## Weighted centroid of the node
            self.xnode[i], self.ynode[i] = np.ma.average(self.xin[listbins], weights=self.weight[listbins]), np.ma.average(self.yin[listbins], weights=self.weight[listbins])
            self.SNnode[i] = self.SN_indexArray(listbins)
            self.listbins.append(listbins)

    ### Assign bins ============================================
    def bin2d_assign_bins(self, sel_pixels=None, scale=None, verbose=0) :
        """ 
        Assign the bins when the nodes are derived With Scaling factor
        """
        if scale is not None: self.scale = scale
        if sel_pixels is None : sel_pixels = range(self.xin.size)
        for i in sel_pixels :
            minind = argmin(dist2(self.xin[i], self.yin[i], self.xnode, self.ynode, scale=self.scale))
            self.status[i] = self.statusnode[minind]
            if verbose :
                print "Pixel ",  self.status[i], self.xin[i], self.yin[i], self.xnode[minind], self.ynode[minind]

        ## reDerive the centroid
        self.bin2d_centroid()

    ### Do  CV tesselation ============================================
    def bin2d_cvt_equal_mass(self, wvt=None, verbose=1) :
        """ 
        Produce a CV Tesselation

        wvt: default is None (will use preset value, see self.wvt)
        """

        ## Reset the status and statusnode for all nodes
        self.status = np.zeros(self.npix, dtype=Nint)
        self.statusnode = np.arange(self.xnode.size) + 1

        if wvt is not None : self.wvt = wvt
        if self.wvt : self.weight = np.ones_like(self.SN)
        else : self.weight = self.SN**4

        self.scale = 1.0

        self.niter = 0
        ## WHILE LOOP: stop when the nodes do not move anymore ============
        Oldxnode, Oldynode = copy.copy(self.xnode[-1]), copy.copy(self.ynode[-1])
        while (not np.array_equiv(self.xnode, Oldxnode)) | (not np.array_equiv(self.ynode, Oldynode)):
            Oldxnode, Oldynode = copy.copy(self.xnode), copy.copy(self.ynode)
            ## Assign the closest centroid to each bin
            self.bin2d_assign_bins()

            ## New nodes weighted centroids
            self.bin2d_weighted_centroid()

            ## Eq. (4) of Diehl & Statler (2006)
            if self.wvt : self.scale = sqrt(self.Areanode/self.SNnode)
            self.niter += 1

    ### Do Voronoi binning ============================================
    def bin_voronoi(self, wvt=None, cvt=None, verbose=1) :
        """ Actually do the Voronoi binning
        
        wvt: default is None (will use preset value, see self.wvt)
        cvt: default is None (will use preset value, see self.cvt)
        """
        if cvt is not None : self.cvt = cvt
        if wvt is not None : self.wvt = wvt

        print "=================="
        print "Accreting Bins... "
        self.bin2d_accretion()
        print "          ...Done"
        print "===================="
        print "Reassigning Bins... "
        self.bin2d_centroid()
        ## Get the bad pixels, not assigned and assign them
        badpixels = np.where(self.status == 0)[0]
        self.bin2d_assign_bins(badpixels)
        print "            ...Done"
        print "===================="
        if self.cvt :
            print "==========================="
            print "Modified Lloyd algorithm..."
            self.bin2d_cvt_equal_mass()
            print "%d iterations Done."%(self.niter)
            print "==========================="
        else : self.scale = 1.0
        ## Final nodes weighted centroids after assigning to the final nodes
        self.bin2d_assign_bins()
        if self.wvt : self.weight = np.ones_like(self.data)
        else : self.weight = self.data
        self.bin2d_weighted_centroid()

    def show_voronoibin(self, datain=None, shownode=1, mycmap=None) :
        """
        Display the voronoi bins on a map
        
        datain: if None (Default), will use random colors to display the bins
                if provided, will display that with a jet (or specified mycmap) cmap 
                   (should be either the length of the voronoi nodes array or the size of the initial pixels)
        shownode: default is 1 -> show the voronoi nodes, otherwise ignore (0)
        mycmap: in case datain is provide, will use that cmpa to display the bins
        """
        from distutils import version
        try:
            import matplotlib
        except ImportError:
            raise Exception("matplotlib 0.99.0 or later is required for this routine")

        if version.LooseVersion(matplotlib.__version__) < version.LooseVersion('0.99.0'):
            raise Exception("matplotlib 0.99.0 or later is required for this routine")

        from matplotlib.collections import PatchCollection
        import matplotlib.patches as mpatches
        import matplotlib.pyplot as plt
        fig = plt.figure(1,figsize=(7,7))
        plt.clf()
        ax = plt.gca()
        patches = []
        binsize = self.pixelsize
        for i in range(len(self.xin)) :
            patches.append(mpatches.Rectangle((self.xin[i],self.yin[i]), binsize, binsize))

        if datain is None :
            dataout = self.status
            mycmap = 'prism'
        else :
            if len(datain) == self.xnode.size :
                dataout = np.zeros(self.xin.size, Nfloat)
                for i in range(self.xnode.size) :
                    listbins = self.listbins[i]
                    dataout[listbins] = [datain[i]]*len(listbins)
            elif len(datain) == self.xin.size :
                dataout = datain
            if mycmap is None : mycmap = 'jet'

        colors = dataout * 100.0 / max(dataout)
        collection = PatchCollection(patches, cmap=mycmap)
        collection.set_array(np.array(colors))
        ax.add_collection(collection)
        if shownode :
            plt.scatter(self.xnode,self.ynode, marker='o', edgecolors='w', facecolors='none')
        plt.axis('image')
        plt.xlabel("X axis")
        plt.ylabel("Y axis")
        plt.title("Voronoi Map")

    def bin_data(self, datain=None, noisein=None) :
        """
        Return a Voronoi adaptive binning of your data.

        datain: if provided, will be used as data input
                if not provided (None = default), will use self.data
        noisein: if provided, will be used as noise input
                if not provided (None = default), will use self.noise

        Output = xnode, ynode, bindata, S/N
        """

        if datain is None: datain = copy.copy(self.data)
        if noisein is None: noisein = copy.copy(self.noise)

        dataout = np.zeros(self.xnode.size, Nfloat)
        xout = np.zeros_like(dataout)
        yout = np.zeros_like(dataout)
        SNout = np.zeros_like(dataout)
        for i in range(self.xnode.size) :
            listbins = self.listbins[i]
            xout[i] = np.ma.average(self.xin.ravel()[listbins], weights=datain[listbins])
            yout[i] = np.ma.average(self.yin.ravel()[listbins], weights=datain[listbins])
            dataout[i] = mean(datain[listbins])
            SNout[i] = self.SN_indexArray(listbins)

        return xout, yout, dataout, SNout
