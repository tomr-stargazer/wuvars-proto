'''
Here's some Robust statistical functions.

Robust Mean, Robust Standard Deviation, and Robust Variance 
borrowed from Jiangang Hao at:
http://code.google.com/p/fmccd/source/browse/trunk/fermiMCCD.py?spec=svn18&r=18

Other versions from Ian Crossfield at:
http://astro.ucla.edu/~ianc/python/analysis.html
Copyright 2010, Ian Crossfield

'''

#import numpy as np
from numpy import *

#-----Robust Mean--------------
#-----from J. Hao--------------
def robust_mean(x):
    y = x.flatten()
    n = len(y)
    y.sort()
    ind_qt1 = round((n+1)/4.)
    ind_qt3 = round((n+1)*3/4.)
    IQR = y[ind_qt3]- y[ind_qt1]
    lowFense = y[ind_qt1] - 1.5*IQR
    highFense = y[ind_qt3] + 1.5*IQR
    ok = (y>lowFense)*(y<highFense)
    yy=y[ok]
    return yy.mean(dtype='double')

#-------Robust Standard Deviation---
#-----from J. Hao--------------
def robust_std(x):
    y = x.flatten()
    n = len(y)
    y.sort()
    ind_qt1 = round((n+1)/4.)
    ind_qt3 = round((n+1)*3/4.)
    IQR = y[ind_qt3]- y[ind_qt1]
    lowFense = y[ind_qt1] - 1.5*IQR
    highFense = y[ind_qt3] + 1.5*IQR
    ok = (y>lowFense)*(y<highFense)
    yy=y[ok]
    return yy.std(dtype='double')

#-------Robust variance---
#-----from J. Hao--------------
def robust_var(x):
    y = x.flatten()
    n = len(y)
    y.sort()
    ind_qt1 = round((n+1)/4.)
    ind_qt3 = round((n+1)*3/4.)
    IQR = y[ind_qt3]- y[ind_qt1]
    lowFense = y[ind_qt1] - 1.5*IQR
    highFense = y[ind_qt3] + 1.5*IQR
    ok = (y>lowFense)*(y<highFense)
    yy=y[ok]
    return yy.var(dtype='double')

#      Stuff from Ian Crossfield (it's better)

def removeoutliers(data, nsigma, remove='both', center='mean', 
                   niter=Inf, retind=False, verbose=False):
    """Strip outliers from a dataset, iterating until converged.

    INPUT:
      data -- 1D numpy array.  data from which to remove outliers.
      nsigma -- positive number.  limit defining outliers: number of
                standard deviations from center of data.

    OPTIONAL INPUTS:               
      remove -- ('min'|'max'|'both') respectively removes outliers
                 below, above, or on both sides of the limits set by
                 nsigma.
      center -- ('mean'|'median'|value) -- set central value, or
                 method to compute it.
      niter -- number of iterations before exit; defaults to Inf,
               which can occasionally result in empty arrays returned
      retind -- (bool) whether to return index of good values as
                second part of a 2-tuple.

    EXAMPLE: 
       from numpy import hist, linspace, randn
       from analysis import removeoutliers
       data = randn(1000)
       hbins = linspace(-5,5,50)
       d2 = removeoutliers(data, 1.5, niter=1)
       hist(data, hbins)
       hist(d2, hbins)
       """
    # 2009-09-04 13:24 IJC: Created
    # 2009-09-24 17:34 IJC: Added 'retind' feature.  Tricky, but nice!
    # 2009-10-01 10:40 IJC: Added check for stdev==0
    # 2009-12-08 15:42 IJC: Added check for isfinite

    from numpy import median, ones, isfinite
    from pylab import find

    def getcen(data, method):
        "Get central value of a 1D array (helper function)"
        if method.__class__==str:
            if method=='median':
                cen = median(data)
            else:
                cen = data.mean()
        else:
            cen = method
        return cen

    def getgoodindex(data, nsigma, center, stdev, remove):
        "Get number of outliers (helper function!)"
        if stdev==0:
            distance = data*0.0
        else:
            distance = (data-center)/stdev
        if remove=='min':
            goodind = distance>-nsigma
        elif remove=='max':
            goodind = distance<nsigma
        else:
            goodind = abs(distance)<=nsigma
        return goodind

    data = data.ravel().copy()

    ndat0 = len(data)
    ndat = len(data)
    iter=0
    goodind = ones(data.shape,bool)
    goodind *= isfinite(data)
    while ((ndat0<>ndat) or (iter==0)) and (iter<niter) and (ndat>0) :
        ndat0 = len(data[goodind])
        cen = getcen(data[goodind], center)
        stdev = data[goodind].std()
        thisgoodind = getgoodindex(data[goodind], nsigma, cen, stdev, remove)
        goodind[find(goodind)] = thisgoodind
        if verbose:
            print "cen>>",cen
            print "std>>",stdev
        ndat = len(data[goodind])
        iter +=1
        if verbose:
            print ndat0, ndat
    if retind:
        ret = data[goodind], goodind
    else:
        ret = data[goodind]
    return ret


def meanr(x, nsigma=3, niter=Inf, finite=True, verbose=False,axis=None):
    """Return the mean of an array after removing outliers.
    
    INPUTS:
      x -- (array) data set to find mean of

    OPTIONAL INPUT:
      nsigma -- (float) number of standard deviations for clipping
      niter -- number of iterations.
      finite -- if True, remove all non-finite elements (e.g. Inf, NaN)
      axis -- (int) axis along which to compute the mean.

    EXAMPLE:
      from numpy import *
      from analysis import meanr
      x = concatenate((randn(200),[1000]))
      print mean(x), meanr(x, nsigma=3)
      x = concatenate((x,[nan,inf]))
      print mean(x), meanr(x, nsigma=3)

    SEE ALSO: medianr, stdr, removeoutliers, numpy.isfinite
    """
    # 2009-10-01 10:44 IJC: Created
    # 2010-07-01 13:52 IJC: Now handles higher dimensions.
    from numpy import array, isfinite, zeros

    x = array(x, copy=True)
    xshape = x.shape
    ndim =  len(xshape)
    if ndim==0:
        return x

    if  ndim==1 or axis is None:
        # "1D" array
        x = x.ravel()
        if finite:
            x = x[isfinite(x)]
        x = removeoutliers(x, nsigma, niter=inf, verbose=verbose)
        return x.mean()
    else:
        newshape = list(xshape)
        oldDimension = newshape.pop(axis) 
        ret = zeros(newshape, float)

        # Prevent us from taking the mean along the axis of primary incidices:
        if axis==0: 
            x=swapaxes(x,0,1)

        if axis>1:
            nextaxis = axis-1
        else:
            nextaxis = 0

        for ii in range(newshape[0]):
            #print 'x.shape>>',x.shape, 'newshape>>',newshape, 'x[ii].shape>>',x[ii].shape, 'ret[ii].shape>>',ret[ii].shape,'ii>>',ii
            ret[ii] = meanr(x[ii], nsigma=nsigma,niter=niter,finite=finite,\
                                 verbose=verbose, axis=nextaxis)
        return ret

def medianr(x, nsigma=3, niter=Inf, finite=True, verbose=False,axis=None):
    """Return the median of an array after removing outliers.
    
    INPUTS:
      x -- (array) data set to find median of

    OPTIONAL INPUT:
      nsigma -- (float) number of standard deviations for clipping
      niter -- number of iterations.
      finite -- if True, remove all non-finite elements (e.g. Inf, NaN)
      axis -- (int) axis along which to compute the mean.

    EXAMPLE:
      from numpy import *
      from analysis import medianr
      x = concatenate((randn(200),[1000]))
      print median(x), medianr(x, nsigma=3)
      x = concatenate((x,[nan,inf]))
      print median(x), medianr(x, nsigma=3)

    SEE ALSO: meanr, stdr, removeoutliers, numpy.isfinite
    """
    # 2009-10-01 10:44 IJC: Created
    #2010-07-01 14:04 IJC: Added support for higher dimensions
    from numpy import array, isfinite, zeros, median

    x = array(x, copy=True)
    xshape = x.shape
    ndim =  len(xshape)
    if ndim==0:
        return x

    if  ndim==1 or axis is None:
        # "1D" array
        x = x.ravel()
        if finite:
            x = x[isfinite(x)]
        x = removeoutliers(x, nsigma, niter=inf, verbose=verbose)
        return median(x)
    else:
        newshape = list(xshape)
        oldDimension = newshape.pop(axis) 
        ret = zeros(newshape, float)

        # Prevent us from taking the action along the axis of primary incidices:
        if axis==0: 
            x=swapaxes(x,0,1)

        if axis>1:
            nextaxis = axis-1
        else:
            nextaxis = 0

        for ii in range(newshape[0]):
            #print 'x.shape>>',x.shape, 'newshape>>',newshape, 'x[ii].shape>>',x[ii].shape, 'ret[ii].shape>>',ret[ii].shape,'ii>>',ii
            ret[ii] = medianr(x[ii], nsigma=nsigma,niter=niter,finite=finite,\
                                 verbose=verbose, axis=nextaxis)
        return ret

def stdr(x, nsigma=3, niter=Inf, finite=True, verbose=False, axis=None):
    """Return the standard deviation of an array after removing outliers.
    
    INPUTS:
      x -- (array) data set to find std of

    OPTIONAL INPUT:
      nsigma -- (float) number of standard deviations for clipping
      niter -- number of iterations.
      finite -- if True, remove all non-finite elements (e.g. Inf, NaN)
      axis -- (int) axis along which to compute the mean.

    EXAMPLE:
      from numpy import *
      from analysis import stdr
      x = concatenate((randn(200),[1000]))
      print std(x), stdr(x, nsigma=3)
      x = concatenate((x,[nan,inf]))
      print std(x), stdr(x, nsigma=3)

    SEE ALSO: meanr, medianr, removeoutliers, numpy.isfinite
    """
    # 2010-02-16 14:57 IJC: Created from mear
    # 2010-07-01 14:06 IJC: ADded support for higher dimensions
    from numpy import array, isfinite, zeros

    x = array(x, copy=True)
    xshape = x.shape
    ndim =  len(xshape)
    if ndim==0:
        return x

    if  ndim==1 or axis is None:
        # "1D" array
        x = x.ravel()
        if finite:
            x = x[isfinite(x)]
        x = removeoutliers(x, nsigma, niter=inf, verbose=verbose)
        return x.std()
    else:
        newshape = list(xshape)
        oldDimension = newshape.pop(axis) 
        ret = zeros(newshape, float)

        # Prevent us from taking the action along the axis of primary incidices:
        if axis==0: 
            x=swapaxes(x,0,1)

        if axis>1:
            nextaxis = axis-1
        else:
            nextaxis = 0

        for ii in range(newshape[0]):
            #print 'x.shape>>',x.shape, 'newshape>>',newshape, 'x[ii].shape>>',x[ii].shape, 'ret[ii].shape>>',ret[ii].shape,'ii>>',ii
            ret[ii] = stdr(x[ii], nsigma=nsigma,niter=niter,finite=finite,\
                                 verbose=verbose, axis=nextaxis)
        return ret
