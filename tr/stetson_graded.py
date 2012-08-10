""" An implementation of Peter Stetson's variability index; see
'On the Automatic Determination of Light-Curve Parameters for Cepheid
Variables', Stetson 1996, PASP..108..851S

This code differs from stetson.py in that it supports "graded" 
Stetson calculations, where the quality of the relevant nights
is factored into the Stetson calculation.

"""

from __future__ import division
import numpy as np


def I(b, sigma_b, grade_b, v, sigma_v, grade_v, min_grade=0.8):
    """
    The Welch/Stetson variability index I, with quality-assigned weights.
    
    Calculates variability index given two magnitude bands and uncertainties.
    Uses quality grades to assign weights.

    Parameters
    ----------
    b, v : array_like
        Two arrays of magnitude values
    sigma_b, sigma_v : array_like
        Two corresponding arrays of uncertainty values.
    grade_b, grade_v : array_like
        Corresponding arrays of quality grades.
    min_grade : float, optional
        The lowest grade to accept (its weight scales to zero, 
        as do any grades lower than it).

    Returns
    -------
    I : float
        The Stetson variability index "I" over the input arrays.

    Notes
    -----
    This code implements the following mathematical formula:
    
    .. math:: I = \sqrt{ \frac{1}{n(n-1)}}\frac{\sum_{i=1}^{n} w_i \left(\frac{b_i-\bar b}{\sigma_{b,i}} \right ) \left(\frac{v_i-\bar v}{\sigma_{v,i}} \right )}{\sum_{i=1}^{n} w_i}

    which is the Welch/Stetson variability index for two color bands, 
    with the slight modification of using weights :math:`w_k` that 
    are a function of the product of the nightly grade in each band.

    """

    if not (b.size == sigma_b.size == grade_b.size == 
            v.size == sigma_v.size == grade_v.size):
        raise Exception("Array dimensions mismatch")
        return

    n = b.size

    # Perhaps hackish
    if n < 2:
        return 0

    # Ensure that no grades are below min_grade
    grade_b[grade_b < min_grade] = min_grade
    grade_v[grade_v < min_grade] = min_grade

    c = 1./(1-min_grade)
    wk = ((grade_b - min_grade) * c) * ((grade_v - min_grade) * c)

    # Ensure that no weights are below zero
    wk[wk < 0] = 0
    

    s = (np.sum( wk * (b - b.mean())/sigma_b * (v - v.mean())/sigma_v ) /
         np.sum( wk ))

    I = np.sqrt(1. / (n*(n-1))) * s

    return I


def delta (m, sigma_m, mean_m, n):
    """ Normalized residual / "relative error" for one observation. 
    Used in Stetson's J variability index.

    INPUTS:
        m: a single magnitude measurement in a certain band
        sigma_m: the uncertainty on that measurement
        mean_m: the mean magnitude in that band
        n: the number of observations in that band

    OUTPUTS:
        d: the "relative error"
    """
    
    d = np.sqrt( n / (n-1) ) * (m - mean_m) / sigma_m
    
    return d

    

def S (j, sigma_j, grade_j, 
       h, sigma_h, grade_h, 
       k, sigma_k, grade_k) :
    """
    Computes the Stetson variability index for one star that has
    3 observations on each night. Uses Carpenter et al.'s notation.
    Simplified from the general expression assuming 3 observations every
    night.

    Parameters
    ----------
    j, h, k : array-like
        Magnitudes in bands J, H, and K.
    sigma_j, sigma_h, sigma_k : array-like
        Three corresponding arrays of uncertainty values.
    grade_j, grade_h, grade_k : array_like
        Corresponding arrays of quality grades.
    min_grade : float, optional
        The lowest grade to accept (its weight scales to zero, 
        as do any grades lower than it).

    Returns
    -------
    S : float
        The Stetson variability index "J" over the 3 input bands.
    
    Notes
    -----
    This code implements the following mathematical formula:
    
    .. math:: J = \frac{\sum_{k=1}^n w_k \textrm{sgn}(P_k)\sqrt{|P_k|}}{\sum_{k=1}^n w_k}
     
    where the meaning of the symbol :math:`P_k` is defined in 
    the reference (Stetson 1996).
    
    """

    n = j.size
    
    # Perhaps hackish
    if n < 2:
        return 0

    # Calculate deltas
    d_j = delta(j, sigma_j, j.mean(), n)
    d_h = delta(h, sigma_h, h.mean(), n)
    d_k = delta(k, sigma_k, k.mean(), n)

    # Ensure that no grades are below min_grade
    grade_j[grade_j < min_grade] = min_grade
    grade_h[grade_h < min_grade] = min_grade
    grade_k[grade_k < min_grade] = min_grade

    c = 1./(1-min_grade)
    wjh = ((grade_j - min_grade) * c) * ((grade_h - min_grade) * c)
    whk = ((grade_h - min_grade) * c) * ((grade_k - min_grade) * c)
    wjk = ((grade_j - min_grade) * c) * ((grade_k - min_grade) * c)

    wk = np.array( [wjh, whk, wjk] )
    
    # Ensure that no weights are below zero
    wk[wk < 0] = 0

    P_i = np.array( [d_j * d_h,
                     d_h * d_k,
                     d_j * d_k] )

    # I originally had two sums going: one over P_i, and one over all the 
    # elements of n, but then I realized that a single sum over all axes
    # did the exact same thing (I tested it) so now it's back to one sum.
    S = (np.sum( wk * np.sign( P_i ) * np.sqrt( np.abs( P_i ))) /
         np.sum( wk ))

    return S

def S_singleton (v, sigma_v, grade_v, min_grade=0.8):
    """
    Computes the 'Stetson' index for a star that has no simultaneous 
    observations (i.e., only observed in a single color).
    """

    n = v.size
    
    # Perhaps hackish
    if n < 2:
        return 0

    d_v = delta(v, sigma_v, v.mean(), n)

    P_i = np.array( [d_v**2 - 1] )

    # Stolen from S(). The sum may be completely unneccessary?
    s = np.sum( np.sign( P_i ) * np.sqrt( np.abs( P_i ))) /(n*1.)

    return s



def S_sid (table, sid, season=123, flags=0) :
    """ Calculates the Stetson J index for a given source.

    Inputs:
      table -- an atpy table with time-series photometry
      sid -- a Source ID from WFCAM (13 digits)
      season -- Which observing season of our dataset (1,2, 3, or all)
    """
    from tr_helpers import season_cut

    s_table = season_cut(table, sid, season, flags=flags)
    jcol = s_table.JAPERMAG3
    hcol = s_table.HAPERMAG3
    kcol = s_table.KAPERMAG3
    jerr = s_table.JAPERMAG3ERR
    herr = s_table.HAPERMAG3ERR
    kerr = s_table.KAPERMAG3ERR


    return S (jcol, jerr, hcol, herr, kcol, kerr)
