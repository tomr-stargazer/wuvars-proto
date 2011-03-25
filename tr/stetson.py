''' An implementation of Peter Stetson's variability index; see
'On the Automatic Determination of Light-Curve Parameters for Cepheid
Variables', Stetson 1996, PASP.

'''

import numpy as np


def I(b, sigma_b, v, sigma_v):
    '''The Welch/Stetson variability index I.
    Calculates variability index given two magnitude bands and uncertainties.

    Inputs: two arrays of magnitudes and two arrays of associated sigmas.
    Returns a number.
    '''

    if not (b.size == sigma_b.size == v.size == sigma_v.size):
        raise Exception("Array dimensions mismatch")
        return

    n = b.size

    s = np.sum( (b - b.mean())/sigma_b * (v - v.mean())/sigma_v )

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

    

def S (j, sigma_j, h, sigma_h, k, sigma_k) :
    """
    Computes the Stetson variability index for one star that has
    3 observations on each night. Uses Carpenter et al.'s notation.
    Simplified from the general expression assuming 3 observations every
      night.
    
    INPUTS:
        j: an array of J-band magnitudes
        sigma_j: an array of corresponding J-band uncertainties
        h: an array of H-band magnitudes
        sigma_h: an array of corresponding H-band uncertainties
        k: an array of K-band magnitudes
        sigma_k: an array of corresponding K-band uncertainties

    OUTPUTS:
        s: the Stetson variability index for 3 bands

    """

    n = j.size
    
    # Perhaps hackish
    if n < 2:
        return 0

    d_j = delta(j, sigma_j, j.mean(), n)
    d_h = delta(h, sigma_h, h.mean(), n)
    d_k = delta(k, sigma_k, k.mean(), n)

    P_i = np.array( [d_j * d_h,
                     d_h * d_k,
                     d_j * d_k] )

    # I originally had two sums going: one over P_i, and one over all the 
    # elements of n, but then I realized that a single sum over all axes
    # did the exact same thing (I tested it) so now it's back to one sum.
    s = np.sum( np.sign( P_i ) * np.sqrt( np.abs( P_i ))) /(n*1.)

    return s

def S_sid (table, sid, season=123) :
    """ Calculates the Stetson J index for a given source.

    Inputs:
      table -- an atpy table with time-series photometry
      sid -- a Source ID from WFCAM (13 digits)
      season -- Which observing season of our dataset (1,2, 3, or all)
    """
    from tr_helpers import season_cut

    s_table = season_cut(table, sid, season, flags=0)
    jcol = s_table.JAPERMAG3
    hcol = s_table.HAPERMAG3
    kcol = s_table.KAPERMAG3
    jerr = s_table.JAPERMAG3ERR
    herr = s_table.HAPERMAG3ERR
    kerr = s_table.KAPERMAG3ERR


    return S (jcol, jerr, hcol, herr, kcol, kerr)
