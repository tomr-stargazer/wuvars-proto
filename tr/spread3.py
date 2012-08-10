""" 
spread3.py

This is the primary statistics-making package within wuvars, 
for "generation 3" of data-handling (i.e. deals with missing 
or error-flagged data sophisticatedly).

Make a spreadsheet of lots of data. 

The purpose of this module is to calculate global statistics for 
the time-series data of many stars, such as mean magnitudes, RMS variability,
the Stetson variability index, (optionally) best-fit periods, etc,
and to bundle all of these values up into a spreadsheet. It's really useful!

Useful functions:
  spreadsheet_write - 
  

Helper functions:
  statcruncher - calculates a large number of stats for one star
                 (formerly arraystat_2)
  

"""

from __future__ import division

import numpy as np

import atpy

import stetson
import stetson_graded
import robust as rb
from helpers3 import data_cut, band_cut
from scargle import fasper as lsp
from timing import lsp_mask
from chi2 import test_analyze
from network2 import get_chip


def reduced_chisq ( m, sigma_m ):
    """ Calculates the reduced chi-squared.

    Inputs:
      m -- an array of photometric magnitudes
      sigma_m -- a corresponding array of photometric uncertainties
      """

    if not (m.size == sigma_m.size):
        raise Exception("Array dimensions mismatch")

    n = m.size
    nu = n - 1

    if n < 2:
        return 0

    return (1./nu) * np.sum( (m - m.mean())**2 / sigma_m**2 )


def Stetson_machine ( s_table, flags=0) :
    """
    Computes the Stetson index on the best combination of bands.

    There's a lot of internal logic here on how to exactly accomplish
    that, and especially which version of the Stetson index to even use.
    
    Parameters
    ----------
    s_table : atpy.Table
        Table with time-series photometry of one star
    flags : int, optional 
        Maximum ppErrBit quality flags to use (default 0)    
    
    Returns
    -------
    Stetson : float
        The Stetson variability index (either "I" or "J" depending 
        on whether 2 or 3 bands were used).
    choice : str {'jhk', 'hk', 'jh', 'jk'}
        Which combination of bands is optimal.
    stetson_nights : int 
        How many nights have all of the optimal combination
        (and therefore, how many nights' worth of data is going into
        the Stetson calculation)

    """
    
    # First, slice the data to find how many nights have
    # a given combination of bands.

    j_table = band_cut( s_table, 'j', max_flag=flags)
    h_table = band_cut( s_table, 'h', max_flag=flags)
    k_table = band_cut( s_table, 'k', max_flag=flags)
    
    jh_table = band_cut( j_table, 'h', max_flag=flags)
    hk_table = band_cut( h_table, 'k', max_flag=flags)
    jk_table = band_cut( j_table, 'k', max_flag=flags)

    jhk_table = band_cut( jh_table, 'k', max_flag=flags)

    # Then we'll measure how many nights are in each combination.

    jh_len = len(jh_table)
    hk_len = len(hk_table)
    jk_len = len(jk_table)
    jhk_len = len(jhk_table)

    # The combination with the most nights (weighted by value^{1})
    # will win. Ties are determined in order: JHK, HK, JH, JK

    max_len = max(jh_len, hk_len, jk_len, jhk_len*2)

    # Now note the winning choice, and compute the relevant index.

    # If there are no simultaneous observations, choose the most-observed band
    # and do a singleband 'Stetson'.

    if max_len == 0:
        j_len = len(j_table)
        h_len = len(h_table)
        k_len = len(k_table)
        max_len_single = max(j_len, h_len, k_len)
        
        if k_len == max_len_single:
            choice = 'k'
            vcol = k_table.KAPERMAG3
            verr = k_table.KAPERMAG3ERR
        elif h_len == max_len_single:
            choice = 'h'
            vcol = h_table.HAPERMAG3
            verr = h_table.HAPERMAG3ERR
        else:
            choice = 'j'
            vcol = j_table.JAPERMAG3
            verr = j_table.JAPERMAG3ERR

        Stetson = stetson.S_singleton(vcol, verr)
        stetson_nights = max_len_single

    elif 2*jhk_len == max_len:
        choice = 'jhk'
        
        jcol = jhk_table.JAPERMAG3; jerr = jhk_table.JAPERMAG3ERR
        hcol = jhk_table.HAPERMAG3; herr = jhk_table.HAPERMAG3ERR
        kcol = jhk_table.KAPERMAG3; kerr = jhk_table.KAPERMAG3ERR

        Stetson = stetson.S(jcol, jerr, hcol, herr, kcol, kerr)

        stetson_nights = jhk_len 
    else:
        if hk_len == max_len:
            choice = 'hk'
            
            bcol = hk_table.HAPERMAG3; berr = hk_table.HAPERMAG3ERR
            vcol = hk_table.KAPERMAG3; verr = hk_table.KAPERMAG3ERR

        elif jh_len == max_len:
            choice = 'jh'

            bcol = jh_table.JAPERMAG3; berr = jh_table.JAPERMAG3ERR
            vcol = jh_table.HAPERMAG3; verr = jh_table.HAPERMAG3ERR

        elif jk_len == max_len:
            choice = 'jk'

            bcol = jk_table.JAPERMAG3; berr = jk_table.JAPERMAG3ERR
            vcol = jk_table.KAPERMAG3; verr = jk_table.KAPERMAG3ERR

        Stetson = stetson.I(bcol, berr, vcol, verr)

        stetson_nights = max_len

    # Finally, return S, the band choice, and how many nights 
    # are going into the Stetson calculation for that choice.

    return (Stetson, choice, stetson_nights)


    # {1}. From Stetson, P. 1996PASP..108..851S, p. 853
    # "... if more than two observations are all obtained 
    # within a time span << the shortest periodicity being 
    # sought, individual observations may be included in 
    # more than one pair. For instance, three closely spaced 
    # observations 'abc' may be distributed into the three 
    # pairs 'ab', 'bc' and 'ac', with two-thirds weight 
    # being assigned to each pair so that in aggregate they 
    # contribute double weight (as a variance determined 
    # from three observations carries twice the weight of a 
    # variance estimated from only two)."  


# I'm going to software hell for copying the following from above...
def graded_Stetson_machine ( s_table, flags=0, min_grade=0.8) :
    """
    Computes the graded Stetson index on the best combination of bands.

    There's a lot of internal logic here on how to exactly accomplish
    that, and especially which version of the Stetson index to even use.
    
    Differs from that "other" Stetson_machine() in that this one uses
    stetson_graded.py.
    
    Parameters
    ----------
    s_table : atpy.Table
        Table with time-series photometry of one star
    flags : int, optional 
        Maximum ppErrBit quality flags to use (default 0)    
    
    Returns
    -------
    Stetson : float
        The Stetson variability index (either "I" or "J" depending 
        on whether 2 or 3 bands were used).
    choice : str {'jhk', 'hk', 'jh', 'jk'}
        Which combination of bands is optimal.
    stetson_nights : int 
        How many nights have all of the optimal combination
        (and therefore, how many nights' worth of data is going into
        the Stetson calculation)
    min_grade : float, optional
        The lowest grade to accept (its weight scales to zero, 
        as do any grades lower than it).

    """
    
    # First, slice the data to find how many nights have
    # a given combination of bands.

    j_table = band_cut( s_table, 'j', max_flag=flags)
    h_table = band_cut( s_table, 'h', max_flag=flags)
    k_table = band_cut( s_table, 'k', max_flag=flags)
    
    jh_table = band_cut( j_table, 'h', max_flag=flags)
    hk_table = band_cut( h_table, 'k', max_flag=flags)
    jk_table = band_cut( j_table, 'k', max_flag=flags)

    jhk_table = band_cut( jh_table, 'k', max_flag=flags)

    # Then we'll measure how many nights are in each combination.

    jh_len = len(jh_table)
    hk_len = len(hk_table)
    jk_len = len(jk_table)
    jhk_len = len(jhk_table)

    # The combination with the most nights (weighted by value^{1})
    # will win. Ties are determined in order: JHK, HK, JH, JK

    max_len = max(jh_len, hk_len, jk_len, jhk_len*2)

    # Now note the winning choice, and compute the relevant index.

    # If there are no simultaneous observations, choose the most-observed band
    # and do a singleband 'Stetson'.

    if max_len == 0:
        j_len = len(j_table)
        h_len = len(h_table)
        k_len = len(k_table)
        max_len_single = max(j_len, h_len, k_len)
        
        if k_len == max_len_single:
            choice = 'k'
            vcol = k_table.KAPERMAG3
            verr = k_table.KAPERMAG3ERR
            vgrade = k_table.KGRADE
        elif h_len == max_len_single:
            choice = 'h'
            vcol = h_table.HAPERMAG3
            verr = h_table.HAPERMAG3ERR
            vgrade = h_table.HGRADE
        else:
            choice = 'j'
            vcol = j_table.JAPERMAG3
            verr = j_table.JAPERMAG3ERR
            vgrade = j_table.JGRADE

        Stetson = stetson_graded.S_singleton(vcol, verr, vgrade)
        stetson_nights = max_len_single

    elif 2*jhk_len == max_len:
        choice = 'jhk'
        
        jcol = jhk_table.JAPERMAG3; jerr = jhk_table.JAPERMAG3ERR
        hcol = jhk_table.HAPERMAG3; herr = jhk_table.HAPERMAG3ERR
        kcol = jhk_table.KAPERMAG3; kerr = jhk_table.KAPERMAG3ERR
        
        jgrade = jhk_table.JGRADE
        hgrade = jhk_table.HGRADE 
        kgrade = jhk_table.KGRADE

        Stetson = stetson_graded.S(jcol, jerr, jgrade, 
                                   hcol, herr, hgrade,
                                   kcol, kerr, kgrade)

        stetson_nights = jhk_len 
    else:
        if hk_len == max_len:
            choice = 'hk'
            
            bcol = hk_table.HAPERMAG3; berr = hk_table.HAPERMAG3ERR
            vcol = hk_table.KAPERMAG3; verr = hk_table.KAPERMAG3ERR
            bgrade = hk_table.HGRADE ; vgrade = hk_table.KGRADE

        elif jh_len == max_len:
            choice = 'jh'

            bcol = jh_table.JAPERMAG3; berr = jh_table.JAPERMAG3ERR
            vcol = jh_table.HAPERMAG3; verr = jh_table.HAPERMAG3ERR
            bgrade = jh_table.JGRADE ; vgrade = jh_table.HGRADE

        elif jk_len == max_len:
            choice = 'jk'

            bcol = jk_table.JAPERMAG3; berr = jk_table.JAPERMAG3ERR
            vcol = jk_table.KAPERMAG3; verr = jk_table.KAPERMAG3ERR
            bgrade = jk_table.JGRADE ; vgrade = jk_table.KGRADE

        Stetson = stetson_graded.I(bcol, berr, bgrade, 
                                   vcol, verr, vgrade)

        stetson_nights = max_len

    # Finally, return S, the band choice, and how many nights 
    # are going into the Stetson calculation for that choice.

    return (Stetson, choice, stetson_nights)


    # {1}. From Stetson, P. 1996PASP..108..851S, p. 853
    # "... if more than two observations are all obtained 
    # within a time span << the shortest periodicity being 
    # sought, individual observations may be included in 
    # more than one pair. For instance, three closely spaced 
    # observations 'abc' may be distributed into the three 
    # pairs 'ab', 'bc' and 'ac', with two-thirds weight 
    # being assigned to each pair so that in aggregate they 
    # contribute double weight (as a variance determined 
    # from three observations carries twice the weight of a 
    # variance estimated from only two)."  


def statcruncher (table, sid, season=0, rob=True, per=True, graded=False, 
                  flags=0) :
    """ 
    Calculates several statistical properties for a given star.

    Will work with "lonely" datapoints (i.e. not all JHK mags are 
    well-defined). Optionally works with graded data, too!

    Parameters
    ----------
    table : atpy.Table
        Table with time-series photometry
    sid : int
        13-digit WFCAM source ID of star to plot
    season : int, optional
        Which observing season of our dataset (1, 2, 3, or all).
        Any value that is not the integers (1, 2, or 3) will be 
        treated as "no season", and no time-cut will be made.
        Note that this is the default behavior.
    rob : bool, optional 
        Use robust statistics, in addition to normal ones?
        (takes longer, default True)
    per : bool, optional 
        Run period-finding? Uses fast chi-squared and lomb-scargle.
        (takes longer, default True)
    graded : bool, optional
        Also calculate Stetson indices using quality grades as weights?
        Uses stetson_graded; requires that the data has been graded by
        night_cleanser.null_cleanser_grader().
    flags : int, optional 
        Maximum ppErrBit quality flags to use (default 0)

    Returns
    -------
    ret : data structure 
        Contains the computed values.
        They can be accessed as attributes 
        (e.g., "ret.j_mean" or "ret.Stetson").

    """
    
    s_table = data_cut ( table, sid, season=season)

    if len(s_table) < 1:
        print "no data for %d!" % sid
        return None

    # First, let's compute single-band statistics. This will require
    # separate data_cuts on each band.

    j_table = band_cut(s_table, 'j', max_flag=flags)
    h_table = band_cut(s_table, 'h', max_flag=flags)
    k_table = band_cut(s_table, 'k', max_flag=flags)

    jmh_table = band_cut(j_table, 'h', max_flag=flags)
    hmk_table = band_cut(h_table, 'k', max_flag=flags)

    # get a date (x-axis) for each 
    jdate = j_table.MEANMJDOBS
    hdate = h_table.MEANMJDOBS
    kdate = k_table.MEANMJDOBS
    jmhdate = jmh_table.MEANMJDOBS
    hmkdate = hmk_table.MEANMJDOBS
#    date = s_table.MEANMJDOBS 
    
    # get a magnitude and magnitude error for each band
    jcol = j_table.JAPERMAG3; jerr = j_table.JAPERMAG3ERR
    hcol = h_table.HAPERMAG3; herr = h_table.HAPERMAG3ERR
    kcol = k_table.KAPERMAG3; kerr = k_table.KAPERMAG3ERR
    jmhcol= jmh_table.JMHPNT; jmherr = jmh_table.JMHPNTERR
    hmkcol= hmk_table.HMKPNT; hmkerr = hmk_table.HMKPNTERR

    # get the RA and DEC columns, checking for sensible values
    racol= s_table.RA[(s_table.RA > 0) & (s_table.RA < 7)]
    decol= s_table.DEC[(s_table.DEC > -4) & (s_table.DEC < 4)]

    # Now let's get some ability to track errorful data.
    # messy_table_j = band_cut( s_table, 'j')
    # messy_table_h = band_cut( s_table, 'h')
    # messy_table_k = band_cut( s_table, 'k')
    # jppcol = messy_table_j.JPPERRBITS
    # hppcol = messy_table_h.HPPERRBITS
    # kppcol = messy_table_k.KPPERRBITS

    # make an empty data structure and just assign it information, then return 
    # the object itself! then there's no more worrying about indices.
    class Empty():
        pass

    ret = Empty()
    
    # How many nights have observations in each band?
    ret.N_j = len(j_table)
    ret.N_h = len(h_table)
    ret.N_k = len(k_table)

    # Mean position of this source
    ret.RA = racol.mean()
    ret.DEC = decol.mean()
    
    # Calculate the Stetson index...
    S, choice, stetson_nights = Stetson_machine (s_table, flags)
    
    ret.Stetson = S
    ret.Stetson_choice = choice
    ret.Stetson_N = stetson_nights

    if graded:
        # Calculate the graded Stetson index...
        g_S, g_choice, g_stetson_nights = Stetson_machine (s_table, flags)
    
        ret.graded_Stetson = g_S
        ret.graded_Stetson_choice = g_choice
        ret.graded_Stetson_N = g_stetson_nights


    # Create parallel data structures for each band, so we can iterate
    ret.j = Empty(); ret.j.data = jcol; ret.j.err = jerr; ret.j.date = jdate   
    ret.h = Empty(); ret.h.data = hcol; ret.h.err = herr; ret.h.date = hdate
    ret.k = Empty(); ret.k.data = kcol; ret.k.err = kerr; ret.k.date = kdate
    ret.jmh = Empty(); ret.jmh.data=jmhcol; ret.jmh.err = jmherr 
    ret.hmk = Empty(); ret.hmk.data=hmkcol; ret.hmk.err = hmkerr
    ret.jmh.date = jmhdate; ret.hmk.date = hmkdate

    ret.j.N = ret.N_j ; ret.h.N = ret.N_h ; ret.k.N = ret.N_k
    ret.jmh.N = len(jmh_table) ; ret.hmk.N = len(hmk_table)

    bands = [ ret.j, ret.h, ret.k, ret.jmh, ret.hmk ]

    for b in bands:
        # use b.data, b.err
        
        # if this band is empty, don't try to do the following assignments
        if b.N == 0: continue

        b.rchi2 = reduced_chisq( b.data, b.err )

        b.mean = b.data.mean()
        b.rms = b.data.std()
        b.min = b.data.min()
        b.max = b.data.max()
        b.range = b.max - b.min

        b.mean_err = b.err.mean()

        # Robust quantifiers simply have an "r" at the end of their names
        if rob:
            b.datar = rb.removeoutliers(b.data, 3, niter=2)
            
            b.meanr = rb.meanr(b.data)
            b.rmsr = rb.stdr(b.data)
            b.minr = b.datar.min()
            b.maxr = b.datar.max()
            b.ranger = b.maxr - b.minr

        # Period finding... is a little dodgy still, and might take forever
        if per==True and b.N > 2:

            b.lsp = lsp(b.date, b.data, 6., 6.) # apologies if this is cluttered
            Jmax = lsp_mask(b.lsp[0], b.lsp[1])
            b.lsp_per = 1./ b.lsp[0][Jmax]
            b.lsp_pow = b.lsp[1][Jmax]
            b.fx2_per = 1./ test_analyze( b.date, b.data, b.err )

    
    # and the pp_max, using the messy table
    # (slated for a re-implementation)
    # ret.jpp_max = jppcol.max()
    # ret.hpp_max = hppcol.max()
    # ret.kpp_max = kppcol.max()

    return ret

def make_sidset ( table ) :
    """ Returns an array of unique source IDs."""

    return np.array(list(set(table.SOURCEID)))

    
def base_lookup (table):
    """ makes a really basic lookup table if you're too lazy"""

    sidarr = np.array( list( set( table.SOURCEID ) ) )
    names = sidarr - 44027700000000L

    Lookup = atpy.Table()
    Lookup.add_column("SOURCEID", sidarr)
    Lookup.add_column("Designation", names )

    return Lookup


def spreadsheet_write (table, lookup, season, outfile, flags=0,
                       Test=False, rob=False, per=False, graded=False):
    """ 
    Makes my spreadsheet! Basically with a big forloop.

    Expect this guy to take a long time.

    Parameters
    ----------
    table : atpy.Table
        Table with time-series photometry
    lookup : atpy.Table
        Table of interesting sources and their names
        (must contain columns "SOURCEID" and "Designation")
    season : int
        Which observing season of our dataset (1, 2, 3, or all).
        Any value that is not the integers (1, 2, or 3) will be 
        treated as "no season", and no time-cut will be made.
        Note that this is the default behavior.
    outfile : str
        What filename to save spreadsheet to.
    flags : int, optional 
        Maximum ppErrBit quality flags to use (default 0)
    Test : bool, optional
        Whether to exit after 30 sources and save to a dummy location
        (default False). Useful for performance/sanity testing.
    rob : bool, optional 
        Use robust statistics, in addition to normal ones?
        (takes longer, default False)
    per : bool, optional 
        Run period-finding? Uses fast chi-squared and lomb-scargle.
        (takes longer, default False)
    graded : bool, optional
        Also calculate Stetson indices using quality grades as weights?
        Uses stetson_graded; requires that the data has been graded by
        night_cleanser.null_cleanser_grader().

      
    Returns
    -------
    None (but writes an output file to `outfile`)

    Note: possibly crashes if given any stars with only 
    1 observation and per=True.
    
    """

    null=np.double(-9.99999488e+08)

    sidarr = lookup.SOURCEID
    names = lookup.Designation
    l = sidarr.size


    Output = atpy.Table()

    N_j = np.ones_like(sidarr)
    N_h = np.copy(N_j)
    N_k = np.copy(N_j)
    RA = np.ones(l)
    DEC= np.ones(l)
#    chip = np.ones_like(sidarr)
#    one_chip = np.ones_like(sidarr==sidarr)
    
    Stetson = np.ones(l)
    Stetson_choice = np.zeros(l, dtype='|S4')
    Stetson_N = np.ones(l, dtype='int')
    
    if graded:
        graded_Stetson = np.ones(l)
        graded_Stetson_choice = np.zeros(l, dtype='|S4')
        graded_Stetson_N = np.ones(l, dtype='int')


    class Band:
        pass

    j = Band()
    h = Band()
    k = Band()
    jmh = Band()
    hmk = Band()

    bands = [ j, h, k, jmh, hmk ]
    band_names = ['j', 'h', 'k', 'jmh', 'hmk']
   
    for b in bands:
        b.mean = np.ones(l) * null
        b.rms =  np.ones(l) * null
        b.rchi2 =np.ones(l) * null
        b.min =  np.ones(l) * null
        b.max =  np.ones(l) * null
        b.range =  np.ones(l) * null
 
        # rob:
        b.meanr =  np.ones(l) * null
        b.rmsr =  np.ones(l) * null
        b.minr =  np.ones(l) * null
        b.maxr =  np.ones(l) * null
        b.ranger =  np.ones(l) * null
        
        # per:
        b.lsp_per =  np.ones(l) * null
        b.lsp_pow =  np.ones(l) * null
        b.fx2_per =  np.ones(l) * null
    
#    color_slope =  np.ones(l)
    # jpp_max = np.ones_like(sidarr)
    # hpp_max = np.ones_like(sidarr)
    # kpp_max = np.ones_like(sidarr)
        

    for sid, i in zip(sidarr, range(len(sidarr)) ):

        # v for values
        v = statcruncher (table, sid, season, rob, per, flags=flags)
        if v == None:
            #skip assigning anything!
            continue
        vbands = [v.j, v.h, v.k, v.jmh, v.hmk]
        vband_names = ['v.j', 'v.h', 'v.k', 'v.jmh', 'v.hmk']

        N_j[i] = v.N_j
        N_h[i] = v.N_h
        N_k[i] = v.N_k

        RA[i] = v.RA
        DEC[i] = v.DEC
#        chip[i] = v.chip
#        one_chip[i] = v.one_chip
        Stetson[i] = v.Stetson
        Stetson_choice[i] = v.Stetson_choice
        Stetson_N[i] = v.Stetson_N

        if graded:
            graded_Stetson[i] = v.graded_Stetson
            graded_Stetson_choice[i] = v.graded_Stetson_choice
            graded_Stetson_N[i] = v.graded_Stetson_N


        for b, bn, vb, vbn in zip(bands, band_names, vbands, vband_names):
            
#            print "writing into band " + repr(b) + " with band " + str(vb)
#            print "writing into band %s with band %s" % (bn, vbn)
            
            if vb.N == 0: continue

            b.mean[i] = vb.mean
#            print "mean: " + str(vb.mean)
#            print "recorded mean: " + str(b.mean[i])
#            return
            b.rms[i] = vb.rms
            b.min[i] = vb.min
            b.max[i] = vb.max
            b.range[i] = vb.range

            b.rchi2[i] = vb.rchi2

#            b.mean_err[i] = vb.mean_err
            
            if rob:
                b.meanr[i] = vb.meanr
                b.rmsr[i] = vb.rmsr
                b.maxr[i] = vb.maxr
                b.minr[i] = vb.minr
                b.ranger[i] = vb.ranger

            if per and vb.N > 2:
                b.lsp_per[i] = vb.lsp_per
                b.lsp_pow[i] = vb.lsp_pow
                b.fx2_per[i] = vb.fx2_per
    
#        color_slope[i] = v.color_slope
        # jpp_max[i] = v.jpp_max
        # hpp_max[i] = v.hpp_max
        # kpp_max[i] = v.kpp_max
            

        if (i > 30 and Test): 
            print i
            print "End of test"
            break

    Output.add_column('SOURCEID',sidarr)
    Output.add_column('Name',names)
    Output.add_column('RA', RA, unit='RADIANS')
    Output.add_column('DEC', DEC, unit='RADIANS')
    Output.add_column('N_j', N_j)
    Output.add_column('N_h', N_h)
    Output.add_column('N_k', N_k)

#    Output.add_column('chip', chip)
#    Output.add_column('one_chip', one_chip)

    Output.add_column('Stetson', Stetson)
    Output.add_column('Stetson_choice', Stetson_choice)
    Output.add_column('Stetson_N', Stetson_N)
    
    if graded:
        Output.add_column('graded_Stetson', graded_Stetson)
        Output.add_column('graded_Stetson_choice', graded_Stetson_choice)
        Output.add_column('graded_Stetson_N', graded_Stetson_N)


    for b, band_name in zip(bands, band_names):
        bn = band_name + '_'
        Output.add_column(bn+'mean', b.mean)
#        print Output[bn+'mean'][0]
        Output.add_column(bn+'rms', b.rms)
        Output.add_column(bn+'min', b.min)
        Output.add_column(bn+'max', b.max)
        Output.add_column(bn+'range', b.range)

        Output.add_column(bn+'rchi2', b.rchi2)

        if rob:
            Output.add_column(bn+'meanr', b.meanr)
            Output.add_column(bn+'rmsr', b.rmsr)
            Output.add_column(bn+'minr', b.minr)
            Output.add_column(bn+'maxr', b.maxr)
            Output.add_column(bn+'ranger', b.ranger)

        if per:
            Output.add_column(bn+'lsp_per', b.lsp_per)
            Output.add_column(bn+'lsp_pow', b.lsp_pow)
            Output.add_column(bn+'fx2_per', b.fx2_per)

            
#    Output.add_column('color_slope', color_slope)
    # Output.add_column('jpp_max', jpp_max)
    # Output.add_column('hpp_max', hpp_max)
    # Output.add_column('kpp_max', kpp_max)

    Output.write(outfile, overwrite=True)
    print "Wrote output to %s" % outfile

    return

def spread_write_test (table, lookup, flags=0) :
    """ Tests spreadsheet_write."""

    test_path = '/home/trice/reu/DATA/Merged_Catalogs/spreadsheet/test.fits'
    spreadsheet_write (table, lookup, -1, test_path, flags=flags, Test=True)
