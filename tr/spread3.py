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
  spread_write - 
  

Helper functions:
  statcruncher - calculates a large number of stats for one star
                 (formerly arraystat_2)
  

"""

import atpy
import numpy as np
import stetson
import robust as rb
from helpers3 import data_cut
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

    if 2*jhk_len == max_len:
        choice = 'jhk'
        
        jcol = jhk_table.JAPERMAG3; jerr = jhk_table.JAPERMAG3ERR
        hcol = jhk_table.HAPERMAG3; herr = jhk_table.HAPERMAG3ERR
        kcol = jhk_table.KAPERMAG3; kerr = jhk_table.KAPERMAG3ERR

        Stetson = stetson.S(jcol, jerr, hcol, herr, kcol, kerr)
 
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


def statcruncher (table, sid, season=0, rob=True, per=True, flags=0) :
    """ Calculates several statistical properties for a given star.

    Will work with "lonely" datapoints (i.e. not all JHK mags are 
    well-defined). 

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
    jmh_date = jmh_table.MEANMJDOBS
    hmk_date = hmk_table.MEANMJDOBS
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
    messy_table_j = band_cut( s_table, 'j')
    messy_table_h = band_cut( s_table, 'h')
    messy_table_k = band_cut( s_table, 'k')
    jppcol = messy_table_j.JPPERRBITS
    hppcol = messy_table_h.HPPERRBITS
    kppcol = messy_table_k.KPPERRBITS

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


    # Create parallel data structures for each band, so we can iterate
    ret.j = Empty();   ret.j.data = jcol;   ret.j.err = jerr
    ret.h = Empty();   ret.h.data = hcol;   ret.h.err = herr
    ret.k = Empty();   ret.k.data = kcol;   ret.k.err = kerr
    ret.jmh = Empty(); ret.jmh.data=jmhcol; ret.jmh.err = jmherr
    ret.hmk = Empty(); ret.hmk.data=hmkcol; ret.hmk.err = hmkerr


    bands = [ ret.j, ret.h, ret.k, ret.jmh, ret.hmk ]

    for b in bands:
        # use b.data, b.err
        
        b.rchi2 = reduced_chisq( b.data, b.err )

        b.mean = b.data.mean()
        b.rms = b.data.std()
        b.min = b.data.min()
        b.max = b.data.max()
        b.peak_trough = b.max - b.min

        b.mean_err = b.err.mean()

        # Robust quantifiers simply have an "r" at the end of their names
        if rob:
            b.datar = rb.removeoutliers(b.data, 3, niter=2)
            
            b.meanr = rb.meanr(b.data)
            b.rmsr = rb.stdr(b.data)
            b.minr = b.datar.min()
            b.maxr = b.datar.max()
            b.peak_troughr = b.maxr - b.minr

        # Period finding... is a little dodgy still, and might take forever
        if per:
            
            b.lsp = lsp(date, b.data, 6., 6.) # apologies if this is cluttered
            Jmax = lsp_mask(b.lsp[0], b.lsp[1])
            b.lsp_per = 1./ b.lsp[0][Jmax]
            b.lsp_pow = b.lsp[1][Jmax]
            b.fx2_per = 1./ test_analyze( date, b.data, b.err )

    
    # and the pp_max, using the messy table
    # (slated for a re-implementation)
    ret.jpp_max = jppcol.max()
    ret.hpp_max = hppcol.max()
    ret.kpp_max = kppcol.max()

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
                       Test=False, rob=False, per=False):
    """ 
    Makes my spreadsheet! Basically with a big forloop.

    Inputs:
      table -- an ATpy table with time-series photometry
      lookup -- an ATpy table of interesting sources and their designations
      season -- which season to select (1,2,3, or other=All)
      outfile -- where to save the output table
      
    Optional inputs:
      flags -- Maximum ppErrBit quality flags to use (default 0)
      Test -- whether to exit after 30 sources (default False)
      rob -- also use Robust statistics? (takes longer, default False)
      per -- run period-finding? (takes longer, default False)

    Returns:
      None (but writes an output file to `outfile`)

    Note: crashes if given any stars with only 1 observation and per=True.
    """

    sidarr = lookup.SOURCEID
    names = lookup.Designation
    l = sidarr.size


    Output = atpy.Table()

    N = np.ones_like(sidarr)
    RA = np.ones(l)
    DEC= np.ones(l)
    chip = np.ones_like(sidarr)
    one_chip = np.ones_like(sidarr==sidarr)
    
    Stetson = np.ones(l)
    
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
        b.mean = np.ones(l)
        b.rms =  np.ones(l)
        b.rchi2 =np.ones(l)
        b.min =  np.ones(l)
        b.max =  np.ones(l)
        b.peak_trough =  np.ones(l)
        
        b.meanr =  np.ones(l)
        b.rmsr =  np.ones(l)
        b.minr =  np.ones(l)
        b.maxr =  np.ones(l)
        b.peak_troughr =  np.ones(l)
        
        # per:
        b.lsp_per =  np.ones(l)
        b.lsp_pow =  np.ones(l)
        b.fx2_per =  np.ones(l)
    
    color_slope =  np.ones(l)
    jpp_max = np.ones_like(sidarr)
    hpp_max = np.ones_like(sidarr)
    kpp_max = np.ones_like(sidarr)
        

    for sid, i in zip(sidarr, range(len(sidarr)) ):

        # v for values
        v = arraystat_2 (table, sid, season, rob, per, flags=flags)
        if v == None:
            #skip assigning anything!
            continue
        vbands = [v.j, v.h, v.k, v.jmh, v.hmk]
        vband_names = ['v.j', 'v.h', 'v.k', 'v.jmh', 'v.hmk']

        N[i] = v.N
        RA[i] = v.RA
        DEC[i] = v.DEC
        chip[i] = v.chip
        one_chip[i] = v.one_chip
        Stetson[i] = v.Stetson

        for b, bn, vb, vbn in zip(bands, band_names, vbands, vband_names):
            
#            print "writing into band " + repr(b) + " with band " + str(vb)
#            print "writing into band %s with band %s" % (bn, vbn)

            b.mean[i] = vb.mean
#            print "mean: " + str(vb.mean)
#            print "recorded mean: " + str(b.mean[i])
#            return
            b.rms[i] = vb.rms
            b.min[i] = vb.min
            b.max[i] = vb.max
            b.peak_trough[i] = vb.peak_trough

            b.rchi2[i] = vb.rchi2

#            b.mean_err[i] = vb.mean_err
            
            if rob:
                b.meanr[i] = vb.meanr
                b.rmsr[i] = vb.rmsr
                b.maxr[i] = vb.maxr
                b.minr[i] = vb.minr
                b.peak_troughr[i] = vb.peak_troughr

            if per:
                b.lsp_per[i] = vb.lsp_per
                b.lsp_pow[i] = vb.lsp_pow
                b.fx2_per[i] = vb.fx2_per
    
        color_slope[i] = v.color_slope
        jpp_max[i] = v.jpp_max
        hpp_max[i] = v.hpp_max
        kpp_max[i] = v.kpp_max
            

        if (i > 30 and Test): 
            print i
            print "End of test"
            break

    Output.add_column('SOURCEID',sidarr)
    Output.add_column('Name',names)
    Output.add_column('RA', RA, unit='RADIANS')
    Output.add_column('DEC', DEC, unit='RADIANS')
    Output.add_column('N', N)
    Output.add_column('chip', chip)
    Output.add_column('one_chip', one_chip)

    Output.add_column('Stetson', Stetson)

    for b, band_name in zip(bands, band_names):
        bn = band_name + '_'
        Output.add_column(bn+'mean', b.mean)
#        print Output[bn+'mean'][0]
        Output.add_column(bn+'rms', b.rms)
        Output.add_column(bn+'min', b.min)
        Output.add_column(bn+'max', b.max)
        Output.add_column(bn+'peak_trough', b.peak_trough)

        Output.add_column(bn+'rchi2', b.rchi2)

        if rob:
            Output.add_column(bn+'meanr', b.meanr)
            Output.add_column(bn+'rmsr', b.rmsr)
            Output.add_column(bn+'minr', b.minr)
            Output.add_column(bn+'maxr', b.maxr)
            Output.add_column(bn+'peak_troughr', b.peak_troughr)

        if per:
            Output.add_column(bn+'lsp_per', b.lsp_per)
            Output.add_column(bn+'lsp_pow', b.lsp_pow)
            Output.add_column(bn+'fx2_per', b.fx2_per)

            
    Output.add_column('color_slope', color_slope)
    Output.add_column('jpp_max', jpp_max)
    Output.add_column('hpp_max', hpp_max)
    Output.add_column('kpp_max', kpp_max)

    Output.write(outfile, overwrite=True)
    print "Wrote output to %s" % outfile

    return

def spread_write_test (table, lookup) :
    """ Tests spreadsheet_write."""

    test_path = '/home/trice/reu/DATA/Merged_Catalogs/spreadsheet/test.fits'
    spreadsheet_write (table, lookup, -1, test_path, Test=True)
