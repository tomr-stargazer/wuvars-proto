"""

NOTHING BELOW THIS IS USEFUL!!! 

This is the primary visualization package within wuvars.

Useful functions:
 # Primary functions
 lc - Plots J, H, K lightcurves
 phase - Plots J, H, K lightcurves folded by some period
 lsp_power - Plots J, H, K periodograms for one star.
 
 # Helper functions
 reddening_vector - Plots a reddening vector in a JHK color-color diagram
 plot_trajectory_vanilla - Plots the main sequence and CTTS locus in JHK
 plot_trajectory_core - Plots the trajectory of some star in JHK space
 plot_trajectory - Takes a source from a table and plots its JHK trajectory
 plot_phase_core - Plots a pretty period-folded lightcurve

"""

import numpy as np
import matplotlib.pyplot as plt

import atpy

from helpers3 import data_cut, band_cut

#import coords
#import stetson
#from chi2 import test_analyze
#from scargle import fasper as lsp
#from timing import lsp_mask
#from tr_helpers import season_cut, data_cut, ensemble_cut



def lc (table, sid, season=0):
    """ 
    Plots J, H, K lightcurves, plus color-color and color-mag
    trajectories, for one star.

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


    Returns
    -------
    None
    
    """

    # Loading data
    s_table = data_cut (table, sid, season)

    if len(s_table) == 0:
        print "no data here"
        return

    # Use band_cut to get relevant thingies
    
    j_table = band_cut(s_table, 'j')
    h_table = band_cut(s_table, 'h')
    k_table = band_cut(s_table, 'k')
    jmh_table = band_cut(s_table, 'jmh')
    hmk_table = band_cut(s_table, 'hmk')
    
    # The Julian date for CE  2000 January  1 00:00:00.0 UT is
    # JD 2451544.500000
    # (i.e. MJD 51544.0)

    # get a date (x-axis) for each plot
    jdate = j_table.MEANMJDOBS - 51544
    hdate = h_table.MEANMJDOBS - 51544
    kdate = k_table.MEANMJDOBS - 51544
    jmhdate = jmh_table.MEANMJDOBS - 51544
    hmkdate = hmk_table.MEANMJDOBS - 51544

    # get a magnitude (y-axis) for each plot
    jcol = j_table.JAPERMAG3
    hcol = h_table.HAPERMAG3
    kcol = k_table.KAPERMAG3
    jmh =  jmh_table.JMHPNT
    hmk =  hmk_table.HMKPNT

    # get a magnitude error (y-error) for each plot
    jerr = j_table.JAPERMAG3ERR
    herr = h_table.HAPERMAG3ERR
    kerr = k_table.KAPERMAG3ERR
    jmherr=jmh_table.JMHPNTERR
    hmkerr=hmk_table.HMKPNTERR

    # get a quality flag for each plot
    jflag = j_table.JPPERRBITS
    hflag = h_table.HPPERRBITS
    kflag = k_table.KPPERRBITS
    jmhflag = jmh_table.JPPERRBITS + jmh_table.HPPERRBITS
    hmkflag = hmk_table.HPPERRBITS + hmk_table.KPPERRBITS
    
    
    ## Define the components of the figure:
    
    fig = plt.figure(figsize = (10, 6), dpi=80,
                     facecolor='w', edgecolor='k')

    bottom = 0.1
    height = .25
    left = 0.075
    width = 0.5

    ax_k = fig.add_axes( (left, bottom, width, height) )
    ax_h = fig.add_axes( (left, bottom+.3, width, height), sharex=ax_k )
    ax_j = fig.add_axes( (left, bottom+.6, width, height), sharex=ax_k )
    
    ax_jhk = fig.add_axes( (.65, bottom, .3, .375) )
    ax_khk = fig.add_axes( (.65, bottom+.475, .3, .375) )

