"""

DOCSTRING INCOMPLETE / OUT OF DATE

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
from plot2 import plot_trajectory_core
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
    fig : plt.Figure 
        The canvas figure that the graphs are plotted onto.
    
    """

    # Loading data
    s_table = data_cut (table, sid, season)

    if len(s_table) == 0:
        print "no data here"
        return

    ### Initializing variables. Most of the following are done thrice:
    # Once with no errors (normal),
    # once with small errors (info),
    # once with large errors (warning).

    ## First: normal

    # Use band_cut to get relevant data chunks.

    j_table = band_cut(s_table, 'j', max_flag=0)
    h_table = band_cut(s_table, 'h', max_flag=0)
    k_table = band_cut(s_table, 'k', max_flag=0)
    hmk_table = band_cut(s_table, 'hmk')    
    # a quick hack to make JHK color-color/mag plots make sense:
    # (since J-H data is only ever used when plotted against H-K, 
    # and we need a second K column to use against H-K)
    jmh_table = band_cut(hmk_table, 'jmh')
    k_table2 = band_cut(k_table, 'hmk')


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
    kcol2 = k_table2.KAPERMAG3 # Needed for the K vs H-K plot

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

    ## Second: info

    # Use band_cut to get relevant data chunks.

    j_table_info = band_cut(s_table, 'j', min_flag=1, max_flag=256)
    h_table_info = band_cut(s_table, 'h', min_flag=1, max_flag=256)
    k_table_info = band_cut(s_table, 'k', min_flag=1, max_flag=256)
    hmk_table_info = band_cut(s_table, 'hmk')    
    # a quick hack to make JHK color-color/mag plots make sense:
    # (since J-H data is only ever used when plotted against H-K, 
    # and we need a second K column to use against H-K)
    jmh_table_info = band_cut(hmk_table, 'jmh')
    k_table_info2 = band_cut(k_table_info, 'hmk', min_flag=1, max_flag=256)


    # The Julian date for CE  2000 January  1 00:00:00.0 UT is
    # JD 2451544.500000
    # (i.e. MJD 51544.0)

    # get a date (x-axis) for each plot
    jdate_info = j_table_info.MEANMJDOBS - 51544
    hdate_info = h_table_info.MEANMJDOBS - 51544
    kdate_info = k_table_info.MEANMJDOBS - 51544
    jmhdate_info = jmh_table_info.MEANMJDOBS - 51544
    hmkdate_info = hmk_table_info.MEANMJDOBS - 51544
    
    # get a magnitude (y-axis) for each plot
    jcol_info = j_table_info.JAPERMAG3
    hcol_info = h_table_info.HAPERMAG3
    kcol_info = k_table_info.KAPERMAG3
    jmh =  jmh_table_info.JMHPNT
    hmk =  hmk_table.HMKPNT
    kcol_info2 = k_table_info2.KAPERMAG3 # Needed for the K vs H-K plot

    # get a magnitude error (y-error) for each plot
    jerr_info = j_table_info.JAPERMAG3ERR
    herr_info = h_table_info.HAPERMAG3ERR
    kerr_info = k_table_info.KAPERMAG3ERR
    jmherr_info=jmh_table_info.JMHPNTERR
    hmkerr_info=hmk_table_info.HMKPNTERR

    # get a quality flag for each plot
    jflag = j_table_info.JPPERRBITS
    hflag = h_table_info.HPPERRBITS
    kflag = k_table_info.KPPERRBITS
    jmhflag = jmh_table_info.JPPERRBITS + jmh_table_info.HPPERRBITS
    hmkflag = hmk_table_info.HPPERRBITS + hmk_table_info.KPPERRBITS
    

    
    ## Define the components and parameters of the figure:
    
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

    ## Start plotting. 
    
    fmt_info = 'D'
    fmt_warn = 'k.'

    # Plot J-band:
    if len(jdate) > 0:
        ax_j.errorbar( jdate, jcol, yerr=jerr, fmt='bo', ecolor='k')
    if len(jdate_info) > 0:
        ax_j.errorbar( jdate_info, jcol_info, yerr=jerr_info, 
                       fmt='b'+fmt_info, ms=4, ecolor='k')
    ax_j.invert_yaxis()

    # Plot H-band:
    if len(hdate) > 0:
        ax_h.errorbar( hdate, hcol, yerr=herr, fmt='go', ecolor='k' )
    if len(hdate_info) > 0:
        ax_h.errorbar( hdate_info, hcol_info, yerr=herr_info, 
                       fmt='g'+fmt_info, ms=4, ecolor='k')
    ax_h.invert_yaxis()

    # Plot K-band:
    if len(kdate) > 0:
        ax_k.errorbar( kdate, kcol, yerr=kerr, fmt='ro', ecolor='k' )
    if len(kdate_info) > 0:
        ax_k.errorbar( kdate_info, kcol_info, yerr=kerr_info, 
                       fmt='r'+fmt_info, ms=4, ecolor='k')
    ax_k.invert_yaxis()

    # Plot J-H vs H-K
    # Note: we use `jmhdate` because of how jmh_table was defined.
    plot_trajectory_core( ax_jhk, hmk, jmh, jmhdate )

    # Plot K vs H-K
    # Note: not sure if this is going to break or not.
    try:
        plot_trajectory_core( ax_khk, hmk, kcol2, hmkdate , ms=False, ctts=False) # gonna update this so that it properly uses K-band main sequence line?
    except:
        pass
        ax_khk.invert_yaxis()

    # Hide the bad labels...
    plt.setp(ax_j.get_xticklabels(), visible=False)
    plt.setp(ax_h.get_xticklabels(), visible=False)

    # Label stuff
    ax_k.set_xlabel( "Time (JD since 01/01/2000)" )

    ax_j.set_ylabel( "J",{'rotation':'horizontal', 'fontsize':'large'} )
    ax_h.set_ylabel( "H",{'rotation':'horizontal', 'fontsize':'large'} )
    ax_k.set_ylabel( "K",{'rotation':'horizontal', 'fontsize':'large'} )

    ax_jhk.set_xlabel( "H-K" )
    ax_jhk.set_ylabel( "J-H")#, {'rotation':'horizontal'})
    ax_khk.set_xlabel( "H-K" )
    ax_khk.set_ylabel( "K")#, {'rotation':'horizontal'})

    return fig
