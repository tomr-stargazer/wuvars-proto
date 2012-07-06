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
from chi2 import test_analyze
from scargle import fasper as lsp
from timing import lsp_mask
from spread3 import Stetson_machine

#import coords
#import stetson

#from tr_helpers import season_cut, data_cut, ensemble_cut



def lc (table, sid, season=0, outfile='', name='', stetson=True, png_too=False):
    """ 
    Plots J, H, K lightcurves, plus color-color and color-mag
    trajectories, for one star.

    Will display "lonely" datapoints (i.e. not all JHK mags are 
    well-defined), and plots error-flagged data as different symbols.
    Compare to plot2.lc() which does neither.

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
    outfile : str, optional
        What filename to save plot to. Default behavior (when 
        `outfile` is an empty string) is to display plot on-screen
        and *not* save to file.
    name : str, optional
        A string to use as a plot header.
    png_too : bool, optional (default: False)
        If `png_too` is True (and `outfile` is not ''), then 
        save the plot in 3 file formats: PDF, PNG, and EPS.
        Do not specify a file extension in `outfile`.

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

    ## Let's do the 3 lightcurve bands.

    ### Initializing variables. Most of the following are done thrice:
    # Once with no errors (normal),
    # once with small errors (info),
    # once with large errors (warning).

    ## First: normal

    # Use band_cut to get relevant data chunks.

    j_table = band_cut(s_table, 'j', max_flag=0)
    h_table = band_cut(s_table, 'h', max_flag=0)
    k_table = band_cut(s_table, 'k', max_flag=0)


    # The Julian date for CE  2000 January  1 00:00:00.0 UT is
    # JD 2451544.500000
    # (i.e. MJD 51544.0)

    # get a date (x-axis) for each plot
    jdate = j_table.MEANMJDOBS - 51544
    hdate = h_table.MEANMJDOBS - 51544
    kdate = k_table.MEANMJDOBS - 51544
    
    # get a magnitude (y-axis) for each plot
    jcol = j_table.JAPERMAG3
    hcol = h_table.HAPERMAG3
    kcol = k_table.KAPERMAG3

    # get a magnitude error (y-error) for each plot
    jerr = j_table.JAPERMAG3ERR
    herr = h_table.HAPERMAG3ERR
    kerr = k_table.KAPERMAG3ERR

    # get a quality flag for each plot
    # (aftercomment: not sure we're going to ever use these)
    jflag = j_table.JPPERRBITS
    hflag = h_table.HPPERRBITS
    kflag = k_table.KPPERRBITS

    ## Second: info

    # Use band_cut to get relevant data chunks.

    j_table_info = band_cut(s_table, 'j', min_flag=1, max_flag=256)
    h_table_info = band_cut(s_table, 'h', min_flag=1, max_flag=256)
    k_table_info = band_cut(s_table, 'k', min_flag=1, max_flag=256)


    # The Julian date for CE  2000 January  1 00:00:00.0 UT is
    # JD 2451544.500000
    # (i.e. MJD 51544.0)

    # get a date (x-axis) for each plot
    jdate_info = j_table_info.MEANMJDOBS - 51544
    hdate_info = h_table_info.MEANMJDOBS - 51544
    kdate_info = k_table_info.MEANMJDOBS - 51544
    
    # get a magnitude (y-axis) for each plot
    jcol_info = j_table_info.JAPERMAG3
    hcol_info = h_table_info.HAPERMAG3
    kcol_info = k_table_info.KAPERMAG3

    # get a magnitude error (y-error) for each plot
    jerr_info = j_table_info.JAPERMAG3ERR
    herr_info = h_table_info.HAPERMAG3ERR
    kerr_info = k_table_info.KAPERMAG3ERR

    # get a quality flag for each plot
    jflag = j_table_info.JPPERRBITS
    hflag = h_table_info.HPPERRBITS
    kflag = k_table_info.KPPERRBITS
    
    ## Third: warning

    # Use band_cut to get relevant data chunks.

    j_table_warn = band_cut(s_table, 'j', min_flag=257)
    h_table_warn = band_cut(s_table, 'h', min_flag=257)
    k_table_warn = band_cut(s_table, 'k', min_flag=257)


    # The Julian date for CE  2000 January  1 00:00:00.0 UT is
    # JD 2451544.500000
    # (i.e. MJD 51544.0)

    # get a date (x-axis) for each plot
    jdate_warn = j_table_warn.MEANMJDOBS - 51544
    hdate_warn = h_table_warn.MEANMJDOBS - 51544
    kdate_warn = k_table_warn.MEANMJDOBS - 51544
    
    # get a magnitude (y-axis) for each plot
    jcol_warn = j_table_warn.JAPERMAG3
    hcol_warn = h_table_warn.HAPERMAG3
    kcol_warn = k_table_warn.KAPERMAG3

    # get a magnitude error (y-error) for each plot
    jerr_warn = j_table_warn.JAPERMAG3ERR
    herr_warn = h_table_warn.HAPERMAG3ERR
    kerr_warn = k_table_warn.KAPERMAG3ERR

    # get a quality flag for each plot
    jflag = j_table_warn.JPPERRBITS
    hflag = h_table_warn.HPPERRBITS
    kflag = k_table_warn.KPPERRBITS

    
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
    fmt_warn = '.'

    # Plot J-band:
    if len(jdate) > 0:
        ax_j.errorbar( jdate, jcol, yerr=jerr, fmt='bo', ecolor='k')
    if len(jdate_info) > 0:
        ax_j.errorbar( jdate_info, jcol_info, yerr=jerr_info, 
                       fmt='b'+fmt_info, ms=4, ecolor='k')
    if len(jdate_warn) > 0:
        ax_j.errorbar( jdate_warn, jcol_warn, yerr=jerr_warn, 
                       fmt='k'+fmt_warn, ecolor='k')
    ax_j.invert_yaxis()

    # Plot H-band:
    if len(hdate) > 0:
        ax_h.errorbar( hdate, hcol, yerr=herr, fmt='go', ecolor='k' )
    if len(hdate_info) > 0:
        ax_h.errorbar( hdate_info, hcol_info, yerr=herr_info, 
                       fmt='g'+fmt_info, ms=4, ecolor='k')
    if len(hdate_warn) > 0:
        ax_h.errorbar( hdate_warn, hcol_warn, yerr=herr_warn, 
                       fmt='k'+fmt_warn, ecolor='k')
    ax_h.invert_yaxis()

    # Plot K-band:
    if len(kdate) > 0:
        ax_k.errorbar( kdate, kcol, yerr=kerr, fmt='ro', ecolor='k' )
    if len(kdate_info) > 0:
        ax_k.errorbar( kdate_info, kcol_info, yerr=kerr_info, 
                       fmt='r'+fmt_info, ms=4, ecolor='k')
    if len(kdate_warn) > 0:
        ax_k.errorbar( kdate_warn, kcol_warn, yerr=kerr_warn, 
                       fmt='k'+fmt_warn, ecolor='k')
    ax_k.invert_yaxis()

    ## Now let's do the 2 color-mag/color-color plots.

    # We'll use different data-cuts for the two different plots.
    # Relevant comment: I made an executive call to include only
    # 'normal' and 'info'-flagged data in the C-C and C-M plots
    # (i.e. max_flag=256 in all relevant bands).
    
    # In the color-mag plot, we need data where H and K are defined 
    # everywhere. That's two cuts.
    khk_table = band_cut( band_cut(s_table, 'k', max_flag=256),
                          'h', max_flag=256)

    khkdate = khk_table.MEANMJDOBS - 51544
    k_khk = khk_table.KAPERMAG3
    hmk_khk = khk_table.HMKPNT

    # In the color-color plot, we need data where J, H, and K are
    # defined everywhere. That's one more cut.
    jhk_table = band_cut(khk_table, 'j', max_flag=256)

    jhkdate = jhk_table.MEANMJDOBS - 51544
    jmh_jhk = jhk_table.JMHPNT
    hmk_jhk = jhk_table.HMKPNT

    # Plot J-H vs H-K using the "jhk_" variables.
    try:
        plot_trajectory_core( ax_jhk, hmk_jhk, jmh_jhk, jhkdate )
    except Exception:
        print "JHK plot broke!"
        pass

    # Plot K vs H-K using the "khk_" variables.
    try:
        plot_trajectory_core( ax_khk, hmk_khk, k_khk, khkdate, 
                              ms=False, ctts=False) 
    except Exception:
        print "KHK plot broke!"
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

    if name != '':
        ax_j.set_title(name)

    if stetson == True:
        S, choice, n = Stetson_machine( s_table, flags=256 )
        stet_string = "S = %.2f" % S
        ax_khk.set_title(stet_string)

    if outfile == '':
        plt.show()
    else:
        if png_too:
            plt.savefig(outfile+".pdf")
            plt.savefig(outfile+".png")
            plt.savefig(outfile+".eps")
            plt.close()
        else:
            plt.savefig(outfile)
            plt.close()


    return fig



# 'hide' parameter is hackish.
def plot_phase_core (ax, t, x, xerr, period, offset=0, 
                     sym='o', color='k', ms=6, hide=False):
    """ 
    Plots a pretty period-folded lightcurve on a given axes object.

    Doesn't assume anything about your data (e.g., that it's in magnitudes)
    
    Parameters
    ----------
    ax : plt.Axes
    t, x, xerr : array_like
    period : float
    offset : float, optional
        How much to shift the phase by. Default is zero.
    sym : str, optional
        Default 'o'. (circles)
    color : str, optional
        Default 'k'. (black)
    ms : float
        Default 6.
        
    Returns
    -------
    period : float
        The input period.
    
    """
    
    phase = ((t % period) / period + offset) % 1.


    if not hide:    ax.errorbar(phase, x, yerr=xerr, fmt= color+sym, ms=ms)
    ax.errorbar(phase-1,x,yerr=xerr,fmt=sym, mfc='0.7',mec='0.7', 
                 ecolor='0.7', ms=ms)
    ax.errorbar(phase+1,x,yerr=xerr,fmt=sym, mfc='0.7',mec='0.7', 
                 ecolor='0.7', ms=ms)
    
    ax.set_xticks( [0, 0.5, 1] )
    ax.set_xticks( np.arange(-.5,1.5,.1), minor=True)

    ax.set_xlim(-0.25, 1.25)

    return period


def phase (table, sid, period='auto', season=0, offset=0, 
           outfile='', name='', stetson=True, png_too=False):
    """ 
    Plots folded J, H, K lightcurves, plus color-color and color-mag
    trajectories, for one star.

    Will display "lonely" datapoints (i.e. not all JHK mags are 
    well-defined), and plots error-flagged data as different symbols.
    Compare to plot2.phase() which does neither.

    If no period is provided, phase() will compute its own guess 
    using the K-band lightcurve.

    Parameters
    ----------
    table : atpy.Table
        Table with time-series photometry
    sid : int
        13-digit WFCAM source ID of star to plot
    period : {'auto', 'lsp', float}, optional
        What period to fold the lightcurves by.
        If 'auto' is provided (default), the fast-chi-squared 
        period of the K-band data will be used.
        If 'lsp' is provided, the Lomb-Scargle Periodogram period
        of the K-band data will be used.
    season : int, optional
        Which observing season of our dataset (1, 2, 3, or all).
        Any value that is not the integers (1, 2, or 3) will be 
        treated as "no season", and no time-cut will be made.
        Note that this is the default behavior.
    offset : float, optional
        How much to shift the phase by. Default is zero.
    outfile : str, optional
        What filename to save plot to. Default behavior (when 
        `outfile` is an empty string) is to display plot on-screen
        and *not* save to file.
    name : str, optional
        A string to use as a plot header.
    png_too : bool, optional (default: False)
        If `png_too` is True (and `outfile` is not ''), then 
        save the plot in 3 file formats: PDF, PNG, and EPS.
        Do not specify a file extension in `outfile`.

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

    ## Let's do the 3 lightcurve bands.

    ### Initializing variables. Most of the following are done thrice:
    # Once with no errors (normal),
    # once with small errors (info),
    # once with large errors (warning).

    ## First: normal

    # Use band_cut to get relevant data chunks.

    j_table = band_cut(s_table, 'j', max_flag=0)
    h_table = band_cut(s_table, 'h', max_flag=0)
    k_table = band_cut(s_table, 'k', max_flag=0)


    # The Julian date for CE  2000 January  1 00:00:00.0 UT is
    # JD 2451544.500000
    # (i.e. MJD 51544.0)

    # get a date (x-axis) for each plot
    jdate = j_table.MEANMJDOBS - 51544
    hdate = h_table.MEANMJDOBS - 51544
    kdate = k_table.MEANMJDOBS - 51544
    
    # get a magnitude (y-axis) for each plot
    jcol = j_table.JAPERMAG3
    hcol = h_table.HAPERMAG3
    kcol = k_table.KAPERMAG3

    # get a magnitude error (y-error) for each plot
    jerr = j_table.JAPERMAG3ERR
    herr = h_table.HAPERMAG3ERR
    kerr = k_table.KAPERMAG3ERR

    # get a quality flag for each plot
    # (aftercomment: not sure we're going to ever use these)
    jflag = j_table.JPPERRBITS
    hflag = h_table.HPPERRBITS
    kflag = k_table.KPPERRBITS

    ## Second: info

    # Use band_cut to get relevant data chunks.

    j_table_info = band_cut(s_table, 'j', min_flag=1, max_flag=256)
    h_table_info = band_cut(s_table, 'h', min_flag=1, max_flag=256)
    k_table_info = band_cut(s_table, 'k', min_flag=1, max_flag=256)


    # The Julian date for CE  2000 January  1 00:00:00.0 UT is
    # JD 2451544.500000
    # (i.e. MJD 51544.0)

    # get a date (x-axis) for each plot
    jdate_info = j_table_info.MEANMJDOBS - 51544
    hdate_info = h_table_info.MEANMJDOBS - 51544
    kdate_info = k_table_info.MEANMJDOBS - 51544
    
    # get a magnitude (y-axis) for each plot
    jcol_info = j_table_info.JAPERMAG3
    hcol_info = h_table_info.HAPERMAG3
    kcol_info = k_table_info.KAPERMAG3

    # get a magnitude error (y-error) for each plot
    jerr_info = j_table_info.JAPERMAG3ERR
    herr_info = h_table_info.HAPERMAG3ERR
    kerr_info = k_table_info.KAPERMAG3ERR

    # get a quality flag for each plot
    jflag = j_table_info.JPPERRBITS
    hflag = h_table_info.HPPERRBITS
    kflag = k_table_info.KPPERRBITS
    
    ## Third: warning

    # Use band_cut to get relevant data chunks.

    j_table_warn = band_cut(s_table, 'j', min_flag=257)
    h_table_warn = band_cut(s_table, 'h', min_flag=257)
    k_table_warn = band_cut(s_table, 'k', min_flag=257)


    # The Julian date for CE  2000 January  1 00:00:00.0 UT is
    # JD 2451544.500000
    # (i.e. MJD 51544.0)

    # get a date (x-axis) for each plot
    jdate_warn = j_table_warn.MEANMJDOBS - 51544
    hdate_warn = h_table_warn.MEANMJDOBS - 51544
    kdate_warn = k_table_warn.MEANMJDOBS - 51544
    
    # get a magnitude (y-axis) for each plot
    jcol_warn = j_table_warn.JAPERMAG3
    hcol_warn = h_table_warn.HAPERMAG3
    kcol_warn = k_table_warn.KAPERMAG3

    # get a magnitude error (y-error) for each plot
    jerr_warn = j_table_warn.JAPERMAG3ERR
    herr_warn = h_table_warn.HAPERMAG3ERR
    kerr_warn = k_table_warn.KAPERMAG3ERR

    # get a quality flag for each plot
    jflag = j_table_warn.JPPERRBITS
    hflag = h_table_warn.HPPERRBITS
    kflag = k_table_warn.KPPERRBITS

    try:
    ## Let's figure out the period.
        if period == 'auto':
            kper_table = band_cut(s_table, 'k', max_flag=256)
            kperdate = kper_table.MEANMJDOBS
            kpercol = kper_table.KAPERMAG3
            kpererr = kper_table.KAPERMAG3ERR
            
            period = 1./test_analyze(kperdate, kpercol, kpererr)
            print period
        elif period == 'lsp':
            kper_table = band_cut(s_table, 'k', max_flag=256)
            kperdate = kper_table.MEANMJDOBS
            kpercol = kper_table.KAPERMAG3
            kpererr = kper_table.KAPERMAG3ERR
            
            lomb = lsp(kperdate,kpercol,6.,6.)
            lsp_freq = lomb[0]
            lsp_power= lomb[1]
            Jmax = lsp_mask( lsp_freq, lsp_power)
            lsp_per = 1./ lomb[0][Jmax]
            period = lsp_per
            print period
    except Exception:
        print 'period-finding failed! returning'
        return

    if np.abs(period) < 1:
        period_string = "%f hours" % (period*24)
        print period_string
    else:
        period_string = "%f days" % period

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
    fmt_warn = '.'

    # Plot J-band:
    if len(jdate) > 0:
#        ax_j.errorbar( jdate, jcol, yerr=jerr, fmt='bo', ecolor='k')
        plot_phase_core( ax_j, jdate, jcol, jerr, period, offset=offset,
                         color='b')
    if len(jdate_info) > 0:
#        ax_j.errorbar( jdate_info, jcol_info, yerr=jerr_info, 
#                       fmt='b'+fmt_info, ms=4, ecolor='k')
        plot_phase_core( ax_j, jdate_info, jcol_info, jerr_info, period, 
                         offset=offset, color='b', sym=fmt_info, ms=4)
    if len(jdate_warn) > 0:
#        ax_j.errorbar( jdate_warn, jcol_warn, yerr=jerr_warn, 
#                       fmt='k'+fmt_warn, ecolor='k')
        plot_phase_core( ax_j, jdate_warn, jcol_warn, jerr_warn, period,
                         offset=offset, sym=fmt_warn)
    ax_j.invert_yaxis()

    # Plot H-band:
    if len(hdate) > 0:
#        ax_h.errorbar( hdate, hcol, yerr=herr, fmt='go', ecolor='k' )
        plot_phase_core( ax_h, hdate, hcol, herr, period, offset=offset,
                         color='g')
    if len(hdate_info) > 0:
#        ax_h.errorbar( hdate_info, hcol_info, yerr=herr_info, 
#                       fmt='g'+fmt_info, ms=4, ecolor='k')
        plot_phase_core( ax_h, hdate_info, hcol_info, herr_info, period, 
                         offset=offset, color='g', sym=fmt_info, ms=4)
    if len(hdate_warn) > 0:
#        ax_h.errorbar( hdate_warn, hcol_warn, yerr=herr_warn, 
#                       fmt='k'+fmt_warn, ecolor='k')
        plot_phase_core( ax_h, hdate_warn, hcol_warn, herr_warn, period,
                         offset=offset, sym=fmt_warn)
    ax_h.invert_yaxis()

    # Plot K-band:
    if len(kdate) > 0:
#        ax_k.errorbar( kdate, kcol, yerr=kerr, fmt='ro', ecolor='k' )
        plot_phase_core( ax_k, kdate, kcol, kerr, period, offset=offset,
                         color='r')

    if len(kdate_info) > 0:
#        ax_k.errorbar( kdate_info, kcol_info, yerr=kerr_info, 
#                       fmt='r'+fmt_info, ms=4, ecolor='k')
        plot_phase_core( ax_k, kdate_info, kcol_info, kerr_info, period, 
                         offset=offset, color='r', sym=fmt_info, ms=4)
    if len(kdate_warn) > 0:
#        ax_k.errorbar( kdate_warn, kcol_warn, yerr=kerr_warn, 
#                       fmt='k'+fmt_warn, ecolor='k')
        plot_phase_core( ax_k, kdate_warn, kcol_warn, kerr_warn, period,
                         offset=offset, sym=fmt_warn)
    ax_k.invert_yaxis()

    ## Now let's do the 2 color-mag/color-color plots.

    # We'll use different data-cuts for the two different plots.
    # Relevant comment: I made an executive call to include only
    # 'normal' and 'info'-flagged data in the C-C and C-M plots
    # (i.e. max_flag=256 in all relevant bands).
    
    # In the color-mag plot, we need data where H and K are defined 
    # everywhere. That's two cuts.
    khk_table = band_cut( band_cut(s_table, 'k', max_flag=256),
                          'h', max_flag=256)

    khkdate = khk_table.MEANMJDOBS - 51544
    k_khk = khk_table.KAPERMAG3
    hmk_khk = khk_table.HMKPNT

    khkphase = ((khkdate % period) / period + offset) % 1.

    # In the color-color plot, we need data where J, H, and K are
    # defined everywhere. That's one more cut.
    jhk_table = band_cut(khk_table, 'j', max_flag=256)

    jhkdate = jhk_table.MEANMJDOBS - 51544
    jmh_jhk = jhk_table.JMHPNT
    hmk_jhk = jhk_table.HMKPNT

    jhkphase = ((jhkdate % period) / period + offset) % 1.

    # Plot J-H vs H-K using the "jhk_" variables.
    try:
        plot_trajectory_core( ax_jhk, hmk_jhk, jmh_jhk, jhkphase )
    except Exception:
        print "JHK plot broke!"
        pass

    # Plot K vs H-K using the "khk_" variables.
    try:
        plot_trajectory_core( ax_khk, hmk_khk, k_khk, khkphase, 
                              ms=False, ctts=False) 
    except Exception:
        print "KHK plot broke!"
        pass
    ax_khk.invert_yaxis()


    # Hide the bad labels...
    plt.setp(ax_j.get_xticklabels(), visible=False)
    plt.setp(ax_h.get_xticklabels(), visible=False)

    # Label stuff
    ax_k.set_xlabel( "Phase (Period = %s)" % period_string )

    ax_j.set_ylabel( "J",{'rotation':'horizontal', 'fontsize':'large'} )
    ax_h.set_ylabel( "H",{'rotation':'horizontal', 'fontsize':'large'} )
    ax_k.set_ylabel( "K",{'rotation':'horizontal', 'fontsize':'large'} )

    ax_jhk.set_xlabel( "H-K" )
    ax_jhk.set_ylabel( "J-H")#, {'rotation':'horizontal'})
    ax_khk.set_xlabel( "H-K" )
    ax_khk.set_ylabel( "K")#, {'rotation':'horizontal'})

    if name != '':
        ax_j.set_title(name)

    if stetson == True:
        S, choice, n = Stetson_machine( s_table, flags=256 )
        stet_string = "S = %.2f" % S
        ax_khk.set_title(stet_string)
        
    if outfile == '':
        plt.show()
    else:
        if png_too:
            plt.savefig(outfile+".pdf")
            plt.savefig(outfile+".png")
            plt.savefig(outfile+".eps")
            plt.close()
        else:
            plt.savefig(outfile)
            plt.close()


    return fig

def lsp_power (table, sid, season=123, outfile='', name='', png_too=False):
    """ 
    Plots J, H, K periodograms for one star.

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
    outfile : str, optional
        What filename to save plot to. Default behavior (when 
        `outfile` is an empty string) is to display plot on-screen
        and *not* save to file.
    png_too : bool, optional (default: False)
        If `png_too` is True (and `outfile` is not ''), then 
        save the plot in 3 file formats: PDF, PNG, and EPS.
        Do not specify a file extension in `outfile`.
        
    Returns
    -------
    fig : plt.Figure 
        The canvas figure that the graphs are plotted onto.
    
    """

    ## Loading data
    s_table = data_cut (table, sid, season)

    if len(s_table) < 2:
        print "no data here"
        return

    ## Set up plot

    fig = plt.figure(figsize = (10, 6), dpi=80,
                     facecolor='w', edgecolor='k')

    ax_j = fig.add_subplot(3,1,1)
    ax_h = fig.add_subplot(3,1,2, sharex=ax_j)
    ax_k = fig.add_subplot(3,1,3, sharex=ax_j)

    ## J
    try:
    # do some band_cutting, with flags = 256 like usual
        j_table = band_cut(s_table, 'j', max_flag=256)
        jdate = j_table.MEANMJDOBS - 51544

    # get a magnitude (y-axis) for each plot
        jcol = j_table.JAPERMAG3

    ## Calculate periodograms
        jlsp = lsp(jdate, jcol, 6., 6.)
        j_lsp_freq = jlsp[0]
        j_lsp_power = jlsp[1]
        
    # best periods, filtered by the lsp_mask
        j_lsp_per = 1./ j_lsp_freq[ lsp_mask( j_lsp_freq, j_lsp_power) ]    

    ## Plot things

        ax_j.plot(1./j_lsp_freq, j_lsp_power, 'b')
    except Exception:
        pass

    ## H
    try:
        h_table = band_cut(s_table, 'h', max_flag=256)
        hdate = h_table.MEANMJDOBS - 51544
        hcol = h_table.HAPERMAG3
        
        hlsp = lsp(hdate, hcol, 6., 6.)
        h_lsp_freq = hlsp[0]
        h_lsp_power = hlsp[1]
        
        h_lsp_per = 1./ h_lsp_freq[ lsp_mask( h_lsp_freq, h_lsp_power) ]
        
        ax_h.plot(1./h_lsp_freq, h_lsp_power, 'g')
        
    except Exception:
        pass


    ## K
    try:
        k_table = band_cut(s_table, 'k', max_flag=256)
        kdate = k_table.MEANMJDOBS - 51544
        kcol = k_table.KAPERMAG3
        
        klsp = lsp(kdate, kcol, 6., 6.)
        k_lsp_freq = klsp[0]
        k_lsp_power = klsp[1]
        
        k_lsp_per = 1./ k_lsp_freq[ lsp_mask( k_lsp_freq, k_lsp_power) ]

        ax_k.plot(1./k_lsp_freq, k_lsp_power, 'r')
    except Exception:
        pass
    
    ax_j.set_xscale('log')
    ax_h.set_xscale('log')
    ax_k.set_xscale('log')

    ax_j.set_title(name)
    ax_k.set_xlabel("Period (days)")
    ax_h.set_ylabel("Periodogram Power")

    ## Save things

    if outfile == '':
        plt.show()
    else:
        if png_too:
            plt.savefig(outfile+".pdf")
            plt.savefig(outfile+".png")
            plt.savefig(outfile+".eps")
            plt.close()
        else:
            plt.savefig(outfile)
            plt.close()

    return fig
