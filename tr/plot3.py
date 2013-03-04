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

from __future__ import division

import numpy as np
import matplotlib.pyplot as plt

import atpy

from helpers3 import data_cut, band_cut
from plot2 import plot_trajectory_core
from chi2 import test_analyze, diagnostic_analyze
from scargle import fasper as lsp
from timing import lsp_mask, lsp_tuning
from spread3 import Stetson_machine
from abridger import abridger
from color_slope import slope

#import coords
#import stetson

#from tr_helpers import season_cut, data_cut, ensemble_cut



def lc (table, sid, season=0, outfile='', name='', 
        stetson=True, png_too=False, 
        date_offset = 51544, color_slope=False, custom_xlabel=None, 
        plot_warn=True):
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
    date_offset : float, optional
        What MJD to use as day "zero". Default 01/01/2000, 
        aka MJD=51544. Date corresponding to the start of Cyg OB7
        observations (04/23/2008) is MJD=54579.
    color_slope : bool, optional (default: False)
        Whether to fit color slope lines to the KvH-K and J-HvH-K plots.
    custom_xlabel : str or None, optional
        What string to use (e.g. "Time (JD since 04/23/2008)") instead of
        the default "Time (MJD - `date_offset`)" for the lightcurve
        x-axis label.
    plot_warn : bool, optional
        Do you want to plot data with "warning" error flags?
        Default True for backwards compatibility, but you probably want
        False for readability.


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
    jdate = j_table.MEANMJDOBS - date_offset
    hdate = h_table.MEANMJDOBS - date_offset
    kdate = k_table.MEANMJDOBS - date_offset
    
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
    jdate_info = j_table_info.MEANMJDOBS - date_offset
    hdate_info = h_table_info.MEANMJDOBS - date_offset
    kdate_info = k_table_info.MEANMJDOBS - date_offset
    
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
    jdate_warn = j_table_warn.MEANMJDOBS - date_offset
    hdate_warn = h_table_warn.MEANMJDOBS - date_offset
    kdate_warn = k_table_warn.MEANMJDOBS - date_offset
    
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
    if (len(jdate_warn) > 0) and plot_warn:
        ax_j.errorbar( jdate_warn, jcol_warn, yerr=jerr_warn, 
                       fmt='k'+fmt_warn, ecolor='k')
    ax_j.invert_yaxis()

    # Plot H-band:
    if len(hdate) > 0:
        ax_h.errorbar( hdate, hcol, yerr=herr, fmt='go', ecolor='k' )
    if len(hdate_info) > 0:
        ax_h.errorbar( hdate_info, hcol_info, yerr=herr_info, 
                       fmt='g'+fmt_info, ms=4, ecolor='k')
    if (len(hdate_warn) > 0) and plot_warn:
        ax_h.errorbar( hdate_warn, hcol_warn, yerr=herr_warn, 
                       fmt='k'+fmt_warn, ecolor='k')
    ax_h.invert_yaxis()

    # Plot K-band:
    if len(kdate) > 0:
        ax_k.errorbar( kdate, kcol, yerr=kerr, fmt='ro', ecolor='k' )
    if len(kdate_info) > 0:
        ax_k.errorbar( kdate_info, kcol_info, yerr=kerr_info, 
                       fmt='r'+fmt_info, ms=4, ecolor='k')
    if (len(kdate_warn) > 0) and plot_warn:
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

    khkdate = khk_table.MEANMJDOBS - date_offset
    k_khk = khk_table.KAPERMAG3
    hmk_khk = khk_table.HMKPNT
    k_khk_err = khk_table.KAPERMAG3ERR
    hmk_khk_err = khk_table.HMKPNTERR

    # In the color-color plot, we need data where J, H, and K are
    # defined everywhere. That's one more cut.
    jhk_table = band_cut(khk_table, 'j', max_flag=256)

    jhkdate = jhk_table.MEANMJDOBS - date_offset
    jmh_jhk = jhk_table.JMHPNT
    hmk_jhk = jhk_table.HMKPNT
    jmh_jhk_err = jhk_table.JMHPNTERR
    hmk_jhk_err = jhk_table.HMKPNTERR

    # Plot J-H vs H-K using the "jhk_" variables.
    try:
        plot_trajectory_core( ax_jhk, hmk_jhk, jmh_jhk, jhkdate )

        if color_slope:
            jhk_slope, jhk_intercept, slope_err = (
                slope(hmk_jhk, jmh_jhk, hmk_jhk_err, jmh_jhk_err,
                      verbose=False) )
            
            ax_jhk.plot([0, 6], [jhk_intercept, jhk_intercept + 6*jhk_slope], 
                        ':', scalex=False, scaley=False)

    except Exception, e:
        print "JHK plot broke!", e
        pass

    # Plot K vs H-K using the "khk_" variables.
    try:
        plot_trajectory_core( ax_khk, hmk_khk, k_khk, khkdate, 
                              ms=False, ctts=False) 

        # plot boundaries are manually set for readability, if necessary
        if len(ax_khk.get_xticks()) > 7:
            khk_xmin = np.floor(hmk_khk.min() * 0.95 * 20)/20.
            khk_xmax = np.ceil( hmk_khk.max() * 1.05 * 20)/20.

            khk_xticks = np.linspace(khk_xmin, khk_xmax, 6)
            ax_khk.set_xticks(khk_xticks)

        if color_slope:
            khk_slope, khk_intercept, slope_err = (
                slope(hmk_khk, k_khk, hmk_khk_err, k_khk_err,
                      verbose=False) )
            
            ax_khk.plot([0, 6], [khk_intercept, khk_intercept + 6*khk_slope],
                        '--', scalex=False, scaley=False)
    
    except Exception, e:
        print "KHK plot broke!", e
        pass
    ax_khk.invert_yaxis()

    # Hide the bad labels...
    plt.setp(ax_j.get_xticklabels(), visible=False)
    plt.setp(ax_h.get_xticklabels(), visible=False)

    # Label stuff
#    ax_k.set_xlabel( "Time (JD since 01/01/2000)" )
    if custom_xlabel:
        ax_k.set_xlabel( custom_xlabel )
    else:
        ax_k.set_xlabel( "Time (MJD - %.1f)" % date_offset )


    ax_j.set_ylabel( "J",{'rotation':'horizontal', 'fontsize':'large'} )
    ax_h.set_ylabel( "H",{'rotation':'horizontal', 'fontsize':'large'} )
    ax_k.set_ylabel( "K",{'rotation':'horizontal', 'fontsize':'large'} )

    ax_jhk.set_xlabel( "H-K" )
    ax_jhk.set_ylabel( "J-H")#, {'rotation':'horizontal'})
    ax_khk.set_xlabel( "H-K" )
    ax_khk.set_ylabel( "K")#, {'rotation':'horizontal'})

    if name != '':
        ax_j.set_title(name)
    else:
        ax_j.set_title(str(sid))

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


    fig.ax_k = ax_k
    fig.ax_h = ax_h
    fig.ax_j = ax_j
    fig.ax_jhk = ax_jhk
    fig.ax_khk = ax_khk

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
           outfile='', name='', stetson=True, png_too=False,
           color_slope=False, period_decimal_places = 6):
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
    color_slope : bool, optional (default: False)
        Whether to fit color slope lines to the KvH-K and J-HvH-K plots.
    period_decimal_places : int, optional (default: 6)
        How many decimal places to print of the period on the
        phase plot's x-label. 

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
            
            hifac = lsp_tuning(kperdate)
            lomb = lsp(kperdate,kpercol,6., hifac)
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
        period_string = "%.*f hours" % (period_decimal_places, period*24)
#        print period_string
    else:
        period_string = "%.*f days" % (period_decimal_places, period)

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
    k_khk_err = khk_table.KAPERMAG3ERR
    hmk_khk_err = khk_table.HMKPNTERR

    khkphase = ((khkdate % period) / period + offset) % 1.

    # In the color-color plot, we need data where J, H, and K are
    # defined everywhere. That's one more cut.
    jhk_table = band_cut(khk_table, 'j', max_flag=256)

    jhkdate = jhk_table.MEANMJDOBS - 51544
    jmh_jhk = jhk_table.JMHPNT
    hmk_jhk = jhk_table.HMKPNT
    jmh_jhk_err = jhk_table.JMHPNTERR
    hmk_jhk_err = jhk_table.HMKPNTERR

    jhkphase = ((jhkdate % period) / period + offset) % 1.

    # Plot J-H vs H-K using the "jhk_" variables.
    try:
        plot_trajectory_core( ax_jhk, hmk_jhk, jmh_jhk, jhkphase )

        if color_slope:
            jhk_slope, jhk_intercept, slope_err = (
                slope(hmk_jhk, jmh_jhk, hmk_jhk_err, jmh_jhk_err,
                      verbose=False) )

            ax_jhk.plot([0, 6], [jhk_intercept, jhk_intercept + 6*jhk_slope], 
                        ':', scalex=False, scaley=False)

    except Exception, e:
        print "JHK plot broke!", e
        pass

    # Plot K vs H-K using the "khk_" variables.
    try:
        plot_trajectory_core( ax_khk, hmk_khk, k_khk, khkphase, 
                              ms=False, ctts=False) 
        # plot boundaries are manually set for readability, if necessary
        if len(ax_khk.get_xticks()) > 7:
            khk_xmin = np.floor(hmk_khk.min() * 0.95 * 20)/20.
            khk_xmax = np.ceil( hmk_khk.max() * 1.05 * 20)/20.

            khk_xticks = np.linspace(khk_xmin, khk_xmax, 6)
            ax_khk.set_xticks(khk_xticks)

        if color_slope:
            khk_slope, khk_intercept, slope_err = (
                slope(hmk_khk, k_khk, hmk_khk_err, k_khk_err,
                      verbose=False) )
            
            ax_khk.plot([0, 6], [khk_intercept, khk_intercept + 6*khk_slope],
                        '--', scalex=False, scaley=False)

    except Exception, e:
        print "KHK plot broke!", e
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
    else:
        ax_j.set_title(str(sid))

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

    fig.ax_k = ax_k
    fig.ax_h = ax_h
    fig.ax_j = ax_j
    fig.ax_jhk = ax_jhk
    fig.ax_khk = ax_khk

    return fig


def jjh (table, sid, season=0, outfile='', name='', 
         stetson=True, png_too=False,         
         date_offset = 51544, color_slope=False):
    """ 
    Plots a J, J-H color-color diagram for a single star, colored by time.

    Parameters
    ----------
    table : atpy.Table
        Table with time-series photometry and grade columns. 
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
    date_offset : float, optional
        What MJD to use as day "zero". Default 01/01/2000, 
        aka MJD=51544, unless `abridged`=True, in which case
        MJD = 54034 is set (initial observations of Orion dataset).
    color_slope : bool, optional (defalt: False)
        Whether to fit a color slope line to the J, J-H plot.
        
    Returns
    -------
    fig : plt.Figure 
        The canvas figure that the graph is plotted onto.
    
    """
    
    # Loading data
    s_table = data_cut (table, sid, season)

    if len(s_table) == 0:
        print "no data here"
        return

    ## Define the components and parameters of the figure:
    
    fig = plt.figure(figsize = (4.286, 3.2143), dpi=80,
                     facecolor='w', edgecolor='k')

    bottom = 0.1
#    height = .25
#    left = 0.075
#    width = 0.5

    ax_jjh = fig.add_axes( (.2, .15, .7, .7) )

    color_vmin = s_table.MEANMJDOBS.min() - date_offset # PROTECT
    color_vmax = s_table.MEANMJDOBS.max() - date_offset # PROTECT
    ## Now let's do the 2 color-mag/color-color plots.

    # Relevant comment: I made an executive call to include only
    # 'normal' and 'info'-flagged data in the C-M plot
    # (i.e. max_flag=256 in all relevant bands).
    
    # In the color-mag plot, we need data where H and J are defined 
    # everywhere. That's two cuts.
    jjh_table = band_cut( band_cut(s_table, 'j', max_flag=256),
                          'h', max_flag=256)

    jjhdate = jjh_table.MEANMJDOBS - date_offset
    j_jjh = jjh_table.JAPERMAG3
    jmh_jjh = jjh_table.JMHPNT
    j_jjh_err = jjh_table.JAPERMAG3ERR
    jmh_jjh_err = jjh_table.JMHPNTERR

        

    # Plot J vs J-H using the "jjh_" variables.
    try:
        plot_trajectory_core( ax_jjh, jmh_jjh, j_jjh, jjhdate,
                              ms=False, ctts=False, 
                              vmin=color_vmin, vmax=color_vmax) 

        # plot boundaries are manually set for readability, if necessary
        if len(ax_jjh.get_xticks()) > 7:
            jjh_xmin = np.floor(jmh_jjh.min() * 0.95 * 20)/20.
            jjh_xmax = np.ceil( jmh_jjh.max() * 1.05 * 20)/20.

            jjh_xticks = np.linspace(jjh_xmin, jjh_xmax, 6)
            ax_jjh.set_xticks(jjh_xticks)

        if color_slope:
            jjh_slope, jjh_intercept, slope_err = (
                slope(jmh_jjh, j_jjh, jmh_jjh_err, j_jjh_err,
                      verbose=False) )
            
            ax_jjh.plot([0, 6], [jjh_intercept, jjh_intercept + 6*jjh_slope],
                        '--', scalex=False, scaley=False)
    
    except Exception:
        print "JJH plot broke!"
        pass
    ax_jjh.invert_yaxis()


    # Label stuff
#    ax_k.set_xlabel( "Time (JD since 01/01/2000)" )
#    ax_j.set_xlabel( "Time (MJD - %.1f)" % date_offset )

    ax_jjh.set_xlabel( "J-H" )
    ax_jjh.set_ylabel( "J", {'rotation':'horizontal'})

    if name != '':
        ax_jjh.set_title(name)
    else:
        ax_jjh.set_title(str(sid))

    if stetson == True:
        S, choice, n = Stetson_machine( s_table, flags=256 )
        stet_string = "S = %.2f" % S
        suptitle = plt.suptitle(ax_jjh.get_title()+", "+stet_string) 
        ax_jjh.set_title('')

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

    fig.ax_jjh = ax_jjh

    return fig



def lsp_power (table, sid, season=0, upper_frequency=0.5,
               outfile='', name='', png_too=False):
    """ 
    Plots J, H, K lomb-scargle periodograms for one star.

    Parameters
    ----------
    table : atpy.Table
        Table with WFCAM time-series photometry
    sid : int
        13-digit WFCAM source ID of star to plot
    season : int, optional
        Which observing season of our dataset (1, 2, 3, or all).
        Any value that is not the integers (1, 2, 3, or 123) will be 
        treated as "no season", and no time-cut will be made.
        Note that this is the default behavior.
    upper_frequency : float, optional
        Highest frequency (1/lowest period) that one desires to scan over
        in the period search. Defaults to 0.5 day^-1 (period of 2 days).
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

    # do some band_cutting, with flags = 256 like usual
    j_table = band_cut(s_table, 'j', max_flag=256)
    if len(j_table) <= 1:
        print "no J data"
    else:
        jdate = j_table.MEANMJDOBS - 51544

    # get a magnitude (y-axis) for each plot
        jcol = j_table.JAPERMAG3

    ## Calculate periodograms
        hifac = lsp_tuning(jdate, upper_frequency=upper_frequency)
        jlsp = lsp(jdate, jcol, 6., hifac)
        j_lsp_freq = jlsp[0]
        j_lsp_power = jlsp[1]
        
    # best periods, filtered by the lsp_mask
        j_lsp_per = 1./ j_lsp_freq[ lsp_mask( j_lsp_freq, j_lsp_power) ]    

    ## Plot things

        ax_j.plot(1./j_lsp_freq, j_lsp_power, 'b')
#    except Exception:
#        print "J power broke!"
#        pass

    ## H
    h_table = band_cut(s_table, 'h', max_flag=256)
    if len(h_table) <= 1:
        print "no H data"
    else:
        hdate = h_table.MEANMJDOBS - 51544
        hcol = h_table.HAPERMAG3
        
        hifac = lsp_tuning(hdate, upper_frequency=upper_frequency)
        hlsp = lsp(hdate, hcol, 6., hifac)
        h_lsp_freq = hlsp[0]
        h_lsp_power = hlsp[1]
        
        h_lsp_per = 1./ h_lsp_freq[ lsp_mask( h_lsp_freq, h_lsp_power) ]
        
        ax_h.plot(1./h_lsp_freq, h_lsp_power, 'g')
        
#    except Exception:
#        print "H power broke!"
#        pass


    ## K
    k_table = band_cut(s_table, 'k', max_flag=256)
    if len(k_table) <= 1:
        print "no K data"
    else:
        kdate = k_table.MEANMJDOBS - 51544
        kcol = k_table.KAPERMAG3

        hifac = lsp_tuning(kdate, upper_frequency=upper_frequency)
        klsp = lsp(kdate, kcol, 6., hifac)
        k_lsp_freq = klsp[0]
        k_lsp_power = klsp[1]
        
        k_lsp_per = 1./ k_lsp_freq[ lsp_mask( k_lsp_freq, k_lsp_power) ]

        ax_k.plot(1./k_lsp_freq, k_lsp_power, 'r')
#    except Exception:
#        print "K power broke!"
#        pass
    
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

    fig.ax_k = ax_k
    fig.ax_h = ax_h
    fig.ax_j = ax_j

    return fig

def fx2_periodogram (table, sid, season=0, 
                     outfile='', name='', png_too=False):
    """ 
    Plots J, H, K fast chi-squared periodograms for one star.

    Parameters
    ----------
    table : atpy.Table
        Table with WFCAM time-series photometry
    sid : int
        13-digit WFCAM source ID of star to plot
    season : int, optional
        Which observing season of our dataset (1, 2, 3, or all).
        Any value that is not the integers (1, 2, 3, or 123) will be 
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

    # do some band_cutting, with flags = 256 like usual
    j_table = band_cut(s_table, 'j', max_flag=256)
    if len(j_table) <= 1:
        print "no J data"
    else:
        jdate = j_table.MEANMJDOBS - 51544

    # get a magnitude (y-axis) for each plot
        jcol = j_table.JAPERMAG3
        jerr = j_table.JAPERMAG3ERR

    ## Calculate periodograms
        j_fx2_freq, j_fx2_csr = diagnostic_analyze(jdate, jcol, jerr)
        
    # best periods, filtered by the lsp_mask
#        j_lsp_per = 1./ j_lsp_freq[ lsp_mask( j_lsp_freq, j_lsp_power) ]    

    ## Plot things

        ax_j.plot(1./j_fx2_freq, j_fx2_csr, 'b')

#    except Exception:
#        print "J power broke!"
#        pass

    ## H
    h_table = band_cut(s_table, 'h', max_flag=256)
    if len(h_table) <= 1:
        print "no H data"
    else:
        hdate = h_table.MEANMJDOBS - 51544
        hcol = h_table.HAPERMAG3
        herr = h_table.HAPERMAG3ERR

    ## Calculate periodograms
        h_fx2_freq, h_fx2_csr = diagnostic_analyze(hdate, hcol, herr)
        
        ax_h.plot(1./h_fx2_freq, h_fx2_csr, 'g')
        
#    except Exception:
#        print "H power broke!"
#        pass


    ## K
    k_table = band_cut(s_table, 'k', max_flag=256)
    if len(k_table) <= 1:
        print "no K data"
    else:
        kdate = k_table.MEANMJDOBS - 51544
        kcol = k_table.KAPERMAG3
        kerr = k_table.KAPERMAG3ERR

    ## Calculate periodograms
        k_fx2_freq, k_fx2_csr = diagnostic_analyze(kdate, kcol, kerr)
        
        ax_k.plot(1./k_fx2_freq, k_fx2_csr, 'r')

#    except Exception:
#        print "K power broke!"
#        pass
    
    ax_j.set_xscale('log')
    ax_h.set_xscale('log')
    ax_k.set_xscale('log')

    ax_j.set_title(name)
    ax_k.set_xlabel("Period (days)")
    ax_h.set_ylabel(r"$\chi^2$ reduction")

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

    fig.ax_k = ax_k
    fig.ax_h = ax_h
    fig.ax_j = ax_j

    return fig


def scatter_phase_core (ax, t, x, xerr, period, offset=0, 
                        ms=6, sym='o', hide=False, **kwargs):
    """ 
    Scatter-plots a period-folded lightcurve on a given axes object.

    Also plots errorbars underneath.
    Doesn't assume anything about your data (e.g., that it's in magnitudes)
    
    Parameters
    ----------
    ax : plt.Axes
        An Axes object to plot onto.
    t, x, xerr : array_like
        The time, value, and uncertainty arrays to plot.
    period : float
        The period to fold the curve by.
    offset : float, optional
        How much to shift the phase by. Default is zero.
    sym : str, optional
        Symbol for gray dudes on the sides. Default 'o'. (circles)
    color : str, optional
        Color of the errorbars. Default 'k'. (black)
    ms : float
        Default 6.
    **kwargs : keyword arguments for plt.scatter()
        Suggested kwargs: `c` for color-array, `cmap` for colormap choice,
        `vmin` and `vmax` for the limits of colormap. Also possible:
        s (size), other things too I guess.
        
    Returns
    -------
    period : float
        The input period, unchanged.
    
    """
    
    # Calculate the "phase" variable as a function of time and the period.
    phase = ((t % period) / period + offset) % 1.


    # Plot the grayed-out guys on the left and right:
    ax.errorbar(phase-1,x,yerr=xerr,fmt=sym, mfc='0.7',mec='0.7', 
                 ecolor='0.7', ms=ms)
    ax.errorbar(phase+1,x,yerr=xerr,fmt=sym, mfc='0.7',mec='0.7', 
                 ecolor='0.7', ms=ms)

    # Now plot our actual scattered dude
    if not hide:    
        # errorbars in the background
        ax.errorbar(phase, x, yerr=xerr, fmt= None, ecolor='k', zorder=0)
        # scatter in the foreground
        ax.scatter(phase, x, zorder=100, **kwargs)
        
    ax.set_xticks( [0, 0.5, 1] )
    ax.set_xticks( np.arange(-.5,1.5,.1), minor=True)

    ax.set_xlim(-0.25, 1.25)

    return period



def graded_lc (table, sid, season=0, outfile='', name='', 
               stetson=True, png_too=False, abridged=False,
               timecolor=False, time_cmap='jet',
               d_cmap={'j':'Blues', 'h': 'Greens', 'k': 'Reds'},
               date_offset = 51544, color_slope=False):
    """ 
    Plots JHK lightcurves of a star, with datapoints colored by grade.

    Also plots color-color and color-mag trajectories, colored by time.

    Lightcurve datapoints can match the "time" coloration by setting
    `timecolor`=True.

    Will display "lonely" datapoints (i.e. not all JHK mags are 
    well-defined), and plots error-flagged data as different symbols.
    Compare to plot2.lc() which does neither, and to plot3.lc() which 
    is unaware of grades.

    Parameters
    ----------
    table : atpy.Table
        Table with time-series photometry and grade columns. 
        "JGRADE", "HGRADE", "KGRADE" must be bestowed by 
        night_cleanser.null_cleanser_grader().
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
    abridged : bool, optional (default: False)
        Create an abridged, panel-like plot?
    timecolor : bool, optional (default: False)
        Color lightcurve datapoints by time?
        Attempts to sync all colored plots.
    time_cmap : str, optional (defalt: 'jet')
        Which colormap to use for `timecolor`; when `timecolor`=True,
        this overrides d_cmap.
        Future possibility: the KHK and JHK plots using this cmap.
    d_cmap : dict or str or tuple
        Which colormaps to use for J, H, and K. If a single string 
        (rather than a dict) is given, all 3 bands will use the 
        same colormap. If a tuple is given, J:0, H:1, K:2.
        Default {'j':'Blues', 'h': 'Greens', 'k': 'Reds'} for now.
    date_offset : float, optional
        What MJD to use as day "zero". Default 01/01/2000, 
        aka MJD=51544, unless `abridged`=True, in which case
        MJD = 54034 is set (initial observations of Orion dataset).
    color_slope : bool, optional (defalt: False)
        Whether to fit color slope lines to the KvH-K and J-HvH-K plots.
        
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

    if abridged: 
        date_offset = 54034
        abridger_stuff = abridger(s_table, 54034, flags=256)

        ab_s2sub = abridger_stuff[0]
        ab_s3sub = abridger_stuff[1]
        ab_s1s2line = abridger_stuff[2]
        ab_s2s3line = abridger_stuff[3]
        ab_xticks = abridger_stuff[4]
        ab_xticklab = abridger_stuff[5]
        ab_xlim = abridger_stuff[6]

    ## Let's do the 3 lightcurve bands.

    ### Initializing variables. Most of the following are done twice:
    # Once with no errors (normal),
    # once with small errors (info).

    ## First: normal

    # Use band_cut to get relevant data chunks.

    j_table = band_cut(s_table, 'j', max_flag=0)
    h_table = band_cut(s_table, 'h', max_flag=0)
    k_table = band_cut(s_table, 'k', max_flag=0)


    # The Julian date for CE  2000 January  1 00:00:00.0 UT is
    # JD 2451544.500000
    # (i.e. MJD 51544.0)

    # get a date (x-axis) for each plot
    jdate = j_table.MEANMJDOBS - date_offset
    hdate = h_table.MEANMJDOBS - date_offset
    kdate = k_table.MEANMJDOBS - date_offset

    # let's save the old versions in case we need em later
    raw_jdate = np.copy(jdate)
    raw_hdate = np.copy(hdate)
    raw_kdate = np.copy(kdate)

    if abridged:
        jdate[jdate > 54300-date_offset] -= ab_s2sub
        jdate[jdate > (54600-date_offset) - ab_s2sub] -= ab_s3sub
        hdate[hdate > 54300-date_offset] -= ab_s2sub 
        hdate[hdate > (54600-date_offset) - ab_s2sub] -= ab_s3sub
        kdate[kdate > 54300-date_offset] -= ab_s2sub 
        kdate[kdate > (54600-date_offset) - ab_s2sub] -= ab_s3sub

    
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

    # Get the grade
    jgrade = j_table.JGRADE
    hgrade = h_table.HGRADE
    kgrade = k_table.KGRADE

    ## Second: info

    # Use band_cut to get relevant data chunks.

    j_table_info = band_cut(s_table, 'j', min_flag=1, max_flag=256)
    h_table_info = band_cut(s_table, 'h', min_flag=1, max_flag=256)
    k_table_info = band_cut(s_table, 'k', min_flag=1, max_flag=256)


    # The Julian date for CE  2000 January  1 00:00:00.0 UT is
    # JD 2451544.500000
    # (i.e. MJD 51544.0)

    # get a date (x-axis) for each plot
    jdate_info = j_table_info.MEANMJDOBS - date_offset
    hdate_info = h_table_info.MEANMJDOBS - date_offset
    kdate_info = k_table_info.MEANMJDOBS - date_offset

    # let's save the old versions in case we need em later
    raw_jdate_info = np.copy(jdate_info)
    raw_hdate_info = np.copy(hdate_info)
    raw_kdate_info = np.copy(kdate_info)

    if abridged:
        jdate_info[jdate_info > 54300-date_offset] -= ab_s2sub
        jdate_info[jdate_info > (54600-date_offset) - ab_s2sub] -= ab_s3sub
        hdate_info[hdate_info > 54300-date_offset] -= ab_s2sub 
        hdate_info[hdate_info > (54600-date_offset) - ab_s2sub] -= ab_s3sub
        kdate_info[kdate_info > 54300-date_offset] -= ab_s2sub 
        kdate_info[kdate_info > (54600-date_offset) - ab_s2sub] -= ab_s3sub
    
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

    # Get the grade
    jgrade_info = j_table_info.JGRADE
    hgrade_info = h_table_info.HGRADE
    kgrade_info = k_table_info.KGRADE

    
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
    # Every band will take two steps: putting the errorbars
    
    fmt_info = 'd'
#    fmt_warn = '.'

    # Let's define dictionaries
    d_ax = {'j': ax_j, 'h': ax_h, 'k': ax_k}
    
    if timecolor:
        d_cmap = {'j': time_cmap, 'h': time_cmap, 'k': time_cmap}
    elif type(d_cmap) is str:
        d_cmap = {'j': d_cmap, 'h': d_cmap, 'k': d_cmap}
    elif type(d_cmap) is not dict:
        d_cmap = {'j': d_cmap[0], 'h': d_cmap[1], 'k': d_cmap[2]}
    d_fmt = {'j': 'bo', 'h': 'go', 'k': 'ro'}
    d_fmt_info = {'j': 'b'+fmt_info, 'h':'g'+fmt_info, 'k': 'o'+fmt_info}

    d_date = {'j': jdate, 'h': hdate, 'k': kdate}
    d_col = {'j': jcol, 'h': hcol, 'k': kcol}
    d_err = {'j': jerr, 'h': herr, 'k': kerr}
    d_grade = {'j': jgrade, 'h': hgrade, 'k': kgrade}

    d_date_info = {'j': jdate_info, 'h': hdate_info, 'k': kdate_info}
    d_col_info = {'j': jcol_info, 'h': hcol_info, 'k': kcol_info}
    d_err_info = {'j': jerr_info, 'h': herr_info, 'k': kerr_info}
    d_grade_info = {'j': jgrade_info, 'h': hgrade_info, 'k': kgrade_info}

    d_rawdate = {'j': raw_jdate, 'h': raw_hdate, 'k': raw_kdate}
    d_rawdate_info = {
        'j': raw_jdate_info, 'h': raw_hdate_info, 'k': raw_kdate_info}

    color_vmin = s_table.MEANMJDOBS.min() - date_offset
    color_vmax = s_table.MEANMJDOBS.max() - date_offset
    # min(
    #     (d_rawdate['j'].min(), d_rawdate['h'].min(), d_rawdate['k'].min(),
    #      d_rawdate_info['j'].min(), d_rawdate_info['h'].min(), d_rawdate_info['k'].min()))
    # color_vmax = max(
    #     (d_rawdate['j'].max(), d_rawdate['h'].max(), d_rawdate['k'].max(),
    #      d_rawdate_info['j'].max(), d_rawdate_info['h'].max(), d_rawdate_info['k'].max()))

    if timecolor:
        d_c = d_rawdate
        d_c_info = d_rawdate_info
        vmin = color_vmin
        vmax = color_vmax
    else:
        d_c = d_grade
        d_c_info = d_grade_info
        vmin = 0.8
        vmax = 1

    for band in ['j', 'h', 'k']:

        # Plot a generic band and reduce the size of the code!

        if len(d_date[band]) > 0:
            # First, plot the errorbars, with no markers, in the background:

            d_ax[band].errorbar( d_date[band], d_col[band], marker=None,
                                 yerr=d_err[band], fmt=None, ecolor='k',
                                 zorder=0)
            
            # Next, scatter the points themselves, colored re:grade :
            d_ax[band].scatter( d_date[band], d_col[band], cmap=d_cmap[band],
                                c=d_c[band], vmin=vmin, vmax=vmax, zorder=100)

            

        if len(d_date_info[band]) > 0:
            # First, plot the errorbars, with no markers, in the background:
            d_ax[band].errorbar( d_date_info[band], d_col_info[band], 
                                 yerr=d_err_info[band], marker=None,
                                 fmt=None, ecolor='k', zorder=0)

            # Next, scatter the points themselves, colored re:grade :
            d_ax[band].scatter( d_date_info[band], d_col_info[band], 
                                marker=fmt_info, 
                                c=d_c_info[band], cmap=d_cmap[band], 
                                vmin=vmin, vmax=vmax, zorder=100)

        # Finally, flip it (magnitudes are backwards).
        d_ax[band].invert_yaxis()

        # And plot the dotted lines, if relevant.
        if abridged:
            d_ax[band].plot([ab_s1s2line, ab_s1s2line], [0,30], "k--",
                            scaley=False, scalex=False)

            d_ax[band].plot([ab_s2s3line, ab_s2s3line], [0,30], "k--",
                            scaley=False, scalex=False)
            


    ## Now let's do the 2 color-mag/color-color plots.

    # We'll use different data-cuts for the two different plots.
    # Relevant comment: I made an executive call to include only
    # 'normal' and 'info'-flagged data in the C-C and C-M plots
    # (i.e. max_flag=256 in all relevant bands).
    
    # In the color-mag plot, we need data where H and K are defined 
    # everywhere. That's two cuts.
    khk_table = band_cut( band_cut(s_table, 'k', max_flag=256),
                          'h', max_flag=256)

    khkdate = khk_table.MEANMJDOBS - date_offset
    k_khk = khk_table.KAPERMAG3
    hmk_khk = khk_table.HMKPNT
    k_khk_err = khk_table.KAPERMAG3ERR
    hmk_khk_err = khk_table.HMKPNTERR

    # In the color-color plot, we need data where J, H, and K are
    # defined everywhere. That's one more cut.
    jhk_table = band_cut(khk_table, 'j', max_flag=256)

    jhkdate = jhk_table.MEANMJDOBS - date_offset
    jmh_jhk = jhk_table.JMHPNT
    hmk_jhk = jhk_table.HMKPNT
    jmh_jhk_err = jhk_table.JMHPNTERR
    hmk_jhk_err = jhk_table.HMKPNTERR

    # Plot J-H vs H-K using the "jhk_" variables.
    try:
        plot_trajectory_core( ax_jhk, hmk_jhk, jmh_jhk, jhkdate,
                              vmin=color_vmin, vmax=color_vmax) 

        if color_slope:
            jhk_slope, jhk_intercept, slope_err = (
                slope(hmk_jhk, jmh_jhk, hmk_jhk_err, jmh_jhk_err,
                      verbose=False) )
            
            ax_jhk.plot([0, 6], [jhk_intercept, jhk_intercept + 6*jhk_slope], 
                        ':', scalex=False, scaley=False)
            
    except Exception:
        print "JHK plot broke!"
        pass
        

    # Plot K vs H-K using the "khk_" variables.
    try:
        plot_trajectory_core( ax_khk, hmk_khk, k_khk, khkdate,
                              ms=False, ctts=False, 
                              vmin=color_vmin, vmax=color_vmax) 

        # plot boundaries are manually set for readability, if necessary
        if len(ax_khk.get_xticks()) > 7:
            khk_xmin = np.floor(hmk_khk.min() * 0.95 * 20)/20.
            khk_xmax = np.ceil( hmk_khk.max() * 1.05 * 20)/20.

            khk_xticks = np.linspace(khk_xmin, khk_xmax, 6)
            ax_khk.set_xticks(khk_xticks)

        if color_slope:
            khk_slope, khk_intercept, slope_err = (
                slope(hmk_khk, k_khk, hmk_khk_err, k_khk_err,
                      verbose=False) )
            
            ax_khk.plot([0, 6], [khk_intercept, khk_intercept + 6*khk_slope],
                        '--', scalex=False, scaley=False)
    
    except Exception:
        print "KHK plot broke!"
        pass
    ax_khk.invert_yaxis()

    # Hide the bad labels...
    plt.setp(ax_j.get_xticklabels(), visible=False)
    plt.setp(ax_h.get_xticklabels(), visible=False)

    # Label stuff
#    ax_k.set_xlabel( "Time (JD since 01/01/2000)" )
    ax_k.set_xlabel( "Time (MJD - %.1f)" % date_offset )

    ax_j.set_ylabel( "J",{'rotation':'horizontal', 'fontsize':'large'} )
    ax_h.set_ylabel( "H",{'rotation':'horizontal', 'fontsize':'large'} )
    ax_k.set_ylabel( "K",{'rotation':'horizontal', 'fontsize':'large'} )

    ax_jhk.set_xlabel( "H-K" )
    ax_jhk.set_ylabel( "J-H")#, {'rotation':'horizontal'})
    ax_khk.set_xlabel( "H-K" )
    ax_khk.set_ylabel( "K")#, {'rotation':'horizontal'})

    # Mess around with the X axis if it's abridged:
    if abridged:
        ax_k.set_xticks(ab_xticks)
        ax_k.set_xticklabels(ab_xticklab)
        ax_k.set_xlim(ab_xlim)

    if name != '':
        ax_j.set_title(name)
    else:
        ax_j.set_title(str(sid))

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

    fig.ax_k = ax_k
    fig.ax_h = ax_h
    fig.ax_j = ax_j
    fig.ax_jhk = ax_jhk
    fig.ax_khk = ax_khk

    return fig


def graded_phase (table, sid, period='auto', season=0, offset=0, 
                  outfile='', name='', stetson=True, png_too=False,
                  timecolor=False, time_cmap='jet',
                  d_cmap={'j':'Blues', 'h': 'Greens', 'k': 'Reds'},
                  date_offset = 54034, color_slope=False):
                  
    """ 
    Plots folded lightcurves of a star, with datapoints colored by grade.

    Also plots color-color and color-mag trajectories, colored by phase.

    Will display "lonely" datapoints (i.e. not all JHK mags are 
    well-defined), and plots error-flagged data as different symbols.
    Compare to plot2.phase() which does neither, and to plot3.phase() which 
    is unaware of grades.

    If no period is provided, phase() will compute its own guess 
    using the K-band lightcurve.

    Parameters
    ----------
    table : atpy.Table
        Table with time-series photometry and grade columns. 
        "JGRADE", "HGRADE", "KGRADE" must be bestowed by 
        night_cleanser.null_cleanser_grader().
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
    timecolor : {False, 'phase', 'time'}, optional (default: False)
        Color lightcurve datapoints by phase or time?
        Attempts to sync all colored plots.
        If `timecolor` evaluates to True and isn't one of the above, 
        defaults to 'phase'.
    time_cmap : str, optional (defalt: 'jet')
        Which colormap to use for `timecolor`; when `timecolor`=True,
        this overrides d_cmap.
        Future possibility: the KHK and JHK plots using this cmap.
    d_cmap : dict or str or tuple
        Which colormaps to use for J, H, and K. If a single string 
        (rather than a dict) is given, all 3 bands will use the 
        same colormap. If a tuple is given, J:0, H:1, K:2.
        Default {'j':'Blues', 'h': 'Greens', 'k': 'Reds'} for now.
    date_offset : float, optional
        What MJD to use as day "zero". Default MJD = 54034 
        (initial observations of Orion dataset).
    color_slope : bool, optional (defalt: False)
        Whether to fit color slope lines to the KvH-K and J-HvH-K plots.

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

    ### Initializing variables. Most of the following are done twice:
    # Once with no errors (normal),
    # once with small errors (info).

    ## First: normal

    # Use band_cut to get relevant data chunks.

    j_table = band_cut(s_table, 'j', max_flag=0)
    h_table = band_cut(s_table, 'h', max_flag=0)
    k_table = band_cut(s_table, 'k', max_flag=0)


    # The Julian date for CE  2000 January  1 00:00:00.0 UT is
    # JD 2451544.500000
    # (i.e. MJD 51544.0)

    # get a date (x-axis) for each plot
    jdate = j_table.MEANMJDOBS - date_offset
    hdate = h_table.MEANMJDOBS - date_offset
    kdate = k_table.MEANMJDOBS - date_offset
    
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

    # Get the grade
    jgrade = j_table.JGRADE
    hgrade = h_table.HGRADE
    kgrade = k_table.KGRADE

    ## Second: info

    # Use band_cut to get relevant data chunks.

    j_table_info = band_cut(s_table, 'j', min_flag=1, max_flag=256)
    h_table_info = band_cut(s_table, 'h', min_flag=1, max_flag=256)
    k_table_info = band_cut(s_table, 'k', min_flag=1, max_flag=256)


    # The Julian date for CE  2000 January  1 00:00:00.0 UT is
    # JD 2451544.500000
    # (i.e. MJD 51544.0)

    # get a date (x-axis) for each plot
    jdate_info = j_table_info.MEANMJDOBS - date_offset
    hdate_info = h_table_info.MEANMJDOBS - date_offset
    kdate_info = k_table_info.MEANMJDOBS - date_offset
    
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

    # Get the grade
    jgrade_info = j_table_info.JGRADE
    hgrade_info = h_table_info.HGRADE
    kgrade_info = k_table_info.KGRADE

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
#        print period_string
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
    # Every band will take two steps: putting the errorbars
    
    fmt_info = 'd'
#    fmt_warn = '.'

    # Let's define dictionaries
    d_ax = {'j': ax_j, 'h': ax_h, 'k': ax_k}
    
    if timecolor:
        d_cmap = {'j': time_cmap, 'h': time_cmap, 'k': time_cmap}
    elif type(d_cmap) is str:
        d_cmap = {'j': d_cmap, 'h': d_cmap, 'k': d_cmap}
    elif type(d_cmap) is not dict:
        d_cmap = {'j': d_cmap[0], 'h': d_cmap[1], 'k': d_cmap[2]}
    d_fmt = {'j': 'bo', 'h': 'go', 'k': 'ro'}
    d_fmt_info = {'j': 'b'+fmt_info, 'h':'g'+fmt_info, 'k': 'o'+fmt_info}

    d_date = {'j': jdate, 'h': hdate, 'k': kdate}
    d_col = {'j': jcol, 'h': hcol, 'k': kcol}
    d_err = {'j': jerr, 'h': herr, 'k': kerr}
    d_grade = {'j': jgrade, 'h': hgrade, 'k': kgrade}

    d_date_info = {'j': jdate_info, 'h': hdate_info, 'k': kdate_info}
    d_col_info = {'j': jcol_info, 'h': hcol_info, 'k': kcol_info}
    d_err_info = {'j': jerr_info, 'h': herr_info, 'k': kerr_info}
    d_grade_info = {'j': jgrade_info, 'h': hgrade_info, 'k': kgrade_info}

    # Do we want to have color-color scatter go with time, or phase?
    if timecolor == 'time':
        color_vmin = s_table.MEANMJDOBS.min() - date_offset
        color_vmax = s_table.MEANMJDOBS.max() - date_offset
    else:
        color_vmin = 0
        color_vmax = 1

    # How do we want to color the datapoints themselves?
    if timecolor == 'time':
        d_c = d_date
        d_c_info = d_date_info
        vmin = color_vmin
        vmax = color_vmax
    elif timecolor:
        # This transforms the dictionary of date arrays to
        # a dictionary of phase arrays: it's a "dict comprehension"
        d_c = {key:(((t % period) / period + offset) % 1.) for 
               key,t in d_date.iteritems()}
        d_c_info = {key:(((t % period) / period + offset) % 1.) for 
               key,t in d_date_info.iteritems()}
        vmin = color_vmin
        vmax = color_vmax
    else:
        d_c = d_grade
        d_c_info = d_grade_info
        vmin = 0.8
        vmax = 1

    for band in ['j', 'h', 'k']:

        # Plot a generic band and reduce the size of the code!

        if len(d_date[band]) > 0:
            scatter_phase_core( d_ax[band], d_date[band], d_col[band],
                                d_err[band], period, c=d_c[band],
                                offset = offset,
                                cmap=d_cmap[band], vmin=vmin, vmax=vmax )

            # First, plot the errorbars, with no markers, in the background:
#            d_ax[band].errorbar( d_date[band], d_col[band], marker=None,
#                                 yerr=d_err[band], fmt=None, ecolor='k',
#                                 zorder=0)
            
            # Next, scatter the points themselves, colored re:grade :
#            d_ax[band].scatter( d_date[band], d_col[band], cmap=d_cmap[band],
#                                c=d_grade[band], vmin=0.8, vmax=1, zorder=100)
            

        if len(d_date_info[band]) > 0:

            scatter_phase_core( d_ax[band], d_date_info[band], d_col_info[band],
                                d_err_info[band], period, c=d_c_info[band],
                                cmap=d_cmap[band], vmin=vmin, vmax=vmax, 
                                offset = offset, marker=fmt_info )

            # First, plot the errorbars, with no markers, in the background:
#            d_ax[band].errorbar( d_date_info[band], d_col_info[band], 
#                                 yerr=d_err_info[band], marker=None,
#                                 fmt=None, ecolor='k', zorder=0)

            # Next, scatter the points themselves, colored re:grade :
#            d_ax[band].scatter( d_date_info[band], d_col_info[band], 
#                                marker=fmt_info, s=4,
#                                c=d_grade_info[band], cmap=d_cmap[band], 
#                                vmin=0.8, vmax=1, zorder=100)

        # Finally, flip it (magnitudes are backwards).
        d_ax[band].invert_yaxis()

    ## Now let's do the 2 color-mag/color-color plots.

    # We'll use different data-cuts for the two different plots.
    # Relevant comment: I made an executive call to include only
    # 'normal' and 'info'-flagged data in the C-C and C-M plots
    # (i.e. max_flag=256 in all relevant bands).
    
    # In the color-mag plot, we need data where H and K are defined 
    # everywhere. That's two cuts.
    khk_table = band_cut( band_cut(s_table, 'k', max_flag=256),
                          'h', max_flag=256)

    khkdate = khk_table.MEANMJDOBS - date_offset
    k_khk = khk_table.KAPERMAG3
    hmk_khk = khk_table.HMKPNT
    k_khk_err = khk_table.KAPERMAG3ERR
    hmk_khk_err = khk_table.HMKPNTERR

    khkphase = ((khkdate % period) / period + offset) % 1.

    # In the color-color plot, we need data where J, H, and K are
    # defined everywhere. That's one more cut.
    jhk_table = band_cut(khk_table, 'j', max_flag=256)

    jhkdate = jhk_table.MEANMJDOBS - date_offset
    jmh_jhk = jhk_table.JMHPNT
    hmk_jhk = jhk_table.HMKPNT
    jmh_jhk_err = jhk_table.JMHPNTERR
    hmk_jhk_err = jhk_table.HMKPNTERR

    jhkphase = ((jhkdate % period) / period + offset) % 1.

    if timecolor != 'time':
        color_label = 'Phase'
        jhktime = jhkphase
        khktime = khkphase
    else:
        color_label = 'Time'
        jhktime = jhkdate
        khktime = khkdate

    # Plot J-H vs H-K using the "jhk_" variables.
    try:
        plot_trajectory_core( ax_jhk, hmk_jhk, jmh_jhk, jhktime,
                              vmin=color_vmin, vmax=color_vmax,
                              label = color_label)

        if color_slope:
            jhk_slope, jhk_intercept, slope_err = (
                slope(hmk_jhk, jmh_jhk, hmk_jhk_err, jmh_jhk_err,
                      verbose=False) )
            
            ax_jhk.plot([0, 6], [jhk_intercept, jhk_intercept + 6*jhk_slope], 
                        ':', scalex=False, scaley=False)

    except Exception:
        print "JHK plot broke!"
        pass

    # Plot K vs H-K using the "khk_" variables.
    try:
        plot_trajectory_core( ax_khk, hmk_khk, k_khk, khktime, 
                              ms=False, ctts=False,
                              vmin=color_vmin, vmax=color_vmax,
                              label = color_label)
 

        # plot boundaries are manually set for readability, if necessary
        if len(ax_khk.get_xticks()) > 7:
            khk_xmin = np.floor(hmk_khk.min() * 0.95 * 20)/20.
            khk_xmax = np.ceil( hmk_khk.max() * 1.05 * 20)/20.

            khk_xticks = np.linspace(khk_xmin, khk_xmax, 6)
            ax_khk.set_xticks(khk_xticks)

        if color_slope:
            khk_slope, khk_intercept, slope_err = (
                slope(hmk_khk, k_khk, hmk_khk_err, k_khk_err,
                      verbose=False) )
            
            ax_khk.plot([0, 6], [khk_intercept, khk_intercept + 6*khk_slope],
                        '--', scalex=False, scaley=False)

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
    else:
        ax_j.set_title(str(sid))

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

    fig.ax_k = ax_k
    fig.ax_h = ax_h
    fig.ax_j = ax_j
    fig.ax_jhk = ax_jhk
    fig.ax_khk = ax_khk

    return fig
