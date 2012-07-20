"""
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

import atpy
import numpy as np
import matplotlib.pyplot as plt
#import coords
#import stetson
from chi2 import test_analyze
from scargle import fasper as lsp
from timing import lsp_mask
from tr_helpers import season_cut, data_cut, ensemble_cut


def reddening_vector(x, y, a_v, 
                     shape='full', lw=2, length_includes_head=True,
                     head_width=.05,**kwargs):
    """ Plots a reddening vector in a JHK color-color diagram
    at (x,y) of A_V length.
    
    Returns the reddening vector as a Patch object.
    """
    
    # Extinction for J-H and H-K: E(J-H), E(H-K)
    ejh = 0.107
    ehk = 0.063

    dx = a_v * ehk
    dy = a_v * ejh

    # I don't want the Arrow patch.
    arr = plt.arrow(x, y, dx, dy, color='k',
                    shape=shape, lw=lw, 
                    length_includes_head=length_includes_head, 
                    head_width=head_width, **kwargs)

    return arr


# Let's define Koornneef's main-sequence colors.
ms_jmh = np.array(
    [-0.16, -0.14, -0.13, -0.12, -0.11, -0.1 , -0.09, -0.08, -0.07,
      -0.06, -0.05, -0.03, -0.03, -0.01,  0.01,  0.02,  0.03,  0.04,
      0.05,  0.07,  0.08,  0.1 ,  0.11,  0.13,  0.16,  0.19,  0.2 ,
      0.23,  0.24,  0.29,  0.32,  0.37,  0.43,  0.49,  0.53,  0.57,
      0.61,  0.65,  0.67,  0.68,  0.66,  0.62,  0.61,  0.58,  0.58,  0.57] )
ms_hmk = np.array(
    [-0.05, -0.05, -0.05, -0.05, -0.04, -0.04, -0.04, -0.03, -0.03,
      -0.02, -0.02, -0.02, -0.01, -0.01,  0.  ,  0.  ,  0.01,  0.01,
      0.02,  0.02,  0.02,  0.02,  0.03,  0.04,  0.04,  0.05,  0.06,
      0.06,  0.07,  0.08,  0.09,  0.1 ,  0.11,  0.13,  0.14,  0.15,
      0.16,  0.18,  0.19,  0.21,  0.26,  0.28,  0.29,  0.3 ,  0.31,  0.33] )

# And Meyer's CTTS locus
tt_hmk = np.array([ 0.2,  1. ])
tt_jmh = np.array([ 0.636,  1.1  ])

def plot_trajectory_vanilla (ax, a_k=1, ctts=True):
    """ Does nothing but plots the main sequence, reddening lines.
    Will be extended to also plot the CTTS locus."""

    global ms_hmk, ms_jmh, tt_jmh, tt_hmk
    ax.plot(ms_hmk, ms_jmh, 'k')
    if ctts:
        ax.plot(tt_hmk, tt_jmh, 'k', linewidth=1.5)

    slope = (3.6/2.1)
    top = np.where( (ms_jmh-slope*ms_hmk) == (ms_jmh - slope*ms_hmk).max())
    # Dotted reddening lines: 1. from the bottom
    ax.plot( [ms_hmk[0], ms_hmk[0] + 1.5*a_k], 
             [ms_jmh[0], ms_jmh[0] + 1.5*slope*a_k], 'k--')

    # 2. from the peak
    ax.plot( [ms_hmk[top], ms_hmk[top] + a_k], 
             [ms_jmh[top], ms_jmh[top] + slope*a_k], 'k--')

    # 3. From the end of the CTTS locus
    if ctts:
        ax.plot( [tt_hmk[1], tt_hmk[1] + a_k], 
                 [tt_jmh[1], tt_jmh[1] + slope*a_k], 'k--')

    return

    

# Change the colorbar call to the figure method (figure it out) [ ]
def plot_trajectory_core (ax, hmk, jmh, c, cmap='jet', label="Time",
                          fmt='k.', ms_hmk=ms_hmk, ms_jmh=ms_jmh, a_k=1,
                          ms=True, ctts=True):
    """ Plots the trajectory of some star in color-color space.
    
    Inputs:
      ax -- a matplotlib Axes object (like a canvas) to draw on
      hmk -- a color index to plot on the X-axis (e.g. H-K)
      jmh -- a color index to plot on the Y-axis (e.g. J-H)
      
    Optional inputs:
      fmt -- a matplotlib plot style. Defaults to black dots.
      ms_hmk -- an array of main sequence colors (e.g. H-K)
      ms_jmh -- an array of main sequence colors (e.g. J-H)
      a_k -- the reddening vector to plot parallel lines for
      ms -- Plot the Main Sequence lines? (defaults to True)
      
    Now with color maps! 
    """

    slope = (3.6/2.1)
    # First, plot the background main-sequence stuff
    if ms:
        ax.plot(ms_hmk, ms_jmh, 'k')


        top = np.where( (ms_jmh-slope*ms_hmk) == (ms_jmh - slope*ms_hmk).max())
    # Dotted reddening lines: 1. from the bottom
        ax.plot( [ms_hmk[0], ms_hmk[0] + 1.5*a_k], 
                 [ms_jmh[0], ms_jmh[0] + 1.5*slope*a_k], 'k--')

    # 2. from the peak
        ax.plot( [ms_hmk[top], ms_hmk[top] + a_k], 
                 [ms_jmh[top], ms_jmh[top] + slope*a_k], 'k--')


    if ctts:
        ax.plot(tt_hmk, tt_jmh, 'k', linewidth=1.5)

    # 3. from the CTTS locus
        ax.plot( [tt_hmk[1], tt_hmk[1] + a_k], 
                 [tt_jmh[1], tt_jmh[1] + slope*a_k], 'k--')


    # Then, plot the actual data we were given
    sc = ax.scatter(hmk, jmh, c=c, marker='o', s=10, 
                    cmap=cmap, edgecolors='none')
    cbar = plt.gcf().colorbar(sc, ax=ax) # This should really be changed to the method
    cbar.set_label(label)
#    ax.plot(hmk, jmh, fmt)
    
    # We may want to consciously scale the viewable limits in a specific way
    # YES THIS MUST HAPPEN
    # I propose: bottom-left is fixed, top-right is determined by the input
    # data and then extend the reddening vectors way out.


    # And label the axes! Wait... I'll leave that to the end-user.

    return

#untested
def plot_trajectory (table, sid, season=123, clear=True, fmt='k.'):
    """ Takes a source from a table and plots its color-color trajectory.

    Inputs:
      table -- ATpy time-series photometry table.
      sid -- WFCAM source ID.
      
    Optional inputs: season -- Season 1,2,3 or all.  clear -- enter
      True to make a new figure when this function calls.  fmt -- a
      matplotlib plot style. Defaults to black dots.

    This is a convenience function that tries to be smart.
      """
    if clear:
        plt.figure()
    ax = plt.gca()
    
    tcut = data_cut(table, [sid], season)
    jmh = tcut.JMHPNT
    hmk = tcut.HMKPNT
    date= tcut.MEANMJDOBS - 54579

    plot_trajectory_core (ax, hmk, jmh, date,  fmt=fmt)
    plt.show()
    return

# 'hide' parameter is hackish.
def plot_phase_core (ax, t, x, xerr, period, offset=0, color='k', hide=False):
    """ Plots a pretty period-folded lightcurve on a given axes object.

    Doesn't assume anything about your data (e.g. that it's in magnitudes)
    """
    # Untested.
    
    phase = ((t % period) / period + offset) % 1.


    if not hide:    ax.errorbar(phase, x, yerr=xerr, fmt= color+'o')
    ax.errorbar(phase-1,x,yerr=xerr,fmt='o', mfc='0.7',mec='0.7', 
                 ecolor='0.7')
    ax.errorbar(phase+1,x,yerr=xerr,fmt='o', mfc='0.7',mec='0.7', 
                 ecolor='0.7')
    
    ax.set_xticks( [0, 0.5, 1] )
    ax.set_xticks( np.arange(-.5,1.5,.1), minor=True)

    ax.set_xlim(-0.25, 1.25)

    return period



def lc (table, sid, outfile='', name='?', season=123, png_too=False, 
        flags=0):
    """ 
    Plots J, H, K lightcurves, as well as JHK color-color and color-mag
    trajectories, for one star.

    Inputs:
      table -- atpy table with time series photometry
      sid -- WFCAM source ID of star to plot
      
    Optional inputs:
      outfile -- a place to save the file 
      name -- a short string to display as a name atop the plot
      season -- the usual
      png_too -- if True, saves both a PDF and a PNG of the outfile
                 (note - iff True, don't give a filename extension)
      flags -- whether to remove bad observations from plotting, 
               and where to draw the cutoff.
                 """
    
    # Loading up the relevant datapoints to plot 
    # (note I set 'flags' as a keyword)
    s_table = season_cut(table, sid, season, flags=flags)

    if len(s_table) == 0:
        print "no data here"
        return

    date = s_table.MEANMJDOBS - 54579

    jcol = s_table.JAPERMAG3
    hcol = s_table.HAPERMAG3
    kcol = s_table.KAPERMAG3
    jmh =  s_table.JMHPNT
    hmk =  s_table.HMKPNT

    jerr = s_table.JAPERMAG3ERR
    herr = s_table.HAPERMAG3ERR
    kerr = s_table.KAPERMAG3ERR
    jmherr=s_table.JMHPNTERR
    hmkerr=s_table.HMKPNTERR

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

    # Plot J-band:
    ax_j.errorbar( date, jcol, yerr=jerr, fmt='bo', ecolor='k')
    ax_j.invert_yaxis()

    # Plot H-band:
    ax_h.errorbar( date, hcol, yerr=herr, fmt='go', ecolor='k' )
    ax_h.invert_yaxis()

    # Plot K-band:
    ax_k.errorbar( date, kcol, yerr=kerr, fmt='ro', ecolor='k' )
    ax_k.invert_yaxis()

    # Plot J-H vs H-K
    plot_trajectory_core( ax_jhk, hmk, jmh, date )

    # Plot K vs H-K
    plot_trajectory_core( ax_khk, hmk, kcol, date , ms=False, ctts=False) # gonna update this so that it properly uses K-band main sequence line
    ax_khk.invert_yaxis()

    # Hide the bad labels...
    plt.setp(ax_j.get_xticklabels(), visible=False)
    plt.setp(ax_h.get_xticklabels(), visible=False)

    # Label stuff
    ax_k.set_xlabel( "Time (JD since 04/23/2008)" )

    ax_j.set_ylabel( "J",{'rotation':'horizontal', 'fontsize':'large'} )
    ax_h.set_ylabel( "H",{'rotation':'horizontal', 'fontsize':'large'} )
    ax_k.set_ylabel( "K",{'rotation':'horizontal', 'fontsize':'large'} )

    ax_jhk.set_xlabel( "H-K" )
    ax_jhk.set_ylabel( "J-H")#, {'rotation':'horizontal'})
    ax_khk.set_xlabel( "H-K" )
    ax_khk.set_ylabel( "K")#, {'rotation':'horizontal'})


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


def phase (table, sid, period='auto', outfile='', season=123, offset=0, 
           flags=0, png_too=False):
    """ 
    Plots J, H, K lightcurves, as well as JHK color-color and color-mag
    trajectories, for one star.

    Inputs:
      table -- atpy table with time series photometry
      sid -- WFCAM source ID of star to plot
      
    Optional inputs:
      outfile -- a place to save the file 
      name -- a short string to display as a name atop the plot
      season -- the usual
      png_too -- if True, saves both a PDF and a PNG of the outfile
                 (note - iff True, don't give a filename extension)
      flags -- whether to remove bad observations from plotting, 
               and where to draw the cutoff.

                 """
    
    # Loading up the relevant datapoints to plot 
    # (note I set 'flags' as a keyword)
    s_table = season_cut(table, sid, season, flags=flags)

    if len(s_table) == 0:
        print "no data here"
        return

    date = s_table.MEANMJDOBS - 54579

    jcol = s_table.JAPERMAG3
    hcol = s_table.HAPERMAG3
    kcol = s_table.KAPERMAG3
    jmh =  s_table.JMHPNT
    hmk =  s_table.HMKPNT

    jerr = s_table.JAPERMAG3ERR
    herr = s_table.HAPERMAG3ERR
    kerr = s_table.KAPERMAG3ERR
    jmherr=s_table.JMHPNTERR
    hmkerr=s_table.HMKPNTERR

# Let's figure out the period.
    if period == 'auto':
        period = 1./test_analyze(date, jcol, jerr)
        print period
    elif period == 'lsp':
        lomb = lsp(date,jcol,6.,6.)
        lsp_freq = lomb[0]
        lsp_power= lomb[1]
        Jmax = lsp_mask( lsp_freq, lsp_power)
        lsp_per = 1./ lomb[0][Jmax]
        period = lsp_per
        print period
        

    if period < 1:
        period_string = "%f hours" % (period*24)
        print period_string
    else:
        period_string = "%f days" % period

    phase = ((date % period) / period + offset) % 1.

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

    # Plot J-band:
    plot_phase_core( ax_j, date, jcol, jerr, period, offset=offset, color='b')
#    ax_j.errorbar( date, jcol, yerr=jerr, fmt='bo', ecolor='k')
    ax_j.invert_yaxis()

    # Plot H-band:
    plot_phase_core( ax_h, date, hcol, herr, period, offset=offset, color='g')
#    ax_h.errorbar( date, hcol, yerr=herr, fmt='go', ecolor='k' )
    ax_h.invert_yaxis()

    # Plot K-band:
    plot_phase_core( ax_k, date, kcol, kerr, period, offset=offset, color='r')
#    ax_k.errorbar( date, kcol, yerr=kerr, fmt='ro', ecolor='k' )
    ax_k.invert_yaxis()

    # Plot J-H vs H-K
    plot_trajectory_core( ax_jhk, hmk, jmh, phase , label='Phase')

    # Plot K vs H-K
    plot_trajectory_core( ax_khk, hmk, kcol, phase, label='Phase', ms=False, ctts=False) # gonna update this so that it properly uses K-band main sequence line
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

    return period


def lsp_power (table, sid, outfile='', name='', season=123, png_too=False):
    """ 
    Plots J, H, K periodograms for one star.

    Inputs:
      table -- atpy table with time series photometry
      sid -- WFCAM source ID of star to plot
      
    Optional inputs:
      outfile -- a place to save the file 
      name -- a short string to display as a name atop the plot
      season -- the usual
      png_too -- if True, saves both a PDF and a PNG of the outfile
                 (note - iff True, don't give a filename extension)
                 """
    
    # Loading up the relevant datapoints to plot (note I set flags to 0)
    s_table = season_cut(table, sid, season, flags=0)

    if len(s_table) < 2:
        print "no data here"
        return

    date = s_table.MEANMJDOBS - 54579

    jcol = s_table.JAPERMAG3
    hcol = s_table.HAPERMAG3
    kcol = s_table.KAPERMAG3
#    jmh =  s_table.JMHPNT
#    hmk =  s_table.HMKPNT

#    jerr = s_table.JAPERMAG3ERR
#    herr = s_table.HAPERMAG3ERR
#    kerr = s_table.KAPERMAG3ERR
#    jmherr=s_table.JMHPNTERR
#    hmkerr=s_table.HMKPNTERR

    jlsp = lsp(date, jcol, 6., 6.)
    hlsp = lsp(date, hcol, 6., 6.)
    klsp = lsp(date, kcol, 6., 6.)

    j_lsp_freq = jlsp[0]
    h_lsp_freq = hlsp[0]
    k_lsp_freq = klsp[0]

    j_lsp_power = jlsp[1]
    h_lsp_power = hlsp[1]
    k_lsp_power = klsp[1]

    # best periods, filtered by the lsp_mask
    j_lsp_per = 1./ j_lsp_freq[ lsp_mask( j_lsp_freq, j_lsp_power) ]    
    h_lsp_per = 1./ h_lsp_freq[ lsp_mask( h_lsp_freq, h_lsp_power) ]
    k_lsp_per = 1./ k_lsp_freq[ lsp_mask( k_lsp_freq, k_lsp_power) ]
    

    fig = plt.figure(figsize = (10, 6), dpi=80,
                     facecolor='w', edgecolor='k')

    ax_j = fig.add_subplot(3,1,1)
    ax_h = fig.add_subplot(3,1,2, sharex=ax_j)
    ax_k = fig.add_subplot(3,1,3, sharex=ax_j)

    ax_j.plot(1./j_lsp_freq, j_lsp_power, 'b')
    ax_h.plot(1./h_lsp_freq, h_lsp_power, 'g')
    ax_k.plot(1./k_lsp_freq, k_lsp_power, 'r')
    
    ax_j.set_xscale('log')
    ax_h.set_xscale('log')
    ax_k.set_xscale('log')

    ax_j.set_title(name)
    ax_k.set_xlabel("Period (days)")
    ax_h.set_ylabel("Periodogram Power")

    # bottom = 0.1
    # height = .25
    # left = 0.075
    # width = 0.5

    # ax_k = fig.add_axes( (left, bottom, width, height) )
    # ax_h = fig.add_axes( (left, bottom+.3, width, height), sharex=ax_k )
    # ax_j = fig.add_axes( (left, bottom+.6, width, height), sharex=ax_k )
    
    # ax_jhk = fig.add_axes( (.65, bottom, .3, .375) )
    # ax_khk = fig.add_axes( (.65, bottom+.475, .3, .375) )

    # # Plot J-band:
    # ax_j.errorbar( date, jcol, yerr=jerr, fmt='bo', ecolor='k')
    # ax_j.invert_yaxis()

    # # Plot H-band:
    # ax_h.errorbar( date, hcol, yerr=herr, fmt='go', ecolor='k' )
    # ax_h.invert_yaxis()

    # # Plot K-band:
    # ax_k.errorbar( date, kcol, yerr=kerr, fmt='ro', ecolor='k' )
    # ax_k.invert_yaxis()

    # # Plot J-H vs H-K
    # plot_trajectory_core( ax_jhk, hmk, jmh, date )

    # # Plot K vs H-K
    # plot_trajectory_core( ax_khk, hmk, kcol, date , ms=False, ctts=False) # gonna update this so that it properly uses K-band main sequence line
    # ax_khk.invert_yaxis()

    # # Hide the bad labels...
    # plt.setp(ax_j.get_xticklabels(), visible=False)
    # plt.setp(ax_h.get_xticklabels(), visible=False)

    # # Label stuff
    # ax_k.set_xlabel( "Time (JD since 04/23/2008)" )

    # ax_j.set_ylabel( "J",{'rotation':'horizontal', 'fontsize':'large'} )
    # ax_h.set_ylabel( "H",{'rotation':'horizontal', 'fontsize':'large'} )
    # ax_k.set_ylabel( "K",{'rotation':'horizontal', 'fontsize':'large'} )

    # ax_jhk.set_xlabel( "H-K" )
    # ax_jhk.set_ylabel( "J-H")#, {'rotation':'horizontal'})
    # ax_khk.set_xlabel( "H-K" )
    # ax_khk.set_ylabel( "K")#, {'rotation':'horizontal'})


    if outfile == '':
        plt.show()
    else:
        if png_too:
            plt.savefig(outfile+".pdf")
            plt.savefig(outfile+".png")
            plt.close()
        else:
            plt.savefig(outfile)
            plt.close()



def phase_trajectory (table, sid, period='auto', outfile='', season=123, 
                      offset=0):
    """ Does just the trajectory window from earlier. """
    # Loading up the relevant datapoints to plot (note I set flags to 0)
    s_table = season_cut(table, sid, season, flags=0)

    if len(s_table) == 0:
        print "no data here"
        return

    date = s_table.MEANMJDOBS - 54579

    jcol = s_table.JAPERMAG3
    hcol = s_table.HAPERMAG3
    kcol = s_table.KAPERMAG3
    jmh =  s_table.JMHPNT
    hmk =  s_table.HMKPNT

    jerr = s_table.JAPERMAG3ERR
    herr = s_table.HAPERMAG3ERR
    kerr = s_table.KAPERMAG3ERR
    jmherr=s_table.JMHPNTERR
    hmkerr=s_table.HMKPNTERR

# Let's figure out the period.
    if period == 'auto':
        period = 1./test_analyze(date, jcol, jerr)
        print period
    elif period == 'lsp':
        lomb = lsp(date,jcol,6.,6.)
        lsp_freq = lomb[0]
        lsp_power= lomb[1]
        Jmax = lsp_mask( lsp_freq, lsp_power)
        lsp_per = 1./ lomb[0][Jmax]
        period = lsp_per
        print period
        

    if period < 1:
        period_string = "%f hours" % (period*24)
        print period_string
    else:
        period_string = "%f days" % period

    phase = ((date % period) / period + offset) % 1.

    fig = plt.figure(figsize = (10, 6), dpi=80,
                     facecolor='w', edgecolor='k')

    ax_jhk = plt.subplot(1,1,1)
    plot_trajectory_core( ax_jhk, hmk, jmh, phase , label='Phase')
    ax_jhk.set_xlabel( "H-K" )
    ax_jhk.set_ylabel( "J-H")#, {'rotation':'horizontal'})

    return period
