''' Functions to plot lightcurves for WFCAM data.

Publicly available: 
plot_lc (make it just lc someday?)
plot_phase (make it just phase someday)
plot_5
plot_ensemble
'''

import atpy
import numpy as np
import matplotlib.pyplot as plt
import coords
import stetson
from chi2 import test_analyze
from scargle import fasper as lsp
from tr_helpers import season_cut, data_cut, ensemble_cut

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


def plot_trajectory_core (ax, hmk, jmh, 
                          fmt='k.', ms_hmk=ms_hmk, ms_jmh=ms_jmh, a_k=1):
    ''' Plots the trajectory of some star in color-color space.
    
    Inputs:
      ax -- a matplotlib Axes object (like a canvas) to draw on
      hmk -- a color index to plot on the X-axis (e.g. H-K)
      jmh -- a color index to plot on the Y-axis (e.g. J-H)
      
    Optional inputs:
      fmt -- a matplotlib plot style. Defaults to black dots.
      ms_hmk -- an array of main sequence colors (e.g. H-K)
      ms_jmh -- an array of main sequence colors (e.g. J-H)
      a_k -- the reddening vector to plot parallel lines for
    '''
    # First, plot the background main-sequence stuff
    ax.plot(ms_hmk, ms_jmh, 'k')

    slope = (3.6/2.1)
    top = np.where( (ms_jmh - slope*ms_hmk) == (ms_jmh - slope*ms_hmk).max())
    # Dotted reddening lines: 1. from the bottom
    ax.plot( [ms_hmk[0], ms_hmk[0] + 1.5*a_k], 
             [ms_jmh[0], ms_jmh[0] + 1.5*slope*a_k], 'k--')

    # 2. from the peak
    ax.plot( [ms_hmk[top], ms_hmk[top] + a_k], 
             [ms_jmh[top], ms_jmh[top] + slope*a_k], 'k--')


    # Then, plot the actual data we were given
    ax.plot(hmk, jmh, fmt)
    
    # We may want to consciously scale the viewable limits in a specific way

    # And label the axes! Wait... I'll leave that to the end-user.

    return

#untested
def plot_trajectory (table, sid, season=123, clear=True, fmt='k.'):
    ''' Takes a source from a table and plots its color-color trajectory.

    Inputs:
      table -- ATpy time-series photometry table.
      sid -- WFCAM source ID.
      
    Optional inputs:
      season -- Season 1,2,3 or all.
      clear -- enter True to make a new figure when this function calls.
      fmt -- a matplotlib plot style. Defaults to black dots.
      '''
    if clear:
        plt.figure()
    ax = plt.gca()
    
    tcut = data_cut(table, [sid], season)
    jmh = tcut.JMHPNT
    hmk = tcut.HMKPNT

    plot_trajectory_core (ax, hmk, jmh, fmt=fmt)
    plt.show()
    return
                         

def plot_jhmk (table,sid, outfile='', sup='', box=True, 
             bands='jhk', season = 123, text=True) :
    ''' Plots J and H-K lightcurves in 3 seasons (separated) over time.
    '''

    s_table = season_cut(table, sid, 123)
    date = s_table.MEANMJDOBS - 54579
    
    jcol = s_table.JAPERMAG3
    hmkcol=s_table.HMKPNT

    jerr = s_table.JAPERMAG3ERR
    hmkerr=s_table.HMKPNTERR

    plt.clf()

    ax1 = plt.subplot(2,3,1)
    ax2 = plt.subplot(2,3,2, sharey=ax1)
    ax3 = plt.subplot(2,3,3, sharey=ax1)

    sx1 = plt.subplot(2,3,4, sharex=ax1)
    sx2 = plt.subplot(2,3,5, sharex=ax2)
    sx3 = plt.subplot(2,3,6, sharex=ax3)

    for ax in (ax1,ax2,ax3):
        ax
        plt.errorbar(date, jcol, yerr=jerr, fmt='k-o', ecolor='k')

    for sx in (sx1,sx2,sx3):
        sx
        plt.errorbar(date, hmkcol, yerr=hmkerr, fmt='k-o', ecolor='k')

    ax1
    plt.gca().invert_yaxis()
    plt.ylabel("WFCAM J Magnitude")
    
    sx1
    d = 54579
    plt.xlim(d, d+100)
    plt.xlabel("Season 1")
    sx2
    plt.xlim(d+100,d+300)
    plt.xlabel("Season 2")
    sx3
    plt.xlim(d+300,d+600)
    plt.xlabel("Season 3")

    plt.show()
    return True

# junk:
'''

    s_table1 = season_cut(table, sid, 1)
    s_table2 = season_cut(table, sid, 2)
    s_table3 = season_cut(table, sid, 3)

    date1 = s_table1.MEANMJDOBS - 54579
    date2 = s_table2.MEANMJDOBS - 54579
    date3 = s_table3.MEANMJDOBS - 54579

    j1 = s_table1.JAPERMAG3
    jmk1=s_table1.JMKPNT
    j2 = s_table2.JAPERMAG3
    jmk2=s_table2.JMKPNT
    j3 = s_table3.JAPERMAG3
    jmk3=s_table3.JMKPNT

'''
#/end junk

def plot_page_periods (table, sid, outfile='', name='?', season = 123):
    """
    Plot one comprehensive page of periodicity information for one source
    in WFCAM time-series JHK data.

     INPUTS:
             table: An atpy table with WFCAM time-series photometry
             sid: A 13-digit WFCAM source ID

     OPTIONAL INPUTS:
             outfile: an output filename, including path and filetype 
                 extension, to save plot to (rather than plot interactively)
             name: a nickname/designation for your source
             season: which season to plot? 1, 2, 3, or all

     OUTPUTS:
             Returns nothing; optionally saves an output plot.

     PROCEDURE:
             This function calls the Lomb Scargle period-finding algorithm and
             the Palmer fast chi-square minimization period-finding algorithm
             to search for periods in the source, in each band.
             Then it plots the periodogram and the lightcurve folded by the
             two best periods (lomb period and palmer period).
             
    """
    
    # 1. Loading up the relevant datapoints to plot (note I set flags to 0)
    s_table = season_cut(table, sid, season, flags=0)

    if len(s_table) < 2:
        print "no data here"
        return

    date = s_table.MEANMJDOBS - 54579

    class Band:
        pass

    j = Band()
    h = Band()
    k = Band()
    jmh = Band()
    hmk = Band()
    kdex = Band()

    bands = [ j, h, k, jmh, kdex, hmk ]

    j.col = s_table.JAPERMAG3
    h.col = s_table.HAPERMAG3
    k.col = s_table.KAPERMAG3
    jmh.col = s_table.JMHPNT
    hmk.col = s_table.HMKPNT

    j.err = s_table.JAPERMAG3ERR
    h.err = s_table.HAPERMAG3ERR
    k.err = s_table.KAPERMAG3ERR
    jmh.err=s_table.JMHPNTERR
    hmk.err=s_table.HMKPNTERR

    kdex.col = 1.71 * hmk.col - jmh.col
    kdex.err = np.sqrt( j.err**2 + h.err**2 + k.err**2 )
    # Done loading up data.

    # 2. Compute periods and basic statistics
    
    # a. Overall stetson variability 

    stet = stetson.S( j.col, j.err, h.col, h.err, k.col, k.err )


    for b in bands:

        # b. Lomb-scargle periodogram for each band        
        b.lsp = lsp(date, b.col, 6., 6.)
        
        b.lsp_freq = b.lsp[0]
        b.lsp_power= b.lsp[1]
        
        b.lsp_per = 1./ b.lsp[0][b.lsp[3]]

        # c. Fast Chi-squared period for each band

        b.fx2_per = 1./ test_analyze( date, b.col, b.err) #confirmed syntax

    # ok, now as a test let's print these quantities that we just calculated

#    print "Stetson index: " + str(stet)
#    for b, n in zip(bands, ('j','h','k', 'j-h','h-k','kdex')):
#        print n.upper() + " band LSP period: " + str(b.lsp_per)
#        print n.upper() + " band fx2 period: " + str(b.fx2_per)
        
    # I gotta silence the output of chi whatever. This may involve some serious popen whatever shit. (or just removing a print statement somewhere...)

    # 3. Create the canvas
    fig = plt.figure(num=None, figsize=(8.5,11), dpi=80, 
                     facecolor='w', edgecolor='k')
    
    # My first approach: using subplots rather than custom coding

    qs = 3 * np.arange(6)
    
    for b, q in zip(bands, qs) : #qs: something about dimension parameters:
#        b.ax1 = fig.add_axes
        b.ax1 = fig.add_subplot(6,3,q+1)
        b.ax2 = fig.add_subplot(6,3,q+2)
        b.ax3 = fig.add_subplot(6,3,q+3)

        b.ax1.plot(1./b.lsp_freq, b.lsp_power)
        b.ax1.set_xscale('log')

    colors = ('b', 'g', 'r')
    for b, c in zip((j,h,k),colors):

        plot_phase_core(b.ax2, date, b.col, b.err, b.lsp_per, color=c)
        plot_phase_core(b.ax3, date, b.col, b.err, b.fx2_per, color=c)

        b.ax2.invert_yaxis()
        b.ax3.invert_yaxis()

    for b in (jmh, kdex, hmk):
        plot_phase_core(b.ax2, date, b.col, b.err, b.lsp_per)
        plot_phase_core(b.ax3, date, b.col, b.err, b.fx2_per)

        if b is kdex:
            # plot a dotted line:
            xs = [-0.25, 1.25]
            ys = [.1, .1]
            b.ax2.plot( xs, ys, 'r--')
            b.ax3.plot( xs, ys, 'r--')
            # plot red dots on disky nights:
            disk = np.where(kdex.col > .1)

            if date[disk].size > 0:
                plot_phase_core(b.ax2, date[disk], b.col[disk], 
                                b.err[disk], b.lsp_per, color='r')
                plot_phase_core(b.ax3, date[disk], b.col[disk], 
                                b.err[disk], b.fx2_per, color='r')


    jmean = j.col.mean()
    hmean = h.col.mean()
    kmean = k.col.mean()
    jrms, hrms, krms = j.col.std(), h.col.std(), k.col.std()

    sra, sdec = s_table.RA[0], s_table.DEC[0]
    sPosition = coords.Position((sra,sdec), units='rad')
    sPositionString = sPosition.hmsdms()

    # I'm testing an invisible big axes thing for my title
#    bigAxes = plt.axes(frameon=False)
#    plt.xticks([])
#    plt.yticks([])

    big_title= (("Object %s.\t" %name)+ 
                (r"$J_{mean} =$ %.2f, $H_{mean} =$ %.2f, $K_{mean} =$ %.2f, " 
                 % (jmean, hmean, kmean) ) + 
                "\nSeason %d \t" % season +
                r"$J_{RMS} =$ %.3f, $H_{RMS} =$ %.3f, $K_{RMS} =$ %.3f" % 
                (jrms, hrms, krms) )

    mean_per = np.mean( [j.lsp_per, h.lsp_per, k.lsp_per,  
                         j.fx2_per, h.fx2_per, k.fx2_per] )
    
    second_title= ( ("Stetson Index: %.1f \n" % stet) +
                    ("Average Period: %.2f days\n" % mean_per)
                    )

#    plt.title(big_title)

    j.ax1.annotate( big_title, xy=(0.025, 0.965), 
                    xycoords='figure fraction',
                    horizontalalignment='left', verticalalignment='top')
    
    j.ax3.annotate( second_title, xy=(0.925, 0.965),
                    xycoords='figure fraction',
                    horizontalalignment='right', verticalalignment='top')
    
# Use text and bbox to draw the fitted periods
    for b in bands:
        if b.lsp_per > .9:
            b.lsp_per_str = "period: %.2f days" % b.lsp_per
        else:
            b.lsp_per_str = "period: %.2f hours" % (b.lsp_per * 24)
        
        if b.fx2_per > .9:
            b.fx2_per_str = "period: %.2f days" % b.fx2_per
        else:
            b.fx2_per_str = "period: %.2f hours" % (b.fx2_per * 24)


        plt.text(0.05, 0.05, b.lsp_per_str,
             horizontalalignment='left', verticalalignment='bottom', 
             bbox=dict(facecolor='white'),transform= b.ax2.transAxes)

        plt.text(0.05, 0.05, b.fx2_per_str,
             horizontalalignment='left', verticalalignment='bottom', 
             bbox=dict(facecolor='white'),transform= b.ax3.transAxes)


    plt.suptitle( "Position: %s,       Source ID %d." % (sPositionString, sid))

    j.ax1.set_title("Lomb-Scargle Periodogram", fontsize=10)
    j.ax2.set_title("Best LSP period", fontsize=10)
    j.ax3.set_title("Best fX2 period", fontsize=10)
    hmk.ax1.set_xlabel("Period (days)")
    hmk.ax2.set_xlabel("Phase")
    hmk.ax3.set_xlabel("Phase")

    if outfile == '':
        plt.show()
    else:
        plt.savefig(outfile)
        plt.close()
            
    return 


# I'll want to write a function that takes in an axes object and plots
# something specific (like the folded phase plot with the different grayscales
# and whatnot) in a consistent way, since I really like that thing.
# Success! The above is now plot_phase_core.

        
    
# Functionality to add to plot_5:
# 1. color-color Trajectory plots (in their own function, to be implemented
# and incorporated) check1[x] check2[ ]
# 2. a map and chip identifier, about the same size as that guy [ ]
# note: perhaps make colors and k-dex smaller vertically 
# to make room for those guys [ ]
# also: rename k-dex to "k excess" [x]
def plot_5 (table,sid, outfile='', name='?', season = 123) :
    ''' Plots all five lightcurves of one star: J, H, K, J-H, H-K, 
    on one page, for one season.
    '''

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

    kdex = 1.71 * hmk - jmh
    kdexerr = np.sqrt( jerr**2 + herr**2 + kerr**2 )
    # Done loading up data.
    
    # Let's create that plot.
    fig = plt.figure(num=None, figsize=(8.5,11), dpi=80, 
                     facecolor='w', edgecolor='k')
    
    left = .125 # don't touch
    width= .775 # these two.

    height = (.8 / 6)
    bottom = np.arange(.1, .9, height) 
    tri_height = (.8/9)
    tri_bottom = np.arange(.1, .9, tri_height)

    ks = .05 # this is the size of the kdex box

    ax_hmk = fig.add_axes( [left, bottom[0], width, tri_height] )
    ax_kdex= fig.add_axes( [left, tri_bottom[1], width, tri_height] 
                           ,sharex=ax_hmk)
    ax_jmh = fig.add_axes( [left, tri_bottom[2], width, tri_height] 
                           ,sharex=ax_hmk)
    ax_k  =  fig.add_axes( [left, bottom[2], width, height] ,sharex=ax_hmk)
    ax_h  =  fig.add_axes( [left, bottom[3], width, height] ,sharex=ax_hmk)
    ax_j  =  fig.add_axes( [left, bottom[4], width, height] ,sharex=ax_hmk)
    
    ax_traj = fig.add_subplot(6,2,1)
    ax_map =  fig.add_subplot(6,2,2)

    # Plot J-band:
    ax_j.errorbar( date, jcol, yerr=jerr, fmt='b-o', ecolor='k' )
    ax_j.invert_yaxis()
#    ax_j.set_xticklabels([])

    # Plot H-band:
    ax_h.errorbar( date, hcol, yerr=herr, fmt='g-o', ecolor='k' )
    ax_h.invert_yaxis()

    # Plot K-band:
    ax_k.errorbar( date, kcol, yerr=kerr, fmt='r-o', ecolor='k' )
    ax_k.invert_yaxis()

    # Plot J-H color:
    ax_jmh.errorbar( date, jmh, yerr=jmherr, fmt='k-o', ecolor='k' )

    # Plot K dex:
    ax_kdex.errorbar( date, kdex, yerr=kdexerr, fmt= 'k-o', ecolor='k' )
    # plot a dotted line:
    xs = [date.min(), date.max()]
    ys = [.1, .1]
    ax_kdex.plot( xs, ys, 'r--')
    # plot red dots on disky nights:
    disk = np.where(kdex > .1)

    ax_kdex.plot( date[disk], kdex[disk], 'ro' )
    

    # Plot H-K color:
    ax_hmk.errorbar( date, hmk, yerr=hmkerr, fmt='k-o', ecolor='k' )

    # Plot the trajectory thingy!
    plot_trajectory_core( ax_traj, hmk, jmh )
    

    # Done plotting.

    plt.setp(ax_j.get_xticklabels(), visible=False)
    plt.setp(ax_h.get_xticklabels(), visible=False)
    plt.setp(ax_k.get_xticklabels(), visible=False)
    plt.setp(ax_jmh.get_xticklabels(), visible=False)
    plt.setp(ax_kdex.get_xticklabels(), visible=False)

    # Now let's create labelling information!

    ax_hmk.set_xlabel( "Time (JD since 04/23/2008)" )
    ax_j.set_ylabel( "J mag" )
    ax_h.set_ylabel( "H mag" )
    ax_k.set_ylabel( "K mag" )
    ax_jmh.set_ylabel( "J-H color" )
    ax_hmk.set_ylabel( "H-K color" )
    ax_kdex.set_ylabel( "K excess" )

    jmean = jcol.mean()
    hmean = hcol.mean()
    kmean = kcol.mean()
    jrms, hrms, krms = jcol.std(), hcol.std(), kcol.std()

    sra, sdec = s_table.RA[0], s_table.DEC[0]
    sPosition = coords.Position((sra,sdec), units='rad')
    sPositionString = sPosition.hmsdms()
    
    # i am aware that the following line is SUCH A MESS
    ax_j.set_title( ("Object %s.\t" %name)+ (r"$J_{mean} =$ %.2f, $H_{mean} =$ %.2f, $K_{mean} =$ %.2f, " % 
                    (jmean, hmean, kmean) ) + "\nSeason %d \t" % season +
                    r"$J_{RMS} =$ %.3f, $H_{RMS} =$ %.3f, $K_{RMS} =$ %.3f" % (jrms, hrms, krms) )
    
    plt.suptitle( "Position: %s,       Source ID %d." % (sPositionString, sid))

    if outfile == '':
        plt.show()
    else:
        plt.savefig(outfile)
        plt.close()



def plot_5_ensemble ( correction_table, chip, outfile='', season=2 ) :
    ''' Plots the "ensemble" lightcurve, i.e. the correction each night.

    Borrows HEAVILY from plot_5
    '''

#    sid = chip
#    table = correction_table
#    table.SOURCEID = table.chip
#    table.MEANMJDOBS = table.date

    # Loading up the relevant datapoints to plot (note I set flags to 0)
    s_table = ensemble_cut( correction_table, chip, season)

    if len(s_table) == 0:
        print "no data here"
        return

    date = s_table.date - 54579

    jcol = -s_table.j_correction  # This is the section where we deviate
    hcol = -s_table.h_correction  # the most from the original - but my 
    kcol = -s_table.k_correction  # goal is to bootstrap as much as I can
    jmh =  jcol - hcol            # so that we can leave the rest of the
    hmk =  hcol - kcol            # function unchanged! 

    jerr = np.zeros_like(jcol)
    herr = np.zeros_like(jcol)
    kerr = np.zeros_like(jcol)
    jmherr=np.zeros_like(jcol)
    hmkerr=np.zeros_like(jcol)

    kdex = 1.71 * hmk - jmh
    kdexerr = np.sqrt( jerr**2 + herr**2 + kerr**2 )
    # Done loading up data.
    
    # Let's create that plot.
    fig = plt.figure(num=None, figsize=(8.5,11), dpi=80, 
                     facecolor='w', edgecolor='k')
    
    left = .125
    width= .775
    bottom = np.arange(.1, .9, .16)
    height = .16

    ks = .05 # this is the size of the kdex box

    ax_hmk = fig.add_axes( [left, bottom[0], width, height-ks] )
    ax_kdex= fig.add_axes( [left, bottom[1]-ks, width, 2*ks] ,sharex=ax_hmk)

    ax_jmh = fig.add_axes( [left, bottom[1]+ks, width, height-ks] ,sharex=ax_hmk)
    ax_k  =  fig.add_axes( [left, bottom[2], width, height] ,sharex=ax_hmk)
    ax_h  =  fig.add_axes( [left, bottom[3], width, height] ,sharex=ax_hmk)
    ax_j  =  fig.add_axes( [left, bottom[4], width, height] ,sharex=ax_hmk)
    

    # Plot J-band:
    ax_j.errorbar( date, jcol, yerr=jerr, fmt='b-o', ecolor='k' )
    ax_j.invert_yaxis()
#    ax_j.set_xticklabels([])

    # Plot H-band:
    ax_h.errorbar( date, hcol, yerr=herr, fmt='g-o', ecolor='k' )
    ax_h.invert_yaxis()

    # Plot K-band:
    ax_k.errorbar( date, kcol, yerr=kerr, fmt='r-o', ecolor='k' )
    ax_k.invert_yaxis()

    # Plot J-H color:
    ax_jmh.errorbar( date, jmh, yerr=jmherr, fmt='k-o', ecolor='k' )

    # Plot K dex:
    ax_kdex.errorbar( date, kdex, yerr=kdexerr, fmt= 'k-o', ecolor='k' )
    # plot a dotted line:
    xs = [date.min(), date.max()]
    ys = [.1, .1]
    ax_kdex.plot( xs, ys, 'r--')
    # plot red dots on disky nights:
    disk = np.where(kdex > .1)

    ax_kdex.plot( date[disk], kdex[disk], 'ro' )
    

    # Plot H-K color:
    ax_hmk.errorbar( date, hmk, yerr=hmkerr, fmt='k-o', ecolor='k' )

    # Done plotting.

    plt.setp(ax_j.get_xticklabels(), visible=False)
    plt.setp(ax_h.get_xticklabels(), visible=False)
    plt.setp(ax_k.get_xticklabels(), visible=False)
    plt.setp(ax_jmh.get_xticklabels(), visible=False)
    plt.setp(ax_kdex.get_xticklabels(), visible=False)

    # Now let's create labelling information!

    ax_hmk.set_xlabel( "Time (JD since 04/23/2008)" )
    ax_j.set_ylabel( "J mag" )
    ax_h.set_ylabel( "H mag" )
    ax_k.set_ylabel( "K mag" )
    ax_jmh.set_ylabel( "J-H color" )
    ax_hmk.set_ylabel( "H-K color" )
    ax_kdex.set_ylabel( "K excess" )

    jmean = jcol.mean()
    hmean = hcol.mean()
    kmean = kcol.mean()
    jrms, hrms, krms = jcol.std(), hcol.std(), kcol.std()

#    sra, sdec = s_table.RA[0], s_table.DEC[0]
#    sPosition = coords.Position((sra,sdec), units='rad')
#    sPositionString = sPosition.hmsdms()
    
    # i am aware that the following line is SUCH A MESS
    ax_j.set_title( ("Ensemble fluctuations for chip %s." % chip ) +
                    ("\nSeason %d " % season ) )

#    ax_j.set_title( ("Object %s.\t" %name)+ (r"$J_{mean} =$ %.2f, $H_{mean} =$ %.2f, $K_{mean} =$ %.2f, " % 
#                    (jmean, hmean, kmean) ) + "\nSeason %d \t" % season +
#                    r"$J_{RMS} =$ %.3f, $H_{RMS} =$ %.3f, $K_{RMS} =$ %.3f" % (jrms, hrms, krms) )
    
#    plt.suptitle( "Position: %s,       Source ID %d." % (sPositionString, sid))

    if outfile == '':
        plt.show()
    else:
        plt.savefig(outfile)
        plt.close()



def plot_lc (table,sid, outfile='', sup='', box=True, 
             bands='jhk', season = 123, text=True) :
    ''' Plots J,H,K lightcurves WITH ERRORBARS for a given input source.

    Written with WFCAM columns in mind, specifically like from WSERV1.
    '''
    
#    w = numpy.where( table.SOURCEID == sid )

    s_table = season_cut(table, sid, season)

    ra1, ra2 = 314.36, 315.77
    dec1,dec2= 52.02, 52.92

    rabox = [ra1, ra1, ra2, ra2, ra1]
    decbox= [dec1,dec2,dec2,dec1,dec1]

#    print w
#    print table.RA[w]
##    print s_table.RA
#    sra, sdec = table.RA[w][0], table.DEC[w][0]
    sra, sdec = s_table.RA[0], s_table.DEC[0]

    date = s_table.MEANMJDOBS - 54579

    jcol = s_table.JAPERMAG3
    hcol = s_table.HAPERMAG3
    kcol = s_table.KAPERMAG3

    jerr = s_table.JAPERMAG3ERR
    herr = s_table.HAPERMAG3ERR
    kerr = s_table.KAPERMAG3ERR

    plt.clf()

    # This is the size of the lightcurve box
    if box:
        plt.axes([.05,.1,.7,.8])
        
    if 'j' in bands:
        plt.errorbar(date,jcol,yerr=jerr,fmt='b-o',ecolor='k',
                     label="J-band, RMS: %f" % jcol.std())
    if 'h' in bands:
        plt.errorbar(date,hcol,yerr=herr,fmt='g-o',ecolor='k',
                     label="H-band, RMS: %f" % hcol.std())
    if 'k' in bands:
        plt.errorbar(date,kcol,yerr=kerr,fmt='r-o',ecolor='k',
                     label="K-band, RMS: %f" % kcol.std())

    plt.gca().invert_yaxis()

    if text:
        plt.legend()

        '''
        if 'j' in bands:
            plt.text(20,jcol.max()+.1,"J-band RMS: %f" % jcol.std())
        if 'h' in bands:
            plt.text(20,hcol.max()+.1,"H-band RMS: %f" % hcol.std())
        if 'k' in bands:
            plt.text(20,kcol.max()+.1,"K-band RMS: %f" % kcol.std())
            '''

    plt.ylabel("WFCAM magnitude")
    plt.xlabel("Julian days since 04/23/2008")
    plt.title("J, H, K with errorbars. Source ID %d." % sid)
    plt.suptitle(sup)

    # This is the size of the position box
    if box:
        plt.axes([.775,.55,.2,.35])
        plt.plot(rabox, decbox)
        plt.plot(np.degrees(sra),np.degrees(sdec),'rD')
    
        dx = 1/8. * (ra2-ra1)
        dy = 1/8. * (dec2-dec1)
        arrx = ra1 + 2*dx
        arry = dec1+6*dy
        plt.arrow(arrx,arry,dx,0)
        plt.arrow(arrx,arry,0,dy)
        
        plt.gca().invert_xaxis()
        
        plt.ylabel("Dec, degrees")
        plt.xlabel("RA, degrees")

        plt.axis('off')

    if outfile == '':
        plt.show()
    else:
        plt.savefig(outfile)
    return


def plot_phase_core (ax, t, x, xerr, period, offset=0, color='k'):
    ''' Plots a pretty period-folded lightcurve on a given axes object.

    Doesn't assume anything about your data (e.g. that it's in magnitudes)
    '''
    # Untested.
    
    phase = ((t % period) / period + offset) % 1.

    ax.errorbar(phase, x, yerr=xerr, fmt= color+'o')
    ax.errorbar(phase-1,x,yerr=xerr,fmt='o', mfc='0.7',mec='0.7', 
                 ecolor='0.7')
    ax.errorbar(phase+1,x,yerr=xerr,fmt='o', mfc='0.7',mec='0.7', 
                 ecolor='0.7')
    
    ax.set_xticks( [0, 0.5, 1] )
    ax.set_xticks( np.arange(-.5,1.5,.1), minor=True)

    ax.set_xlim(-0.25, 1.25)

    return period


def plot_phase (table, sid, period, band='j', outfile='', season=123,
                offset=0,clear=True):
    ''' Plots magnitude as a function of phase for one source in one band. '''

    # Yo! This needs one more argument! A constant term to shift over by...
    # I'll make it scaled to phase (i.e. between 0 and 1)

    if band.lower() not in ['j','h','k','jmh','hmk']:
        print "Error: keyword 'band' must be 'j','h', 'k', 'jmh', or 'hmk'."
        print "Keyword 'band' defaulting to 'j'."
        band = 'j'
    
    if 'm' in band.lower():
        bandname = band.upper() + "PNT"
        thing = 'color'
    else:
        bandname = band.upper() + "APERMAG3"
        thing = 'magnitude'

#    w = numpy.where( table.SOURCEID == sid )

    s_table = season_cut(table, sid, season)

    # ra1, ra2 = 314.36, 315.77
    # dec1,dec2= 52.02, 52.92

    # rabox = [ra1, ra1, ra2, ra2, ra1]
    # decbox= [dec1,dec2,dec2,dec1,dec1]

    # sra, sdec = table.RA[w][0], table.DEC[w][0]

    date = s_table.MEANMJDOBS #- 54579
    phase = ((date % period) / period + offset) % 1.

    mag = s_table.data[bandname]
    err = s_table.data[bandname+"ERR"]

    if clear:
        plt.clf()

    ax = plt.gca()

    plt.errorbar(phase,mag,yerr=err,fmt='ko')#,ecolor='k')
#    plt.errorbar(phase-1,mag,yerr=err,fmt='ko',ecolor='0.7', alpha=0.3) 
    plt.errorbar(phase-1,mag,yerr=err,fmt='o', mfc='0.7',mec='0.7', 
                 ecolor='0.7')
    plt.errorbar(phase+1,mag,yerr=err,fmt='o', mfc='0.7',mec='0.7', 
                 ecolor='0.7')
#    plt.errorbar(phase+1,mag,yerr=err,fmt='ko',ecolor='0.7', alpha=0.3) 


    plt.xticks( [0, 0.5, 1] )
    ax.set_xticks( np.arange(-.5,1.5,.1), minor=True)

    plt.xlim(-0.25, 1.25)

    if len(band) == 1:
        plt.gca().invert_yaxis()

    plt.ylabel("WFCAM %s %s" % (band.upper(), thing))
    plt.xlabel("Phase")

    if period < 1:
        period_string = "%f hours" % (period*24)
        print period_string
    else:
        period_string = "%f days" % period

    plt.title ("Phase-folded lightcurve. Source ID %d. Period: %s." %
               (sid, period_string))

    if outfile == '':
        plt.show()
    else:
        plt.savefig(outfile)
    return

    
