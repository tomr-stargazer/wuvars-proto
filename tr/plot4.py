"""
plot4.py

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


class StarData(object):

    def __init__(self, table, sid):

        self.table = table
        self.sid = sid

        # Loading data
        self.s_table = data_cut (table, sid)

        if len(self.s_table) == 0:
            raise ValueError("no data here")

        # Once with no errors (normal),
        # once with small errors (info).

        ## First: normal
        # Use band_cut to get relevant data chunks.
        j_table = band_cut(self.s_table, 'j', max_flag=0)
        h_table = band_cut(self.s_table, 'h', max_flag=0)
        k_table = band_cut(self.s_table, 'k', max_flag=0)

        # get a date (x-axis) for each plot
        self.jdate = j_table.MEANMJDOBS
        self.hdate = h_table.MEANMJDOBS
        self.kdate = k_table.MEANMJDOBS
        
        # get a magnitude (y-axis) for each plot
        self.jcol = j_table.JAPERMAG3
        self.hcol = h_table.HAPERMAG3
        self.kcol = k_table.KAPERMAG3

        # get a magnitude error (y-error) for each plot
        self.jerr = j_table.JAPERMAG3ERR
        self.herr = h_table.HAPERMAG3ERR
        self.kerr = k_table.KAPERMAG3ERR

        # get a quality flag for each plot
        # (aftercomment: not sure we're going to ever use these)
        self.jflag = j_table.JPPERRBITS
        self.hflag = h_table.HPPERRBITS
        self.kflag = k_table.KPPERRBITS

        ## Second: info
        j_table_info = band_cut(self.s_table, 'j', min_flag=1, max_flag=256)
        h_table_info = band_cut(self.s_table, 'h', min_flag=1, max_flag=256)
        k_table_info = band_cut(self.s_table, 'k', min_flag=1, max_flag=256)

        # get a date (x-axis) for each plot
        self.jdate_info = j_table_info.MEANMJDOBS
        self.hdate_info = h_table_info.MEANMJDOBS
        self.kdate_info = k_table_info.MEANMJDOBS

        # get a magnitude (y-axis) for each plot
        self.jcol_info = j_table_info.JAPERMAG3
        self.hcol_info = h_table_info.HAPERMAG3
        self.kcol_info = k_table_info.KAPERMAG3

        # get a magnitude error (y-error) for each plot
        self.jerr_info = j_table_info.JAPERMAG3ERR
        self.herr_info = h_table_info.HAPERMAG3ERR
        self.kerr_info = k_table_info.KAPERMAG3ERR

        # get a quality flag for each plot
        self.jflag = j_table_info.JPPERRBITS
        self.hflag = h_table_info.HPPERRBITS
        self.kflag = k_table_info.KPPERRBITS

        # We'll use different data-cuts for the two different plots.
        # Relevant comment: I made an executive call to include only
        # 'normal' and 'info'-flagged data in the C-C and C-M plots
        # (i.e. max_flag=256 in all relevant bands).
        
        # In the color-mag plot, we need data where H and K are defined 
        # everywhere. That's two cuts.
        khk_table = band_cut( band_cut(self.s_table, 'k', max_flag=256),
                              'h', max_flag=256)

        self.khkdate = khk_table.MEANMJDOBS #- date_offset
        self.k_khk = khk_table.KAPERMAG3
        self.hmk_khk = khk_table.HMKPNT
        self.k_khk_err = khk_table.KAPERMAG3ERR
        self.hmk_khk_err = khk_table.HMKPNTERR

        # In the color-color plot, we need data where J, H, and K are
        # defined everywhere. That's one more cut.
        jhk_table = band_cut(khk_table, 'j', max_flag=256)

        self.jhkdate = jhk_table.MEANMJDOBS #- date_offset
        self.jmh_jhk = jhk_table.JMHPNT
        self.hmk_jhk = jhk_table.HMKPNT
        self.jmh_jhk_err = jhk_table.JMHPNTERR
        self.hmk_jhk_err = jhk_table.HMKPNTERR

    def get_columns(self, band, flags=0):

        b_table = band_cut(self.s_table, band, max_flag=flags)

        columns = {}

        columns['date'] = b_table['MEANMJDOBS']
        columns['mag'] = b_table['{0}APERMAG3'.format(band.upper())]
        columns['err'] = b_table['{0}APERMAG3ERR'.format(band.upper())]
        columns['flag'] = b_table['{0}PPERRBITS'.format(band.upper())]
        try:
            columns['grade'] = b_table['{0}GRADE'.format(band.upper())]
        except:
            columns['grade'] = np.zeros_like(columns['mag'])

        return columns


def lightcurve_axes_with_info(stardata, band, axes, colorscale, 
                              cmap, vmin, vmax):


        columns = stardata.get_columns(band, flags=0)
        columns_info = stardata.get_columns(band, flags=256)

        if len(columns['date']) > 0:
            # First, plot the errorbars, with no markers, in the background:
            axes.errorbar( columns['date'], columns['mag'], marker=None,
                                 yerr=columns['err'], fmt=None, ecolor='k',
                                 zorder=0)
            
            # Next, scatter the points themselves, colored re:grade :
            axes.scatter( columns['date'], columns['mag'], cmap=cmap,
                                c=columns[colorscale], vmin=vmin, vmax=vmax, zorder=100)

        if len(columns_info['date']) > 0:
            # First, plot the errorbars, with no markers, in the background:
            axes.errorbar( columns_info['date'], columns_info['mag'], 
                                 yerr=columns_info['err'], marker=None,
                                 fmt=None, ecolor='k', zorder=0)

            # Next, scatter the points themselves, colored re:grade :
            axes.scatter( columns_info['date'], columns_info['mag'], 
                                marker='d', 
                                c=columns_info[colorscale], cmap=cmap, 
                                vmin=vmin, vmax=vmax, zorder=100)

        # Finally, flip it (magnitudes are backwards).
        axes.invert_yaxis()

        # And plot the dotted lines, if relevant.
        # if abridged:
        #     d_ax[band].plot([ab_s1s2line, ab_s1s2line], [0,30], "k--",
        #                     scaley=False, scalex=False)

        #     d_ax[band].plot([ab_s2s3line, ab_s2s3line], [0,30], "k--",
        #                     scaley=False, scalex=False)

def basic_lc(stardata):

    # kwargs defaulting over
    timecolor = True
    time_cmap = 'jet'
    color_slope = False

    if timecolor is True:
        colorscale='date'

    fig = plt.figure(figsize = (10, 6), dpi=80, facecolor='w', edgecolor='k')

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

    d_date = {'j': stardata.jdate, 'h': stardata.hdate, 'k': stardata.kdate}
    d_date_info = {'j': stardata.jdate_info, 'h': stardata.hdate_info, 'k': stardata.kdate_info}

    # d_rawdate = {'j': stardata.raw_jdate, 'h': stardata.raw_hdate, 'k': stardata.raw_kdate}
    # d_rawdate_info = {
    #     'j': stardata.raw_jdate_info, 'h': stardata.raw_hdate_info, 'k': stardata.raw_kdate_info}

    color_vmin = stardata.s_table.MEANMJDOBS.min() 
    color_vmax = stardata.s_table.MEANMJDOBS.max() 
    # min(
    #     (d_rawdate['j'].min(), d_rawdate['h'].min(), d_rawdate['k'].min(),
    #      d_rawdate_info['j'].min(), d_rawdate_info['h'].min(), d_rawdate_info['k'].min()))
    # color_vmax = max(
    #     (d_rawdate['j'].max(), d_rawdate['h'].max(), d_rawdate['k'].max(),
    #      d_rawdate_info['j'].max(), d_rawdate_info['h'].max(), d_rawdate_info['k'].max()))

    if timecolor:
        d_c = d_date
        d_c_info = d_date_info
        vmin = color_vmin
        vmax = color_vmax
    else:
        d_c = d_grade
        d_c_info = d_grade_info
        vmin = 0.8
        vmax = 1

    for band in ['j', 'h', 'k']:

        lightcurve_axes_with_info(stardata, band, d_ax[band], colorscale, 
                                  cmap=d_cmap[band], vmin=vmin, vmax=vmax)

            


    ## Now let's do the 2 color-mag/color-color plots.

    # Plot J-H vs H-K using the "jhk_" variables.
    try:
        plot_trajectory_core( ax_jhk, stardata.hmk_jhk, stardata.jmh_jhk, stardata.jhkdate,
                              vmin=color_vmin, vmax=color_vmax) 

        if color_slope:
            jhk_slope, jhk_intercept, slope_err = (
                slope(stardata.hmk_jhk, stardata.jmh_jhk, stardata.hmk_jhk_err, stardata.jmh_jhk_err,
                      verbose=False) )
            
            ax_jhk.plot([0, 6], [jhk_intercept, jhk_intercept + 6*jhk_slope], 
                        ':', scalex=False, scaley=False)
            
    except Exception as e:
        print "JHK plot broke: {0}".format(e)
        pass
        

    # Plot K vs H-K using the "khk_" variables.
    try:
        plot_trajectory_core( ax_khk, stardata.hmk_khk, stardata.k_khk, stardata.khkdate,
                              ms=False, ctts=False, 
                              vmin=color_vmin, vmax=color_vmax) 

        # plot boundaries are manually set for readability, if necessary
        if len(ax_khk.get_xticks()) > 7:
            khk_xmin = np.floor(stardata.hmk_khk.min() * 0.95 * 20)/20.
            khk_xmax = np.ceil( stardata.hmk_khk.max() * 1.05 * 20)/20.

            khk_xticks = np.linspace(khk_xmin, khk_xmax, 6)
            ax_khk.set_xticks(khk_xticks)

        if color_slope:
            khk_slope, khk_intercept, slope_err = (
                slope(hmk_khk, k_khk, hmk_khk_err, k_khk_err,
                      verbose=False) )
            
            ax_khk.plot([0, 6], [khk_intercept, khk_intercept + 6*khk_slope],
                        '--', scalex=False, scaley=False)
    
    except Exception as e:
        print "KHK plot broke: {0}".format(e)
        pass
    ax_khk.invert_yaxis()

    # Hide the bad labels...
    plt.setp(ax_j.get_xticklabels(), visible=False)
    plt.setp(ax_h.get_xticklabels(), visible=False)

    # Label stuff
#    ax_k.set_xlabel( "Time (JD since 01/01/2000)" )

    ax_j.set_ylabel( "J",{'rotation':'horizontal', 'fontsize':'large'} )
    ax_h.set_ylabel( "H",{'rotation':'horizontal', 'fontsize':'large'} )
    ax_k.set_ylabel( "K",{'rotation':'horizontal', 'fontsize':'large'} )

    ax_jhk.set_xlabel( "H-K" )
    ax_jhk.set_ylabel( "J-H")#, {'rotation':'horizontal'})
    ax_khk.set_xlabel( "H-K" )
    ax_khk.set_ylabel( "K")#, {'rotation':'horizontal'})

    fig.ax_k = ax_k
    fig.ax_h = ax_h
    fig.ax_j = ax_j
    fig.ax_jhk = ax_jhk
    fig.ax_khk = ax_khk

    return fig
