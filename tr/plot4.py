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

    def get_columns(self, band, max_flag=0, min_flag=0):

        b_table = band_cut(self.s_table, band, max_flag=max_flag, min_flag=min_flag)

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

    def get_colormag_columns(self, band, max_flag=256, min_flag=0):

        # So we have some options: we might want 
        # J, J-H,
        # K, H-K 
        if band.lower() not in ('jjh', 'khk'):
            raise ValueError("Invalid color-mag combination: {0} not in ('jjh', 'khk')".format(band.lower()))

        mag = band.lower()[0]
        red = band.lower()[-1]
        blue = band.lower()[-2]

        colormag_table = band_cut(band_cut(self.s_table, red, max_flag=max_flag),
                                  blue, max_flag=max_flag)

        columns = {}

        columns['date'] = colormag_table['MEANMJDOBS'] #- date_offset
        columns['mag'] = colormag_table['{0}APERMAG3'.format(mag.upper())]
        columns['color'] = colormag_table['{0}M{1}PNT'.format(blue.upper(), red.upper())]
        columns['mag_err'] = colormag_table['{0}APERMAG3ERR'.format(mag.upper())]
        columns['color_err'] = colormag_table['{0}M{1}PNTERR'.format(blue.upper(), red.upper())]

        return columns


    def get_colorcolor_columns(self, max_flag=256, min_flag=0):

        colorcolor_table = band_cut(band_cut(band_cut(self.s_table, 'k', max_flag=max_flag),
                                  'h', max_flag=max_flag), 'j', max_flag=max_flag)

        columns = {}

        columns['date'] = colorcolor_table['MEANMJDOBS'] #- date_offset
        columns['jmh'] = colorcolor_table['JMHPNT']
        columns['hmk'] = colorcolor_table['HMKPNT']
        columns['jmh_err'] = colorcolor_table['JMHPNTERR']
        columns['hmk_err'] = colorcolor_table['HMKPNTERR']

        return columns


def lightcurve_axes_with_info(stardata, band, axes, colorscale, 
                              cmap, vmin, vmax):


        columns = stardata.get_columns(band, max_flag=0)
        columns_info = stardata.get_columns(band, min_flag=1, max_flag=256)

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


def colormag_axes(stardata, band, axes, colorscale, cmap, vmin, vmax, color_slope=False):

    colormag_columns = stardata.get_colormag_columns(band, max_flag=256)

    try:
        plot_trajectory_core(axes, colormag_columns['color'], colormag_columns['mag'], colormag_columns['date'],
                             ms=False, ctts=False, 
                             vmin=vmin, vmax=vmax)

        # plot boundaries are manually set for readability, if necessary
        if len(axes.get_xticks()) > 7:
            xmin = np.floor(colormag_columns['color'].min() * 0.95 * 20)/20.
            xmax = np.ceil( colormag_columns['color'].max() * 1.05 * 20)/20.

            xticks = np.linspace(xmin, xmax, 6)
            axes.set_xticks(xticks)

        if color_slope:
            slope_color, color_intercept, slope_err = (
                slope(colormag_columns['color'], colormag_columns['mag'], 
                      colormag_columns['color_err'], colormag_columns['mag_err'],
                      verbose=False) )
            
            axes.plot([0, 6], [color_intercept, color_intercept + 6*slope_color],
                        '--', scalex=False, scaley=False)
    
    except Exception as e:
        print "Color-mag plot broke: {0}".format(e)
        pass
    axes.invert_yaxis()


def colorcolor_axes(stardata, axes, colorscale, cmap, vmin, vmax, color_slope=False):

    colorcolor_columns = stardata.get_colorcolor_columns(max_flag=256)

    try:
        plot_trajectory_core(axes, colorcolor_columns['hmk'], colorcolor_columns['jmh'], colorcolor_columns['date'],
                             vmin=vmin, vmax=vmax)

        # plot boundaries are manually set for readability, if necessary
        if len(axes.get_xticks()) > 7:
            xmin = np.floor(colorcolor_columns['hmk'].min() * 0.95 * 20)/20.
            xmax = np.ceil( colorcolor_columns['hmk'].max() * 1.05 * 20)/20.

            xticks = np.linspace(xmin, xmax, 6)
            axes.set_xticks(xticks)

        if color_slope:
            colorcolor_slope, jmh_intercept, slope_err = (
                slope(colorcolor_columns['hmk'], colorcolor_columns['jmh'], 
                      colorcolor_columns['hmk_err'], colorcolor_columns['jmh_err'],
                      verbose=False) )
            
            axes.plot([0, 6], [jmh_intercept, jmh_intercept + 6*colorcolor_slope],
                        '--', scalex=False, scaley=False)
    
    except Exception as e:
        print "Color-color plot broke: {0}".format(e)
        pass

def basic_lc(stardata, timecolor=True):

    # kwargs defaulting over
    time_cmap = 'jet'
    color_slope = False
    d_cmap={'j':'Blues', 'h': 'Greens', 'k': 'Reds'}

    if timecolor is True:
        colorscale='date'
    else:
        colorscale='grade'

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
        vmin = color_vmin
        vmax = color_vmax
    else:
        vmin = 0.8
        vmax = 1

    for band in ['j', 'h', 'k']:

        lightcurve_axes_with_info(stardata, band, d_ax[band], colorscale, 
                                  cmap=d_cmap[band], vmin=vmin, vmax=vmax)

    ## Now let's do the 2 color-mag/color-color plots.

    colorcolor_axes(stardata, ax_jhk, colorscale, cmap='jet', vmin=vmin, vmax=vmax)
    colormag_axes(stardata, 'khk', ax_khk, colorscale, cmap='jet', vmin=vmin, vmax=vmax)


    # # Plot J-H vs H-K using the "jhk_" variables.
    # try:
    #     plot_trajectory_core( ax_jhk, stardata.hmk_jhk, stardata.jmh_jhk, stardata.jhkdate,
    #                           vmin=color_vmin, vmax=color_vmax) 

    #     if color_slope:
    #         jhk_slope, jhk_intercept, slope_err = (
    #             slope(stardata.hmk_jhk, stardata.jmh_jhk, stardata.hmk_jhk_err, stardata.jmh_jhk_err,
    #                   verbose=False) )
            
    #         ax_jhk.plot([0, 6], [jhk_intercept, jhk_intercept + 6*jhk_slope], 
    #                     ':', scalex=False, scaley=False)
            
    # except Exception as e:
    #     print "JHK plot broke: {0}".format(e)
    #     pass
        

    # # Plot K vs H-K using the "khk_" variables.
    # try:
    #     plot_trajectory_core( ax_khk, stardata.hmk_khk, stardata.k_khk, stardata.khkdate,
    #                           ms=False, ctts=False, 
    #                           vmin=color_vmin, vmax=color_vmax) 

    #     # plot boundaries are manually set for readability, if necessary
    #     if len(ax_khk.get_xticks()) > 7:
    #         khk_xmin = np.floor(stardata.hmk_khk.min() * 0.95 * 20)/20.
    #         khk_xmax = np.ceil( stardata.hmk_khk.max() * 1.05 * 20)/20.

    #         khk_xticks = np.linspace(khk_xmin, khk_xmax, 6)
    #         ax_khk.set_xticks(khk_xticks)

    #     if color_slope:
    #         khk_slope, khk_intercept, slope_err = (
    #             slope(hmk_khk, k_khk, hmk_khk_err, k_khk_err,
    #                   verbose=False) )
            
    #         ax_khk.plot([0, 6], [khk_intercept, khk_intercept + 6*khk_slope],
    #                     '--', scalex=False, scaley=False)
    
    # except Exception as e:
    #     print "KHK plot broke: {0}".format(e)
    #     pass
    # ax_khk.invert_yaxis()

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
