"""
plot4.py

"""

from __future__ import division

import numpy as np
import matplotlib.pyplot as plt

from helpers3 import data_cut, band_cut
from plot2 import plot_trajectory_core
from chi2 import test_analyze, diagnostic_analyze
from scargle import fasper as lsp
from timing import lsp_mask, lsp_tuning
from spread3 import Stetson_machine
from abridger import abridger
from color_slope import slope


class StarData(object):

    """
    StarData contains the data for a single star,
    plus methods to access its relevant columns.
    Used primarily for making light curves -- probably
    too slow for running computations.

    Example use:

        >>> from plot4 import StarData. basic_lc
        >>> sd = StarData(variables_photometry, 44199508443333, date_offset=54034)
        >>> light_curve = basic_lc(sd)

    """

    def __init__(self, table, sid, date_offset=0, name=None, abridger=None):

        self.sid = sid
        self.date_offset = date_offset

        # Loading data
        self.s_table = data_cut (table, sid)

        self.min_date = self.s_table.MEANMJDOBS.min() - date_offset
        self.max_date = self.s_table.MEANMJDOBS.max() - date_offset

        self.abridger = abridger

        if len(self.s_table) == 0:
            raise ValueError("no data here")

        if name is not None:
            self.name = name
        else:
            self.name = ''

    def get_columns(self, band, max_flag=0, min_flag=0):
        """
        Returns relevant columns for a given photometry band.

        `band` must be 'j', 'h', or 'k'.

        """

        if band.lower() not in ('j', 'h', 'k'):
            raise ValueError("Invalid band: {0} not in ('j', 'h', 'k')".format(band.lower()))

        b_table = band_cut(self.s_table, band, max_flag=max_flag, min_flag=min_flag)

        columns = {}

        columns['date'] = b_table['MEANMJDOBS'] - self.date_offset
        columns['mag'] = b_table['{0}APERMAG3'.format(band.upper())]
        columns['err'] = b_table['{0}APERMAG3ERR'.format(band.upper())]
        columns['flag'] = b_table['{0}PPERRBITS'.format(band.upper())]
        try:
            columns['grade'] = b_table['{0}GRADE'.format(band.upper())]
        except:
            columns['grade'] = np.zeros_like(columns['mag'])

        return columns

    def get_colormag_columns(self, band, max_flag=256, min_flag=0):
        """
        Returns relevant columns for color+magnitude pair.

        """

        if band.lower() not in ('jjh', 'khk'):
            raise ValueError("Invalid color-mag combination: {0} not in ('jjh', 'khk')".format(band.lower()))

        mag, blue, red = band.lower()

        colormag_table = band_cut(band_cut(self.s_table, red, max_flag=max_flag),
                                  blue, max_flag=max_flag)

        columns = {}

        columns['date'] = colormag_table['MEANMJDOBS'] - self.date_offset
        columns['mag'] = colormag_table['{0}APERMAG3'.format(mag.upper())]
        columns['color'] = colormag_table['{0}M{1}PNT'.format(blue.upper(), red.upper())]
        columns['mag_err'] = colormag_table['{0}APERMAG3ERR'.format(mag.upper())]
        columns['color_err'] = colormag_table['{0}M{1}PNTERR'.format(blue.upper(), red.upper())]

        return columns


    def get_colorcolor_columns(self, max_flag=256, min_flag=0):
        """
        Returns relevant columns for the J-H, H-K color+color pair.

        """

        colorcolor_table = band_cut(band_cut(band_cut(self.s_table, 'k', max_flag=max_flag),
                                    'h', max_flag=max_flag), 'j', max_flag=max_flag)

        columns = {}

        columns['date'] = colorcolor_table['MEANMJDOBS'] - self.date_offset
        columns['jmh'] = colorcolor_table['JMHPNT']
        columns['hmk'] = colorcolor_table['HMKPNT']
        columns['jmh_err'] = colorcolor_table['JMHPNTERR']
        columns['hmk_err'] = colorcolor_table['HMKPNTERR']

        return columns


def lightcurve_axes_with_info(stardata, band, axes, colorscale, cmap, vmin, vmax):

    columns = stardata.get_columns(band, max_flag=0)
    columns_info = stardata.get_columns(band, min_flag=1, max_flag=256)

    date = np.copy(columns['date'])
    date_info = np.copy(columns_info['date'])

    if stardata.abridger:
        bridge = stardata.abridger(stardata, flags=256)
        # this logic should get moved into StarData...
        # The following uses the signature of wuvars-proto/tr/abridger.py
        date[date > bridge['s1_s2_bound']] -= bridge['s2_subtraction_factor']
        date[date > bridge['s2_s3_bound'] - bridge['s2_subtraction_factor']] -= bridge['s3_subtraction_factor']
        date_info[date_info > bridge['s1_s2_bound']] -= bridge['s2_subtraction_factor']
        date_info[date_info > bridge['s2_s3_bound'] - bridge['s2_subtraction_factor']] -= bridge['s3_subtraction_factor']

    if len(columns['date']) > 0:
        # First, plot the errorbars, with no markers, in the background:
        axes.errorbar( date, columns['mag'], marker=None,
                             yerr=columns['err'], fmt=None, ecolor='k',
                             zorder=0)
        
        # Next, scatter the points themselves, colored re:colorscale :
        axes.scatter( date, columns['mag'], cmap=cmap,
                            c=columns[colorscale], vmin=vmin, vmax=vmax, zorder=100)

    if len(columns_info['date']) > 0:
        # First, plot the errorbars, with no markers, in the background:
        axes.errorbar( date_info, columns_info['mag'], 
                             yerr=columns_info['err'], marker=None,
                             fmt=None, ecolor='k', zorder=0)

        # Next, scatter the points themselves, colored re:colorscale :
        axes.scatter( date_info, columns_info['mag'], 
                            marker='d', 
                            c=columns_info[colorscale], cmap=cmap, 
                            vmin=vmin, vmax=vmax, zorder=100)

    # Finally, flip it (magnitudes are backwards).
    axes.invert_yaxis()

    # don't go negative on the X axis ever
    if stardata.min_date >= 0:
        xlims = axes.get_xlim()
        axes.set_xlim( max(xlims[0], 0), xlims[1] )

    if stardata.abridger:
        axes.plot([bridge['s1_s2_line'], bridge['s1_s2_line']], [0,30], "k--",
                        scaley=False, scalex=False)

        axes.plot([bridge['s2_s3_line'], bridge['s2_s3_line']], [0,30], "k--",
                        scaley=False, scalex=False)

        axes.set_xticks(bridge['xticks'])
        axes.set_xticklabels(bridge['xticklabels'])
        axes.set_xlim(bridge['xlim_bounds'])

    axes.get_figure().canvas.draw()


def phase_axes_with_info(stardata, band, period, axes, colorscale, cmap, vmin, vmax, offset=0):

    columns = stardata.get_columns(band, max_flag=0)
    columns_info = stardata.get_columns(band, min_flag=1, max_flag=256)

    date = np.copy(columns['date'])
    date_info = np.copy(columns_info['date'])

    phase = ((date % period) / period + offset) % 1
    phase_info = ((date_info % period) / period + offset) % 1

    if len(columns['date']) > 0:
        # plot the greyed-out versions on left and right
        axes.errorbar(phase-1,columns['mag'],yerr=columns['err'],fmt='o', mfc='0.7',mec='0.7', 
                     ecolor='0.7', ms=6, zorder=-5)
        axes.errorbar(phase+1,columns['mag'],yerr=columns['err'],fmt='o', mfc='0.7',mec='0.7', 
                     ecolor='0.7', ms=6, zorder=-5)

        # First, plot the errorbars, with no markers, in the background:
        axes.errorbar(phase, columns['mag'], marker=None,
                      yerr=columns['err'], fmt=None, ecolor='k',
                      zorder=0)
        # Next, scatter the points themselves, colored re:colorscale :
        axes.scatter(phase, columns['mag'], cmap=cmap,
                     c=columns[colorscale], vmin=vmin, vmax=vmax, zorder=100)

    if len(columns_info['date']) > 0:
        # plot the greyed-out versions on left and right
        axes.errorbar(phase_info-1,columns_info['mag'],yerr=columns_info['err'],fmt='d', mfc='0.7',mec='0.7', 
                     ecolor='0.7', ms=6, zorder=-5)
        axes.errorbar(phase_info+1,columns_info['mag'],yerr=columns_info['err'],fmt='d', mfc='0.7',mec='0.7', 
                     ecolor='0.7', ms=6, zorder=-5)

        # First, plot the errorbars, with no markers, in the background:
        axes.errorbar(phase_info, columns_info['mag'], 
                      yerr=columns_info['err'], marker=None,
                      fmt=None, ecolor='k', zorder=0)
        # Next, scatter the points themselves, colored re:colorscale :
        axes.scatter(phase_info, columns_info['mag'], 
                     marker='d', 
                     c=columns_info[colorscale], cmap=cmap, 
                     vmin=vmin, vmax=vmax, zorder=100)

    # Finally, flip it (magnitudes are backwards).
    axes.invert_yaxis()

    axes.set_xticks( [0, 0.5, 1] )
    axes.set_xticks( np.arange(-.5,1.5,.1), minor=True)

    axes.set_xlim(-0.25, 1.25)

    axes.get_figure().canvas.draw()


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
    axes.get_figure().canvas.draw()    


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

        axes.get_figure().canvas.draw()            
    
    except Exception as e:
        print "Color-color plot broke: {0}".format(e)
        pass

def basic_lc(stardata, timecolor=True, custom_xlabel=False):
    """
    Proof-of-concept reimplementation of plot3.graded_lc.

    Fewer bells and whistles, but looks perfect, and is much cleaner.

    """

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

    d_ax = {'j': ax_j, 'h': ax_h, 'k': ax_k}
    
    if timecolor:
        d_cmap = {'j': time_cmap, 'h': time_cmap, 'k': time_cmap}
    elif type(d_cmap) is str:
        d_cmap = {'j': d_cmap, 'h': d_cmap, 'k': d_cmap}
    elif type(d_cmap) is not dict:
        d_cmap = {'j': d_cmap[0], 'h': d_cmap[1], 'k': d_cmap[2]}

    color_vmin = stardata.min_date
    color_vmax = stardata.max_date 

    if timecolor:
        vmin = color_vmin
        vmax = color_vmax
    else:
        vmin = 0.8
        vmax = 1

    for band in ['j', 'h', 'k']:
        lightcurve_axes_with_info(stardata, band, d_ax[band], colorscale, 
                                  cmap=d_cmap[band], vmin=vmin, vmax=vmax)

    colorcolor_axes(stardata, ax_jhk, colorscale, cmap='jet', vmin=vmin, vmax=vmax,
                    color_slope=color_slope)
    colormag_axes(stardata, 'khk', ax_khk, colorscale, cmap='jet', vmin=vmin, vmax=vmax,
                  color_slope=color_slope)

    # Hide the bad labels...
    plt.setp(ax_j.get_xticklabels(), visible=False)
    plt.setp(ax_h.get_xticklabels(), visible=False)

    # Label stuff
#    ax_k.set_xlabel( "Time (JD since 01/01/2000)" )
    if custom_xlabel:
        ax_k.set_xlabel( custom_xlabel )
    else:
        ax_k.set_xlabel( "Time (MJD - %.1f)" % stardata.date_offset )

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

def basic_phase(stardata, period, timecolor=True, offset=0):
    """
    Proof-of-concept reimplementation of plot3.graded_phase.

    This is begging to be refactored, with basic_lc, into a single function.
    Proposed idea: Make the single function pretty general and take a couple extra "messy"
    kwargs, and then define basic_phase and basic_lc as partials of the generic version.

    """

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

    d_ax = {'j': ax_j, 'h': ax_h, 'k': ax_k}
    
    if timecolor:
        d_cmap = {'j': time_cmap, 'h': time_cmap, 'k': time_cmap}
    elif type(d_cmap) is str:
        d_cmap = {'j': d_cmap, 'h': d_cmap, 'k': d_cmap}
    elif type(d_cmap) is not dict:
        d_cmap = {'j': d_cmap[0], 'h': d_cmap[1], 'k': d_cmap[2]}

    color_vmin = stardata.min_date
    color_vmax = stardata.max_date

    if timecolor:
        vmin = color_vmin
        vmax = color_vmax
    else:
        vmin = 0.8
        vmax = 1

    for band in ['j', 'h', 'k']:
        phase_axes_with_info(stardata, band, period, d_ax[band], colorscale, 
                             cmap=d_cmap[band], vmin=vmin, vmax=vmax, offset=offset)

    colorcolor_axes(stardata, ax_jhk, colorscale, cmap='jet', vmin=vmin, vmax=vmax,
                    color_slope=color_slope)
    colormag_axes(stardata, 'khk', ax_khk, colorscale, cmap='jet', vmin=vmin, vmax=vmax,
                  color_slope=color_slope)

    # Hide the bad labels...
    plt.setp(ax_j.get_xticklabels(), visible=False)
    plt.setp(ax_h.get_xticklabels(), visible=False)

    ax_j.set_ylabel( "J",{'rotation':'horizontal', 'fontsize':'large'} )
    ax_h.set_ylabel( "H",{'rotation':'horizontal', 'fontsize':'large'} )
    ax_k.set_ylabel( "K",{'rotation':'horizontal', 'fontsize':'large'} )

    ax_k.set_xlabel("Phase (Period = {0:.4} days)".format(float(period)))

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


def multi_lightcurve(stardatas, dimensions, bands, cmap='jet', colorscale='date'):

    xdim, ydim = dimensions

    if len(stardatas) > (xdim * ydim):
        raise ValueError("Number of input stars should be less than or equal to product of dimensions")
    elif len(stardatas) != len(bands):
        raise ValueError("List of bands should be same length as list of input stars")

    fig = plt.figure(figsize = (1.5+xdim*5, 0.6+ydim*1.8), 
                     dpi=80, facecolor='w', edgecolor='k')

    # single colorscale across all light curves
    if colorscale == 'date':
        vmin = min([stardata.min_date for stardata in stardatas])
        vmax = max([stardata.max_date for stardata in stardatas])
    elif colorscale == 'grade':
        vmin = 0.8
        vmax = 1.0

    fig.xlim = (0,0)
    fig.xticks = []
    fig.xticklabels = []

    for stardata, band, i in zip(stardatas, bands, range(1, 1+len(stardatas))):

        if i == 1: sharex = None
        else: sharex = fig.ax1

        ax = fig.add_subplot(ydim, xdim, i, sharex=sharex)

        lightcurve_axes_with_info(stardata, band, ax, 'date', 
                                  cmap=cmap, vmin=vmin, vmax=vmax)

        ax.set_ylabel( band.upper(),{'rotation':'horizontal', 'fontsize':'large'} )

        if i <= len(bands) - xdim:
            plt.setp(ax.get_xticklabels(), visible=False)

        fig.__setattr__('ax{0}'.format(i), ax)

        fig.xlim = (min(fig.xlim[0], ax.get_xlim()[0]), max(fig.xlim[1], ax.get_xlim()[1]))
        if len(ax.get_xticks()) > len(fig.xticks):
            fig.xticks = ax.get_xticks()
            fig.xticklabels = [x.get_text() for x in ax.get_xticklabels()]

        ax.set_xlim(fig.xlim)
        ax.set_xticks(fig.xticks)
        ax.set_xticklabels(fig.xticklabels)

        ax.text(0.7, 0.1, stardata.name, transform=ax.transAxes, fontsize='small')

    fig.canvas.draw()
    return fig


