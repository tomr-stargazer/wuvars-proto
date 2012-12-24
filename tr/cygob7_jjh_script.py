"""
This is a script for generating J, J-H plots in Cyg OB7 in 
a repeatable way.

4 December 2012 - Tom Rice t.rice90@gmail.com - created

"""

import numpy as np
import matplotlib.pyplot as plt

import atpy

import plot3, spread3

from cygob7_jjh_script_helper import *

# where the data is stored
path = '/home/tom/Dropbox/Cyg_OB7/paper2/'
path2 ='/home/tom/Dropbox/Cyg_OB7/paper2/book/December2012/' 
path3 ='/home/tom/reu/DATA/Merged_Catalogs/September_2011/'

# The data file that contains data for all of the wise disks, aspin sources,
# and wise transition disks
watso = atpy.Table(path+"Watso_ALLDATA_cleaned_errorcorrected_ce.fits")

# removing a ~single date datapoint
watso = watso.where(watso.MEANMJDOBS > 54570)
 
# general data file that has most of the wise_extras
data = atpy.Table(path+"cleaned_errorcorrected_cecorrected_data.fits")

# special data file for 2 super faint wise_extras that were missed in the
# previous passes
welo = atpy.Table(
    path+"WISE_extras_lonelytwo_only_ALLDATA_errorcorrected_ce.fits")

# RWA data
rwa_stats = atpy.Table(path3+"disk_stats_radec_withRWA_extended17.fits")

rwa_sources = rwa_stats.SOURCEID
rwa_names = ["RWA "+str(x) for x in rwa_stats.Designation]

rwa16 = 44027710020132
rwa_data = data.where( 
    np.in1d(data.SOURCEID,rwa_sources) & 
    (((data.JPPERRBITS + data.HPPERRBITS + data.KPPERRBITS) == 0) |
     (data.SOURCEID == rwa16)) )

# This is no good because of errorbars not being calibrated:
#rwa_data = atpy.Table(path3+"disk_data_extended17.fits")


# We have created lists for Wise disks, Wise extras, Aspin sources, and 
# Transition disks. variables:
# `wise_disks`, `wise_extras`, `aspin_sources`, `wise_trans`
# with identifiers as the variable name + "_names"

# "real" wise disks - to filter out non-matches
rwd = wise_disks > 1
ras = aspin_sources > 1
rwt = wise_trans > 1
rwe = wise_extras > 1

def wise_guys():
    """
    Creates JJH plots for WISE and Aspin sources.
    
    """

    # wise disks
    for s, n in zip(wise_disks[rwd], wise_disks_names[rwd]):
        fig = plot3.jjh(watso, s, name=str(n), color_slope=True, 
                        date_offset=54579,
                        outfile=path2+'Wise/'+str(n), png_too=True)

        if fig == None:
            print "dude %s failed to plot right" % str(s)

    # wise special case
    for s, n in zip(wise_v, wise_v_names):
        fig = plot3.jjh(data, s, name=str(n), color_slope=True, 
                        date_offset=54579,
                        outfile=path2+'Wise/'+str(n), png_too=True)


    # aspins
    for s, n in zip(aspin_sources[ras], aspin_sources_names[ras]):
        fig = plot3.jjh(watso, s, name=str(n), color_slope=True, 
                        date_offset=54579,
                        outfile=path2+'Aspin/'+str(n), png_too=True)

        if fig == None:
            print "dude %s failed to plot right" % str(s)


    # wise trans
    for s, n in zip(wise_trans[rwt], wise_trans_names[rwt]):
        fig = plot3.jjh(watso, s, name=str(n), color_slope=True, 
                        date_offset=54579,
                        outfile=path2+'transition/'+str(n), png_too=True)

        if fig == None:
            print "dude %s failed to plot right" % str(s)

    # wise extras
    for s, n in zip(wise_extras[rwe], wise_extras_names[rwe]):
        fig = plot3.jjh(data, s, name=str(n), color_slope=True, 
                        date_offset=54579,
                        outfile=path2+'Wise_extras/'+str(n), png_too=True)

        if fig == None:
            fig2 = plot3.jjh(welo, s, name=str(n), color_slope=True, 
                             date_offset=54579,
                             outfile=path2+'Wise_extras/'+str(n), png_too=True)
            if fig2 == None:
                print "dude %s failed to plot right" % str(s)

    print "did that."

def rwa_guys():
    """
    Creates JJH plots for RWA sources.

    """
    
    print "Now for RWA sources:"

    for s, n in zip(rwa_sources, rwa_names):
        fig = plot3.jjh(rwa_data, s, name=str(n), color_slope=True, 
                        date_offset=54579,
                        outfile=path2+'RWA_sources/'+str(n), png_too=True)

        if fig == None:
            print "dude %s failed to plot right" % str(s)


def tables():
    """
    Creates spreadsheets for Wise, Aspin, wise_trans, wise_extras, 
    and rwa_sources that include color-slope terms and robust stats.

    """

    # wise disks
    
    wise_disks_lookup = atpy.Table()
    wise_disks_lookup.add_column("SOURCEID", wise_disks[rwd])
    wise_disks_lookup.add_column("Designation", wise_disks_names[rwd])

    spread3.spreadsheet_write(
        watso, wise_disks_lookup, 0, 
        path2+"Wise/spreadsheet_withcolors.fits",
        flags=256, rob=True, colorslope=True, per=True)

    # wise special case
    wise_v_lookup = atpy.Table()
    wise_v_lookup.add_column("SOURCEID", wise_v)
    wise_v_lookup.add_column("Designation", wise_v_names)

    spread3.spreadsheet_write(
        data, wise_v_lookup, 0, 
        path2+"Wise/spreadsheet_withcolors_wise29434.fits",
        flags=256, rob=True, colorslope=True, per=True)
    
    # aspins
    
    aspin_sources_lookup = atpy.Table()
    aspin_sources_lookup.add_column("SOURCEID", aspin_sources[ras])
    aspin_sources_lookup.add_column("Designation", aspin_sources_names[ras])

    spread3.spreadsheet_write(
        watso, aspin_sources_lookup, 0, 
        path2+"Aspin/spreadsheet_withcolors.fits",
        flags=256, rob=True, colorslope=True, per=True)
    
    # wise trans

    wise_trans_lookup = atpy.Table()
    wise_trans_lookup.add_column("SOURCEID", wise_trans[rwt])
    wise_trans_lookup.add_column("Designation", wise_trans_names[rwt])

    spread3.spreadsheet_write(
        watso, wise_trans_lookup, 0,
        path2+"transition/spreadsheet_withcolors.fits",
        flags=256, rob=True, colorslope=True, per=True)

    # wise extras

    wise_extras_lookup = atpy.Table()
    wise_extras_lookup.add_column("SOURCEID", wise_extras[rwe])
    wise_extras_lookup.add_column("Designation", wise_extras_names[rwe])

    spread3.spreadsheet_write(
        data, wise_extras_lookup, 0,
        path2+"Wise_extras/spreadsheet_withcolors.fits",
        flags=256, rob=True, colorslope=True)

    # special case of whatevers - this is gonna be annoying to disentangle.
    # I'll have to see what the error behavior of statcruncher is when
    # there's no data for a star. Hopefully null values and no crashing.

    welo_lookup = atpy.Table()
    welo_lookup.add_column("SOURCEID", welo_sources)
    welo_lookup.add_column("Designation", welo_names)
    
    spread3.spreadsheet_write(
        welo, welo_lookup, 0,
        path2+"Wise_extras/spreadsheet_withcolors_lonely.fits",
        flags=256, rob=True, colorslope=True)
    
    # rwa dudes

    rwa_lookup = atpy.Table()
    rwa_lookup.add_column("SOURCEID", rwa_sources)
    rwa_lookup.add_column("Designation", rwa_names)

    spread3.spreadsheet_write(
        rwa_data, rwa_lookup, 0,
        path2+"RWA_sources/spreadsheet_withcolors.fits",
        flags=256, rob=True, colorslope=True, per=True)
    
    print "Did the spreadsheets!"
