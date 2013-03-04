"""

This is a script for re-generating five-panel lightcurves 
(J, H, K, KvH-K, J-HvH-K) in Cyg OB7 with the following adjustments/upgrades:

1. All plots have fitted color slope lines
2. A consistent "zero date" of 04/23/2008 across all plots
3. Period-folded plots only print two decimal places in the x-label

for the following groups of stars:

1. all 30 RWA stars
2. Wise disk sources
3. Wise transition disk sources
4. Wise "extra" sources
5. Aspin sources

"""

import numpy as np
import matplotlib.pyplot as plt

import atpy

import plot3, spread3

from cygob7_jjh_script_helper import *

# where the data is stored
path = '/home/tom/Dropbox/Cyg_OB7/paper2/'
path2 ='/home/tom/Dropbox/Cyg_OB7/paper2/book/March2013/' 
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

def do_it_all_cygob7():
    """
    Creates Lightcurve plots for WISE, Aspin, and RWA sources.

    For each (season), (period-type), (subcategory of stars), do stuff!
    
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

    print "Now for RWA sources:"

    for s, n in zip(rwa_sources, rwa_names):
        fig = plot3.jjh(rwa_data, s, name=str(n), color_slope=True, 
                        date_offset=54579,
                        outfile=path2+'RWA_sources/'+str(n), png_too=True)

        if fig == None:
            print "dude %s failed to plot right" % str(s)


    print "did that."

    


