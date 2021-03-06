import atpy
import numpy as np
import matplotlib.pyplot as plt
import tr

datapath = '/home/tom/reu/DATA/Merged_Catalogs/'
savepath = '/home/tom/reu/logs/feb_15_2011/huge/'

# Let's load up our data...

wc = atpy.Table(datapath + 'wserv_errbits_combined.fits')

goods = atpy.Table(datapath + '5_sigma_disks_good1.fits')

for row in goods:
    desig = row['Designation']
    name = "p%d.%d" % ( row[4], round(row[3]*100, 0) )
    sid = row[0]
    
    for season in [1,2,3]:
        outname = savepath+ ("%.3d" % desig) + "." + str(season) + ".pdf"
        outname_png = ( savepath+"png/"+ ("%.3d" % desig) + 
                        "." + str(season) + ".png" )
        tr.plot_5(wc, sid, season=season, name=name, outfile=outname)
        tr.plot_5(wc, sid, season=season, name=name, outfile=outname_png)

    print name + " is done"


