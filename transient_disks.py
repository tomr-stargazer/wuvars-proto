''' Now I'm gonna search for these disk guys '''

import numpy as np
import atpy
import matplotlib.pyplot as plt
from tr import make_sidset

datapath = '/home/tom/reu/DATA/Merged_Catalogs/'
savepath = '/home/tom/reu/logs/dec_10_2010/movies/'

# Let's load up our data...

#wc = atpy.Table(datapath + 'wserv_errbits_combined.fits')

def trand(wc):

# First I'm gonna generate a list of dates:


    dates = np.array( sorted( list( set( wc.MEANMJDOBS ) ) ) )

    print dates.size
    
    sids = make_sidset(wc)

    disk_counter = atpy.Table()
    disk_counter.add_column("SOURCEID", sids)
    disk_counter.add_column("n_disks", np.zeros_like(sids))
    
    disk_set = set()



    for date, i in zip(dates, range(dates.size)):
        
        print date, i
        # Let's do our awesome filtering
        
        tab = wc.where(wc.MEANMJDOBS == date)
        
        tab = tab.where(tab.PSTAR > 0.8)
        
        tab = tab.where(tab.JAPERMAG3ERR < .02)
        tab = tab.where(tab.JAPERMAG3ERR > 0)
        tab = tab.where(tab.HAPERMAG3ERR < .02)
        tab = tab.where(tab.HAPERMAG3ERR > 0)
        tab = tab.where(tab.KAPERMAG3ERR < .02)
        tab = tab.where(tab.KAPERMAG3ERR > 0)

        tab = tab.where(tab.JPPERRBITS < 5)
        tab = tab.where(tab.HPPERRBITS < 5)
        tab = tab.where(tab.KPPERRBITS < 5)
        
        disks = tab.where(tab.JMHPNT < 1.71*tab.HMKPNT - .08)
        
        for d in disks.SOURCEID:
            disk_set.add(d)
            
            inc_col = disk_counter.where(disk_counter.SOURCEID == d).n_disks 
            inc_col += 1
            
        plt.clf()
        plt.plot( tab.HMKPNT, tab.JMHPNT, 'k.')
        plt.plot( disks.HMKPNT, disks.JMHPNT, 'ro')

        plt.title( str(date) )
        
        plt.xlim(0,2.5)
        plt.ylim(0,3.5)
        plt.savefig( savepath + str(i % 4) + ("/%.4d" % i) + ".png" )

    return disk_set, disk_counter
