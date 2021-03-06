''' Now I'm gonna search for these disk guys '''

import numpy as np
import atpy
import matplotlib.pyplot as plt
from tr import make_sidset

datapath = '/home/tom/reu/DATA/Merged_Catalogs/'
#savepath = '/home/tom/reu/logs/dec_10_2010/movies/'
savepath = '/home/tom/reu/logs/jan_5_2011/movies/'

# Let's load up our data...

wc = atpy.Table(datapath + 'wserv_errbits_combined.fits')


# First I'm gonna generate a list of dates:

dates = np.array( list( set( wc.MEANMJDOBS ) ) )

sids = make_sidset(wc)

disk_counter = atpy.Table()
disk_counter.add_column("SOURCEID", sids)
disk_counter.add_column("n_disks", np.zeros_like(sids))
disk_counter.add_column("n_nodisks", np.zeros_like(sids))

disk_set = set()



for date, i in zip(dates, range(dates.size)):
    
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

    disks = tab.where(tab.JMHPNT < 1.71*tab.HMKPNT - .08).SOURCEID

    for d in disks:
        disk_set.add(d)

# This makes copies of the table! Erggghhh        
#        disk_counter.where(disk_counter.SOURCEID == d).n_disks += 1
        disk_counter.n_disks[np.where (disk_counter.SOURCEID == d)] += 1


    plt.plot( tab.HMKPNT, tab.JMHPNT, 'k.')
    plt.savefig( savepath + ("%.4d" % i) + ".png" )

# What should I do next? For all the stars that DO ever have disks,
# how many nights do they NOT have disks.

plt.clf()
plt.hist(disk_counter.n_disks, 30)
plt.show()
