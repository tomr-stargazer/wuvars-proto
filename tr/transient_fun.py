import numpy as np
import atpy
import matplotlib.pyplot as plt
from spread3 import make_sidset


def disk_analyzer ( wc, cutoff ):
    ''' 
    A function that analyzes the transient-ness of disks.
    Returns a table of disk statistics


    Inputs: 

    wc: An atpy table
    cutoff: disk criteria cutoff (e.g. 0.1)
    '''

    dates = np.array( list( set( wc.MEANMJDOBS ) ) )

    sids = make_sidset(wc)

    disk_counter = atpy.Table()
    disk_counter.add_column("SOURCEID", sids)
    disk_counter.add_column("n_disks", np.zeros_like(sids))
    disk_counter.add_column("n_nodisks", np.zeros_like(sids))
    
    disk_set = set()

# First, we're gonna see how many stars even have disks and how many

    for date, i in zip(dates, range(dates.size)):
        
        # Let's do our awesome filtering

        tab = wc.where(wc.MEANMJDOBS == date)

        tab = tab.where(tab.PSTAR > 0.8)

        tab = tab.where(tab.JAPERMAG3ERR < .05)
        tab = tab.where(tab.JAPERMAG3ERR > 0)
        tab = tab.where(tab.HAPERMAG3ERR < .05)
        tab = tab.where(tab.HAPERMAG3ERR > 0)
        tab = tab.where(tab.KAPERMAG3ERR < .05)
        tab = tab.where(tab.KAPERMAG3ERR > 0)
        
        tab = tab.where(tab.JPPERRBITS < 5)
        tab = tab.where(tab.HPPERRBITS < 5)
        tab = tab.where(tab.KPPERRBITS < 5)

        disks = tab.where(tab.JMHPNT < 1.71*tab.HMKPNT - .08).SOURCEID

        for d in disks:
            disk_set.add(d)
            # therefore disk_set is a set of disk SOURCEIDs
            disk_counter.n_disks[np.where (disk_counter.SOURCEID == d)] += 1


##    plt.plot( tab.HMKPNT, tab.JMHPNT, 'k.')
##    plt.savefig( savepath + ("%.4d" % i) + ".png" )



# Second, we're gonna see (for stars that DO have disks) how often they don't 
# have a disk. This is done after we populate our disk_set.

# (The reason we're not doing this at the same time as the other one... 
# it would probably take too long to execute)

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

        nodisks = tab.where(tab.JMHPNT > 1.71*tab.HMKPNT - .08).SOURCEID

        for d in disk_set:
            if d in nodisks:
                disk_counter.n_nodisks[np.where(disk_counter.SOURCEID==d)] += 1


            place = np.where(disk_counter.SOURCEID == d)


    useful_disk_counter = disk_counter.where( disk_counter.n_disks > 0 )

# Finally, we'll return something we care about.

    return disk_set, useful_disk_counter



# What should I do next? For all the stars that DO ever have disks,
# how many nights do they NOT have disks.

##    plt.clf()
##    plt.hist(disk_counter.n_disks, 30)
##    plt.show()

