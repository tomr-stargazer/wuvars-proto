""" 
This is the second iteration of transient_fun, and it uses the second
method that I think will work - where it tabulates the number of "i'm a disk"
and "i'm not a disk" for each source in the same run.

It's more elegant but I expect it to take much, much longer, unfortunately.

Godspeed yo.
"""

import numpy as np
import atpy
import matplotlib.pyplot as plt
from spread3 import make_sidset

def date_analyzer  ( table, cutoff=0.08 ):
    ''' 
    A function that analyzes the transient-ness of disks.
    Returns a table of date statistics: date, n_disks


    Inputs: 

    table: An atpy table
    cutoff: disk criteria cutoff (e.g. 0.1)
    '''

    dates = np.array( list( set( table.MEANMJDOBS ) ) )

    date_counter = atpy.Table()
    date_counter.add_column("date", dates)
    date_counter.add_column("n_disks", np.zeros_like(dates))

    for date, i in zip(dates, range(dates.size)):
        
        # Let's do our awesome filtering

        tab = table.where(table.MEANMJDOBS == date)

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

        disks = tab.where(tab.JMHPNT < 1.71*tab.HMKPNT - cutoff).SOURCEID
#        nodisks = tab.where(tab.JMHPNT > 1.71*tab.HMKPNT - cutoff).SOURCEID

        date_counter.n_disks[i] = disks.size

    return date_counter

################################


def disk_analyzer2 ( table, cutoff=0.08 ):
    ''' 
    A function that analyzes the transient-ness of disks.
    Returns a table of disk statistics


    Inputs: 

    table: An atpy table
    cutoff: disk criteria cutoff (e.g. 0.1)
    '''

    dates = np.array( list( set( table.MEANMJDOBS ) ) )

    sids = make_sidset(table)

    disk_counter = atpy.Table()
    disk_counter.add_column("SOURCEID", sids)
    disk_counter.add_column("n_disks", np.zeros_like(sids))
    disk_counter.add_column("n_nodisks", np.zeros_like(sids))
    
    disk_set = set()


    for date, i in zip(dates, range(dates.size)):
        
        # Let's do our awesome filtering

        tab = table.where(table.MEANMJDOBS == date)

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

        disks = tab.where(tab.JMHPNT < 1.71*tab.HMKPNT - cutoff).SOURCEID
        nodisks = tab.where(tab.JMHPNT > 1.71*tab.HMKPNT - cutoff).SOURCEID

# First, we're gonna see how many stars even have disks and how many

        for d in disks:
            disk_set.add(d)
            # therefore disk_set is a set of disk SOURCEIDs
            disk_counter.n_disks[np.where (disk_counter.SOURCEID == d)] += 1

# Second, we're gonna count up how often stars DON'T have
# a disk. This is done while we populate our disk_set.

        for nd in nodisks:
            disk_counter.n_nodisks[np.where(disk_counter.SOURCEID == nd)] += 1



# Second, we're gonna see (for stars that DO have disks) how often they don't 
# have a disk. This is done after we populate our disk_set.

# (The reason we're not doing this at the same time as the other one... 
# it would probably take too long to execute)



    useful_disk_counter = disk_counter.where( disk_counter.n_disks > 0 )

# Finally, we'll return something we care about.

    return disk_set, useful_disk_counter



# What should I do next? For all the stars that DO ever have disks,
# how many nights do they NOT have disks.

##    plt.clf()
##    plt.hist(disk_counter.n_disks, 30)
##    plt.show()

