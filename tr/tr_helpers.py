'''
Here's some Helper functions.

This module does not call any of my other modules, for dependency reasons.
'''
import atpy
import numpy as np


def data_cut (table, sid_list, season, flags=-1):
    ''' Returns a subset of a table that corresponds to 
    one (or more!) source(s) for one (or more!) season(s). 

    Also filters datapoints for quality!

    Inputs:
      table -- an atpy table with time-series photometry
      sid_list -- a list of Source IDs from WFCAM (13 digits each)
      season -- Which observing season of our dataset (1, 2, 3, or all)

    Optional inputs:
      flags -- whether to remove bad observations (default: no)
               and where to draw the cutoff.
    '''

    # First, select these sources' data from the table
    source = table.where( 
        np.array( [sid in sid_list for sid in table.SOURCEID] )
        )
#        table.SOURCEID == sid)

    # If we had Numpy 1.4 or higher, we could use the following:
    # source = table.where( np.in1d(table.SOURCEID, sid_list) )


    # Next, figure out where to slice the data (in terms of dates)
    offset = 54579
    cut1 = 100
    cut2 = 300
    cut3 = 600

    if season == 1:
        low = offset
        high = offset+cut1
    elif season == 2:
        low = offset+cut1
        high = offset+cut2
    elif season == 3:
        low = offset+cut2
        high = offset+cut3
    else:
        low = offset
        high = offset+cut3
    
    source = source.where( (source.MEANMJDOBS < high) & 
                           (source.MEANMJDOBS > low))

    if flags >= 0:
        source = source.where( (source.JPPERRBITS <= flags) &
                               (source.HPPERRBITS <= flags) &
                               (source.KPPERRBITS <= flags) )
    
    return source





def season_cut (table, sid, season, flags=-1) :
    ''' Returns a subset of a table that corresponds to 
    one source for one (or more!) season(s). 

    Also filters datapoints for quality!

    Arguments:
      table: an atpy table with time-series photometry
      sid: a Source ID from WFCAM (13 digits)
      season: Which observing season of our dataset (1,2, or 3)

    Keywords:
      flags: whether to remove bad observations (default: no)
             and where to draw the cutoff.
    '''

    # First, select this source's data from the table
    source = table.where(table.SOURCEID == sid)

    # Next, figure out where to slice the data (in terms of dates)
    offset = 54579
    cut1 = 100
    cut2 = 300
    cut3 = 600

    if season == 1:
        low = offset
        high = offset+cut1
    elif season == 2:
        low = offset+cut1
        high = offset+cut2
    elif season == 3:
        low = offset+cut2
        high = offset+cut3
    else:
        low = offset
        high = offset+cut3
    
    source = source.where( (source.MEANMJDOBS < high) & 
                           (source.MEANMJDOBS > low))

    if flags >= 0:
        source = source.where( (source.JPPERRBITS <= flags) &
                               (source.HPPERRBITS <= flags) &
                               (source.KPPERRBITS <= flags) )
    
    return source

# I hope to deprecate this function someday.
def ensemble_cut (corr_table, chip, season, flags=-1) :
    ''' Returns a subset of a correction table that corresponds to 
    one chip for one (or more!) season(s). 


    Arguments:
      corr_table: an atpy table with CORRECTIONS to time-series photometry
      chip: a WFCAM chip (1-16)
      season: Which observing season of our dataset (1,2, or 3)

    '''

    # First, select this source's data from the table
    source = corr_table.where(corr_table.chip == chip)

    # Next, figure out where to slice the data (in terms of dates)
    offset = 54579
    cut1 = 100
    cut2 = 300
    cut3 = 600

    if season == 1:
        low = offset
        high = offset+cut1
    elif season == 2:
        low = offset+cut1
        high = offset+cut2
    elif season == 3:
        low = offset+cut2
        high = offset+cut3
    else:
        low = offset
        high = offset+cut3
    
    source = source.where( (source.date < high) & 
                           (source.date > low))
    
    return source

