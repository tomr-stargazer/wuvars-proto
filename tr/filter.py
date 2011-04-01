'''
A module to filter out data in various ways.
Currently in draft form.
'''

import numpy as np
import atpy
from tr_helpers import data_cut
import network2 as n2

# 1. A function that trims a table for "happy" sources, based on SOURCEID

def happy_table ( big_table, sid_list, flags=-1 ):
    ''' Basically calls data_cut in a special syntax that makes life easy. 
    Inputs:
      big_table -- a WFCAM time-series photometry table.
      sid_list -- a list or array of SOURCEIDs to cut along.
      
      '''

    return data_cut (big_table, sid_list, season=123, flags=flags)



# 2. A function that eliminates bad chip-nights from a dataset
# It works! (For a table of one sourceID) 
# Now, how long does it take for longer than that...
# Started at 4:37 ended at 4:39 - not that bad!! (for hh, the 473-sid table)
# How long does the whole wserv table take?? well ima go to sleep
# Almost exactly 5 hours, and uses up all my RAM so if I try to save
# anything then it crashes T_T
def happy_chipnights ( data_table, corrections_table, max_correction=.0125 ):
    ''' Removes bad chip-nights from a data table.
    
    Inputs:
      data_table -- an ATpy table with WFCAM time-series photometry.
      corrections_table -- an ATpy table with per-night per-chip corrections.
      
    Optional inputs:
      max_correction -- the minimum "chip deviation" to flag a chip_night as
                        needing removal.
    '''

    # So, we'll load up the table, and... for every datapoint we'll 
    # filter stuff? That might take forever!
    # Alternatively, we could check every source for whether it's on only
    # one chip

    # How do I even delete a row from an atpy table? 
    # Ah... I have a clever way to do it I think. 
    # The best way to do this seems to be to add the j,h,k corrections columns to this table, then populate them row-by-row, and filter it all at the end! Okay, I'll do that.
    
    ctable = corrections_table
    

    # I think I just learned how to copy tables! 
    # it uses really silly "where" calls like where something is positive or negative or zero
    if atpy.__version__ == '0.9.3':
        try:
            # deleting the columns we need to create in a sec
            data_table.remove_columns( ['j_correction',
                                        'h_correction',
                                        'k_correction'])
        except: 
            pass
        table = data_table
    else:
        # making a friggin copy so we can use it, and calling it table
        table = data_table.copy()
    
    table.add_empty_column('j_correction', np.float64)
    table.add_empty_column('h_correction', np.float64)
    table.add_empty_column('k_correction', np.float64)
   
    for row in data_table:
        chip = n2.get_chip( row['MEANMJDOBS'], 
                            np.degrees(row['RA']),
                            np.degrees(row['DEC']) )
        
        local_corrections = ctable.where( (ctable.date == row['MEANMJDOBS']) &
                                          (ctable.chip == chip) )
        
        row['j_correction'] = local_corrections.j_correction
        row['h_correction'] = local_corrections.h_correction
        row['k_correction'] = local_corrections.k_correction

    
    print "Made corrections rows! They look like this (j, h, k):"
    print table.j_correction[0]
    print table.h_correction[0]
    print table.k_correction[0]

    # Now I guess it's time to filter stuff... though the above seems 
    # pretty useful on its own... like, that's basically the implementation
    # for apply_corrections...
    
    f_table = table.where( (np.abs(table.j_correction) < max_correction) &
                           (np.abs(table.h_correction) < max_correction) &
                           (np.abs(table.k_correction) < max_correction) )

    print "done filtering!"
    
    return f_table
