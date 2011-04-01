'''
A module to create my lightcurve and periodogram books in a repeatable way.

Replaces stuff like plot_1200.py etc.
May be extended to do other things too.
'''

import atpy
import numpy as np
import matplotlib.pyplot as plt
import plot as tplot
import tr

# Some ideas for things I'll want to do with Publisher:
# 1. Make my official, combined period+lightcurve books, with 
#    pages side-to-side so you can view it all two-up
#    (this one was obvious)
#
# 2. Make lightcurve books for my constant stars

# 3. Make comparison books for before-and-after data filtering.
#    This one might actually be kind of urgent.

# I should definitely try to automate the system calls to pdfjoin and such,
# though I should be careful to keep a low number of things joined together at 
# a time

# 3. The comparison before-and-after book:

def before_after ( before_table, after_table, lookup_table, savepath,
                   prefix = ''):
    ''' A function to visually compare many lightcurves 
    before and after corrections. 

    Inputs:
      before_table -- ATpy table with unfiltered photometry.
      after_table -- ATpy table with filtered photometry.
      lookup_table -- table with source and name information
      savepath -- directory to write output plots to.

    Optional inputs:
      prefix -- a small string to put in front of every name in the plots

    Plots all 3 seasons' worth of data, no segregating into seasons.
    '''

    for row in lookup_table:
        desig = row['Designation'] 
        # This only works when desig is a number!
        name_b = prefix + "%d_before" % desig 
        name_a = prefix + "%d_after" % desig
        sid = row['SOURCEID']

        # perhaps see the implementation below
        outname_a = savepath+ ("%.3d" % desig) + "." + 'a'
        outname_b = savepath+ ("%.3d" % desig) + "." + 'b'
        
        tplot.plot_5 ()
        return


def book_maker ( data_table, lookup_table, savepath, sort_column, 
                 max_i=-1, reverse=True ):
    ''' A function to make lightcurve-and-period books.
    
    Inputs:
      data_table -- an ATpy table with time-series photometry.
      lookup_table -- an ATpy table with star names, IDs, and other parameters
      savepath -- Where to save output files to
      sort_column -- (str) a column in lookup_table to sort the book by.

    Optional inputs:
      max_i -- how many stars to go through before stopping
               (defaults to all of them)
      reverse -- whether to iterate over the sorted column backwards
                 (defaults to True)
    '''
   
    # First, let's sort our lookup table by the sort column
    lookup_table.sort("%s" % sort_column)

    # This is SO HACKISH
    if reverse:
        look_table = reversed(lookup_table)
    else: 
        look_table = lookup_table

    # For each star
    for row, i in zip(look_table, np.arange(len(lookup_table))):
        
        desig = row['Designation']
        name = "p%d.%d" % ( desig, round(row['fraction']*100, 0) )
        sid = row['SOURCEID']
        # For each season

        for season in [0,1,2,3]:
            outname_a = ( savepath + ("%.3d_%.3d" % (i, desig)) + "." +
                          str(season) +"a")
            
            outname_b = ( savepath + ("%.3d_%.3d" % (i, desig)) + "." +
                          str(season) +"b")

            tr.plot_5 ( data_table, sid, outfile=outname_a, name=name, 
                        season=season, png_too=True )
            tr.plot_page_periods ( data_table, sid, outfile=outname_b, 
                                   name=name, season=season, png_too=True )
            
        print name + " is done"

        if i == max_i:
            break
