'''
Generic, extensible code to create and use the photometric corrections
for the constant network in my WFCAM data.
Does so on a per-chip basis, i.e. each of the 16 chips has its 
own constant network.

Replaces network.py and will eventually be renamed network.py 
after the smoke clears.
'''

import os

import numpy as np
import robust as rb

import atpy
import pickle
from tr_helpers import data_cut, season_cut

# Loading up some good default data, which can be completely ignored
# with proper keyword usage. (I don't anticipate this happening, but...)

script_location = os.path.dirname(os.path.realpath(__file__))
path = script_location+'/../network/'

dates = np.loadtxt( path + "sorted_date_list.txt")

f = open(path+"Centers")
center_map = pickle.load(f)
f.close()

f2 = open(path+"Quadrants")
chip_quadrant_map = pickle.load(f2)
f2.close()


# This function works!
def get_chip( date, ra, dec, 
              dates = dates, center_map = center_map, 
              chip_quadrant_map = chip_quadrant_map ) :
    ''' Gets a star's chip based only on the fundamental parameters.
    
    Inputs:
      date -- timestamp of observation.
      ra -- position in degrees.
      dec -- position in degrees.
      
    Keywords:
      dates -- a numpy array of valid dates to use.
      center_map -- a dictionary mapping chipsets to coordinates.
      chip_quadrant_map -- a dictionary mapping chipset + location information
                           to specific chips.

    '''
    # I might move these back outside the function definition (for speed),
    # but then it breaks the portability of this code.
    which_chipset = np.arange(dates.size) % 4
    chipset_map = dict( zip( dates.astype('|S14'), which_chipset ) )

    if date not in dates:
        print "Faulty date"
        return -1

    # First, determine which chipset we're on based on the timestamp
    chipset = chipset_map[ str(date) ]
    
    # Then, use the chipset and the star's position to pinpoint 
    # which chip we're on
    center = center_map[chipset]
    
    chip = chip_quadrant_map[ ( chipset, 
                                ra > center[0], 
                                dec > center[1] ) ]
    
    return chip


# This function works
def feed_chip_one ( table, sid ):
    ''' Gets a star's chip intelligently based on its SOURCEID and table.

    Inputs:
      table -- an ATpy table with time-series data and positions.
      sid -- a 13-digit WFCAM source ID.

    Calls get_chip. A convenience function. Uses default get_chip keywords.
    Only looks up a single event's chip. 
    SEE ALSO: feed_chip_four if you want more.
    '''
    s_table = season_cut(table, sid, 123)
    date = s_table.MEANMJDOBS
    RA = np.degrees(s_table.RA)
    DEC= np.degrees(s_table.DEC)
    
    return get_chip( date[0], RA[0], DEC[0] )

# Untested but probably works
def feed_chip_four ( table, sid ):
    ''' Gets a star's chip intelligently based on its SOURCEID and table.
    Checks to see if the star is on multiple chips.

    Inputs:
      table -- an ATpy table with time-series data and positions.
      sid -- a 13-digit WFCAM source ID.

    Calls get_chip. A convenience function. Uses default get_chip keywords.
    Looks up a star's chip for its first four detections,
    to see if the star is on multiple chips.

    SEE ALSO: feed_chip_one for less.

    I'll have to think carefully about what to return here. 
    (or whether this is even useful or necessary)
    '''
    s_table = season_cut(table, sid, 123)
    date = s_table.MEANMJDOBS
    RA = np.degrees(s_table.RA)
    DEC= np.degrees(s_table.DEC)
    
    chips = [ get_chip( date[x], RA[x], DEC[x] ) for x in range(4) ]

    # chips = [ get_chip( date[0], RA[0], DEC[0] ),
    #           get_chip( date[1], RA[1], DEC[1] ),
    #           get_chip( date[2], RA[2], DEC[2] ),
    #           get_chip( date[3], RA[3], DEC[3] ) ]

    # For now I'll just return the four chips in a list, rather than some kind
    # of flag indicating whether the star was on multiple chips. This seems the
    # most naive way to handle it.
    return chips


# Untested but based on what I did interactively that works
def make_perfect_withchips_table ( stats_table ):
    ''' Filters a stats table for perfect sources and adds a chip column.'''

    ws = stats_table
    wsh =  ws.where( (ws.n_detect > 110) & (ws.n_detect < 130) & 
                     (ws.jpp_max == 0) & (ws.kpp_max ==0) & (ws.hpp_max == 0) )

    chip_arr = np.zeros_like( wsh.j_rms )
    
    for sid, i in zip(wsh.SOURCEID, arange(chip_arr.size)):
        chip_arr[i] = feed_chip_one(w, sid)


# This function works
def make_constants_table ( stats_table, selection_column='j_rms', num=10 ):
    ''' Makes a table of the best constant stars per chip.

    Inputs:
      stats_table -- an ATpy table with RMS and chip information for each star.

    Optional inputs:
      selection_column -- the criteria to choose the "best" stars from.
      num -- how many stars per chip to choose.

    Returns:
      An ATpy table of 'num' most constant stars per chip, with the same 
      columns as stats_table.
      '''
    
    # For each chip, make a copy of the j_rms column, sort it, and 
    # do a "where" query

    # no, do a lexsort! it makes more sense
    
    # sc is the data in the chosen selection column (usually j_rms)
    sc = stats_table[selection_column]
    ids= stats_table.SOURCEID

    constants_list = []
    
    chip_list = list( set( stats_table.chip ) )

    for chip in chip_list:
        # Make a copy of the j_rms and SOURCEID
        local_sc = sc[stats_table.chip == chip].copy()
        local_ids= ids[stats_table.chip == chip].copy()

        # then sort the two together
        ind = np.lexsort( ( local_ids, local_sc ) )
        
        # and save the 'num' best from this chip
        bests = local_ids[ind][:num]
        constants_list.append( bests )

    # finally, concatenate them all and lookup the proper rows in stats_table
    constant_ids = np.concatenate( constants_list )
    constant_table = stats_table.where(
        np.array( [sid in constant_ids for sid in stats_table.SOURCEID] )
        )

    print "constant table has %d rows!" % constant_table.shape[0]
    return constant_table


# This function works!
def make_corrections_table ( constants, table ):
    ''' Creates a table of photometric corrections per chip per night.

    Inputs:
      constants -- an ATpy table which gives 10 constant stars per chip.
                   Columns: "SOURCEID" (13-digit int), "chip" (1-16 int)
      table -- an ATpy table with time-series photometry

    Returns:
      an ATpy table with the following format:
      THE CORRECTIONS TABLE:
      night            chip   correction_J  corr_H   corr_K
      '54582.6251067'  3      +0.13         +0.07    -0.03

    '''
    # rb.meanr(x) is the robust mean

    # First - let's compute every constant star's robust mean in each band.
    # And keep track of them.

    j_meanr = np.zeros(constants.SOURCEID.size) * 1.
    h_meanr = np.zeros(constants.SOURCEID.size) * 1.
    k_meanr = np.zeros(constants.SOURCEID.size) * 1.
    
    for sid, i in zip(constants.SOURCEID, range(constants.SOURCEID.size)):
        stable = season_cut(table, sid, 123, flags=0)
        
        j_meanr[i] = rb.meanr( stable.JAPERMAG3 )
        h_meanr[i] = rb.meanr( stable.HAPERMAG3 )
        k_meanr[i] = rb.meanr( stable.KAPERMAG3 )
        del stable
    
    try:
        constants.add_column('j_meanr', j_meanr)
        constants.add_column('h_meanr', h_meanr)
        constants.add_column('k_meanr', k_meanr)

        print "computed robust-mean magnitude for each constant star"
    except:
        print "looks like you already computed robust-mean magnitudes"
        
    # Second - Calculate mean(r) deviations for each chip for each night
    chip_list = list( set( constants.chip ) )
    
    corrections_list = [] # add tables to this list, join them up at the end

    for chip in chip_list:
        
        local_network = constants.where(constants.chip == chip)
        
        # so I'm taking all of the local stars, and for every night...
        # do i calculate it by each star first, or each night first?
        # each night first I think would work better.

        #let's grab a slice of the big table corresponding only to our 
        # favorite sources' photometry. 

        # by joining together the season_cuts from all 10 sources!
        # no that's stupid. Use an "or |" operator! this is gonna be painful

        cids= local_network.SOURCEID
        ids = table.SOURCEID
        #local_table = season_cut( table, local_network[0], 123 )

        local_table = data_cut( table, cids, 123, flags=0 ) # aww yeah

        # local_table = table.where( ( (ids == cids[0]) | # I really, really wish
        #                              (ids == cids[1]) | # I knew how to make
        #                              (ids == cids[2]) | # this more elegant.
        #                              (ids == cids[3]) | 
        #                              (ids == cids[4]) |
        #                              (ids == cids[5]) |
        #                              (ids == cids[6]) |
        #                              (ids == cids[7]) |
        #                              (ids == cids[8]) |
        #                              (ids == cids[9]) ) &
        #                            (table.JPPERRBITS <= 0) &
        #                            (table.HPPERRBITS <= 0) &
        #                            (table.KPPERRBITS <= 0) )

        # okay, now that i've got the local table... let's get each night's
        # meanr deviation.
        # first, let's make some dates to iterate through
        date_list = list( set( local_table.MEANMJDOBS ) )
        
        # at some point i need to make a structure to save the corrections to

        ld = len(date_list)
        
        date_arr = np.zeros(ld)
        j_correction = np.zeros(ld)
        h_correction = np.zeros(ld)
        k_correction = np.zeros(ld)
        chip_arr = chip * np.ones(ld, dtype=int)

        # get each night's correction!
        for date, j in zip( date_list, range(ld) ):
            
            # a temporary place to keep the individual deviations
            j_deviation = np.zeros_like(cids) * 1.
            h_deviation = np.zeros_like(cids) * 1.
            k_deviation = np.zeros_like(cids) * 1.

            for star, i in zip(cids, range(cids.size) ):
                star_night_row = local_table.where( 
                    (local_table.SOURCEID == star) &
                    (local_table.MEANMJDOBS == date) )

                # deviation: the meanr minus that night's magnitude.
                j_deviation[i] = (constants.j_meanr[constants.SOURCEID==star]- 
                                  star_night_row.JAPERMAG3 )
                h_deviation[i] = (constants.h_meanr[constants.SOURCEID==star]- 
                                  star_night_row.HAPERMAG3 )
                k_deviation[i] = (constants.k_meanr[constants.SOURCEID==star]- 
                                  star_night_row.KAPERMAG3 )

            date_arr[j] = date
            j_correction[j] = -rb.meanr(j_deviation)
            h_correction[j] = -rb.meanr(h_deviation)
            k_correction[j] = -rb.meanr(k_deviation)


        # make a table for each chip, and (at the end) 
        # add it to the corrections_list. We'll join them up at the end.

        correction_subtable = atpy.Table(name="The Corrections Table")
        # add_column( 'name', data )
        correction_subtable.add_column('date', date_arr)
        correction_subtable.add_column('chip', chip_arr)
        correction_subtable.add_column('j_correction', j_correction)
        correction_subtable.add_column('h_correction', h_correction)
        correction_subtable.add_column('k_correction', k_correction)

        corrections_list.append( correction_subtable )

    #whoo, finally!

    correction_table = corrections_list[0]
    for subtable in corrections_list[1:] :
        correction_table.append( subtable )

    return correction_table

    ''' 
    To compute the corrections table:
      -Define 10 "good" stars per chip based on lowest j_rms 
      (from a sample of perfect sources)
      -Calculate each goodstar's mean magnitude. Keep track of it.
      (Probably the robust mean?)
     FOR EACH chip, FOR EACH night,
      -Calculate each of the 10 stars' deviation, average them (robust mean!) 
      and save its negative as the "correction".

    THE CORRECTIONS TABLE:
    night            chip   correction_J  corr_H   corr_K
    '54582.6251067'  3      +0.13         +0.07    -0.03
                         .
                         .
                         .
      '''
    



def apply_corrections( data_table, corrections_table ):
    ''' the real docstring comes later '''

    ''' 
    To apply corrections:
    FOR EACH star, FOR EACH night (timestamp),
      -Lookup chip for that timestamp (using get_chip),
      -Lookup appropriate correction, 
       add it to photometry, 
       and save the new magnitude.
      -Include a correction to the uncertainty:
       sigma_new = sqrt( sigma_old**2 + correction**2 )
       "that's not exactly right, but it's good enough"

     this might take forever.
    '''

    pass
