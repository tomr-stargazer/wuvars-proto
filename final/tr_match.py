''' A collection of functions for matching tables based on source positions.'''

import atpy 
import numpy, math
import matplotlib.pyplot as plt
import coords

where = numpy.where
sect = numpy.intersect1d

v = False

def vprint ( string ):
    '''Prints a string only if verbose is true'''
    if v:
        print string

def core_match ( radd1, dedd1, radd2, dedd2, max_match, verbose = True ) :
    ''' Matches two sets of position arrays and returns (two numpy arrays).

    Assumes all positions are decimal degrees.
    '''
    
    global v
    v = verbose

    delta = math.cos(math.radians(numpy.abs(dedd1).max()))
    boxsize = max_match / 3600.

    min_offset = -0.1 * numpy.ones_like(radd1)
    match      = -1   * numpy.ones_like(radd1).astype(int)

    counter = 1

    for s1 in range(radd1.size):
        p1 = coords.Position( (radd1[s1], dedd1[s1]), \
                                 units = 'deg')

        w1 = where(radd2 < radd1[s1] + boxsize/delta)[0]
        w2 = where(radd2 > radd1[s1] - boxsize/delta)[0]
        w3 = where(dedd2 < dedd1[s1] + boxsize)[0]
        w4 = where(dedd2 > dedd1[s1] - boxsize)[0]

        # Let's slice a box around our source
        box = sect(sect(w1,w2),sect(w3,w4))

        # And calculate offsets to all sources inside that box
        offset = -1. * numpy.ones_like(radd2[box])
        if offset.size != 0:
            for s2 in range(offset.size):
                p2 = coords.Position( (radd2[box][s2], dedd2[box][s2]),  units = 'deg')
                offset[s2] = p1.angsep(p2).arcsec()

            # If the closest match is within our matching circle
            if offset.min() < max_match:
                min_offset[s1] = offset.min()
                match[s1] = box[ where(offset == offset.min()) ][0]
                vprint( "Source %d: Matched with %f arcsec" \
                        % (counter, offset.min() ) )
            else:
                vprint( "Source %d: Failed to match" % counter)
        else:
            vprint( "Source %d: Failed to match" % counter)
        counter += 1

    return (match, min_offset)



def gen_match ( table1, table2, ra1, dec1, ra2, dec2, max_match, verbose=True):
    ''' Matches two catalogs and returns a tuple of two numpy arrays.

    Arguments:
    (atpy.Table, atpy.Table, string, string, string, string, float)
    
    Two ATpy tables to cross-match, the string names (respectively) of 
    Table 1's RA, Table 1's Dec, Table 2's RA, Table 2's Dec,
    the upper bound on source-matches in arcseconds 
    
    
    Keywords:
    verbose -- print running progress? (default True)
    '''

    global v
    v = verbose

    # Let's cast all position arrays to decimal degrees
    if ('deg' or 'DEG') in table1.columns[ra1].unit:
        table1.radd = table1.data[ra1]
        table1.dedd = table1.data[dec1]
    elif ('rad' or 'RAD') in table1.columns[ra1].unit:
        table1.radd = numpy.degrees(table1.data[ra1])
        table1.dedd = numpy.degrees(table1.data[dec1])

    if ('deg' or 'DEG') in table2.columns[ra2].unit:
        table2.radd = table2.data[ra2]
        table2.dedd = table2.data[dec2]
    elif ('rad' or 'RAD') in table2.columns[ra2].unit:
        table2.radd = numpy.degrees(table2.data[ra2])
        table2.dedd = numpy.degrees(table2.data[dec2])

    vprint('Matching two tables: ')
    vprint( table1.columns )
    vprint( table2.columns )
    
    return core_match (table1.radd, table1.dedd, table2.radd, table2.dedd,
                       max_match, verbose=verbose)

def smart_match(table1, table2, boxsize, verbose=True) : 
    ''' 
    Matches two catalogs and returns a tuple of two numpy arrays.

    Usage:
    (atpy.Table, atpy.Table, float)
    
    Two ATpy tables to cross-match, and the upper bound 
    on source-matches in arcseconds.
    This function guesses what strings to use for the column names.
    
    Keywords:
    verbose -- print running progress? (default True)
    '''

    global v
    v = verbose

    if 'RA' in table1.columns: 
        ra1 = 'RA'
        dec1= 'DEC'
    elif 'ra' in table1.columns:
        ra1 = 'ra'
        dec1= 'dec'

    if 'RA' in table2.columns: 
        ra2 = 'RA'
        dec2= 'DEC'
    elif 'ra' in table2.columns:
        ra2 = 'ra'
        dec2= 'dec'

    vprint("Table 1: %s, %s" % (ra1, dec1))
    vprint("Table 2: %s, %s" % (ra2, dec2))

    return gen_match (table1, table2, ra1, dec1, ra2, dec2, boxsize, 
                      verbose=verbose)

def small_match ( ra, dec, radd2, dedd2, max_match, verbose=True ) :
    ''' Matches one position in decimal degrees to a source in a table. 

    Negative return values indicate failure to match.
    '''

    radd1 = numpy.array([ra])
    dedd1 = numpy.array([dec])

    match, min_offset = core_match(radd1,dedd1,radd2,dedd2,max_match,verbose)
    return (match[0], min_offset[0])

def coords_match ( position, table, max_match = 10, verbose=True) :
    ''' Uses a coords.Position object and a table'''
    ra,dec = position.dd()
    radd2 = numpy.degrees(table.RA)
    dedd2 = numpy.degrees(table.DEC)
    return small_match(ra,dec,radd2,dedd2,max_match,verbose=True)
