''' 
A collection of functions for matching tables based on source positions.

Available functions:


Modification History:

June 2010: Created. Tom Rice (t.rice90@gmail.com)
17 March 2012: Added 'units' keyword to coords_match.
15 May 2012: Updating code.
17 May 2012: Still updating code, mainly by standardizing documentation.

'''

import atpy 
import numpy as np
#import math
import matplotlib.pyplot as plt
import coords

where = np.where
sect = np.intersect1d

v = False

def vprint ( string ):
    """
    Prints `string` only if the global variable `verbose` is True.
    """
    if v:
        print string

def core_match ( radd1, dedd1, radd2, dedd2, max_match, verbose = True ) :
    """ 
    Matches two sets of position arrays and returns (two numpy arrays).

    Assumes all positions are decimal degrees.

    Parameters
    ----------
    radd1, dedd1 : numpy arrays
        R.A. and Decl. arrays for the first table. In decimal degrees.
    radd2, dedd2 : numpy arrays
        R.A. and Decl. arrays for the second table. In decimal degrees.
    max_match : float
        Largest acceptable offset (in arcsec) for two sources to match.

    Returns
    -------
    match : numpy array
        Gives the index in table 2 corresponding to each item 
        in table 1. A value of -1 is a failed match.
    min_offset : numpy array
        Gives the distance between the best-matched items in 
        table 1 vs table 2.
    """
    
    global v
    v = verbose

    delta = np.cos(np.radians(np.abs(dedd1).max()))
    boxsize = max_match / 3600.

    min_offset = -0.1 * np.ones_like(radd1)
    match      = -1   * np.ones_like(radd1).astype(int)

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

#       HEY. dear anybody who ever edits this code:
#       DO NOT change the implementation of how "box" works. It works.
#       Leave it be. Do not attempt anything silly like the following.

#        box = np.array( (radd2 < radd1[s1] + boxsize/delta) &
#                        (radd2 > radd1[s1] - boxsize/delta) &
#                        (dedd2 < dedd1[s1] + boxsize) &
#                        (dedd2 > dedd1[s1] - boxsize))

        # And calculate offsets to all sources inside that box
        offset = -1. * np.ones_like(radd2[box])
        if offset.size != 0:
            for s2 in range(offset.size):
                p2 = coords.Position( (radd2[box][s2], dedd2[box][s2]),  units = 'deg')
                offset[s2] = p1.angsep(p2).arcsec()

            # If the closest match is within our matching circle
            if offset.min() < max_match:
                min_offset[s1] = offset.min()
                match[s1] = box[where(offset == offset.min() )][0]
                vprint( "Source %d: Matched with %f arcsec" \
                        % (counter, offset.min() ) )
            else:
                vprint( "Source %d: Failed to match" % counter)
        else:
            vprint( "Source %d: Failed to match" % counter)
        counter += 1

    return (match, min_offset)



def gen_match ( table1, table2, ra1, dec1, ra2, dec2, max_match, verbose=True):
    """ 
    Matches two catalogs and returns a tuple of two numpy arrays.

    Arguments:
    (atpy.Table, atpy.Table, string, string, string, string, float)
    
    Two ATpy tables to cross-match, the string names (respectively) of 
    Table 1's RA, Table 1's Dec, Table 2's RA, Table 2's Dec,
    the upper bound on source-matches in arcseconds 
    
    
    Keywords:
    verbose -- print running progress? (default True)
    """

    global v
    v = verbose

    # Let's cast all position arrays to decimal degrees.
    # Assume degrees by default.
    if ('deg' or 'DEG') in table1.columns[ra1].unit:
        table1.radd = table1.data[ra1]
        table1.dedd = table1.data[dec1]
    elif ('rad' or 'RAD') in table1.columns[ra1].unit:
        table1.radd = np.degrees(table1.data[ra1])
        table1.dedd = np.degrees(table1.data[dec1])
    else:
        table1.radd = table1.data[ra1]
        table1.dedd = table1.data[dec1]

    if ('deg' or 'DEG') in table2.columns[ra2].unit:
        table2.radd = table2.data[ra2]
        table2.dedd = table2.data[dec2]
    elif ('rad' or 'RAD') in table2.columns[ra2].unit:
        table2.radd = np.degrees(table2.data[ra2])
        table2.dedd = np.degrees(table2.data[dec2])
    else:
        table2.radd = table2.data[ra2]
        table2.dedd = table2.data[dec2]


    vprint('Matching two tables: ')
    vprint( table1.columns )
    vprint( table2.columns )
    
    return core_match (table1.radd, table1.dedd, table2.radd, table2.dedd,
                       max_match, verbose=verbose)

def smart_match(table1, table2, boxsize, verbose=True) : 
    """ 
    Matches two catalogs and returns a tuple of two numpy arrays.

    Usage:
    (atpy.Table, atpy.Table, float)
    
    Two ATpy tables to cross-match, and the upper bound 
    on source-matches in arcseconds.
    This function guesses what strings to use for the column names.
    
    Keywords:
    verbose -- print running progress? (default True)
    """

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
    """ 
    Matches one position in decimal degrees to a source in a table. 

    Negative return values indicate failure to match.
    """

    radd1 = np.array([ra])
    dedd1 = np.array([dec])

    match, min_offset = core_match(radd1,dedd1,radd2,dedd2,max_match,verbose)
    return (match[0], min_offset[0])

def coords_match ( position, table, max_match = 10, verbose=True, units='rad') :
    """ 
    Matches a coords.Position object to a table with coordinates.

    Parameters
    ----------
    position : coords.Position object
    table : ATpy table 
        Assumes the table has columns named "RA" and "DEC".
    max_match : int, optional
        The upper limit on match size, in arcseconds
        (defaults to 10).
    verbose : bool, optional
        Whether to print internal statements (defaults to True).
    units : str, optional
        What units to interpret the input Table in
        (defaults to 'rad' for radians).

    Returns
    -------

    """
    ra,dec = position.dd()
    if units=='rad':
        radd2 = np.degrees(table.RA)
        dedd2 = np.degrees(table.DEC)
    else:
        radd2 = table.RA
        dedd2 = table.DEC

        
    return small_match(ra,dec,radd2,dedd2,max_match,verbose=verbose)
