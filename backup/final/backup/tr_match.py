import atpy 
import numpy, math
import matplotlib.pyplot as plt
import coords

where = numpy.where
sect = numpy.intersect1d

def vprint ( string, v ):
    '''Prints a string only if verbose is true'''
    if v:
        print string

def gen_match ( table1, table2, ra1, dec1, ra2, dec2, max_match, verbose=True):
    ''' Matches two catalogs and returns a tuple of two numpy arrays.'''

    # Let's cast all position arrays to decimal degrees
    if table1.columns[ra1].unit == 'deg':
        table1.radd = table1.data[ra1]
        table1.dedd = table1.data[dec1]
    elif table1.columns[ra1].unit == 'RADIANS':
        table1.radd = numpy.degrees(table1.data[ra1])
        table1.dedd = numpy.degrees(table1.data[dec1])

    if table2.columns[ra2].unit == 'deg':
        table2.radd = table2.data[ra2]
        table2.dedd = table2.data[dec2]
    elif table2.columns[ra2].unit == 'RADIANS':
        table2.radd = numpy.degrees(table2.data[ra2])
        table2.dedd = numpy.degrees(table2.data[dec2])
    

    delta = math.cos(math.radians(numpy.abs(table1.dedd).max()))
    boxsize = max_match / 3600.

    min_offset = -0.1 * numpy.ones_like(table1.radd)
    match      = -1   * numpy.ones_like(table1.radd).astype(int)

    for s1 in range(table1.radd.size):
        p1 = coords.Position( (table1.radd[s1], table1.dedd[s1]), \
                                 units = 'deg')

        w1 = where(table2.radd < table1.radd + boxsize/delta)[0]
        w2 = where(table2.radd > table1.radd - boxsize/delta)[0]
        w3 = where(table2.dedd < table1.dedd + boxsize)[0]
        w4 = where(table2.dedd > table1.dedd - boxsize)[0]

        # Let's slice a box around our source
        box = sect(sect(w1,w2),sect(w3,w4))

        # And calculate offsets to all sources inside that box
        offset = -1. * numpy.ones_like(table2.radd[box])
        if offset.size != 0:
            for s2 in range(offset.size):
                p2 = coords.Position( (table2.radd[box][s2], table2.dedd[box][s2]),  units = 'deg')
                offset[s2] = p1.angsep(p2).arcsec()

            # If the closest match is within our matching circle
            if offset.min() < max_match:
                min_offset[s1] = offset.min()
                match[s1] = box[ where(offset == offset.min()) ][0]

    return (match, min_offset)


def smart_match(table1, table2, boxsize) : 
    ''' Matches two catalogs and returns a tuple of two numpy arrays.'''
    if 'RA' in table1.columns: 
        ra1 = 'RA'
        dec1= 'DEC'
    elif 'ra' in table1.columns:
        ra1 = 'ra'
        dec1= 'dec'

    if 'RA' in table2.columns: 
        ra1 = 'RA'
        dec1= 'DEC'
    elif 'ra' in table2.columns:
        ra1 = 'ra'
        dec1= 'dec'

    return gen_match (table1, table2, ra1, dec1, ra2, dec2, boxsize)
    # try:
    #     gen_match(table1, table2, 'RA', 'DEC', 'RA', 'DEC', boxsize)
    # except ValueError:
    #     gen_match(table1, table2, 'ra', 'dec', 'RA', 'DEC', boxsize)
