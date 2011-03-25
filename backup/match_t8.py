''' Here is my first attempt at consecutive matching, to make my first 
lightcurves. Let's stick to maybe 100 sources of reliable magnitude...
and think of creative ways to store the data that glues the catalogs together.

I am a fan of simply marking the index where the match is, and then grabbing 
the data at that index point when we need it later.

Oh! I've suddenly been enlightened. Make the matching procedure a simple 
function! Let's attempt to define the function here and then, if it's a
success, we'll make it official somewhere.
'''

import atpy 
import numpy, math
import matplotlib.pyplot as plt
import coords

where = numpy.where
sect = numpy.intersect1d

def gen_match ( table1, table2, ra1, dec1, ra2, dec2, min_match ):
    ''' Matches two catalogs and returns a tuple of two numpy arrays.'''
    
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
    boxsize = min_match / 3600.

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

        #
        offset = -1. * numpy.ones_like(table2.radd[box])
        if offset.size != 0:
            for s2 in range(offset.size):
                p2 = coords.Position( (table2.radd[box][s2], table2.dedd[box][s2]),  units = 'deg')
                offset[s2] = p1.angsep(p2).arcsec()

            min_offset[s1] = offset.min()
            match[s1] = box[ where(offset == offset.min()) ][0]


    return (match, min_offset)


def smart_match(table1, table2, boxsize) : 
    try:
        gen_match(table1, table2, 'RA', 'DEC', 'RA', 'DEC', boxsize)
    except ValueError:
        gen_match(table1, table2, 'ra', 'dec', 'RA', 'DEC', boxsize)

where = numpy.where

boxsize_as = 1. #size, in arcseconds, that the best matches will come from
boxsize = boxsize_as/3600. # boxsize in degrees

print "box size for selection is %f arcsec" % boxsize_as

data = '/home/trice/reu/DATA/Merged_Catalogs/'
data2= '/home/trice/reu/DATA/2MASS/'

# Does ATpy have options to suppress output for these fits files?

Userv= atpy.Table(data+'USERV1678_error_filtered_0.1.fits',verbose=False)
#u08b = atpy.Table(data+'U08BH2_error_filtered_0.1.fits',verbose=False)
#u09b = atpy.Table(data+'U09BH14_error_filtered_0.1.fits',verbose=False)
Tmass= atpy.Table(data2+'fp_2mass.fp_psc26068.tbl',verbose=False)
#iras = atpy.Table(data2+'iras.iraspsc26148.tbl',verbose=False)

Output = atpy.Table() #Output is my output table!

delta = math.cos(math.radians(Tmass.dec.mean()))

#brights = numpy.where(Tmass.j_m < 10)[0]
#print brights

min_offset = -0.1* numpy.ones_like(Tmass.ra)
match =  -1* numpy.ones_like(Tmass.err_ang) #the only reason I use err_ang is because I want "match" to be an integer array, this might be irrelevant

counter = 1

for s1 in range(Tmass.ra.shape[0]):
    print "================ Source number %d: ================" % counter

    
    c1 = coords.Position((Tmass.ra[s1],Tmass.dec[s1]))
    print "Coordinates: ", c1
    #WFCAM is still in radians!!!
    # let's take the intersection of these four where queries
    w1= where(numpy.degrees(Userv.DEC) < Tmass.dec[s1] + boxsize)[0]
    w2= where(numpy.degrees(Userv.DEC) > Tmass.dec[s1] - boxsize)[0]
    w3= where(numpy.degrees(Userv.RA)  < Tmass.ra[s1] + boxsize/delta)[0]
    w4= where(numpy.degrees(Userv.RA)  > Tmass.ra[s1] - boxsize/delta)[0]

    rabox = numpy.intersect1d(w3,w4)
    decbox= numpy.intersect1d(w1,w2)

    box = numpy.intersect1d(rabox,decbox)
    
#    print box, "these guys got matched"

    offset = -1.* numpy.ones_like(Userv.RA[box])
    if offset.size != 0:
        for s2 in range(Userv.RA[box].shape[0]):
            # Don't forget! WFCAM data comes in radians!
            c2 = coords.Position((Userv.RA[box][s2],Userv.DEC[box][s2]),units='radians')
            offset[s2] = c1.angsep(c2).arcsec()
        
    #print "The min,max,mean offset: ", str(offset.min()), str(offset.max()), str(offset.mean()) #all in arcsec

        min_offset[s1] = offset.min()

        print "Offset to best match: %f arcsec" % min_offset[s1]
        match[s1] = box[numpy.where(offset == offset.min())][0]
        print "Index number of catalog 2 source match: %s" % str(match[s1])
#numpy.where(offset.min())
    #print Userv.data[match]
    #try also Userv.row(match)

        print "2MASS jmag: "+str(Tmass.j_m[s1])
        print "WFCAM jmag: "+str(Userv.JAPERMAG3[match[s1]])
    else:
        print "No match found within %d arcsec" % boxsize_as
    counter += 1
    #if counter == 100: break

print "finished successfully"

Output.add_column('match',match)
Output.add_column('offset',min_offset,unit='arcsec')


Output.write('test_output.fits',overwrite=True)

plt.hist(numpy.histogram(min_offset)[0])
plt.show()
