''' This is the one in which I match the whole table!'''

import atpy 
import numpy, math
import matplotlib.pyplot as plt
import coords

where = numpy.where

boxsize_as = 5. #size, in arcseconds, that the best matches will come from
boxsize = boxsize_as/3600. # boxsize in degrees

print "box size for selection is %f arcsec" % boxsize_as

data = '/home/trice/reu/DATA/Merged_Catalogs/'
data2= '/home/trice/reu/DATA/2MASS/'

# Does ATpy have options to suppress output for these fits files?

userv= atpy.Table(data+'USERV1678_error_filtered_0.1.fits',verbose=False)
#u08b = atpy.Table(data+'U08BH2_error_filtered_0.1.fits',verbose=False)
#u09b = atpy.Table(data+'U09BH14_error_filtered_0.1.fits',verbose=False)
tmass= atpy.Table(data2+'fp_2mass.fp_psc26068.tbl',verbose=False)
#iras = atpy.Table(data2+'iras.iraspsc26148.tbl',verbose=False)

# I am going to brute this now and later figure out the most elegant 
# array-based solution; for now I just want working code.

delta = math.cos(math.radians(tmass.dec.mean()))

#brights = numpy.where(tmass.j_m < 10)[0]
#print brights

min_offset = -0.1* numpy.ones_like(tmass.ra)


counter = 1

for s1 in range(tmass.ra.shape[0]):
    print "================ Source number %d: ================" % counter

    
    c1 = coords.Position((tmass.ra[s1],tmass.dec[s1]))
    print "Coordinates: ", c1
    #WFCAM is still in radians!!!
    # let's take the intersection of these four where queries
    w1= where(numpy.degrees(userv.DEC) < tmass.dec[s1] + boxsize)[0]
    w2= where(numpy.degrees(userv.DEC) > tmass.dec[s1] - boxsize)[0]
    w3= where(numpy.degrees(userv.RA)  < tmass.ra[s1] + boxsize/delta)[0]
    w4= where(numpy.degrees(userv.RA)  > tmass.ra[s1] - boxsize/delta)[0]

    rabox = numpy.intersect1d(w3,w4)
    decbox= numpy.intersect1d(w1,w2)

    box = numpy.intersect1d(rabox,decbox)
    
#    print box, "these guys got matched"

    offset = -1.* numpy.ones_like(userv.RA[box])
    if offset.size != 0:
        for s2 in range(userv.RA[box].shape[0]):
            # Don't forget! WFCAM data comes in radians!
            c2 = coords.Position((userv.RA[box][s2],userv.DEC[box][s2]),units='radians')
            offset[s2] = c1.angsep(c2).arcsec()
        
    #print "The min,max,mean offset: ", str(offset.min()), str(offset.max()), str(offset.mean()) #all in arcsec

        min_offset[s1] = offset.min()

        print "Offset to best match: %f arcsec" % offset.min()
        match = box[numpy.where(offset == offset.min())][0]
        print "Index number of catalog 2 source match: %s" % str(match)
#numpy.where(offset.min())
    #print userv.data[match]
    #try also userv.row(match)

        print "2MASS jmag: "+str(tmass.j_m[s1])
        print "WFCAM jmag: "+str(userv.JAPERMAG3[match])
    else:
        print "No match found within %d arcsec" % boxsize_as
    counter += 1
    #if counter == 100: break

print "finished successfully"

plt.plot(min_offset,'r.')
plt.show()
