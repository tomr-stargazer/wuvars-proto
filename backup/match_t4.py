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

userv= atpy.Table(data+'USERV1678_error_filtered_0.1.fits')
u08b = atpy.Table(data+'U08BH2_error_filtered_0.1.fits')
u09b = atpy.Table(data+'U09BH14_error_filtered_0.1.fits')
tmass= atpy.Table(data2+'fp_2mass.fp_psc26068.tbl')
iras = atpy.Table(data2+'iras.iraspsc26148.tbl')

# I am going to brute this now and later figure out the most elegant 
# array-based solution; for now I just want working code.

delta = math.cos(math.radians(tmass.dec.mean()))

brights = numpy.where(tmass.j_m < 10)[0]
print brights

counter = 1

for s1 in brights:
    print "================ Source number %d: ================" % counter

    
    c1 = coords.Position((tmass.ra[s1],tmass.dec[s1]))
    #WFCAM is still in radians!!!
    w1= where(numpy.degrees(userv.DEC) < tmass.dec[s1] + boxsize)
    w2= where(numpy.degrees(userv.DEC[w1]) > tmass.dec[s1] - boxsize)
    w3= where(numpy.degrees(userv.RA[w2])  < tmass.ra[s1] + boxsize/delta)
    w4= where(numpy.degrees(userv.RA[w3])  > tmass.ra[s1] - boxsize/delta)

    box = w4
    
#    print box

    offset = -1.* numpy.ones_like(userv.RA[box])
    if offset.size != 0:
        for s2 in range(userv.RA[box].shape[0]):
            # Don't forget! WFCAM data comes in radians!
            c2 = coords.Position((userv.RA[box][s2],userv.DEC[box][s2]),units='radians')
            offset[s2] = c1.angsep(c2).arcsec()
        
    #print "The min,max,mean offset: ", str(offset.min()), str(offset.max()), str(offset.mean()) #all in arcsec

        print "Offset to best match: %f arcsec" % offset.min()
        match = numpy.where(offset == offset.min())
        print "Index number of catalog 2 source match: %s" % str(match[0][0])
#numpy.where(offset.min())
    #print userv.data[match]

        print "2MASS jmag: "+str(tmass.j_m[s1])
        print "WFCAM jmag: "+str(userv.JAPERMAG3[box][match[0]][0])
    else:
        print "No match found within %d arcsec" % boxsize
    counter += 1
    if counter == 10: break

print "finished successfully"
