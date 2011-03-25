import atpy 
import numpy
import matplotlib.pyplot as plt
import coords

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

# Next time I'll say "for s1 in range(tmass.ra.shape[0])" or something

brights = numpy.where(tmass.j_m < 10)[0]
print brights

counter = 1

for s1 in brights:
    print "================ Source number %d: ================" % counter

    offset = numpy.zeros_like(userv.RA)
    for s2 in range(userv.RA.shape[0]):
        c1 = coords.Position((tmass.ra[s1],tmass.dec[s1]))
    # Don't forget! WFCAM data comes in radians!
        c2 = coords.Position((userv.RA[s2],userv.DEC[s2]),units='radians')
        offset[s2] = c1.angsep(c2).arcsec()
        
    #print "The min,max,mean offset: ", str(offset.min()), str(offset.max()), str(offset.mean()) #all in arcsec
    print "Offset to best match: %f arcsec" % offset.min()
    match = numpy.where(offset == offset.min())
    print "Index number of catalog 2 source match: %s" % str(match[0][0])
#numpy.where(offset.min())
    #print userv.data[match]

    print "2MASS jmag: "+str(tmass.j_m[s1])
    print "WFCAM jmag: "+str(userv.JAPERMAG3[match[0]][0])
    counter += 1
    if counter == 10: break

print "finished successfully"
