import atpy 
import numpy
import matplotlib.pyplot as plt
import coords

''' Here's a test of my matching algorithm on some dummy data'''

tra = numpy.array([ 10.] )
tdec= numpy.array([ 20.] )

ura = numpy.array([10.,10.01,10,10.02])
udec= numpy.array([20.01,20.,20.01,19.98])

s1 = 0

offset = numpy.zeros_like(ura)
for s2 in range(ura.shape[0]):
    c1 = coords.Position((tra[s1],tdec[s1]))
    # Don't forget! WFCAM data comes in radians!
    c2 = coords.Position((ura[s2],udec[s2]))#,units='radians')
    offset[s2] = c1.angsep(c2).arcsec()

print "The min,max,mean offset: ", str(offset.min()), str(offset.max()), str(offset.mean()) #all in arcsec
match = numpy.where(offset == offset.min())
print match
#numpy.where(offset.min())
#print userv.data[match]

#print "2MASS jmag: "+str(tmass.j_m[s1])
#print "WFCAM jmag: "+str(userv.JAPERMAG3[match])
