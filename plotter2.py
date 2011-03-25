import atpy
import numpy
import matplotlib.pyplot as plt

data = '/home/trice/reu/DATA/'

# don't forget: 0 mag is brighter than 30 mag

t = atpy.Table(data+'results14_22_20_4_15864.fits')

#print t.columns

mag = {}

mag['j'] = t.JAPERMAG3
mag['h'] = t.HAPERMAG3
mag['k'] = t.KAPERMAG3

# All of the undefined magnitudes have been set to -9.9999949e+08 but I prefer
# them to be 0 for these purposes.

for i in mag:
    #print "Min is "+str(mag[i].min())
    mag[i][mag[i]<0] = 0
    #print "Now it's "+str(mag[i].min())

plt.plot(mag['j'],mag['k'],'ro')
plt.xlabel("J magnitude")
plt.ylabel("K magnitude")
plt.title("Braid nebula region J vs K, U08B")
plt.axis([35,5,28,8])

plt.show()
