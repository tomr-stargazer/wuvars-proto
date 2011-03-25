#/usr/bin/python

# This is a test of the neat matplotlib and ATpy stuff

import atpy
import numpy
import matplotlib.pyplot as plt

data = '/home/trice/reu/DATA/U08BH2v20091106/'

lt = []
for i in range(1,5):
    lt.append(atpy.Table(data+'w20080919_00162_sf_st_cat.fits', hdu=i))

print 'list created'

for t in lt:
    plt.plot(numpy.degrees(t.RA)/15,numpy.degrees(t.DEC),'ko')

plt.show()

print 'plot plotted'
