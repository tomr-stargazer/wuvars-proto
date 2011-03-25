import atpy,numpy

data = '/home/trice/reu/DATA/'

t = atpy.Table(data+'results15_21_51_0_15909.fits')

print t.columns

j = t.JAPERMAG3
h = t.HAPERMAG3
k = t.KAPERMAG3
jmh = t.JMHPNT
hmk = t.HMKPNT

print j.max(),h.max(),k.max()

a =  numpy.where(j==j.max())
print t.data[a]
print h[a],k[a]
