import atpy
import numpy
import matplotlib.pyplot as plt

data = '/home/trice/reu/DATA/'

# don't forget: 0 mag is brighter than 30 mag

t = atpy.Table(data+'results15_2_42_24_15868.fits')

j = t.JAPERMAG3
h = t.HAPERMAG3
k = t.KAPERMAG3
jmh = t.JMHPNT
hmk = t.HMKPNT

mag = [j,h,k]

#for i in mag:
#    i[i<0] = -100
j[j<0] = -100
#h[h<0] = -200
#k[k<0] = -300
hmk[hmk<-50] = -50
jmh[hmk<-70] = -70

'''
plt.plot(hmk,j,'bs')
plt.ylabel("J magnitude")
plt.xlabel("H-K color")
plt.title("Braid nebula region H-K vs J, U08B")
plt.axis([-7,7,35,5])
'''

plt.plot(hmk,jmh,'g.')
plt.ylabel("J-H color")
plt.xlabel("H-K color")
plt.title("Braid nebula region H-K vs J-H, U08B")
plt.axis([-7,7,-7,7])

plt.show()
