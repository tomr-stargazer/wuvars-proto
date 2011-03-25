import atpy
import numpy
import matplotlib.pyplot as plt

data = '/home/trice/reu/DATA/Merged_Catalogs/'

# don't forget: 0 mag is brighter than 30 mag

userv = atpy.Table(data+'USERV1678_error_filtered_0.1.fits')
u08b = atpy.Table(data+'U08BH2_error_filtered_0.1.fits')
u09b = atpy.Table(data+'U09BH14_error_filtered_0.1.fits')

l = [userv,u08b,u09b]

for t in l:
    t.j = t.JAPERMAG3
    t.h = t.HAPERMAG3
    t.k = t.KAPERMAG3
    t.jmh = t.JMHPNT
    t.hmk = t.HMKPNT
    
    print t.j.min(),t.h.min(),t.k.min()

    t.brights = numpy.where(t.j < 11)
    print t.brights
    t.RAb = numpy.degrees(t.RA[t.brights])
    t.DECb= numpy.degrees(t.DEC[t.brights])

plt.figure(1)

plt.plot(numpy.degrees(userv.RA)/15,numpy.degrees(userv.DEC),'r.',alpha=0.01)
plt.plot(numpy.degrees(u08b.RA)/15,numpy.degrees(u08b.DEC),'g.',alpha=0.01)
plt.plot(numpy.degrees(u09b.RA)/15,numpy.degrees(u09b.DEC),'b.',alpha=0.01)


plt.plot(userv.RAb/15, userv.DECb,'r^')
plt.plot(u08b.RAb/15, u08b.DECb,'gs',alpha=0.5)
plt.plot(u09b.RAb/15, u09b.DECb,'b.')

plt.ylabel("Dec, degrees")
plt.xlabel("RA, hours")
plt.title("Braid nebula brightest J-band sources, WFCAM")

plt.text(20.94,52.2,"Red = USERV\nGreen = U08B\nBlue = U09B")


plt.show()
