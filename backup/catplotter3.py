import atpy
import numpy
import matplotlib.pyplot as plt

data = '/home/trice/reu/DATA/Merged_Catalogs/'
data2= '/home/trice/reu/DATA/2MASS/'

# don't forget: 0 mag is brighter than 30 mag

userv = atpy.Table(data+'USERV1678_error_filtered_0.1.fits')
u08b = atpy.Table(data+'U08BH2_error_filtered_0.1.fits')
u09b = atpy.Table(data+'U09BH14_error_filtered_0.1.fits')
tm   = atpy.Table(data2+'fp_2mass.fp_psc26068.tbl')
iras = atpy.Table(data2+'iras.iraspsc26148.tbl')

l = [userv,u08b,u09b]

for t in l:
    t.j = t.JAPERMAG3
    t.h = t.HAPERMAG3
    t.k = t.KAPERMAG3
    t.jmh = t.JMHPNT
    t.hmk = t.HMKPNT
    
    print t.j.min(),t.h.min(),t.k.min()
'''
# FIGURE 1: Color-magnitude
plt.figure(1)
plt.plot(userv.hmk,userv.j,'ro',alpha=0.5)

plt.plot(u08b.hmk,u08b.j,'go',alpha=0.5)

plt.plot(u09b.hmk,u09b.j,'bo',alpha=0.5)


plt.ylabel("J magnitude")
plt.xlabel("H-K color")
plt.title("Braid nebula region H-K vs J, WFCAM data")
#plt.axis([-2,5,22,8])
plt.gca().axes.invert_yaxis()

plt.text(-6,20,"Red = USERV\nGreen = U08B\nBlue = U09B")

plt.show()

## FIGURE 2: Color-color
plt.figure(2)
plt.plot(userv.hmk,userv.jmh,'ro', alpha=0.5)

plt.plot(u08b.hmk,u08b.jmh,'go',alpha=0.5)

plt.plot(u09b.hmk,u09b.jmh,'bo',alpha=0.5)

plt.ylabel("J-H color")
plt.xlabel("H-K color")
plt.title("Braid nebula region H-K vs J-H, WFCAM")
#plt.axis([-7,7,-7,7])

plt.text(-6,-4,"Red = USERV\nGreen = U08B\nBlue = U09B")

plt.show()
'''
## FIGURE 3: Positions

plt.figure(3)
plt.plot(numpy.degrees(userv.RA)/15,numpy.degrees(userv.DEC),'r.',alpha=0.1)
#plt.plot(numpy.degrees(u08b.RA)/15,numpy.degrees(u08b.DEC),'g.',alpha=0.1)
plt.plot(numpy.degrees(u09b.RA)/15,numpy.degrees(u09b.DEC),'b.',alpha=0.05)
plt.plot(tm.ra/15,tm.dec,'k.',alpha=0.2)
plt.plot(iras.ra/15,iras.dec,'g*',ms=15)

plt.ylabel("Dec, degrees")
plt.xlabel("RA, hours")
plt.title("Braid nebula region coverage, WFCAM+2MASS+IRAS")

plt.text(20.93,52.2,"Red = USERV\nBlack = 2MASS\nBlue = U09B\nGreen = IRAS")

plt.show()
