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

plt.figure(1)

#plt.subplot(131)

plt.plot(userv.hmk,userv.j,'ro',alpha=0.8)
plt.ylabel("J magnitude")
plt.xlabel("H-K color")
plt.title("Braid nebula region H-K vs J, USERV")
plt.axis([-2,5,22,8])

#plt.show()

#plt.subplot(132)

plt.plot(u08b.hmk,u08b.j,'go',alpha=0.8)
#plt.ylabel("J magnitude")
#plt.xlabel("H-K color")
#plt.title("Braid nebula region H-K vs J, U08B")
#plt.axis([-2,5,22,8])

plt.plot(u09b.hmk,u09b.j,'bo',alpha=0.8)

plt.show()
'''
plt.figure(2)
plt.plot(userv.hmk,userv.jmh,'g.')
plt.ylabel("J-H color")
plt.xlabel("H-K color")
plt.title("Braid nebula region H-K vs J-H, USERV")
#plt.axis([-7,7,-7,7])



plt.show()
'''
