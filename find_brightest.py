import atpy,numpy

data = '/home/trice/reu/DATA/Merged_Catalogs/'

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


a =  numpy.where(userv.j==userv.j.min())
print userv.data[a]
print userv.h[a],userv.k[a]

print "The brightest source in userv is located at " + str(numpy.degrees(userv.RA[a])/15) +', ' + str(numpy.degrees(userv.DEC[a]))

b =  numpy.where(u08b.j==u08b.j.min())
print "The brightest source in u08b is located at " + str(numpy.degrees(u08b.RA[b])/15) +', ' + str(numpy.degrees(u08b.DEC[b]))

c =  numpy.where(u09b.j==u09b.j.min())
print "The brightest source in u09b is located at " + str(numpy.degrees(u09b.RA[c])/15) +', ' + str(numpy.degrees(u09b.DEC[c]))
