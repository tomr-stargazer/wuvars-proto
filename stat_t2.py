import atpy 
import numpy as np
import matplotlib.pyplot as plt

where = np.where

data = '/home/trice/reu/'
data2= '/home/trice/reu/DATA/2MASS/'
data3 = '/home/trice/reu/DATA/Merged_Catalogs/'
data4= '/home/trice/reu/DATA/'

Tmass = atpy.Table(data2+'fp_2mass.fp_psc26068.tbl',verbose=False)
Userv = atpy.Table(data3+'USERV1678_error_filtered_0.1.fits',verbose=False)
Match = atpy.Table(data + 'test_output.fits',verbose=False)

# Let's try to do this on the whole datachunk before slicing it up

col_1 = Tmass.j_m[where(Match.match > 0)]
col_2 = Userv.JAPERMAG3[Match.match[where(Match.match > 0)]]

#plt.plot(col_1,col_2,'ro')
#plt.show()

#exit()

print "Here are their sizes: %d, %d" % (col_1.shape[0],col_2.shape[0])

#calculating 

jmin, jmax = 18, 10
jrange = np.arange(jmin,jmax,-2) #let's start with big interval size

#pieces = jrange.shape[0]

means = np.zeros_like(jrange)
medians=np.zeros_like(jrange)
stddevs=np.zeros_like(jrange)
lower = np.zeros_like(jrange)
upper = np.zeros_like(jrange)

for i in range(1,jrange.shape[0]):

    #print jrange[i],jrange[i-1]
    lower[i],upper[i] = jrange[i-1],jrange[i]

    # make a slice
    interval = col_2[ where (col_1 < jrange[i]) and \
                          where (col_1 > jrange[i-1])]
    #print "The interval is this big: ", interval.shape[0]
    print "The interval looks like this: ", interval

    means[i] = interval.mean()
    medians[i] = np.median(interval)
    stddevs[i] = interval.std()


Output = atpy.Table()
Output.add_column('2mass_lower',lower)
Output.add_column('2mass_upper',upper)
Output.add_column('means',means)
Output.add_column('medians',medians)
Output.add_column('stds',stddevs)

Output.write(data4+'test_stat.fits',overwrite=True)

print Output.data
