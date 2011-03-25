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

col_1 = Tmass.j_m[where(Match.match != 0)]
col_2 = Userv.JAPERMAG3[Match.match[where(Match.match != 0)]]

#calculating 

jmin, jmax = 18, 10
jrange = np.arange(jmin,jmax,0.2)

pieces = 1

means = np.zeros(pieces)
medians=np.zeros(pieces)
stddevs=np.zeros(pieces)

for dude in range(pieces):

    # make a slice
    
    interval = col_2 #next this will be an interval
    
    means[dude] = interval.mean()
    medians[dude] = np.median(interval)
    stddevs[dude] = interval.std()
    

Output = atpy.Table()
Output.add_column('means',means)
Output.add_column('medians',medians)
Output.add_column('stds',stddevs)

Output.write(data4+'test_stat.fits',overwrite=True)

print Output.data
