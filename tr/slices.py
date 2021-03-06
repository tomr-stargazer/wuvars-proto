import atpy
import numpy

import coords

where = numpy.where
sect = numpy.intersect1d

def slices (table, n_slices, m_slices) :
    '''
    Slices a position table into an n x m list of sub-tables; returns a list of
    those tables so you can iterate over them and do whatever you want, 
    like plot K-band histograms or find the most constant source.
    '''

    
    ramin = table.RA.min()
    ramax = table.RA.max()
    demin = table.DEC.min()
    demax = table.DEC.max()

    raside = (ramax - ramin) / n_slices
    deside = (demax - demin) / m_slices

    raspace = numpy.linspace(ramin,ramax,n_slices, endpoint=False)
    despace = numpy.linspace(demin,demax,m_slices, endpoint=False)

    slicelist = []

    for left in raspace:
        
        for bottom in despace:
            w1 = where(table.RA > left)[0]
            w2 = where(table.RA < left+raside)[0]
            w3 = where(table.DEC > bottom)[0]
            w4 = where(table.DEC < bottom+deside)[0]
            
            w = sect(sect(w1,w2),sect(w3,w4))

            u = table.where(w)

            slicelist.append(u)

    return slicelist
