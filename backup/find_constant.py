import atpy
import numpy
import matplotlib.pyplot as plt
import coords

from stat_fun3 import *
import stat_fun3

where = numpy.where
sect = numpy.intersect1d

def get_list (data_table, stat_table, slices) :
    ''' Finds constant H-band sources and returns a list of source IDs. 

    It slices the field into (slices**2) regions and picks the most constant
    source from each of them.

    This procedure is agnostic to position units.
    '''

    plots = '/home/trice/reu/PLOTS/constants/'

    ramin = data_table.RA.min()
    ramax = data_table.RA.max()
    demin = data_table.DEC.min()
    demax = data_table.DEC.max()

    raside = (ramax - ramin) / slices
    deside = (demax - demin) / slices

    raspace = numpy.linspace(ramin,ramax,slices, endpoint=False)
    despace = numpy.linspace(demin,demax,slices, endpoint=False)

    sidlist = []

    ws1 = where(stat_table.n_detect > 100)[0]
    ws2 = where(stat_table.n_detect < 150)[0]
    
    ws = sect(ws1,ws2)

    good_sources = stat_table.SOURCEID[ws]

    print "I think there are %d good sources total" % good_sources.size

    counter = 1

    for left in raspace:
        
        for bottom in despace:
            w1 = where(data_table.RA > left)[0]
            w2 = where(data_table.RA < left+raside)[0]
            w3 = where(data_table.DEC > bottom)[0]
            w4 = where(data_table.DEC < bottom+deside)[0]
            
            w = sect(sect(w1,w2),sect(w3,w4))

            # a set of source IDs of the proper sources
            right_place = numpy.array(list(set(data_table.SOURCEID[w])))

            #print "I think there are %d sources in region %d" %  \
            #    (right_place.size, counter)

            # The sids of all the good stars in this neck of the woods
            doubly_good_sources = sect(good_sources,right_place)

            print "When I intersect goods with rights I find %d sources" % \
                doubly_good_sources.size
            print doubly_good_sources

            # This line looks fishy. I want to turn my list of SIDs to a list of indices...

            wdg = numpy.empty_like ( doubly_good_sources )
            i = 0
            # I die a little bit whenever I use a forloop.
            for sid in doubly_good_sources:
                wdg[i] = where ( stat_table.SOURCEID == sid )[0]
                i += 1
            print wdg

            plt.clf()
            # Let's plot a histogram of the RMS of this region
            plt.hist(stat_table.j_rms[wdg],bins=50)
            plt.savefig(plots+'region_%d_hist' % counter)

            local_constant = stat_table.SOURCEID[where (stat_table.j_rms == stat_table.j_rms[wdg].min())]
            sidlist.append(local_constant[0])

            # Also let's plot and save the most-constant lightcurve
            stat_fun3.plot_lc(data_table,local_constant[0], sup= \
                        "This is REALLY the most constant source in region %d"\
                        % counter, outfile = plots+"region_%d" % counter)

            print "Region %d: seems like there are %d 'local constants'" % (counter, local_constant.size)
            print "Local constant is SID %d" % local_constant[0]

            counter +=1
            
            if counter == 2: return

    return sidlist
