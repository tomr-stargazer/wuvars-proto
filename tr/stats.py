''' A collection of statistical methods for the WSERV data from WFCAM. 

Publicly available:

make_sidset (for now)

'''

import atpy
import numpy
import matplotlib.pyplot as plt
import coords
import stetson

tmass = atpy.Table('/home/trice/reu/DATA/2MASS/fp_2mass.fp_psc26068.tbl', 
                   verbose=False)

def make_sidset ( table ) :
    ''' Returns an array of unique source IDs.'''

    return numpy.array(list(set(table.SOURCEID)))

# def new_arraystat (table, sid) :
#     ''' Inputs: source ID integer, atpy.Table. Outputs: Dictionary!

#     Calculates statistics on a given source ID from values in a table.
#     '''
    
#     w = numpy.where( table.SOURCEID == sid )

#     jcol = table.JAPERMAG3[w]
#     hcol = table.HAPERMAG3[w]
#     kcol = table.KAPERMAG3[w]
#     jmhcol=table.JMHPNT[w]
#     hmkcol=table.HMKPNT[w]
#     racol= table.RA[w]
#     decol= table.DEC[w]

#     #print "running stats on source %s" % sid


#     return (w[0].size,
#             jcol.min(),jcol.max(),jcol.mean(),jcol.std(),
#             hcol.min(),hcol.max(),hcol.mean(),hcol.std(),
#             kcol.min(),kcol.max(),kcol.mean(),kcol.std(),
#             jmhcol.min(),jmhcol.max(),jmhcol.mean(),jmhcol.std(),
#             hmkcol.min(),hmkcol.max(),hmkcol.mean(),hmkcol.std(),
#             racol.mean(), decol.mean()
#             )


# #old tuple-based code
def arraystat (table, sid) :
    ''' Inputs: source ID integer, atpy.Table. Outputs: tuple of 26 numbers.

    Calculates statistics on a given source ID from values in a table.
    '''
    
    w = numpy.where( table.SOURCEID == sid )

    jcol = table.JAPERMAG3[w]
    hcol = table.HAPERMAG3[w]
    kcol = table.KAPERMAG3[w]
    jmhcol=table.JMHPNT[w]
    hmkcol=table.HMKPNT[w]
    racol= table.RA[w]
    decol= table.DEC[w]
    jppcol=table.JPPERRBITS[w]
    hppcol=table.HPPERRBITS[w]
    kppcol=table.KPPERRBITS[w]

    #print "running stats on source %s" % sid

    return (w[0].size,
            jcol.min(),jcol.max(),jcol.mean(),jcol.std(),
            hcol.min(),hcol.max(),hcol.mean(),hcol.std(),
            kcol.min(),kcol.max(),kcol.mean(),kcol.std(),
            jmhcol.min(),jmhcol.max(),jmhcol.mean(),jmhcol.std(),
            hmkcol.min(),hmkcol.max(),hmkcol.mean(),hmkcol.std(),
            jppcol.max(),hppcol.max(),kppcol.max(),
            racol.mean(), decol.mean()
            )

#Here comes my table-making function that probably will have to run overnight

def stat_write (table, outfile, Test = False) :
    ''' Writes a table with statistical information about every source.

    This function is VERY SLOW for large numbers of sources (~100k).
    It took my computer five hours.

    Keyword arguments:
    Test: if true, writes only 100 lines and then leaves the rest empty.

    This function is strongly dependent on the behavior of arraystat.
    '''

    sidarr = make_sidset(table)

    #let's create a bunch of empty/zeroed columns, do a for loop to 
    # use arraystat to fill the columns row by row, and then add them columns
    # to the output table which we can then write.

    Output = atpy.Table()

    # wonder if there's a better way to do this...
    num = numpy.ones_like(sidarr)
    j_min = numpy.ones_like(sidarr)*1.0
    j_max = numpy.ones_like(sidarr)*1.0
    j_mean = numpy.ones_like(sidarr)*1.0
    j_std = numpy.ones_like(sidarr)*1.0
    h_min = numpy.ones_like(sidarr)*1.0
    h_max = numpy.ones_like(sidarr)*1.0
    h_mean = numpy.ones_like(sidarr)*1.0
    h_std = numpy.ones_like(sidarr)*1.0
    k_min = numpy.ones_like(sidarr)*1.0
    k_max = numpy.ones_like(sidarr)*1.0
    k_mean = numpy.ones_like(sidarr)*1.0
    k_std = numpy.ones_like(sidarr)*1.0
    jmh_min = numpy.ones_like(sidarr)*1.0
    jmh_max = numpy.ones_like(sidarr)*1.0
    jmh_mean = numpy.ones_like(sidarr)*1.0
    jmh_std = numpy.ones_like(sidarr)*1.0
    hmk_min = numpy.ones_like(sidarr)*1.0
    hmk_max = numpy.ones_like(sidarr)*1.0
    hmk_mean = numpy.ones_like(sidarr)*1.0
    hmk_std = numpy.ones_like(sidarr)*1.0
    jpp_max = numpy.ones_like(sidarr)
    hpp_max = numpy.ones_like(sidarr)
    kpp_max = numpy.ones_like(sidarr)
    ra  = numpy.ones_like(sidarr)*1.0
    dec = numpy.ones_like(sidarr)*1.0


    # The following forloop is VERY slow but I can't think of a way to make it
    # array-based.
    for sid in sidarr:
        i = numpy.where(sidarr == sid)[0]
        

        num[i], \
            j_min[i], j_max[i], j_mean[i], j_std[i], \
            h_min[i], h_max[i], h_mean[i], h_std[i], \
            k_min[i], k_max[i], k_mean[i], k_std[i],  \
            jmh_min[i], jmh_max[i], jmh_mean[i], jmh_std[i],  \
            hmk_min[i], hmk_max[i], hmk_mean[i], hmk_std[i],  \
            jpp_max[i], hpp_max[i], kpp_max[i], \
            ra[i], dec[i] \
         = arraystat(table, sid)
        
        if (i[0] > 100 and Test): 
            print i
            print "End of test"
            break

    # The following could probably be achieved much more easily 
    # with a dictionary and a for loop.
         # (except I just realized a dictionary is unsorted)


    # And with .add_empty_column ...
    Output.add_column('SOURCEID',sidarr)
    Output.add_column('RA', ra, unit='RADIANS')
    Output.add_column('DEC', dec, unit='RADIANS')
    Output.add_column('n_detect', num)
    Output.add_column('j_min', j_min)
    Output.add_column('j_max', j_max)
    Output.add_column('j_mean', j_mean)
    Output.add_column('j_rms', j_std)
    Output.add_column('h_min', h_min)
    Output.add_column('h_max', h_max)
    Output.add_column('h_mean', h_mean)
    Output.add_column('h_rms', h_std)
    Output.add_column('k_min', k_min)
    Output.add_column('k_max', k_max)
    Output.add_column('k_mean', k_mean)
    Output.add_column('k_rms', k_std)
    Output.add_column('jpp_max', jpp_max)
    Output.add_column('hpp_max', hpp_max)
    Output.add_column('kpp_max', kpp_max)

    Output.write(outfile, overwrite=True)

    print "Wrote output to %s" % outfile

def stat_write_test (table) :
    ''' A method to test the output of stat_write. 

    Writes output to /home/trice/reu/DATA/Merged_Catalogs/stat/test.fits .'''

    stat_write(table,'/home/trice/reu/DATA/Merged_Catalogs/stat/test.fits', 
               Test = True)

def stat_write_smart (table, outfile) :
    ''' Calls stat_write with the output into a specific directory. 

    That directory is /home/trice/reu/DATA/Merged_Catalogs/stat/ .'''

    stat_write(table,'/home/trice/reu/DATA/Merged_Catalogs/stat/' + outfile)

def sidI(table, sid):
    ''' Calculates the Welch/Stetson variability index I for a given source
    in a table. Uses J and K mags'''
    ws = numpy.where(table.SOURCEID == sid)
    
    if ws[0].size == 0:
        print "Source not found"
        return
    
    b = table.JAPERMAG3[ws]
    bs= table.JAPERMAG3ERR[ws]
    v = table.KAPERMAG3[ws]
    vs= table.KAPERMAG3ERR[ws]

    return stetson.I(b,bs,v,vs)
