''' A collection of statistical methods for the WSERV data from WFCAM. '''

import atpy
import numpy
import matplotlib.pyplot as plt
import coords

tmass = atpy.Table('/home/trice/reu/DATA/2MASS/fp_2mass.fp_psc26068.tbl', \
                       verbose=False)

def make_sidset ( table ) :
    ''' Returns an array of unique source IDs'''
    # sidcol = table.SOURCEID
    # sidset = set(sidcol)
    # sidlist= list(sidset)
    # sidarr = numpy.array(sidlist)

    # I think it's hilarious I can get the desired effect via a triple-cast.
    # Also: this is exactly the kind of thing that IDL can never ever do.
    # Boom! One-liner.
    return numpy.array(list(set(table.SOURCEID)))


''' 
Let's take one source in a table (known by its Source ID) and spit out some
statistics on it as return values. 
Returns a tuple of min, max, mean, stddev
'''

def arraystat (table, sid) :
    ''' Inputs: source ID integer, atpy.Table. Outputs: tuple of 13 numbers'''
    
    w = numpy.where( table.SOURCEID == sid )

    jcol = table.JAPERMAG3[w]
    hcol = table.HAPERMAG3[w]
    kcol = table.KAPERMAG3[w]
    jmhcol=table.JMHPNT[w]
    hmkcol=table.HMKPNT[w]
    racol= table.RA[w]
    decol= table.DEC[w]

    #print "running stats on source %s" % sid

    return (w[0].size,
            jcol.min(),jcol.max(),jcol.mean(),jcol.std(),
            hcol.min(),hcol.max(),hcol.mean(),hcol.std(),
            kcol.min(),kcol.max(),kcol.mean(),kcol.std(),
            jmhcol.min(),jmhcol.max(),jmhcol.mean(),jmhcol.std(),
            hmkcol.min(),hmkcol.max(),hmkcol.mean(),hmkcol.std(),
            racol.mean(), decol.mean()
            )

#Here comes my table-making function that probably will have to run overnight

def stat_write (table, outfile, Test = False) :
    ''' Writes a table with statistical information about every source.

    This function is VERY SLOW for large numbers of sources (~100k).
    It took my computer five hours.

    Keyword arguments:
    Test: if true, writes only 100 lines and then leaves the rest empty.
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
            ra[i], dec[i] \
         = arraystat(table, sid)
        
        if (i[0] > 10 and Test): 
            print i
            break

    # The following could probably be achieved much more easily 
    # with a dictionary and a for loop.

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


# def plot_lc (table,sid) :
#     ''' Plots J,H,K lightcurves for a given input source.

#     Written with WFCAM columns in mind, specifically like from WSERV1.
#     '''
    
#     w = numpy.where( table.SOURCEID == sid )

#     date = table.MEANMJDOBS[w] - table.MEANMJDOBS.min()

#     jcol = table.JAPERMAG3[w]
#     hcol = table.HAPERMAG3[w]
#     kcol = table.KAPERMAG3[w]

#     plt.plot(date,jcol,'b-o')
#     plt.plot(date,hcol,'g-o')
#     plt.plot(date,kcol,'r-o')

#     plt.gca().invert_yaxis()

#     plt.show()
#     return

# def small_match (table, sid) :
#     ''' Returns a 2MASS identifier so you can plot 2MASS mags for a source.
    
#     Assumes WFCAM-style inputs, specifically RA/Dec in radians.
#     Tries to match within 1 arcsec
#     '''

#     # Let's grab the position from the first night, that ought to be good nuff
#     ra = table.RA [numpy.where(table.SOURCEID == sid)][0] 
#     dec= table.DEC[numpy.where(table.SOURCEID == sid)][0] 

#     delta = math.cos(abs(dec))
#     boxsize = 1 / 3600.

#     p1 = coords.Position( (ra,dec), units = 'rad')
    
    
#     w1 = where(tmass.ra < ra + boxsize/delta)[0]
#     w2 = where(tmass.ra > ra - boxsize/delta)[0]
#     w3 = where(tmass.dec < dec + boxsize)[0]
#     w4 = where(tmass.dec > dec - boxsize)[0]
# '''
#         # Let's slice a box around our source
#     box = sect(sect(w1,w2),sect(w3,w4))

#         # And calculate offsets to all sources inside that box
#     offset = -1. * numpy.ones_like(tmass.ra[box])
#     if offset.size != 0:
#         for s2 in range(offset.size):
#             p2 = coords.Position( (tmass.ra[box][s2], tmass.dec[box][s2]),  units = 'deg')
#             offset[s2] = p1.angsep(p2).arcsec()

#             # If the closest match is within our matching circle
#         if offset.min() < 1.:
#             min_offset[s1] = offset.min()
#             match[s1] = box[ where(offset == offset.min()) ][0]
#             vprint( "Source %d: Matched with %f arcsec" \
#                         % (counter, offset.min() ) )
#         else:
#             vprint( "Source %d: Failed to match" % counter)
#     else:
#         vprint( "Source %d: Failed to match" % counter)
# '''
# the above is totally incomplete



def plot_lc (table,sid, text=True, sup='', outfile='') :
    ''' Plots J,H,K lightcurves WITH ERRORBARS for a given input source.

    Written with WFCAM columns in mind, specifically like from WSERV1.
    '''
    
    w = numpy.where( table.SOURCEID == sid )

    ra1, ra2 = 314.36, 315.77
    dec1,dec2= 52.02, 52.92

    rabox = [ra1, ra1, ra2, ra2, ra1]
    decbox= [dec1,dec2,dec2,dec1,dec1]

    sra, sdec = table.RA[w][0], table.DEC[w][0]

    date = table.MEANMJDOBS[w] - 54579

    jcol = table.JAPERMAG3[w]
    hcol = table.HAPERMAG3[w]
    kcol = table.KAPERMAG3[w]

    jerr = table.JAPERMAG3ERR[w]
    herr = table.HAPERMAG3ERR[w]
    kerr = table.KAPERMAG3ERR[w]

    plt.clf()

    # This is the size of the lightcurve box
    plt.axes([.05,.1,.7,.8])

    plt.errorbar(date,jcol,yerr=jerr,fmt='b-o',ecolor='k')
    plt.errorbar(date,hcol,yerr=herr,fmt='g-o',ecolor='k')
    plt.errorbar(date,kcol,yerr=kerr,fmt='r-o',ecolor='k')

    plt.gca().invert_yaxis()

    if text:
        plt.text(20,jcol.max()+.1,"J-band RMS: %f" % jcol.std())
        plt.text(20,hcol.max()+.1,"H-band RMS: %f" % hcol.std())
        plt.text(20,kcol.max()+.1,"K-band RMS: %f" % kcol.std())

    plt.ylabel("WFCAM magnitude")
    plt.xlabel("Julian days since 04/23/2008")
    plt.title("J, H, K with errorbars. Source ID %d." % sid)
    plt.suptitle(sup)

    # This is the size of the position box
    plt.axes([.775,.55,.2,.35])
    plt.plot(rabox, decbox)
    plt.plot(numpy.degrees(sra),numpy.degrees(sdec),'rD')
    
    plt.gca().invert_xaxis()

    dx = 1/8. * (ra2-ra1)
    dy = 1/8. * (dec2-dec1)
    arrx = ra1 + 2*dx
    arry = dec1+6*dy
    plt.arrow(arrx,arry,dx,0)
    plt.arrow(arrx,arry,0,dy)

    plt.ylabel("Dec, degrees")
    plt.xlabel("RA, degrees")

    plt.axis('off')

    if outfile == '':
        plt.show()
    else:
        plt.savefig(outfile)
    return

plot_lc_err = plot_lc
    

# this don't work
"""
def plot_lc_list (table, sidlist ):
    ''' Plots, sequentially, a bunch of lightcurves interactively.'''

    plt.clf()
    for sid in sidlist:
        plot_lc (table, sid)
        plt.show()
"""
