import atpy
import numpy
import matplotlib.pyplot as plt

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
#is this how i want to do it? maybe a dictionary would be better,
# or something conducive to storing in 
    
    w = numpy.where( table.SOURCEID == sid )

#    date = table.MEANMJDOBS[w] - table.MEANMJDOBS.min()

    jcol = table.JAPERMAG3[w]
    hcol = table.HAPERMAG3[w]
    kcol = table.KAPERMAG3[w]

# Let's refactor this code so it uses tuple assignments and is less ugly.
# And is there any reason why I have all these assignments and not just return 
# the results of the functions? Save so much space.

    num = w[0].size #jcol.size #(try w.size?)
    j_min = jcol.min()
    j_max = jcol.max()
    j_mean = jcol.mean()
    j_rms = jcol.std()
    h_min = hcol.min()
    h_max = hcol.max()
    h_mean = hcol.mean()
    h_rms = hcol.std()
    k_min = kcol.min()
    k_max = kcol.max()
    k_mean = kcol.mean()
    k_rms = kcol.std()

    return (num, 
            j_min, j_max, j_mean, j_rms, 
            h_min, h_max, h_mean, h_rms, 
            k_min, k_max, k_mean, k_rms)

#Here comes my table-making function that probably will have to run overnight

def stat_write (table, outfile) :
    ''' Writes a table with statistical information about every source.
    '''

    sidarr = make_sidset(table)

    #let's create a bunch of empty/zeroed columns, do a for loop to 
    # use arraystat to fill the columns row by row, and then add them columns
    # to the output table which we can then write.

    Output = atpy.Table()
    



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

def plot_lc (table,sid, text=True) :
    ''' Plots J,H,K lightcurves WITH ERRORBARS for a given input source.

    Written with WFCAM columns in mind, specifically like from WSERV1.
    '''
    
    w = numpy.where( table.SOURCEID == sid )

    date = table.MEANMJDOBS[w] - table.MEANMJDOBS.min()

    jcol = table.JAPERMAG3[w]
    hcol = table.HAPERMAG3[w]
    kcol = table.KAPERMAG3[w]

    jerr = table.JAPERMAG3ERR[w]
    herr = table.HAPERMAG3ERR[w]
    kerr = table.KAPERMAG3ERR[w]

    plt.errorbar(date,jcol,yerr=jerr,fmt='b-o',ecolor='k')
    plt.errorbar(date,hcol,yerr=herr,fmt='g-o',ecolor='k')
    plt.errorbar(date,kcol,yerr=kerr,fmt='r-o',ecolor='k')

    plt.gca().invert_yaxis()

    if text:
        plt.text(20,jcol.max()+.1,"J-band RMS: %f" % jcol.std())
        plt.text(20,hcol.max()+.1,"H-band RMS: %f" % hcol.std())
        plt.text(20,kcol.max()+.1,"K-band RMS: %f" % kcol.std())

    plt.ylabel("WFCAM magnitude")
    plt.xlabel("Julian days since 04/26/2008")
    plt.title("J, H, K with errorbars. Source ID %d." % sid)

    plt.show()
    return

plot_lc_err = plot_lc
    
