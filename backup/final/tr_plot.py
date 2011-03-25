''' Functions to plot lightcurves for WFCAM data.'''

import atpy
import numpy
import matplotlib.pyplot as plt
import coords

def plot_lc (table,sid, outfile='', sup='', box=True, bands='jhk', text=True) :
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
    if box:
        plt.axes([.05,.1,.7,.8])
        
    if 'j' in bands:
        plt.errorbar(date,jcol,yerr=jerr,fmt='b-o',ecolor='k')
    if 'h' in bands:
        plt.errorbar(date,hcol,yerr=herr,fmt='g-o',ecolor='k')
    if 'k' in bands:
        plt.errorbar(date,kcol,yerr=kerr,fmt='r-o',ecolor='k')

    plt.gca().invert_yaxis()

    if text:
        if 'j' in bands:
            plt.text(20,jcol.max()+.1,"J-band RMS: %f" % jcol.std())
        if 'h' in bands:
            plt.text(20,hcol.max()+.1,"H-band RMS: %f" % hcol.std())
        if 'k' in bands:
            plt.text(20,kcol.max()+.1,"K-band RMS: %f" % kcol.std())

    plt.ylabel("WFCAM magnitude")
    plt.xlabel("Julian days since 04/23/2008")
    plt.title("J, H, K with errorbars. Source ID %d." % sid)
    plt.suptitle(sup)

    # This is the size of the position box
    if box:
        plt.axes([.775,.55,.2,.35])
        plt.plot(rabox, decbox)
        plt.plot(numpy.degrees(sra),numpy.degrees(sdec),'rD')
    
        dx = 1/8. * (ra2-ra1)
        dy = 1/8. * (dec2-dec1)
        arrx = ra1 + 2*dx
        arry = dec1+6*dy
        plt.arrow(arrx,arry,dx,0)
        plt.arrow(arrx,arry,0,dy)
        
        plt.gca().invert_xaxis()
        
        plt.ylabel("Dec, degrees")
        plt.xlabel("RA, degrees")

        plt.axis('off')

    if outfile == '':
        plt.show()
    else:
        plt.savefig(outfile)
    return


def plot_phase (table, sid, period, band='j', clear=True):
    ''' Plots magnitude as a function of phase for one source in one band. '''

    # Yo! This needs one more argument! A constant term to shift over by...

    if band.lower() not in ['j','h','k','jmh','hmk']:
        print "Error: keyword 'band' must be 'j','h', 'k', 'jmh', or 'hmk'."
        print "Keyword 'band' defaulting to 'j'."
        band = 'j'
    
    if 'm' in band.lower():
        bandname = band.upper() + "PNT"
        thing = 'color'
    else:
        bandname = band.upper() + "APERMAG3"
        thing = 'magnitude'

    w = numpy.where( table.SOURCEID == sid )

    # ra1, ra2 = 314.36, 315.77
    # dec1,dec2= 52.02, 52.92

    # rabox = [ra1, ra1, ra2, ra2, ra1]
    # decbox= [dec1,dec2,dec2,dec1,dec1]

    # sra, sdec = table.RA[w][0], table.DEC[w][0]

    date = table.MEANMJDOBS[w] #- 54579
    phase = (date % period) / period #Uh, I should probably divide this by the period

    mag = table.data[bandname][w]
    err = table.data[bandname+"ERR"][w]

    if clear:
        plt.clf()

    plt.errorbar(phase,mag,yerr=err,fmt='ko',ecolor='k')

    if len(band) == 1:
        plt.gca().invert_yaxis()

    plt.ylabel("WFCAM %s %s" % (band.upper(), thing))
    plt.xlabel("Phase")

    plt.title ("Phase-shifted lightcurve. Source ID %d. Period: %f days." %
               (sid,period))

    plt.show()
    
