"""
This is the package with super-functions to generate 
hundreds of figures at a time.
It calls upon most of the other wuvars packages.

Useful functions:
 # Primary functions
 
 # Helper functions
"""

import atpy
import numpy as np
import matplotlib.pyplot as plt
import spreadsheet
import plot3 as tplot

import os, errno

def mkdir_p(path):
    """ A helper function, copied from the following URL:
    http://stackoverflow.com/questions/600268/mkdir-p-functionality-in-python
    """
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST:
            pass
        else: raise


def do_it_all( table, sid_list, name_list, path='', 
               option=['lc','tables','phase'] ):
    """ 
    Does some stuff. Not sure exactly what yet, but I'll 
    want it to make tons of plots and tables for a list of 
    SOURCEIDs and their corresponding data.
    
    Parameters
    ----------
    table : atpy.Table 
        The WFCAM time-series data to be extracted
    sid_list : (list or array) of int
        SOURCEIDs of stars to be analyzed
    name_list : list of str
        Names that correspond to each SOURCEID: will be used in
        the filename and title section of each plot, and in the tables
    path : str
        The parent file path to save the output to. 
        Please make sure it has a "/" at the end of it.
    option : list of str, optional
        Which sub-components of this function you'd like to actually
        call. Default is all (lc, tables, phase).

    Returns
    -------
    I don't know yet. Probably nothing.
    Produces a bunch of files, though!

    
    Notes:
    The input data must already be cleaned via some other method.
    
    """

    if path=='' or type(path) is not str:
        print "`path` must be a string. Exiting without action."
        return
    
    # Force `path` to have a trailing forward slash
    path = path.rstrip('/')+'/'
    
    ## First, make directories for everything (or try, at least).

    # Structure of directories:
    #  path/lc
    #  path/phase
    #  path/tables
    # within each of those (now it's a *),
    #  path/*/(s1 | s2 | s3 | s123)
    # Within lc/s*, we put the files directly in.
    # Within phase/s*, we do
    #  path/phase/s*/(h_fx2 | h_lsp | j_fx2 | j_lsp | k_fx2 | k_lsp | lsp_power)
    # Within tables/s*, we put the files directly in.

    mkdir_p(path+"lc/")
    mkdir_p(path+"phase/")
    mkdir_p(path+"tables/")

    ss = ['s1', 's2', 's3', 's123']
    types = ['h_fx2', 'h_lsp', 'j_fx2', 'j_lsp', 'k_fx2', 'k_lsp', 'lsp_power']

    for s in ss:
        mkdir_p(path+"lc/"+s)
        mkdir_p(path+"tables/"+s)
        for t in types:
            mkdir_p("%sphase/%s/%s"%(path,s,t))
            
    # We should now be done making directories. Let's test this.
    # Tested! Woo.

    ## Second, let's make tables.

    if 'tables' in option:

        tables = path+"tables/"

        # Make a lookup table. It should have a column with the `name_list`
        # parameter we fed into this function, as well as a column with the 
        # SOURCEIDs.
        # They should be named "SOURCEID" and "Designation", respectively.
    
        lookup = atpy.Table()
        lookup.add_column("SOURCEID", sid_list)
        lookup.add_column("Designation", name_list)
        
        for season, s in zip([1,2,3,123], ss):
            
            # Write the spreadsheet and save it to the relevant directory.
            spreadsheet.spreadsheet_write(table, lookup, season, 
                                          tables+s+'/spreadsheet.fits', 
                                          flags=16, per=True)
            
        # Tested!
        #return

    # What command do we want to make plots?
    # Probably plot3.lc and plot3.phase, which are going to be almost 
    # identical to plot2 equivalents except that they can handle missing data 
    # and flags and stuff.

    ## Third, make lightcurves.
    # And put the gorram Stetson index in the title!

        
    for season, s in zip([1,2,3,123], ss):
        for name, sid in zip(name_list, sid_list):
            # The specific plot command we use here depends a lot
            # on what functions are available.
            tplot.lc(table, sid, season=season, name=name, #flags=16,
                     outfile=path+"lc/"+s+"/"+name, png_too=True) #png, eps, pdf


        # A bunch of the following code will be substantially rewritten
        # once "spread3" is functional.
            for t in types:
                if t == 'lsp_power':
                    tplot.lsp_power(table, sid, season=season, name=name,
                                    outfile=path+"phase/"+s+"/lsp_power/"+name, 
                                    png_too=True) 
                elif t == 'k_fx2':
                    tplot.phase(table, sid, season=season, name=name,
                                outfile=path+"phase/"+s+"/k_fx2/"+name, 
                                png_too=True) 
                elif t == 'k_lsp':
                    tplot.phase(table, sid, season=season, name=name,
                                period='lsp',
                                outfile=path+"phase/"+s+"/k_lsp/"+name, 
                                png_too=True) 

                else:
                    pass


    return


 
