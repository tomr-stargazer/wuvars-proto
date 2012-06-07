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


def do_it_all( table, sid_list, name_list, path ):
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

    Returns
    -------
    I don't know yet. Probably nothing.
    Produces a bunch of files, though!

    
    Notes:
    The input data must already be cleaned via some other method.
    
    """
    
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
        spreadsheet_write(table, lookup, season, tables+ss+'/spreadsheet.fits',
                          per=True)


    # What command do we want to make plots?

    return


 
