'''
This is the package with super-functions to generate 
hundreds of figures at a time.
It calls upon most of the other wuvars packages.

Useful functions:
 # Primary functions
 
 # Helper functions
'''

import atpy
import numpy as np
import matplotlib.pyplot as plt


import os, errno

def mkdir_p(path):
    ''' A helper function, copied from the following URL:
    http://stackoverflow.com/questions/600268/mkdir-p-functionality-in-python
    '''
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST:
            pass
        else: raise

def do_it_all( table, sid, path ):
    ''' 
    Does some stuff. Not sure exactly what yet, but I'll 
    want it to make tons of plots and tables for a list of 
    SOURCEIDs and their corresponding data.
    '''

    pass
