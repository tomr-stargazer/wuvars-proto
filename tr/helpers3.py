"""
This is a module that contains 'helper' functions that are called 
by the other '3'-series packages in wuvars. It does not import any 
of my other modules, for dependency reasons.

Useful functions:
  data_cut - Cuts a table for a selection of sources and seasons
  band_cut - Selects only data where a certain band (J,H,K) is well-defined.

"""

import numpy as np
import atpy

def data_cut (table, sid_list, season=0, null=-9.99999e+08):

    pass
    
