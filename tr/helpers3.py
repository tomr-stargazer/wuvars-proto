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

def data_cut (table, sid_list, season=0):
    """
    Selects data corresponding to specified source(s).

    Returns a subset of `table` containing only data for sources
    specified in `sid_list`, for dates within `season`.

    Parameters
    ----------
    table : atpy.Table
        An ATpy Table with time-series photometry from the WFCAM 
        Science Archive. It does not need to necessarily have 
        well-defined J,H,K data for each timestamp (unlike tables 
        fed into tr_helpers.data_cut, a similar function).
    sid_list : array_like
        A list of 13-digit WFCAM Source IDs whose data to extract.
    season : int, optional
        Which observing season of our dataset (1, 2, 3, or all).
        Any value that is not the integers (1, 2, or 3) will be 
        treated as "no season", and no time-cut will be made.
        Note that this is the default behavior.

    Returns
    -------
    cut_table : atpy.Table
        Subset of `table` containing only data for sources
        specified in `sid_list`, for dates within `season`.

    """

    
    

    
