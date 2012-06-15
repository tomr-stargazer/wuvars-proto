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
        fed into tr_helpers.data_cut(), a similar function).
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

    # 0. make sure input data is a list
    if type(sid_list) is not list:
        sid_list = [sid_list]
    
    # First, select these sources' data from the table.
    # This code taken from tr_helpers.data_cut.
    source = table.where( 
        np.array( [sid in sid_list for sid in table.SOURCEID] )
        )
    
    # Second, slice the data by season.
    # These seasons defined by the Cyg OB7 variability study.
    offset = 54579
    cut1 = 100
    cut2 = 300
    cut3 = 600

    if season == 1:
        low = offset
        high = offset+cut1
    elif season == 2:
        low = offset+cut1
        high = offset+cut2
    elif season == 3:
        low = offset+cut2
        high = offset+cut3
    else:
        low = 0
        high = 1e6
    
    cut_table = source.where( (source.MEANMJDOBS < high) & 
                              (source.MEANMJDOBS > low))
    
    return cut_table


def band_cut (table, band, min_flag=0, max_flag=2147483648,
              null=np.double(-9.99999488e+08)): 
    """
    Selects data corresponding to a specified filter.

    Returns a subset of `table` containing only data where photometry in 
    filter `band` is well-defined and the ppErrBits are between `min_flag`
    and `max_flag`. Assumes an input table that has already 
    been processed by data_cut().
    
    Note that `min_flag` and `max_flag` do NOTHING currently if `band` is 
    'jmh' or 'hmk', because the colors don't have their own flags.
    (The right solution is probably to take the relevant J and H or H and K
    flags and, say, add them, but that's extra work.)

    Parameters
    ----------
    table : atpy.Table
        An ATpy Table containing time-series photometry from the 
        WFCAM Science Archive corresponding to a single source.
        (i.e., a table that has been processed by data_cut.)
    band : str ('j'|'h'|'k')
        Which band to select: J, H, or K.
    min_flag : int, optional
        The lowest ppErrBits flag to accept. Default 0.
    max_flag : int, optional
        The highest ppErrBits flag to accept. Default 2147483648 (2**31).
    null : float, optional
        What value to use as a 'null' when filtering data.
        Default value -9.99999e+08 (as used by WSA).

    Returns
    -------
    cut_table : atpy.Table
        Subset of `table` containing only data where photometry in 
        filter `band` is well-defined. 

    """

    # Input check.
    if band.lower() not in ['j', 'h', 'k', 'jmh', 'hmk']:
        raise Exception(
            "Invalid `band`: must use 'j', 'h', 'k', 'jmh', or 'hmk'.")

    # Distinguish between magnitude or color.
    if len(band) == 1:
        metric = "APERMAG3"
    else:
        metric = "PNT"

    band_name = band.upper() + metric
    banderr_name = band_name + "ERR"
    pperrbits_name = band.upper() + "PPERRBITS"

    # Now let's select only data where these guys aren't null
    # and the ppErrBits flags are within the acceptable range.
    # (The ppErrBits part is only currently applicable to J, H, or K.)
    
    if len(band) == 1:
        cut_table = table.where( (table.data[band_name] != null) &
                                 (table.data[banderr_name] != null) &
                                 (table.data[pperrbits_name] > min_flag) &
                                 (table.data[pperrbits_name] < max_flag) )
    else:
        cut_table = table.where( (table.data[band_name] != null) &
                                 (table.data[banderr_name] != null) )
        

    return cut_table

    
