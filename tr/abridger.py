""" 
abridger.py : a module to handle how to make abridged lightcurves

"""

import numpy as np

def abridger( s_table, date_offset, flags=256 ):
    """
    A function that "intelligently" calculates 'abridging' parameters.

    Used as a subroutine for plot3.py lightcurve functions.
    Does these things:

    1. Tells you exactly how much to subtract (and where), and 
       where to put dotted lines, given season breaks.
    2. Gives you arrays to use for xticks and xticklabels.
    3. Tells you where to put your xlim()s.

    Parameters
    ----------
    s_table : atpy.Table
        Table with time-series photometry of ONE STAR please
    date_offset : atpy.Table
        MJD value to use as "zero" date.
    flags : int, optional 
        Maximum ppErrBit quality flags to use (default 0)    
        
    Returns
    -------
    s2_subtraction_factor, s3_subtraction_factor : float
        What values to subtract for (1) all data after the end of season
        1, and (2) all data after the end of season 2. Note that season 3 
        will need to have both subtractions performed on it.
    s1_s2_line, s2_s3_line : float
        The x-values to plot vertical dotted separator lines
    xticks : array-like
        X-values to place xticks at. 
        Generally separated by increments of 50; season boundaries are funny.
    xticklabels : list of str
        Labels for the above xticks
    xlim_bounds : tuple of float
        Values to use for the min and max xlim().
    
    """
    
    # how much buffer space, in nights, between season separators 
    # and on the edges?
    spacing = 5

    s1_start = 54034 - date_offset
    s1_s2_bound = 54300 - date_offset
    s2_s3_bound = 54600 - date_offset

    # extract the data that beat "flags" from s_table; 
    # that's our "working table"
    
    w_table = s_table.where( (s_table.JPPERRBITS < flags) &
                             (s_table.HPPERRBITS < flags) &
                             (s_table.KPPERRBITS < flags) )

    dates = w_table.MEANMJDOBS - date_offset

    # define the precise boundaries of each season
    s1_xmin = dates.min() 
    s1_xmax = dates[dates < s1_s2_bound].max() 
    
    s2_xmin = dates[dates > s1_s2_bound].min() 
    s2_xmax = dates[dates < s2_s3_bound].max() 

    try:
        s3_xmin = dates[dates > s2_s3_bound].min() 
    except ValueError:
        s3_xmin = dates.max()
    s3_xmax = dates.max() 

    s2_subtraction_factor = s2_xmin - s1_xmax - 2*spacing
    s3_subtraction_factor = s3_xmin - s2_xmax - 2*spacing

    # now determine where to place xticks and what to call them
    s1_xticks = np.arange(int(s1_xmin), int(s1_xmax), 50)
    s1_xticklabels = [ str(x) for x in s1_xticks ]
    
    s2_xticks_raw = np.arange(int(s2_xmin), int(s2_xmax), 50)
    s2_xticks = s2_xticks_raw - int(s2_subtraction_factor)
    s2_xticklabels = [ str(x) for x in s2_xticks_raw ]
 
    s3_xticks_raw = np.arange(int(s3_xmin+10), int(s3_xmax), 50)
    s3_xticks = s3_xticks_raw-int(s3_subtraction_factor+s2_subtraction_factor)
    s3_xticklabels = [ str(x) for x in s3_xticks_raw ]

    # combine the respective xtick and xticklabel arrays together!

    xticks = np.concatenate( (s1_xticks, s2_xticks, s3_xticks) )
    xticklabels = s1_xticklabels + s2_xticklabels + s3_xticklabels

    # now determine where to put the dotted line spacings.
    s1_s2_line = s2_xmin - spacing - s2_subtraction_factor
    s2_s3_line = s2_xmax + spacing - s2_subtraction_factor
    
    # and how to size the plot, overall
    xlim_bounds = (s1_xmin-spacing, 
                   s3_xmax+spacing-s2_subtraction_factor-s3_subtraction_factor)

    # set the locations of the xticks
    # xticks( arange(6) )
    
    # set the locations and labels of the xticks
    # xticks( arange(5), ('Tom', 'Dick', 'Harry', 'Sally', 'Sue') )
  
    return (s2_subtraction_factor, s3_subtraction_factor,
            s1_s2_line, s2_s3_line, xticks, xticklabels, xlim_bounds)
