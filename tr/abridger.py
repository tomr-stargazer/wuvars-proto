""" 
abridger.py : a module to handle how to make abridged lightcurves

"""

def abridger( s_table, flags=256 ):
    """
    A function that "intelligently" calculates 'abridging' parameters.

    Used as a subroutine for plot3.py lightcurve functions.
    Does these things:

    1. Tells you exactly how much to subtract (and where), and 
       where to put dotted lines, given season breaks.
    2. Gives you arrays to use for xticks and xticklabels.

    Parameters
    ----------
    s_table : atpy.Table
        Table with time-series photometry of ONE STAR please
    flags : int, optional 
        Maximum ppErrBit quality flags to use (default 0)    

        
    Returns
    -------
    
    """
    
    s1_start = 54034
    s1_s2_bound = 54300
    s2_s3_bound = 54600

    # extract the data that beat "flags" from s_table; 
    # that's our "working table"
    
    w_table = s_table.where( (s_table.JPPERRBITS < flags) &
                             (s_table.HPPERRBITS < flags) &
                             (s_table.KPPERRBITS < flags) )

    dates = w_table.MEANMJDOBS

    s1_xmin = dates.min()
    s1_xmax = dates[dates < s1_s2_bound].max()
    
    s2_xmin = dates[dates > s1_s2_bound].min()
    s2_xmax = dates[dates < s2_s3_bound].max()


    # set the locations of the xticks
    # xticks( arange(6) )
    
    # set the locations and labels of the xticks
    # xticks( arange(5), ('Tom', 'Dick', 'Harry', 'Sally', 'Sue') )
  
