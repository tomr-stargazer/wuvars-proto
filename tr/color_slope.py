"""
color_slope.py : a module to fit slopes to color trajectories.

"""

from scipy.odr import RealData, Model, ODR

from helpers3 import data_cut, band_cut

def f(B, x):
    ''' Linear function y = m*x + b '''
    return B[0]*x + B[1]

    # B is a vector of the parameters.
    # x is an array of the current x values.
    # x is same format as the x passed to Data or RealData.

    # Return an array in the same format as y passed to Data or RealData.

linear = Model(f)

def slope(x, y, xerr, yerr, verbose=True):
    """
    Calculates the slope of a color trajectory.

    Uses the scipy ODRpack and total least squares algorithm.

    Parameters
    ----------
    x, y : array-like
        Colors or magnitudes to calculate slopes of.
    xerr, yerr : array-like
        Corresponding errors (stdev) on `x` and `y`.
        Be sure to make sure that the HMKPNT errors are properly adjusted!
    verbose : bool. optional
        Whether to print a verbose output. Default true.

    Returns
    -------
    slope : float
        Slope (in rise/run) of the linear fit.
    intercept : float
        Y-value where the linear fit intercepts the Y-axis.
    slope_error : float
        The standard error on the fitted slope: an indication of fit quality.

    Notes
    -----
    Much of this code is borrowed from the dosctring example found
    at https://github.com/scipy/scipy/blob/master/scipy/odr/odrpack.py#L27 .
    
    See http://docs.scipy.org/doc/scipy/reference/odr.html and
    http://stackoverflow.com/questions/9376886/orthogonal-regression-fitting-in-scipy-least-squares-method
    for further discussion.

    """
    
    mydata = RealData(x, y, sx=xerr, sy=yerr)

    # Someday, we may want to improve the "initial guess" with 
    # a leastsq first-pass.
    myodr = ODR(mydata, linear, beta0=[1.,2.])

    myoutput = myodr.run()

    print myoutput.pprint()

    return myoutput.beta[0], myoutput.beta[1], myoutput.sd_beta[0]

def star_slope( table, sid, xband='hmk', yband='k', flags=0, verbose=True):
    """
    Calculates the color slope, given an input table and ID.

    Parameters
    ----------
    table : atpy.Table
        Table with time-series photometry
    sid : int
        13-digit WFCAM source ID of star to plot
    xband : {'jmh', 'hmk'}
        The x-axis array to use for the slope. Default 'k'.
    yband : {'j', 'jmh', 'hmk'}
        The y-axis array to use for the slope. Default 'hmk'.
    flags : int, optional 
        Maximum ppErrBit quality flags to use (default 0).
    verbose : bool. optional
        Whether to print a verbose output. Default true.
   
    Returns
    -------
    slope : float
        Slope (in rise/run) of the linear fit.
    intercept : float
        Y-value where the linear fit intercepts the Y-axis.
    slope_error : float
        The standard error on the fitted slope: an indication of fit quality.

    """

    if (xband not in ['jmh', 'hmk']) or (yband not in ['j', 'jmh', 'k']):
        raise ValueError("Incorrect argument to `xband` or `yband`")

    # define this dict thing
    band_dict = {'j':'JAPERMAG3', 'k':'KAPERMAG3', 
                 'jmh':'JMHPNT', 'hmk':'HMKPNT'}


    # Loading data
    s_table = data_cut (table, sid, season=0)

    if len(s_table) == 0:
        print "no data here"
        return


    j_table = band_cut( s_table, 'j', max_flag=flags)
    h_table = band_cut( s_table, 'h', max_flag=flags)
#    k_table = band_cut( s_table, 'k', max_flag=flags)
    
    jh_table = band_cut( j_table, 'h', max_flag=flags)
    hk_table = band_cut( h_table, 'k', max_flag=flags)
 #   jk_table = band_cut( j_table, 'k', max_flag=flags)

    jhk_table = band_cut( jh_table, 'k', max_flag=flags)
    

    if (xband, yband) == ('hmk', 'k'):
        data = hk_table
    elif (xband, yband) == ('jmh', 'j'):
        data = jh_table
    elif (xband, yband) == ('hmk', 'jmh'):
        data = jhk_table
    else:
        data = jhk_table
        print "Incorrect combination of `xband`, `yband`."

    x_array = data.data[band_dict[xband]]
    xerr_array = data.data[band_dict[xband]+"ERR"]
    y_array = data.data[band_dict[yband]]
    yerr_array = data.data[band_dict[yband]+"ERR"]
    
    return slope(x_array, y_array, xerr_array, yerr_array, verbose)
