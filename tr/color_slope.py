"""
color_slope.py : a module to fit slopes to color trajectories.

"""

from scipy.odr import RealData, Model, ODR

def f(B, x):
    ''' Linear function y = m*x + b '''
    return B[0]*x + B[1]

    # B is a vector of the parameters.
    # x is an array of the current x values.
    # x is same format as the x passed to Data or RealData.

    # Return an array in the same format as y passed to Data or RealData.

linear = Model(f)

def slope(x, y, xerr, yerr):
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

    Returns
    -------
    slope : float
        Slope (in rise/run) of the linear fit.
    intercept : float
        Y-value where the linear fit intercepts the Y-axis.
    residual : float
        The residual variance of the fit: an indication of fit quality.

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

    return myoutput.beta[0], myoutput.beta[1], myoutput.res_var
