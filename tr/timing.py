'''
This package has code to filter out bad periods from periodogram 
analysis, and to return the best filtered period.

Useful functions:
  lsp_mask - Trims unreliable frequencies from a periodogram
  lsp_tuning - Selects the proper `hifac` value for a desired highest frequency.
'''

from __future__ import division
import numpy as np


''' 
Okay, what periods (frequencies) are we removing on our first run?

Exclude periods above 500 days,
# exclude < 0.002
between .98, 1.02 days
# exclude  < 1.0526315789473684
# and > 0.95238095238095233
Use 13 hours min
# exclude > 1.8461538461538463

'''

upper_f = 1 / (13/24) # 13 hours is the highest frequency (shortest period)
lower_f = 1 / (500) # 500 days is the lowest frequency (longest period)
midrange = [ (1/1.10, 1/.90) ] # a list of mid-ranges to exclude


def lsp_mask ( Wk1, Wk2, upper_f=upper_f, lower_f=lower_f, midrange=midrange ):
    ''' Trims unreliable frequencies from a periodogram and returns index of
    the reliably best frequency.

    Inputs:
      Wk1 -- An array of Lomb periodogram frequencies.
      Wk2 -- An array of corresponding values of the Lomb periodogram.

    Optional inputs:
      upper_f -- The upper frequency to cut (shortest period)
      lower_f -- The lower frequency to cut (longest period)
      midrange -- A list of pairs of frequencies to exclude

    Returns:
      Jmax -- The array index corresponding to the MAX( Wk2 ), but trimmed!

# deprecated:
#      fWk1 -- The trimmed frequencies.
#      fWk2 -- The corresponding trimmed periodogram values.

      '''


    trim = np.where( (Wk1 < upper_f) & (Wk1 > lower_f) &
                     ( (Wk1 < midrange[0][0]) | (Wk1 > midrange[0][1]) ) )

    print "LSP mask executed properly."
    Jmax = np.where( Wk2 == Wk2[trim].max() )[0][0]
    return Jmax



def lsp_tuning(t, upper_frequency=0.5):
    """
    Tunes the `hifac` value in a periodogram to a desired upper frequency.

    By default, it'll tell you what `hifac` value to use so that your 
    periodogram starts at a 2 day period on the lower end.
    
    Note: 2 days is the empirically-derived (and unsurprising) effective 
    Nyquist frequency (0.5 day^-1) for our dataset, based on looking at lots
    of aliased periodograms. (and thinking about it.)

    Parameters
    ----------
    t : array-like
        Array of timestamps going into the periodogram
    upper_frequency : float, optional
        Highest frequency (1/lowest period) that one desires to scan over
        in the period search.

    Returns
    -------
    hifac : float
        Suggested value of `hifac` to feed into scargle.fasper().
        Automatically rounded up (np.ceil) for your convenience.

    """

    n = len(t)

    tmin = t.min()
    tmax = t.max()
    tdif = tmax-tmin

    # The "average nyqust frequency" is given by
    # 0.5 * (n / tdif), so we invert that to find the proper hifac:
    hifac = np.ceil(2. * tdif / n * upper_frequency)

    return hifac
