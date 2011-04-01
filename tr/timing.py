'''
Code to filter out bad periods from the periodograms, 
and then return the best period!
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
midrange = [ (1/1.05, 1/.95) ] # a list of mid-ranges to exclude


def lsp_mask ( Wk1, Wk2, upper_f=upper_f, lower_f=lower_f, midrange=midrange ):
    ''' Trims unreliable frequencies from a periodogram and returns the
    reliably best frequency.

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
                     ( (Wk1 > midrange[0][0]) | (Wk1 < midrange[0][1]) ) )

    Jmax = np.where( Wk2 == Wk2[trim].max() )[0]
    return Jmax

