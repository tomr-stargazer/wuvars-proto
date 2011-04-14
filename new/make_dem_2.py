''' A one-time script to make a stats table (no period information) from 
my WSERV more-or-less complete catalog! it should probably run overnight.'''

import atpy
import filter as ft
from datetime import datetime
import spreadsheet

path = '/media/storage/Documents/Research/reu/DATA/Merged_Catalogs/'

w_clipped = atpy.Table(path+'wserv_happy.fits')


print "loaded stuff at " + datetime.now().strftime("%a, %d %b %Y %H:%M:%S")


mock_lookup = spreadsheet.base_lookup(w_clipped)

print "mock lookup at " + datetime.now().strftime("%a, %d %b %Y %H:%M:%S")

global_stats = spreadsheet.spreadsheet_write ( w_clipped, mock_lookup, 0, 
                                               path+"global_stats1.fits")

print "made stats at " + datetime.now().strftime("%a, %d %b %Y %H:%M:%S")+"!!!"


