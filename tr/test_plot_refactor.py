"""
This is a module in which I'll write "visual tests" for my
plot3.py refactoring.

"""

from variables_data_filterer import variables_photometry, ukvar_spread
import plot3


plot3.lc(variables_photometry, ukvar_spread.SOURCEID[0], 
         name='lc: default settings')

plot3.graded_lc(variables_photometry, ukvar_spread.SOURCEID[0], 
                name='graded_lc: default settings')

plot3.lc(variables_photometry, ukvar_spread.SOURCEID[0], color_slope=True,
         name='lc: color_slope=True')

plot3.graded_lc(variables_photometry, ukvar_spread.SOURCEID[0], 
                color_slope=True,
                name='graded_lc: color_slope=True')
