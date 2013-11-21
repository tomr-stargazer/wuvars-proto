"""
This is a module in which I'll write "visual tests" for my
plot3.py refactoring.

"""

from variables_data_filterer import variables_photometry, ukvar_spread
import plot3

def test_1():
    plot3.lc(variables_photometry, ukvar_spread.SOURCEID[0], 
             name='lc: default settings, ONCvar 1')

    pass # will be "new plot"

def test_2():
    plot3.graded_lc(variables_photometry, ukvar_spread.SOURCEID[0], 
                    name='graded_lc: default settings, ONCvar 1')

    pass # will be "new plot"

def test_3():
    plot3.lc(variables_photometry, ukvar_spread.SOURCEID[0], color_slope=True,
             name='lc: color_slope=True, ONCvar 1')

    pass # will be "new plot"    

def test_4():
    plot3.graded_lc(variables_photometry, ukvar_spread.SOURCEID[0], 
                    color_slope=True,
                    name='graded_lc: color_slope=True, ONCvar 1')

    pass # will be "new plot"

def test_5():
    plot3.graded_lc(variables_photometry, ukvar_spread.SOURCEID[7], 
                    color_slope=True,
                    name='graded_lc: color_slope=True, ONCvar 8')
    
    pass # will be "new plot"

def test_6():
    plot3.graded_phase(variables_photometry, ukvar_spread.SOURCEID[6], 
                       color_slope=True, period=7.244539,
                       name='graded_phase: color_slope=True, ONCvar 7')

    pass # will be "new plot"

def test_7():
    plot3.graded_phase(variables_photometry, ukvar_spread.SOURCEID[6], 
                       timecolor='time', period=7.244539,
                       name='graded_phase: timecolor="time", ONCvar 7')

    pass # will be "new plot"
    

tests = [test_1, test_2, test_3, test_4, test_5, test_6, test_7]

for t in tests:
    t()