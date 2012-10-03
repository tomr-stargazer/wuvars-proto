''' 
This package interfaces with Palmer's runchi2 program, which
is found at http://public.lanl.gov/palmer/fastchi.html and described in
2009ApJ...695..496P.

Useful functions:
  chi_analyze - returns best frequency for one source in a table
  test_analyze - returns best frequency for raw input arrays t, x, xerr

'''

# This is the correct one.

# Current major issues: This code does not recognize error messages given by 
# runchi2! If something goes wrong everything crashes and burns. 
# Update (2 Oct '12) : The above issue never comes up in practice.

import subprocess

import numpy

from tr_helpers import season_cut


def chi_input_writer (name, t, x, err, outfile):
    """ 
    Writes data to a table suitable for runchi2's input. Returns outfile.
    
    Parameters
    ----------
    name : str
        The name of the source.
    t, x, err : array-like
        Arrays for the time values, x-values, and error-bar values to 
        use in the period-finding.
    outfile : str
        Path to write an output file to.

    Returns
    -------
    outfile : str
        Same value as input `outfile`.

    """

    if not (t.size == x.size == err.size):
        print "Input arrays must be the same size!"
        return None

    f = open(outfile,'w')
    f.write(name+"\n")
    f.write(str(t.size)+"\n") 
    
    for i in range(t.size):
        f.write("%f \t %f \t %f \n" % (t[i], x[i], err[i]))

    f.close()

    return outfile


def smart_chi_writer (table, sid, band = 'j', season=123) :
    """
    Writes one source's data into runchi2's format.

    Parameters
    ----------
    table : atpy.Table
        Table of time-series data from WFCAM.
    sid : int
        13-digit WFCAM source ID.
    band : {'j', 'h', 'k', 'jmh', 'hmk'}
        Which band to run period-finding over.
    season : {0, 1, 2, 3, 123}
        which season to select (1, 2, 3, 123, or other=no cut)

    Returns
    -------
    outfile : str
        An outfile path based on the name of the source.

    """

    print band

    name = "." + band + str(sid)
#    w = numpy.where(table.SOURCEID == sid) # I may alter this line to filter error 
                                     # bits and such
    s_table = season_cut(table, sid, season)

    if 'm' in band:
        tail = "PNT"
    else:
        tail = "APERMAG3"

    t = s_table.MEANMJDOBS
    x = s_table.data[band.upper()+tail]
    err=s_table.data[band.upper()+tail+"ERR"]

    # do i want to use a data path for the outfile? ... yes

    outfile = '/home/trice/reu/DATA/chi2/' + name
    return chi_input_writer (name, t, x, err, outfile) 

def run_chi (infile, diagnostic=False) :
    """ 
    Runs Palmer's runchi2 program on a given input file. 

    Parameters
    ----------
    infile : str
        An input file for runchi2, formatted to the standards 
        specified in /home/tom/reu/software/FastChi2-1.03/README.
        Created by chi_input_writer() or, if you enjoy punshment,
        by hand.
    diagnostic : bool or str, optional (default False)
        "diagnosticfile contains the chi-squared reduction 
        (larger is better) for each frequency"

    Returns
    -------
    result : str
        The "stdout" output from runchi2, which needs to be parsed by
        parse_chi().
       
    """
    # These arguments, after "runchi2", correspond to 
    # `nharmonics`, `freqmax`, and `infile` (with its "-i" flag) respectively.
    # See /home/tom/reu/software/FastChi2-1.03/README for more details.
    args = ["runchi2","3","12","-i",infile]

    # If we want to get a diagnosticfile, we'll add an argument:
    if diagnostic and type(diagnostic) is str:
        args.extend(["-D", diagnostic])

    # Here I  pass the output of runchi2 to stdout, 
    # then python's subprocess call returns that stdout 
    # string, which can later be parsed by parse_chi(). 

    runchi = subprocess.Popen(args, stdout=subprocess.PIPE)

    result = runchi.communicate()[0]

    return result


def parse_chi (string, ret_chimin=False) :
    """
    Parses the output of runchi2 for one source, returns frequency 

    Parameters
    ----------
    string : str
        A single string containing the text output of runchi2,
        called with only the -i option.
    ret_chimin : bool, optional (default False)
        Return a chimin value?
        If True, then the return value is a tuple of two floats.
        If False, then the return value is a float.

    Returns
    -------
    fbest : float
        The best-fit period returned by runchi2. If the parsing fails, 
        then a value of -1 is returned.
    chimin : float, only if ret_chimin=True
        The minimum chisquared value at `fbest`
    
    """
    # maybe it would be a better idea to define this earlier and have
    # run_chi call this function before returning. The big sloppy string is
    # never ever useful until it's been parsed...

    # and then parse the string from there, and thus return
    # a tuple or _dictionary_ of important values such as periods.

    # Okay, we're gonna have to do some stripping and some parsing (splitting)

    # The following assumes that the data is on the second-to-last line
    # (the last line is el blanco) ... (I mean blank)
    results = string.split('\n')[-2].split('\t')

#    print results
    
    name = results[0]

    # THE FOLLOWING WAS A 6AM HACK
    try:
        fbest= float(results[1])
        chimin = float(results[3]) #tom go check which one is chimin
    except:
        fbest= -1
        chimin = -1

    if ret_chimin:
        return fbest, chimin
    else:
        return fbest
    
def chi_analyze (table, sid, band = 'j', season=123) :
    """
    Does everything for one source and returns the frequency information.

    Parameters
    ----------
    table : atpy.Table
        Table of WFCAM time-series photometry
    sid : int
        13-digit WFCAM source ID of star to analyze
    band : {'j', 'h', 'k', 'jmh', 'hmk'}, optional
        Which band to run period-finding over.
    season : {0, 1, 2, 3, 123}, optional
        which season to select (1, 2, 3, 123, or other=no cut)

    Returns
    -------
    fbest : float
        The star's best frequency at that band (or -1 if failure)

    """

    datafile = smart_chi_writer(table,sid,band, season)
    
    return parse_chi ( run_chi ( datafile ) )

def test_analyze (t, x, err, ret_chimin=False):
    """
    Takes in test data, returns best frequency 

    Parameters
    ----------
    t, x, err : array-like
        Arrays for the time values, x-values, and error-bar values to 
        use in the period-finding.
    ret_chimin : bool, optional (default False)
        Return a chimin value?

    Returns
    -------
    fbest : float
        The best-fit period returned by runchi2. If the parsing fails, 
        then a value of -1 is returned.
    chimin : float, only if ret_chimin=True
        The minimum chisquared value at `fbest`
        
    """

    datafile = chi_input_writer("test", t, x, err, 
                                '/home/trice/reu/DATA/chi2/chitest.in')

    return parse_chi( run_chi(datafile), ret_chimin=ret_chimin )

# Next step... combine stat and chi2 to make a table of best-freq and chi^2 
# values for every source. I've been assuming that chi^2 can be used as a 
# measure of how good the periodicity actually is.

def parse_diagnostic( diagnostic ) :
    """
    Loads in the diagnostic file and returns freq. and chisq arrays.

    Parameters
    ----------
    diagnostic : str
        Filepath of the input diagnostic file.
    
    Returns
    -------
    freq : array-like
        An array of frequencies tested by Fast chi-squared.
    chisq_red : array-like
        An array of the chi-squared reduction (larger is better) 
        for each frequency.

    """
    
    # 9 rows are skipped because of the form of the input file.
    # If this varies, I may need to do this more intelligently.
    contents = np.loadtxt(diagnostic, skiprows=9)

    freq = contents[:,0]
    chisq_red = contents[:,1]

    return freq, chisq_red


def diagnostic_analyze(t, x, err):
    """
    Takes in test data, returns fx2 periodogram. 
    
    Note: much of this is of the same form as test_analyze(), and 
    could in principle be easily combined with that function.

    Parameters
    ----------
    t, x, err : array-like
        Arrays for the time values, x-values, and error-bar values to 
        use in the period-finding.

    Returns
    -------
    freq : array-like
        An array of frequencies tested by Fast chi-squared.
    chisq_red : array-like
        An array of the chi-squared reduction (larger is better) 
        for each frequency.

    """
    
    input_file = '/home/trice/reu/DATA/chi2/chitest.in'
    diag_file = '/home/trice/reu/DATA/chi2/diagnostic.txt'

    # This line is the same as for test_analyze()
    datafile = chi_input_writer("test", t, x, err, input_file)

    # do a runchi2, but ignore the piped output; we just want diag_file
    run_chi(datafile, diagnostic=diag_file)
    
    return parse_diagnostic( diag_file )
