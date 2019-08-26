#!/usr/bin/env python3

# -----------------------------------------------------------------------------
# This file is part of the BHC Complexity Toolkit.
#
# The BHC Complexity Toolkit is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# The BHC Complexity Toolkit is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with the BHC Complexity Toolkit.  If not, 
# see <https://www.gnu.org/licenses/>.
# -----------------------------------------------------------------------------
# Copyright 2019, Mark D. Flood
#
# Author: Mark D. Flood
# Last revision: 22-Jun-2019
# -----------------------------------------------------------------------------

__all__ = ['tic',
           'toc',
           'currtime',
           'parse_command_line', 
           'read_config', 
           'print_config',
           'make_asof',
           'assemble_asofs',
           'stringify_qtrend',
           'next_qtrend',
           'rcnt_qtrend',
           ]
__version__ = '0.3'
__author__ = 'Mark D. Flood'


import time
import getopt
import sys
import os
import logging
import logging.config as logcfg
import re
import configparser as cp

LOG = logging.getLogger(__file__.split(os.path.sep)[-1].split('.')[0])



lasttime = time.time()
def currtime():
    """Get the current system time
    
    Returns
    -------
    currtime : float
        The new benchmark system time, from time.time()
    """
    currtime = time.time()
    return currtime

def tic():
    """Reset the timer to a new benchmark time
    
    Returns
    -------
    lastime : float
        The new benchmark system time, from time.time()
    """
    global lasttime
    lasttime = time.time()
    return lasttime

def toc():
    """Reset the benchmark time, and report time elapsed since last reset
    
    Returns
    -------
    delta : float
        Time change since last call to tic() or toc(), in seconds
    """
    global lasttime
    newtime = time.time()
    delta = newtime - lasttime
    lasttime = newtime
    return delta



def parse_command_line(argv, modulefile):
    """Parses command-line arguments and overrides items in the config.
    
    This method assumes that the Python ConfigParser has already read in
    a config object (during module import), and that additional arguments
    (argv) are available as command-line parameters. The following 
    parameters (all optional) are recognized:
        
        * -c   Print the config dictionary for this module to stdout
        * -l <loglevel_file>      Set the logging threshold for file output
        * -L <loglevel_console>   Set the logging threshold for console output
        * -C <configfile>   Read custom configuration from a separate file
        * -h | --help   Print usage help and quit
        * -p <paramkey>:<paramval>  Override/set individual config parameters
    """
    usagestring = ('python '+modulefile+
                  ' [-c]'+
                  ' [-C <configfile>]'+
                  ' [-l <loglevel_file>]'+
                  ' [-L <loglevel_console>]'+
                  ' [-h | --help]'+
                  ' [-p <paramkey>:<paramval>]')
    section = os.path.splitext(os.path.basename(modulefile))[0]
    config = None
    showconfig = False
    if argv is None:
        argv = sys.argv
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "hcl:L:C:p:", ["help"])
        except getopt.error as msg:
            raise Usage(msg)
        for o, a in opts:
            # Scan through once, to see if help was requested
            if o in ("-h", "--help"):
                print(usagestring)
                sys.exit()
        for o, a in opts:
            # Scan through once, to see if a special config is named
            if "-C"==o:
                cfgfile = a
                config = read_config(cfgfile)
        # If we still have no config, attempt to read the default
        if (None==config):
            config = read_config()
        for o, a in opts:
            # Now we have the right config file, get the parameters
            if "-p"==o:
                [paramkey, paramval] = a.split(':')
                config[section][paramkey] = paramval
            if "-C"==o:
                pass
            elif "-c"==o:
                showconfig = True
            elif "-l"==o:
                config['handler_file']['level'] = a
            elif "-L"==o:
                config['handler_console']['level'] = a
            elif o in ("-h", "--help"):
                print(usagestring)
                sys.exit()
            else:
                assert False, "unhandled option: "+o
    except Usage as err:
        print(err.msg, file=sys.stderr)
        print("for help use --help", file=sys.stderr)
        return 2
    if (showconfig):
        print_config(config, modulefile)
    return config
class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg



def read_config(config_file='bhc_complex.ini'):
    """Reads the application parameters from a local configuration file
    
    Parameters
    ----------
    config_file : filename
        The name of a local configuration file. The default value is
        bhc_complex.ini, which should reside in the same directory with
        the Python source code.
    """
    config = cp.ConfigParser(interpolation=cp.ExtendedInterpolation())
    config.read(config_file)
    log_config()
    return config



def log_config(logger_name=None, config_file='logging.ini'):
    """Reads the logging parameters from a local configuration file
    
    Parameters
    ----------
    logger_name : string
        The name of a logger to create. If this parameter is set to None,
        no logger is configured and the function returns None. 
        Else a logger of the given name is returned.
    config_file : filename
        The name of a local configuration file. The default value is
        logging.ini, which should reside in the local (.) directory.
        
    Returns
    -------
    log : logging.Logger
        The logger for the given logger_name, if specified. If 
        logger_name==None, then this function returns None. 
    """
    log = None
    if not(logging.getLogger().hasHandlers()):
        # Logging has not been configured yet
        log_config = cp.ConfigParser(interpolation=cp.ExtendedInterpolation())
        log_config.read(config_file)
        log_dir = log_config['handler_file']['args']
        log_dir = log_dir.split(sep="'")[1]
        log_dir = os.path.split(log_dir)[0]
        os.makedirs(log_dir, exist_ok=True)
        logcfg.fileConfig(config_file)
        # Set level on the root logger to the smallest of the handler levels
        LEVELS = {'NOTSET':logging.NOTSET, 'DEBUG':logging.DEBUG,
                  'INFO':logging.INFO, 'WARNING':logging.WARNING,
                  'ERROR':logging.ERROR, 'CRITICAL':logging.CRITICAL}
        loglevel_file = log_config['handler_file']['level']
        loglevel_console = log_config['handler_console']['level']
        loglevel_min = min(LEVELS[loglevel_file], LEVELS[loglevel_console])
        logging.getLogger().setLevel(loglevel_min)
    if not(logger_name is None):
        # Logging is configured, but we need to add this logger
        log = logging.getLogger(logger_name)
    return log



def print_config(config, modulefile):
    """Prints the parameters for a specific configuration section
    
    Simple formatted dump of the config parameters relevant for a given 
    configuration section. This is useful to confirm that the app has 
    ingested the correct configuation.
    
    Parameters
    ----------
    config : configparser.ConfigParser
        The configuration object to display
    modulefile : str
        The module name whose configuration section should be displayed

    """
    section = os.path.splitext(os.path.basename(modulefile))[0]
    print('-------------------------- CONFIG ----------------------------')
    print('Current working directory:', os.getcwd())
    print('[DEFAULT]')
    config_dict = dict(config['DEFAULT'])
    for k,v in sorted(config_dict.items()):
        print(k, '=', v)
    print('['+section+']')
    config_dict = dict(config[section])
    for k,v in sorted(config_dict.items()):
        print(k, '=', v)
    print('--------------------------------------------------------------')



def delim_norm(raw_delim):
    """Normalizes a delimiter specifier to convert tab characters
    
    The Python ConfigParser does not handle tab characters well, so 
    we specify those with the string '<TAB>' instead. This function
    converts any instances of '<TAB>' in a raw_delim string to '\t'
    instead. If the raw_delim is a comma (or any other character or
    string without a '<TAB>' substring), this function has no effect.
    
    Parameters
    ----------
    raw_delim : str
        The candidate delimiter string to convert
    
    Returns
    -------
    delim : str
        A string normalized to convert '<TAB>' to '\t'
    """
    delim = re.sub('<TAB>', '\t', raw_delim)
    return delim




def make_asof(YYYYQQ):
    """Converts a string in YYYYQQ form to an integer in yyyymmdd form

    The function parses the string and returns a tuple of date components.
    
    Parameters
    ----------
    YYYYQQ : string
        A string of form YYYYQQ, for example, '1995Q3'
        
    Returns
    -------
    tuple
       *  yyyymmdd:  An int variable indicating year/mo/day, 19950930
       *  y:         An int variable indicating the year, 1995
       *  q:         An int variable indicating the quarter, 3
       *  Q:         A string variable indicating the quarter, 'Q3'
       
    Examples
    --------
    >>> KeyDate = '1999Q2'
    >>> make_asof(KeyDate)
    19990630
    """
    MMDDs = [331, 630, 930, 1231]
    y = int(YYYYQQ[0:4])
    q = int(YYYYQQ[5:6])
    Q = YYYYQQ[4:6]
    yyyymmdd = y*10000+MMDDs[q-1]
    return (yyyymmdd, y, q, Q)
 


def assemble_asofs(YQ0, YQ1):
    """Converts quarter (YYYYQQ) start/end dates to a list of ints (YYYYMMDD)
    
    Parses the strings YQ0 and YQ1, each of the form YYYYQQ (e.g., "1986Q2")
    into asofdate variables (of type int), each of the form YYYYMMDD 
    (e.g., 19860630). Every quarter-end asofdate between YQ0 and YQ1 
    (inclusive) is added to the asofs list, which is returned. 
    """
    asofs = []
    (yyyymmdd0, Y0, q0, Q0) =  make_asof(YQ0)
    (yyyymmdd1, Y1, q1, Q1) = make_asof(YQ1)
    if (yyyymmdd0 > yyyymmdd1):
        print('ERROR: End date,', yyyymmdd1, 'precedes start date,', yyyymmdd0)
    if (Y0 == Y1):
        # Full range is within one year
        for q in range(q0,q1+1):
            asofs.append(make_asof(str(Y0)+'Q'+str(q))[0])
    else:
        # For the (possibly partial) first year in the range
        for q in range(q0,5):
            asofs.append(make_asof(str(Y0)+'Q'+str(q))[0])
        # For the interior (full) years in the range
        for y in range(Y0+1, Y1):
            for q in range(1,5):
                asofs.append(make_asof(str(y)+'Q'+str(q))[0])
        # For the (possibly partial) last year in the range
        for q in range(1,q1+1):
            asofs.append(make_asof(str(Y1)+'Q'+str(q))[0])
    return asofs


def stringify_qtrend(asofdate):
    """Converts an as-of date to a YYYYQQ string for the next quarter end
    """
    yyyy = int(asofdate/10000)
    mmdd = int(asofdate) -yyyy*10000
    if (930 < mmdd):
        Nqtr = str(yyyy)+'Q4'
    elif (630 < mmdd):
        Nqtr = str(yyyy)+'Q3'
    elif (331 < mmdd):
        Nqtr = str(yyyy)+'Q2'
    elif (100 < mmdd):
        Nqtr = str(yyyy)+'Q1'
    return Nqtr
        
    
def next_qtrend(asofdate):
    """Constructs the next quarter-end date for a given as-of date
    """
    yyyy = int(asofdate/10000)
    mmdd = int(asofdate) -yyyy*10000
    if (1231 == mmdd):
        MRqtr = asofdate
    elif (930 < mmdd):
        MRqtr = yyyy*10000 + 1231
    elif (630 < mmdd):
        MRqtr = yyyy*10000 + 930
    elif (331 < mmdd):
        MRqtr = yyyy*10000 + 630
    elif (100 < mmdd):
        MRqtr = yyyy*10000 + 331
    return MRqtr

    
def rcnt_qtrend(asofdate):
    """Constructs the most recent past quarter-end date for a given as-of date
    """
    yyyy = int(asofdate/10000)
    mmdd = int(asofdate) -yyyy*10000
    if (1231 == mmdd):
        MRqtr = asofdate
    elif (930 <= mmdd):
        MRqtr = yyyy*10000 + 930
    elif (630 <= mmdd):
        MRqtr = yyyy*10000 + 630
    elif (331 <= mmdd):
        MRqtr = yyyy*10000 + 331
    elif (100 < mmdd):
        MRqtr = (yyyy-1)*10000 + 1231
    return MRqtr

    
def rcnt_midyear(asofdate):
    """Constructs the most recent prior mid-year date for a given as-of date
    """
    yyyy = int(asofdate/10000)
    mmdd = int(asofdate) -yyyy*10000
    if (mmdd > 630):
        midyear = yyyy
    else:
        midyear = yyyy-1
    return midyear



if __name__ == "__main__":
    import doctest
    doctest.testmod()