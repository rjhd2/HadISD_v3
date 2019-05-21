#!/usr/local/sci/bin/python
#*****************************
#
# Precipitation cross check (PCC)
#   Ensure that shorter accumulation periods have lower amounts than longer ones
#
#************************************************************************
#                    SVN Info
#$Rev:: 160                                           $:  Revision of last commit
#$Author:: rdunn                                      $:  Author of last commit
#$Date:: 2018-06-13 18:28:02 +0100 (Wed, 13 Jun 2018) $:  Date of last commit
#************************************************************************
import numpy as np
import scipy as sp
import datetime as dt
# RJHD routines
import qc_utils as utils

#************************************************************************
def pcc_accumulations(times, precips, start, logfile, plots = False, diagnostics = False):
    '''
    Accumulation cross check.  Check that e.g. 1h accumul <= 24h accumul
    
    :param array times: timestamps
    :param array precips: precipitation depth arrays for 1,3,6,9,12,18,24 hourly accumul.
    :param datetime start: DATASTART (for plotting)
    :param file logfile: logfile to store outputs
    :param bool plots: do plots or not
    :param bool diagnostics: extra verbose output

    :returns: flags - locations where flags have been set
    '''
 
    flags = np.zeros(len(times))

    for t,tt in enumerate(times):
        
        if (tt > 0) and (tt < times[-1]):
            # at a coincident timestamp, presume that all accumulations reporting then
            # were for the previous N hours.  Hence longer accumulations should
            # have larger values.

            good_precips = precips[:, t].compressed()

            if len(good_precips) > 1:
                # if sorting into ascending order makes a difference
                #   then one value is out of order.  Flag all.
                for p, pptn in enumerate(good_precips):
                    if sorted(good_precips)[p] != pptn:
                        flags[t] = 1

    nflags = len(np.where(flags != 0)[0])
    utils.print_flagged_obs_number(logfile, "Precipitation Cross Check", "precipitation", nflags, noWrite = diagnostics)

    return flags # pcc_accumulations

#************************************************************************
def pcc(station, flag_col, start, end, logfile, diagnostics = False, plots = False):
    '''
    Run humidity cross checks on precipitation obs

    :param object station: station object
    :param array flag_col: which columns to fill in flag_array
    :param datetime start: start of dataset
    :param datetime end: end of dataset
    :param file logfile: logfile to store outputs
    :param bool plots: do plots or not
    :param bool diagnostics: extra verbose information

    :returns:    
    '''

    precip1 = getattr(station, 'precip1_depth') 
    precip2 = getattr(station, 'precip2_depth') 
    precip3 = getattr(station, 'precip3_depth') 
    precip6 = getattr(station, 'precip6_depth') 
    precip9 = getattr(station, 'precip9_depth') 
    precip12 = getattr(station, 'precip12_depth') 
    precip15 = getattr(station, 'precip15_depth') 
    precip18 = getattr(station, 'precip18_depth') 
    precip24 = getattr(station, 'precip24_depth') 

    times = station.time.data

    # combine all the precips together
    precips = np.ma.array([precip1.data, precip2.data, precip3.data, precip6.data, precip9.data, precip12.data, precip15.data, precip18.data, precip24.data])

    station.qc_flags[:,flag_col[0]] = pcc_accumulations(times, precips, start, logfile, plots = plots, diagnostics = diagnostics) 
    
    station = utils.append_history(station, "Precipitation Cross Check")  

    return # pcc
    
#************************************************************************
if __name__ == "__main__":
    
    print "precipitation cross checks" 
    
    
    
    
