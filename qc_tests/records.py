#!/usr/local/sci/bin/python
#*****************************
#
# Known Records Check (KRC)
#
#   Check for exceedence of world records
#
#
#************************************************************************
#                    SVN Info
#$Rev:: 219                                           $:  Revision of last commit
#$Author:: rdunn                                      $:  Author of last commit
#$Date:: 2019-05-20 16:56:47 +0100 (Mon, 20 May 2019) $:  Date of last commit
#************************************************************************
import numpy as np
import scipy as sp
import datetime as dt
# RJHD routines
import qc_utils as utils

# updated max dewpoint after David's paper.
T_X = {"Africa":55.0,"Asia":53.9,"South_America":48.9,"North_America":56.7,"Europe":48.0,"Pacific":50.7,"Antarctica":15.0,"ROW":56.7}
T_N = {"Africa":-23.9,"Asia":-67.8,"South_America":-32.8,"North_America":-63.0,"Europe":-58.1,"Pacific":-23.0,"Antarctica":-89.2,"ROW":-89.2}
D_X = {"Africa":55.0,"Asia":53.9,"South_America":48.9,"North_America":56.7,"Europe":48.0,"Pacific":50.7,"Antarctica":15.0,"ROW":56.7}
D_N = {"Africa":-50.,"Asia":-100.,"South_America":-60.,"North_America":-100.,"Europe":-100.,"Pacific":-50.,"Antarctica":-100.,"ROW":-100.}
W_X = {"Africa":113.2,"Asia":113.2,"South_America":113.2,"North_America":113.2,"Europe":113.2,"Pacific":113.2,"Antarctica":113.2,"ROW":113.2}
W_N = {"Africa":0.,"Asia":0.,"South_America":0.,"North_America":0.,"Europe":0.,"Pacific":0.,"Antarctica":0.,"ROW":0.}
S_X = {"Africa":1083.3,"Asia":1083.3,"South_America":1083.3,"North_America":1083.3,"Europe":1083.3,"Pacific":1083.3,"Antarctica":1083.3,"ROW":1083.3}
S_N = {"Africa":870.,"Asia":870.,"South_America":870.,"North_America":870.,"Europe":870.,"Pacific":870.,"Antarctica":870.,"ROW":870.}


maxes = {"temperatures": T_X, "dewpoints": D_X, "windspeeds": W_X, "slp": S_X}
mins = {"temperatures": T_N, "dewpoints": D_N, "windspeeds": W_N, "slp": S_N}


#************************************************************************
def krc_get_wmo_region(stnid):
    ''' Get the WMO region from the station id '''
    
    region = "ROW" # rest of world
    if 600000 <= stnid <= 699999:
        region = "Africa"
    if 200000 <= stnid <= 200999:
        region = "Asia"
    if 202000 <= stnid <= 219999:
        region = "Asia"
    if 230000 <= stnid <= 259999:
        region = "Asia"
    if 280000 <= stnid <= 329999:
        region = "Asia"
    if 350000 <= stnid <= 369999:
        region = "Asia"
    if 380000 <= stnid <= 399999:
        region = "Asia"
    if 403500 <= stnid <= 485999:
        region = "Asia"
    if 488000 <= stnid <= 499999:
        region = "Asia"
    if 500000 <= stnid <= 599999:
        region = "Asia"
    if 800000 <= stnid <= 889999:
        region = "South_America"
    if 700000 <= stnid <= 799999:
        region = "North_America"
    if 486000 <= stnid <= 487999:
        region = "Pacific"
    if 900000 <= stnid <= 989999:
        region = "Pacific"
    if stnid <= 199999:
        region = "Europe"
    if 201000 <= stnid <= 201999:
        region = "Europe" 
    if 220000 <= stnid <= 229999:
        region = "Europe"
    if 260000 <= stnid <= 279999:
        region = "Europe"
    if 330000 <= stnid <= 349999:
        region = "Europe"
    if 370000 <= stnid <= 379999:
        region = "Europe"
    if 400000 <= stnid <= 403499:
        region = "Europe"
    if 890000 <= stnid <= 899999:
        region = "Antarctica" 

    return region # krc_get_wmo_region


#************************************************************************
def krc_set_flags(locs, flags, col):
    '''
    Set the flags to 1 for correct column if data exists
    
    :param array locs: locations for flags
    :param array flags: the flags
    :param int col: column to use
    '''
    
    if len(locs[0]) > 0:
                
        flags[locs, col] = 1  
    
    return # krc_set_flags

#************************************************************************
def krc(station, var_list, flag_col, logfile, diagnostics = False, plots = False):
    '''
    Run the known records check for each variable in list
    
    :param object station: station to process
    :param list var_list: list of variables to process
    :param list flag_col: which columns to use for which variable 
    :param file logfile: logfile to store output 
    :param bool diagnostics: diagnostic output (unused)
    :param bool plots: do the plots (unused)
    
    '''
    
    
    for v, variable in enumerate(var_list):
        
        st_var = getattr(station, variable)
    
        st_region = krc_get_wmo_region(station.id)
        
        all_filtered = utils.apply_filter_flags(st_var)
        
        too_high = np.where(all_filtered > maxes[variable][st_region])
        
        krc_set_flags(too_high, station.qc_flags, flag_col[v])
        
        # make sure that don't flag the missing values!
        too_low = np.where(np.logical_and(all_filtered < mins[variable][st_region], all_filtered.mask == False ))
        
        krc_set_flags(too_low, station.qc_flags, flag_col[v])
        
        flag_locs = np.where(station.qc_flags[:, flag_col[v]] != 0)

        utils.print_flagged_obs_number(logfile, "World Record", variable, len(flag_locs[0]), noWrite = diagnostics)

        # copy flags into attribute
        st_var.flags[flag_locs] = 1
	    
    station = utils.append_history(station, "World Record Check")  

    return # krc

#************************************************************************
if __name__ == "__main__":
    
    print "checking for exceedence of world records"
