#!/usr/local/sci/bin/python
#*****************************
#
# Cloud Coverage Logical Check (CCC)
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



#************************************************************************
def unobservable(station, flag_col, logfile, diagnostics = False, plots = False):
    '''
    Cloud observation code given as unobservable (==9 or 10)

    :param obj station: station object
    :param list flag_col: flag columns to use
    :param file logfile: logfile to store outpu
    :param bool plots: to do any plots
    :param bool diagnostics: to do any extra diagnostic output

    :returns:
    '''

    # for each cloud variable, find bad locations and flag
    for c, cloud in enumerate(['total_cloud_cover','low_cloud_cover','mid_cloud_cover','high_cloud_cover']):

        cloud_obs = getattr(station, cloud)

        bad_locs = np.ma.where(np.logical_or(cloud_obs.data == 9, cloud_obs.data == 10))
        
        station.qc_flags[bad_locs, flag_col[c]] = 1

        flag_locs = np.where(station.qc_flags[:, flag_col[c]] != 0)
        utils.print_flagged_obs_number(logfile, "Unobservable cloud", cloud, len(flag_locs[0]), noWrite = diagnostics)

        # copy flags into attribute
        cloud_obs.flags[flag_locs] = 1

    return # unobservable


#************************************************************************
def total_lt_max(station, flag_col, logfile, diagnostics = False, plots = False):
    '''
    Total cloud cover less than maximum of low, mid and high
    
    :param obj station: station object
    :param list flag_col: flag columns to use
    :param file logfile: logfile to store outpu
    :param bool plots: to do any plots
    :param bool diagnostics: to do any extra diagnostic output

    :returns:
    '''
    
    total = getattr(station, "total_cloud_cover")
    low   = getattr(station, "low_cloud_cover")
    mid   = getattr(station, "mid_cloud_cover")
    high  = getattr(station, "high_cloud_cover")
    
    maximum = np.ma.max([low.data, mid.data, high.data], axis = 0)
    
    bad_locs = np.ma.where(maximum > total.data)
    
    station.qc_flags[bad_locs, flag_col] = 1

    flag_locs = np.where(station.qc_flags[:, flag_col] != 0)
    utils.print_flagged_obs_number(logfile, "Total < Max cloud", "cloud", len(flag_locs[0]), noWrite = diagnostics)

    # copy flags into attribute
    total.flags[flag_locs] = 1
    low.flags[flag_locs] = 1
    mid.flags[flag_locs] = 1
    high.flags[flag_locs] = 1

    return # total_lt_max

#************************************************************************
def low_full(station, flag_col, logfile, diagnostics = False, plots = False):
    '''
    Low cloud full, but values in mid or high
    
    :param obj station: station object
    :param list flag_col: flag columns to use
    :param file logfile: logfile to store outpu
    :param bool plots: to do any plots
    :param bool diagnostics: to do any extra diagnostic output

    :returns:
    '''
    
    low   = getattr(station, "low_cloud_cover")
    mid   = getattr(station, "mid_cloud_cover")
    high  = getattr(station, "high_cloud_cover")
   
    low_full_locs = np.ma.where(low.data == 8)
    
    bad_mid = np.where(mid.data.mask[low_full_locs] != True)
    
    station.qc_flags[low_full_locs[0][bad_mid[0]], flag_col] = 1    
    
    bad_high = np.where(high.data.mask[low_full_locs] != True)
 
    station.qc_flags[low_full_locs[0][bad_high[0]], flag_col] = 1
    
    flag_locs = np.where(station.qc_flags[:, flag_col] != 0)
    utils.print_flagged_obs_number(logfile, "Low full cloud", "cloud", len(flag_locs[0]), noWrite = diagnostics)

    # copy flags into attribute
    mid.flags[flag_locs] = 1
    high.flags[flag_locs] = 1

    return # low_full

 
#************************************************************************
def mid_full(station, flag_col, logfile, diagnostics = False, plots = False):
    '''
    Mid cloud full, but values in high
    
    :param obj station: station object
    :param list flag_col: flag columns to use
    :param file logfile: logfile to store outpu
    :param bool plots: to do any plots
    :param bool diagnostics: to do any extra diagnostic output

    :returns:
    '''
    
    mid   = getattr(station, "mid_cloud_cover")
    high  = getattr(station, "high_cloud_cover")
   
    mid_full_locs = np.ma.where(mid.data == 8)
       
    bad_high = np.where(high.data.mask[mid_full_locs] != True)
 
    station.qc_flags[mid_full_locs[0][bad_high[0]], flag_col] = 1
    
    flag_locs = np.where(station.qc_flags[:, flag_col] != 0)
    utils.print_flagged_obs_number(logfile, "Mid full cloud", "cloud", len(flag_locs[0]), noWrite = diagnostics)

    # copy flags into attribute
    high.flags[flag_locs] = 1

    return # mid_full

#************************************************************************
def fix_cloud_base(station):
    '''
    If cloud base is 22000ft, then set to missing
    
    :param obj station: station object

    :returns:
    '''

    cloud_base = getattr(station, "cloud_base")
 
    bad_cb = np.where(cloud_base.data == 22000)
 
    # no flag set on purpose - just set to missing (unobservable)
    cloud_base.data[bad_cb] = cloud_base.mdi
    cloud_base.data.mask[bad_cb] = True
    
    return

#************************************************************************
def negative_cloud(station, flag_col, logfile, diagnostics = False, plots = False):
    '''
    Non-sensical cloud value
    
    :param obj station: station object
    :param list flag_col: flag columns to use
    :param file logfile: logfile to store outpu
    :param bool plots: to do any plots
    :param bool diagnostics: to do any extra diagnostic output

    :returns:
    '''

    # go through each cloud varaible and flag bad locations
    for c, cloud in enumerate(['total_cloud_cover','low_cloud_cover','mid_cloud_cover','high_cloud_cover']):

        cloud_obs = getattr(station, cloud)

        bad_locs = np.ma.where(cloud_obs.data < 0)

        station.qc_flags[bad_locs, flag_col] = 1

        # copy flags into attribute
        cloud_obs.flags[bad_locs] = 1


    flag_locs = np.where(station.qc_flags[:, flag_col] != 0)
    utils.print_flagged_obs_number(logfile, "Negative Cloud", "cloud", len(flag_locs[0]), noWrite = diagnostics)

    return

#************************************************************************
def ccc(station, flag_col, logfile, diagnostics = False, plots = False):
    '''
    Call the logical cloud checks
    
    :param obj station: station object
    :param list flag_col: flag columns to use
    :param file logfile: logfile to store output
    :param bool diagnostics: diagnostic output (unused)
    :param bool plots: do the plots (unused)

    :returns:
    '''

    if len(flag_col) != 8:
        print "insufficient flag columns given"
        return

    unobservable(station, flag_col[0:4], logfile, plots = plots, diagnostics = diagnostics)
    
    total_lt_max(station, flag_col[4], logfile, plots = plots, diagnostics = diagnostics) 
    
    low_full(station, flag_col[5], logfile, plots = plots, diagnostics = diagnostics) 
    
    mid_full(station, flag_col[6], logfile, plots = plots, diagnostics = diagnostics) 
    
    fix_cloud_base(station)
    
    negative_cloud(station, flag_col[7], logfile, plots = plots, diagnostics = diagnostics)
    
    station = utils.append_history(station, "Cloud - Logical Cross Check")  
                     
    return # ccc

#************************************************************************
if __name__ == "__main__":
    
    print "cloud level cross checks" 
    

  
