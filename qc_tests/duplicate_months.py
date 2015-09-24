#!/usr/local/sci/bin/python
#*****************************
#
# Duplicate Months Check (DMC)
#
#
#************************************************************************
#                    SVN Info
#$Rev:: 55                                            $:  Revision of last commit
#$Author:: rdunn                                      $:  Author of last commit
#$Date:: 2015-02-06 16:38:46 +0000 (Fri, 06 Feb 2015) $:  Date of last commit
#************************************************************************

import numpy as np
import scipy as sp
import datetime as dt

# RJHD routines
import qc_utils as utils

#************************************************************************
def dmc_plot(source_data, target_data,start, source_offset, target_offset, ylabel):
    '''
    Plot the two timeseries which are duplicated
    :param array source_data: input source month (masked array)
    :param array target_data: input comparison month (masked array)
    :param datetime start: start of dataseries
    :param int source_offset: # hours after start for beginning of source data
    :param int target_offset: # hours after start for beginning of target data
    :param str ylabel: y-label of plot    
    '''
    
    import matplotlib.pyplot as plt
    plot_date = start + dt.timedelta(hours = source_offset)                            
    plot_times = []                           
    while len(plot_times) < len(source_data):
        plot_times += [plot_date]
        plot_date = plot_date + dt.timedelta(hours = 1)      

    # plot overlapping data
    plt.clf()
    plt.plot(np.ma.array(plot_times, mask = source_data.mask).compressed(), source_data.compressed(), \
                 'bo', ls='-', label = "source %s" % (dt.datetime.strftime(plot_times[0], "%Y-%m")))

    target_date = start + dt.timedelta(hours = target_offset)
    plt.plot(np.ma.array(plot_times, mask = target_data.mask).compressed(), target_data.compressed(), \
                 'ro', ls='-', label = "target %s" % (dt.datetime.strftime(target_date, "%Y-%m")))

    plt.legend(loc=8)
    plt.ylabel(ylabel)
    plt.show()

    return # dmc_plot

#************************************************************************
def is_month_duplicated(source_data, target_data, valid, sm, tm, duplicated):
    '''
    Perform test whether month is duplicated
    :param array source_data: source data
    :param array target_data: target data
    :param array valid: valid values
    :param int sm: source month counter
    :param int tm: target month counter
    :param array duplicated: True/False array for duplicated months
    '''
    
    match = np.where(source_data.compressed()[valid] == target_data.compressed()[valid])
                            
    if len(match[0]) == len(valid[0]):
        duplicated[sm] = 1
        duplicated[sm + 1 + tm] = 1
                            
    return duplicated # test_if_duplicated


 #************************************************************************
def duplication_test(source_data, target_data, valid, sm, tm, source_month, target_month, duplicated, diagnostics, flags, flag_col):
    '''
    Pass into test whether month is duplicated and set flags
    :param array source_data: source data
    :param array target_data: target data
    :param array valid: valid values
    :param int sm: source month counter
    :param int tm: target month counter
    :param array duplicated: True/False array for duplicated months
    '''
    
    
    duplicated = is_month_duplicated(source_data, target_data, valid, sm, tm, duplicated)
                                                                    
    if duplicated[sm] == 1:
        # make flags if this month is duplicated somewhere
        flags[source_month[0]:source_month[1], flag_col] = 1
        flags[target_month[0]:target_month[1], flag_col] = 1
        if diagnostics: print "   Flagging"
                            
        
    return duplicated # duplication_test


#************************************************************************
def dmc(station, variable_list, full_variable_list, flag_col, start, end, logfile, diagnostics=False, plots=False):
    '''
    Method copied from check_duplicates.pro
    :param obj station: station object with suitable attributes (see netcdf_procs.py)
    :param list variable_list: list of netcdf variables to process
    :param list full_variable_list: the variables for flags to be applied to
    :param list flag_col: which column to set in flag array
    :param datetime start: data start
    :param datetime end: data end
    :param file logfile: logfile to store outputs
    :param bool diagnostics: extra verbosity
    :param bool plots: do plots
    '''
    MIN_DATA_REQUIRED = 20 # obs per month
    
    # get array of Nx2 start/end pairs
    month_ranges = utils.month_starts_in_pairs(start, end)
    
    for v, variable in enumerate(variable_list):
    
        st_var = getattr(station, variable)
        
        # double loop structure - not ideal
        duplicated = np.zeros(len(month_ranges))
        
        for sm, source_month in enumerate(month_ranges):
            if diagnostics: print "Month %i of %i" % (sm + 1, len(month_ranges))
    
            source_data = st_var.data[source_month[0]:source_month[1]]
    
            if duplicated[sm] == 0:
                # don't repeat if already a duplicated
                for tm, target_month in enumerate(month_ranges[sm + 1:]):
                    target_data = st_var.data[target_month[0]:target_month[1]]
               
                    # match the data periods
                    overlap = np.min([len(source_data), len(target_data)])
                    s_data, t_data = source_data[:overlap], target_data[:overlap]
                    s_valid, t_valid = np.where(s_data.compressed() != st_var.fdi), \
                        np.where(t_data.compressed() != st_var.fdi)
                
                    # if enough of an overlap
                    if (len(s_valid[0]) >= MIN_DATA_REQUIRED) and \
                        (len(t_valid[0]) >= MIN_DATA_REQUIRED):                        
                        
                        if len(s_valid[0]) < len(t_valid[0]):   
                            
                            duplicated = duplication_test(source_data, target_data, s_valid, sm, tm, source_month, target_month, duplicated, diagnostics, station.qc_flags, flag_col[v])
                                                           
                        else:     
                            # swap the list of valid points                       
                            duplicated = duplication_test(source_data, target_data, t_valid, sm, tm, source_month, target_month, duplicated, diagnostics, station.qc_flags, flag_col[v])

                        if plots:
                            dmc_plot(s_data, t_data, start, source_month[0], target_month[0], st_var.name)
        
                    # target month
            # source month
        # variable list
        flag_locs = np.where(station.qc_flags[:, flag_col[v]] != 0)
        if plots or diagnostics:
            utils.print_flagged_obs_number(logfile, "Duplicate Month", variable, len(flag_locs[0]), noWrite = True)
        else:
            utils.print_flagged_obs_number(logfile, "Duplicate Month", variable, len(flag_locs[0]))

        # copy flags into attribute
        st_var.flags[flag_locs] = 1
	
    utils.apply_flags_all_variables(station, full_variable_list, flag_col[variable_list == "temperatures"], logfile, "Duplicate Months", plots = plots, diagnostics = diagnostics)

    station = utils.append_history(station, "Duplicate Months Check")  
                     
    return # dmc

#************************************************************************
if __name__ == "__main__":

    print "Checking for Duplicate Months"

#************************************************************************
