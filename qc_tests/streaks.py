#!/usr/local/sci/bin/python
#*****************************
#
# Repeated Streaks Check (RSC)
#
#   Checks for replication of 
#     1) checks for consecutive repeating values
#     2) checks if one year has more repeating strings than expected
#     3) checks for repeats at a given hour across a number of days
#     4) checks for repeats for whole days - all 24 hourly values
#
#
#   Some thresholds now determined dynamically
#
#************************************************************************
#                    SVN Info
#$Rev:: 114                                           $:  Revision of last commit
#$Author:: rdunn                                      $:  Author of last commit
#$Date:: 2017-01-17 17:26:42 +0000 (Tue, 17 Jan 2017) $:  Date of last commit
#************************************************************************
import numpy as np
import scipy as sp
import datetime as dt
import copy

# RJHD routines
import qc_utils as utils


# threshold values for low, mid and high resolution for each of the variables tested
# consecutive repeats (ignoring missing data); spread over N days; same at each hour of day for N days; full day repeats

T = {1: [40, 14, 25, 10], 0.5: [30, 10, 20, 7], 0.1: [24, 7, 15, 5]}
D = {1: [80, 14, 25, 10], 0.5: [60, 10, 20, 7], 0.1: [48, 7, 15, 2]}
S = {1: [120, 28, 25, 10], 0.5:[100, 21, 20, 7], 0.1:[72, 14, 15, 5]}
WS = {1: [40, 14, 25, 10], 0.5: [30, 10, 20, 7], 0.1: [24, 7, 15, 5]}
WD = {90: [120, 28, 28, 10], 45: [96, 28, 28, 10], 22: [72, 21, 21, 7], 10: [48, 14, 14, 7], 1: [24, 7, 14, 5]}

limits_dict = {"temperatures": T, "dewpoints":  D, "slp": S, "windspeeds": WS, "winddirs": WD}

WIND_MIN_VALUE = {1:0.5, 0.5:1.0, 0.1:0.5}

#*********************************************
def linear(X,p):
    '''
    decay function for line fitting
    p[0]=intercept
    p[1]=slope
    '''
    return p[1]*X + p[0] # linear

#*********************************************
def residuals_LS(p, Y, X):
    '''
    Least squared residuals from linear trend
    '''
    err = ((Y-linear(X,p))**2.0)

    return err # ResidualsLS

#************************************************************************
def rsc_get_straight_string_threshold(st_var, start, end, reporting = 0., diagnostics = False, plots = False, old_threshold = 0):
    '''
    Derive threshold number for strings/streaks of repeating values
    
    :param object st_var: station variable object
    :param datetime start: start of data
    :param datetime end: end of data    
    :param float reporting: reporting accuracy
    :param bool diagnostics: do diagnostic output
    :param bool plots: do plots
    :param float old_threshold: old threshold to use as comparison
    '''
    all_filtered = utils.apply_filter_flags(st_var)
    
   
    # find and count the length of all repeating strings
    
    prev_value = st_var.mdi
    this_string = []
    
    string_lengths =[]
    
    # run through all obs, the inefficient (non-pythonic) way
    for o, obs in enumerate(all_filtered):
        
        if all_filtered.mask[o] == False:
            
            if obs != prev_value:
                # if different value to before
                string_lengths += [len(this_string)]
                                       
                this_string = [o]
            else:
                # if same value as before, note and continue
                this_string += [o]
            prev_value = obs

    if plots:
        import calendar
        title = "Straight String Distribution"                  
        line_label = st_var.name
        xlabel = "String length"
    else:
        title, line_label, xlabel = "","",""
        
    threshold = utils.get_critical_values(string_lengths, binmin = 1, binwidth = 1, plots = plots, diagnostics = diagnostics, title = title, line_label = line_label, xlabel = xlabel, old_threshold = old_threshold)
 
    return threshold # rsc_get_straight_string_threshold


#************************************************************************
def rsc_diagnostics_and_plot(time, data, flags, title, start, plots = False):
    ''' plots time series of data with flagged streaks highlighted
    
    :param array time: time stamps in hours since
    :param array data: data to be plotted
    :param list flags: locations of obs to be flagged
    :param string title: title of plot (parameter)
    :param datetime start: dataset start date
    :param bool plots: do the plot
    '''


    YLABELS = {"temperatures":"Temperature (C)", "dewpoints":"Dewpoints (C)", "slp":"SLP (hPa)", "windspeeds":"Wind Speed (m/s)", "winddirs":"Degrees"}

    # get period to plot and convert times
    extra = 48
    min_t = flags[0] - extra
    max_t = flags[-1] + extra

    if min_t < 0: min_t = 0

    time = utils.times_hours_to_datetime(time[min_t:max_t], start)

    print "Streak at %s, %i observations" % (dt.datetime.strftime(time[extra], "%Y %m %d %H:%M"), len(flags))


    if plots:
        import matplotlib.pyplot as plt
        plt.clf()
        plt.plot(time, data[min_t:max_t], 'bo', ls = '-')
        
        flag_time = np.array(flags) - min_t
        plt.plot(time[flag_time], data[flags], 'ro', markersize = 10)
        plt.title(title.capitalize())
        plt.ylabel(YLABELS[title])
        
        plt.show()

    return # rsc_plots

#************************************************************************
def rsc_annual_string_expectance(all_filtered, value_starts, value_lengths, flags, start, end, st_var, times, diagnostics = False, plots = False):
    '''
    Find years where have more strings than expected, but not long enough to set off test
    
    :param array all_filtered: data filtered by all flags set so far
    :param array value_starts: locations of start of strings/streaks of data
    :param array value_lengths: lengths of each streak
    :param array flags: array of flags to be set
    :param datetime start: start of data
    :param datetime end: end of data    
    :param bool diagnostics: do diagnostic output
    :param bool plots: do plots
    '''

    month_starts = utils.month_starts(start,end)

    month_starts = np.array(month_starts).reshape(-1,12)

    year_proportions = np.zeros(month_starts.shape[0])
    year_proportions.fill(st_var.mdi)
    
    # churn through each year in turn
    for y in range(month_starts.shape[0]):
        
        if y != month_starts.shape[0] -1:
            year = all_filtered[month_starts[y,0] : month_starts[y+1,0]]
        else:
            year = all_filtered[month_starts[y,0] :]
        
        if len(year.compressed()) >= 200:
            # if there are strings (streaks of same value) in this year
            if y != month_starts.shape[0] -1:
                string_starts = np.where(np.logical_and((value_starts >= month_starts[y,0]),(value_starts < month_starts[y+1,0])))
            else:
                string_starts = np.where(value_starts >= month_starts[y,0])
            
            year_proportions[y] = 0
            if len(string_starts[0]) >= 1:
                # work out the proportion of the amount of data
                year_proportions[y] = np.sum(value_lengths[string_starts[0]])/float(len(year.compressed()))
            
    # if enough dirty years

    good_years = np.where(year_proportions != st_var.mdi)

    if len(good_years[0]) >= 10:
        
        median = np.median(year_proportions[good_years])
        
        if median < 0.005 : median = 0.005
        
        # find the number which have proportions > 5 x median
        bad_years = np.where(year_proportions > 5.*median)
        
        if len(bad_years[0]) >= 1:
            
            for bad in bad_years[0]:
                # and flag
                if bad == month_starts.shape[0]-1:
                    # if last year, just select all
                    locs, = np.where(value_starts >= month_starts[bad,0])
                else:
                    locs, = np.where((value_starts >= month_starts[bad,0]) & (value_starts <= month_starts[bad+1,0]))
                
                for loc in locs:
                    # need to account for missing values here  26/9/2014

                    goods, = np.where(all_filtered.mask[value_starts[loc]:] == False)

                    flags[value_starts[loc]+goods[:value_lengths[loc]]] = 1

                if plots or diagnostics:
                    plot_year = all_filtered[month_starts[bad,0]:month_starts[bad+1,0]]
                    plot_time = times[month_starts[bad,0]:month_starts[bad+1,0]]
                    plot_flags = np.where(flags[month_starts[bad,0]:month_starts[bad+1,0]] == 1)[0]

                    rsc_diagnostics_and_plot(plot_time, plot_year, plot_flags, st_var.name, start, plots = plots)           


    return flags # rsc_annual_string_expectance


#************************************************************************
def rsc_straight_strings(st_var, times, n_obs, n_days, start, end, wind = False, reporting = 0., diagnostics = False, plots = False, dynamic = True):
    '''
    Check for strings/streaks of repeating values
    
    :param object st_var: station variable object
    :param int n_days: number of days to exceed
    :param int n_obs: number of observations to exceed
    :param datetime start: start of data
    :param datetime end: end of data    
    :param float reporting: reporting accuracy
    :param bool wind: whether there is wind data to account for - extra minimum value
    :param bool diagnostics: do diagnostic output
    :param bool plots: do plots
    :param bool dynamic: calculate threshold of number of observations dynamically rather than using n_obs
    '''

    # January 2015 - changed to dynamically calculating the thresholds, but only use if less than current ^RJHD


    if st_var.name == "winddirs":
        # remove calm periods for this check.
        wd_st_var = copy.deepcopy(st_var)
        calms, = np.ma.where(st_var.data == 0) # True calms have direction set to 0, northerlies to 360
        wd_st_var.data[calms] = wd_st_var.mdi

        if dynamic:
            threshold = rsc_get_straight_string_threshold(wd_st_var, start, end, reporting = reporting, diagnostics = diagnostics, plots = plots, old_threshold = n_obs)          

            if threshold < n_obs: n_obs = threshold

        all_filtered = utils.apply_filter_flags(wd_st_var) # calms have been removed

    else:
        if dynamic:
            threshold = rsc_get_straight_string_threshold(st_var, start, end, reporting = reporting, diagnostics = diagnostics, plots = plots, old_threshold = n_obs)          
            
            if threshold < n_obs: n_obs = threshold

        all_filtered = utils.apply_filter_flags(st_var)
    
    flags = np.zeros(len(all_filtered))
    
    
    ''' Look for continuous straight strings '''
    
    prev_value = st_var.mdi
    string_points = []
    
    # storage for excess over years
    value_starts = []
    value_lengths =[]
    
    
    for o, obs in enumerate(all_filtered):
        
        if all_filtered.mask[o] == False:
            
            if obs != prev_value:
                if (st_var.name == "winddirs") and (prev_value == 0):
                    # this was a calm as a string of zeros.
                    # shouldn't be necessary - but just in case!
                    pass

                else:
                    # if different value to before, which is long enough (and large enough for Wind)
                    if len(string_points) >= 10:
                        if wind == False or (wind == True and prev_value > WIND_MIN_VALUE[reporting]):
                            # note start and length for the annual excess test
                            value_starts += [string_points[0]]
                            value_lengths += [len(string_points)]

                            time_diff = times[string_points[-1]] - times[string_points[0]]

                            # if length above threshold and spread over sufficient time frame, flag
                            if (len(string_points) >= n_obs) or (time_diff >= (n_days * 24)): # measuring time in hours 
                                flags[string_points] = 1
                                if plots or diagnostics:
                                    rsc_diagnostics_and_plot(times, all_filtered, string_points, st_var.name, start, plots = plots)           
                                
                string_points = [o]
            else:
                # if same value as before, note and continue
                string_points += [o]
            prev_value = obs

    # matches value_lengths 030660-99999, 1/7/2014 - seems to flag more though - compressed vs full time?

    flags = rsc_annual_string_expectance(all_filtered, np.array(value_starts), np.array(value_lengths), flags, start, end, st_var, times, diagnostics = diagnostics, plots = plots)
    
    
    return flags # rsc_straight_strings

#************************************************************************
def rsc_whole_day_repeats(data, n_wday, st_var, diagnostics = False, plots = False):
    '''
    Look through for repeats of full 24h of data
    
    :param array data: data to process
    :param int n_wday: number of whole days to exceed
    :param MetVar st_var: station variable - mainly for plotting
    :param bool diagnostics: do diagnostic output
    :param bool plots: do plots
    '''
    
    flags = np.zeros(data.shape)
    
    ndays = [] # to get distribution

    for day in range(data.shape[1]):
        if day == 0:
            prev_day = data[day, :]
            nday = 1
            continue
        else:
        
            matches = np.where(prev_day == data[day, :])

            # if this day matches previous one (not IDL wording, but matches logic)
            if len(matches[0]) == 24:
                if nday >= n_wday:

                    if len(data[day, :].compressed()) > 0:
                        # if above threshold and not because all values are missing, flag
                        flags[day-nday:day, :] = 1

                        if plots:
                            rsc_diagnostics_and_plot(st_var.time.data, st_var.data, np.arange((day-nday)*24, day*24), st_var.name, start, plots = plots)           
                            
            else:
                ndays += [nday]
                prev_day = data[day, :]
                nday = 1

            nday += 1
 
    return flags.reshape(-1) # rsc_whole_day_repeats


#************************************************************************
def rsc_hourly_repeats(st_var, times, n_hrs, n_wday, diagnostics = False, plots = False):
    '''
    Repeat of same value at given hour for >N days 
    
    :param object st_var: station variable object
    :param array times: timestamps
    :param int n_hrs: number of hours to exceed
    :param int n_wday: number of whole days to exceed (passed on)
    :param bool diagnostics: do diagnostic output
    :param bool plots: do plots
    '''


    flags = np.zeros(len(st_var.data))
    hourly_data = utils.apply_filter_flags(st_var)
    
    hourly_data = hourly_data.reshape(-1,24)
    hourly_times = times.reshape(hourly_data.shape)
    
    for hour in range(24):
        
        match_values = -999.
        match_times = [] # assumes start at time zero

        len_matches = [] # for distribution

        for day in range(hourly_data.shape[0]):
            # for each day at each given hour
            if hourly_data.mask[day, hour] == False:
                if hourly_data[day, hour] != match_values:
                    # if different value, check if string/streak above threshold
                    if len(match_times) > n_hrs:
                        
                        bad = np.where(match_times == times)
                        flags[bad] = 1

                        if plots:
                            rsc_diagnostics_and_plot(st_var.time.data, st_var.data, bad, st_var.name, start, plots = plots)           
                            
                    
                    len_matches += [len(match_times)]
                    match_values = hourly_data[day, hour]
                    match_times = [hourly_times[day, hour]]
                else:
                    # if same value
                    match_times +=[hourly_times[day, hour]]
                    
                    
    day_flags = rsc_whole_day_repeats(hourly_data, n_wday, st_var, diagnostics = diagnostics, plots = plots)
 
    return flags, day_flags # rsc_hourly_repeats


#************************************************************************
def rsc(station, var_list, flag_col, start, end, logfile, diagnostics = False, plots = False):
    ''' Wrapper for the four individual repeating streak check tests '''
    
    times = station.time.data
    
    for v, variable in enumerate(var_list):
        
        st_var = getattr(station, variable)
                
        if len(utils.apply_filter_flags(st_var).compressed()) > 0:
        
            wind = False
            if variable == "windspeeds": wind = True
            winddir= False
            if variable == "winddirs": winddir = True

            reporting_resolution = utils.reporting_accuracy(utils.apply_filter_flags(st_var), winddir = winddir, plots = plots)

            limits = limits_dict[variable][reporting_resolution]  

            # need to apply flags to st_var.flags each time for filtering
            station.qc_flags[:,flag_col[v][0]] = rsc_straight_strings(st_var, times, limits[0], limits[1], start, end, reporting = reporting_resolution, wind = wind, diagnostics = diagnostics, plots = plots, dynamic = True)

            station.qc_flags[:, flag_col[v][1]], station.qc_flags[:, flag_col[v][2]]= rsc_hourly_repeats(st_var, times, limits[2], limits[3], diagnostics = diagnostics, plots = plots)


            for streak_type in range(3):
                flag_locs = np.where(station.qc_flags[:, flag_col[v][streak_type]] != 0)
                if plots or diagnostics:
                    utils.print_flagged_obs_number(logfile, "Streak Check", variable, len(flag_locs[0]), noWrite = True)
                else:
                    utils.print_flagged_obs_number(logfile, "Streak Check", variable, len(flag_locs[0]))


                # copy flags into attribute
                st_var.flags[flag_locs] = 1

    station = utils.append_history(station, "Streak Check")
    
    return


#************************************************************************
if __name__ == "__main__":
    
    print "checking repeated strings"
