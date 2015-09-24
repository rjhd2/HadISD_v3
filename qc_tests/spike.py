#!/usr/local/sci/bin/python
#*****************************
#
# Spike Check (SC)
#
#
#************************************************************************
#                    SVN Info
#$Rev:: 55                                            $:  Revision of last commit
#$Author:: rdunn                                      $:  Author of last commit
#$Date:: 2015-02-06 16:38:46 +0000 (Fri, 06 Feb 2015) $:  Date of last commit
#************************************************************************
import numpy as np
import math
import scipy as sp
import datetime as dt

# RJHD routines
import qc_utils as utils


#************************************************************************
def sc_diagnostics_and_plots(times, indata, start, end, datastart, title, plots = True):
    '''
    Plot each set of values highlighted by Spike Check

    :param array times: array of times (hours since)
    :param array indata: data to plot
    :param int start: start of flagging period
    :param int end : end of flagging period
    :param datetime datastart: start of dataset
    :param string title: title of plot
    :param bool plots: do the plot

    :returns:
    '''

    YLABELS = {"temperatures":"Temperature (C)", "dewpoints":"Dewpoints (C)", "slp":"SLP (hPa)", "windspeeds":"Wind Speed (m/s)"}

    extra = 48

    plot_times = utils.times_hours_to_datetime(times[start-extra:end+extra], datastart)

    print "Spike at %s" % dt.datetime.strftime(plot_times[extra], "%Y %m %d %H:%M")

    if plots:

        import matplotlib.pyplot as plt
        plt.clf()
        plt.plot(plot_times, indata[start-extra:end+extra], 'bo', ls='-')
        plt.plot(plot_times[extra:-extra], indata[start:end], 'ro', markersize=10)
        plt.title("Spike "+title.capitalize())
        plt.ylabel(YLABELS[title])
        plt.show()

    return # sc_plots


#************************************************************************
def sc(station, variable_list, flag_col, start, end, logfile, diagnostics = False, plots = False, second = False):
    '''
    Spike Check, looks for spikes up to 3 observations long, using thresholds
    calculated from the data itself.

    :param MetVar station: the station object
    :param list variable_list: list of observational variables to process
    :param list flag_col: the columns to set on the QC flag array
    :param datetime start: dataset start time
    :param datetime end: dataset end time
    :param file logfile: logfile to store outputs
    :param bool plots: do plots
    :param bool second: run for second time

    :returns:    
    '''
    print "refactor"
    
    for v, variable in enumerate(variable_list):

        flags = station.qc_flags[:, flag_col[v]]

        prev_flag_number = 0
        if second:
            # count currently existing flags:
            prev_flag_number = len(flags[flags != 0])
    
        st_var = getattr(station, variable)
    
        all_filtered = utils.apply_filter_flags(st_var)
      
        reporting_resolution = utils.reporting_accuracy(utils.apply_filter_flags(st_var))
        # to match IDL system - should never be called as would mean no data
        if reporting_resolution == -1: reporting_resolution = 1 

        month_ranges = utils.month_starts_in_pairs(start, end)
        month_ranges = month_ranges.reshape(-1,12,2)
        
        good = np.where(all_filtered.mask == False)
        
        full_time_diffs = np.ma.zeros(len(all_filtered))
        full_time_diffs.mask = all_filtered.mask
        full_time_diffs[good] = station.time.data[good][1:] - station.time.data[good][:-1]
        
        # develop critical values using clean values
        # NOTE 4/7/14 - make sure that Missing and Flagged values treated appropriately
        print "sort the differencing if values were flagged rather than missing"

        full_filtered_diffs = np.ma.zeros(len(all_filtered))
        full_filtered_diffs.mask = all_filtered.mask
        full_filtered_diffs[good] = all_filtered.compressed()[1:] - all_filtered.compressed()[:-1]
        
        # test all values
        good_to_uncompress = np.where(st_var.data.mask == False)
        full_value_diffs = np.ma.zeros(len(st_var.data))
        full_value_diffs.mask = st_var.data.mask
        full_value_diffs[good_to_uncompress] = st_var.data.compressed()[1:] - st_var.data.compressed()[:-1]

        # convert to compressed time to match IDL
        value_diffs = full_value_diffs.compressed()
        time_diffs = full_time_diffs.compressed()
        filtered_diffs = full_filtered_diffs.compressed()
        flags = flags[good_to_uncompress]
        

        critical_values = np.zeros([9,12])
        critical_values.fill(st_var.mdi)
        
        # link observation to calendar month
        month_locs = np.zeros(full_time_diffs.shape)
                
        for month in range(12):
            for year in range(month_ranges.shape[0]):
                
                if year == 0:
                    this_month_time_diff = full_time_diffs[month_ranges[year,month,0]:month_ranges[year,month,1]]
                    this_month_filtered_diff = full_filtered_diffs[month_ranges[year,month,0]:month_ranges[year,month,1]]
                else:
                    this_month_time_diff = np.ma.concatenate([this_month_time_diff, full_time_diffs[month_ranges[year,month,0]:month_ranges[year,month,1]]])
                    this_month_filtered_diff = np.ma.concatenate([this_month_filtered_diff, full_filtered_diffs[month_ranges[year,month,0]:month_ranges[year,month,1]]])


                month_locs[month_ranges[year,month,0]:month_ranges[year,month,1]] = month
      
            for delta in range(1,9):
                
                locs = np.ma.where(this_month_time_diff == delta)
        
                if len(locs[0]) >= 100:
                    
                    iqr = utils.IQR(this_month_filtered_diff[locs])

                    if iqr == 0. and delta == 1:
                        critical_values[delta-1,month] = 6.
                    elif iqr == 0: 
                        critical_values[delta-1,month] = st_var.mdi
                    else:
                        critical_values[delta-1,month] = 6. * iqr      

                    # January 2015 - changed to dynamically calculating the thresholds if less than IQR method ^RJHD

                    if plots:
                        import calendar
                        title = "{}, {}-hr differences".format(calendar.month_name[month+1], delta)                  
                        line_label = st_var.name
                        xlabel = "First Difference Magnitudes"
                    else:
                        title, line_label, xlabel = "","",""

                    threshold = utils.get_critical_values(this_month_filtered_diff[locs], binmin = 0, binwidth = 0.5, plots = plots, diagnostics = diagnostics, title = title, line_label = line_label, xlabel = xlabel, old_threshold = critical_values[delta-1,month])

                    if threshold < critical_values[delta-1,month]: critical_values[delta-1,month] = threshold

                    if plots or diagnostics:

                        print critical_values[delta-1,month] , iqr, 6 * iqr
           

        month_locs = month_locs[good_to_uncompress]
        if diagnostics:
            print critical_values[0,:]
                
        # not less than 5x reporting accuracy
        good_critical_values = np.where(critical_values != st_var.mdi)
        low_critical_values = np.where(critical_values[good_critical_values] <= 5.*reporting_resolution)
        temporary = critical_values[good_critical_values]
        temporary[low_critical_values] = 5.*reporting_resolution
        critical_values[good_critical_values] = temporary
        
        
        if diagnostics:
            print critical_values[0,:], 5.*reporting_resolution

        # check hourly against 2 hourly, if <2/3 the increase to avoid crazy rejection rate
        for month in range(12):
            if critical_values[0,month] != st_var.mdi and critical_values[1,month] != st_var.mdi:
                if critical_values[0,month]/critical_values[1,month] <= 0.66:
                    critical_values[0,month] = 0.66 * critical_values[1,month]
        
        if diagnostics:
            print critical_values[0,:]


        # get time differences for unfiltered data

        full_time_diffs = np.ma.zeros(len(st_var.data))
        full_time_diffs.mask = st_var.data.mask
        full_time_diffs[good_to_uncompress] = station.time.data[good_to_uncompress][1:] - station.time.data[good_to_uncompress][:-1]
        time_diffs = full_time_diffs.compressed()

        # go through each difference, identify which month it is in if passes spike thresholds 
    
        # spikes at the beginning or ends of sections
        for t in np.arange(len(time_diffs)):
            if (np.abs(time_diffs[t - 1]) > 240) and (np.abs(time_diffs[t]) < 3):
                # 10 days before but short gap thereafter
                
                next_values = st_var.data[good_to_uncompress[0][t + 1:]] 
                good, = np.where(next_values.mask == False)
        
                next_median = np.ma.median(next_values[good[:10]])
        
                next_diff = np.abs(value_diffs[t]) # out of spike
                median_diff = np.abs(next_median - st_var.data[good_to_uncompress[0][t]]) # are the remaining onees
                       
                if (critical_values[time_diffs[t] - 1, month_locs[t]] != st_var.mdi):
                    
                    # jump from spike > critical but average after < critical / 2
                    if (np.abs(median_diff) < critical_values[time_diffs[t] - 1, month_locs[t]] / 2.) and\
                        (np.abs(next_diff) > critical_values[time_diffs[t] - 1, month_locs[t]]) :
                    
                        flags[t] = 1
                        if plots or diagnostics:
                            sc_diagnostics_and_plots(station.time.data, st_var.data, good_to_uncompress[0][t], good_to_uncompress[0][t+1], start, variable, plots = plots)
                        
                        
            elif (np.abs(time_diffs[t - 1]) < 3) and (np.abs(time_diffs[t]) > 240):
                # 10 days after but short gap before
                
                prev_values = st_var.data[good_to_uncompress[0][:t - 1]]
                good, = np.where(prev_values.mask == False)
        
                prev_median = np.ma.median(prev_values[good[-10:]])
        
                prev_diff = np.abs(value_diffs[t - 1])
                median_diff = np.abs(prev_median - st_var.data[good_to_uncompress[0][t]])
        
                if (critical_values[time_diffs[t - 1] - 1, month_locs[t]] != st_var.mdi):
                    
                    # jump into spike > critical but average before < critical / 2
                    if (np.abs(median_diff) < critical_values[time_diffs[t - 1] - 1, month_locs[t]] / 2.) and\
                        (np.abs(prev_diff) > critical_values[time_diffs[t - 1] - 1, month_locs[t]]) :
                    
                        flags[t] = 1
                        if plots or diagnostics:
                            sc_diagnostics_and_plots(station.time.data, st_var.data, good_to_uncompress[0][t], good_to_uncompress[0][t+1], start, variable, plots = plots)
                        
        
        
        
        ''' this isn't the nicest way, but a direct copy from IDL
            masked arrays might help remove some of the lines

            Also, this is relatively slow'''
            
        for t in np.arange(len(time_diffs)):
            for spk_len in [1,2,3]:
                if t >= spk_len and t < len(time_diffs) - spk_len:
                    
                    # check if time differences are appropriate, for multi-point spikes
                    if (np.abs(time_diffs[t - spk_len]) <= spk_len * 3) and\
                    (np.abs(time_diffs[t]) <= spk_len * 3) and\
                    (time_diffs[t - spk_len - 1] - 1 < spk_len * 3) and\
                    (time_diffs[t + 1] - 1 < spk_len * 3) and \
                    ((spk_len == 1) or \
                    ((spk_len == 2) and (np.abs(time_diffs[t - spk_len + 1]) <= spk_len * 3)) or \
                    ((spk_len == 3) and (np.abs(time_diffs[t - spk_len + 1]) <= spk_len * 3) and (np.abs(time_diffs[t - spk_len + 2]) <= spk_len * 3))):
                        
                        # check if differences are valid                        
                        if (value_diffs[t - spk_len] != st_var.mdi) and \
                        (value_diffs[t - spk_len] != st_var.fdi) and \
                        (critical_values[time_diffs[t - spk_len] - 1, month_locs[t]] != st_var.mdi):
                        
                            # if exceed critical values
                            if (np.abs(value_diffs[t - spk_len]) >= critical_values[time_diffs[t - spk_len] - 1, month_locs[t]]):

                                # are signs of two differences different
                                if (math.copysign(1, value_diffs[t]) != math.copysign(1, value_diffs[t - spk_len])):
                                    
                                    # are within spike differences small
                                    if (spk_len == 1) or\
                                    ((spk_len == 2) and (np.abs(value_diffs[t - spk_len + 1]) < critical_values[time_diffs[t - spk_len + 1] -1, month_locs[t]] / 2.)) or \
                                    ((spk_len == 3) and (np.abs(value_diffs[t - spk_len + 1]) < critical_values[time_diffs[t - spk_len + 1] -1, month_locs[t]] / 2.) and\
                                      (np.abs(value_diffs[t - spk_len + 2]) < critical_values[time_diffs[t - spk_len + 2] -1, month_locs[t]] / 2.)):
                                    
                                        # check if following value is valid
                                        if (value_diffs[t] != st_var.mdi) and (critical_values[time_diffs[t] - 1, month_locs[t]] != st_var.mdi) and\
                                            (value_diffs[t] != st_var.fdi):
                                            
                                            # and if at least critical value                                            
                                            if (np.abs(value_diffs[t]) >= critical_values[time_diffs[t] - 1, month_locs[t]]):
                                                
                                                # test if surrounding differences below 1/2 critical value
                                                if (np.abs(value_diffs[t - spk_len - 1]) <= critical_values[time_diffs[t - spk_len - 1] - 1, month_locs[t]] / 2.): 
                                                    if (np.abs(value_diffs[t + 1]) <= critical_values[time_diffs[t + 1] - 1, month_locs[t]] / 2.): 
                                                    
                                                        # set the flags
                                                        flags[ t - spk_len + 1 : t +1] = 1   

                                                        if plots or diagnostics:
                                                            
                                                            sc_diagnostics_and_plots(station.time.data, st_var.data, good_to_uncompress[0][t-spk_len+1], good_to_uncompress[0][t+1], start, variable, plots = plots)
                                                           

        station.qc_flags[good_to_uncompress, flag_col[v]] = flags
                                    
        flag_locs = np.where(station.qc_flags[:, flag_col[v]] != 0)

        if plots or diagnostics:
            utils.print_flagged_obs_number(logfile, "Spike", variable, len(flag_locs[0]) - prev_flag_number, noWrite = True) # additional flags
        else:
            utils.print_flagged_obs_number(logfile, "Spike", variable, len(flag_locs[0]) - prev_flag_number) # additional flags

        # copy flags into attribute
        st_var.flags[flag_locs] = 1
 
        # matches 030660 - but with adapted IDL
        # matches 030220 OK, but finds more but all are reasonable 1/9/14

        do_interactive = False
        if plots and do_interactive == True:
            import matplotlib.pyplot as plt
        
            plot_times = utils.times_hours_to_datetime(station.time.data, start)
            
            plt.clf()
            plt.plot(plot_times, all_filtered, 'bo', ls='-')
            flg = np.where(flags[:, flag_col[v]] == 1)
            plt.plot(plot_times[flg], all_filtered[flg], 'ro', markersize=10)
            plt.show()
	    
    station = utils.append_history(station, "Spike Check")  

    return # sc



#************************************************************************
if __name__ == "__main__":
    
    print "checking for short period spikes"
