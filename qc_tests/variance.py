#!/usr/local/sci/bin/python
#*****************************
#
# Excess Variance Check (EVC)
#   
#
#************************************************************************
#                    SVN Info
#$Rev:: 67                                            $:  Revision of last commit
#$Author:: rdunn                                      $:  Author of last commit
#$Date:: 2015-05-01 16:18:52 +0100 (Fri, 01 May 2015) $:  Date of last commit
#************************************************************************
import numpy as np
import scipy as sp
import datetime as dt
# RJHD routines

import qc_utils as utils

MAD_THRESHOLD = 4
DATA_COUNT_THRESHOLD = 120

#************************************************************************
def evc_plot_slp_wind(times, slp, slp_diffs, median_slp, slp_MAD, winds, median_wind, wind_MAD):
    # plot showing the pressure, pressure first differences and the wind speeds
    
    import matplotlib.pyplot as plt
    plt.clf()
    
    # pressure for this month
    plt.subplot(311)
    plt.plot(times, slp, 'bo')
    plt.axhline(median_slp, color = 'r')
    plt.axhline(median_slp - (MAD_THRESHOLD * slp_MAD), color = 'r', ls= '--')
    plt.ylabel("SLP (hPa)")
    
    # pressure first differences
    plt.subplot(312)

    plot_diffs = np.ma.zeros(len(slp))
    plot_diffs.mask = slp.mask
    good = np.where(slp.mask == False)
    plot_diffs[good] = slp_diffs

    plt.plot(times, plot_diffs,'go')
    plt.ylabel("1st diffs")
    
    # windspeeds
    plt.subplot(313)
    plt.plot(times, winds, 'bo')
    plt.axhline(median_wind, color = 'r')
    plt.axhline(median_wind + (MAD_THRESHOLD * wind_MAD), color = 'r', ls= '--')
    plt.ylabel("Wind Speed (m/s)")
    
    # showing median values (red solid) and the MAD offset (red dashed)
    plt.show()                             

    return # plot_slp_wind

#************************************************************************
def evc_plot_hist(plot_variances, iqr_threshold, title):
    '''
    Plot the histogram, with removed observations highlighted
    
    :param array plot_variances: values to be shown on histogram
    :param array iqr_threshold: threshold for removal
    :param str title: title of plot
    
    :returns:
    '''

    import matplotlib.pyplot as plt
    # set up the bins
    bins, bincenters = utils.create_bins(plot_variances, 1.0)

    # make the histogram
    hist, binEdges = np.histogram(plot_variances, bins = bins)
                
    plot_hist = np.array([0.01 if h == 0 else h for h in hist]) # allow for log y-scale

    plt.clf()
    plt.step(bincenters, plot_hist, 'k-', label = 'standardised months', where='mid')
    
    # sort the labels etc
    plt.xlabel("variance offset (IQR)")                    
    plt.ylabel("Frequency")
    plt.gca().set_yscale('log')

    plt.axvline(-iqr_threshold,c='r')
    plt.axvline(iqr_threshold,c='r')
    plt.step(bincenters[bincenters < -iqr_threshold], plot_hist[bincenters < -iqr_threshold], 'r-', where='mid')
    plt.step(bincenters[bincenters > iqr_threshold], plot_hist[bincenters > iqr_threshold], 'r-', where='mid')

    plt.ylim(ymin=0.1)
    plt.title(title)
                
    plt.show()

    return # plot_hist


#************************************************************************
def evc(station, variable_list, flag_col, start, end, logfile, diagnostics = False, plots = False, idl = False):
    
    if plots or diagnostics:
        import matplotlib.pyplot as plt
        import calendar

    
    # very similar to climatological check - ensure that not duplicating
    
    for v, variable in enumerate(variable_list):
    
        st_var = getattr(station, variable)
    
        reporting_resolution = utils.reporting_accuracy(utils.apply_filter_flags(st_var))
        reporting_freq = utils.reporting_frequency(utils.apply_filter_flags(st_var))
   
        month_ranges = utils.month_starts_in_pairs(start, end)
        month_ranges = month_ranges.reshape(-1,12,2)

        month_data_count = np.zeros(month_ranges.shape[0:2])

        # for each month
        for month in range(12):

            # set up hourly climatologies
            hourly_clims = np.zeros(24)
            hourly_clims.fill(st_var.data.fill_value)

            this_month, year_ids, month_data_count[:,month] = utils.concatenate_months(month_ranges[:,month,:], st_var.data, hours = True)

            
            # # extract each year and append together
            # year_ids = [] # counter to determine which year each day corresponds to
            # for year in range(month_ranges.shape[0]):
                
            #     this_year = st_var.data[month_ranges[year,month][0]:month_ranges[year,month][1]]
            #     if year == 0:
            #         # store so can access each hour of day separately
            #         this_month = this_year.reshape(-1,24)
                    
            #         year_ids = [year for x in range(this_month.shape[0])]
                    
            #         month_data_count[year,month] = len(this_year.compressed())
                    
            #     else:
            #         this_year = this_year.reshape(-1,24)
                       
            #         this_month = np.ma.concatenate((this_month, this_year), axis = 0)
                    
            #         year_ids.extend([year for x in range(this_year.shape[0])])
                    
            #         month_data_count[year,month] = len(this_year.compressed())

                
                  
            # winsorize and get hourly climatology 
            for h in range(24):
                
                this_hour = this_month[:,h]
                
                if len(this_hour.compressed()) > 100:

                    
                    # winsorize & climatologies - done to match IDL
                    if idl:
                        this_hour_winsorized = utils.winsorize(np.append(this_hour.compressed(), -999999), 0.05, idl = idl)
                        hourly_clims[h] = np.ma.sum(this_hour_winsorized)/(len(this_hour_winsorized) - 1)
                        
                    else:
                        this_hour_winsorized = utils.winsorize(this_hour.compressed(), 0.05, idl = idl)
                        hourly_clims[h] = np.ma.mean(this_hour_winsorized)
                    
            
            hourly_clims = np.ma.masked_where(hourly_clims == st_var.data.fill_value, hourly_clims)           
            anomalies = this_month - np.tile(hourly_clims, (this_month.shape[0], 1))
            
            # extract IQR of anomalies (using 1/2 value to match IDL)
            if len(anomalies.compressed()) >= 10:
                
                iqr = utils.IQR(anomalies.compressed().reshape(-1)) / 2. # to match IDL
                if iqr < 1.5: iqr = 1.5

            else:
                iqr = st_var.mdi
            
            normed_anomalies = anomalies / iqr
            

            variances = np.ma.zeros(month_ranges.shape[0])
            variances.mask = [False for i in range(month_ranges.shape[0])]
            rep_accuracies = np.zeros(month_ranges.shape[0])
            rep_freqs = np.zeros(month_ranges.shape[0])
            
            variances.fill(st_var.mdi)
            rep_accuracies.fill(st_var.mdi)
            rep_freqs.fill(st_var.mdi)
                
            year_ids = np.array(year_ids)
            
            # extract variance of normalised anomalies for each year
            for y, year in enumerate(range(month_ranges.shape[0])):
            
                year_locs = np.where(year_ids == y)
            
                this_year = normed_anomalies[year_locs,:]
                this_year = this_year.reshape(-1)
                
            # end of similarity with Climatological check
            
                if len(this_year.compressed()) >= 30:
            
                    variances[y] = utils.mean_absolute_deviation(this_year, median = True)
                    
                    rep_accuracies[y] = utils.reporting_accuracy(this_year)
                    rep_freqs[y] = utils.reporting_frequency(this_year)

                else:
                    variances.mask[y] = True

            good = np.where(month_data_count[:,month] >= 100)
            
            # get median and IQR of variance for all years for this month
            if len(good[0]) >= 10:
                
                median_variance = np.median(variances[good])
                
                iqr_variance = utils.IQR(variances[good]) / 2. # to match IDL
                
                if iqr_variance < 0.01: iqr_variance = 0.01
            else:
                
                median_variance = st_var.mdi
                iqr_variance = st_var.mdi

                
            # if SLP, then get median and MAD of SLP and windspeed for month
            if variable in ["slp", "windspeeds"]:
                
                winds = getattr(station, "windspeeds")
                slp = getattr(station, "slp")
        
                # refactor this as similar in style to how target data extracted  
                for y, year in enumerate(range(month_ranges.shape[0])):
                    
                    if y == 0:
                        winds_year = winds.data[month_ranges[year,month][0]:month_ranges[year,month][1]]
                        winds_month = winds_year.reshape(-1,24)
                                            
                        slp_year = slp.data[month_ranges[year,month][0]:month_ranges[year,month][1]]
                        slp_month = slp_year.reshape(-1,24)
                                            
                    else:
                        winds_year = winds.data[month_ranges[year,month][0]:month_ranges[year,month][1]]                  
                        winds_year = winds_year.reshape(-1,24)
                        winds_month = np.ma.concatenate((winds_month, winds_year), axis = 0)
                        
                        slp_year = slp.data[month_ranges[year,month][0]:month_ranges[year,month][1]]                  
                        slp_year =  slp_year.reshape(-1,24)
                        slp_month = np.ma.concatenate((slp_month, slp_year), axis = 0)
                        
                median_wind = np.ma.median(winds_month)
                median_slp  = np.ma.median(slp_month)
                
                wind_MAD = utils.mean_absolute_deviation(winds_month.compressed())
                slp_MAD = utils.mean_absolute_deviation(slp_month.compressed())
                
                if diagnostics:
                    print "median windspeed {} m/s, MAD = {}".format(median_wind, wind_MAD)
                    print "median slp {} hPa, MAD = {}".format(median_slp, slp_MAD)

            # now test to see if variance exceeds expected range
            for y, year in enumerate(range(month_ranges.shape[0])):


                if (variances[y] != st_var.mdi) and (iqr_variance != st_var.mdi) and \
                    (median_variance != st_var.mdi) and (month_data_count[y,month] >= DATA_COUNT_THRESHOLD):
                    
                    # if SLP, then need to test if deep low pressure ("hurricane/storm") present
                    #   as this will increase the variance for this month + year
                    if variable in ["slp", "windspeeds"]:
                        
                        iqr_threshold = 6.
                        
                        # increase threshold if reporting frequency and resolution of this
                        #   year doesn't match average
                        if (rep_accuracies[y] != reporting_resolution) and \
                            (rep_freqs[y] != reporting_freq):
                            iqr_threshold = 8.
                       
                        if diagnostics:
                            print np.abs(variances[y] - median_variance) / iqr_variance, variances[y] , median_variance , iqr_variance , iqr_threshold, month+1, year+start.year
                        
                        if np.abs((variances[y] - median_variance) / iqr_variance) > iqr_threshold:
                        
                            # check for storms     
                            winds_month = winds.data[month_ranges[year,month][0]:month_ranges[year,month][1]]                  
                            slp_month = slp.data[month_ranges[year,month][0]:month_ranges[year,month][1]]                  
                   
                            storm = False
                            if (len(winds_month.compressed()) >= 1) and (len(slp_month.compressed()) >= 1):
                                # find max wind & min SLP
                                # max_wind_loc = np.where(winds_month == np.max(winds_month))[0][0]
                                # min_slp_loc = np.where(slp_month == np.min(slp_month))[0][0]

                                # if these are above thresholds and within one day of each other,
                                #    then it likely was a storm
                                # print "fix this in case of multiple max/min locations"
                                # if (np.abs(max_wind_loc - min_slp_loc) <= 24) and \ 
                                #     (((np.max(winds_month) - median_wind) / wind_MAD) > MAD_THRESHOLD) and \
                                #     (((median_slp - np.min(slp_month)) / slp_MAD) > MAD_THRESHOLD): 

                                # locations where winds greater than threshold
                                high_winds, = np.where((winds_month - median_wind)/wind_MAD > MAD_THRESHOLD)
                                # and where SLP less than threshold
                                low_slps, = np.where((median_slp - slp_month)/slp_MAD > MAD_THRESHOLD)

                                # if any locations match, then it's a storm
                                match_loc = high_winds[np.in1d(high_winds, low_slps)]
                                    
                                if len(match_loc) > 0:
                                    storm = True
                            else:
                                print "write spurious"
                                
                            # check the SLP first difference series
                            #   to ensure a drop down and climb out of minimum SLP/or climb up and down from maximum wind speed
                            if variable == "slp":
                                diffs = np.diff(slp_month.compressed())
                            elif variable == "windspeeds":
                                diffs = np.diff(winds_month.compressed())
                            
                            negs, poss = 0,0
                            biggest_neg, biggest_pos = 0,0
                            
                            for diff in diffs:
                                
                                if diff > 0:
                                    if negs > biggest_neg: biggest_neg = negs
                                    negs = 0
                                    poss += 1
                                else:
                                    if poss > biggest_pos: biggest_pos = poss
                                    poss = 0
                                    negs += 1
                                
                            if (biggest_neg < 10) and (biggest_pos < 10) and not storm:
                                
                                # not a hurricane, so mask
                                station.qc_flags[month_ranges[year,month,0]:month_ranges[year,month,1], flag_col[v]] = 1
                                if plots or diagnostics:
                                    print "No Storm or Hurricane in %i %i - flagging\n" % (month+1, y+start.year)
                                else:
                                    logfile.write("No Storm or Hurricane in %i %i - flagging\n" % (month+1, y+start.year))
                                
                            else:
                                # hurricane
                                if plots or diagnostics:
                                    print "Storm or Hurricane in %i %i - not flagging\n" % (month+1, y+start.year)
                                else:
                                    logfile.write("Storm or Hurricane in %i %i - not flagging\n" % (month+1, y+start.year))
                        
                            if plots:
                                # plot showing the pressure, pressure first differences and the wind speeds
                                plot_times = utils.times_hours_to_datetime(station.time.data[month_ranges[year,month][0]:month_ranges[year,month][1]], start)

                                evc_plot_slp_wind(plot_times, slp_month, diffs, median_slp, slp_MAD, winds_month, median_wind, wind_MAD)

                    else:
                        
                        iqr_threshold = 8.
                        
                        if (rep_accuracies[y] != reporting_resolution) and \
                            (rep_freqs[y] != reporting_freq):
                            iqr_threshold = 10.
                            

                        if np.abs(variances[y] - median_variance) / iqr_variance > iqr_threshold:
                                
                            if diagnostics:
                                print "flagging {} {}".format(year+start.year,calendar.month_name[month+1])
                            # remove the data 
                            station.qc_flags[month_ranges[year,month,0]:month_ranges[year,month,1], flag_col[v]] = 1


            if plots:
                plot_variances = (variances - median_variance) / iqr_variance

                plot_variances = np.ma.masked_where(month_data_count[:,month] < DATA_COUNT_THRESHOLD,plot_variances)
                
                evc_plot_hist(plot_variances, iqr_threshold, "Variance Check - %s - %s" % (variable, calendar.month_name[month+1]))
 
        flag_locs = np.where(station.qc_flags[:, flag_col[v]] != 0)
        if plots or diagnostics:
            utils.print_flagged_obs_number(logfile, "Variance", variable, len(flag_locs[0]), noWrite = True)
        else:
            utils.print_flagged_obs_number(logfile, "Variance", variable, len(flag_locs[0]))
            
        # copy flags into attribute
        st_var.flags[flag_locs] = 1

    # matches 030660 for T, D and SLP 21/8/2014

    station = utils.append_history(station, "Excess Variance Check")

    return # evc

#************************************************************************
if __name__ == "__main__":

    print "Checking for High/Low Variance Months"
