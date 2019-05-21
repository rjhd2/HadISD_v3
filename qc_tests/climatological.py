#!/usr/local/sci/bin/python
#*****************************
#
# Climatological Outlier Check (COC)
#
#   Outliers against monthly hourly climatologies
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
import copy

# RJHD routines
import qc_utils as utils
import distributional_gap as dgc # very similar in parts, so avoid duplication of code


FREQUENCY_THRESHOLD = 0.1

IMAGELOCATION = "blank"

#************************************************************************
def coc_set_up_plot(bincenters, hist, gaussian, variable, threshold = 0, sub_par = ""):
    '''
    Set up the plotting space for the Climatological Outlier Check

    :param array bincenters: bin centres of histogram
    :param array hist: histogram values
    :param array gaussian: parameters of gaussian fit [m, s, n]
    :param str variable: name of variable for title
    :param int threshold: threshold to plot
    :param str sub_par: sub parameter for axis label
    '''   
    import matplotlib.pyplot as plt
    
    plt.clf()
    plt.axes([0.1,0.15,0.85,0.75])
    plot_hist = np.array([0.01 if h == 0 else h for h in hist])  
    plt.step(bincenters, plot_hist, 'k-', label = 'standardised months', where='mid')

    # plot fitted Gaussian
    plot_gaussian = utils.gaussian(bincenters, gaussian)
    plt.plot(bincenters, plot_gaussian, 'b-', label = 'Gaussian fit')

    # sort the labels etc
    plt.xlabel("%s offset (IQR)" % variable)                    
    plt.ylabel("Frequency (%s)" % sub_par)
    plt.gca().set_yscale('log')
    plt.axvline(-threshold-1,c='r')
    plt.axvline(threshold+1,c='r')
    plt.axvline(-threshold,c='orange')
    plt.axvline(threshold,c='orange')
    plt.ylim(ymin=0.1)
    plt.title("Climatological Gap Check - %s - %s" % (sub_par, variable) )        

    return  # coc_set_up_plot

#************************************************************************
def coc_get_weights(monthly_vqvs, monthly_subset, filter_subset):
    '''
    Get the weights for the low pass filter.

    :param array monthly_vqvs: monthly anomalies
    :param array monthly_subset: which values to take
    :param array filter_subset: which values to take
    :returns:
        weights
    '''

    filterweights = np.array([1.,2.,3.,2.,1.])

    if np.sum(filterweights[filter_subset] * np.ceil(monthly_vqvs[monthly_subset] - np.floor(monthly_vqvs[monthly_subset]))) == 0:
        weights = 0
    else:
        weights = np.sum(filterweights[filter_subset] * monthly_vqvs[monthly_subset]) / \
            np.sum(filterweights[filter_subset] * np.ceil(monthly_vqvs[monthly_subset] - np.floor(monthly_vqvs[monthly_subset])))

    return weights # coc_get_weights

#************************************************************************
def coc_low_pass_filter(normed_anomalies, year_ids, monthly_vqvs, years):
    '''
    Run the low pass filter - get suitable ranges, get weights, and apply

    :param array normed_anomalies: input normalised anomalies
    :param array year_ids: look up array for which row is which year
    :param array monthly_vqvs: monthly average anomalies
    :param int years: number of years in data
    :returns:
       normed_anomalies - with weights applied
    '''

    for year in range(years):
        
        if year == 0:
            monthly_range = np.arange(0,3)
            filter_range = np.arange(2,5)
        elif year == 1:
            monthly_range = np.arange(0,4)
            filter_range = np.arange(1,5)
        elif year == years - 2:
            monthly_range = np.arange(-4,0,-1)
            filter_range = np.arange(0,4)
        elif year == years - 1:
            monthly_range = np.arange(-3,0,-1)
            filter_range = np.arange(0,3)
        else:
            monthly_range = np.arange(year-2, year+3)
            filter_range = np.arange(5)

            
        if np.ma.sum(np.ma.abs(monthly_vqvs[monthly_range])) != 0:
                
            weights = coc_get_weights(monthly_vqvs, monthly_range, filter_range)
            
            year_locs = np.where(year_ids == year)
            normed_anomalies[year_locs,:] = normed_anomalies[year_locs,:] - weights

    return normed_anomalies # coc_low_pass_filter

#************************************************************************
def coc_find_and_apply_flags(month_ranges, normed_anomalies, flags, year_ids, threshold, gap_start, upper = True, plots = False, gpv = [], tpv = []):
    '''
    Find which values need flagging, do so, and extract for plotting

    :param array month_ranges: start and end locations for each month
    :param array normed_anomalies: data to flag
    :param array flags: flag array
    :param array year_ids: year look-up array
    :param int threshold: critical value to start flagging
    :param int gap_start: location of gap if exists, else 0
    :param bool upper: upper or lower part of distribution
    :param bool plots: do plots
    :param array gpv: gap plot values
    :param array tpb: tentative plot values

    :returns:
        flags - updated flag array
        gpv - updated values which have been flagged - for plotting
        tpv - updated values which have been tentatively flagged - for plotting
    '''

    for y, year in enumerate(month_ranges):

        year_locs = np.where(year_ids == y)
        this_year_data = normed_anomalies[year_locs[0],:]

        this_year_flags = np.array(flags[year[0]:year[1]])
        this_year_flags = this_year_flags.reshape(-1,24)

        if upper:
            tentative_locations = np.ma.where(this_year_data > threshold)
        else:
            tentative_locations = np.ma.where(this_year_data < -threshold)
        this_year_flags[tentative_locations] = 2

        if gap_start != 0:
            if upper:
                gap_cleaned_locations = np.ma.where(this_year_data > gap_start)
            else:
                gap_cleaned_locations = np.ma.where(this_year_data < gap_start)

            this_year_flags[gap_cleaned_locations] = 1

        this_year_flags = this_year_flags.reshape(-1)

        flags[year[0]:year[1]] = this_year_flags

        if plots: 
            if gap_start != 0:
                gpv = np.append(gpv, this_year_data[gap_cleaned_locations].compressed())
            tpv = np.append(tpv, this_year_data[tentative_locations].compressed())


    return flags, gpv, tpv

#************************************************************************
def coc(station, variable_list, flag_col, start, end, logfile, diagnostics = False, plots = False, idl = False, doMonth = False):
    
    
    for v, variable in enumerate(variable_list):
        
        st_var = getattr(station, variable)
        all_filtered = utils.apply_filter_flags(st_var, doMonth = doMonth, start = start, end = end)
        
        st_var_complete_year = copy.deepcopy(st_var)
        if doMonth:
            # restrict the incomplete year if appropriate - keep other flagged obs.
            full_year_end = utils.get_first_hour_this_year(start, end)
            st_var_complete_year.data.mask[full_year_end :] = True


        # is this needed 13th Nov 2014 RJHD
        #reporting_resolution = utils.reporting_accuracy(utils.apply_filter_flags(st_var))
        
        month_ranges = utils.month_starts_in_pairs(start, end)
        month_ranges = month_ranges.reshape(-1,12,2)
    
        for month in range(12):
            
            hourly_climatologies = np.zeros(24)
            hourly_climatologies.fill(st_var.mdi)
            
            # append all e.g. Januaries together

            this_month, year_ids, dummy = utils.concatenate_months(month_ranges[:,month,:], st_var.data, hours = True)
            this_month_complete, dummy, dummy = utils.concatenate_months(month_ranges[:,month,:], st_var_complete_year.data, hours = True)
            this_month_filtered, dummy, dummy = utils.concatenate_months(month_ranges[:,month,:], all_filtered, hours = True)

            # if fixed climatology period, sort this here
            
            # get as array of 24 hrs.  
            this_month = np.ma.array(this_month)
            this_month = this_month.reshape(-1,24)

            this_month_complete = np.ma.array(this_month_complete)
            this_month_complete = this_month_complete.reshape(-1,24)

            this_month_filtered = np.ma.array(this_month_filtered)
            this_month_filtered = this_month_filtered.reshape(-1,24)

            # get hourly climatology for each month
            for hour in range(24):
                
                this_hour = this_month_complete[:,hour]

                # need to have data if this is going to work!
                if len(this_hour.compressed()) > 0:

                    # winsorize & climatologies - done to match IDL
                    if idl:
                        this_hour = utils.winsorize(np.append(this_hour.compressed(), -999999), 0.05, idl = idl)
                        hourly_climatologies[hour] = np.ma.sum(this_hour)/(len(this_hour) - 1)

                    else:
                        this_hour = utils.winsorize(this_hour.compressed(), 0.05, idl = idl)
                        hourly_climatologies[hour] = np.ma.mean(this_hour)

            if diagnostics:
                print "hourly clims", hourly_climatologies

            if len(this_month.compressed()) > 0 and len(this_month_complete.compressed()) > 0:
                # can get stations with few obs in a particular variable.
                # or, with monthly running, if this variable has just started reporting
                #  only get data from the recent month, but not from previous year.

                # anomalise each hour over month appropriately

                anomalies = this_month - np.tile(hourly_climatologies, (this_month.shape[0],1))
                anomalies_complete = this_month_complete - np.tile(hourly_climatologies, (this_month_complete.shape[0],1))
                anomalies_filtered = this_month_filtered - np.tile(hourly_climatologies, (this_month_filtered.shape[0],1))

                if len(anomalies_complete.compressed()) >= 10:
                    iqr = utils.IQR(anomalies_complete.compressed().reshape(-1))/2.  # to match IDL
                    if iqr < 1.5: iqr = 1.5
                else:
                    iqr = st_var.mdi

                normed_anomalies = anomalies / iqr
                normed_anomalies_complete = anomalies_complete / iqr
                normed_anomalies_filtered = anomalies_filtered / iqr

                if diagnostics:
                    print np.ma.mean(this_month), np.ma.mean(this_month_complete), np.ma.mean(this_month_filtered)
                    print np.ma.mean(anomalies), np.ma.mean(anomalies_complete), np.ma.mean(anomalies_filtered)
                    print np.ma.mean(normed_anomalies), np.ma.mean(normed_anomalies_complete), np.ma.mean(normed_anomalies_filtered)


                # get average anomaly for year
                year_ids = np.array(year_ids)
                monthly_vqvs = np.ma.zeros(month_ranges.shape[0])
                monthly_vqvs.mask = [False for x in range(month_ranges.shape[0])]
                for year in range(month_ranges.shape[0]):
                    year_locs = np.where(year_ids == year)
                    this_year = normed_anomalies_filtered[year_locs,:]

                    if len(this_year.compressed()) > 0:
                        # need to have data for this to work!
                        if idl:
                            monthly_vqvs[year] = utils.idl_median(this_year.compressed().reshape(-1))
                        else:
                            monthly_vqvs[year] = np.ma.median(this_year)
                    else:
                        monthly_vqvs.mask[year] = True

                if diagnostics:
                    print "monthly vqvs", monthly_vqvs

                # low pass filter
                normed_anomalies = coc_low_pass_filter(normed_anomalies, year_ids, monthly_vqvs, month_ranges.shape[0])

                if doMonth:
                    # run low pass filter, ignoring the final incomplete year.
                    not_final_year_locs, = np.where(year_ids != year_ids[-1])

                    normed_anomalies_complete[not_final_year_locs] = coc_low_pass_filter(normed_anomalies_complete[not_final_year_locs], year_ids[not_final_year_locs], monthly_vqvs[:-1], month_ranges.shape[0]-1)
                else:
                    normed_anomalies_complete = coc_low_pass_filter(normed_anomalies_complete, year_ids, monthly_vqvs, month_ranges.shape[0])


                # copy from distributional_gap.py - refactor!
                # get the threshold value using complete values
                bins, bincenters = utils.create_bins(normed_anomalies_complete, 1.)
                
                hist, binEdges = np.histogram(normed_anomalies_complete.compressed(), bins = bins)

                gaussian = utils.fit_gaussian(bincenters, hist, max(hist), mu=np.ma.mean(normed_anomalies_complete), sig = np.ma.std(normed_anomalies_complete))
                minimum_threshold = round(1. + utils.invert_gaussian(FREQUENCY_THRESHOLD, gaussian))

                if diagnostics:
                    print iqr, minimum_threshold, 1. + utils.invert_gaussian(FREQUENCY_THRESHOLD, gaussian)
                    print gaussian
                    print hist
                    print bins

                if plots:
                    coc_set_up_plot(bincenters, hist, gaussian, variable, threshold = minimum_threshold, sub_par = "observations")

                # apply to uncomplete values
                uppercount = len(np.where(normed_anomalies > minimum_threshold)[0])
                lowercount = len(np.where(normed_anomalies < -minimum_threshold)[0])

                these_flags = station.qc_flags[:, flag_col[v]]
                gap_plot_values, tentative_plot_values = [], []

                # find the gaps and apply the flags

                gap_start = dgc.dgc_find_gap(hist, binEdges, minimum_threshold, gap_size = 1) # in DGC it is 2.
                these_flags, gap_plot_values, tentative_plot_values =\
                    coc_find_and_apply_flags(month_ranges[:,month,:],normed_anomalies, these_flags, year_ids, minimum_threshold, gap_start, \
                                                           upper = True, plots = plots, gpv = gap_plot_values, tpv = tentative_plot_values)

                gap_start = dgc.dgc_find_gap(hist, binEdges, -minimum_threshold, gap_size = 1) # in DGC it is 2.
                these_flags, gap_plot_values, tentative_plot_values =\
                    coc_find_and_apply_flags(month_ranges[:,month,:],normed_anomalies, these_flags, year_ids, minimum_threshold, gap_start, \
                                                           upper = False, plots = plots, gpv = gap_plot_values, tpv = tentative_plot_values)

                station.qc_flags[:, flag_col[v]] = these_flags

                if uppercount + lowercount > 1000:
                    #print "not sorted spurious stations yet"
                    pass
                if plots:
                    import matplotlib.pyplot as plt
                    hist, binEdges = np.histogram(tentative_plot_values, bins = bins)
                    plot_hist = np.array([0.01 if h == 0 else h for h in hist])
                    plt.step(bincenters, plot_hist, c='orange', ls='-', label = 'tentative', where='mid')

                    hist, binEdges = np.histogram(gap_plot_values, bins = bins)
                    plot_hist = np.array([0.01 if h == 0 else h for h in hist])
                    plt.step(bincenters, plot_hist, 'r-', label = 'flagged', where='mid')
                    import calendar
                    plt.text(0.1,0.9,calendar.month_name[month+1], transform = plt.gca().transAxes)
                    leg=plt.legend(loc='lower center',ncol=4, bbox_to_anchor=(0.5,-0.2),frameon=False,prop={'size':13},labelspacing=0.15,columnspacing=0.5)
                    plt.setp(leg.get_title(), fontsize=14)
                    plt.show()
                    #plt.savefig(IMAGELOCATION+'/'+station.id+'_ClimatologicalGap_'+str(month+1)+'.png')


        
        flag_locs = np.where(station.qc_flags[:, flag_col[v]] != 0)

        # copy flags into attribute
        st_var.flags[flag_locs] = 1

        utils.print_flagged_obs_number(logfile, "Climatological", variable, len(flag_locs[0]), noWrite = diagnostics)
        if diagnostics: print "where\n"
        logfile.write("where\n")
        nflags = len(np.where(station.qc_flags[:, flag_col[v]] == 1)[0])
        utils.print_flagged_obs_number(logfile, "  Firm Clim", variable, nflags, noWrite = diagnostics)
        nflags = len(np.where(station.qc_flags[:, flag_col[v]] == 2)[0])
        utils.print_flagged_obs_number(logfile, "  Tentative Clim", variable, nflags, noWrite = diagnostics)
 
        # firm flags match 030220
    station = utils.append_history(station, "Climatological Check")  
                     
    return

#************************************************************************
if __name__ == "__main__":

    print "Checking for climatological outliers"
