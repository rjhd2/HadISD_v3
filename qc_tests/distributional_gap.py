#!/usr/local/sci/bin/python
#*****************************
#
# Distributional Gap Check (DGC)
#
#   At times this is a direct translation from IDL
#     Could be made more pythonic, but need to match outputs!
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
from scipy.optimize import leastsq
import scipy.stats as stats
import copy

# RJHD routines
import qc_utils as utils


# mean vs median
# SD vs IQR
IMAGELOCATION = "blank"
# need to explain all of these
OBS_LIMIT = 50 
VALID_MONTHS = 5
SPREAD_LIMIT = 2
MONTH_LIMIT = 60
LARGE_LIMIT = 5 # IQR
MEAN = False
NON_GAP_FLAG = False
FREQUENCY_THRESHOLD = 0.1
GAP_SIZE = 2

BIN_SIZE = 0.25

MAD_THRESHOLD = 4.5

#************************************************************************
def dgc_get_monthly_averages(data, limit, mdi, MEAN = False):

    if len(data.compressed()) >= limit:
        if MEAN:
            return np.ma.mean(data)
        else:
            return np.ma.median(data)
        
    else:
        return mdi # dgc_get_monthly_averages
                
#************************************************************************
def dgc_expand_storms(storms, maximum):

    for i in range(6):

        if storms[0]-1 >=0:
            storms = np.insert(storms, 0, storms[0]-1)
        if storms[-1]+1 < maximum:
            storms = np.append(storms, storms[-1]+1)

    return np.unique(storms) # dgc_expand_storms
                
#************************************************************************
def dgc_find_gap(hist, bins, threshold, gap_size = GAP_SIZE):
    '''
    Walk the bins of the distribution to find a gap and return where it starts
   
    :param array hist: histogram values
    :param array bins: bin values
    :param flt threshold: limiting value
    :param int gap_size: gap size to record
    :returns:
        flt: gap_start
    '''

    start = np.argmax(hist)

    if bins[start] < threshold:
        positive = True
    else:
        positive = False
        
    n = 0
    gap_length = 0
    gap_start = 0
    while True:
        if hist[start + n] == 0:
            gap_length += 1
            if gap_start == 0:
                # plus 1 to get upper bin boundary
                if (positive and bins[start + n + 1] >= threshold):
                    gap_start = bins[start + n + 1]
                elif (not positive and bins[start + n] <= threshold):
                    gap_start = bins[start + n]
                
        else:
            if gap_length < gap_size:
                gap_length = 0
                
            elif gap_length >= gap_size and gap_start != 0:
                break
            
        if (start + n == len(hist) - 1) or (start + n == 0):
            break
        
        if positive:
            n += 1
        else:
            n -= 1

    return gap_start # dgc_find_gap

#************************************************************************
def dgc_set_up_plot(plot_gaussian, standardised_months, variable, threshold = (1.5,-1.5), sub_par = "", GH = False):
    '''
    Set up the histogram plot and the Gaussian Fit.

    :param array standardised_months: input array of months standardised by IQR
    :param str variable: label for title and axes
    :param int threshold: x values to draw vertical lines
    :param str sub_par: sub-parameter for labels
    :returns:
    '''

    # set up the bins
    bins, bincenters = utils.create_bins(standardised_months, BIN_SIZE)
    dummy, plot_bincenters = utils.create_bins(standardised_months, BIN_SIZE/10.)

    # make the histogram
    hist, binEdges = np.histogram(standardised_months, bins = bins)
    plot_hist = np.array([0.01 if h == 0 else h for h in hist]) # allow for log y-scale
    
    import matplotlib.pyplot as plt

    plt.clf()
    plt.axes([0.1,0.15,0.8,0.7])
    plt.step(bincenters, plot_hist, 'k-', label = 'standardised months', where='mid')

    # # plot fitted Gaussian
    # if GH:
    #     initial_values = [np.max(hist), np.mean(standardised_months), np.std(standardised_months), stats.skew(standardised_months), stats.kurtosis(standardised_months)] # norm, mean, std, skew, kurtosis
        
    #     fit = leastsq(utils.residualsGH, initial_values, [bincenters, hist, np.ones(len(hist))])
    #     res = utils.hermite2gauss(fit[0])
        
    #     bins, bincenters = utils.create_bins(standardised_months, 0.025)
    #     plot_gaussian = utils.funcGH(fit[0], bincenters)

    # else:

    #     fit = utils.fit_gaussian(bincenters, hist, max(hist), mu = np.mean(standardised_months), sig = np.std(standardised_months))
    #     bins, bincenters = utils.create_bins(standardised_months, 0.025)
    #     plot_gaussian = utils.gaussian(bincenters, fit)


    plt.plot(plot_bincenters, plot_gaussian, 'b-', label = 'Gaussian fit')

    # sort the labels etc
    plt.xlabel("%s offset (IQR)" % variable)                    
    plt.ylabel("Frequency (%s)" % sub_par)
    plt.gca().set_yscale('log')
    plt.axvline(threshold[0],c='r')
    plt.axvline(threshold[1],c='r')
    plt.ylim(ymin=0.1)
    plt.title("Distributional Gap Check - %s - %s" % (sub_par, variable) )        

    return  # dgc_set_up_plot

#************************************************************************
def dgc_monthly(station, variable, flags, start, end, logfile, plots=False, diagnostics=False, idl = False, doMonth = False):
    '''
    Original Distributional Gap Check

    :param obj station: station object
    :param str variable: variable to act on
    :param array flags: flags array
    :param datetime start: data start
    :param datetime end: data end
    :param file logfile: output logfile
    :param bool plots: run plots
    :param bool diagnostics: run diagnostics
    :param bool idl: run IDL equivalent routines for median
    :returns: 
       flags - updated flag array
    '''

    if plots:
        import matplotlib.pyplot as plt
    
    st_var = getattr(station, variable)
    
    month_ranges = utils.month_starts_in_pairs(start, end)
    
    # get monthly averages
    month_average = np.empty(month_ranges.shape[0])
    month_average.fill(st_var.mdi)
    month_average_filtered = np.empty(month_ranges.shape[0])
    month_average_filtered.fill(st_var.mdi)
    
    all_filtered = utils.apply_filter_flags(st_var, doMonth = doMonth, start = start, end = end)

    for m, month in enumerate(month_ranges):
        
        data = st_var.data[month[0]:month[1]]
        
        filtered = all_filtered[month[0]:month[1]]
        
        month_average[m] = dgc_get_monthly_averages(data, OBS_LIMIT, st_var.mdi, MEAN)
        month_average_filtered[m] = dgc_get_monthly_averages(filtered, OBS_LIMIT, st_var.mdi, MEAN)
            
    # get overall monthly climatologies - use filtered data
    
    month_average = month_average.reshape(-1,12)
    month_average_filtered = month_average_filtered.reshape(-1,12)
    
    standardised_months = np.empty(month_average.shape)
    standardised_months.fill(st_var.mdi)
    
    for m in range(12):
        
        valid_filtered = np.where(month_average_filtered[:,m] != st_var.mdi)
        
        if len(valid_filtered[0]) >= VALID_MONTHS:
            
            valid_data = month_average_filtered[valid_filtered,m][0]
            
            if MEAN:
                clim = np.mean(valid_data)
                spread = np.stdev(valid_data)
                
            else:        
                if idl:
                    clim = utils.idl_median(valid_data.compressed().reshape(-1))
                else:
                    clim = np.median(valid_data)
                spread = utils.IQR(valid_data)
                if spread <= SPREAD_LIMIT:
                    spread = SPREAD_LIMIT
                    
            standardised_months[valid_filtered,m] = (month_average[valid_filtered,m] - clim) / spread 
                    
    standardised_months = standardised_months.reshape(month_ranges.shape[0]) 
    
    good_months = np.where(standardised_months != st_var.mdi)

    # must be able to do this with masked arrays
    if plots:
        bins, bincenters = utils.create_bins(standardised_months[good_months], BIN_SIZE)
        dummy, plot_bincenters = utils.create_bins(standardised_months[good_months], BIN_SIZE/10.)

        hist, binEdges = np.histogram(standardised_months[good_months], bins = bins)   

        fit = utils.fit_gaussian(bincenters, hist, max(hist), mu = np.mean(standardised_months[good_months]), sig = np.std(standardised_months[good_months]))
        plot_gaussian = utils.gaussian(plot_bincenters, fit)

        dgc_set_up_plot(plot_gaussian, standardised_months[good_months], variable, sub_par = "Months")
        
    # remove all months with a large standardised offset
        
    if len(good_months[0]) >= MONTH_LIMIT:
                
        standardised_months = np.ma.masked_values(standardised_months, st_var.mdi)
        large_offsets = np.where(standardised_months >= LARGE_LIMIT)

        if len(large_offsets[0]) > 0:
            
            for lo in large_offsets[0]:
                flags[month_ranges[lo,0]:month_ranges[lo,1]] = 1
                
            if plots:
                
                hist, binEdges = np.histogram(standardised_months[large_offsets], bins = bins)
                plot_hist = np.array([0.01 if h == 0 else h for h in hist])
                plt.step(bincenters, plot_hist, 'g-', label = '> %i' % LARGE_LIMIT, where = 'mid', zorder = 5)
                
                plt.axvline(5,c='g')
                plt.axvline(-5,c='g')



        # walk distribution from centre and see if any assymetry
        sort_order = standardised_months[good_months].argsort()

        mid_point = len(good_months[0]) / 2
        
        good = True
        iter = 1
        while good:
            
            if standardised_months[good_months][sort_order][mid_point - iter] != standardised_months[good_months][sort_order][mid_point + iter]:
                # using IDL notation
                tempvals = [np.abs(standardised_months[good_months][sort_order][mid_point - iter]),np.abs(standardised_months[good_months][sort_order][mid_point + iter])]
                
                if min(tempvals) != 0:
                    if max(tempvals)/min(tempvals) >= 2. and min(tempvals) >= 1.5:
                        # substantial asymmetry in distribution - at least 1.5 from centre and difference of 2.
                        
                        if tempvals[0] == max(tempvals):
                            # LHS
                            bad = good_months[0][sort_order][:mid_point - iter]
                            if plots: badplot = standardised_months[good_months][sort_order][:mid_point - iter]
                        elif tempvals[1] == max(tempvals):
                            #RHS
                            bad = good_months[0][sort_order][mid_point + iter:]
                            if plots: badplot = standardised_months[good_months][sort_order][mid_point + iter:]
                            
                        for b in bad:
                            flags[month_ranges[b,0]:month_ranges[b,1]] = 1
                
                        if plots:
                            
                            hist, binEdges = np.histogram(badplot, bins = bins)
                            plot_hist = np.array([0.01 if h == 0 else h for h in hist])
                            plt.step(bincenters, plot_hist, 'r-', label = 'Gap', where = 'mid', zorder = 4)
                
                        good = False        
                            
                
            iter += 1
            if iter == mid_point: break
                
                          
        if plots: 
            plt.legend(loc='lower center',ncol=4, bbox_to_anchor=(0.5,-0.2),frameon=False,prop={'size':13})
            plt.show()
            #plt.savefig(IMAGELOCATION+'/'+station.id+'_DistributionalGap.png')

    nflags, = np.where(flags != 0)
    utils.print_flagged_obs_number(logfile, "Distributional Gap Month", variable, len(nflags), noWrite=diagnostics)
                   
    return flags # dgc_monthly

#************************************************************************
def dgc_all_obs(station, variable, flags, start, end, logfile, plots = False, diagnostics = False, idl = False, windspeeds = False, GH = False, doMonth = False):
    '''RJHD addition working on all observations'''
    
    if plots:
        import matplotlib.pyplot as plt

    month_ranges = utils.month_starts_in_pairs(start, end)
    month_ranges = month_ranges.reshape(-1,12,2)

    # extract variable
    st_var = getattr(station, variable)
    # apply flags (and mask incomplete year if appropriate)
    all_filtered = utils.apply_filter_flags(st_var, doMonth = doMonth, start = start, end = end)
    
    st_var_complete_year = copy.deepcopy(st_var)
    if doMonth:
        # restrict the incomplete year if appropriate - keep other flagged obs.
        full_year_end = utils.get_first_hour_this_year(start, end)
        st_var_complete_year.data.mask[full_year_end :] = True
 

    for month in range(12):
    
        # if requiring wind data, extract data and find monthly averages
        if windspeeds == True:
            st_var_wind = getattr(station, "windspeeds")

            if doMonth:
                # restrict the incomplete year if appropriate               
                st_var_wind.data.mask[full_year_end :] = True

            # get monthly averages
            windspeeds_month = np.empty([])
            for y, year in enumerate(month_ranges[:,month,:]):
            
                if y == 0:
                    windspeeds_month = np.ma.array(st_var_wind.data[year[0]:year[1]])
                else:
                    windspeeds_month = np.ma.concatenate([windspeeds_month, st_var_wind.data[year[0]:year[1]]])
                  
            windspeeds_month_average = dgc_get_monthly_averages(windspeeds_month, OBS_LIMIT, st_var_wind.mdi, MEAN)
            windspeeds_month_mad = utils.mean_absolute_deviation(windspeeds_month, median=True)
    
                
        # pull data from each calendar month together
        this_month_data, dummy, dummy = utils.concatenate_months(month_ranges[:,month,:], st_var.data, hours = False)
        this_month_filtered, dummy, dummy = utils.concatenate_months(month_ranges[:,month,:], all_filtered, hours = False)
        this_month_complete, dummy, dummy = utils.concatenate_months(month_ranges[:,month,:], st_var_complete_year.data, hours = False)
                
        # if enough clean and complete data for this calendar month find the median and IQR
        if len(this_month_filtered.compressed()) > OBS_LIMIT:
            
            if idl:
                monthly_median = utils.idl_median(this_month_filtered.compressed().reshape(-1))
            else:
                monthly_median = np.ma.median(this_month_filtered)
                  
            iqr = utils.IQR(this_month_filtered.compressed())
            
            
            if iqr == 0.0:
                # to get some spread if IQR too small                   
                iqr = utils.IQR(this_month_filtered.compressed(), percentile = 0.05)  
                print "Spurious_stations file not yet sorted"
    
                
            # if have an IQR, anomalise using median and standardise using IQR
            if iqr != 0.0:               

                monthly_values = np.ma.array((this_month_data.compressed() - monthly_median) / iqr)
                complete_values = np.ma.array((this_month_complete.compressed() - monthly_median) / iqr)

                # use complete years only for the histogram - aiming to find outliers.
                bins, bincenters = utils.create_bins(complete_values, BIN_SIZE)
                dummy, plot_bincenters = utils.create_bins(complete_values, BIN_SIZE/10.)
                hist, binEdges = np.histogram(complete_values, bins = bins)
                """
                Change to monthly updates Oct 2017
                Thought about changing distribution to use filtered values
                But this changes the test beyond just dealing with additional months
                Commented out lines below would be alternative.
                """
                # bins, bincenters = utils.create_bins(filtered_values, BIN_SIZE)
                # dummy, plot_bincenters = utils.create_bins(filtered_values, BIN_SIZE/10.)
                # hist, binEdges = np.histogram(filtered_values, bins = bins)
                
                # used filtered (incl. incomplete year mask) to determine the distribution.                
                if GH:
                    # Use Gauss-Hermite polynomials to add skew and kurtosis to Gaussian fit - January 2015 ^RJHD

                    # Feb 2019 - if large amounts off centre, can affect initial values
                    # switched to median and MAD
                    initial_values = [np.max(hist), np.median(complete_values), utils.mean_absolute_deviation(complete_values, median=True), stats.skew(complete_values), stats.kurtosis(complete_values)] # norm, mean, std, skew, kurtosis
                    
                    fit = leastsq(utils.residualsGH, initial_values, [bincenters, hist, np.ones(len(hist))])
                    res = utils.hermite2gauss(fit[0], diagnostics = diagnostics)
                    
                    plot_gaussian = utils.funcGH(fit[0], plot_bincenters)

                    # adjust to remove the rising bumps seen in some fits - artefacts of GH fitting?
                    mid_point = np.argmax(plot_gaussian)
                    bad, = np.where(plot_gaussian[mid_point:] < FREQUENCY_THRESHOLD/10.)
                    if len(bad) > 0: plot_gaussian[mid_point:][bad[0]:] = FREQUENCY_THRESHOLD/10.

                    bad, = np.where(plot_gaussian[:mid_point] < FREQUENCY_THRESHOLD/10.)
                    if len(bad) > 0: plot_gaussian[:mid_point][:bad[-1]] = FREQUENCY_THRESHOLD/10.                   

                    # extract threshold values
                    good_values = np.argwhere(plot_gaussian > FREQUENCY_THRESHOLD)

                    l_minimum_threshold = round(plot_bincenters[good_values[0]]) - 1
                    u_minimum_threshold = 1 + round(plot_bincenters[good_values[-1]])
                                      
                    if diagnostics:
                        print hist
                        print res
                        print iqr, l_minimum_threshold, u_minimum_threshold


                # or just a standard Gaussian
                else:
                    gaussian = utils.fit_gaussian(bincenters, hist, max(hist), mu = np.median(complete_values), sig = utils.mean_absolute_value(complete_values))

                    # assume the same threshold value
                    u_minimum_threshold = 1 + round(utils.invert_gaussian(FREQUENCY_THRESHOLD, gaussian))
                    l_minimum_threshold = -u_minimum_threshold


                    plot_gaussian = utils.gaussian(plot_bincenters, gaussian)

                    if diagnostics:
                        print hist
                        print gaussian
                        print iqr, u_minimum_threshold, 1. + utils.invert_gaussian(FREQUENCY_THRESHOLD, gaussian)

                if plots:
                    dgc_set_up_plot(plot_gaussian, complete_values, variable, threshold = (u_minimum_threshold, l_minimum_threshold), sub_par = "observations", GH = GH)
                     
                    if GH:
                        plt.figtext(0.15, 0.67, 'Mean %.2f, S.d. %.2f,\nSkew %.2f, Kurtosis %.2f' %(res['mean'], res['dispersion'], res['skewness'], res['kurtosis']), color='k', size='small')

                    
                # now trying to find gaps in the distribution
                uppercount = len(np.where(monthly_values > u_minimum_threshold)[0])
                lowercount = len(np.where(monthly_values < l_minimum_threshold)[0])
                
                # this needs refactoring - but lots of variables to pass in
                if plots or diagnostics: gap_plot_values = np.array([])

                # do one side of distribution and then other
                if uppercount > 0:
                    gap_start = dgc_find_gap(hist, binEdges, u_minimum_threshold)
                        
                    if gap_start != 0:
                        
                        # if found a gap, then go through each year for this calendar month
                        #  and flag observations further from middle
                        for y, year in enumerate(month_ranges[:,month,:]):
                
                            # not using filtered - checking all available data
                            this_year_data = np.ma.array(st_var.data[year[0]:year[1]])
                            this_year_flags = np.array(flags[year[0]:year[1]])
                            gap_cleaned_locations = np.ma.where(((this_year_data - monthly_median) / iqr) > gap_start)

                            this_year_flags[gap_cleaned_locations] = 1
                            flags[year[0]:year[1]] = this_year_flags

                            if plots or diagnostics: 
                                gap_plot_values = np.append(gap_plot_values, (this_year_data[gap_cleaned_locations].compressed() - monthly_median)/iqr)

                                if len(gap_cleaned_locations[0]) > 0:
                                    print "Upper {}-{} - {} obs flagged".format(y+start.year, month, len(gap_cleaned_locations[0]))
                                    print gap_cleaned_locations, this_year_data[gap_cleaned_locations]

                if lowercount > 0:
                    gap_start = dgc_find_gap(hist, binEdges, l_minimum_threshold)

                    if gap_start != 0:

                        # if found a gap, then go through each year for this calendar month
                        #  and flag observations further from middle
                        for y, year in enumerate(month_ranges[:,month,:]):
                
                            this_year_data = np.ma.array(st_var.data[year[0]:year[1]])
                            this_year_flags = np.array(flags[year[0]:year[1]])
                            gap_cleaned_locations = np.ma.where(np.logical_and(((this_year_data - monthly_median) / iqr) < gap_start, this_year_data.mask != True))
                            # add flag requirement for low pressure bit if appropriate

                            this_year_flags[gap_cleaned_locations] = 1
                            flags[year[0]:year[1]] = this_year_flags

                            if plots or diagnostics: 
                                gap_plot_values = np.append(gap_plot_values, (this_year_data[gap_cleaned_locations].compressed() - monthly_median)/iqr)

                                if len(gap_cleaned_locations[0]) > 0:
                                    print "Lower {}-{} - {} obs flagged".format(y+start.year, month, len(gap_cleaned_locations[0]))
                                    print gap_cleaned_locations, this_year_data[gap_cleaned_locations]

                            # if doing SLP then do extra checks for storms
                            if windspeeds:
                                windspeeds_year = np.ma.array(st_var_wind.data[year[0]:year[1]])

                                this_year_flags[gap_cleaned_locations] = 2 # tentative flags
                                
                                slp_average = dgc_get_monthly_averages(this_month_data, OBS_LIMIT, st_var.mdi, MEAN)
                                slp_mad = utils.mean_absolute_deviation(this_month_data, median=True)
                                
                                # need to ensure that this_year_data is less than slp_average, hence order of test
                                storms, = np.ma.where((((windspeeds_year - windspeeds_month_average) / windspeeds_month_mad) > MAD_THRESHOLD) &\
                                                   (((slp_average - this_year_data) / slp_mad) > MAD_THRESHOLD))
                                
                                # using IDL terminology
                                if len(storms) >= 2:
                                    # use the first difference series to find when there are gaps in 
                                    # contiguous sequences of storm observations - want to split up into
                                    # separate storm events                                    
                                    storm_1diffs = np.diff(storms)                                    
                                    separations, = np.where(storm_1diffs != 1)

                                    # expand around storm signal so that all low SLP values covered, and unflagged
                                    if len(separations) >= 1:
                                        print "  multiple storms in {} {}".format(y+start.year, month)

                                        # if more than one storm signal that month, then use intervals
                                        #    in the first difference series to expand around the first interval alone
                                        storm_start = 0
                                        storm_finish = separations[0] + 1                                           
                                        first_storm = dgc_expand_storms(storms[storm_start: storm_finish], len(this_year_data))
                                        final_storms = copy.deepcopy(first_storm)

                                        for j in range(len(separations)):
                                            # then do the rest in a loop
                                                
                                            if j+1 == len(separations):
                                                # final one
                                                this_storm = dgc_expand_storms(storms[separations[j]+1: ], len(this_year_data))
                                            else:
                                                this_storm = dgc_expand_storms(storms[separations[j]+1: separations[j+1]+1], len(this_year_data))
                                                
                                            final_storms = np.append(final_storms, this_storm)
                                        
                                    else:
                                        # else just expand around the signal by 6 hours either way
                                        final_storms = dgc_expand_storms(storms, len(this_year_data))

                                else:
                                    final_storms = storms

                                if len(storms) >= 1:
                                    print "Tropical Storm signal in {} {}".format(y+start.year, month)
                                    this_year_flags[final_storms] = 0

                            # and write flags back into array
                            flags[year[0]:year[1]] = this_year_flags

                if plots:
                    hist, binEdges = np.histogram(gap_plot_values, bins = bins)
                    plot_hist = np.array([0.01 if h == 0 else h for h in hist])
                    plt.step(bincenters, plot_hist, 'r-', label = 'flagged', where='mid')
                    import calendar
                    plt.text(0.1,0.9,calendar.month_name[month+1], transform = plt.gca().transAxes)
                    plt.legend(loc='lower center',ncol=3, bbox_to_anchor=(0.5,-0.2),frameon=False,prop={'size':13})
                    plt.show()
                    #plt.savefig(IMAGELOCATION+'/'+station.id+'_DistributionalGap_'+str(month+1)+'.png')


    nflags, = np.where(flags != 0)
    utils.print_flagged_obs_number(logfile, "Distributional Gap All", variable, len(nflags), noWrite=diagnostics)

    return flags # dgc_all_obs

#************************************************************************
def dgc(station, variable_list, flag_col, start, end, logfile, plots=False, diagnostics = False, idl = False, GH = False, doMonth = False):
    '''Controller for two individual tests'''

    if plots:
        import matplotlib.pyplot as plt
        
    
    for v, variable in enumerate(variable_list):    
        station.qc_flags[:,flag_col[v]] = dgc_monthly(station, variable, station.qc_flags[:,flag_col[v]], start, end, logfile, plots=plots, diagnostics=diagnostics, idl = idl)
        
        if variable == "slp":
            # need to send in windspeeds too        
            station.qc_flags[:,flag_col[v]] = dgc_all_obs(station, variable, station.qc_flags[:,flag_col[v]], start, end, logfile, plots=plots, diagnostics=diagnostics, idl = idl, windspeeds = True, GH = GH, doMonth = doMonth)

        else:
            station.qc_flags[:,flag_col[v]] = dgc_all_obs(station, variable, station.qc_flags[:,flag_col[v]], start, end, logfile, plots=plots, diagnostics=diagnostics, idl = idl, GH = GH, doMonth = doMonth)

    
        flag_locs = np.where(station.qc_flags[:, flag_col[v]] != 0)

        
        utils.print_flagged_obs_number(logfile, "Distributional Gap", variable, len(flag_locs[0]), noWrite=diagnostics)


        # copy flags into attribute
        st_var = getattr(station, variable)
        st_var.flags[flag_locs] = 1
 
        # MATCHES IDL for 030660-99999, 2 flags in T, 30-06-2014

    station = utils.append_history(station, "Distributional Gap Check")  
                     
    return # dgc

#************************************************************************

if __name__ == "__main__":

    print "running distributional gap check"
