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
#$Rev:: 67                                            $:  Revision of last commit
#$Author:: rdunn                                      $:  Author of last commit
#$Date:: 2015-05-01 16:18:52 +0100 (Fri, 01 May 2015) $:  Date of last commit
#************************************************************************

import numpy as np
import scipy as sp
import datetime as dt
from scipy.optimize import leastsq
import scipy.stats as stats

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
def dgc_monthly(station, variable, flags, start, end, plots=False, diagnostics=False, idl = False):
    '''
    Original Distributional Gap Check

    :param obj station: station object
    :param str variable: variable to act on
    :param array flags: flags array
    :param datetime start: data start
    :param datetime end: data end
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
    
    all_filtered = utils.apply_filter_flags(st_var)
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
                   
    return flags # dgc_monthly

#************************************************************************
def dgc_all_obs(station, variable, flags, start, end, plots = False, diagnostics = False, idl = False, windspeeds = False, GH = False):
    '''RJHD addition working on all observations'''
    
    if plots:
        import matplotlib.pyplot as plt

    st_var = getattr(station, variable)
    
    month_ranges = utils.month_starts_in_pairs(start, end)
    month_ranges = month_ranges.reshape(-1,12,2)
    
    all_filtered = utils.apply_filter_flags(st_var)

 
    for month in range(12):
    
        if windspeeds == True:
            st_var_wind = getattr(station, "windspeeds")
            
            # get monthly averages
            windspeeds_month = np.empty([])
            for y, year in enumerate(month_ranges[:,month,:]):
            
                if y == 0:
                    windspeeds_month = np.ma.array(st_var_wind.data[year[0]:year[1]])
                else:
                    windspeeds_month = np.ma.concatenate([windspeeds_month, st_var_wind.data[year[0]:year[1]]])
                  
            windspeeds_month_average = dgc_get_monthly_averages(windspeeds_month, OBS_LIMIT, st_var_wind.mdi, MEAN)
            windspeeds_month_mad = utils.mean_absolute_deviation(windspeeds_month, median=True)
    
        
        this_month_data = np.array([])
        this_month_filtered = np.array([])
        
        this_month_data, dummy, dummy = utils.concatenate_months(month_ranges[:,month,:], st_var.data, hours = False)
        this_month_filtered, dummy, dummy = utils.concatenate_months(month_ranges[:,month,:], all_filtered, hours = False)
                
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
    

            if iqr != 0.0:               
                monthly_values = np.ma.array((this_month_data.compressed() - monthly_median) / iqr)

                bins, bincenters = utils.create_bins(monthly_values, BIN_SIZE)
                dummy, plot_bincenters = utils.create_bins(monthly_values, BIN_SIZE/10.)
        
                hist, binEdges = np.histogram(monthly_values, bins = bins)
                                               
                if GH:
                    # Use Gauss-Hermite polynomials to add skew and kurtosis to Gaussian fit - January 2015 ^RJHD

                    initial_values = [np.max(hist), np.mean(monthly_values), np.std(monthly_values), stats.skew(monthly_values), stats.kurtosis(monthly_values)] # norm, mean, std, skew, kurtosis
                    
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
                                      

                else:
                    gaussian = utils.fit_gaussian(bincenters, hist, max(hist), mu = np.mean(monthly_values), sig = np.std(monthly_values))

                    # assume the same threshold value
                    u_minimum_threshold = 1 + round(utils.invert_gaussian(FREQUENCY_THRESHOLD, gaussian))
                    l_minimum_threshold = -u_minimum_threshold


                    plot_gaussian = utils.gaussian(plot_bincenters, gaussian)

                if diagnostics:
                    if GH:
                        print hist
                        print res
                        print iqr, l_minimum_threshold, u_minimum_threshold

                    else:
                        print hist
                        print gaussian
                        print iqr, u_minimum_threshold, 1. + utils.invert_gaussian(FREQUENCY_THRESHOLD, gaussian)

                if plots:
                    dgc_set_up_plot(plot_gaussian, monthly_values, variable, threshold = (u_minimum_threshold, l_minimum_threshold), sub_par = "observations", GH = GH)
                     
                    if GH:
                        plt.figtext(0.15, 0.67, 'Mean %.2f, S.d. %.2f,\nSkew %.2f, Kurtosis %.2f' %(res['mean'], res['dispersion'], res['skewness'], res['kurtosis']), color='k', size='small')

                    

                uppercount = len(np.where(monthly_values > u_minimum_threshold)[0])
                lowercount = len(np.where(monthly_values < l_minimum_threshold)[0])
                
                # this needs refactoring - but lots of variables to pass in
                if plots or diagnostics: gap_plot_values = np.array([])

                if uppercount > 0:
                    gap_start = dgc_find_gap(hist, binEdges, u_minimum_threshold)
                        
                    if gap_start != 0:
                        
                        for y, year in enumerate(month_ranges[:,month,:]):
                
                            this_year_data = np.ma.array(all_filtered[year[0]:year[1]])
                            this_year_flags = np.array(flags[year[0]:year[1]])
                            gap_cleaned_locations = np.where(((this_year_data - monthly_median) / iqr) > gap_start)

                            this_year_flags[gap_cleaned_locations] = 1
                            flags[year[0]:year[1]] = this_year_flags

                            if plots or diagnostics: gap_plot_values = np.append(gap_plot_values, (this_year_data[gap_cleaned_locations].compressed() - monthly_median)/iqr)


                if lowercount > 0:
                    gap_start = dgc_find_gap(hist, binEdges, l_minimum_threshold)
                        
                    if gap_start != 0:

                        for y, year in enumerate(month_ranges[:,month,:]):
                
                            this_year_data = np.ma.array(all_filtered[year[0]:year[1]])
                            this_year_flags = np.array(flags[year[0]:year[1]])
                            gap_cleaned_locations = np.where(np.logical_and(((this_year_data - monthly_median) / iqr) < gap_start, this_year_data.mask != True))

                            this_year_flags[gap_cleaned_locations] = 1
                            flags[year[0]:year[1]] = this_year_flags

                            if plots or diagnostics: gap_plot_values = np.append(gap_plot_values, (this_year_data[gap_cleaned_locations].compressed() - monthly_median)/iqr)
                    

                            if windspeeds:
                                this_year_flags[gap_cleaned_locations] = 2 # tentative flags
                                
                                slp_average = dgc_get_monthly_averages(this_month_data, OBS_LIMIT, st_var.mdi, MEAN)
                                slp_mad = utils.mean_absolute_deviation(this_month_data, median=True)
                                storms = np.where((((windspeeds_month - windspeeds_month_average) / windspeeds_month_mad) > 4.5) &\
                                                   (((this_month_data - slp_average) / slp_mad) > 4.5))
                                
                                if len(storms[0]) >= 2:
                                    
                                    storm_1diffs = np.diff(storms)
                                    
                                    separations = np.where(storm_1diffs != 1)

                                    #for sep in separations:


                if plots:
                    hist, binEdges = np.histogram(gap_plot_values, bins = bins)
                    plot_hist = np.array([0.01 if h == 0 else h for h in hist])
                    plt.step(bincenters, plot_hist, 'r-', label = 'flagged', where='mid')
                    import calendar
                    plt.text(0.1,0.9,calendar.month_name[month+1], transform = plt.gca().transAxes)
                    plt.legend(loc='lower center',ncol=3, bbox_to_anchor=(0.5,-0.2),frameon=False,prop={'size':13})
                    plt.show()
                    #plt.savefig(IMAGELOCATION+'/'+station.id+'_DistributionalGap_'+str(month+1)+'.png')
    if diagnostics:
        utils.print_flagged_obs_number("", "Distributional Gap", variable, len(gap_plot_values), noWrite=True)

    return flags # dgc_all_obs

#************************************************************************
def dgc(station, variable_list, flag_col, start, end, logfile, plots=False, diagnostics = False, idl = False, GH = False):
    '''Controller for two individual tests'''

    if plots:
        import matplotlib.pyplot as plt
        
    
    for v, variable in enumerate(variable_list):    
        station.qc_flags[:,flag_col[v]] = dgc_monthly(station, variable, station.qc_flags[:,flag_col[v]], start, end, plots=plots, diagnostics=diagnostics, idl = idl)
        
        if variable == "slp":
            # need to send in windspeeds too        
            station.qc_flags[:,flag_col[v]] = dgc_all_obs(station, variable, station.qc_flags[:,flag_col[v]], start, end, plots=plots, diagnostics=diagnostics, idl = idl, windspeeds = True, GH = GH)
        else:
            station.qc_flags[:,flag_col[v]] = dgc_all_obs(station, variable, station.qc_flags[:,flag_col[v]], start, end, plots=plots, diagnostics=diagnostics, idl = idl, GH = GH)
    
        flag_locs = np.where(station.qc_flags[:, flag_col[v]] != 0)

        if plots or diagnostics:
            utils.print_flagged_obs_number(logfile, "Distributional Gap", variable, len(flag_locs[0]), noWrite=True)
        else:
            utils.print_flagged_obs_number(logfile, "Distributional Gap", variable, len(flag_locs[0]))

        # copy flags into attribute
        st_var = getattr(station, variable)
        st_var.flags[flag_locs] = 1
 
        # MATCHES IDL for 030660-99999, 2 flags in T, 30-06-2014

    station = utils.append_history(station, "Distributional Gap Check")  
                     
    return # dgc

#************************************************************************

if __name__ == "__main__":

    print "running distributional gap check"
