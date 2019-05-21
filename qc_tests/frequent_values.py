#!/usr/local/sci/bin/python
#*****************************
#
# Frequent Value Check (FVC)
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

SEASONS = ['Ann','MAM','JJA','SON','D+JF']



#************************************************************************
def fvc_plot_setup(hist_data, hist, binEdges, xlabel, title = ""):
    '''
    Plot the histogram, with removed observations highlighted
    
    :param array hist_data: raw values which have been binned to create hist
    :param array hist: values of histogram
    :param array binEdges: location of LH bin edge
    :param str xlabel: label for x-axis
    :param str title: title of plot
    
    :returns:
        plot-hist - useful histogram data to plot in log-scale
        bincenters - locations of centres of bins
    '''
    import matplotlib.pyplot as plt

    plot_hist = np.array([float(x) if x != 0 else 1e-1 for x in hist])
    plt.clf()
    bincenters = 0.5 * (binEdges[1:] + binEdges[:-1])
    plt.step(bincenters, plot_hist, 'b-', label = 'observations', where='mid')
            
    fit = utils.fit_gaussian(bincenters, hist, max(hist), mu = np.mean(hist_data), sig = np.std(hist_data))
    plot_gaussian = utils.gaussian(bincenters, fit)
    plt.plot(bincenters, plot_gaussian, 'r-', label = 'Gaussian fit')
    # sort labels and prettify
    plt.xlabel(xlabel)                    
    plt.ylabel("Frequency")
    plt.gca().set_yscale('log')
    plt.ylim([0.1,10000])
    plt.title(title)

    return plot_hist, bincenters # fvc_plot_setup

#************************************************************************
def fvc(station, variable_list, flag_col, start, end, logfile, diagnostics = False, plots = False, doMonth = False):
    '''
    Check for certain values occurring more frequently than would be expected
    
    :param object station: station object to process
    :param list variable_list: list of variables to process
    :param list flag_col: columns to fill in flag array
    :param datetime start: datetime object of start of data
    :param datetime end: datetime object of end of data
    :param file logfile: logfile to store outputs
    :param bool diagnostics: produce extra diagnostic output
    :param bool plots: produce plots
    :param bool month: ignore months after last complete year/season for distribution
    '''
    
    MIN_DATA_REQUIRED = 500 # to create histogram for complete record
    MIN_DATA_REQUIRED_YEAR = 100 # to create histogram

    month_ranges = utils.month_starts_in_pairs(start, end)

    month_ranges_years = month_ranges.reshape(-1,12,2)

    for v,variable in enumerate(variable_list):
    
        st_var = getattr(station, variable)
        
        reporting_accuracy = utils.reporting_accuracy(utils.apply_filter_flags(st_var))
        
        # apply flags - for detection only
        filtered_data = utils.apply_filter_flags(st_var, doMonth = doMonth, start = start, end = end)

        for season in range(5): # Year,MAM,JJA,SON,JF+D
 
            if season == 0:
                # all year
                season_data = np.ma.masked_values(filtered_data.compressed(), st_var.fdi)
                thresholds = [30,20,10]

            else:
                thresholds = [20,15,10]
                season_data = np.ma.array([])
                
                for y,year in enumerate(month_ranges_years):
                    # churn through months extracting data, accounting for fdi and concatenating together
                    if season == 1:
                        #mam
                        season_data = np.ma.concatenate([season_data, np.ma.masked_values(filtered_data[year[2][0]:year[4][-1]], st_var.fdi)])
                    elif season == 2:
                        #jja
                        season_data = np.ma.concatenate([season_data, np.ma.masked_values(filtered_data[year[5][0]:year[7][-1]], st_var.fdi)])
                    elif season == 3:
                        #son
                        season_data = np.ma.concatenate([season_data, np.ma.masked_values(filtered_data[year[8][0]:year[10][-1]], st_var.fdi)])
                    elif season == 4:
                        #d+jf
                        season_data = np.ma.concatenate([season_data, np.ma.masked_values(filtered_data[year[0][0]:year[1][-1]], st_var.fdi)])
                        season_data = np.ma.concatenate([season_data, np.ma.masked_values(filtered_data[year[-1][0]:year[-1][-1]], st_var.fdi)])
                    

            season_data = season_data.compressed()

            if len(season_data) > MIN_DATA_REQUIRED:    

                if 0 < reporting_accuracy <= 0.5: # -1 used as missing value
                    bins, bincenters = utils.create_bins(season_data, 0.5)
                else:
                    bins, bincenters = utils.create_bins(season_data, 1.0)

                hist, binEdges = np.histogram(season_data, bins = bins)
            
                if plots:
                    plot_hist, bincenters = fvc_plot_setup(season_data, hist, binEdges, st_var.name, title = "%s" % (SEASONS[season]))

                bad_bin = np.zeros(len(hist))

                # scan through bin values and identify bad ones
                for e, element in enumerate(hist):                  
                    if e > 3 and e <= (len(hist) - 3):
                        # don't bother with first three or last three bins
                        seven_bins = hist[e-3:e+3+1]
                        if (seven_bins[3] == seven_bins.max()) and (seven_bins[3] != 0):
                            # is local maximum and != zero
                            if (seven_bins[3]/float(seven_bins.sum()) >= 0.5) and (seven_bins[3] >= thresholds[0]):
                                # contains >50% of data and is greater than threshold
                                bad_bin[e] = 1

                            # for plotting remove good bins
                            else:
                                if plots: plot_hist[e]=1e-1
                        else:
                            if plots: plot_hist[e]=1e-1
                    else:
                        if plots: plot_hist[e]=1e-1



                if plots:
                    import matplotlib.pyplot as plt
                    plt.step(bincenters, plot_hist, 'r-', where='mid')
                    plt.show()
            
                # having identified possible bad bins, check each year in turn, on unfiltered data
                for y,year in enumerate(month_ranges_years):

                    if season == 0:
                        # year
                        year_data = np.ma.masked_values(st_var.data[year[0][0]:year[-1][-1]], st_var.fdi)
                        year_flags = station.qc_flags[year[0][0]:year[-1][-1],flag_col[v]]
                    elif season == 1:
                        #mam
                        year_data = np.ma.masked_values(st_var.data[year[2][0]:year[4][-1]], st_var.fdi)
                        year_flags = station.qc_flags[year[2][0]:year[4][-1],flag_col[v]]
                    elif season == 2:
                        #jja
                        year_data = np.ma.masked_values(st_var.data[year[5][0]:year[7][-1]], st_var.fdi)
                        year_flags = station.qc_flags[year[5][0]:year[7][-1],flag_col[v]]
                    elif season == 3:
                        #son
                        year_data = np.ma.masked_values(st_var.data[year[8][0]:year[10][-1]], st_var.fdi)
                        year_flags = station.qc_flags[year[8][0]:year[10][-1],flag_col[v]]
                    elif season == 4:
                        #d+jf
                        year_data = np.ma.concatenate([np.ma.masked_values(st_var.data[year[0][0]:year[1][-1]], st_var.fdi),\
                                                       np.ma.masked_values(st_var.data[year[-1][0]:year[-1][-1]], st_var.fdi)])
                        year_flags = np.append(station.qc_flags[year[0][0]:year[1][-1],flag_col[v]],station.qc_flags[year[-1][0]:year[-1][-1],flag_col[v]])


                    if len(year_data.compressed()) > MIN_DATA_REQUIRED_YEAR:    

                        hist, binEdges = np.histogram(year_data.compressed(), bins = bins)

                        if plots:
                            plot_hist, bincenters = fvc_plot_setup(year_data.compressed(), hist, binEdges, st_var.name, title = "%s - %s" % (y+start.year, SEASONS[season]))

                        for e, element in enumerate(hist):

                            if bad_bin[e] == 1:
                                # only look at pre-identified bins

                                if e >= 3 and e <= (len(hist) - 3):
                                    # don't bother with first three or last three bins
                                    seven_bins = hist[e-3:e+3+1].astype('float')
                                    if (seven_bins[3] == seven_bins.max()) and (seven_bins[3] != 0):
                                        # is local maximum and != zero
                                        if (seven_bins[3]/seven_bins.sum() >= 0.5 and seven_bins[3] >= thresholds[1]) \
                                            or (seven_bins[3]/seven_bins.sum() >= 0.9 and seven_bins[3] >= thresholds[2]):
                                            # contains >50% or >90% of data and is greater than appropriate threshold

                                            # Flag these data
                                            bad_points = np.where((year_data >= binEdges[e]) & (year_data < binEdges[e+1]))
                                            year_flags[bad_points] = 1

                                        # for plotting remove good bins
                                        else:
                                            if plots: plot_hist[e]=1e-1
                                    else:
                                        if plots: plot_hist[e]=1e-1
                                else:
                                    if plots: plot_hist[e]=1e-1
                            else:
                                if plots: plot_hist[e]=1e-1

                        if diagnostics or plots:
                            nflags = len(np.where(year_flags != 0)[0])
                            print "{} {}".format(y + start.year, nflags)

                        if plots:
                            if nflags > 0:
                                plt.step(bincenters, plot_hist, 'r-', where='mid')
                                plt.show()
                            else:
                                plt.clf()

                    # copy flags back

                    if season == 0:
                        station.qc_flags[year[0][0]:year[-1][-1], flag_col[v]] = year_flags   
                    elif season == 1:
                        station.qc_flags[year[2][0]:year[4][-1], flag_col[v]] = year_flags   
                    elif season == 2:
                        station.qc_flags[year[5][0]:year[7][-1], flag_col[v]] = year_flags   
                    elif season == 3:
                        station.qc_flags[year[8][0]:year[10][-1], flag_col[v]] = year_flags   
                    elif season == 4:
                        split = len(station.qc_flags[year[0][0]:year[1][-1], flag_col[v]])
                        station.qc_flags[year[0][0]:year[1][-1], flag_col[v]] = year_flags[:split]
                        station.qc_flags[year[-1][0]:year[-1][-1], flag_col[v]] = year_flags[split:]
 
        flag_locs = np.where(station.qc_flags[:, flag_col[v]] != 0)
        utils.print_flagged_obs_number(logfile, "Frequent Value", variable, len(flag_locs[0]), noWrite = diagnostics)

        # copy flags into attribute
        st_var.flags[flag_locs] = 1
        
    station = utils.append_history(station, "Frequent Values Check")  
                     
    return # fvc


#************************************************************************
if __name__ == "__main__":

    print "Checking for Frequent Values"
