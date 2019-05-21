#!/usr/local/sci/bin/python
#*****************************
#
# Detect outliers from station and neighbour
#
#
#************************************************************************
#                    SVN Info
#$Rev:: 129                                           $:  Revision of last commit
#$Author:: rdunn                                      $:  Author of last commit
#$Date:: 2017-08-15 10:02:07 +0100 (Tue, 15 Aug 2017) $:  Date of last commit
#************************************************************************


import numpy as np
import scipy as sp
import math
import os
import gc

# RJHD routines
import qc_utils as utils
import netcdf_procs as ncdfp
from set_paths_and_vars import *


N_NEIGHBOURS = 10

UNFLAG_COL_DICT = {"spike":{"temperatures":27 ,"dewpoints":28, "slp":29},\
                       "climatological":{"temperatures":24 ,"dewpoints":25, "slp":26},\
                       "odd":{"temperatures":54 ,"dewpoints":55, "slp":57},\
                       "gap":{"slp":7},\
                       "dpd":{"dewpoints":31}}

#*******************************************************
def get_neighbour_parameters(station_info):
    '''
    Return two big symmetrical arrays of the station separations and bearings
    '''

    distances = np.zeros([len(station_info),len(station_info)])
    angles = np.zeros([len(station_info),len(station_info)])
    
    all_lats = np.array(station_info[:,1], dtype=float)
    all_lons = np.array(station_info[:,2], dtype=float)


    for st, stat in enumerate(station_info):
        distances[st,:], angles[st,:] = utils.get_dist_and_bearing([all_lats[st],all_lons[st]],[all_lats,all_lons])
        distances[st,st] = 0. # by hand set the diagonal to zero
        angles[st,st] = 0. # by hand set the diagonal to zero

    return distances, angles # get_neighbour_parameters

#*******************************************************
def get_neighbours(station_loc, st_elev, distances, bearings, elevations, sep_limit = 500., elev_limit = 300.):
    '''
    Get the neighbours for the station using distance, angles and elevations
    Will return 10 or less neighbours.

    :param int station_loc: sequence number of station
    :param float st_elev: station elevation
    :param array distances: distances of all other stations
    :param array bearings: bearings of all other stations
    :param array elevations: elevations of all other stations
    :param float sep_limit: separation limit (default 500km)
    :param float elev_limit: elevation limit (default 300m)

    :returns: final set of neighbours as numbers in station list
    '''


    # get bearings into quadrants to make searching easier
    quadrants = np.zeros(bearings.shape)
    for a,angle in enumerate([270,180,90,0]):
        locs = np.where((bearings >= angle) & (bearings < angle+90.))
        quadrants[locs] = a+1

    # start by sorting by distance
    ordering = np.argsort(distances)

    neighbour_locations = np.array([], dtype = int)
    neighbour_quadrants = np.array([])
    for index in ordering:

        if index == station_loc:
            # don't choose itself!
            continue

        # if within range of distance and elevation
        if distances[index] <= sep_limit:
            if np.abs(elevations[index] - st_elev) <= elev_limit:

                neighbour_locations=np.append(neighbour_locations, index)
                neighbour_quadrants = np.append(neighbour_quadrants, quadrants[index]) # combine into one data structure
                               
                # if have more than 10, see if angularly spread
                if len(neighbour_locations) >= N_NEIGHBOURS:
                    
                    locs1 = neighbour_locations[neighbour_quadrants == 1]
                    locs2 = neighbour_locations[neighbour_quadrants == 2]
                    locs3 = neighbour_locations[neighbour_quadrants == 3]
                    locs4 = neighbour_locations[neighbour_quadrants == 4]

                    # if at least two per quadrant then stop looking
                    if len(locs1) >= 2 and len(locs2) >= 2 and len(locs3) >= 2 and len(locs4) >= 2:
                        break


    # if have more than 10
    if len(neighbour_locations) >= N_NEIGHBOURS:

        # now go through to ensure have only 10, and two per quadrant
        final_locs = np.concatenate((locs1[:2], locs2[:2], locs3[:2], locs4[:2]), axis = 0).reshape(-1)

        # and add the rest in order of distance
        for index in neighbour_locations:
            if index not in final_locs:
                final_locs = np.append(final_locs, index)

            if len(final_locs) == N_NEIGHBOURS:
                break

        return final_locs
    # else just return the ones found
    else:
        return neighbour_locations # get_neighbours

#*******************************************************
def get_all_neighbours(station_loc, st_elev, distances, bearings, elevations, sep_limit = 500., elev_limit = 300., max_neighbours = 20):
    '''
    Get the neighbours for the station using distance, angles and elevations
    Returns all neighbours within distance and elevation ranges 

    :param int station_loc: sequence number of station
    :param float st_elev: station elevation
    :param array distances: distances of all other stations
    :param array bearings: bearings of all other stations
    :param array elevations: elevations of all other stations
    :param float sep_limit: separation limit (default 500km)
    :param float elev_limit: elevation limit (default 300m)
    :param int max_neighbours: maximum number to return (default 50)

    :returns: final set of neighbours as numbers in station list
    '''


    # get bearings into quadrants to make searching easier
    quadrants = np.zeros(bearings.shape)
    for a,angle in enumerate([270,180,90,0]):
        locs = np.where((bearings >= angle) & (bearings < angle+90.))
        quadrants[locs] = a+1

    # start by sorting by distance
    ordering = np.argsort(distances)

    neighbour_locations = np.array([], dtype = int)
    neighbour_quadrants = np.array([])
    for index in ordering:

        if index == station_loc:
            # don't choose itself!
            continue

        # if within range of distance and elevation
        if distances[index] <= sep_limit:
            if np.abs(elevations[index] - st_elev) <= elev_limit:

                neighbour_locations=np.append(neighbour_locations, index)
                neighbour_quadrants = np.append(neighbour_quadrants, quadrants[index]) # combine into one data structure
                               

    return neighbour_locations[:max_neighbours], neighbour_quadrants[:max_neighbours] # get_all_neighbours

#*******************************************************
def hourly_daily_anomalies(timeseries, obs_per_day = 6):
    '''
    Format the time series to get allow for sensible correlations

     - Process into 24h x N_days
     - Obtain daily average and hence hourly anomalies from daily average (removes annual cycle)
     - Obtain hourly average of anomalies to get double-anomalies (removes diurnal cycle)

    :param array timeseries: data to be processed in 1-D array
    :param int obs_per_day: number of observations per 24hr to get a daily average to process (default = 6)

    :returns: anomalies - timeseries with removed annual and diurnal cycles.
    '''

    timeseries = timeseries.reshape((-1,24))

    daily_mean = np.ma.mean(timeseries, axis = 1)

    not_mask_count = timeseries.count(axis = 1)

    daily_mean = np.ma.masked_where(not_mask_count < obs_per_day, daily_mean) # completeness check - only for correlations

    hourly_anomalies = timeseries - np.ma.repeat(daily_mean, 24).reshape((-1,24)) # removed annual cycle

    hourly_mean = np.ma.mean(hourly_anomalies, axis = 0)

    anomalies = hourly_anomalies - np.tile(hourly_mean, timeseries.shape[0]).reshape((-1,24)) # removed diurnal cycle

    return anomalies.ravel() # hourly_daily_anomalies

#*******************************************************
def select_neighbours(station, variable, neighbour_info, neighbours, neighbour_distances, neighbour_quadrants, data_locs, datastart, dataend, logfile, diagnostics = False, plots = False):
    '''
    From the list of nearby stations select the ones which will be good neighours for the test.
    Select on basis of correlation, overlap of data points and bearing (quadrants)
    
    :param object station: station object
    :param str variable: which variable to proces
    :param array neighbour_info: array of ID, lat, lon and elev
    :param array neighbours: which station sequence numbers are the nearby stations
    :param array neighbour_distances: distances to nearby stations
    :param array neighbour_quadrants: bearings to nearby stations (in 90deg bins)
    :param array data_locs: path to data files
    :param datetime datastart: start of data set
    :param datetime dataend: end of data set
    :param file logfile: logfile to store outputs
    :param boolean diagnostics: output diagnostic information
    :param boolean plots: make a plot

    :returns: final_locs - array of station sequence numbers to use.
    '''

    # set up storage arrays
    n_correlations = np.zeros(len(neighbours))
    n_distances = np.zeros(len(neighbours))
    n_quadrants = np.zeros(len(neighbours))
    n_overlaps = np.zeros(len(neighbours))
    combined_score = np.zeros(len(neighbours))

    # get station data
    st_var = getattr(station, variable)
    st_anomalies = hourly_daily_anomalies(st_var.data[:])

    # go through initial list and extract correlations and overlaps
    for nn, nn_loc in enumerate(neighbours):

        n_details = neighbour_info[nn]
        neigh = utils.Station(n_details[0], float(n_details[1]), float(n_details[2]), float(n_details[3]))

        ncdfp.read(os.path.join(NETCDF_DATA_LOCS, "hadisd.{}_19310101-{}_{}_internal.nc".format(LONG_VERSION, END_TIME, station.id)), neigh, [variable], diagnostics = diagnostics, read_input_station_id = False)
       

        dummy = utils.create_fulltimes(neigh, [variable], datastart, dataend, [], do_input_station_id = False)

        # get the correlations of data to this neighbour
        neigh_var = getattr(neigh, variable)
        neigh_anomalies = hourly_daily_anomalies(neigh_var.data[:])
        # correlation = np.ma.corrcoef(neigh_var.data, st_var.data)[1,0]
        correlation = np.ma.corrcoef(neigh_anomalies, st_anomalies)[1,0]

        overlap = len(np.where(np.logical_or(neigh_var.data.mask, st_var.data.mask) == False)[0])/float(len(st_var.data.compressed()))

        if not math.isnan(correlation):
            n_correlations[nn] = correlation
            n_overlaps[nn] = overlap
            combined_score[nn] = correlation + overlap
            n_distances[nn] = neighbour_distances[nn]
            n_quadrants[nn] = neighbour_quadrants[nn]
            
        # clear up to save memory
        del dummy
        del neigh_var
        del neigh_anomalies
        gc.collect()
    # sort in order of the combination of correlation and overlap
    sort_order = np.argsort(combined_score)[::-1]

    # and select the best 10
    # final_selection = neighbours[sort_order][:10]

    # sort out the quadrants
                    
    locs1 = neighbours[sort_order][n_quadrants[sort_order] == 1]
    locs2 = neighbours[sort_order][n_quadrants[sort_order] == 2]
    locs3 = neighbours[sort_order][n_quadrants[sort_order] == 3]
    locs4 = neighbours[sort_order][n_quadrants[sort_order] == 4]

    final_locs = np.concatenate((locs1[:2], locs2[:2], locs3[:2], locs4[:2]), axis = 0).reshape(-1)

    # and add the rest in order of combined score
    for index in neighbours[sort_order]:
        if index not in final_locs:
            final_locs = np.append(final_locs, index)
            
        if len(final_locs) == N_NEIGHBOURS:
            break

    # output table showing distances, correlations, overlaps, the combined score and which ones were selected
    if plots or diagnostics:
        print "{:14s} {:10s} {:10s} {:10s} {:10s} {:10s} {:10s} {:10s}".format("Neighbour","Distance","Elevation", "Correl'n", "Overlap", "Combined", "Quadrant","Selected")
    else:
        logfile.write("{:14s} {:10s} {:10s} {:10s} {:10s} {:10s} {:10s} {:10s}\n".format("Neighbour","Distance","Elevation", "Correl'n", "Overlap", "Combined", "Quadrant","Selected")
)

    selected_correlations = []
    selected_overlaps = []
    for nn, nn_loc in enumerate(neighbours[sort_order]):

        selected = ""
        if nn_loc in final_locs: 
            selected = "Y"
            if plots:
                selected_correlations += [n_correlations[sort_order[nn]]]
                selected_overlaps += [n_overlaps[sort_order[nn]]]


        neigh_details = neighbour_info[sort_order][nn]
        if plots or diagnostics:
            print "{:14s} {:10.1f} {:10.1f} {:10.5f} {:10.3f} {:10.3f} {:10.0f} {:10s}".format(neigh_details[0], n_distances[sort_order][nn], float(neigh_details[3]), n_correlations[sort_order][nn], n_overlaps[sort_order][nn], combined_score[sort_order][nn], n_quadrants[sort_order][nn], selected)
        else:
            logfile.write("{:14s} {:10.1f} {:10.1f} {:10.5f} {:10.3f} {:10.3f} {:10.0f} {:10s}\n".format(neigh_details[0], n_distances[sort_order][nn], float(neigh_details[3]), n_correlations[sort_order][nn], n_overlaps[sort_order][nn], combined_score[sort_order][nn], n_quadrants[sort_order][nn], selected))

            
    # plot of correlations and overlaps, with selected stations highlighted
    if plots:
        import matplotlib.pyplot as plt

        plt.clf()
        plt.plot(n_correlations, n_overlaps, 'bo')
        plt.plot(selected_correlations, selected_overlaps, 'ro')
        plt.xlabel("correlations")
        plt.ylabel("data overlap")
        plt.title("{} - {}".format(station.id, variable))
        plt.show()

    return final_locs # select_neighbours

#*******************************************************
def plot_target_neigh_diffs_dist(differences, iqr):
    '''
    Plot the distribution of target-neighbour differences
    
    :param array differences: masked difference array
    :param float iqr: inter quartile range of differences

    :returns: 
    '''
    import matplotlib.pyplot as plt
    
    plt.clf()
    
    bins, bincenters = utils.create_bins(differences.compressed(), 1.0)
    
    hist, binEdges = np.histogram(differences.compressed(), bins=bins)
    plot_hist = np.array([float(x) if x != 0 else 1e-1 for x in hist])
    plt.step(bincenters, plot_hist, 'k-', label = 'observations', where='mid')
    
    fit = utils.fit_gaussian(bincenters, hist, max(hist), mu=np.mean(differences.compressed()), sig = np.std(differences.compressed()))
    plot_gaussian = utils.gaussian(bincenters, fit)
    plt.plot(bincenters, plot_gaussian, 'b-', label = 'Gaussian fit')
    
    plt.axvline(5.*iqr, c = 'r')
    plt.axvline(-5.*iqr, c = 'r')
    
    print "only shows lowest of monthly IQRs"

    plt.ylabel("Frequency")
    plt.gca().set_yscale('log')
    plt.ylim([0.1,2*max(hist)])
    
    plt.show()      

    return # plot_target_neigh_diffs_dist


#*******************************************************
def detect(station, neighbour, variable, flags, neighbour_count, start, end, distance = 0, diagnostics = False, plots = False):
    '''
    Detect which observations are outliers

    :param MetVar station: station object (target)
    :param MetVar neighbour: station object (neighbour)
    :param string variable: which variable to process
    :param array flags: array to store how many neighbours thing each obs is bad
    :param array neighbour_count: how many neighbours present at each obs
    :param datetime start: start of dataset
    :param datetime end: end of dataset
    :param int distance: separation of target and neighbour
    :param bool diagnostics: extra output
    :param bool plots: make figures

    :returns: None
    '''

    FILTERING_FLAG_COL = {"temperatures":[0,1,4,5,8,12,16,20,27,41,44,58],
                          "dewpoints":[0,2,4,6,8,9,13,17,21,28,30,31,32,42,45,59],
                          "slp":[0,3,4,7,11,15,19,23,29,43,46,60],
                          "windspeeds":[0,4,10,14,18,22,56,62,63,64]} # not used, but ready for it.
    
    
    st_var = getattr(station, variable)
    neigh_var = getattr(neighbour, variable)

    # filter by flags - not all (no Climatological [24,25], or Odd cluster [54,55,56,57]), T record check not in D, 
    total_flags = np.sum(station.qc_flags[:,FILTERING_FLAG_COL[variable]], axis = 1)
    st_filtered = np.ma.masked_where(total_flags == 1, st_var.data)
    neigh_filtered = np.ma.masked_where(total_flags == 1, neigh_var.data)

    # match the observation times
    match = np.where(np.logical_and((st_filtered.data != st_var.mdi), (neigh_filtered.data != neigh_var.mdi)))

    month_ranges = utils.month_starts_in_pairs(start, end).reshape(-1,12,2) # in year-long sets of pairs.        

    if len(match[0]) >= 100:

        neighbour_count[match] += 1 # number of neighbours with data present

        differences = np.ma.zeros(len(st_filtered))
        differences.fill(st_var.mdi)
        differences.mask = True

        differences[match] = st_filtered.data[match] - neigh_filtered.data[match]
        differences.mask[match] = False

        all_iqrs = np.zeros(len(differences))
        # get monthly IQR values
        for month in range(12):
  
            this_month, dummy1, dummy2 = utils.concatenate_months(month_ranges[:,month,:], differences, hours = False)

            if len(this_month.compressed()) > 4:

                iqr = utils.IQR(this_month.compressed())
                if iqr <= 2.: iqr = 2.
            else:
                iqr = 2.

            # and copy back into the array
            for year in month_ranges[:,month,:]:
                all_iqrs[year[0]:year[1]] = iqr

        if plots:
            plot_target_neigh_diffs_dist(differences, min(all_iqrs))

        dubious = np.ma.where(np.ma.abs(differences) > 5. * all_iqrs)

        if len(dubious[0]) >= 1.:

            if variable == "slp":
                # check if they are storms
                positive = np.ma.where(differences > 5. * iqr)
                negative = np.ma.where(differences < -5. * iqr)
                
                # if majority negative (2/3) and separation > 100

                if (distance > 100.) and (float(len(positive[0]))/len(dubious[0]) < 0.333):

                    if len(positive[0]) > 0:
                        flags[positive] += 1
                    if len(negative[0]) > 0:
                        neighbour_count[match] -= 1
                    
                else:
                    flags[dubious] += 1
            else:

                flags[dubious] += 1

    return # detect

#*******************************************************
def plot_outlier(station, variable, outlier_locs, all_data, datastart):
    '''
    Plot the outlier location (either to flag or unflag) with the target and all neighbours

    :param MetVar station: station object
    :param str variable: variable to process
    :param array outlier_locs: locations which are outliers to plot
    :param array all_data: all data from neighbours for plotting
    :param datetime datastart: start of dataset

    :returns: None
    '''

    import matplotlib.pyplot as plt
    import datetime as dt
    
    YLABELS = {"temperatures":"Temperature (C)", "dewpoints":"Dewpoints (C)", "slp":"SLP (hPa)", "windspeeds":"Wind Speed (m/s)"}
    extra = 48 # hours

    indata = getattr(station, variable).data
    
    for location in outlier_locs:
        
        plot_times = utils.times_hours_to_datetime(station.time.data[location-extra: location+extra], datastart)
        
        plt.clf()
        plt.plot(plot_times, indata[location-extra:location+extra], 'bo', ls='-')
        plt.plot(plot_times[extra], indata[location], 'ro', markersize=10)
        
        for nn in range(all_data.shape[0]):
            plt.plot(plot_times,all_data[nn,location-extra:location+extra] , c='0.5', ls='-')
            
        plt.ylabel(YLABELS[variable])
        plt.title("{:s} {:s}".format(station.id,dt.datetime.strftime(plot_times[extra], "%d/%m/%Y")))
        plt.show()
                    
    return # plot_outlier

#*******************************************************
def unflagging_locs(differences, flags, neigh_count, dpd_count = [], flag_value = 1):
    '''
    Return locations where flags to be set to zero
    Currently has deliberate bug, but this is to match IDL

    :param array differences: normalised differences of target - median of neighbours
    :param array flags: flag array to check where to unset
    :param array neigh_count: number of neighbour obs at each timestamp
    :param list dpd_count: number of DPD flags set at each timestamp
    :param int flag_value: can use to look for tentative flags

    :returns:
       unset_locs - list of locations to unset
    '''

    sufficient_neighbours, = np.where(neigh_count >= 3)

    if len(sufficient_neighbours) > 0:

        flag_locs, = np.where(flags[sufficient_neighbours] == flag_value)

        if len(flag_locs) > 0:

            if dpd_count != []:
                # if at least 2/3 of neighbours also have DPD flag set, then unset these
                dpd_proportion = dpd_count/neigh_count
                
                unset_locs = np.where(dpd_proportion[sufficient_neighbours][flag_locs] >= 2./3.)
                
            else:
                
                unset_locs = np.where(differences[sufficient_neighbours][flag_locs] <= 4.5)
                
            if len(unset_locs) > 0:
                # have some values that fall within 4.5MAD of median or where DPD is
                #  present in sufficiently many neighbours, remove these flags

                return sufficient_neighbours[flag_locs][unset_locs]
            

    return [] # unflagging

#*******************************************************
def bn_median(masked_array, axis=None):
    """
    https://github.com/astropy/ccdproc/blob/122cdbd5713140174f057eaa8fdb6f9ce03312df/docs/ccdproc/bottleneck_example.rst
    Perform fast median on masked array

    Parameters

    masked_array : `numpy.ma.masked_array`
        Array of which to find the median.

    axis : int, optional
        Axis along which to perform the median. Default is to find the median of
        the flattened array.
    """
    import numpy as np
    import bottleneck as bn
    data = masked_array.filled(fill_value=np.NaN)
    med = bn.nanmedian(data, axis=axis)
    # construct a masked array result, setting the mask from any NaN entries
    return np.ma.array(med, mask=np.isnan(med))

#*******************************************************
def median_absolute_deviation(a, axis=None):
    """from astropy.stats, converted to use np.ma - 13-Oct-2014 RJHD

    Compute the median absolute deviation.

    Returns the median absolute deviation (MAD) of the array elements.
    The MAD is defined as ``median(abs(a - median(a)))``.

    Parameters
    ----------
    a : array_like
        Input array or object that can be converted to an array.
    axis : int, optional
        Axis along which the medians are computed. The default (axis=None)
        is to compute the median along a flattened version of the array.

    Returns
    -------
    median_absolute_deviation : ndarray
        A new array holding the result. If the input contains
        integers, or floats of smaller precision than 64, then the output
        data-type is float64.  Otherwise, the output data-type is the same
        as that of the input.

    Examples
    --------

    This will generate random variates from a Gaussian distribution and return
    the median absolute deviation for that distribution::

        >>> from astropy.stats import median_absolute_deviation
        >>> from numpy.random import randn
        >>> randvar = randn(10000)
        >>> mad = median_absolute_deviation(randvar)

    See Also
    --------
    numpy.median

    """

    a = np.ma.array(a, copy=False)
    a_median = bn_median(a, axis=axis)

    # re-broadcast the output median array to subtract it
    if axis is not None:
        a_median = np.ma.expand_dims(a_median, axis=axis)

    # calculated the median average deviation
    return bn_median(np.ma.abs(a - a_median), axis=axis) # median_absolute_deviation

#*******************************************************
def do_unflagging(station, variable, all_data, reporting_accuracies, neigh_count, dpd_flags, FLAG_COL_DICT, start, logfile, plots = False, diagnostics = False):
    '''
    Set up and run the unflagging process for the specified tests

    :param MetVar station: station object
    :param string variable: variable to process
    :param array all_data: array containing all neighbour obs for full time period
    :param array reporting accuracies: reporting accuracy for each neighbour
    :param array neigh_count: number of neighbours with data at each time stamp
    :param array dpd_flags: number of neighbours that have DPD set at each time stamp
    :param dict FLAG_COL_DICT: look up dictionary to 
    :param datetime start: start of dataset
    :param file logfile: logfile to store outputs
    :param bool plots: do plots
    '''


    # unflagging using neighbours
    '''This is slow - np.ma.median is known to be slow
    https://github.com/astropy/ccdproc/issues/74
    https://github.com/astropy/ccdproc/blob/122cdbd5713140174f057eaa8fdb6f9ce03312df/docs/ccdproc/bottleneck_example.rst'''
    mean_of_neighbours = bn_median(all_data, axis = 0)
    std_of_neighbours = median_absolute_deviation(all_data, axis = 0)

    # find where the spread of neighbour observations is less than 1/2
    #    of maximum reporting accuracy
    std_of_neighbours[std_of_neighbours < 0.5*max(reporting_accuracies)] = 0.5*max(reporting_accuracies)

    # create series of normalised differences of obs from neighbour mean
    st_var = getattr(station, variable)
    normalised_differences = np.ma.abs(st_var.data - mean_of_neighbours)/std_of_neighbours
            
    for qc_test in ["climatological","gap","odd","dpd"]:
        
        if qc_test == "dpd" and variable == "dewpoints":
            flags = station.qc_flags[:, UNFLAG_COL_DICT[qc_test][variable]] 
            unset_locs = unflagging_locs(normalised_differences, flags, neigh_count, dpd_count = dpd_flags)
        elif qc_test == "dpd":
            # only unflag DPD on dewpoints
            continue

        elif qc_test == "gap" and variable != "slp":
            # only unflag gap check on slp observations
            continue

        else:
            flags = station.qc_flags[:, UNFLAG_COL_DICT[qc_test][variable]] 
            if qc_test == "gap" or qc_test == "climatological":
                # only tentative flags
                unset_locs = unflagging_locs(normalised_differences, flags, neigh_count, flag_value = 2)
            else:
                unset_locs = unflagging_locs(normalised_differences, flags, neigh_count)

        if len(unset_locs) > 0:
            station.qc_flags[unset_locs, UNFLAG_COL_DICT[qc_test][variable]] = 0
            
            # need to unflag attribute if and only if no other flags are set
            subset_flags = station.qc_flags[:, FLAG_COL_DICT[variable]]
            total_flags = np.sum(subset_flags[unset_locs, :], axis = 1)
            clean_locs = np.where(total_flags == 0)
            st_var.flags[unset_locs[clean_locs]] = 0

        # and print result
        if plots or diagnostics:
            utils.print_flagged_obs_number(logfile, "Unflagging "+qc_test, variable, len(unset_locs), noWrite = True)
        else:
            utils.print_flagged_obs_number(logfile, "Unflagging "+qc_test, variable, len(unset_locs))
            

        if plots:
            if len(unset_locs) > 0:
                plot_outlier(station, variable, unset_locs, all_data, start)

    station = utils.append_history(station, "Unflagging - "+variable)  
    
    return # do_unflagging


#************************************************************************
if __name__ == "__main__":
    
    print "Utilities for neighbour check"
