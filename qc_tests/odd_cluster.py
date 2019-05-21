#!/usr/local/sci/bin/python
#*****************************
#
# QC checks.
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


#*********************************************
class OddCluster:
    '''
    Class for odd cluster information
    '''
    
    def __init__(self, start, end, length, locations, data_mdi, last_data):
        self.start = start
        self.end = end
        self.length = length
        self.locations = locations
        self.data_mdi = data_mdi
        self.last_data = last_data
        
    def __str__(self):
        return "odd cluster, starting {}, ending {}, length {}".format(self.start, self.end, self.length)
    
    __repr__ = __str__
  
#*********************************************
def oc_plots(station, cluster, time, start, indata, variable, oc_details):
    '''
    Plot each odd cluster highlighted against surrounding data

    :param MetVar station: station object
    :param OddCluster cluster: cluster object
    :param int time: timestamp
    :param datetime start: start of dataseries
    :param masked array indata: input data
    :param string variable: variable name
    
    :returns:
    '''


    import matplotlib.pyplot as plt
    YLABELS = {"temperatures":"Temperature (C)", "dewpoints":"Dewpoints (C)", "slp":"SLP (hPa)", "windspeeds":"Wind Speed (m/s)"}
    
    plot_start, plot_end = cluster.locations[0] - 10*24 , time + 10*24
    if plot_start < 0 : plot_start = 0
    
    plot_times = utils.times_hours_to_datetime(station.time.data[plot_start: plot_end], start)
    
    plt.clf()
    plt.plot(plot_times, indata[plot_start: plot_end], 'bo')
    plt.plot(plot_times[np.array(oc_details.locations) - plot_start], indata[oc_details.locations], 'ro')

    plt.ylim(utils.sort_ts_ylim(indata[plot_start: plot_end]))
    plt.ylabel(YLABELS[variable])
    plt.show()

    return # oc_plots

#*********************************************
def occ_normal(cluster, obs_type, time, flags):
    '''
    Just a normal observation. (Not all inputs used here, 
       but required for consistency)

    :param OddCluster cluster: cluster object
    :param int obs_type: type determinator of observation
    :param int time: timestamp
    :param array flags: flag array 
    
    :returns:
       cluster - updated cluster information
       obs_type - updated observation type
    '''
    
    cluster.last_data = time

    return cluster, obs_type # occ_normal


#*********************************************
def occ_in_cluster(cluster, obs_type, time, flags):
    '''currently in a potential cluster.  (Not all inputs used here, 
       but required for consistency)

    :param OddCluster cluster: cluster object
    :param int obs_type: type determinator of observation
    :param int time: timestamp
    :param array flags: flag array 
    
    :returns:
       cluster - updated cluster information
       obs_type - updated observation type

    '''
    if (cluster.length == 6) or (time - cluster.start > 24.):
        '''longer than 6 hours or span over 24hrs --> not a cluster --> reset'''
        
        mdi = cluster.data_mdi
        # set all specifically
        cluster.start = mdi
        cluster.end   = mdi
        cluster.length = 0
        cluster.locations = mdi
        cluster.last_data = mdi

        obs_type = 0

    else:
        '''in a cluster, less than 6hr and not over 24hr - increment '''
        cluster.length += 1
        cluster.locations += [time]
        cluster.end = time

    return cluster, obs_type # occ_in_cluster

#*********************************************
def occ_start_cluster(cluster, obs_type, time, flags):
    '''
    There has been a gap in the data, check if long enough to start cluster
       (Not all inputs used here, but required for consistency)

    :param OddCluster cluster: cluster object
    :param int obs_type: type determinator of observation
    :param int time: timestamp
    :param array flags: flag array 
    
    :returns:
       cluster - updated cluster information
       obs_type - updated observation type
    '''

    if time - cluster.last_data >= 48:
        '''If gap of 48hr, then start cluster increments '''
        
        obs_type = 2
        cluster.length += 1
        cluster.start = time
        cluster.end   = time
        cluster.locations = [time]
        
    else:
        '''Gap in data not sufficiently large '''
        obs_type = 1
        cluster.last_data = time
        cluster.start = cluster.data_mdi
        cluster.end   = cluster.data_mdi
        cluster.length = 0
        cluster.locations = 0


    return cluster, obs_type # occ_start_cluster


#*********************************************
def occ_after_cluster(cluster, obs_type, time, flags):
    '''
    There has been a gap in the data after a cluster;
      check if long enough to mean cluster is sufficiently isolated.
      If so, flag else reset

    :param OddCluster cluster: cluster object
    :param int obs_type: type determinator of observation
    :param int time: timestamp
    :param array flags: flag array 
    
    :returns:
       cluster - updated cluster information
       obs_type - updated observation type
    '''

    if time - cluster.end >= 48:
        '''isolated cluster with 48hr gap either side'''

        # plotting done outside of this def.
        flags[cluster.locations] = 1
   
        # as have had a 48hr gap, start a new cluster 
        cluster.last_data = cluster.end
        cluster.start = time
        cluster.end = time
        cluster.locations = [time]
        cluster.length = 1
        obs_type = 2


    elif (time - cluster.start <= 24) and (cluster.length < 6):
        '''have data, but cluster is small and within thresholds --> increment'''
        obs_type = 2
        cluster.length += 1
        cluster.locations += [time]
        cluster.end = time

    else:
        '''actually it is now normal data, so reset'''

        obs_type = 0
        
        cluster.last_data = time
        cluster.start = cluster.data_mdi
        cluster.end   = cluster.data_mdi
        cluster.length = 0
        cluster.locations = 0
            
    return cluster, obs_type # occ_after_cluster


#*********************************************
def occ(station, variable_list, flag_col, datastart, logfile, diagnostics = False, plots = False):
    '''
    Check for odd clusters of data surrounded by missing 
        up to 6hr/24hr surrounded by at least 48 on each side

    :param MetVar station: the station object
    :param list variable_list: list of observational variables to process
    :param list flag_col: the columns to set on the QC flag array
    :param datetime datastart: dataset start time
    :param file logfile: logfile to store outputs
    :param bool diagnostics: do extra verbose output
    :param bool plots: do plots

    :returns:    
    '''

    # the four options of what to do with each observation
    #   the keys give values which are subroutines, and can be called
    #   all subroutines have to take the same set of inputs
    options = {0 : occ_normal, 1 : occ_start_cluster, 2 : occ_in_cluster, 3 : occ_after_cluster}

    for v,variable in enumerate(variable_list):
    
        st_var = getattr(station, variable)

        filtered_data = utils.apply_filter_flags(st_var)

        var_flags = station.qc_flags[:,flag_col[v]]
	
        prev_flag_number = 0

        # using IDL copy as method to ensure reproducibility (initially)
        
        oc_details = OddCluster(st_var.mdi, st_var.mdi, 0, st_var.mdi, st_var.mdi, -1)

        obs_type = 1

        for time in station.time.data:

            if filtered_data.mask[time] == False:
                # process observation point using subroutines, called from named tuple

                if plots and (obs_type == 3) and (time - oc_details.end >= 48):
                    # do plotting if matches flagging criteria
                    oc_plots(station, oc_details, time, datastart, filtered_data, variable, oc_details)

                oc_details, obs_type = options[obs_type](oc_details, obs_type, time, var_flags)

            else:
                # have missing data, 
        
                if obs_type  == 2:
                    obs_type = 3
                elif obs_type == 0:
                    obs_type = 1

        station.qc_flags[:,flag_col[v]] = var_flags

        flag_locs = np.where(station.qc_flags[:, flag_col[v]] != 0)
        utils.print_flagged_obs_number(logfile, "Odd Cluster", variable, len(flag_locs[0]) - prev_flag_number, noWrite = diagnostics)
        
        
        # copy flags into attribute
        st_var.flags[flag_locs] = 1

        # matches 032070 temperature 26/8/2014
    station = utils.append_history(station, "Isolated Odd Cluster Check")  

    return # occ


#************************************************************************
if __name__ == "__main__":

    print "Checking for Odd Clusters"
