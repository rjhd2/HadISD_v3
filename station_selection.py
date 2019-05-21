#!/usr/local/sci/bin/python2.7
#------------------------------------------------------------
#                    SVN Info
#$Rev:: 219                                           $:  Revision of last commit
#$Author:: rdunn                                      $:  Author of last commit
#$Date:: 2019-05-20 16:56:47 +0100 (Mon, 20 May 2019) $:  Date of last commit
#------------------------------------------------------------
# START
#------------------------------------------------------------
"""
station_selection.py invoked by typing::

  python2.7 station_selection.py --datastart 1931 --data_end 2018 --plots --update_files [--uk]

Downloads the ISD inventories and constructs a set of stations to use.

Input arguments:

--data_start        [1931]           First year of data

--data_end          [current year]   Last year of data

--years_present     [15]             Number of years between required between listed start and end date to select station

--update_files      [False]          Download updated files from ISD server

--plots             [False]          Create plots into the image directory

--uk                [False]          Run only for UK stations (03*)
"""

import struct
import os
import datetime as dt
import numpy as np
import sys
import argparse
import matplotlib
# use the Agg environment to generate an image rather than outputting to screen
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

# RJHD Utilities
import sort_canada as s_canada
import qc_utils as utils
import selection_utils as sel_utils
from set_paths_and_vars import *
import compare_station_lists as compare_lists


ISD_LISTING = "isd-history.txt"
ISD_INVENTORY = "isd-inventory.txt"

# what is available
START_YEAR = 1901
END_YEAR = dt.date.today().year

DAYS_IN_AVERAGE_MONTH=[31,28,31,30,31,30,31,31,30,31,30,31]


#*********************************************
class Station(object):
    """
    Class for station

    At bare minimum contains .id .lat .lon .elev attributes, set in that order
    """
    
    def __init__(self, stn_id, lat, lon, elev):
        self.id = stn_id
        self.lat = lat
        self.lon = lon
        self.elev = elev
        
    def __str__(self):
        return "station {}, lat {}, lon {}, elevation {}".format(self.id, self.lat, self.lon, self.elev)
    
    __repr__ = __str__
 
#*********************************************
def read_stations(limit, start_year,  uk = False, germany = True, canada = True):
    """
    Process the ISD history file
    
    Select stations which have defined lat, lon and elevation;
    at least N years of data (using start/end dates).  Also do additional
    processing for Canadian (71*) and German (09* and 10*) stations.

    :param int limit: number of years of data required
    :param int start_year: year of data start
    :param bool uk: only run for stations starting 03* 
    :param bool germany: do extra selection for 09* and 10* stations
    :param bool canada: do extra selection and processing for 71* stations

    :returns: list of selected station objects
    """

    fieldwidths = (7,6,30,5,3,5,8,9,8,9,9)
    fmtstring = ''.join('%ds' % f for f in fieldwidths)
    parse = struct.Struct(fmtstring).unpack_from

    all_stations = []

    try:
        with open(INPUT_FILE_LOCS+ISD_LISTING, 'r') as infile:

            lc = 0
            for line in infile:
                lc += 1
                if lc <= 22:

                    continue

                fields = parse(line)

                # unpack
                st_id = fields[0].strip()+'-'+fields[1].strip()
                name = fields[2].rstrip()
                country = fields[3].strip()
                state = fields[4].strip()
                call = fields[5].strip()
                lat = fields[6].strip()
                lon = fields[7].strip()
                elev = fields[8].strip()
                start = dt.datetime.strptime(fields[9].strip(),"%Y%m%d").date()
                end = dt.datetime.strptime(fields[10].strip(),"%Y%m%d").date()

                if uk:
                    if st_id[:2] != "03":
                        continue


                # need to have lat, lon and elev and the elev defined
                if lat != "" and lon != "" and elev != "" and float(elev) != -999.9:

                    # test if elevation is known (above Dead-Sea shore in Jordan at -423m)
                    if float(elev) > -430.0:

                        # create station object
                        station = Station(st_id, float(lat), float(lon), float(elev))
                        station.name = name
                        station.state = state
                        station.call = call
                        station.start = start
                        station.end = end
                        station.mergers = []
                        all_stations += [station]


    except OSError:
        print "{:s} does not exist.  Check download".format(isd_listing)
        raise OSError
    
    print "{} stations in full ISD".format(lc)

    # go through German stations and pick out those which match in last 4 places of ID
    if germany:
        
        all_stations = process_germany(all_stations)

    candidate_stations = []

    # go through Canadian stations to check for undocumented moves/mergers
    if canada:

        s_canada.process_EC_file()

        all_stations = s_canada.process_canadian_stations(all_stations, start_year)

    for station in all_stations:

        # need to have sufficient number of years of data
        if station.end.year - station.start.year >= limit:
            
            # store in list
            candidate_stations += [station]


    print "{} candidate stations (lat/lon/elev and timespan limits)".format(len(candidate_stations))
  
    return all_stations, candidate_stations # read_stations


#*********************************************
def extract_inventory(station, inventory, data_start, data_end, do_mergers = True):
    """
    Extract the information from the ISD inventory file

    :param station station: station object to extract data for
    :param array inventory: the full ISD inventory to extract information from
    :param int data_start: first year of data
    :param int data_end: last year of data
    :param bool do_mergers: include the information from merged stations when doing selection cuts
    :returns: array of Nyears x 12 months containing number of observation records in each month
    """

    monthly_obs = np.zeros([END_YEAR - START_YEAR + 1, 12])

    # are WMO IDs sufficient?
    locs, = np.where(np.logical_and(inventory[:,0] == int(station.id[:6]), inventory[:,1] == int(station.id[7:])))
    this_station = inventory[locs,:]

    for year in this_station:
        monthly_obs[year[2] - START_YEAR, :] = year[3:]

    if do_mergers:
        if station.mergers != []:

            for mstation in station.mergers:

                merger_station = inventory[inventory[:,0] == int(mstation[:6]),:]

                for year in merger_station:
                    # only replace if have more observations
                    # can't tell how the merge will work with interleaving obs, so just take the max
                    locs = np.where(year[3:] > monthly_obs[year[2] - START_YEAR, :])
                    monthly_obs[year[2] - START_YEAR, locs] = year[3:][locs]

    # if restricted data period
    if data_start != START_YEAR:
        monthly_obs = monthly_obs[(data_start - START_YEAR):, :]

    if data_end != END_YEAR:
        monthly_obs = monthly_obs[:(data_end - END_YEAR), :]

    return monthly_obs # extract_inventory


#*********************************************
def process_inventory(candidate_stations, data_start, data_end, min_years_present, do_mergers = True):

    """
    Process the ISD inventory file

    Use input list to select those which have at least 6 hourly reporting over min_years_present years.

    :param list candidate_stations: list of station objects to match up
    :param int data_start: first year of data
    :param int data_end: last year of data
    :param int min_years_present: how many years need to be there for station to be selected
    :param bool do_mergers: include the information from merged stations when doing selection cuts
    :returns: list of selected station objects
    """
    

    print "filtering ISD inventory"

    # read in the inventory
    try:
        inventory = np.genfromtxt(INPUT_FILE_LOCS+ISD_INVENTORY, skip_header = 8, dtype = int)
    except OSError:
        pass

    final_stations = []
    # spin through each station
    for s, station in enumerate(candidate_stations):
#        print "{}/{}".format(s, len(candidate_stations))

        # extract the observations in each month for this station
        monthly_obs = extract_inventory(station, inventory, data_start, data_end, do_mergers = True)

        # mask where no observations (should keep array structure)
        obs = np.ma.masked_where(monthly_obs <= 0, monthly_obs)

        # only store if have > 6 hourly on average (4 * 30) 
        if np.ma.median(obs) > 120:

            # now mask where insufficient observations (should keep array structure)
            obs = np.ma.masked_where(monthly_obs <= 120, monthly_obs)

            # and have N years of 12 months data at 6 hourly intervals
            if len(obs.compressed()) > (12*min_years_present):

                # test for strange data structures - so do a monthly check for median of 6 hourly data.
                monthly_medians = np.ma.median(obs, axis = 0).astype("int")

                if all(mm >= 4 * DAYS_IN_AVERAGE_MONTH[i] for i, mm in enumerate(monthly_medians)):

                    station.obs = monthly_obs
                    final_stations += [station]

    print "{} candidate stations (6 hourly reporting) and data in at least {} months".format(len(final_stations),12*min_years_present)

    inventory = [] #  clear memory    

    return final_stations # process_inventory

#*********************************************
def stations_per_year(candidate_stations, data_start):
    """
    Return list of which years each station is present in

    :param list candidate_station: list of station objects
    :param str data_start: year of truncated record start

    :returns: stations_active_in_years - list of lists - which years a station is active in.   

    """

    # set up blank lists
    stations_active_in_years=[[] for stn in candidate_stations]

    # spin through stations
    for s, station in enumerate(candidate_stations):

        # process each year at a time - using .obs attribute (monthly_obs from extract_inventory)
        for y, year in enumerate(station.obs):

            if np.sum(year) != 0:

                stations_active_in_years[s] += [y + data_start]

    return stations_active_in_years # stations_per_year

#*********************************************
def process_germany(all_stations):

    """go through German stations and pick out those which match in last 4 places of ID"""

    g_count = 0
    station_ids = np.array([stn.id for stn in all_stations])

    german_ids = np.array([s for s, stn in enumerate(station_ids) if stn[:2] == "10"])

    for gid in german_ids:

        old_id_loc, = np.where(station_ids == "09"+station_ids[gid][2:])

        if len(old_id_loc) != 0:
            # there is a match - check distance, elev and date.
            old_stn = all_stations[old_id_loc[0]]
            probs = sel_utils.do_match(old_stn, all_stations[gid], sel_utils.LATITUDE_THRESHOLD, sel_utils.ELEVATION_THRESHOLD, sel_utils.DISTANCE_THRESHOLD)

            if np.product(probs) > sel_utils.PROB_THRESHOLD:
                g_count += 1

                print "Merging {} and {}".format(all_stations[gid].id, old_stn.id)

                all_stations[gid].mergers = [old_stn.id]
                if old_stn.start < all_stations[gid].start:
                    all_stations[gid].start = old_stn.start
                if old_stn.end > all_stations[gid].end:
                    all_stations[gid].end = old_stn.end

                # ensure this station won't be selected
                all_stations[old_id_loc[0]].start = dt.datetime.today()
                all_stations[old_id_loc[0]].end = dt.datetime.today()

    print "{} German stations merged".format(g_count)

    return all_stations # process_germany


#*********************************************
def doftp_transfer(infile, diagnostics = False):
    """
    Do the FTP transfer of the files

    :param infile outfile: name of file to retrieve
    :param bool diagnostics: extra verbosity
    :returns:
    """

    HOST=r'ftp.ncdc.noaa.gov'
    ISD_LOC=r'/pub/data/noaa/'

    if diagnostics: print "getting {:s}".format(infile)

    os.system('doftp -host '+HOST+' -cwd /pub/data/noaa/ -get '+infile+'='+INPUT_FILE_LOCS+'/'+infile)

    if diagnostics: print "got {:s}".format(infile)

    return # doftp

#*********************************************
def plot_stations(station_list, outfile, merges = False, title = ""):
    """
    Plot the stations on a global map

    :param list station_list: list of station objects
    :param str outfile: name of output file
    :param bool merges: plot merged stations in different colour
    :param str title: plot title
    :returns:
    """


    plt.figure(figsize=(8,5))
    plt.clf()
    ax = plt.axes([0,0,1,1],projection = ccrs.Robinson())
    ax.coastlines('50m')
    try:
        ax.gridlines(draw_labels = True)
    except TypeError:
        ax.gridlines()

    lats, lons = [], []
    mlats, mlons =  [], []

    for stn in station_list:

        lats += [stn.lat]
        lons += [stn.lon]

        if stn.mergers != []:

            mlats += [stn.lat]
            mlons += [stn.lon]

    ax.scatter(lons, lats, transform = ccrs.Geodetic(), s = 15, c = 'b', edgecolor = 'b', label = 'HadISD stations')

    if merges:
        ax.scatter(mlons, mlats, transform = ccrs.Geodetic(), c = 'r', edgecolor = 'r', s = 15, label = "mergers")

    # for dummy points at max E,W,S,N to get full globe (set_extent didn't work well)
    ax.scatter([-179.9,179.9,0,0], [0,0,-90,90], transform = ccrs.Geodetic(), s = 1, c = 'w', edgecolor = 'w')

    plt.legend(loc='lower center',ncol=2, bbox_to_anchor=(0.5,-0.1),frameon=False,prop={'size':13})
    if title!= "":
        plt.suptitle("{} - {} stations".format(title, len(lats)))
    else:
        plt.suptitle("{} stations".format(len(lats)))

    watermarkstring="/".join(os.getcwd().split('/')[4:])+'/'+os.path.basename( __file__ )+"   "+dt.datetime.strftime(dt.datetime.now(), "%d-%b-%Y %H:%M")
    plt.figtext(0.01,0.01,watermarkstring,size=5)

    plt.savefig(IMAGE_LOCS+outfile)
    plt.close()

    return # plot_stations


#****************************************************
def plot_station_number_over_time(candidate_stations, data_start, outfile, after_merges = False, unmerged_station_numbers = []):
    """
    Plot the number of stations in each year

    :param list station_list: list of station objects
    :param int data_start: year of truncated record start
    :param str outfile: name of output file
    :param bool after_merges: plot stations after merger as well as before
    :returns:
    """

    # flatten list
    stations_active_in_years = stations_per_year(candidate_stations, data_start)

    station_years = np.array([item for sublist in stations_active_in_years for item in list(set(sublist))] )

    plt.clf()

    # if this is post merger
    if after_merges:
        assert unmerged_station_numbers != [] # should have been set on call
        
        # plot the unmerged
        plt.plot(range(START_YEAR, END_YEAR+1),unmerged_station_numbers, 'co', ls='-', label = "raw")

        # calculate the merged
        station_numbers_merged = []

        for y in range(START_YEAR, END_YEAR+1):
            station_numbers_merged += [len(station_years[station_years == y])]

        plt.plot(range(START_YEAR, END_YEAR+1),station_numbers_merged, 'rs', ls='-', label = "merged")

        # do a legend in this case
        plt.legend(loc = "upper left", frameon=False)

    else:
        # if before the merging, then make the raw station number list.
        station_numbers = []

        for y in range(START_YEAR, END_YEAR+1):
            station_numbers += [len(station_years[station_years == y])]

        plt.plot(range(START_YEAR, END_YEAR+1),station_numbers, 'co', ls='-', label = "raw")

    # prettify the plot
    plt.ylim([0,8000])
    plt.ylabel("Stations with data")
    plt.xlim([START_YEAR,END_YEAR+1])
    plt.xticks(range(1900,2020,10))
    watermarkstring="/".join(os.getcwd().split('/')[4:])+'/'+os.path.basename( __file__ )+"   "+dt.datetime.strftime(dt.datetime.now(), "%d-%b-%Y %H:%M")
    plt.figtext(0.01,0.01,watermarkstring,size=5)


    plt.savefig(IMAGE_LOCS+outfile, transparent=True)

    if after_merges:
        return station_numbers_merged # plot_station_number_over_time
    else:
        return station_numbers # plot_station_number_over_time

#****************************************************
def find_matches(candidate_stations):
    """
    Find the matches within the station list.  

    Returns a cross-matched list

    :param list candidate_stations: list of station objects to test
    :returns:
        list of locations of matches with  match[i]=j & match[j]=i 
    """

    print "finding merger candidate matches"
    match = [[] for i in range(len(candidate_stations))]

    outfile = file(INPUT_FILE_LOCS+'internal_merging_log.txt','w')
    outfile.write("Station1 Station2, Probs, Prob-product\n\n")
    for i1,stn1 in enumerate(candidate_stations):

        for i2,stn2 in enumerate(candidate_stations):
            if i2 <= i1: continue

            probs = sel_utils.do_match(stn1, stn2, sel_utils.LATITUDE_THRESHOLD, sel_utils.ELEVATION_THRESHOLD, sel_utils.DISTANCE_THRESHOLD)

            if np.product(probs) > sel_utils.PROB_THRESHOLD:

                match[i1] += [i2] 
                match[i2] += [i1] 

                outfile.write("{:30s} {:13s} {:30s} {:13s} {} {:8.4f}\n".format(stn1.name, stn1.id, stn2.name, stn2.id, " ".join(["{:8.4f}".format(p) for p in probs]), np.product(probs)))
    outfile.close()

    return match # find_matches

#****************************************************
def process_matches(match, candidate_stations):
    """
    Sort through the list of matches to determine which is the primary station

    :param list match: list of cross matches
    :param list candidate_stations: list of station objects
    
    :returns:
       ordered_merges = list of ordered merge locations (first entry is primary station in each case)
    """


    ordered_merges = []

    for mc,merge_candidates in enumerate(match):

        if merge_candidates != []:

            # first item is the station that's being assessed
            these_stns = [mc]

            prev_len = len(these_stns)

            # add on the other stations which have been matched to the first one
            these_stns.extend(merge_candidates)
            match[mc] = [] # and wipe this record in matched to avoid double counting

            # go through each of these subsequent station records
            #   want to aggregate all possible matches together.
            while prev_len != len(these_stns):
                prev_len = len(these_stns)

                # check if the stn has been processed
                for stn in these_stns:
                    if match[stn] != []:
                        # those not dealt with
                        these_stns.extend(match[stn])
                        match[stn] = []

                these_stns = list(set(these_stns))

            # have list of stations to merge.  Extract order - just in terms of station length

            length = []
            for stn in these_stns:

                length += [candidate_stations[stn].end.year - candidate_stations[stn].start.year]

            sort_order = np.argsort(np.array(length))

            ordered_merges += [np.array(these_stns)[sort_order[::-1]]]

    return ordered_merges # process_matches

#****************************************************
def internal_merges(ordered_merges, candidate_stations):
    """
    Take the list of ordered stations to be merged, and store in station object.

    Remove these stations from the station list and return the new list.

    :param list ordered_merges: list of lists with ordered merge station IDs
    :param list candidate_stations: list of stations to process

    :returns:
       new_candidates - updated list of stations
    """


    # reduce the list and put the candidate IDs into the objects
    for merges in ordered_merges:
        # merges contains parent station followed by secondaries
        candidate_stations[merges[0]].mergers = [candidate_stations[loc].id for loc in merges[1:]]  

        # sort the .obs attribute
        for loc in merges[1:]:
            # np.maximum returns the max of each array, element wise
            candidate_stations[merges[0]].obs = np.maximum(candidate_stations[merges[0]].obs, candidate_stations[loc].obs)

        for secondary_stn in merges[1:]:
            candidate_stations[secondary_stn] = 0

    # and remove the stations which will be merged as secondaries
    new_candidates = []
    for stn in candidate_stations:
        if stn != 0:
            new_candidates += [stn]

    print "{} original candidates, {} merged, leaving {}".format(len(candidate_stations), len(candidate_stations) - len(new_candidates), len(new_candidates))

    return new_candidates # internal_merges

#****************************************************
def external_merges(candidate_stations, all_stations, data_start, data_end):
    """
    Run the merging candidate selection against the full ISD stations

    Unconfound where possible to ensure each station only included once!
    Update attributes to store the potential matches
    
    :param list candidate_stations: the stations selected so far
    :param list all_stations: the full ISD list
    :param int data_start: first year of data
    :param int data_end: last year of data

    :returns:
        candidate_stations - updated to include new merger candidates
    """

    outfile = file(INPUT_FILE_LOCS+'external_merging_log.txt','w')
    outfile.write("Station1, Station2, Probs, Prob-product\n\n")

    # now check matches against full 29k ISD stations
    match = [[] for i in range(len(candidate_stations))]
    reverse_match = [[] for i in range(len(all_stations))] # so that candidates are not merged into two stations!!

    for c_stn, candidate in enumerate(candidate_stations):
        for a_stn, target in enumerate(all_stations):
            if candidate.id != target.id:
                probs = sel_utils.do_match(candidate, target, sel_utils.LATITUDE_THRESHOLD, sel_utils.ELEVATION_THRESHOLD, sel_utils.DISTANCE_THRESHOLD)
                if np.product(probs) > sel_utils.PROB_THRESHOLD:
                   match[c_stn] += [a_stn]
                   reverse_match[a_stn] +=[[c_stn, np.product(probs)]]
                   outfile.write("{:30s} {:13s} {:30s} {:13s} {} {:8.4f}\n".format(candidate.name, candidate.id, target.name, target.id, " ".join(["{:8.4f}".format(p) for p in probs]), np.product(probs)))


    # check no stations being merged into two primaries
    for s, stn in enumerate(reverse_match):
        if len(stn) > 1:
            print stn
            # find maximum likelihood match
            max_prob = np.argmax(np.array(stn)[:,1])
            # and remove all the others
            for p, potential in enumerate(stn):
                if p != max_prob:
                    match[potential[0]] = [stnid for stnid in match[potential[0]] if stnid != s ]


    try:
        inventory = np.genfromtxt(INPUT_FILE_LOCS+ISD_INVENTORY, skip_header = 8, dtype = int)
    except OSError:
        pass

    # write into the attribute
    for l,locs in enumerate(match):
        if locs != []:
            candidate_stations[l].mergers += [all_stations[loc].id for loc in locs]
            candidate_stations[l].mergers = list(set(candidate_stations[l].mergers)) # as will have picked up internal set again!
            # sort the .obs attribute
            # np.maximum returns the max of each array, element wise
            monthly_obs = extract_inventory(all_stations[loc], inventory, data_start, data_end, do_mergers = False)

            candidate_stations[l].obs = np.maximum(candidate_stations[l].obs, monthly_obs)

    inventory = [] #  clear memory    

    return candidate_stations # external_merges


#*********************************************
def run_selection(data_start = 1931, 
                  data_end = dt.datetime.now().year, 
                  min_years_present = 15, 
                  updateFiles = False, 
                  plots = True, 
                  uk = False):

    """
    Main selection and merger script

    Processes ISD history and inventory files to extract list of possible stations.
    
    German (09* and 10*) and Canadian (71*) stations have additional processing

    Initially select those stations which have defined latitude, longitude and elevation.
    Then those which have 15 years between listed start and end date.  Finally use the 
    inventory file to select those which have 6 hourly observations on average.  

    Check which stations within this final list can be merged together to obtain
    final station number.  Then compare to full ISD and extract further files which can be merged.

    :param int data_start: first year of final dataset
    :param int data_end: last year of final dataset
    :param int min_years_present: number of years of data present
    :param bool updateFiles: update the ISD input files
    :param bool plots: make plots
    :param bool uk: run for UK stations only
    """

    # get the isd_history file to do a gross filtering
    if updateFiles:
        doftp_transfer(ISD_LISTING)


    # parse text file into candidate list, with lat, lon, elev and time span limits applied
    all_stations, candidate_stations = read_stations(min_years_present, start_year = data_start, uk = uk, canada = True, germany = True)      

    if plots: plot_stations(candidate_stations, 'isd_lat_lon_elev_{}.png'.format(HADISD_VERSION))
    # write candidate files
    outfile = file(INPUT_FILE_LOCS+"candidate_stations_lat_lon_elev.txt",'w')
    for stn in candidate_stations:
        outfile.write("{:12s} {:7.3f} {:8.3f} {:7.1f}\n".format(stn.id, stn.lat, stn.lon, stn.elev))
    outfile.close()


    # get the isd_inventory file to do a more detailed selection
    if updateFiles:
        doftp_transfer(ISD_INVENTORY)


    candidate_stations = process_inventory(candidate_stations, data_start, data_end, min_years_present)

    if plots: plot_stations(candidate_stations, 'isd_6hourly_{}.png'.format(HADISD_VERSION))

    # write candidate files
    outfile = file(INPUT_FILE_LOCS+"candidate_stations_6hourly.txt",'w')
    for stn in candidate_stations:
        outfile.write("{:12s} {:7.3f} {:8.3f} {:7.1f}\n".format(stn.id, stn.lat, stn.lon, stn.elev))
    outfile.close()

    # plot distribution of stations
    if plots: 
        unmerged_station_numbers = plot_station_number_over_time(candidate_stations, data_start, 'hadisd_station_number_raw_{}.png'.format(HADISD_VERSION)) 

        

    # do matching within this list
    match = find_matches(candidate_stations)

    # sort station lists - take longest of all matched pairs, and set merging order

    ordered_merges = process_matches(match, candidate_stations)


    # write out internal matches
    outfile = file(INPUT_FILE_LOCS+"internal_mergers.txt", 'w')
    for merges in ordered_merges:
        stn_ids = [candidate_stations[stn].id for stn in merges]  
        outfile.write(" ".join(stn_ids)+"\n")
    outfile.close()

    # store the internal merges in the station objects and reform the list
    candidate_stations = internal_merges(ordered_merges, candidate_stations)

    # now check matches against full 29k ISD stations
    candidate_stations = external_merges(candidate_stations, all_stations, data_start, data_end)

    # write final set of merges
    outfile = file(INPUT_FILE_LOCS + MERGER_LIST, 'w')                
    for stn in candidate_stations:
        if stn.mergers != []:
            outfile.write(stn.id + " "+" ".join(stn.mergers)+"\n")
    outfile.close()

    # write candidate files
    outfile = file(INPUT_FILE_LOCS + STATION_LIST,'w')
    for stn in candidate_stations:
        outfile.write("{:12s} {:7.3f} {:8.3f} {:7.1f}\n".format(stn.id, stn.lat, stn.lon, stn.elev))
    outfile.close()

    # write output files - full details
    outfile = file(INPUT_FILE_LOCS+"candidate_stations_details.txt",'w')
    for stn in candidate_stations:
        outfile.write("{:12s} {:30s} {:7.3f} {:8.3f} {:7.1f} {} {}\n".format(stn.id, stn.name, stn.lat, stn.lon, stn.elev, stn.start, stn.end))
    outfile.close()


    if plots: 
        # plot the station locations
        plot_stations(candidate_stations, 'isd_6hourly_mergers_{}.png'.format(HADISD_VERSION), merges = True)

        # plot of stations over time
        plot_station_number_over_time(candidate_stations, data_start, 'hadisd_station_number_merged_{}.png'.format(HADISD_VERSION), after_merges = True, unmerged_station_numbers = unmerged_station_numbers)
    
        #*************************
        # plot a gridded map

        # gridcell size
        delta_lon, delta_lat = 1.5, 1.5
        gridmax = 20.

        # create array of cell boundaries and make grid
        rawlons=np.arange(-180.,180.+delta_lon,delta_lon)
        rawlats=np.arange(-90.,90.+delta_lat,delta_lat)
        gridlon, gridlat=np.meshgrid(rawlons,rawlats)
        # set up empty array for gridded data
        griddata=np.zeros(list(gridlon.shape))

        UsedStation=np.zeros(len(candidate_stations))

        lats, lons = [], []
        for stn in candidate_stations:
            lats += [stn.lat]
            lons += [stn.lon]

        nstations = len(candidate_stations)

        StationNumbers=np.zeros(gridlon.shape)
        start=0
        for tlons,longitude in enumerate(gridlon[0,1:]):

            for tlats,latitude in enumerate(gridlat[1:,0]):

                for st in range(start,nstations):
                     if UsedStation[st]!=1.:
                        # if station not already in a grid box
                        if lats[st] < latitude and \
                                lons[st] < longitude:
                            """counts from bottom LH corner upwards
                            so starts at -180, -90.
                            Grid box values are to the top right of the
                            coordinates, with the final set ignored
                            This counts the stations to the bottom left
                            so do an offset to "place" result in "correct"
                            grid box.  No stations are "<" -180,-90
                            """
                            StationNumbers[tlats,tlons]+=1
                            UsedStation[st]=1.

        print "Total number of stations passing criteria ",np.sum(StationNumbers)
        print "Grid maximum set to "+str(gridmax)

        import cartopy.crs as ccrs
        plt.figure(figsize=(8,5.5))
        plt.clf()
        ax = plt.axes([0.05,0.05,.9,.9],projection = ccrs.Robinson())
        ax.coastlines('50m')
        try:
            ax.gridlines(draw_labels = True)
        except TypeError:
            ax.gridlines()

        # mask the grid for zeros.
        StationNumbers[StationNumbers > gridmax] = gridmax
        MaskedGriddata=np.ma.masked_where(StationNumbers==0,StationNumbers)

        # plot the grid - set max for colourbar at 50
        cs=plt.pcolormesh(gridlon, gridlat, MaskedGriddata, cmap=plt.cm.jet_r, alpha=0.8, vmax=gridmax, vmin=1, transform = ccrs.PlateCarree())

        # colour bar
        cb=plt.colorbar(cs,orientation='horizontal', pad=0.04, fraction=0.04, ticks=[1,5,10,15,20], shrink = 0.8)
        cb.set_label('Number of Stations')

        plt.suptitle('Station distribution for 6 hourly reporting frequency')
        watermarkstring="/".join(os.getcwd().split('/')[4:])+'/'+os.path.basename( __file__ )+"   "+dt.datetime.strftime(dt.datetime.now(), "%d-%b-%Y %H:%M")
        plt.figtext(0.01,0.01,watermarkstring,size=5)
        #plt.show()
        plt.savefig(IMAGE_LOCS+'hadisd_gridded_station_distribution_{}.png'.format(HADISD_VERSION), dpi=200)


        # plot station numbers in all years
        for year in range(data_start, data_end + 1):

            plot_list = []

            # spin through all stations
            for stn in candidate_stations:

                # if at least one observation in that year, select it.
                if np.sum(stn.obs[year - data_start]) > 0:
                    plot_list += [stn]

            plot_stations(plot_list, "hadisd_station_number_in_{}_{}.png".format(year, HADISD_VERSION), title = str(year)),
            print year, len(plot_list)


    # and compare the station lists from previous and this version
    compare_lists.compare()

    return # run_selection

#*******************
if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--start', dest='start', action='store', default = 1931,
                        help='Start year, default = 1931')
    parser.add_argument('--end', dest='end', action='store', default = dt.datetime.now().year,
                        help='End year, default = current')
    parser.add_argument('--years_present', dest='yearsPresent', action='store', default = 15,
                        help='Number of years of data required')
    parser.add_argument('--update_files', dest='updateFiles', action='store_true', default = False,
                        help='Update the input files from the ISD server')
    parser.add_argument('--plots', dest='plots', action='store_true', default = False,
                        help='Do all the output plots')
    parser.add_argument('--uk', dest='uk', action='store_true', default = False,
                        help='UK stations only')

    args = parser.parse_args()

    print "Available data from {} to {}".format(START_YEAR, END_YEAR)
    print "Extracting and testing from {} to {}".format(args.start, args.end)


    run_selection(int(args.start), int(args.end), int(args.yearsPresent), args.updateFiles, args.plots, args.uk)

#------------------------------------------------------------
# END
#------------------------------------------------------------
