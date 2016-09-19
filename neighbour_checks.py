#!/usr/local/sci/bin/python
#*****************************
#
# controller for external QC checks.
#
#
#************************************************************************
#                    SVN Info
#$Rev:: 108                                           $:  Revision of last commit
#$Author:: rdunn                                      $:  Author of last commit
#$Date:: 2016-09-19 13:56:07 +0100 (Mon, 19 Sep 2016) $:  Date of last commit
#************************************************************************


import numpy as np
import scipy as sp
import os
import sys
import datetime as dt
import subprocess
import time

# RJHD utilities
import netcdf_procs as ncdfp
import qc_utils as utils
import neighbour_utils as n_utils
from set_paths_and_vars import *

import qc_tests



# look-up dictionaries for flag columns for different parts of this test
FLAG_OUTLIER_DICT = {"temperatures": 41, "dewpoints": 42, "slp": 43}
FLAG_COL_DICT = {"temperatures":np.array([0,1,4,5,8,12,16,20,24,27,41,54,58]),
                  "dewpoints":np.array([0,2,4,6,8,9,13,17,21,25,28,30,31,32,42,55,59]),
                  "slp":np.array([0,3,4,7,11,15,19,23,26,29,43,57,60]), 
                  "windspeeds":np.array([0,4,10,14,18,22,56,62,63,64]), 
                  "winddirs":np.array([0,4,10,14,18,22,56,62,63,64])}

N_NEIGHBOURS = 10



#*********************************************
def get_distances_angles(station_info):
    """
    Wrapper to calculate the distances and angles between all neighbours

    May seem superfluous, but means don't need to know details of n_utils if 
         just using this script

    :param list station_info: list of lists with strings of [ID, lat, lon, elev]
    :returns:
        distances - array of distances between station pairs
        angles - array of angles between station pairs
    """

    # calculate distances before station loop!
    distances, angles = n_utils.get_neighbour_parameters(station_info)
    
    return distances, angles # get_distances_angles


#*********************************************
def neighbour_checks(station_info, restart_id = "", end_id = "", distances=np.array([]), angles=np.array([]), second = False, masking = False, doZip=False, plots = False, diagnostics = False):
    """
    Run through neighbour checks on list of stations passed
    
    :param list station_info: list of lists - [[ID, lat, lon, elev]] - strings
    :param array distances: array of distances between station pairs
    :param array angles: array of angles between station pairs
    :param bool second: do the second run
    :param bool masking: apply the flags to the data to mask the observations.

    """
    first = not second

    qc_code_version = subprocess.check_output(['svnversion']).strip()

    # if distances and angles not calculated, then do so
    if (len(distances) == 0) or (len(angles) == 0):
        print "calculating distances and bearings matrix"
        distances, angles = get_distances_angles(station_info)

    # extract before truncate the array
    neighbour_elevations = np.array(station_info[:,3], dtype=float) 
    neighbour_ids        = np.array(station_info[:,0])
    neighbour_info       = np.array(station_info[:,:])

    # sort truncated run
    startindex = 0
    if restart_id != "":
        startindex, = np.where(station_info[:,0] == restart_id)


    if end_id != "":
        endindex, = np.where(station_info[:,0] == end_id)
        if endindex != len(station_info) -1:
            station_info = station_info[startindex: endindex+1]
            distances = distances[startindex:endindex+1,:]
            angles = angles[startindex:endindex+1,:]
        else:
            station_info = station_info[startindex:]
            distances = distances[startindex:,:]
            angles = angles[startindex:,:]
    else:
        station_info = station_info[startindex:]
        distances = distances[startindex:,:]
        angles = angles[startindex:,:]
        

    # process each neighbour
    for st, stat in enumerate(station_info):       

        print dt.datetime.strftime(dt.datetime.now(), "%A, %d %B %Y, %H:%M:%S")
        print "Neighbour Check"
        print "{:35s} {}".format("Station Identifier :", stat[0])

        if not plots and not diagnostics:
            logfile = file(LOG_OUTFILE_LOCS+stat[0]+'.log','a') # append to file if second iteration.
            logfile.write(dt.datetime.strftime(dt.datetime.now(), "%A, %d %B %Y, %H:%M:%S\n"))
            logfile.write("Neighbour Check\n")
            logfile.write("{:35s} {}\n".format("Station Identifier :", stat[0]))
        else:
            logfile = ""

        process_start_time = time.time()

        station = utils.Station(stat[0], float(stat[1]), float(stat[2]), float(stat[3]))

        # if running through the first time
        if first:

            if os.path.exists(os.path.join(NETCDF_DATA_LOCS, station.id + "_internal.nc.gz")):
                # if gzip file, unzip here
                subprocess.call(["gunzip",os.path.join(NETCDF_DATA_LOCS, station.id + "_internal.nc.gz")])
                time.sleep(5) # make sure it is unzipped before proceeding

            # read in the data
            ncdfp.read(os.path.join(NETCDF_DATA_LOCS, station.id + "_internal.nc"), station, process_vars, carry_thru_vars, diagnostics = diagnostics)

            if plots or diagnostics:
                print "{:35s}  {}\n".format("Total station record size :",len(station.time.data))
            else:
                logfile.write("{:35s}  {}\n".format("Total station record size :",len(station.time.data)))

            match_to_compress = utils.create_fulltimes(station, process_vars, DATASTART, DATAEND, carry_thru_vars)

        # or if second pass through?
        elif second:
            if os.path.exists(os.path.join(NETCDF_DATA_LOCS, station.id + "internal2.nc.gz")):
                # if gzip file, unzip here
                subprocess.call(["gunzip",os.path.join(NETCDF_DATA_LOCS, station.id + "_internal2.nc.gz")])
                time.sleep(5) # make sure it is unzipped before proceeding

            ncdfp.read(os.path.join(NETCDF_DATA_LOCS, station.id + "_internal2.nc"), station, process_vars, carry_thru_vars, diagnostics = diagnostics)
            if plots or diagnostics:
                print "{:35s}  {}\n".format("Total station record size :",len(station.time.data))
            else:
                logfile.write("{:35s}  {}\n".format("Total station record size :",len(station.time.data)))

            match_to_compress = utils.create_fulltimes(station, process_vars, DATASTART, DATAEND, carry_thru_vars)


        # select neighbours
        neighbour_distances  = distances[st,:]
        neighbour_bearings   = angles[st,:]

        # have to add in start index so that can use location in distance file.
        # neighbours = n_utils.get_neighbours(st+startindex, np.float(stat[3]), neighbour_distances, neighbour_bearings, neighbour_elevations)

        # return all neighbours up to a limit from the distance and elevation offsets (500km and 300m respectively)
        neighbours, neighbour_quadrants = n_utils.get_all_neighbours(st+startindex, np.float(stat[3]), neighbour_distances, neighbour_bearings, neighbour_elevations)

        if plots or diagnostics:
            print "{:14s} {:10s} {:10s}".format("Neighbour","Distance","Elevation")
            for n in neighbours:
                print "{:14s} {:10.1f} {:10.1f}".format(neighbour_ids[n],neighbour_distances[n],neighbour_elevations[n])

        else:
            logfile.write("{:14s} {:10s} {:10s}\n".format("Neighbour","Distance","Elevation"))
            for n in neighbours:
                logfile.write("{:14s} {:10.1f} {:10.1f}\n".format(neighbour_ids[n],neighbour_distances[n],neighbour_elevations[n]))

        # if sufficient neighbours
        if len(neighbours) >= 3:

            for variable, col in FLAG_OUTLIER_DICT.items():
                # NOTE - this requires multiple reads of the same file
                #      but does make it easier to understand and code

                st_var = getattr(station, variable)

                if plots or diagnostics:
                    print "Length of {} record: {}".format(variable, len(st_var.data.compressed()))
                else:
                    logfile.write("Length of {} record: {}\n".format(variable, len(st_var.data.compressed())))

                
                if len(st_var.data.compressed()) > 0:

                    final_neighbours = n_utils.select_neighbours(station, variable, neighbour_info[neighbours], neighbours, neighbour_distances[neighbours], neighbour_quadrants, NETCDF_DATA_LOCS, DATASTART, DATAEND, logfile, second = second, diagnostics = diagnostics, plots = plots)


                    # now read in final set of neighbours and process

                    neigh_flags = np.zeros(len(station.time.data)) # count up how many neighbours think this obs is bad
                    neigh_count = np.zeros(len(station.time.data)) # number of neighbours at each time stamp
                    dpd_flags = np.zeros(len(station.time.data)) # number of neighbours at each time stamp
                    reporting_accuracies = np.zeros(len(neighbours)) # reporting accuracy of each neighbour

                    all_data = np.ma.zeros([len(final_neighbours), len(station.time.data)]) # store all the neighbour values

                    for nn, nn_loc in enumerate(final_neighbours):

                        neigh_details = neighbour_info[nn_loc]
                        neigh = utils.Station(neigh_details[0], float(neigh_details[1]), float(neigh_details[2]), float(neigh_details[3]))

                        if first:
                            ncdfp.read(os.path.join(NETCDF_DATA_LOCS, neigh.id + "_internal.nc"), neigh, [variable], diagnostics = diagnostics, read_input_station_id = False)
                        elif second:
                            ncdfp.read(os.path.join(NETCDF_DATA_LOCS, neigh.id + "_internal2.nc"), neigh, [variable], diagnostics = diagnostics, read_input_station_id = False)

                        dummy = utils.create_fulltimes(neigh, [variable], DATASTART, DATAEND, [], do_input_station_id = False)

                        all_data[nn, :] = utils.apply_filter_flags(getattr(neigh, variable))

                        if diagnostics:
                            print neigh_details

                        n_utils.detect(station, neigh, variable, neigh_flags, neigh_count, DATASTART, DATAEND, distance = neighbour_distances[nn_loc], diagnostics = diagnostics, plots = plots)

                        reporting_accuracies[nn] = utils.reporting_accuracy(getattr(neigh,variable).data)

                        dpd_flags += neigh.qc_flags[:,31]
                    # gone through all neighbours


                    # if at least 2/3 of neighbours have flagged this point (and at least 3 neighbours)
                    some_flags, = np.where(neigh_flags > 0)            
                    outlier_locs, = np.where(np.logical_and((neigh_count[some_flags] >= 3),(neigh_flags[some_flags].astype("float")/neigh_count[some_flags] > 2./3.)))

                    # flag where < 3 neighbours
                    locs = np.where(neigh_count[some_flags] < 3)
                    station.qc_flags[some_flags[locs], col] = -1

                    if len(outlier_locs) >= 1:
                        station.qc_flags[some_flags[outlier_locs], col] = 1

                        # print number flagged and copy into attribute
                        if plots or diagnostics:
                            utils.print_flagged_obs_number(logfile, "Neighbour", variable, len(outlier_locs), noWrite = True)
                        else:
                            utils.print_flagged_obs_number(logfile, "Neighbour", variable, len(outlier_locs))
                        st_var = getattr(station, variable)
                        st_var.flags[some_flags[outlier_locs]] = 1

                    else:
                        if plots or diagnostics:
                            utils.print_flagged_obs_number(logfile, "Neighbour", variable, len(outlier_locs), noWrite = True)
                        else:
                            utils.print_flagged_obs_number(logfile, "Neighbour", variable, len(outlier_locs))


                    if plots:
                        n_utils.plot_outlier(station, variable, some_flags[outlier_locs], all_data, DATASTART)

                    # unflagging using neighbours
                    n_utils.do_unflagging(station, variable, all_data, reporting_accuracies, neigh_count, dpd_flags, FLAG_COL_DICT, DATASTART, logfile, plots = plots, diagnostics = diagnostics)

                else:
                    if plots or diagnostics:
                        print "No observations to assess for {}".format(variable)
                    else:
                        logfile.write("No observations to assess for {}\n".format(variable))
                    


            # variable loop
        else:
            if plots or diagnostics:
                print "Fewer than 3 neighbours"
            else:
                logfile.write("Fewer than 3 neighbours\n")

        print dt.datetime.strftime(dt.datetime.now(), "%A, %d %B %Y, %H:%M:%S")
        print "processing took {:4.0f}s\n\n".format(time.time() - process_start_time)

        # end of neighbour check
	utils.append_history(station, "Neighbour Outlier Check")

        # clean up months 

        qc_tests.clean_up.clu(station, ["temperatures","dewpoints","windspeeds","winddirs","slp"], [44,45,46,47,48], FLAG_COL_DICT, DATASTART, DATAEND, logfile, plots = plots)


        if diagnostics or plots: raw_input("stop")

        # masking (at least call from here - optional call from internal?)

        # write to file
        if first:
            ncdfp.write(os.path.join(NETCDF_DATA_LOCS, station.id + "_external.nc"), station, process_vars, os.path.join(INPUT_FILE_LOCS,'attributes.dat'), opt_var_list = carry_thru_vars, compressed = match_to_compress, processing_date = '', qc_code_version = qc_code_version)
            # gzip the raw file
        elif second:
            ncdfp.write(os.path.join(NETCDF_DATA_LOCS, station.id + "_external2.nc"), station, process_vars, os.path.join(INPUT_FILE_LOCS,'attributes.dat'), opt_var_list = carry_thru_vars, compressed = match_to_compress, processing_date = '', qc_code_version = qc_code_version)
            # gzip the raw file
 

        # masking - apply the flags and copy masked data to flagged_obs attribute
        if masking:

            station = utils.mask(station, process_vars, logfile)

        # write to file
            if first:
                ncdfp.write(os.path.join(NETCDF_DATA_LOCS, station.id + "_mask.nc"), station, process_vars, os.path.join(INPUT_FILE_LOCS,'attributes.dat'), opt_var_list = carry_thru_vars, compressed = match_to_compress, processing_date = '', qc_code_version = qc_code_version)
            elif second:
                ncdfp.write(os.path.join(NETCDF_DATA_LOCS, station.id + "_mask2.nc"), station, process_vars, os.path.join(INPUT_FILE_LOCS,'attributes.dat'), opt_var_list = carry_thru_vars, compressed = match_to_compress, processing_date = '', qc_code_version = qc_code_version)

        if plots or diagnostics:
            print "Masking completed\n"
            print dt.datetime.strftime(dt.datetime.now(), "%A, %d %B %Y, %H:%M:%S\n")
            print "processing took {:4.0f}s\n\n".format(time.time() - process_start_time)
        else:
            logfile.write("Masking completed\n")
            logfile.write(dt.datetime.strftime(dt.datetime.now(), "%A, %d %B %Y, %H:%M:%S\n"))
            logfile.write("processing took {:4.0f}s\n\n".format(time.time() - process_start_time))
            logfile.close()
            
    # gzip up all the raw files
    if doZip:
        for st, stat in enumerate(station_info):       
            if first:
                subprocess.call(["gzip",os.path.join(NETCDF_DATA_LOCS, stat[0]+"_internal.nc")])
                if masking:
                    subprocess.call(["gzip",os.path.join(NETCDF_DATA_LOCS, stat[0]+"_external.nc")])
                    subprocess.call(["gzip",os.path.join(NETCDF_DATA_LOCS, stat[0]+"_mask.nc")])

            elif second:
                subprocess.call(["gzip",os.path.join(NETCDF_DATA_LOCS, stat[0]+"_internal2.nc")])
                if masking:
                    subprocess.call(["gzip",os.path.join(NETCDF_DATA_LOCS, stat[0]+"_external2.nc")])
                    subprocess.call(["gzip",os.path.join(NETCDF_DATA_LOCS, stat[0]+"_mask2.nc")])

    print "Neighbour Checks completed\n"

    return # neighbour_checks 

#************************************************************************
if __name__=="__main__":

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--second', dest='second', action='store_true', default = False,
                        help='Second run through')
    parser.add_argument('--masking', dest='masking', action='store_true', default = False,
                        help='Apply flags')
    parser.add_argument('--do_zip', dest='doZip', action='store_true', default = False,
                        help='Zip files at end')
    parser.add_argument('--restart_id', dest='restart_id', action='store', default = "",
                        help='Restart ID for truncated run, default = ""')
    parser.add_argument('--end_id', dest='end_id', action='store', default = "",
                        help='End ID for truncated run, default = ""')
    parser.add_argument('--diagnostics', dest='diagnostics', action='store_true', default = False,
                        help='Run diagnostics (will not write out file)')
    parser.add_argument('--plots', dest='plots', action='store_true', default = False,
                        help='Run plots (will not write out file)')
    args = parser.parse_args()

    """To run as stand alone, process the file and obtain station list"""

    station_list = "candidate_stations.txt"

    try:
        station_info = np.genfromtxt(os.path.join(INPUT_FILE_LOCS, station_list), dtype=(str))
    except IOError:
        print "station list not found"
        sys.exit()

    uk = False
    if uk:
        uk_locs = []
        for s,station in enumerate(station_info[:,0]):
            if station[:2] == "03":
                uk_locs += [s]
                
        station_info = station_info[uk_locs]

    # station_info = [["030220-99999","57.4670","-7.3670","6.0000"]]

    
    neighbour_checks(station_info, restart_id = args.restart_id, end_id = args.end_id, second = args.second, masking = args.masking, doZip = args.doZip, diagnostics = args.diagnostics, plots = args.plots)

#    import cProfile

#    cProfile.run("neighbour_checks(station_info, restart_id = args.restart_id, end_id = args.end_id, second = args.second, masking = args.masking, doZip = args.doZip, diagnostics = args.diagnostics, plots = args.plots)")

