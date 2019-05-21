#!/usr/local/sci/bin/python
#*****************************
#
# controller for external QC checks.
#
#
#************************************************************************
#                    SVN Info
#$Rev:: 219                                           $:  Revision of last commit
#$Author:: rdunn                                      $:  Author of last commit
#$Date:: 2019-05-20 16:56:47 +0100 (Mon, 20 May 2019) $:  Date of last commit
#************************************************************************
'''
neighbour_checks.py invoked by typing::

  python2.7 neighbour_checks.py --restart_id 000000-99999 --end_id 999999-99999 --masking --do_zip  [--plots] [--diagnostics]

Input arguments:

--restart_id        First station to process

--end_id            Last station to process

--masking           Apply the results of the QC tests and move these observations to new fields

--plots             [False] Create plots into the image directory

--diagnostics       [False] Verbose output
'''





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
FLAG_COL_DICT = { "temperatures":np.array([0,1,4,5,8,12,16,20,24,27,41,44,54,58]),
                  "dewpoints":np.array([0,2,4,6,8,9,13,17,21,25,28,30,31,32,42,45,55,59]),
                  "slp":np.array([0,3,4,7,11,15,19,23,26,29,43,46,57,60]), # 26 should be empty
                  "windspeeds":np.array([0,4,10,14,18,22,47,56,61,62,63,64,65]), 
                  "winddirs":np.array([0,4,10,14,18,22,47,48,56,61,62,63,64,65,66,67,68]),
                  "clouds": np.array([33,34,35,36,37,38,39,40]), # for completeness, but unused
                  "total_cloud_cover":[33,37,40],
                  "low_cloud_cover":[34,38,40],
                  "mid_cloud_cover":[35,38,39,40],
                  "high_cloud_cover":[36,38,39,40],
                  "stnlp" : [69],
                  "precip1_depth" : [70],
                  "precip2_depth" : [70],
                  "precip3_depth" : [70],
                  "precip6_depth" : [70],
                  "precip9_depth" : [70],
                  "precip12_depth" : [70],
                  "precip15_depth" : [70],
                  "precip18_depth" : [70],
                  "precip24_depth" : [70]}

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
def neighbour_checks(restart_id = "", end_id = "", distances=np.array([]), angles=np.array([]), masking = False, doZip=False, plots = False, diagnostics = False):
    """
    Run through neighbour checks on list of stations passed
    
    :param str restart_id: which station to start on
    :param str end_id: which station to end on
    :param array distances: array of distances between station pairs
    :param array angles: array of angles between station pairs
    :param bool masking: apply the flags to the data to mask the observations.
    :param bool doZip: use netCDF4 compression
    :param bool plots: produce plots on the fly
    :param bool diagnostics: verbose output

    """

#    qc_code_version = subprocess.check_output(['svnversion']).strip()
    qc_code_version = subprocess.check_output(['svn', 'info', 'file:///home/h05/rdunn/svn/hadisd_py_qc/branches/monthly/'])
    for line in qc_code_version.split("\n"):
        if line.split(":")[0] == "Revision":
            qc_code_version = line.split(":")[1]
            break

    # get station information
    try:
        station_info = np.genfromtxt(os.path.join(INPUT_FILE_LOCS, STATION_LIST), dtype=(str))
    except IOError:
        print "station list not found"
        sys.exit()

    # if distances and angles not calculated, then do so
    if (len(distances) == 0) or (len(angles) == 0):
        print "calculating distances and bearings matrix"
        distances, angles = get_distances_angles(station_info)

    # extract before truncate the array
    neighbour_elevations = np.array(station_info[:,3], dtype=float) 
    neighbour_ids        = np.array(station_info[:,0])
    neighbour_info       = np.array(station_info[:,:])

    # sort truncated run
    startindex = [0]
    if restart_id != "":
        startindex, = np.where(station_info[:,0] == restart_id)


    if end_id != "":
        endindex, = np.where(station_info[:,0] == end_id)
        if endindex != len(station_info) -1:
            station_info = station_info[startindex[0]: endindex[0]+1]
            distances = distances[startindex[0]:endindex[0]+1,:]
            angles = angles[startindex[0]:endindex[0]+1,:]
        else:
            station_info = station_info[startindex[0]:]
            distances = distances[startindex[0]:,:]
            angles = angles[startindex[0]:,:]
    else:
        station_info = station_info[startindex[0]:]
        distances = distances[startindex[0]:,:]
        angles = angles[startindex[0]:,:]
        

    # process each neighbour
    for st, stat in enumerate(station_info):       

        print dt.datetime.strftime(dt.datetime.now(), "%A, %d %B %Y, %H:%M:%S")
        print "Neighbour Check"
        print "{:35s} {}".format("Station Identifier :", stat[0])

        if not plots and not diagnostics:
            logfile = file(LOG_OUTFILE_LOCS+stat[0]+'.log','a') 
            logfile.write(dt.datetime.strftime(dt.datetime.now(), "%A, %d %B %Y, %H:%M:%S\n"))
            logfile.write("Neighbour Check\n")
            logfile.write("{:35s} {}\n".format("Station Identifier :", stat[0]))
        else:
            logfile = ""

        process_start_time = time.time()

        station = utils.Station(stat[0], float(stat[1]), float(stat[2]), float(stat[3]))

        if os.path.exists(os.path.join(NETCDF_DATA_LOCS, station.id + "_internal.nc.gz")):
            # if gzip file, unzip here
            subprocess.call(["gunzip",os.path.join(NETCDF_DATA_LOCS, station.id + "_internal.nc.gz")])
            time.sleep(5) # make sure it is unzipped before proceeding

        # read in the data
        ncdfp.read(os.path.join(NETCDF_DATA_LOCS, "hadisd.{}_19310101-{}_{}_internal.nc".format(LONG_VERSION, END_TIME, station.id)), station, process_vars, carry_thru_vars, diagnostics = diagnostics)

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

                    final_neighbours = n_utils.select_neighbours(station, variable, neighbour_info[neighbours], neighbours, neighbour_distances[neighbours], neighbour_quadrants, NETCDF_DATA_LOCS, DATASTART, DATAEND, logfile, diagnostics = diagnostics, plots = plots)


                    # now read in final set of neighbours and process

                    neigh_flags = np.zeros(len(station.time.data)) # count up how many neighbours think this obs is bad
                    neigh_count = np.zeros(len(station.time.data)) # number of neighbours at each time stamp
                    dpd_flags = np.zeros(len(station.time.data)) # number of neighbours at each time stamp
                    reporting_accuracies = np.zeros(len(neighbours)) # reporting accuracy of each neighbour

                    all_data = np.ma.zeros([len(final_neighbours), len(station.time.data)]) # store all the neighbour values

                    for nn, nn_loc in enumerate(final_neighbours):

                        neigh_details = neighbour_info[nn_loc]
                        neigh = utils.Station(neigh_details[0], float(neigh_details[1]), float(neigh_details[2]), float(neigh_details[3]))

                        ncdfp.read(os.path.join(NETCDF_DATA_LOCS, "hadisd.{}_19310101-{}_{}_internal.nc".format(LONG_VERSION, END_TIME, station.id)), neigh, [variable], diagnostics = diagnostics, read_input_station_id = False)

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
        
        # clean up months - no difference for odd months - should just run through.

        qc_tests.clean_up.clu(station, ["temperatures","dewpoints","slp","windspeeds","winddirs"], [44,45,46,47,48], FLAG_COL_DICT, DATASTART, DATAEND, logfile, plots = plots, diagnostics = diagnostics)
        utils.apply_flags_from_A_to_B(station, "slp", "stnlp", diagnostics = diagnostics)

        if diagnostics or plots: raw_input("stop")

        # masking (at least call from here - optional call from internal?)

        # write to file
        ncdfp.write(os.path.join(NETCDF_DATA_LOCS, "hadisd.{}_19310101-{}_{}_external.nc".format(LONG_VERSION, END_TIME, station.id)), station, process_vars, os.path.join(INPUT_FILE_LOCS,'attributes.dat'), opt_var_list = carry_thru_vars, compressed = match_to_compress, processing_date = dt.datetime.strftime(dt.datetime.now(), "%d-%b-%Y"), qc_code_version = qc_code_version)
 

        # masking - apply the flags and copy masked data to flagged_obs attribute
        if masking:

            station = utils.mask(station, process_vars, logfile, FLAG_COL_DICT)

            # write to file
            ncdfp.write(os.path.join(NETCDF_DATA_LOCS, "hadisd.{}_19310101-{}_{}.nc".format(LONG_VERSION, END_TIME, station.id)), station, process_vars, os.path.join(INPUT_FILE_LOCS,'attributes.dat'), opt_var_list = carry_thru_vars, compressed = match_to_compress, processing_date = dt.datetime.strftime(dt.datetime.now(), "%d-%b-%Y"), qc_code_version = qc_code_version)

        if plots or diagnostics:
            print "Masking completed\n"
            print dt.datetime.strftime(dt.datetime.now(), "%A, %d %B %Y, %H:%M:%S\n")
            print "processing took {:4.0f}s\n\n".format(time.time() - process_start_time)
        else:
            logfile.write("Masking completed\n")
            logfile.write(dt.datetime.strftime(dt.datetime.now(), "%A, %d %B %Y, %H:%M:%S\n"))
            logfile.write("processing took {:4.0f}s\n\n".format(time.time() - process_start_time))
            logfile.close()
            
    # looped through all stations

    # gzip up all the raw files
    if doZip:
        for st, stat in enumerate(station_info):       
            subprocess.call(["gzip",os.path.join(NETCDF_DATA_LOCS, "hadisd.{}_19310101-{}_{}_internal.nc".format(LONG_VERSION, END_TIME, station.id))])
            if masking:
                subprocess.call(["gzip",os.path.join(NETCDF_DATA_LOCS, "hadisd.{}_19310101-{}_{}_external.nc".format(LONG_VERSION, END_TIME, station.id))])
                subprocess.call(["gzip",os.path.join(NETCDF_DATA_LOCS, "hadisd.{}_19310101-{}_{}.nc".format(LONG_VERSION, END_TIME, station.id))])


    print "Neighbour Checks completed\n"

    return # neighbour_checks 

#************************************************************************
if __name__=="__main__":

    import argparse
    parser = argparse.ArgumentParser()
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
    
    neighbour_checks(restart_id = args.restart_id, end_id = args.end_id, masking = args.masking, doZip = args.doZip, diagnostics = args.diagnostics, plots = args.plots)

#    import cProfile

#    cProfile.run("neighbour_checks(station_info, restart_id = args.restart_id, end_id = args.end_id, second = args.second, masking = args.masking, doZip = args.doZip, diagnostics = args.diagnostics, plots = args.plots)")

