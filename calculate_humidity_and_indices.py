#!/usr/local/sci/bin/python
#*****************************
#
# Calculates humidity quantities and heat-stress indices for new output file
#
#
#************************************************************************
#                    SVN Info
#$Rev:: 219                                           $:  Revision of last commit
#$Author:: rdunn                                      $:  Author of last commit
#$Date:: 2019-05-20 16:56:47 +0100 (Mon, 20 May 2019) $:  Date of last commit
#************************************************************************
'''
calculate_humidity_and_indices.py invoked by typing::

  python2.7 calculate_humidity_and_indices.py --restart_id 000000-99999 --end_id 999999-99999

Input arguments:

--restart_id        First station to process

--end_id            Last station to process
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
import humidity_vars as humidity
import heat_stress as heat_stress
from set_paths_and_vars import *


#************************************************************************
def make_hum_heat_vars(restart_id = "", end_id = "", diagnostics = False, plots = False):
    """
    Make the humidity and heat-stress variable netCDF files

    Make two sets of output files containing the humidity and heat-stress
    parameters calculated on an hourly basis from the QC'd HadISD data

    :param str restart_id: first station to process
    :param str end_id: last station to process
    :param bool diagnostics: verbose output to screen
    :param bool plots: make plots (placeholder)
    """

    # get station information
    try:
        station_info = np.genfromtxt(os.path.join(INPUT_FILE_LOCS, STATION_LIST), dtype=(str))
    except IOError:
        print "station list not found"
        sys.exit()

    # sort truncated run
    startindex = 0
    if restart_id != "":
        startindex, = np.where(station_info[:,0] == restart_id)


    if end_id != "":
        endindex, = np.where(station_info[:,0] == end_id)
        if endindex != len(station_info) -1:
            station_info = station_info[startindex[0]: endindex[0]+1]
        else:
            station_info = station_info[startindex[0]:]
    else:
        station_info = station_info[startindex[0]:]
        

    for st,stat in enumerate(station_info):       
        print dt.datetime.strftime(dt.datetime.now(), "%A, %d %B %Y, %H:%M:%S")
        print "{:35s} {:d}/{:d}".format("Station Number : ", st + 1, len(station_info))
        print "{:35s} {}".format("Station Identifier :", stat[0])

        if plots or diagnostics:
            logfile = ""
        else:
            logfile = file(LOG_OUTFILE_LOCS+stat[0]+'.log','a') 
            logfile.write(dt.datetime.strftime(dt.datetime.now(), "%A, %d %B %Y, %H:%M:%S\n"))
            logfile.write("Calculating Humidity and Heat Stress variables\n")
            logfile.write("{:35s} {}\n".format("Station Identifier :", stat[0]))
        process_start_time = time.time()


        station = utils.Station(stat[0], float(stat[1]), float(stat[2]), float(stat[3]))
        if os.path.exists(os.path.join(NETCDF_DATA_LOCS, "hadisd.{}_19310101-{}_{}.nc.gz".format(LONG_VERSION, END_TIME, station.id))):
            # if gzip file, unzip here
            subprocess.call(["gunzip",os.path.join(NETCDF_DATA_LOCS, "hadisd.{}_19310101-{}_{}.nc.gz".format(LONG_VERSION, END_TIME, station.id))])
            time.sleep(5) # make sure it is unzipped before proceeding

        # read in the data
        ncdfp.read(os.path.join(NETCDF_DATA_LOCS, "hadisd.{}_19310101-{}_{}.nc".format(LONG_VERSION, END_TIME, station.id)), station, process_vars, diagnostics = diagnostics, read_qc_flags = False, read_flagged_obs = False)

        match_to_compress = utils.create_fulltimes(station, process_vars, DATASTART, DATAEND, do_qc_flags = False, do_flagged_obs = False)

        # run through calculations, each one should add a new variable to object.


        """
        1) Use T and P to get e  [to get es, use Td]
        2) Use e, P, Td and T to get Tw
        3) If Tw < 0C, recalculate e w.r.t ice, and re-obtain Tw - keep both!
        4) Use e and P to calculate q
        5) Use e and es to get rh (use appropriate es too) - or q and qs

        what P to use if no measurement - using monthly mean probably isn't appropriate in this instance??
        """

        station = humidity.run_calcs(station, logfile)
        
        # run through heat stress calculations

        station = heat_stress.run_calcs(station, logfile)

        if diagnostics or plots: raw_input("stop")


        # adjust this to work with the desired output file - will need a separate write function - output humidity in one set, heat indices in another?
        humidity_vars = ["temperatures","dewpoints","slp","vapour_pressure","saturation_vapour_pressure","wetbulb_temperature","specific_humidity","relative_humidity"]
        ncdfp.write(os.path.join(NETCDF_DATA_LOCS, "hadisd.{}_19310101-{}_{}_humidity.nc".format(LONG_VERSION, END_TIME, station.id)), station, humidity_vars, os.path.join(INPUT_FILE_LOCS,'attributes.dat'), compressed = match_to_compress, processing_date = '', qc_code_version = '', write_QC_flags = False, write_flagged_obs = False, least_significant_digit = 5)

        heat_stress_vars = ["temperatures","dewpoints","windspeeds","THI","WBGT","humidex","apparent_t","heat_index"]
        ncdfp.write(os.path.join(NETCDF_DATA_LOCS, "hadisd.{}_19310101-{}_{}_heat_stress.nc".format(LONG_VERSION, END_TIME, station.id)), station, heat_stress_vars, os.path.join(INPUT_FILE_LOCS,'attributes.dat'), compressed = match_to_compress, processing_date = '', qc_code_version = '', write_QC_flags = False, write_flagged_obs = False, least_significant_digit = 5)

        # gzip the raw file
        # subprocess.call(["gzip","-f",os.path.join(NETCDF_DATA_LOCS, station.id + "_humidity.nc")])
        # subprocess.call(["gzip","-f",os.path.join(NETCDF_DATA_LOCS, station.id + "_heat_stress.nc")])
        # subprocess.call(["gzip",os.path.join(NETCDF_DATA_LOCS, station.id + "_mask.nc")])

        logfile.write(dt.datetime.strftime(dt.datetime.now(), "%A, %d %B %Y, %H:%M:%S\n"))
        logfile.write("processing took {:4.0f}s\n\n".format(time.time() - process_start_time))
        logfile.close()

        print "Humidity and Heat Stress Indices calculated"

    return # make_hum_heat_vars

#************************************************************************
if __name__=="__main__":

    import argparse

    # set up keyword arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--restart_id', dest='restart_id', action='store', default = "",
                        help='Restart ID for truncated run, default = ""')
    parser.add_argument('--end_id', dest='end_id', action='store', default = "",
                        help='End ID for truncated run, default = ""')

    args = parser.parse_args()


    make_hum_heat_vars(restart_id = args.restart_id, end_id = args.end_id)

#************************************************************************
