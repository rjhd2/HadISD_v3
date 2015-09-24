#!/usr/local/sci/bin/python
#*****************************
#
# Calculates humidity quantities and heat-stress indices for new output file
#
#
#************************************************************************
#                    SVN Info
#$Rev:: 76                                            $:  Revision of last commit
#$Author:: rdunn                                      $:  Author of last commit
#$Date:: 2015-06-15 19:03:14 +0100 (Mon, 15 Jun 2015) $:  Date of last commit
#************************************************************************

import numpy as np
import scipy as sp
import os
import sys
import datetime as dt
import subprocess
import time

# RJHD utilities
import netcdf_procs as ncdf
import qc_utils as utils
import humidity_vars as humidity
import heat_stress as heat_stress


DATASTART = dt.datetime(1931,1,1,0,0)
DATAEND = dt.datetime(2015,1,1,0,0)

FILE_LOCS = "/project/hadobs2/hadisd/v200_2014/code_v200_2014/input_files/"
DATA_LOCS = "/project/hadobs2/hadisd/v200_2014/netcdf_files_v200_2014/"
OUTFILE_LOCS = "/project/hadobs2/hadisd/v200_2014/suppl_files_v200_2014/"


process_vars = ["temperatures","dewpoints","slp","windspeeds","input_station_id"]



def make_hum_heat_vars(station_info, restart_id = "", end_id = "", diagnostics = False, plots = False):
    '''
    


    '''


    # sort truncated run
    startindex = 0
    if restart_id != "":
        startindex, = np.where(station_info[:,0] == restart_id)


    if end_id != "":
        endindex, = np.where(station_info[:,0] == end_id)
        if endindex != len(station_info) -1:
            station_info = station_info[startindex: endindex+1]
        else:
            station_info = station_info[startindex:]
    else:
        station_info = station_info[startindex:]
        

    for st,stat in enumerate(station_info):       
        print dt.datetime.strftime(dt.datetime.now(), "%A, %d %B %Y, %H:%M:%S")
        print "{:35s} {:d}/{:d}".format("Station Number : ", st + 1, len(station_info))
        print "{:35s} {}".format("Station Identifier :", stat[0])



        station = utils.Station(stat[0], float(stat[1]), float(stat[2]), float(stat[3]))
        if os.path.exists(os.path.join(DATA_LOCS, station.id + "_mask.nc.gz")):
            # if gzip file, unzip here
            subprocess.call(["gunzip",os.path.join(DATA_LOCS, station.id + "_mask.nc.gz")])
            time.sleep(5) # make sure it is unzipped before proceeding

        # read in the data
        ncdf.read(os.path.join(DATA_LOCS, station.id + "_mask.nc"), station, process_vars, diagnostics = diagnostics, read_qc_flags = False, read_flagged_obs = False)

        match_to_compress = utils.create_fulltimes(station, process_vars, DATASTART, DATAEND, do_qc_flags = False, do_flagged_obs = False)

        # run through calculations, each one should add a new variable to object.


        '''
        1) Use T and P to get e  [to get es, use Td]
        2) Use e, P, Td and T to get Tw
        3) If Tw < 0C, recalculate e w.r.t ice, and re-obtain Tw - keep both!
        4) Use e and P to calculate q
        5) Use e and es to get rh (use appropriate es too) - or q and qs

        what P to use if no measurement - using monthly mean probably isn't appropriate in this instance??
        '''

        station = humidity.run_calcs(station)
        
        # run through heat stress calculations

        station = heat_stress.run_calcs(station)

        if diagnostics or plots: raw_input("stop")


        # adjust this to work with the desired output file - will need a separate write function - output humidity in one set, heat indices in another?
        humidity_vars = ["temperatures","dewpoints","slp","vapour_pressure","saturation_vapour_pressure","wetbulb_temperature","specific_humidity","relative_humidity"]
        ncdf.write(os.path.join(DATA_LOCS, station.id + "_humidity.nc"), station, humidity_vars, os.path.join(FILE_LOCS,'attributes.dat'), compressed = match_to_compress, processing_date = '', qc_code_version = '', write_QC_flags = False, write_flagged_obs = False, least_significant_digit = 5)

        heat_stress_vars = ["temperatures","dewpoints","windspeeds","THI","WBGT","humidex","apparent_t","heat_index"]
        ncdf.write(os.path.join(DATA_LOCS, station.id + "_heat_stress.nc"), station, heat_stress_vars, os.path.join(FILE_LOCS,'attributes.dat'), compressed = match_to_compress, processing_date = '', qc_code_version = '', write_QC_flags = False, write_flagged_obs = False, least_significant_digit = 5)
        # gzip the raw file
        subprocess.call(["gzip","-f",os.path.join(DATA_LOCS, station.id + "_humidity.nc")])
        subprocess.call(["gzip","-f",os.path.join(DATA_LOCS, station.id + "_heat_stress.nc")])
        # subprocess.call(["gzip",os.path.join(DATA_LOCS, station.id + "_mask.nc")])

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

    '''To run as stand alone, process the file and obtain station list'''
    station_list = "candidate_stations.txt"

    try:
        station_info = np.genfromtxt(os.path.join(FILE_LOCS, station_list), dtype=(str))
    except IOError:
        print "station list not found"
        sys.exit()


    make_hum_heat_vars(station_info, restart_id = args.restart_id, end_id = args.end_id)
