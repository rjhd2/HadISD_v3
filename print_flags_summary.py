#!/usr/local/sci/bin/python
#*****************************
#
# prints out the totals of the big array of flags for each station
#
#************************************************************************
#                    SVN Info
#$Rev:: 52                                            $:  Revision of last commit
#$Author:: rdunn                                      $:  Author of last commit
#$Date:: 2015-01-28 15:01:15 +0000 (Wed, 28 Jan 2015) $:  Date of last commit
#************************************************************************

import numpy as np
import scipy as sp
import os
import sys
import datetime as dt

# RJHD utilities
import netcdf_procs as ncdf
import qc_utils as utils


FILE_LOCS = "/project/hadobs2/hadisd/v200_2014/code_v200_2014/input_files/"
DATA_LOCS = "/project/hadobs2/hadisd/v200_2014/netcdf_files_v200_2014/"
OUTFILE_LOCS = "/project/hadobs2/hadisd/v200_2014/suppl_files_v200_2014/"

#FILE_LOCS = "/home/h05/rdunn/LocalData/HadISD/QC/uk_testing_hadobs_copy/infiles/"
#DATA_LOCS = "/home/h05/rdunn/LocalData/HadISD/QC/uk_testing_hadobs_copy/netcdffiles_v000_2013x/"
#OUTFILE_LOCS = "/home/h05/rdunn/LocalData/HadISD/QC/uk_testing_hadobs_copy/suppl_files_v000_2013x/"

station_list = "candidate_stations.txt"



try:
    station_info = np.genfromtxt(os.path.join(FILE_LOCS, station_list), dtype=(str))
except IOError:
    print "station list not found"
    sys.exit()



with open(os.path.join(OUTFILE_LOCS,'qc_summary.txt'),'w') as outfile:

    outfile.write("     station   times   fails     DUP     TFV     DFV     SFV     DNL     TGP     DGP     SGP     TRC     DRC     WRC     PRC     TSS     DSS     WSS     PSS     HTS     HDS     HWS     HPS     DTS     DDS     DWS     DPS     TCM     DCM     PCM     TSP     DSP     PSP     SSS     DPD     DCF    CUOT    CUOL    CUOM    CUOH     CST     FLW     FMC     NGC     TOT     DOT     SOT     TMB     DMB     SMB     WMB     BBB     CMB     LMB     MMB     HMB     BMB     OCT     OCD     OCW     OCS     TVR     DVR     SVR\n")

    for st,stat in enumerate(station_info):  

        print "{:35s} {:d}/{:d}".format("Station Number : ", st + 1, len(station_info))

        station = utils.Station(stat[0], float(stat[1]), float(stat[2]), float(stat[3]))

        ncdf.read(os.path.join(DATA_LOCS, station.id + "_mask.nc"), station, ['temperatures'])

        qc_flags = station.qc_flags[:]

        # create array of "ones" and mask all locations where no flag was set
        totals = np.ma.masked_where(qc_flags <= 0, np.ma.ones(qc_flags.shape))

        # sum down columns to get total for each test
        totals = np.ma.sum(totals, axis = 0)
        totals.fill_value = 0 

        # make string of all these numbers
        flag_string = "".join(["{:8.0f}".format(flag) for flag in totals.data])

        # and print to file
        print_string = "{:12s}{:8d}{:8.0f}{:s}\n".format(stat[0], len(station.time.data), np.ma.sum(totals), flag_string)

        outfile.write(print_string)


    outfile.write("     station   times   fails     DUP     TFV     DFV     SFV     DNL     TGP     DGP     SGP     TRC     DRC     WRC     PRC     TSS     DSS     WSS     PSS     HTS     HDS     HWS     HPS     DTS     DDS     DWS     DPS     TCM     DCM     PCM     TSP     DSP     PSP     SSS     DPD     DCF    CUOT    CUOL    CUOM    CUOH     CST     FLW     FMC     NGC     TOT     DOT     SOT     TMB     DMB     SMB     WMB     BBB     CMB     LMB     MMB     HMB     BMB     OCT     OCD     OCW     OCS     TVR     DVR     SVR\n")


    
