#!/usr/local/sci/bin/python
#*****************************
#
# print failure rate for each test for each station
#
#
#************************************************************************
#                    SVN Info
#$Rev:: 77                                            $:  Revision of last commit
#$Author:: rdunn                                      $:  Author of last commit
#$Date:: 2015-07-06 16:23:09 +0100 (Mon, 06 Jul 2015) $:  Date of last commit
#************************************************************************

import numpy as np
import datetime as dt
import os

import netcdf_procs as ncdf
import qc_utils as utils

#*******************************************************

DATA_LOCS = "/project/hadobs2/hadisd/v200_2014/netcdf_files_v200_2014/"
FILE_LOCS = "/project/hadobs2/hadisd/v200_2014/code_v200_2014/input_files/"
IMG_LOCS = "/project/hadobs2/hadisd/v200_2014/code_v200_2014/images/"

qc_test=['DUP','TFV','DFV','SFV','DNL','TGP','DGP','SGP','TRC','DRC',\
	  'WRC','PRC','TSS','DSS','WSS','PSS','HTS','HDS','HWS','HPS',\
	  'DTS','DDS','DWS','DPS','TCM','DCM','PCM','TSP','DSP','PSP',\
	  'SSS','DPD','DCF','CUOT','CUOL','CUOM','CUOH','CST','FLW','FMC',\
	  'NGC','TOT','DOT','SOT','TMB','DMB','SMB','WMB','BBB','CMB',\
	  'LMB','MMB','HMB','BMB','OCT','OCD','OCW','OCS','TVR','DVR',\
	  'SVR','WVR','WSL','WDL','WRS']

station_list = "candidate_stations.txt"
process_vars = ["temperatures","dewpoints","slp","windspeeds","total_cloud_cover"]
diagnostics = False
start_time_string = dt.datetime.strftime(dt.datetime.now(), "%Y%m%d")

try:
    station_info = np.genfromtxt(os.path.join(FILE_LOCS, station_list), dtype=(str))
except IOError:
    print "station list not found"
    sys.exit()


Lons = []
Lats = []

outfile = file(FILE_LOCS+"qc_summary_{}.dat".format(start_time_string),'w')

for st,stat in enumerate(station_info):  


    # set up station
    station = utils.Station(stat[0], float(stat[1]), float(stat[2]), float(stat[3]))

#    if station.id[:2] != "03":
#        continue
    print st, station.id
    

    # read attributes and qc_flags
    ncdf.read(os.path.join(DATA_LOCS, station.id + "_mask.nc"), station, process_vars, [], diagnostics = diagnostics)

    # sum qc_flags:
    # remove multi-level flagging
    qc_flags = station.qc_flags[:]

    qc_flags[qc_flags[:] > 1] = 1

    # remove multi-level flagging - neighbour flags
    no_neighbours = qc_flags[qc_flags[:] == -1].size
    qc_flags[qc_flags[:] < 0] = 0

    total_flags = qc_flags[qc_flags[:] != 0].size

    sum_flags = np.sum(qc_flags[:], axis = 0) # 61 column array

    outfile.write("{:12s} {:7d} {:7.0f}".format(stat[0], len(station.time.data), np.sum(sum_flags))+''.join(['%8i' % n for n in sum_flags])+"\n")

outfile.write("{:12s} {:7s} {:7s}".format("station", "times", "flags")+''.join(['%8s' % n for n in qc_test])+"\n")
outfile.close()
