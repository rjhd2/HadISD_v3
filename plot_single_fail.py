#!/usr/local/sci/bin/python
#*****************************
#
# plot data and flags for each year in turn
#
#
#************************************************************************
#                    SVN Info
#$Rev:: 92                                            $:  Revision of last commit
#$Author:: rdunn                                      $:  Author of last commit
#$Date:: 2016-02-08 12:00:09 +0000 (Mon, 08 Feb 2016) $:  Date of last commit
#************************************************************************

import numpy as np
import datetime as dt
import os
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import sys

# RJHD utils
import qc_utils as utils
import netcdf_procs as ncdfp
from set_paths_and_vars import *



input_id = sys.argv[1]
year = int(sys.argv[2])
variable = sys.argv[3]
test = sys.argv[4]


qc_test=np.array(['DUP','TFV','DFV','SFV','DNL','TGP','DGP','SGP','TRC','DRC',\
	  'WRC','PRC','TSS','DSS','WSS','PSS','HTS','HDS','HWS','HPS',\
	  'DTS','DDS','DWS','DPS','TCM','DCM','PCM','TSP','DSP','PSP',\
	  'SSS','DPD','DCF','CUOT','CUOL','CUOM','CUOH','CST','FLW','FMC',\
	  'NGC','TOT','DOT','SOT','TMB','DMB','SMB','WMB','BBB','CMB',\
	  'LMB','MMB','HMB','BMB','OCT','OCD','OCW','OCS','TVR','DVR',\
	  'SVR','WVR','WSL','WDL','WRS','STR_T','STR_D','STR_w','STR_S','ALL_T','ALL_Td','ALL_SLP','ACL'])

ms=20


try:
    station_info = np.genfromtxt(os.path.join(INPUT_FILE_LOCS, STATION_LIST), dtype=(str))
except IOError:
    print "station list not found"
    sys.exit()


for st,stat in enumerate(station_info):  

    if stat[0] == input_id:

        # set up station
        station = utils.Station(stat[0], float(stat[1]), float(stat[2]), float(stat[3]))
        break
else:
    sys.exit(0)

# read attributes and qc_flags
ncdfp.read(os.path.join(NETCDF_DATA_LOCS, station.id + "_external.nc"), station, process_vars, [], read_qc_flags = True)

match_to_compress = utils.create_fulltimes(station, process_vars, DATASTART, DATAEND, [])

# nyears x 12 months
month_start_locs = np.array(utils.month_starts(DATASTART, DATAEND)).reshape(-1,12)

# which years
years = DATASTART.year + np.arange(month_start_locs.shape[0])

# find which year and test to plot
year_loc, = np.where(years == year)    
test_loc, = np.where(qc_test == test)[0]

# and get the plot range
if year != DATAEND.year - 1:
    plot_range = (month_start_locs[year_loc,0], month_start_locs[year_loc+1,0])
else:
    plot_range = (month_start_locs[year_loc,0], -1) # misses last hour

# convert to useful numbers
plot_times = utils.times_hours_to_datetime(station.time.data[plot_range[0]:plot_range[1]], DATASTART)

# get all the QC flags
plot_qc_flags = station.qc_flags[plot_range[0]:plot_range[1],:]

plot_var = getattr(station, variable)
plot_data = plot_var.data[plot_range[0]:plot_range[1]]

plot_test_loc, = np.where(plot_qc_flags[:, test_loc] != 0)
plot_all_test_loc, = np.where(np.sum(plot_qc_flags[:,:], axis = 1) != 0)

plt.clf()


plt.scatter(plot_times, plot_data, c = 'k', marker = 'o', s = ms, edgecolor = 'k')
plt.scatter(plot_times[plot_all_test_loc], plot_data[plot_all_test_loc], c = 'c', marker = 'o', s = ms, edgecolor = '0.5')
plt.scatter(plot_times[plot_test_loc], plot_data[plot_test_loc], c = 'r', marker = 'o', s = ms, edgecolor = '0.5')


ymin, ymax = plt.gca().get_ylim()
if variable == "windspeeds":
    plt.ylim([-0.1*ymax, ymax])
elif variable in ["total_cloud_cover","low_cloud_cover","mid_cloud_cover","high_cloud_cover"]:
    plt.ylim([0, ymax])
    

plt.ylabel(plot_var.name.capitalize())
plt.title(input_id)

plt.show()
