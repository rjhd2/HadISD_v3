#!/usr/local/sci/bin/python
#*****************************
#
# plot data and flags for each year in turn
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
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import sys

import qc_utils as utils
import netcdf_procs as ncdf

class QCtest(object):
    '''
    Class to hold QC test parameters
    '''

    def __init__(self, name, color, ymin, ymax):
        self.name = name
        self.color = color
        self.ymin = ymin
        self.ymax = ymax     


#*******************************************************
FILE_LOCS = "/project/hadobs2/hadisd/v200_2014/code_v200_2014/input_files/"
DATA_LOCS = "/project/hadobs2/hadisd/v200_2014/netcdf_files_v200_2014/"
IMG_LOCS = "/project/hadobs2/hadisd/v200_2014/img_files_v200_2014/"


DATASTART = dt.datetime(1931,1,1,0,0)
DATAEND = dt.datetime(2015,1,1,0,0)

station_list = "candidate_stations.txt"
process_vars = ["temperatures","dewpoints","slp","windspeeds","total_cloud_cover"]

input_id = "035715-99999"

# set up plot objects

Frequent=QCtest('Frequent','k',0,0.1)
Diurnal=QCtest('Diurnal','purple',0.1,0.2)
Gap=QCtest('Gap','b',0.1,0.2)
Record=QCtest('Record','g',0.2,0.3)
SString=QCtest('Straight Strings','g',0.7,0.8)
HString=QCtest('Hour Strings','g',0.8,0.9)
DString=QCtest('Day Strings','g',0.9,1.0)
Climatological=QCtest('Climatological','lime',0.3,0.4)
Spike=QCtest('Spike','y',0.4,0.5)
SuperSat=QCtest('Supersaturation','orange',0.5,0.6)
DPD=QCtest('Dewpoint Depression','r',0.6,0.7)
DewpCutOff=QCtest('Dewpoint Cut Off','purple',0.7,0.8)
Neighbour=QCtest('Neighbour','b',0.9,1.0)
MonthClean=QCtest('Month Removed','r',0.0,0.1)
OddCluster=QCtest('Odd Cluster','orange',0.1,0.2)
Variance=QCtest('Variance','y',0.2,0.3)
WindLogic=QCtest('Wind','r',0.3,0.4)
WindRose=QCtest('Wind','b',0.4,0.5)
# Clouds
UnobsT=QCtest('Unobservable - total','b',0.,0.25)
UnobsL=QCtest('Unobservable - low','c',0.25,0.5)
UnobsM=QCtest('Unobservable - medium','c',0.5,0.75)
UnobsH=QCtest('Unobservable - high','c',0.75,1.0)
MonthT=QCtest('Month Removed - total','g',0.,0.1)
MonthL=QCtest('Month Removed - low','r',0.1,0.2)
MonthM=QCtest('Month Removed - medium','r',0.2,0.3)
MonthH=QCtest('Month Removed - high','r',0.3,0.4)
MonthC=QCtest('Month Removed - CBase','orange',0.4,0.5)

qc_test=['DUP','TFV','DFV','SFV','DNL','TGP','DGP','SGP','TRC','DRC',\
	  'WRC','PRC','TSS','DSS','WSS','PSS','HTS','HDS','HWS','HPS',\
	  'DTS','DDS','DWS','DPS','TCM','DCM','PCM','TSP','DSP','PSP',\
	  'SSS','DPD','DCF','CUOT','CUOL','CUOM','CUOH','CST','FLW','FMC',\
	  'NGC','TOT','DOT','SOT','TMB','DMB','SMB','WMB','BBB','CMB',\
	  'LMB','MMB','HMB','BMB','OCT','OCD','OCW','OCS','TVR','DVR',\
	  'SVR','WVR','WSL','WDL','WRS','STR_T','STR_D','STR_w','STR_S','ALL_T','ALL_Td','ALL_SLP','ACL']
# dictionary of tests
qc_dict={'DUP':QCtest("Duplicate", '0.5',0,1),
         'TFV':Frequent,
         'DFV':Frequent,
         'SFV':Frequent,
         'DNL':Diurnal,
         'TGP':Gap,
         'DGP':Gap,
         'SGP':Gap,
         'TRC':Record,
         'DRC':Record,
         'WRC':Record,
         'PRC':Record,
         'TSS':SString,
         'DSS':SString,
         'WSS':SString,
         'PSS':SString,
         'HTS':HString,
         'HDS':HString,
         'HWS':HString,
         'HPS':HString,
         'DTS':DString,
         'DDS':DString,
         'DWS':DString,
         'DPS':DString,
         'TCM':Climatological,
         'DCM':Climatological,
         'PCM':Climatological,
         'TSP':Spike,
         'DSP':Spike,
         'PSP':Spike,
         'SSS':SuperSat,
         'DPD':DPD,
         'DCF':DewpCutOff,
         'CUOT':UnobsT,
         'CUOL':UnobsL,
         'CUOM':UnobsM,
         'CUOH':UnobsH,
         'CST':QCtest("Small Total", '0.5',0,1), 
         'FLW':QCtest("Full Low", '0.5',0,1), 
         'FMC':QCtest("Full Mid", '0.5',0,1), 
         'NGC':QCtest("Negative Value", '0.5',0,1),
         'TOT':Neighbour,
         'DOT':Neighbour,
         'SOT':Neighbour,
         'TMB':MonthClean,
         'DMB':MonthClean,
         'SMB':MonthClean,
         'WMB':MonthClean,
         'BBB':MonthClean,
         'CMB':MonthT,
         'LMB':MonthL,
         'MMB':MonthM,
         'HMB':MonthH,
         'BMB':MonthC,
         'OCT':OddCluster,
         'OCD':OddCluster,
         'OCW':OddCluster,
         'OCS':OddCluster,
         'TVR':Variance,
         'DVR':Variance,
         'SVR':Variance,
         'WVR':Variance,
         'WSL':WindLogic,
         'WDL':WindLogic,
         'WRS':WindRose
         }


# tests to use for each variable
tests = {"temperatures" : [0,1,4,5,8,12,16,20,24,27,41,44,54,58],
         "dewpoints":[0,2,4,6,8,9,13,17,21,25,28,30,31,32,42,45,55,59],
         "slp":[0,3,4,7,11,15,19,23,26,29,43,46,57,60],
         "windspeeds":[10,14,18,22,47,56,61,62,63,64],
         'total_cloud_cover':range(33,41)}




try:
    station_info = np.genfromtxt(os.path.join(FILE_LOCS, station_list), dtype=(str))
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
ncdf.read(os.path.join(DATA_LOCS, station.id + "_external.nc"), station, process_vars, [])

match_to_compress = utils.create_fulltimes(station, process_vars, DATASTART, DATAEND, [])

# nyears x 12 months
month_start_locs = np.array(utils.month_starts(DATASTART, DATAEND)).reshape(-1,12)

# which years
years = DATASTART.year + np.arange(month_start_locs.shape[0])

for year in range(DATASTART.year, DATAEND.year):

    year_loc, = np.where(years == year)

    if year != DATAEND.year - 1:
        plot_range = (month_start_locs[year_loc,0], month_start_locs[year_loc+1,0])
    else:
        plot_range = (month_start_locs[year_loc,0], -1) # misses last hour

    plot_times = utils.times_hours_to_datetime(station.time.data[plot_range[0]:plot_range[1]], DATASTART)

    plot_qc_flags = station.qc_flags[plot_range[0]:plot_range[1],:]

    plt.figure(figsize=(12, 12), dpi=100)
    plt.clf()

    MakePlot = False
    for v,var in enumerate(process_vars):

        ax = plt.subplot(5, 1, v+1)

        plot_var = getattr(station, var)
        plot_data = plot_var.data[plot_range[0]:plot_range[1]]

        if len(plot_data.compressed()) > 0:
            MakePlot = True

            plt.plot(plot_times, plot_data, 'k.', zorder=3, alpha = 0.5)  # 'k,' gives pixel marker

            ymin, ymax = plt.gca().get_ylim()
            if var == "windspeeds":
                plt.ylim([-0.1*ymax, ymax])
            elif var in ["total_cloud_cover","low_cloud_cover","mid_cloud_cover","high_cloud_cover"]:
                plt.ylim([0, ymax])
                

            plt.ylabel(plot_var.name.capitalize())

            if v != len(process_vars)-1:
                plt.setp(ax.get_xticklabels(), visible=False)

            # run through each test
            for test_number in tests[plot_var.name]:

                these_flags = plot_qc_flags[:,test_number]

                locs, = np.where(these_flags != 0)

                test = qc_dict[qc_test[test_number]] # horrible coding

                # if something to plot, then do
                if len(locs) > 0:

                    for loc in locs:

                        plt.axvline(plot_times[loc], ymin = test.ymin, ymax = test.ymax, color = test.color, alpha=0.5)

    if MakePlot:
        # only show if there is something to show!
        plt.gcf().subplots_adjust(hspace=0.0, wspace=0.05, top=0.95, bottom=0.05, left=0.075, right=0.925)
        plt.suptitle("{} - {}".format(station.id, DATASTART.year+year_loc[0]))
        watermarkstring="/".join(os.getcwd().split('/')[4:])+'/'+os.path.basename( __file__ )+"   "+dt.datetime.strftime(dt.datetime.now(), "%d-%b-%Y %H:%M")
        plt.figtext(0.01,0.01,watermarkstring,size=5)

        plt.savefig(IMG_LOCS+"{}_{}_fails.png".format(station.id, DATASTART.year+year_loc[0]))
        plt.close()
    else:
        plt.close()
 
