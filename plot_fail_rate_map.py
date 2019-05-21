#!/usr/local/sci/bin/python
#*****************************
#
# plot failure rate for each test and overall for variables
#
#
#************************************************************************
#                    SVN Info
#$Rev:: 219                                           $:  Revision of last commit
#$Author:: rdunn                                      $:  Author of last commit
#$Date:: 2019-05-20 16:56:47 +0100 (Mon, 20 May 2019) $:  Date of last commit
#************************************************************************
'''
plot_fail_rate_map.py invoked by typing::

  python2.7 plot_fail_rate_map.py

Runs automatically, spinning through each variable and QC test
'''

import numpy as np
import datetime as dt
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

# RJHD utils
import qc_utils as utils
import netcdf_procs as ncdfp
from set_paths_and_vars import *


#*******************************************************
qc_test=['DUP','TFV','DFV','SFV','DNL','TGP','DGP','SGP','TRC','DRC',\
	  'WRC','PRC','TSS','DSS','WSS','PSS','HTS','HDS','HWS','HPS',\
	  'DTS','DDS','DWS','DPS','TCM','DCM','PCM','TSP','DSP','PSP',\
	  'SSS','DPD','DCF','CUOT','CUOL','CUOM','CUOH','CST','FLW','FMC',\
	  'NGC','TOT','DOT','SOT','TMB','DMB','SMB','WMB','BBB','CMB',\
	  'LMB','MMB','HMB','BMB','OCT','OCD','OCW','OCS','TVR','DVR',\
	  'SVR','WVR','WSL','WDL','WRS','WSP','RSS','HRS','DRS','STNLP','PPTN',\
          'STR_T','STR_D','STR_W','STR_WD','STR_S','ALL_T','ALL_Td','ALL_SLP','ALL_W','ALL_WD','ACL']

# list of tests not the same as where they were applied to!
T_QC=[0,1,4,5,8,12,16,20,24,27,41,44,54,58]
D_QC=[0,2,4,6,8,9,13,17,21,25,28,30,31,32,42,45,55,59]
S_QC=[0,3,4,7,11,15,19,23,26,29,43,46,57,60,69] # 26 should be empty (added SLP/STNLP)
WS_QC=[0,10,14,18,22,47,56,61,62,63,64,65]
WD_QC=[0,10,14,18,22,47,48,56,61,62,63,64,65,66,67,68]
C_QC=range(33,41)

strT_QC=[12,16,20]
strD_QC=[13,17,21]
strWS_QC=[14,18,22]
strWD_QC=[66,67,68]
strS_QC=[15,19,23]

diagnostics = False
start_time_string = dt.datetime.strftime(dt.datetime.now(), "%Y%m%d")


#************************************************************************
def make_test_dictionary():

    with open(os.path.join(INPUT_FILE_LOCS, "tests_codes.txt"), "r") as infile:

        qc_test_dict = {}

        for line in infile:

            (key, val) = line.split("=")
            qc_test_dict[key.strip()] = " ".join([v.capitalize() for v in val.split()])

    return qc_test_dict # make_test_dictionary

#************************************************************************
def main():
    """
    Main plot function - no inputs.  Runs from settings in set_paths_and_vars.
    """

    qc_test_names = make_test_dictionary()

    try:
        station_info = np.genfromtxt(os.path.join(INPUT_FILE_LOCS, STATION_LIST), dtype=(str))
    except IOError:
        print "station list not found"
        sys.exit()


    all_flag_sums = np.zeros([len(station_info), len(qc_test)])
    all_flag_pct  = np.zeros([len(station_info), len(qc_test)])

    Lons = []
    Lats = []

    uk_stns = []

    for st,stat in enumerate(station_info):  

        # set up station
        station = utils.Station(stat[0], float(stat[1]), float(stat[2]), float(stat[3]))

#        if station.id[:2] != "03":
#            continue
        print st, station.id


        # read attributes and qc_flags
        try:
            ncdfp.read(os.path.join(NETCDF_DATA_LOCS, "hadisd.{}_19310101-{}_{}.nc".format(LONG_VERSION, END_TIME, station.id)), station, process_vars, diagnostics = diagnostics)

            # sum qc_flags:
            # remove multi-level flagging
            qc_flags = station.qc_flags[:]

            qc_flags[qc_flags[:] > 1] = 1

            # remove multi-level flagging - neighbour flags
            no_neighbours = qc_flags[qc_flags[:] == -1].size
            qc_flags[qc_flags[:] < 0] = 0

            total_flags = qc_flags[qc_flags[:] != 0].size

            sum_flags = np.sum(qc_flags[:], axis = 0) # 71 column array

            for cols in [strT_QC, strD_QC, strWS_QC, strWD_QC, strS_QC, T_QC, D_QC, S_QC, WS_QC, WD_QC, C_QC]:

                # to prevent double counting of flags on individual time stamps, re-sum
                combined_flags = np.sum(np.max(qc_flags[:,cols], axis = 1))

                sum_flags = np.append(sum_flags, combined_flags)

            all_flag_sums[st] = sum_flags

            # now do percentage flagged of total obs

            pct_flag = np.zeros(len(qc_test), dtype = float)

            for t,test in enumerate(qc_test):

                if t in T_QC:
                    if station.temperatures.data.compressed().size > 0:
                        pct_flag[t] = sum_flags[t] / station.temperatures.data.compressed().size
                elif t in D_QC:
                    if station.dewpoints.data.compressed().size > 0:
                        pct_flag[t] = sum_flags[t] / station.dewpoints.data.compressed().size
                elif t in S_QC:
                    if station.slp.data.compressed().size > 0:
                        pct_flag[t] = sum_flags[t] / station.slp.data.compressed().size
                elif t in WS_QC:
                    if station.windspeeds.data.compressed().size > 0:
                        pct_flag[t] = sum_flags[t] / station.windspeeds.data.compressed().size
                elif t in WD_QC:
                    if station.winddirs.data.compressed().size > 0:
                        pct_flag[t] = sum_flags[t] / station.winddirs.data.compressed().size
                elif t in C_QC:
                    if station.total_cloud_cover.data.compressed().size > 0:
                        pct_flag[t] = sum_flags[t] / station.total_cloud_cover.data.size
                else:
                    if station.temperatures.data.compressed().size > 0:
                        pct_flag[t] = sum_flags[t] / station.temperatures.data.size

            all_flag_pct[st] = 100. * pct_flag

            # get occasions when more locations are flagged than have data.

            over_100, = np.where(all_flag_pct[st] > 100.)
            all_flag_pct[st][over_100] = 100.

            Lons += [station.lon]
            Lats += [station.lat]
            uk_stns += [st]
        except RuntimeError:
            # file doesn't exist
            pass



    Lats = np.array(Lats)
    Lons = np.array(Lons)

    outfile = file(INPUT_FILE_LOCS+"all_fails_summary_{}.dat".format(start_time_string),'w')

    for t,test in enumerate(qc_test):

        plt.figure(figsize=(8,6))
        plt.clf()
        ax = plt.axes([0,0,1,1],projection = ccrs.Robinson())
        ax.set_global()
        ax.coastlines('50m')
        try:
            ax.gridlines(draw_labels = True)
        except TypeError:
            ax.gridlines()

        # colors are the exact same RBG codes as in IDL
        colors = [(150,150,150),(41,10,216),(63,160,255),(170,247,255),(255,224,153),(247,109,94),(165,0,33),(0,0,0)]
        limits = [0.0, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 100.]

        all_locs = []

        for u, upper in enumerate(limits):

            if u == 0:
                locs, = np.where(all_flag_pct[uk_stns,t] == 0)
                label = "{}%: {}".format(upper, len(locs))
            else:
                locs, = np.where(np.logical_and(all_flag_pct[uk_stns,t] <= upper, all_flag_pct[uk_stns,t] > limits[u-1]))
                label = ">{} to {}%: {}".format(limits[u-1], upper, len(locs))
                if upper == limits[-1]:
                    label = ">{}%: {}".format(limits[u-1], len(locs))


            if len(locs) > 0:
                ax.scatter(Lons[locs], Lats[locs], transform = ccrs.Geodetic(), s = 15, c = tuple([float(c)/255 for c in colors[u]]), edgecolors="none", label = label)

            else:
                ax.scatter([0], [-90], transform = ccrs.Geodetic(), s = 15, c = tuple([float(c)/255 for c in colors[u]]), edgecolors="none", label = label)

            all_locs += [len(locs)]

        plt.title(qc_test_names[test])
        watermarkstring="/".join(os.getcwd().split('/')[4:])+'/'+os.path.basename( __file__ )+"   "+dt.datetime.strftime(dt.datetime.now(), "%d-%b-%Y %H:%M")
        plt.figtext(0.01,0.01,watermarkstring,size=5)

        leg=plt.legend(loc='lower center',ncol=4, bbox_to_anchor = (0.5,-0.2), frameon=False, title='',prop={'size':11},labelspacing=0.15,columnspacing=0.5, numpoints=1)

        plt.savefig(IMAGE_LOCS+"All_fails_{}_{}.png".format(test, start_time_string))
        plt.close()


        outfile.write("{:10s}".format(test)+''.join(['%7i' % n for n in all_locs])+''.join(["%7.1f" % n for n in [100.*n/len(Lats) for n in all_locs]])+"\n")

    outfile.close()
              
    return # main

#************************************************************************
if __name__=="__main__":
    main()
#*******************************************************
