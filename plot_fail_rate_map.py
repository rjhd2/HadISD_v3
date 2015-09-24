#!/usr/local/sci/bin/python
#*****************************
#
# plot failure rate for each test and overall for variables
#
#
#************************************************************************
#                    SVN Info
#$Rev:: 78                                            $:  Revision of last commit
#$Author:: rdunn                                      $:  Author of last commit
#$Date:: 2015-07-15 16:49:32 +0100 (Wed, 15 Jul 2015) $:  Date of last commit
#************************************************************************

import numpy as np
import datetime as dt
import os
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

import qc_utils as utils
import netcdf_procs as ncdf


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
	  'SVR','WVR','WSL','WDL','WRS','STR_T','STR_D','STR_w','STR_S','ALL_T','ALL_Td','ALL_SLP','ALL_W','ACL']

T_QC=[0,1,4,5,8,12,16,20,24,27,41,44,54,58]
D_QC=[0,2,4,6,8,9,13,17,21,25,28,30,31,32,42,45,55,59]
S_QC=[0,3,4,7,11,15,19,23,26,29,43,46,57,60]
W_QC=[10,14,18,22,47,56,61,62,63,64]
C_QC=range(33,41)

strT_QC=[12,16,20]
strD_QC=[13,17,21]
strW_QC=[14,18,22]
strS_QC=[15,19,23]

station_list = "candidate_stations.txt"
process_vars = ["temperatures","dewpoints","slp","windspeeds","total_cloud_cover"]
diagnostics = False
start_time_string = dt.datetime.strftime(dt.datetime.now(), "%Y%m%d")

try:
    station_info = np.genfromtxt(os.path.join(FILE_LOCS, station_list), dtype=(str))
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

    for cols in [strT_QC, strD_QC, strW_QC, strS_QC, T_QC, D_QC, S_QC, W_QC, C_QC]:

        sum_flags = np.append(sum_flags, np.sum(sum_flags[cols]))

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
        elif t in W_QC:
            if station.windspeeds.data.compressed().size > 0:
                pct_flag[t] = sum_flags[t] / station.windspeeds.data.compressed().size
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

    

Lats = np.array(Lats)
Lons = np.array(Lons)

outfile = file(FILE_LOCS+"all_fails_summary_{}.dat".format(start_time_string),'w')

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

    if test != "ALL_SLP":
        plt.title(" ".join([s.capitalize() for s in test.split("_")]))
    else:
        plt.title("All SLP")
    watermarkstring="/".join(os.getcwd().split('/')[4:])+'/'+os.path.basename( __file__ )+"   "+dt.datetime.strftime(dt.datetime.now(), "%d-%b-%Y %H:%M")
    plt.figtext(0.01,0.01,watermarkstring,size=5)

    leg=plt.legend(loc='lower center',ncol=4, bbox_to_anchor = (0.5,-0.2), frameon=False, title='',prop={'size':11},labelspacing=0.15,columnspacing=0.5, numpoints=1)

    plt.savefig(IMG_LOCS+"All_fails_{}_{}.png".format(test, start_time_string))
    plt.close()


    outfile.write("{:10s}".format(test)+''.join(['%7i' % n for n in all_locs])+''.join(["%7.1f" % n for n in [100.*n/len(Lats) for n in all_locs]])+"\n")

outfile.close()
              
