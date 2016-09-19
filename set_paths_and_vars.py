#!/usr/local/sci/bin/python2.7
#------------------------------------------------------------
#                    SVN Info
#$Rev:: 84                                            $:  Revision of last commit
#$Author:: rdunn                                      $:  Author of last commit
#$Date:: 2015-12-18 16:35:07 +0000 (Fri, 18 Dec 2015) $:  Date of last commit
#------------------------------------------------------------
# START
#------------------------------------------------------------
"""
Sets a load of paths and defaults.

Should be the only place to edit these each time around.
"""
import datetime as dt

# File paths to read and use

HADISD_VERSION = "v200_2015p"
PREVIOUS_VERSION = "v200_2015p"


ROOT_LOC = "/project/hadobs2/hadisd/{}".format(HADISD_VERSION)

INPUT_FILE_LOCS = "{}/code_{}/input_files/".format(ROOT_LOC, HADISD_VERSION)
IMAGE_LOCS = "{}/img_files_{}/".format(ROOT_LOC, HADISD_VERSION)
NETCDF_DATA_LOCS = "{}/netcdf_files_{}/".format(ROOT_LOC, HADISD_VERSION)
ISD_DATA_LOCS = "{}/isd_files_{}/".format(ROOT_LOC, HADISD_VERSION)
LOG_OUTFILE_LOCS = "{}/suppl_files_{}/".format(ROOT_LOC, HADISD_VERSION)

OLD_ISD_DATA_LOCS = "/project/hadobs2/hadisd/{}/isd_files_{}/".format(PREVIOUS_VERSION, PREVIOUS_VERSION)
OLD_INPUT_FILE_LOCS = "/project/hadobs2/hadisd/{}/code_{}/input_files/".format(PREVIOUS_VERSION, PREVIOUS_VERSION)


# Other settings

DATASTART = dt.datetime(1931,1,1,0,0)
DATAEND = dt.datetime(2016,1,1,0,0)

process_vars = ["temperatures","dewpoints","slp","windspeeds", "winddirs", "total_cloud_cover","low_cloud_cover","mid_cloud_cover","high_cloud_cover"]
carry_thru_vars = ["cloud_base","precip1_depth","precip1_period","wind_gust", "past_sigwx1"]


# For SPICE/Slurm
# NETCDF_DATA_LOCS = "/data/local/hadobs/hadisdv2/"

#------------------------------------------------------------
# END
#------------------------------------------------------------
