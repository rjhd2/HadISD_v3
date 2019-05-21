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
import os
import shutil
import glob
import ConfigParser
import sys
import json

#************************************************************************
#
# Config file held in hadobs homespace for the moment
#   + it's group writable
#   + it saves a major recode of all scripts to take a --config option
#   - there's a hard coded file path
#
# Similarly with json file for parameters to run
#
# Later it may be reasonable to add the --config option to all scripts
#   and re-order this script to work with that
#
# RJHD October 2018
#
#************************************************************************

CONFIG_FILE = "~hadobs/hadisd/input_files_v3/configuration.txt"

if not os.path.exists(os.path.expanduser(CONFIG_FILE)):
    print "Configuration file missing - {}".format(os.path.expanduser(CONFIG_FILE))
    sys.exit

#************************************************************************
def copyfiles(src, dest):
    print r'{}*'.format(os.path.expanduser(src))
    for filename in glob.glob(r'{}*'.format(os.path.expanduser(src))):
        shutil.copy(filename, dest)
        print filename, dest
    return # copyfiles

#************************************************************************
def expand_version(ver, end_time, ceda = False):
    '''
    Process the short version (e.g. v100_2011f) into the longer form and an end date

    :param str ver: short version
    
    :returns: long_ver and end_str
    '''

    last = ver.split("_")[-1]

    long_ver = "{}.{}.{}.{}".format(ver[1], ver[2], ver[3], last)
    
    # need to be cleverer when go to monthly!
    end_str = dt.datetime.strftime(end_time, "%Y%m%d")
 
    return long_ver, end_str # expand_version
#************************************************************************

# read in configuration file
config = ConfigParser.ConfigParser()
config.read(os.path.expanduser(CONFIG_FILE))

# set Dates
DATASTART = dt.datetime(config.getint("Dates","startyear"), config.getint("Dates", "startmonth"), config.getint("Dates", "startday"), 0, 0)
DATAEND = dt.datetime(config.getint("Dates","endyear"), config.getint("Dates", "endmonth"), config.getint("Dates", "endday"), 0, 0)

# set station lists
STATION_LIST = config.get("Files", "station_list")
MERGER_LIST = config.get("Files", "merger_list")
PARAM_FILE = config.get("Files", "parameters")

# set versions
HADISD_VERSION = config.get("Versions", "hadisd_version")
PREVIOUS_VERSION = config.get("Versions", "previous_version")
LONG_VERSION, END_TIME = expand_version(HADISD_VERSION, DATAEND)

# set paths
ROOT = config.get("Paths", "root")
ANCILS_LOC = config.get("Paths", "ancillaries")


#************************************************************************
# make all the other paths, and files
if not os.path.exists(ROOT): os.mkdir(ROOT)
ROOT_LOC = "/{}/{}/".format(ROOT, HADISD_VERSION)

# Set up the path variables and create if necessary

# set of files for running the suite - some static, others dynamic
INPUT_FILE_LOCS = "{}/input_files_{}/".format(ROOT_LOC, HADISD_VERSION)
if not os.path.exists(INPUT_FILE_LOCS): 
    os.mkdir(INPUT_FILE_LOCS)

# Copy files if not already there.
# for attributes.dat, test_codes.txt, Canadian files etc
copyfiles(ANCILS_LOC, INPUT_FILE_LOCS)
# remainder are created during the run

# the images produced during the suite
IMAGE_LOCS = "{}/img_files_{}/".format(ROOT_LOC, HADISD_VERSION)
if not os.path.exists(IMAGE_LOCS): os.mkdir(IMAGE_LOCS)

# the processed netcdf data
NETCDF_DATA_LOCS = "{}/netcdf_files_{}/".format(ROOT_LOC, HADISD_VERSION)
if not os.path.exists(NETCDF_DATA_LOCS): os.mkdir(NETCDF_DATA_LOCS)

# and the data for CEDA
CLIPC_DATA_LOCS = "{}/clipc_files_{}/".format(ROOT_LOC, HADISD_VERSION)
if not os.path.exists(CLIPC_DATA_LOCS): os.mkdir(CLIPC_DATA_LOCS)
if not os.path.exists(CLIPC_DATA_LOCS+"/esgf"): os.mkdir(CLIPC_DATA_LOCS+"/esgf")
if not os.path.exists(CLIPC_DATA_LOCS+"/ceda"): os.mkdir(CLIPC_DATA_LOCS+"/ceda")

# the raw ISD data
ISD_DATA_LOCS = "{}/isd_files_{}/".format(ROOT_LOC, HADISD_VERSION)
if not os.path.exists(ISD_DATA_LOCS): os.mkdir(ISD_DATA_LOCS)

# logs from the QC suite
LOG_OUTFILE_LOCS = "{}/suppl_files_{}/".format(ROOT_LOC, HADISD_VERSION)
if not os.path.exists(LOG_OUTFILE_LOCS): os.mkdir(LOG_OUTFILE_LOCS)

# A copy of the QC suite and associated files for this version
CODE_LOCS = "{}/code_{}/".format(ROOT_LOC, HADISD_VERSION) # to store working copy of files
if not os.path.exists(CODE_LOCS): os.mkdir(CODE_LOCS)

# For finding already downloaded files
OLD_ISD_DATA_LOCS = "/scratch/rdunn/hadisd/{}/isd_files_{}/".format(PREVIOUS_VERSION, PREVIOUS_VERSION)
OLD_INPUT_FILE_LOCS = "/scratch/rdunn/hadisd/{}/input_files_{}/".format(PREVIOUS_VERSION, PREVIOUS_VERSION)

# copy candidate_stations.txt to new version.  Will be overwritten if updated station selection
for cfile in [STATION_LIST, MERGER_LIST]:
    if not os.path.exists("{}/{}".format(INPUT_FILE_LOCS, cfile)): 
        print "copying {}/{} to {}".format(OLD_INPUT_FILE_LOCS, cfile, INPUT_FILE_LOCS)
        shutil.copy("{}/{}".format(OLD_INPUT_FILE_LOCS, cfile), INPUT_FILE_LOCS )

# copy Canada files to new version.  Again, overwritten if Canada re-run
for cfile in glob.glob("{}/Canada*".format(OLD_INPUT_FILE_LOCS)):
    cfilename = cfile.split("/")[-1]
    if not os.path.exists("{}/{}".format(INPUT_FILE_LOCS, cfilename)): 
        print "copying {}/{} to {}".format(OLD_INPUT_FILE_LOCS, cfilename, INPUT_FILE_LOCS)
        shutil.copy("{}/{}".format(OLD_INPUT_FILE_LOCS, cfilename), INPUT_FILE_LOCS )
    
#************************************************************************
# variables
with open(os.path.join(os.path.expanduser(ANCILS_LOC),PARAM_FILE), "r") as pf:
    parameters = json.load(pf)

process_vars = parameters["variables"]["process_vars"]
carry_thru_vars = parameters["variables"]["carry_thru_vars"]


#************************************************************************
# print for information each time - enables clearer checking
print "HadISD version: {}".format(HADISD_VERSION)
print "Data location : {}".format(NETCDF_DATA_LOCS)
print "Data range    : {} - {}\n".format(dt.datetime.strftime(DATASTART, "%Y-%m-%d"), dt.datetime.strftime(DATAEND, "%Y-%m-%d"))


#------------------------------------------------------------
# END
#------------------------------------------------------------
