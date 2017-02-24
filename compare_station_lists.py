#!/usr/local/sci/bin/python2.7
#*****************************
#
# Compares current and previous station lists to give summary of changes.
#
#
#************************************************************************
#                    SVN Info
#$Rev:: 92                                            $:  Revision of last commit
#$Author:: rdunn                                      $:  Author of last commit
#$Date:: 2016-02-08 12:00:09 +0000 (Mon, 08 Feb 2016) $:  Date of last commit
#************************************************************************


import numpy as np
import os

# RJHD utils
from set_paths_and_vars import *


#**************************************
def WriteFile(data, outfilename):
    """
    Write the output files

    :param array data: data array of station IDs
    :param str outfilename: output file name
    :returns:
    """

    with file(INPUT_FILE_LOCS + outfilename, "w") as outfile:
        for stn in data:
            outfile.write("{}\n".format(stn))

    return # WriteFile

#**************************************
#**************************************

station_list = "candidate_stations.txt"

try:
    new_station_info = np.genfromtxt(os.path.join(INPUT_FILE_LOCS, station_list), dtype=(str))
except IOError:
    print "new station list not found"
    sys.exit()

try:
    old_station_info = np.genfromtxt(os.path.join(OLD_INPUT_FILE_LOCS, station_list), dtype=(str))
except IOError:
    print "new station list not found"
    sys.exit()


new_ids = new_station_info[:,0]
old_ids = old_station_info[:,0]


in_both = new_ids[np.in1d(new_ids, old_ids)]
in_new = new_ids[np.in1d(new_ids, old_ids, invert = True)]
in_old = old_ids[np.in1d(old_ids, new_ids, invert = True)]


WriteFile(in_new, "candidate_stations_additions.dat")
WriteFile(in_old, "candidate_stations_removals.dat")
WriteFile(in_both, "candidate_stations_continued.dat")

#************************************************************************
