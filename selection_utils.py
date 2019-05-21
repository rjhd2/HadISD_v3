#!/usr/local/sci/bin/python2.7
#*****************************
#
# general utilities for station selection scripts
#
#
#************************************************************************
#                    SVN Info
#$Rev:: 117                                           $:  Revision of last commit
#$Author:: rdunn                                      $:  Author of last commit
#$Date:: 2017-01-30 15:33:46 +0000 (Mon, 30 Jan 2017) $:  Date of last commit
#************************************************************************

import numpy as np

# RJHD utilities
import qc_utils as utils

# merging limits
DISTANCE_THRESHOLD = 25. # km - 1/e drop in probability by 25km sep
ELEVATION_THRESHOLD = 100. # m - 1/e drop in probability by 100m sep
LATITUDE_THRESHOLD = 2 # deg
PROB_THRESHOLD = 0.5



#****************************************************
def jaccard(seq1, seq2):
    """
    Compute the Jaccard distance between the two sequences `seq1` and `seq2`.
    
    They should contain hashable items.
	
    The return value is a float between 0 and 1, where 1 means equal, and 0 totally different.
    """
    set1, set2 = set(seq1), set(seq2)
    return len(set1 & set2) / float(len(set1 | set2)) # jaccard


#****************************************************
def do_match(station1, station2, latitude, elevation, distance):
    """
    Perform the match between two stations.  

    Do initial latitude check to speed up the test
    (not longitude as this isn't a constant distance)

    Return probabilities for elevation, separation and Jaccard Index

    :param Station Class station1: 
    :param Station Class station2: 
    :returns:
       list of 3 probabilities [elev, dist, jaccard]
    """


    # latitude - pre check to make quicker
    if np.abs(station1.lat - station2.lat) > LATITUDE_THRESHOLD:
        return False

    # elevation
    height = np.abs(station1.elev - station2.elev)
    if height < (ELEVATION_THRESHOLD*4):
        height_Pr = np.exp(-1.0 * height / ELEVATION_THRESHOLD)
    else:
        height_Pr = 0

    # latitude & longitude
    distance, bearing = utils.get_dist_and_bearing([station1.lat, station1.lon],[station2.lat, station2.lon])
    if distance < (DISTANCE_THRESHOLD*4):
        dist_Pr = np.exp(-1.0 * distance / DISTANCE_THRESHOLD)
    else:
        dist_Pr = 0.

    # Jaccard Index on name - remove all whitespace
    jac_Pr = jaccard(station1.name.strip(), station2.name.strip())

    # Jaccard Index on METAR call sign
    if station1.call != "" and station2.call != "":
        jac_Pr_metar = jaccard(station1.call, station2.call)
    

    # name matching
    return [height_Pr, dist_Pr, jac_Pr] # do_match

#-------------------------------------------------
#                END
#-------------------------------------------------
