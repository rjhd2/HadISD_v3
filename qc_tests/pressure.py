#!/usr/local/sci/bin/python
#*****************************
#
# Station Pressure Check (SPC)
#   Ensure SLP and STNLP are reasonable (SLP taken as truth)
#
#************************************************************************
#                    SVN Info
#$Rev:: 219                                           $:  Revision of last commit
#$Author:: rdunn                                      $:  Author of last commit
#$Date:: 2019-05-20 16:56:47 +0100 (Mon, 20 May 2019) $:  Date of last commit
#************************************************************************
import numpy as np
import scipy as sp
import datetime as dt
# RJHD routines
import qc_utils as utils

MAD_THRESHOLD = 4.5


#*********************************************
def spc_diff(sfc, stn, flags, month_ranges, start, end, logfile, plots = False, diagnostics = False, doMonth = False):
    '''
    Pressure difference check, on individual obs.  Remove very silly stnlp
    
    :param array sfc: SLP
    :param array stn: STNLP 
    :param array flags: flags_array
    :param array month_ranges: array of month start and end times
    :param datetime start: DATASTART
    :param datetime end: DATAEND
    :param file logfile: logfile to store outputs
    :param bool plots: do plots or not
    :param bool diagnostics: extra verbose output
    :param bool doMonth: account for spare month

    :returns: flags - locations where flags have been set
    '''
    
    month_ranges = utils.month_starts_in_pairs(start, end)
    month_ranges = month_ranges.reshape(-1,12,2)

    # apply flags (and mask incomplete year if appropriate)
    sfc_filtered = utils.apply_filter_flags(sfc, doMonth = doMonth, start = start, end = end)
    stn_filtered = utils.apply_filter_flags(stn, doMonth = doMonth, start = start, end = end)

    # get the differences
    diffs = sfc.data - stn.data
    diffs_filtered = sfc_filtered - stn_filtered

    # robust statistics
    median_difference = np.ma.median(diffs)
    mad_difference = utils.mean_absolute_deviation(diffs, median = True)

    # where exceed
    high, = np.ma.where(diffs > (median_difference + MAD_THRESHOLD*mad_difference))
    low, = np.ma.where(diffs < (median_difference - MAD_THRESHOLD*mad_difference))

    # set flags
    if len(high) != 0:
        if diagnostics: print "Number of high differences {}".format(len(high))
        flags[high] = 1
    if len(low) != 0:
        if diagnostics: print "Number of low differences {}".format(len(low))
        flags[low] = 1

    if plots:
        import matplotlib.pyplot as plt
        plt.clf()
        plt.hist(diffs.compressed(), bins = np.arange(np.round(median_difference)-10, np.round(median_difference)+10, 0.1))
        plt.axvline(x = (median_difference + 4*mad_difference), ls = "--", c = "r")
        plt.axvline(x = (median_difference - 4*mad_difference), ls = "--", c = "r")
        plt.xlim([median_difference - 11, median_difference + 11])
        plt.ylabel("Observations")
        plt.xlabel("Difference (hPa)")
        plt.show()

    # How to set the range of allowable values.
    nflags, = np.where(flags != 0)
    utils.print_flagged_obs_number(logfile, "Station Level Pressure", "stnlp", len(nflags), noWrite=diagnostics)

    return flags # spc_diff

#************************************************************************
def spc(station, flag_col, start, end, logfile, diagnostics = False, plots = False, doMonth = False):
    '''
    Run pressure cross checks on SLP and STNLP

    :param object station: station object
    :param array flag_col: which columns to fill in flag_array
    :param datetime start: start of dataset
    :param datetime end: end of dataset
    :param file logfile: logfile to store outputs
    :param bool plots: do plots or not
    :param bool diagnostics: extra verbose information
    :param bool doMonth: account for incomplete months

    :returns:    
    '''

    slp = getattr(station, 'slp')
    stnlp = getattr(station, 'stnlp')

    month_ranges = utils.month_starts_in_pairs(start, end)

    station.qc_flags[:,flag_col[0]] = spc_diff(slp, stnlp, station.qc_flags[:,flag_col[0]], month_ranges, start, end, logfile, plots = plots, diagnostics = diagnostics, doMonth = doMonth)

    station = utils.append_history(station, "Pressure Cross Check")  

    return # spc

#************************************************************************
if __name__ == "__main__":
    
    print "station level pressure check" 
    
