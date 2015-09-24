#!/usr/local/sci/bin/python
#*****************************
#
# Clean Up (CLU)
#
#   Check for excessive flagging or very low obs numbers on a monthly basis
#
#
#************************************************************************
#                    SVN Info
#$Rev:: 62                                            $:  Revision of last commit
#$Author:: rdunn                                      $:  Author of last commit
#$Date:: 2015-04-30 17:03:29 +0100 (Thu, 30 Apr 2015) $:  Date of last commit
#************************************************************************
import numpy as np
import scipy as sp
import datetime as dt
# RJHD routines
import qc_utils as utils

#************************************************************************
def cu_plots(times, indata, start, end, datastart, title, extra_text = ""):
    '''
    Plot each set of values highlighted by Clean Up Check

    :param array times: array of times (hours since)
    :param array indata: data to plot
    :param int start: start of flagging period
    :param int end : end of flagging period
    :param datetime datastart: start of dataset
    :param string title: title of plot
    :param string extra_text: more text for title
    :returns:
    '''

    YLABELS = {"temperatures":"Temperature (C)", "dewpoints":"Dewpoints (C)", "slp":"SLP (hPa)", "windspeeds":"Wind Speed (m/s)"}

    extra = 480

    plot_times = utils.times_hours_to_datetime(times[start-extra:end+extra], datastart)

    import matplotlib.pyplot as plt
    plt.clf()
    plt.plot(plot_times, indata[start-extra:end+extra], 'bo', ls='-')
    plt.plot(plot_times[extra:-extra], indata[start:end], 'ro', markersize=10)
    plt.xlim([plot_times[0], plot_times[-1]])
    plt.title("Clean Up - "+title.capitalize()+" - "+extra_text)
    plt.ylabel(YLABELS[title])
    plt.show()

    return # sc_plots

#*******************************************************
def clean_up(st_var, flags, input_flag_cols, out_flag_col, start, end, times, plots = False):
    '''
    Clean up the remaining observations if many flagged or few left in a month

    :param MetVar st_var: input station object
    :param array flags: QC flags array
    :param array input_flag_cols: which columns to check over
    :param int out_flag_col: in which column to set the flags
    :param datetime start: start of dataset
    :param datetime end: end of dataset
    :param array times: hourly time stamps
    :param bool plots: show plots
    '''
    

    month_ranges = utils.month_starts_in_pairs(start, end)

    total_flags = np.sum(flags[:, input_flag_cols], axis = 1)

    filtered = utils.apply_filter_flags(st_var)

    # test each month
    for month in month_ranges:

        this_month = filtered[month[0]:month[1]]


        # if less than 20 obs, then flag remaining
        if len(this_month.compressed()) < 20:
            
            locs = np.where(this_month.mask == False)[0]

            # only flag those observations that actually exist.
            if len(locs) > 0:

                month_range = np.arange(month[0],month[1], dtype=("int"))

                flags[month_range[locs], out_flag_col] = 1
                
                if plots:
                    cu_plots(times, st_var.data, month_range[0], month_range[-1], start, st_var.name, extra_text = "few_obs")

        
        # if 40% of obs flagged, then flag remainig
        else:

            good_locs = np.where(this_month.mask == False)[0]
            flag_locs = np.where(total_flags[month[0]:month[1]] > 0)[0]

            # good_locs - internal and neighbour flags already applied, so need to add
            #      these back in.
            proportion = float(len(flag_locs))/float(len(good_locs)+len(flag_locs))

            if proportion > 0.4:

                month_range = np.arange(month[0],month[1], dtype=("int"))
                flags[month_range[good_locs], out_flag_col] = 1 

                if plots:
                    cu_plots(times, st_var.data, month_range[0], month_range[-1], start, st_var.name, extra_text = "lots flagged")

                
    return # clean_up

#*******************************************************
def clu(station, var_list, flag_cols, FLAG_COL_DICT, start, end, logfile, plots = False):
    '''
    Run the clean up for each variable

    :param file logfile: logfile to store outputs
    '''


    for v,variable in enumerate(var_list):

        st_var = getattr(station, variable)

        clean_up(st_var, station.qc_flags, FLAG_COL_DICT[variable], flag_cols[v], start, end, station.time.data, plots = plots)

        flag_locs = np.where(station.qc_flags[:, flag_cols[v]] != 0)

        if plots:
            utils.print_flagged_obs_number(logfile, "Clean Up Months", variable, len(flag_locs[0]), noWrite = True)
        else:
            utils.print_flagged_obs_number(logfile, "Clean Up Months", variable, len(flag_locs[0]))


        # copy flags into attribute
        st_var.flags[flag_locs] = 1

    station = utils.append_history(station, "Clean Up Months")
                     
    return # clu

#************************************************************************
if __name__ == "__main__":

    print "Checking for sparsely observed or highly flagged months"
