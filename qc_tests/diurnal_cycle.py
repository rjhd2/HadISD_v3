#!/usr/local/sci/bin/python
#*****************************
#
# Diurnal Cycle Check (DCC)
#
#   At times this is a direct translation from IDL
#     Could be made more pythonic, but need to match outputs!
#
#
#************************************************************************
#                    SVN Info
#$Rev:: 67                                            $:  Revision of last commit
#$Author:: rdunn                                      $:  Author of last commit
#$Date:: 2015-05-01 16:18:52 +0100 (Fri, 01 May 2015) $:  Date of last commit
#************************************************************************

import numpy as np
import scipy as sp
import datetime as dt

# RJHD routines
import qc_utils as utils


OBS_PER_DAY = 4.
DAILY_RANGE = 5.
HOURS = 24
INTMDI = -999

DYNAMIC_DIURNAL = True


#************************************************************************
def dcc_quartile_check(day):
    '''
    Check if >=3 quartiles of the day have data
    
    :param array day: 24 hours of data in array
    :returns: boolean
    '''
    
    quartile_has_data = np.zeros(4)

    quart_start = np.arange(0, HOURS, HOURS / 4)
    quart_end   = quart_start + 6

    for q in range(4):
        if len(day[quart_start[q]:quart_end[q]].compressed()) > 0:
            quartile_has_data[q] = 1

    if quartile_has_data.sum() >= 3:
        return True

    else:
        return False # dcc_quartile_check

#************************************************************************
def dcc_make_sine():
    '''Return sine curve over 24 points '''
    return np.sin(2.*np.pi* np.arange(HOURS,dtype=(np.float64)) / HOURS) # dcc_make_sine



#************************************************************************
def dcc(station, variable_list, full_variable_list, flag_col, logfile, plots = False, diagnostics = False):
    '''
    The diurnal cycle check.
    
    :param object station: the station object to be processed
    :param list variable_list: the variables to be processed
    :param list full_variable_list: the variables for flags to be applied to
    :param list flag_col: which column in the qc_flags array to work on
    :param file logfile: logfile to store outputs
    :param bool plots: to do any plots
    :param bool diagnostics: to do any extra diagnostic output
    :returns:
    '''

    # list of flags for each variable
    diurnal_flags = []

    for v,variable in enumerate(variable_list):
    
        st_var = getattr(station, variable)
 
 	# is this needed 21/08/2014        
#        reporting_accuracy = utils.reporting_accuracy(utils.apply_filter_flags(st_var))
        
        # apply flags - for detection only
        filtered_data = utils.apply_filter_flags(st_var)

        filtered_data = filtered_data.reshape(-1,24) # working in fulltimes.
        number_of_days = filtered_data.shape[0]

        if plots:
            import matplotlib.pyplot as plt
            plt.clf()
            plot_data = np.ma.zeros(filtered_data.shape)
            plot_data.mask = True
#            best_estimate_counter = np.zeros(HOURS)

        diurnal_best_fits     = np.zeros(filtered_data.shape[0], dtype = (int))
        diurnal_best_fits.fill(INTMDI)
        diurnal_uncertainties = np.zeros(filtered_data.shape[0])
        diurnal_uncertainties.fill(INTMDI)

        for d,day in enumerate(filtered_data):

            '''enough observations and have large enough diurnal range '''
            if len(day.compressed()) >= OBS_PER_DAY:

                obs_daily_range = max(day.compressed()) - min(day.compressed())
                if obs_daily_range >= DAILY_RANGE:
                    
                    if dcc_quartile_check(day):
                        scaled_sine = ((dcc_make_sine() + 1.) / 2. * obs_daily_range) + min(day.compressed())
                        diffs = np.zeros(HOURS)

                        '''Find differences for each shifted sine --> cost function'''
                        for h in range(HOURS):
                            diffs[h] = np.sum(np.abs(day - scaled_sine).compressed())
                            scaled_sine = np.roll(scaled_sine, 1) # matched to IDL SHIFT()
                                                                           
                        diurnal_best_fits[d] = np.argmin(diffs)

                        # default uncertainty is the average time resolution of the data
                        diurnal_uncertainties[d] = round(float(HOURS) / len(day.compressed()))

                        if DYNAMIC_DIURNAL:
                            critical_value = min(diffs) + ((max(diffs) - min(diffs)) * 0.33)

                            # centre so minimum in middle
                            diffs = np.roll(diffs, 11 - diurnal_best_fits[d])
                            
                            uncertainty = 1
                            while uncertainty < 11:
                                if (diffs[11 - uncertainty] > critical_value) and\
                                        (diffs[11 + uncertainty] > critical_value):
                                    # break if both sides greater than critical difference
                                    # when counting outwards
                                    #    see diurnal_example.py
                                    break

                                uncertainty += 1

                            # check if uncertainty greater than time resolution for day
                            if uncertainty > diurnal_uncertainties[d] :
                                diurnal_uncertainties[d] = uncertainty

                        if plots:
#                            best_estimate_counter[np.argmin(diffs)] += 1
                            # scale daily data to range -1 -> 1, plot with random scatter for clarity
                            plot_data[d] = ((2 * (day - min(day.compressed())) / obs_daily_range) - 1.)
                            plt.plot(np.arange(24)+np.random.randn(24)*0.25, plot_data[d]+np.random.randn(24)*0.05, 'k,')
                            

        if plots:
            plt.plot(np.arange(24),np.roll(dcc_make_sine(), np.argmax(np.bincount(diurnal_best_fits[np.where(diurnal_best_fits != INTMDI)]))),'r-')
            plt.xlim([-1,25])
            plt.ylim([-1.2,1.2])
            plt.show()

        # dumb copy of IDL

        '''For each uncertainty range (1-6h) find median of cycle offset'''
        best_fits = np.zeros(6)
        best_fits.fill(-9)
        for h in range(6):
            locs = np.where(diurnal_uncertainties == h+1)

            if len(locs[0]) > 300:
                # best_fits[h] = int(np.median(diurnal_best_fits[locs])) 
                # Numpy median gives average of central two values which may not be integer
                # 25/11/2014 use IDL style which gives lower value
                best_fits[h] = utils.idl_median(diurnal_best_fits[locs])
 
        '''Build up range of cycles incl, uncertainty to find where best of best located'''

        hours = np.arange(24)
        hour_matches=np.zeros(24)
        diurnal_peak = -9
        number_estimates = 0
        for h in range(6):
            if best_fits[h] != -9:

                '''Store lowest uncertainty best fit as first guess'''
                if diurnal_peak == -9: 
                    diurnal_peak = best_fits[h]
                    hours = np.roll(hours,11-int(diurnal_peak))
                    hour_matches[11-(h+1):11+(h+2)] = 1
                    number_estimates += 1
                 
                centre = np.where(hours == best_fits[h])
 
                if (centre[0] - h + 1) >= 0:
                    if (centre[0] + h + 1 ) <=23:
                        hour_matches[centre[0] - (h + 1) : centre[0] + h + 2] += 1
                    else:
                        hour_matches[centre[0] - (h + 1) : ] += 1
                        hour_matches[ : centre[0] + h + 2- 24] += 1                                        
                else:
                    hour_matches[: centre[0] + h + 2] += 1
                    hour_matches[centre[0] - (h + 1) :] += 1

                number_estimates += 1

        
        '''If value at lowest uncertainty not found in all others, then see what value is found by all others '''
        if hour_matches[11] != number_estimates:  # central estimate at 12 o'clock
            all_match = np.where(hour_matches == number_estimates)

            # if one is, then use it
            if len(all_match[0]) > 0:
                diurnal_peak = all_match[0][0]
            else:
                diurnal_peak = -9
            
 
        '''Now have value for best fit diurnal offset'''

        potentially_spurious = np.zeros(number_of_days)
        potentially_spurious.fill(INTMDI)

        if diurnal_peak != -9:
            hours = np.arange(24)
            hours = np.roll(hours,11-int(diurnal_peak))
            for d in range(number_of_days):
                if diurnal_best_fits[d] != INTMDI:

                    '''Checks if global falls inside daily value+/-range
                    rather than seeing if each day falls in global value+/-range'''

                    

                    min_range = 11 - diurnal_uncertainties[d]
                    max_range = 11 + diurnal_uncertainties[d]
                    maxloc = np.where(hours == diurnal_best_fits[d])[0][0]


                    if maxloc < min_range or maxloc > max_range:
                        potentially_spurious[d] = 1
                    else:
                        potentially_spurious[d] = 0
                    

 

            # count number of good, missing and not-bad days
            n_good = 0
            n_miss = 0
            n_not_bad = 0
            total_points = 0
            total_not_miss = 0
            to_flag = np.zeros(number_of_days)

            for d in range(number_of_days):

                if potentially_spurious[d] == 1:
                    
                    n_good = 0
                    n_miss = 0
                    n_not_bad = 0
                    total_points += 1
                    total_not_miss +=1
                    
                else:

                    if potentially_spurious[d] == 0:

                        n_good += 1
                        n_not_bad += 1
                        if n_miss != 0:
                            n_miss = 0
                        total_not_miss += 1

                    if potentially_spurious[d] == -999:

                        n_miss += 1
                        n_not_bad += 1
                        if n_good != 0:
                            n_good = 0

                    total_points += 1


                    if (n_good == 3) or (n_miss == 3) or (n_not_bad >=6):

                        
                        if total_points >= 30:
                            if float(total_not_miss)/total_points >= 0.5:
                                to_flag[d - total_points : d ] = 1
                        
                        n_good = 0
                        n_miss = 0
                        n_not_bad = 0
                        total_points = 0 
                        total_not_miss = 0

            dcc_flags = np.zeros(filtered_data.shape)

            for d in range(number_of_days):

                if to_flag[d] == 1:
                    good = np.where(filtered_data.mask[d,:] == False)
                    if len(good[0]) >= 1:
                        dcc_flags[d,good]=1

            if diagnostics:
                print len(np.where(dcc_flags == 1)[0])
                print "currently matches IDL, but should all hours in days have flags set, not just the missing/flagged ones?"

            diurnal_flags += [dcc_flags]
        else: 
            diurnal_flags += [np.zeros(filtered_data.shape)]

        station.qc_flags[:, flag_col[v]] = np.array(diurnal_flags).reshape(-1)

        flag_locs = np.where(station.qc_flags[:, flag_col[v]] != 0)
        if plots or diagnostics:
            utils.print_flagged_obs_number(logfile, "Diurnal Cycle", variable, len(flag_locs[0]), noWrite = True)
        else:
            utils.print_flagged_obs_number(logfile, "Diurnal Cycle", variable, len(flag_locs[0]))

        # copy flags into attribute
        st_var.flags[flag_locs] = 1

        # CHECKED 030660-99999, 30-06-2014, 855 flagged RJHD
    
    utils.apply_flags_all_variables(station, full_variable_list, flag_col[variable_list == "temperatures"], logfile, "Diurnal Cycle", plots = plots, diagnostics = diagnostics)

    station = utils.append_history(station, "Diurnal Cycle Check")  
                     
    return # dcc


#************************************************************************

if __name__ == "__main__":
    
    print "checking diurnal cycle"
