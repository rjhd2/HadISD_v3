#!/usr/local/sci/bin/python
#*****************************
#
# Humidity Cross Check (HCC)
#   Temperature and Dewpoints consistency check
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

SSS_MONTH_FRACTION = 0.2

#************************************************************************
def hcc_time_plot(Ts, Ds, start, end, datastart):
    '''
    Plot time series of Temperature and Dewpoints 
    showing which points have been flagged from the tests

    :param array Ts: Temperatures
    :param array Ds: Dewpoints
    :param int start: start of flag period
    :param int end: end of flag period
    :param datetime datastart: start of data set

    :returns:
    '''

    extra = 48

    times = utils.times_hours_to_datetime(np.arange(start - extra, end + extra, 1), datastart)

    import matplotlib.pyplot as plt
    plt.clf()

    plt.plot(times, Ts[start - extra: end + extra], 'ko', ls = '-', label = "Temperature")
    plt.plot(times, Ds[start - extra: end + extra], 'bs', ls = '-', label = "Dewpoints")

    plt.plot(times[extra:-extra], Ts[start: end], 'ro')
    plt.plot(times[extra:-extra], Ds[start: end], 'rs')
    plt.ylabel("(Dewpoint) Temperature (C)")
    plt.legend(loc='lower center', ncol=2, frameon=False,prop={'size':13})
    plt.show()

    return # hcc_time_plot

#************************************************************************
def hcc_sss(T, D, month_ranges, start, logfile, plots = False, diagnostics = False):
    '''
    Supersaturation check, on individual obs, and then if >20% of month affected
    
    :param array T: temperatures
    :param array D: dewpoint temperatures
    :param array month_ranges: array of month start and end times
    :param datetime start: DATASTART (for plotting)
    :param file logfile: logfile to store outputs
    :param bool plots: do plots or not
    :param bool diagnostics: extra verbose output

    :returns: flags - locations where flags have been set
    '''
    
    flags = np.zeros(len(T))

    # flag each location where D > T
    for m,month in enumerate(month_ranges):

        data_count = 0.
        sss_count = 0.

        try:
            for t in np.arange(month[0],month[1]):

                data_count += 1

                if D[t] > T[t]:
                    sss_count += 1

                    flags[t] = 1

                    if plots:
                        hcc_time_plot(T, D, t-1, t, start)
        except IndexError:
            # no data for that month - incomplete year
            pass

        # test whole month
        # if more than 20% flagged, flag whole month
        if sss_count / data_count >= SSS_MONTH_FRACTION:
            
            flags[month[0]:month[1]] = 1

            if plots:
                hcc_time_plot(T, D, month[0], month[1], start)

    nflags = len(np.where(flags != 0)[0])
    utils.print_flagged_obs_number(logfile, "Supersaturation", "temperature", nflags, noWrite = diagnostics)

    # not yet tested.
    return flags # hcc_sss


#************************************************************************
def hcc_dpd(times, T, D, P, C, SX, start, logfile, plots = False, diagnostics = False):
    '''
    Dew point Depression check.  If long string of DPD = 0, then flag
    
    :param array times: timestamps
    :param array T: temperatures
    :param array D: dewpoint temperatures
    :param array P: precipitation depth arrays for 1,3,6,12,24 hourly accumul.
    :param array C: cloud base
    :param array SX: past significant weather
    :param datetime start: DATASTART (for plotting)
    :param file logfile: logfile to store outputs
    :param bool plots: do plots or not
    :param bool diagnostics: extra verbose output

    :returns: flags - locations where flags have been set
    '''
 
    flags = np.zeros(len(T))

    dpds = T - D
    
    last_dpds = -9999.
    
    string_start_time = times[0]
    start_loc = 0

   
    for t,tt in enumerate(times):
        
        if (tt > 0) and (tt < times[-1]):
            
            if (dpds.mask[t] == False):

                # if change in DPD, examine previous string
                if (dpds[t] != last_dpds):

                    # if long enough
                    if (times[t-1]-string_start_time >= 24):
                        these_dpds = dpds[start_loc : t]
                        good = np.where(these_dpds.mask == False)

                        if T[t] >= 0:
                            abs_diff = 0.25
                        else:
                            abs_diff = 1.

                        # has enough data and is small enough
                        if (len(good[0]) >=4) and (abs(last_dpds) <= abs_diff):
                            
                            # check if weather event could explain it.
                            these_sigwx = SX[start_loc : t]
                            these_P = P[:, start_loc : t] # not 1-D anymore
                            these_CB = C[start_loc : t]
                            
                            # use past significant weather, precipitation or low cloud base
                            fog = np.where(np.logical_or.reduce((these_sigwx[good] >= 4, \
                                               these_P[0][good] > 0., \
                                               these_P[1][good] > 0., \
                                               these_P[2][good] > 0., \
                                               these_P[3][good] > 0., \
                                               these_P[4][good] > 0., \
                                               these_P[5][good] > 0., \
                                               these_P[6][good] > 0., \
                                               these_P[7][good] > 0., \
                                               these_P[8][good] > 0., \
                                               np.logical_and(these_CB[good] > 0., these_CB[good] < 1000.))))
                                               # these_P[0] to [8] are the 1, 2, 3, 6, 9, 12, 15, 18 and 24 hourly accumulations

                            if len(fog[0]) >= 1:
 
                                if len(fog[0])/float(len(good[0])) < 0.333:
                                    flags[start_loc : t][good] = 1

                                    if plots:
                                        hcc_time_plot(T, D, start_loc, t, start)

                            else:
                                flags[start_loc : t][good] = 1
                                if plots:
                                    hcc_time_plot(T, D, start_loc, t, start)

                    string_start_time = tt
                    start_loc = t
                    last_dpds = dpds[t]
                      

    nflags = len(np.where(flags != 0)[0])
    utils.print_flagged_obs_number(logfile, "Dewpoint Depression", "temperature", nflags, noWrite = diagnostics)

    # checked on 032220 on 19/8/2014 and matches identically
    return flags # hcc_dpd

#************************************************************************
def hcc_cutoffs(T, D, month_ranges, logfile, start, plots = False, diagnostics = False):
    '''
    Check each month to see if most T obs have corresponding D obs
    do in bins of 10C - if one bin has <50% of match, then remove month
    
    :param array T: temperatures
    :param array D: dewpoint temperatures
    :param array month_ranges: array of month start and end times
    :param file logfile: logfile to store outputs
    :param datetime start: start of dataset for labelling of plots if required
    :param bool plots: do plots
    :param bool diagnostics: output extra verbose information

    :returns: flags - locations where flags have been set
    '''
    
    flags = np.zeros(len(T))

    binwidth = 10
    bins = np.arange(-90,70,binwidth)
    
    for m,month in enumerate(month_ranges):
        
        this_month_T = T[month[0]:month[1]]
        this_month_D = D[month[0]:month[1]]
    
        goodT = np.where(this_month_T.mask == False)
        goodD = np.where(this_month_D.mask == False)
        
        # check if more than 112 obs (4/d * 28 d)
        if len(goodT[0]) > 112:
            
            # run through each bin
            for bin in bins:
                
                bin_locs = np.where((this_month_T >= bin) & (this_month_T < bin+binwidth))
                
                # if data in this bin
                if len(bin_locs[0]) > 20:
                  
                    # get the data for this bin    
                    binT = this_month_T[bin_locs]                       
                    binD = this_month_D[bin_locs]

                    # find good temperatures
                    good_binT = np.where(binT.mask == False)
                    good_binD = np.where(binD.mask == False)

                    # and the bad dewpoints coincident with the good temperatures
                    bad_binD_T = np.where(binD.mask[good_binT] == False)

                    if len(bad_binD_T[0]) != 0:

                        bad_D_fraction = 1.-len(bad_binD_T[0])/float(len(good_binT[0]))

                        # if more than 50% missing
                        if bad_D_fraction >= 0.5:

                            # check the temporal resolution - if OK, then flag
                            #  This is a better temporal resolution calculation than IDL - will pick up more months
                            T_resoln = np.median(np.diff(goodT))
                            D_resoln = np.median(np.diff(goodD))

                            # only flag if resolutions the same and number of observations in total are similar
                            if (T_resoln != D_resoln) and (float(len(goodD[0]))/len(goodT[0]) < 0.666):
                                continue
                            else:
                                
                                flags[month[0]:month[1]][goodD] = 1
                               
                                ''' break this loop testing bins as whole month flagged
                                no point doing any more testing
                                return to month loop and move to next month'''
                                dt_month = start + dt.timedelta(hours = month[0])
                                if diagnostics:
                                    print dt.datetime.strftime(dt_month, "%B %Y"), bin, len(goodD[0])

                                if plots:
                                    # show the histogram
                                    import matplotlib.pyplot as plt
                                    
                                    plt.clf()
                                    plt.hist(this_month_T, bins=bins, color='r', label='Temperature', alpha=0.5)
                                    plt.hist(this_month_D, bins=bins, color='b', label='Dewpoint', alpha=0.5)
                                    
                                    plt.title(dt.datetime.strftime(dt_month, "%B %Y"))
                                    plt.legend()
                                    plt.show()

                                break


    # checked on 032220-99999 on 20/8/2014 and matches mostly
    # one extra month found - issues with IDL month start/end times.

    nflags = len(np.where(flags != 0)[0])
    utils.print_flagged_obs_number(logfile, "Dewpoint Cut-off", "temperature", nflags, noWrite = diagnostics)
        
    return flags # hcc_cutoffs
                           
#************************************************************************
def hcc(station, flag_col, start, end, logfile, diagnostics = False, plots = False):
    '''
    Run humidity cross checks on temperature and dewpoint temperature obs

    :param object station: station object
    :param array flag_col: which columns to fill in flag_array
    :param datetime start: start of dataset
    :param datetime end: end of dataset
    :param file logfile: logfile to store outputs
    :param bool plots: do plots or not
    :param bool diagnostics: extra verbose information

    :returns:    
    '''

    
    temperatures = getattr(station, 'temperatures')
    dewpoints = getattr(station, 'dewpoints')
    
    month_ranges = utils.month_starts_in_pairs(start, end)

    # Supersaturation
    station.qc_flags[:,flag_col[0]] = hcc_sss(temperatures.data, dewpoints.data, month_ranges, start, logfile, plots = plots, diagnostics = diagnostics) 
            
    # Dew point depression  
    precip1 = getattr(station, 'precip1_depth') 
    precip2 = getattr(station, 'precip2_depth') 
    precip3 = getattr(station, 'precip3_depth') 
    precip6 = getattr(station, 'precip6_depth') 
    precip9 = getattr(station, 'precip9_depth') 
    precip12 = getattr(station, 'precip12_depth') 
    precip15 = getattr(station, 'precip15_depth') 
    precip18 = getattr(station, 'precip18_depth') 
    precip24 = getattr(station, 'precip24_depth') 

    cloudbase = getattr(station, 'cloud_base')
    past_sigwx = getattr(station, 'past_sigwx1')
    times = station.time.data

    # combine all the precips together
    precips = np.array([precip1.data, precip2.data, precip3.data, precip6.data, precip9.data, precip12.data, precip15.data, precip18.data, precip24.data])
    
    station.qc_flags[:,flag_col[1]] = hcc_dpd(times, temperatures.data, dewpoints.data, precips, cloudbase.data, past_sigwx.data, start, logfile, plots = plots, diagnostics = diagnostics)    
   
    # Dew point cutoffs
    station.qc_flags[:,flag_col[2]] = hcc_cutoffs(temperatures.data, dewpoints.data, month_ranges, logfile, start, plots = plots, diagnostics = diagnostics) 
    
    for col in range(3):
        flag_locs = np.where(station.qc_flags[:, flag_col[col]] != 0)
        station.dewpoints.flags[flag_locs] = 1
    
    station = utils.append_history(station, "Temperature-Humidity Cross Check")  

    return # hcc
    
#************************************************************************
if __name__ == "__main__":
    
    print "temperature - humidity cross checks" 
    
    
    
    
