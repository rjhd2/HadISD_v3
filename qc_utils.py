#!/usr/local/sci/bin/python
#*****************************
#
# general utilities & classes for Python QC.
#
#
#************************************************************************
#                    SVN Info
#$Rev:: 117                                           $:  Revision of last commit
#$Author:: rdunn                                      $:  Author of last commit
#$Date:: 2017-01-30 15:33:46 +0000 (Mon, 30 Jan 2017) $:  Date of last commit
#************************************************************************

import numpy as np
import datetime as dt
from scipy.optimize import leastsq,fsolve

#*********************************************
class MetVar(object):
    '''
    Class for meteorological variable
    '''
    
    def __init__(self, name, long_name):
        self.name = name
        self.long_name = long_name
        

    def __str__(self):     
        return "variable: {}, long_name: {}".format(self.name, self.long_name)

    __repr__ = __str__
   
   
   
#*********************************************
class Station(object):
    '''
    Class for station
    '''
    
    def __init__(self, stn_id, lat, lon, elev):
        self.id = stn_id
        self.lat = lat
        self.lon = lon
        self.elev = elev
        
    def __str__(self):
        return "station {}, lat {}, lon {}, elevation {}".format(self.id, self.lat, self.lon, self.elev)
    
    __repr__ = __str__
    
    
 
#*********************************************
def set_MetVar_attributes(name, long_name, standard_name, units, mdi, dtype):
    '''
    Wrapper to set up a new MetVar object and populate some of the attibute fields

    :param str name: name for variable - ideally CF compliant
    :param str long_name: longname for variable - ideally CF compliant
    :param str standard_name: standard_name for variable - ideally CF compliant
    :param str units: units for variable - ideally CF compliant
    :param float/int mdi: missing data indicator
    :param type dtype: dtype for variable

    :returns MetVar new_var: new MetVar variable.
    '''

    new_var = MetVar(name, long_name)
    new_var.units = units
    new_var.dtype = dtype
    new_var.mdi = mdi
    new_var.standard_name = standard_name
    
    return new_var # set_MetVar_attributes
   
#*********************************************
def create_fulltimes(station, var_list, start, end, opt_var_list = [], do_input_station_id = True, do_qc_flags = True, do_flagged_obs = True):
    '''
    expand the time axis of the variables.
    ''' 
    
    # print "Expanding time axis"

    time_range = end - start
    fulltimes = np.arange(time_range.days * 24)

    # adjust if input netCDF file has different start date to desired
    netcdf_start = dt.datetime.strptime(station.time.units.split()[2], '%Y-%m-%d')

    offset = start - netcdf_start
    offset = offset.days *24
    fulltimes = fulltimes + offset

    match = np.in1d(fulltimes, station.time.data)
    match_reverse = np.in1d(station.time.data, fulltimes)

    # if optional/carry through variables given, then set to extract these too
    if opt_var_list != []:
        full_var_list = np.append(var_list, opt_var_list)
    else:
        full_var_list = var_list
  
    if do_input_station_id:
        final_var_list = np.append(full_var_list, ["input_station_id"])
    else:
        final_var_list = full_var_list

    
    for variable in final_var_list:
        st_var = getattr(station, variable)
               
        # use masked arrays for ease of filtering later
        if variable in ["input_station_id"]:
            new = np.ma.array(["            " for i in range(len(fulltimes))], fill_value = st_var.mdi, dtype=(str))
        elif variable in ["precip1_condition","windtypes","precip2_condition","precip3_condition","precip4_condition"]:
            new = np.ma.array([" " for i in range(len(fulltimes))], fill_value = st_var.mdi, dtype=(str))
        else:
            new = np.ma.zeros(len(fulltimes), fill_value = st_var.mdi)

        new.fill(st_var.mdi)
        new.mask = True
        
        new[match] = st_var.data[match_reverse]
        new.mask[match] = False # unmask the filled timestamps
        
        # but re-mask those filled timestamps which have missing data
        st_var.data = np.ma.masked_where(new == st_var.mdi, new)
        
        if variable in var_list and do_flagged_obs == True:
            # flagged values
            new = np.zeros(len(fulltimes))
            new.fill(st_var.mdi)
            new[match] = st_var.flagged_obs[match_reverse]
            st_var.flagged_obs = new
        
            # flags - for filtering
            new = np.zeros(len(fulltimes))
            new.fill(st_var.mdi)
            new[match] = st_var.flags[match_reverse]
            st_var.flags = new

   

    # do the QC flags, using try/except
    if do_qc_flags == True:
        try:
            qc_var = getattr(station, "qc_flags")

            # use masked arrays for ease of filtering later
            new = np.zeros([len(fulltimes), qc_var.shape[1]])

            new[match, :] = qc_var[match_reverse, :]

            station.qc_flags = new

        except AttributeError:
            pass

    # working in fulltimes throughout and filter by missing
    if offset != 0:
        fulltimes = np.arange(time_range.days * 24)

    station.time.data = fulltimes
    
    return match

#*********************************************
def month_starts(start, end):
    '''
    Returns locations of month starts (using hours as index)
    '''
    
    month_locs = []
    
    date = start
    
    while date < end:
        
        difference = date - start
        
        month_locs += [difference.days*24]
    
        # increment counter
        if date.month < 12:
            date = dt.datetime(date.year, date.month+1, 1)
        else:
            date = dt.datetime(date.year+1, 1, 1)
    
    return month_locs # month_starts

#*********************************************
def month_starts_in_pairs(start, end):
    '''
    Create array of month start/end pairs
    :param datetime start: start of data 
    :param datetime end: end of data
    :returns: month_ranges: Nx2 array    
    '''
    # set up the arrays of month start locations
    m_starts = month_starts(start, end)
    
    month_ranges = np.zeros((len(m_starts),2))
    month_ranges[:-1,0] = m_starts[:-1]
    month_ranges[:-1,1] = m_starts[1:]

    difference = end - start

    month_ranges[-1,:] = [m_starts[-1], difference.days * 24.]
    
    return month_ranges # month_starts_in_pairs

#*********************************************
def reporting_accuracy(indata, winddir = False, plots = False):
    '''
    Following reporting_accuracy.pro method.
    Uses histogram of remainders to look for special values

    :param array indata: masked array
    :param bool winddir: true if processing wind directions
    :param bool plots: make plots (winddir only)

    :output: resolution - reporting accuracy (resolution) of data
    '''
    
    
    good_values = indata.compressed()

    resolution = -1
    if winddir:
        # 360/36/16/8/ compass points ==> 1/10/22.5/45/90 deg resolution
        if len(good_values) > 0:

            hist, binEdges = np.histogram(good_values, bins = np.arange(0,362,1))

            # normalise
            hist = hist / float(sum(hist))

            #
            if sum(hist[np.arange(90,360+90,90)]) >= 0.6:
                resolution = 90
            elif sum(hist[np.arange(45,360+45,45)]) >= 0.6:
                resolution = 45
            elif sum(hist[np.round(0.1 + np.arange(22.5,360+22.5,22.5)).astype("int")]) >= 0.6:
                # added 0.1 because of floating point errors!
                resolution = 22
            elif sum(hist[np.arange(10,360+10,10)]) >= 0.6:
                resolution = 10
            else:
                resolution = 1

            print "Wind dir resolution = {} degrees".format(resolution)
            if plots:
                import matplotlib.pyplot as plt
                plt.clf()
                plt.hist(good_values, bins = np.arange(0,362,1))
                plt.show()
        
    else:
        if len(good_values) > 0:

            remainders = np.abs(good_values) - np.floor(np.abs(good_values))

            hist, binEdges = np.histogram(remainders, bins = np.arange(-0.05,1.05,0.1))

            # normalise
            hist = hist / float(sum(hist))

            if hist[0] >= 0.3:
                if hist[5] >= 0.15:
                    resolution = 0.5
                else:
                    resolution = 1.0
            else:
                resolution = 0.1

    return resolution # reporting_accuracy

#*********************************************
def reporting_frequency(indata):
    '''
    Following reporting_accuracy.pro method.
    Uses histogram of remainders to look for special values

    :param array indata: masked array
    :output: frequency - reporting frequency of data
    '''
        
    masked_locs, = np.where(indata.mask == False)    

    frequency = -1
    if len(masked_locs) > 0:

        difference_series = np.diff(masked_locs)

        hist, binEdges = np.histogram(difference_series, bins = np.arange(1,25,1), density=True)
        # 1,2,3,6
        if hist[0] >= 0.5:
            frequency = 1
        elif hist[1] >= 0.5:
            frequency = 2
        elif hist[2] >= 0.5:
            frequency = 3
        elif hist[3] >= 0.5:
            frequency = 4
        elif hist[5] >= 0.5:
            frequency = 6
        else:
            frequency = 24
       
    return frequency # reporting_frequency

#*********************************************
def gaussian(X,p):
    '''
    Gaussian function for line fitting
    p[0]=mean
    p[1]=sigma
    p[2]=normalisation
    '''
    return (p[2]*(np.exp(-((X-p[0])*(X-p[0]))/(2.0*p[1]*p[1])))) # gaussian

#*********************************************
def invert_gaussian(Y,p):
    '''
    X value of Gaussian at given Y
    p[0]=mean
    p[1]=sigma
    p[2]=normalisation
    '''
    return p[0] + (p[1]*np.sqrt(-2*np.log(Y/p[2]))) # invert_gaussian

#*********************************************
def residuals_gaussian(p, Y, X):
    '''
    Least squared residuals from linear trend
    '''
    err = ((Y-gaussian(X,p))**2.0)

    return err # residuals_gaussian

#*********************************************
def fit_gaussian(x,y,norm, mu=False, sig=False):
    '''
    Fit a straight line to the data provided
    Inputs:
      x - x-data
      y - y-data
      norm - norm
    Outputs:
      fit - array of [mu,sigma,norm]
    '''
    if not mu:
        mu=np.mean(x)
    if not sig:
        sig=np.std(x)
    p0 = np.array([mu,sig,norm])
    fit,success=leastsq(residuals_gaussian, p0, args=(y,x), maxfev=10000,full_output=False)

    return fit # fit_gaussian

#*********************************************
def apply_filter_flags(st_var):
    '''
    Return the data masked by the flags
    '''

    return  np.ma.masked_where(st_var.flags == 1, st_var.data) # apply_filter_flags

#*********************************************
def apply_flags_to_mask(station, variable):
    '''
    Return the data as masked array including the flagged values
    '''
    st_var = getattr(station, variable)

    st_var.data.fill_value = st_var.mdi

    st_var.data = np.ma.masked_where(st_var.data == st_var.fdi, st_var.data)

    return  st_var # apply_flags_to_mask

#*********************************************
def IQR(data, percentile = 0.25):
    ''' Calculate the IQR of the data '''
    # perhaps combine with percentile - but this may be more efficient
    sorted_data = sorted(data)
    
    n_data = len(sorted_data)

    quartile = int(round(percentile * n_data))
                       
    return sorted_data[n_data - quartile] - sorted_data[quartile] # IQR
    
#*********************************************
def mean_absolute_deviation(data, median = False):    
    ''' Calculate the MAD of the data '''
    
    if median:
        mad = np.ma.mean(np.ma.abs(data - np.ma.median(data)))
        
    else:        
        mad = np.ma.mean(np.ma.abs(data - np.ma.mean(data)))

    return mad # mean_absolute_deviation

#*********************************************
def percentiles(data, percent, idl = False):
    ''' Calculate the percentile of data '''
    
    sorted_data = sorted(data)
    
    if idl:
        n_data = len(data)-1
        percentile = sorted_data[int(np.ceil(n_data * percent))] # matches IDL formulation
    else:
        n_data = len(data)
        percentile = sorted_data[int(n_data * percent)]
        
    
    return percentile # percentile

#*********************************************
def winsorize(data, percent, idl = False):
    

    for pct in [percent, 1-percent]:
        
        
        if pct < 0.5:
            percentile = percentiles(data, pct, idl = idl)
            locs = np.where(data < percentile)
        else:
            percentile = percentiles(data, pct, idl = idl)
            locs = np.where(data > percentile)
        
        data[locs] = percentile

    return data # winsorize

#*********************************************
def times_hours_to_datetime(times, start):
    ''' convert the hours since into datetime objects
    makes for more intelligible plots'''

    import matplotlib.dates as mdt

    offset = mdt.date2num(start)

    return np.array(mdt.num2date(offset + (times/24.)))

#************************************************************************
def create_bins(indata, binwidth):
    ''' create bins and bin centres from data 
    given bin width covers entire range'''

    # set up the bins
    bmins = np.floor(np.min(indata))
    bmaxs = np.ceil(np.max(indata))
    bins = np.arange(bmins - binwidth, bmaxs + (3. * binwidth), binwidth)
    bincenters = 0.5 * (bins[1:] + bins[:-1])

    return bins, bincenters # create_bins

#************************************************************************
def print_flagged_obs_number(logfile, test, variable, nflags, noWrite=False):

    if noWrite:
        print "{:50s} {:20s}{:5d}\n".format(test+" Check Flags :", variable.capitalize()+" :", nflags)
    else:
        logfile.write("{:50s} {:20s}{:5d}\n".format(test+" Check Flags :", variable.capitalize()+" :", nflags))
    return # print_flagged_obs_number

#************************************************************************
def sort_ts_ylim(data):
    ''' sort the y-limit of timeseries plots '''

    dmax = np.max(data)
    dmin = np.min(data)
    
    if dmax > 0:
        dmax = dmax * 2
    else:
        dmax = dmax / 2

    if dmin > 0:
        dmin = dmin / 2
    else:
        dmin = dmin * 2

    return [dmax, dmin] # print_flagged_obs_number

#************************************************************************
def idl_median(indata):
    ''' matches IDL version of median '''

    if len(indata)/2. ==  len(indata)/2:
        # even
        return sorted(indata)[len(indata)/2]

    else:
        return np.median(indata) # idl_median

#************************************************************************
def get_dist_and_bearing(coord1,coord2):
    '''
    Get the distance between two points long Earth's surface
    '''
    def get_phi(lat):
        return np.deg2rad(90. - lat)
    

    lat1, lon1 = coord1
    lat2, lon2 = coord2
    R = 6371.229 # km

    phi1 = get_phi(lat1)
    phi2 = get_phi(lat2)

    theta1 = np.deg2rad(lon1)
    theta2 = np.deg2rad(lon2)

    cos = (np.sin(phi1) * np.sin(phi2) * np.cos(theta1 - theta2) + np.cos(phi1) * np.cos(phi2))
    arc = np.arccos( cos )

    lat1, lon1 = np.deg2rad((lat1, lon1))
    lat2, lon2 = np.deg2rad((lat2, lon2))

    bearing = np.rad2deg(np.arctan2(np.sin(lon2-lon1)*np.cos(lat2), np.cos(lat1)*np.sin(lat2)-np.sin(lat1)*np.cos(lat2)*np.cos(lon2-lon1)))

    try:
        if bearing < 0: bearing += 360.
    except ValueError:
        # then a list or array
        bearing[bearing < 0] = bearing[bearing < 0] + 360.

    distance = arc * R

    return distance.astype("int"), bearing.astype("int") # get_dist_and_bearing

#************************************************************************
def concatenate_months(month_ranges, data, hours = True):
    '''
    Sum up a single month across all years (e.g. all Januaries)
    '''
    

    datacount = np.zeros(month_ranges.shape[0])
    year_ids = [] 

    for y, year in enumerate(month_ranges):
        
        this_year = data[month_ranges[y][0]:month_ranges[y][1]]

        datacount[y] = len(this_year.compressed())

        if y == 0:
            
            # store so can access each hour of day separately
            if hours:
                this_month = this_year.reshape(-1,24)           
            else:
                this_month = this_year[:] # to ensure a copy not a view
            year_ids = [y for x in range(this_month.shape[0])]
            
        else:
            if hours:
                this_year = this_year.reshape(-1,24)            
                this_month = np.ma.concatenate((this_month, this_year), axis = 0)            
            else:
                this_month = np.ma.concatenate((this_month, this_year))
            
            year_ids.extend([y for x in range(this_year.shape[0])])
            

    return this_month, year_ids, datacount # concatenate_months

#************************************************************************
def mask_old(station, var_list):
    '''
    Apply the flags to the data and copy across to storage attribute

    :param object station: station object
    :param list var_list: list of variables to process
    :returns:
        station - updated station
    '''


    for variable in var_list:
        st_var = getattr(station, variable)

        flags = np.ma.where(st_var.flags != 0)

        st_var.flagged_obs[flags] = st_var.data[flags]

        st_var.data[flags] = st_var.fdi

    station = append_history(station, "Masking")

    return station # mask_old

#************************************************************************
def mask(station, var_list, logfile, FLAG_COL_DICT):
    '''
    Apply the flags to the data and copy across to storage attribute
    Uses flag array rather than built in system

    :param object station: station object
    :param list var_list: list of variables to process
    :param file logfile: log file to write to
    :param dict FLAG_COL_DICT: dictionary of flag columns to apply

    :returns:
        station - updated station
    '''

    for variable in var_list:
        st_var = getattr(station, variable)

        # winds logical test notes recovered directions with -1 flag
        temp_flags = station.qc_flags[:, FLAG_COL_DICT[variable]]
        neg_locs = np.where(temp_flags < 0)
        temp_flags[neg_locs] = 0

        flags = np.sum(temp_flags, axis = 1)

        flag_locs = np.ma.where(flags != 0)

        st_var.flagged_obs[flag_locs] = st_var.data[flag_locs]

        st_var.data[flag_locs] = st_var.fdi

        if logfile == "":
            print "Mask applied to {}".format(variable)
        else:
            logfile.write("Mask applied to {}\n".format(variable))

    station = append_history(station, "Masking")

    return station # mask

#************************************************************************
def append_history(station, text):
    '''
    Append text to the station history attribute

    :param object station: station object
    :param str text: text to append with date.

    '''

    station.history = station.history + text + dt.datetime.strftime(dt.datetime.now(), " %Y-%m-%d, %H:%M \n")
    
    return station # append_history

#************************************************************************
def monthly_reporting_statistics(st_var, start, end):
    '''
    Return reporting accuracy & reporting frequency for variable

    :param obj st_var: station variable object
    :param datetime start: start of data series
    :param datatime end: end of data series

    :returns:
       reporting_stats - Nx2 array, one pair for each month
    '''



    monthly_ranges = month_starts_in_pairs(start, end)

    reporting_stats = -np.ones(monthly_ranges.shape)

    for m, month in enumerate(monthly_ranges):

        reporting_stats[m] = [reporting_frequency(st_var.data[month[0]:month[1]]),reporting_accuracy(st_var.data[month[0]:month[1]])]

    return reporting_stats # monthly_reporting_statistics

#***************************************
def gausshermiteh3h4(x, A, x0, s, h3, h4):
    '''
    The Gauss-Hermite function is a superposition of functions of the form
    F = (x-xc)/s                                            
    E =  A.Exp[-1/2.F^2] * {1 + h3[c1.F+c3.F^3] + h4[c5+c2.F^2+c4.F^4]} 
    From http://www.astro.rug.nl/software/kapteyn-beta/plot_directive/EXAMPLES/kmpfit_gausshermite.py
    '''
    c0 =     np.sqrt(6.0)/4.0
    c1 =    -np.sqrt(3.0)
    c2 =    -np.sqrt(6.0)
    c3 = 2.0*np.sqrt(3.0)/3.0
    c4 =     np.sqrt(6.0)/3.0
    
    F = (x-x0)/s
    E = A*np.exp(-0.5*F*F)*( 1.0 + h3*F*(c3*F*F+c1) + h4*(c0+F*F*(c2+c4*F*F)) )

    return E # gausshermiteh3h4


#***************************************
def hermite2gauss(par, diagnostics = False):
    '''
    Convert Gauss-Hermite parameters to Gauss(like)parameters.
                                                                
    We use the first derivative of the Gauss-Hermite function   
    to find the maximum, usually around 'x0' which is the center
    of the (pure) Gaussian part of the function.                          
    If F = (x-x0)/s then the function for which we want the     
    the zero's is A0+A1*F+A2*F^2+A3*F^3+A4*F^4+A5*F^5 = 0       
    c0 = 1/4sqrt(6) c1 = -sqrt(3) c2 = -sqrt(6)                 
    c3 = 2/3sqrt(3) c4 = 1/3sqrt(6)                             
    From http://www.astro.rug.nl/software/kapteyn-beta/plot_directive/EXAMPLES/kmpfit_gausshermite.py
    removed the error calculation as not using the Kapteyn "fitter" object/function
    '''
    sqrt2pi = np.sqrt(2.0*np.pi)
    amp, x0, s, h3, h4 = par
    c0 = np.sqrt(6.0)/4.0
    c1 = -np.sqrt(3.0)
    c2 = -np.sqrt(6.0)
    c3 = 2.0*np.sqrt(3.0)/3.0
    c4 = np.sqrt(6.0)/3.0

    A = np.zeros(6)
    A[0] = -c1*h3
    A[1] = h4*(c0-2.0*c2) + 1.0
    A[2] = h3*(c1-3.0*c3)
    A[3] = h4*(c2 - 4.0*c4)
    A[4] = c3*h3
    A[5] = c4*h4
    
    # Define the function that represents the derivative of
    # the GH function. You need it to find the position of the maximum.
    fx = lambda x: A[0] + x*(A[1]+x*(A[2]+x*(A[3]+x*(A[4]+x*A[5]))))
    xr = fsolve(fx, 0, full_output=True)
    xm = s*xr[0] + x0
    ampmax = gausshermiteh3h4(xm, amp, x0, s, h3, h4)
    
    # Get line strength
    f = 1.0 + h4 * np.sqrt(6.0) / 4.0
    area  = amp * s * f * sqrt2pi

    # Get mean
    mean  = x0 + np.sqrt(3.0)*h3*s
    
    # Get dispersion
    f = 1.0 + h4*np.sqrt(6.0)
    dispersion = abs(s * f)

    # Skewness
    f = 4.0 * np.sqrt(3.0)
    skewness = f * h3
    
    # Kurtosis
    f = 8.0 * np.sqrt(6.0)
    kurtosis = f * h4
    
    res = dict(xmax=xm, amplitude=ampmax, area=area, mean=mean, dispersion=dispersion,\
                   skewness=skewness, kurtosis=kurtosis)

    if diagnostics:
        print "Gauss-Hermite max=%g at x=%g"%(res['amplitude'], res['xmax'])
        print "Area      :", res['area']
        print "Mean (X0) :", res['mean']
        print "Dispersion:", res['dispersion']
        print "Skewness  :", res['skewness']
        print "Kurtosis  :", res['kurtosis']

    return res # hermite2gauss

#***************************************
def funcGH(p, x):
    # Model is a Gauss-Hermite function
    A, xo, s, h3, h4 = p
    return gausshermiteh3h4(x, A, xo, s, h3, h4) # funcGH

#***************************************
def residualsGH(p, data):
    # Return weighted residuals of Gauss-Hermite
    x, y, err = data
    return (y-funcGH(p,x)) / err # residualsGH

#*********************************************
def linear(X,p):
    '''
    decay function for line fitting
    p[0]=intercept
    p[1]=slope
    '''
    return p[1]*X + p[0] # linear

#*********************************************
def residuals_linear(p, Y, X):
    '''
    Least squared residuals from linear trend
    '''
    err = ((Y-linear(X,p))**2.0)

    return err # residuals_linear


#*********************************************
def plot_log_distribution(edges, hist, fit, threshold, line_label, xlabel, title, old_threshold = 0):

    import matplotlib.pyplot as plt
    
    plt.clf()
    # stretch bars, so can run off below 0
    plot_hist = np.array([np.log10(x) if x != 0 else -1 for x in hist])
    plt.step(edges[1:], plot_hist, color = 'k', label = line_label)
    plt.plot(edges, fit, 'b-', label = "best fit")          
    
    plt.xlabel(xlabel)
    plt.ylabel("log10(Frequency)")
    
    # set y-lim to something sensible
    plt.ylim([-0.3, max(plot_hist)])
    plt.xlim([0, max(edges)])
    
    plt.axvline(threshold, c = 'r', label = "threshold = {}".format(threshold))
    if old_threshold != 0:
        plt.axvline(old_threshold, c = 'g', label = "old threshold = {}".format(old_threshold))
    
    plt.legend(loc = "upper right")
    plt.title(title)
       
    plt.show()

    return # plot_log_distribution

#*********************************************
def get_critical_values(indata, binmin = 0, binwidth = 1, plots = False, diagnostics = False, line_label = "", xlabel = "", title = "", old_threshold = 0):
    """
    Plot histogram on log-y scale and fit 1/x decay curve to set threshold

    :param array indata: input data to bin up
    :param int binmin: minimum bin value
    :param int binwidth: bin width
    :param bool plots: do the plots
    :param bool diagnostics : do diagnostic outputs
    :param str line_label: label for plotted histogram
    :param str xlabel: label for x axis
    :param str title: plot title
    :param float old_threshold: (spike) plot the old threshold from IQR as well

    :returns:
       threshold value

    """
    
    if len(set(indata)) > 1:

        bins = np.arange(binmin, 3 * max(np.ceil(np.abs(indata))), binwidth)
        full_hist, full_edges = np.histogram(np.abs(indata), bins = bins)

        if len(full_hist) > 1:

            # use only the central section (as long as it's not just 2 bins)
            i = 0
            limit = 0
            while limit < 2:
                try:
                    limit = np.argwhere(full_hist == 0)[i][0]
                    i += 1
                except IndexError:
                    # no zero bins in this histogram
                    limit = len(full_hist)
                    break          

            edges = full_edges[:limit]
            hist  = np.log10(full_hist[:limit])

            # Working in log-yscale from hereon

            # a 10^-bx
            a = hist[np.argmax(hist)]
            b = 1

            p0 = np.array([a,b])
            fit,success=leastsq(residuals_linear, p0, args=(hist, edges), maxfev=10000,full_output=False)

            fit_curve = linear(full_edges, fit)

            if fit[1] < 0:
                # in case the fit has a positive slope

                # where does fit fall below log10(-0.1)
                try:
                    fit_below_point1, = np.argwhere(fit_curve < -1)[0]
                    first_zero_bin, = np.argwhere(full_hist[fit_below_point1:] == 0)[0] + 1
                    threshold = binwidth * (fit_below_point1 + first_zero_bin)
                except IndexError:
                    # too shallow a decay - use default maximum
                    threshold = len(full_hist)

                # find first empty bin after that

            else:
                threshold = len(full_hist)


            if plots:
                plot_log_distribution(full_edges, full_hist, fit_curve, threshold, line_label, xlabel, title, old_threshold = old_threshold)
                
        else:
            threshold = max(indata) + binwidth

    else:
        threshold = max(indata) + binwidth
 
    return threshold # get_critical_values

#************************************************************************
def apply_flags_all_variables(station, all_variables, flag_col, logfile, test_name, plots = False, diagnostics = False):
    """
    Apply these flags to all variables

    :param object station: the station object to be processed
    :param list all_variables: the variables where the flags are to be applied
    :param list flag_col: which column in the qc_flags array to work on
    :param file logfile: logfile to store outputs
    :param str test_name: test name for printing/loggin
    :param bool plots: to do any plots
    :param bool diagnostics: do any extra diagnostic output
    :returns:
    """
    
    flag_locs, = np.where(station.qc_flags[:, flag_col] != 0)

    for var in all_variables:
        st_var = getattr(station, var)
        
        # copy flags into attribute
        st_var.flags[flag_locs] = 1
   
        if plots or diagnostics:
            print "Applying {} flags to {}".format(test_name, var)
        else:
            logfile.write("Applying {} flags to {}\n".format(test_name, var))

    return # apply_flags_all_variables

#************************************************************************
def apply_windspeed_flags_to_winddir(station, diagnostics = False):
    """
    Applying windspeed flags to wind directions synergistically
    Called after every test which assess windspeeds

    :param object station: the station object to be processed
    :param bool diagnostics: do any extra diagnostic output
    """


    windspeeds = getattr(station, "windspeeds")
    winddirs = getattr(station, "winddirs")
    
    winddirs.flags = windspeeds.flags

    if diagnostics:

        old_flags, = np.where(winddirs.flags != 0)
        new_flags, = np.where(windspeeds.flags != 0)

        print "{} flags copied from windspeeds to winddirs".format(len(new_flags) - len(old_flags))

    return # apply_windspeed_flags_to_winddir

#************************************************************************
def nearly_equal(a,b,sig_fig=5):
    """
    Returns it two numbers are nearly equal within sig_fig decimal places

http://stackoverflow.com/questions/558216/function-to-determine-if-two-numbers-are-nearly-equal-when-rounded-to-n-signific

    :param flt a: number 1
    :param flt b: number 2
    :param int sig_fig: number of decimal places to check agreement to

    :returns bool:    
    """

    return ( a==b or 
             int(a*10**sig_fig) == int(b*10**sig_fig)
           ) # nearly_equal
