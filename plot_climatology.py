#!/usr/local/sci/bin/python
#***************************************
#
#   Copy of IDL code to produce climatology plots
#   25 September 2012 RJHD
#************************************************************************
#                    SVN Info
# $Rev:: 77                                            $:  Revision of last commit
# $Author:: rdunn                                      $:  Author of last commit
# $Date:: 2016-11-21 17:51:26 +0000 (Mon, 21 Nov 2016) $:  Date of last commit
#************************************************************************
#                                 START
#************************************************************************


'''
Takes one station, and plots the climatology and a specific year
'''

import numpy as np
import datetime as dt
import netCDF4 as ncdf
import matplotlib.pyplot as plt
import calendar
import scipy.stats
import scipy.signal
from PIL import Image
import os

import utils

DAYSPERYEAR=365.25
HOURSPERDAY=24
HADISD_VERSION = "v200_2015p"
DATALOCATION='/project/hadobs2/hadisd/{}/netcdf_files_{}/'.format(HADISD_VERSION, HADISD_VERSION)
DATASUFFIX='_mask.nc'

DATASTART=dt.datetime(1931,1,1,0,0)
DATAEND=dt.datetime(2016,1,1,0,0)

NORMSTART=dt.datetime(1975,1,1,0,0)
NORMEND=dt.datetime(2005,12,31,0,0)

OBSPERDAY=4
OBSDAYSPAN=12

month_starts = [0,31,59,90,120,151,181,212,243,273,304,334]

#***************************************
# Definitions/Subroutines
#***************************************
#***************************************
def DataRead(filename,variable, mask=False):    
    '''
    Read in data from the netcdf file given
      Extract the time and the variable specified

    :param string filename: full path to file
    :param string variable: variable as defined in netcdf file
    :param boolean mask: run mask
    :returns: np array of times, np array of data, float of missing value
    '''
    try:
        ncfile=ncdf.Dataset(filename,'r')

        timevar=ncfile.variables['time']

        var=ncfile.variables[variable] # this is a masked array

        try:
            dummy = var[:].mask
            masked_data=var[:].data # extract data with mdis where mask is True

        except AttributeError:
            masked_data=var[:] # no mask, so this is the data

        time_data=timevar[:] # this makes final filtering faster

        mdi=var.missing_value # find the fill value of the array

        if mask:
            # this is the slow step
            good_data=np.where(masked_data > mdi/2.) # filter masked and missing data
            data=masked_data[good_data]
            times=time_data[good_data]

        else:
            data=masked_data # convert masked array to normal with mdis
            times=time_data
    

        ncfile.close()
    except IOError,RuntimeError:
        print "File not found: ",filename
        raise IOError,RuntimeError

    return times,data,mdi # DataRead

#***************************************
def CalculateFullTimes(times,data,start,end,mdi):
    '''
    Match the compressed time axis onto a complete time axis

    :param array times: times
    :param array data: data
    :returns: fulltimes,fulldata
    '''

    dataspan=end-start
    fulltimes=np.arange(dataspan.days*HOURSPERDAY)
    inttimes=times.astype(int)
    mask=np.in1d(fulltimes,inttimes,assume_unique=True)

    # use mask to expand times and data
    fulltimes[mask]=inttimes
    fulldata=np.zeros(dataspan.days*HOURSPERDAY)
    fulldata.fill(mdi)
    fulldata[mask]=data

    return fulltimes,fulldata # CalculateFullTimes


#***************************************
def CalculateDailyMeans(times,data,start,end,mdi,OBSPERDAY = 4, OBSDAYSPAN = 12):
    '''
    Calculate daily mean value - 4 obs per day spread over at least 12 hours
    
    :param array times: np.array of timestamps
    :param array data: np.array of data
    :param float mdi: missing data indicator
    '''

    dataspan=end-start

    fulltimes,fulldata=CalculateFullTimes(times,data,start,end,mdi)

    # Ndays x 24hrs 
    fulldata=fulldata.reshape(dataspan.days,HOURSPERDAY)

    dailymeans=np.zeros(dataspan.days)
    dailymeans.fill(mdi)

    for day,thisdaydata in enumerate(fulldata):

        good_data=np.where(thisdaydata > mdi/2)

        if len(good_data[0])>=OBSPERDAY:
            # at least 4 observations in a day

            if good_data[0][-1]-good_data[0][0] >=OBSDAYSPAN:
                # at least 12 hours out of the day

                dailymeans[day]=np.mean(thisdaydata[good_data])

    return dailymeans, fulldata.reshape(dataspan.days*HOURSPERDAY) # CalculateDailyMeans

#***************************************
def binomial(n,k): 
    '''
    calculate binomial function
    '''
    bc = [1 for i in range(0,k+1)] 
    for j in range(1,n-k+1): 
        for i in range(1,k+1): 
            bc[i] = bc[i-1]+bc[i] 
    return bc[k] #binomial


#***************************************
def binomialfilter(data,mdi,n):
    '''
    Create binomial filter
    data - data array
    mdi - missing data indicator
    n - filter size
    '''
    

    tosmooth=np.concatenate((data,data,data))
    smoothed=np.zeros(len(tosmooth))
    smoothed.fill(mdi)

    weights=[]
    for k in range(n):
        weights+=[binomial(n-1,k)]

    weights=np.array(weights)

    for o,obs in enumerate(tosmooth):
        if (o >= n/2) and (o <= len(tosmooth)-n/2):
            
            chunk=tosmooth[o-(n/2):o+(n/2)+1]
            good=np.where(chunk != mdi)
            
            if len(good[0]) >= (n/3):
                
                norm=sum(weights[good])
                weighted=sum(chunk[good]*weights[good]/norm)

                smoothed[o]=weighted

    return smoothed[len(data):2*len(data)] #binomialfilter

#***************************************
def GetPercentiles(data, limit):
    '''
    Return a tuple of the upper and lower percentile values

    Inputs:
      data - flattened numpy array of the data
      limit - percentile to caclulate

    Outputs:
      p,p - tuple of upper and lower percentile
    '''

    return scipy.stats.scoreatpercentile(data,per=limit),scipy.stats.scoreatpercentile(data,per=100-limit) #GetPercentiles

#***************************************
def CalculateClimatologies(data, start, end, normstart, normend, mdi, percentile=95, hours=False):
    '''
    Calculate the climatological means, and upper and lower percentiles

    data - data (dailies/hourlies) in Nyears x Nperiods array
    mdi - missing data indicator
    percentiles - value of upper percentile
    hours - hourly data as opposed to daily
    '''

    reference_period=data[normstart.year-start.year:normend.year-end.year,:]

    if hours:
        xlen=365*HOURSPERDAY
    else:
        xlen=365
    
    # two versions of the climatologies
    Cmedian=np.zeros(xlen)
    Cmean=np.zeros(xlen)
    Cupper=np.zeros(xlen)
    Clower=np.zeros(xlen)
    
    # make the climatology for each day, and apply to create anomalies
    for obs in range(xlen):
        this_obs_data=reference_period[:,obs]
        good_data=np.where(this_obs_data > mdi/2)
        if len(good_data[0])>10:
            Cmedian[obs]=np.median(this_obs_data[good_data])
            Cmean[obs]=np.mean(this_obs_data[good_data])
            Cupper[obs],Clower[obs]=GetPercentiles(this_obs_data[good_data],percentile)
        else:
            Cmedian[obs]=Cmean[obs]=Cupper[obs]=Clower[obs]=mdi

    return Cmean, Cupper, Clower #CalculateClimatologies

#***************************************
def gauss_kern(size):
    """ Returns a normalized 1D gauss kernel array for convolutions """
    size = int(size)
    x=np.arange(-size,size+1)
    g = np.exp(-(x**2/float(size)))
    return g / g.sum() #gauss_kern


#***************************************
def gauss_smooth(data, n) :
    """ smooths the timeseries by convolving with a gaussian kernel of typical
        size n. The optional keyword argument ny allows for a different
        size in the y direction.
    """
    # as yearly climatologies, can wrap
    tosmooth=np.concatenate((data,data,data))

    g = gauss_kern(n)
    improc = scipy.signal.convolve(tosmooth, g, mode='valid')
    return improc[len(data)-n:2*len(data)-n] #gauss_smooth

#***************************************
def plot_figure(Cmean,Cupper,Clower,data, selyear, kernelsize, variable, stn, hours=False, selmonth='Jan',fill=False):
    '''
    Plot the climatological range, and the selected year

    Cmean,Cupper,Clower - the climatological mean, upper and lower percentiles
    data - the data (dailies/hourlies) - to get the selected year
    selyear - the selected year
    kernelsize - for the smoothing
    variable - what are we plotting?
    stn - station ID
    hours - to plot hourly data
    selmonth - which months to plot for hourly data
    '''
    
    plt.clf()
    fig=plt.figure()
    
    x_range=365

    # process the data
    smoothed_mean=binomialfilter(Cmean,mdi,kernelsize)
    good=np.where(smoothed_mean != mdi)
    
    # plot the smoothed mean
    plt.plot(np.arange(x_range)[good], smoothed_mean[good], 'k-')
        
    # smoothed percentiles - and fill between
    smoothed_upper,smoothed_lower=utils.binomialfilter(Cupper,mdi,kernelsize), utils.binomialfilter(Clower,mdi,kernelsize)
    good=np.where((smoothed_upper != mdi) & (smoothed_lower!=mdi))
    plt.fill_between(np.arange(x_range)[good],smoothed_upper[good],smoothed_lower[good] , facecolor='yellow', alpha=0.5, edgecolor='red')

    # the year in question
    thisdata=data[selyear-DATASTART.year]
    good=np.where(thisdata > mdi/2)
        
    if len(good[0])>0:

        plt.plot(np.arange(x_range)[good], thisdata[good], 'g-', lw=1.)

        if fill:
            # highlight hot/cold bits
            good=np.where((smoothed_upper != mdi) & (smoothed_lower!=mdi) & (thisdata > mdi/2))

            plt.fill_between(np.arange(x_range)[good], thisdata[good], smoothed_upper[good], where=thisdata[good]>smoothed_upper[good],facecolor='red')
            plt.fill_between(np.arange(x_range)[good], thisdata[good], smoothed_lower[good], where=thisdata[good]<smoothed_lower[good],facecolor='cyan')
        
    else:
        print "no data to plot for selected year "+str(selyear)
        
    # prettify
    ymin,ymax=plt.gca().get_ylim()
    plt.ylabel(variable.capitalize())

    # sort text on figure
    plt.text(10, ymin+(ymax-ymin)*0.95, 'Climatology period is '+str(NORMSTART.year)+' - '+str(NORMEND.year-1))
    plt.text(10, ymin+(ymax-ymin)*0.90, 'Station is '+str(stn))
    plt.text(10, ymin+(ymax-ymin)*0.85, 'Selected year is '+str(selyear))
    
    plt.xticks(month_starts,['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'], ha='left')

    # add logo to image
    im = Image.open('hc_logo_sm.jpg')
    height = im.size[1]

    # We need a float array between 0-1, rather than
    # a uint8 array between 0-255
    im = np.array(im).astype(np.float) / 255
    plt.figimage(im, 0, fig.bbox.ymax, zorder=0)

    # add text to show what code created this and when
    watermarkstring="/".join(os.getcwd().split('/')[4:])+'/'+os.path.basename( __file__ )+"   "+dt.datetime.strftime(dt.datetime.now(), "%d-%b-%Y %H:%M")
    plt.figtext(0.01,0.01,watermarkstring,size=6)

    # save the figure
    plt.xlim(0,x_range)
    plt.savefig(stn+'_'+variable[0:4]+'.png')
    plt.savefig(stn+'_'+variable[0:4]+'.eps')

    return # plot_figure

#***************************************
def plot_figure_hours(Cmean,Cupper,Clower,data, selyear, kernelsize, variable, stn, selmonth='Jan',fill=False):
    '''
    Plot the climatological range, and the selected year

    Cmean,Cupper,Clower - the climatological mean, upper and lower percentiles
    data - the data (dailies/hourlies) - to get the selected year
    selyear - the selected year
    kernelsize - for the smoothing
    variable - what are we plotting?
    stn - station ID
    selmonth - which months to plot for hourly data
    '''
    
    plt.clf()
    fig=plt.figure()

    x_range=365*HOURSPERDAY

    # set the start and end month
    start_month=np.where(np.array(calendar.month_abbr[:]) == selmonth)[0][0]
    start_day=dt.datetime(2009,start_month,1).timetuple().tm_yday
    start_hour=start_day*24-24

    # Terminate at end of year
    if start_month >= 11:
        end_hour=365*24
    else:
        end_day=dt.datetime(2010,start_month+2,1).timetuple().tm_yday
        end_hour=end_day*24
        
    # plot smoothed mean
    smoothed_mean=binomialfilter(Cmean,mdi,kernelsize)
    good=np.where(smoothed_mean != mdi)
    plt.plot(np.arange(x_range)[good], smoothed_mean[good], 'k-')
    
    # plot smoothed percentiles
    smoothed_upper,smoothed_lower=binomialfilter(Cupper,mdi,kernelsize), binomialfilter(Clower,mdi,kernelsize)
    good=np.where((smoothed_upper != mdi) & (smoothed_lower!=mdi))
    plt.fill_between(np.arange(x_range)[good],smoothed_upper[good],smoothed_lower[good] , facecolor='yellow', alpha=0.5, edgecolor='red')

    # plot year in question
    thisdata=data[selyear-DATASTART.year]
    good=np.where(thisdata > mdi/2)
        
    if len(good[0])>0:

        plt.plot(np.arange(x_range)[good], thisdata[good], 'g-', lw=1.)

        if fill:

            good=np.where((smoothed_upper != mdi) & (smoothed_lower!=mdi) & (thisdata > mdi/2))

            plt.fill_between(np.arange(x_range)[good], thisdata[good], smoothed_upper[good], where=thisdata[good]>smoothed_upper[good],facecolor='red')
            plt.fill_between(np.arange(x_range)[good], thisdata[good], smoothed_lower[good], where=thisdata[good]<smoothed_lower[good],facecolor='cyan')

        
    else:
        print "no data to plot for selected year "+str(selyear)
        
    # prettify
    ymin,ymax=plt.gca().get_ylim()   
    plt.text(start_hour+10, ymin+(ymax-ymin)*0.15, 'Climatology period is '+str(NORMSTART.year)+' - '+str(NORMEND.year-1))
    plt.text(start_hour+10, ymin+(ymax-ymin)*0.10, 'Station is '+str(stn))
    plt.text(start_hour+10, ymin+(ymax-ymin)*0.05, 'Selected year is '+str(selyear))

    # tick the days and months 
    months = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
    plt.xticks(np.arange(0,365*24,24),[months[x/24/29] if x in [i*24 for i in month_starts] else "" for x in np.arange(0,365*24,24)])

    # add logo to image
    im = Image.open('hc_logo_sm.jpg')
    height = im.size[1]

    # We need a float array between 0-1, rather than
    # a uint8 array between 0-255
    im = np.array(im).astype(np.float) / 255
    plt.figimage(im, 0, fig.bbox.ymax, zorder=0)

    # add text to show what code created this and when
    watermarkstring="/".join(os.getcwd().split('/')[4:])+'/'+os.path.basename( __file__ )+"   "+dt.datetime.strftime(dt.datetime.now(), "%d-%b-%Y %H:%M")
    plt.figtext(0.01,0.01,watermarkstring,size=6)

    plt.xlim(start_hour,end_hour)
    plt.savefig(stn+'_'+variable[0:4]+'_hours.png')
    plt.savefig(stn+'_'+variable[0:4]+'_hours.eps')

    return # plot_figure_hours

#***************************************
def main(stn, selyear, selmonth):

    kernelsize=11 # has to be integer

    truetimes, data, mdi = DataRead(DATALOCATION+stn+DATASUFFIX,variable)

    dataspan = DATAEND-DATASTART

    dailymeans, fulldata = CalculateDailyMeans(truetimes, data, DATASTART, DATAEND, mdi)

    dailymeans_year = np.zeros([(DATAEND.year - DATASTART.year), 365])

    hours_year = np.zeros([(DATAEND.year - DATASTART.year), 365 * HOURSPERDAY])

    yearstart = 0
    hourstart = 0

    for year in range(DATASTART.year, DATAEND.year):
        if calendar.isleap(year):

            this_year = np.concatenate((dailymeans[yearstart: yearstart+59], dailymeans[yearstart+60: yearstart+366]))
            yearstart += 366

            this_hours = np.concatenate((fulldata[hourstart: hourstart+(59*HOURSPERDAY)], fulldata[hourstart+(60*HOURSPERDAY): hourstart+366*HOURSPERDAY]))
            hourstart += 366*HOURSPERDAY

        else:
            this_year = dailymeans[yearstart: yearstart+365]
            yearstart += 365

            this_hours = fulldata[hourstart: hourstart+(365*HOURSPERDAY)]
            hourstart += 365*HOURSPERDAY

        dailymeans_year[year-DATASTART.year] = this_year

        hours_year[year-DATASTART.year] = this_hours

    climatologies_mean, climatologies_upper, climatologies_lower = CalculateClimatologies(dailymeans_year, DATASTART, DATAEND, NORMSTART, NORMEND, mdi)

    climatologies_mean_h, climatologies_upper_h, climatologies_lower_h = CalculateClimatologies(hours_year, DATASTART, DATAEND, NORMSTART, NORMEND, mdi, hours=True)

    plot_figure(climatologies_mean, climatologies_upper, climatologies_lower, dailymeans_year, selyear, 11, variable, stn, fill=True)

    plot_figure_hours(climatologies_mean_h, climatologies_upper_h, climatologies_lower_h, hours_year, selyear, 11, variable, stn, fill=True, selmonth = selmonth)

    return # main

#************************************************************************
if __name__=="__main__":

    import argparse

    # set up keyword arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--stn', dest='stn', action='store', default = "",
                        help='Station ID, default = ""')
    parser.add_argument('--year', dest='year', action='store', default = "",
                        help='Year, default = ""')
    parser.add_argument('--month', dest='month', action='store', default = "",
                        help='Month (Jan, Feb... etc), default = ""')
    parser.add_argument('--variable', dest='variable', action='store', default = "temperatures",
                        help='Variable, default = "temperatures"')
 
    args = parser.parse_args()

    main(stn, year, month, variable)


#************************************************************************
#                                 END
#************************************************************************

