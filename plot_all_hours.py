#!/usr/local/sci/bin/python
############################COPYRIGHT###############################
# British Crown Copyright (c) 2013, the Met Office All rights reserved.
#
# Redistribution and use in source and binary forms, with or without 
# modification, are permitted provided that the following conditions 
# are met:
#
# Redistributions of source code must retain the above copyright notice, 
# this list of conditions and the following disclaimer. 
#
# Redistributions in binary form must reproduce the above copyright 
# notice, this list of conditions and the following disclaimer in the 
# documentation and/or other materials provided with the distribution. 
#
# Neither the name of the Met Office nor the names of its contributors 
# may be used to endorse or promote products derived from this software 
# without specific prior written permission. 
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR 
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
############################COPYRIGHT###############################
#***************************************
#   Quickly reads file and extracts one variable to plot interactively
#  
#   23 Jan 2013 RJHD
#***************************************

import sys
import numpy as np
import netCDF4 as ncdf
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt 
import os

import datetime as dt
import matplotlib.patches as patches
import matplotlib.dates as mdates
import matplotlib.ticker as ticker
from matplotlib import rcParams

import utils
from set_paths_and_vars import *

if len(sys.argv) < 2:
    print "usage: plot_all_hours stationID variable"
    sys.exit()

stn=sys.argv[1]
if len(sys.argv) == 3:
    var=sys.argv[2]
else:
    var='temperatures'


rcParams["xtick.direction"] = "out"
rcParams["ytick.direction"] = "out"

HOURSPERDAY=24
DAYSPERYEAR=365


#***************************************
def ReadStations(filename):
    ''' Use numpy genfromtxt reading to read in all station data in ID,Lat,Lon,Elev list'''

    return np.genfromtxt(filename, dtype=(str))
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
        input_stn_var=ncfile.variables['input_station_id'] # this is a masked array

        masked_data=var[:].data # extract data with mdis where mask is True
        input_stn_id=input_stn_var[:]
        
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

    return times,data,mdi,input_stn_id # DataRead
#***************************************
def DataRead(filename, variable, mask=False):    
    '''
    Read in data from the netcdf file given
      Extract the time and the variable specified

    :param string filename: full path to file
    :param string variable: variable as defined in netcdf file
    :param boolean mask: run mask
    :returns: np array of times, np array of data, float of missing value
    '''
    try:
        ncfile = ncdf.Dataset(filename,'r')

        timevar = ncfile.variables['time']

        var = ncfile.variables[variable] # this is a masked array
        input_stn_var = ncfile.variables['input_station_id'] # this is a masked array

        masked_data = var[:].data # extract data with mdis where mask is True
        input_stn_id = input_stn_var[:]
        
        time_data = timevar[:] # this makes final filtering faster

        mdi = var.missing_value # find the fill value of the array

        if mask:
            # this is the slow step
            good_data = np.where(masked_data > mdi/2.) # filter masked and missing data
            data = masked_data[good_data]
            times = time_data[good_data]

        else:
            data = masked_data # convert masked array to normal with mdis
            times = time_data
    

        ncfile.close()
    except IOError,RuntimeError:
        print "File not found: ", filename
        raise IOError,RuntimeError

    return times, data, mdi, input_stn_id # DataRead

#***************************************
def CalculateFullTimes(times,data,mdi):
    '''
    Match the compressed time axis onto a complete time axis

    :param array times: times
    :param array data: data
    :returns: fulltimes,fulldata
    '''


    dataspan=DATAEND-DATASTART
    fulltimes=np.arange(dataspan.days*HOURSPERDAY)
    inttimes=times.astype(int)
    mask=np.in1d(fulltimes,inttimes,assume_unique=True)

    # use mask to expand times and data
    fulltimes[mask]=inttimes
    fulldata=np.zeros(dataspan.days*HOURSPERDAY)
    fulldata.fill(mdi)
    fulldata[mask]=data

    return fulltimes,fulldata,mask # CalculateFullTimes

#***************************************
def CalculateFullTimes(times, data, mdi):
    '''
    Match the compressed time axis onto a complete time axis

    :param array times: times
    :param array data: data
    :returns: fulltimes,fulldata
    '''

    dataspan = DATAEND - DATASTART
    fulltimes = np.arange(dataspan.days * HOURSPERDAY)
    inttimes = times.astype(int)
    mask = np.in1d(fulltimes, inttimes, assume_unique=True)

    # use mask to expand times and data
    fulltimes[mask] = inttimes
    fulldata = np.zeros(dataspan.days*HOURSPERDAY)
    fulldata.fill(mdi)
    fulldata[mask] = data

    return fulltimes, fulldata, mask # CalculateFullTimes

#***************************************
def fix_limits(value, upper = True):

    if value > 0:
        if upper: 
            value = value * 0.9
        else:
            value = value * 1.1
    if value < 0:
        if upper: 
            value = value * 1.1
        else:
            value = value * 0.9

    return value.round() # fix_limits
 
#***************************************
HOMOGLOCATION='/data/local/rdunn/Homogenisation/code/trunk/'

if var in ["relative_humidity", "specific_humidity", "vapor_pressure", "saturation_vapor_pressure", "wet_bulb_temperature"]:
    DATASUFFIX = "_humidity.nc"
else:
    DATASUFFIX='.nc'
    
RawStationData = ReadStations('{}/{}'.format(INPUT_FILE_LOCS, STATION_LIST))

# Read in the data
truetimes, data, mdi, input_stn = DataRead("{}/hadisd.{}_19310101-{}_{}{}".format(NETCDF_DATA_LOCS, LONG_VERSION, END_TIME, stn, DATASUFFIX), var)

# Sort the input stations to show merged stations
# this bit not used at the moment 27-Nov-2017
input_stns = []
for ts in range(input_stn.shape[0]):
    input_stns += ["".join(input_stn[ts,:])]
unique = set(input_stns)

# Uncompress the time axis
fulltimes, fulldata, mask = CalculateFullTimes(truetimes, data, mdi)
fulldata = np.ma.masked_where(fulldata <= -1.e30, fulldata)

# uncompress the input stations
source_stns = np.chararray(fulltimes.shape, itemsize = 12)
source_stns[mask] = input_stns

# convert to datetimes
plottimes = np.array([DATASTART + dt.timedelta(hours = n) for n in range(len(fulltimes))])

# choose a leap year
dummy_year = np.array([dt.datetime(2016,1,1,0,0) + dt.timedelta(hours = n) for n in range(366*24)])

# plot the data
fig = plt.figure(figsize = (20, 6))
ax = plt.axes([0.05,0.1,0.93,0.85])

thirty_mins = dt.timedelta(minutes = 30)
cmap = plt.cm.spring

maximum = fix_limits(np.ma.max(fulldata))
minimum = fix_limits(np.ma.min(fulldata), upper = False)

minimum = -10
maximum = 25

normal = plt.normalize(minimum, maximum)
cmap = plt.cm.gnuplot
colours = cmap(normal(fulldata))

for year in range(DATASTART.year, DATAEND.year):
    print year

    this_year_locs, = np.where(np.logical_and(plottimes >= dt.datetime(year,1,1,0,0), plottimes < dt.datetime(year+1, 1, 1, 0, 0)))

    if this_year_locs.shape[0] == 0: continue

    for o, obs in enumerate(this_year_locs):

        if fulldata.mask[obs] == True: continue

        left = mdates.date2num(dummy_year[o]-thirty_mins)
        right = mdates.date2num(dummy_year[o]+thirty_mins)
        width = right - left
        
        rect = patches.Rectangle((left, year-0.4), width, 0.8, color = colours[obs])
        ax.add_patch(rect)
        
#    if year == 1980: break

# fake a colourbar
sm = plt.cm.ScalarMappable(cmap=cmap, norm=normal)
sm._A = []
plt.colorbar(sm, label="Temperatures (C)", orientation = "vertical", pad=0.03, fraction=0.03, ticks=[-10,5,0,5,10,15,20,25], shrink = 0.8, extend = "both")

# sort the axes
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
ax.yaxis.set_ticks_position("left")
ax.xaxis.set_ticks_position("bottom")
ax.spines["bottom"].set_position(("axes",-0.01))
ax.spines["left"].set_position(("axes", -0.01))

# sort labels
monthLocs = mdates.MonthLocator()  # every month
nameLocs = mdates.DayLocator(bymonthday=[15])  # every month
monthFmt = mdates.DateFormatter('%b')

# tick locations
locator = mdates.AutoDateLocator(minticks=3)
ax.xaxis.set_major_locator(monthLocs) # major ticks to separate months
ax.xaxis.set_minor_locator(nameLocs) # minor ticks to name months

# tick label formats
formatter = mdates.AutoDateFormatter(locator)
ax.xaxis.set_major_formatter(ticker.NullFormatter()) # no label on majors
ax.xaxis.set_minor_formatter(monthFmt) # label minors
plt.tick_params(which='minor', length=0)

# set the limits
start = mdates.date2num(dummy_year[0]-thirty_mins)
end = mdates.date2num(dummy_year[-1]-thirty_mins)
plt.xlim([start-width, end+width])
plt.ylim([1972, DATAEND.year + 1])

# and labels
plt.title("Station ID - {}".format(stn))
watermarkstring="/".join(os.getcwd().split('/')[4:])+'/'+os.path.basename( __file__ )+"   "+dt.datetime.strftime(dt.datetime.now(), "%d-%b-%Y %H:%M")
plt.figtext(0.01,0.01,watermarkstring,size=7)
sourcestring = "HadISD {}".format(version)
plt.figtext(0.9,0.01,sourcestring,size=7)

plt.savefig("{}/all_hours_{}.png".format(IMAGE_LOCS, stn))

#************************************************************************
#                                 END
#************************************************************************
