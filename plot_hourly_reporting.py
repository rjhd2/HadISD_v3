#!/usr/local/sci/bin/python
#***************************************
#     Calculates average reporting rate
#     Per station per year and plots results
#     May 2013 RJHD
#***************************************
#************************************************************************
#                    SVN Info
# $Rev:: 77                                            $:  Revision of last commit
# $Author:: rdunn                                      $:  Author of last commit
# $Date:: 2016-11-21 17:51:26 +0000 (Mon, 21 Nov 2016) $:  Date of last commit
#************************************************************************
#                                 START
#************************************************************************
import numpy as np
import datetime as dt
import netCDF4 as ncdf
import matplotlib.pyplot as plt
import os

# RJHD utils
from set_paths_and_vars import *

DAYSPERYEAR=365.25
MONTHSPERYEAR=12
HOURSPERDAY=24


#***************************************
# Definitions/Subroutines
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

        masked_data=var[:].data
        time_data=timevar[:]

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
    except IOError:
        print "File not found: ",filename
        raise IOError

    return times,data,mdi # DataRead

#***************************************
def YearBreaks():
    '''
    Find hour number for at each month change over the years.

    Inputs: 
    :returns: 12xNyears array giving start of each month in days since startyear
    :requires: DATASTART and DATAEND - datetimes of start and end of dataset
    '''
 

    Months=range(1,13) # MAM,JJA,SON,DJF & year start

    YearBreaks=[]
    year=DATASTART.year
    while year < DATAEND.year+1:
        difference=dt.datetime(year,1,1,0,0)-DATASTART
        
        YearBreaks+=[float(difference.days)*HOURSPERDAY]

        year+=1

    return YearBreaks # MonthBreaks

#------------------------------------------------------------

var = "temperatures"

StationData = ReadStations('{}/{}'.format(INPUT_FILE_LOCS, STATION_LIST))
StationIDs=np.array(RawStationData[:,0])

StationData=np.array(StationData)
StationLat=np.array([float(x) for x in StationData[:,1]])
StationLon=np.array([float(x) for x in StationData[:,2]])
StationElv=np.array([float(x) for x in StationData[:,3]])

nstations=len(StationIDs)
nyears=DATAEND.year-DATASTART.year

YearStarts=YearBreaks()

Resolution=np.zeros([nstations,nyears+1])

for s,stn in enumerate(StationIDs):
    print stn

    # read data - mask on read
    try:
        truetimes, data, mdi=DataRead("{}/hadisd.{}_19310101-{}_{}{}".format(NETCDF_DATA_LOCS, LONG_VERSION, END_TIME, stn, DATASUFFIX), var, mask = True)
        
    except IOError:
        print "File for station not found.  Skipping station ", StationIDs[st]
        continue

    # convert the year breaks onto truetimes axis
    YearIndices=[]         
    for ys in YearStarts:
        match=np.where(truetimes>ys)
        
        if len(match[0]) != 0:
            YearIndices+=[match[0][0]]
        else:
            YearIndices+=[len(truetimes)-1]

    for y,year in enumerate(YearIndices):

        if y == len(YearIndices)-1:
            this_years_times=truetimes[year:]
        else:
            this_years_times=truetimes[year:YearIndices[y+1]]

            
        if len(this_years_times) > 2:

            first_diff=this_years_times[1:]-this_years_times[:-1]

            Resolution[s,y]=np.median(first_diff)

            
results=np.zeros([5,nyears])

for y in range(nyears):
    hourly=len(np.where(Resolution[:,y] == 1)[0])
    threes=len(np.where(Resolution[:,y] == 3)[0])
    sixes=len(np.where(Resolution[:,y] == 6)[0])
    others=len(np.where(Resolution[:,y] != 0)[0])-hourly-threes-sixes
    off=len(np.where(Resolution[:,y] == 0)[0])
    results[:,y]=[hourly,threes,sixes,others,off]

plt.clf()

plt.plot(np.arange(nyears)+DATASTART.year,results[4,:],'b-',label='none')
plt.plot(np.arange(nyears)+DATASTART.year,results[3,:],'c-',label='unknown')
plt.plot(np.arange(nyears)+DATASTART.year,results[2,:],c='lime',ls='-',label='6 hourly')
plt.plot(np.arange(nyears)+DATASTART.year,results[1,:],c='orange',ls='-',label='3 hourly')
plt.plot(np.arange(nyears)+DATASTART.year,results[0,:],'r-',label='1 hourly')

plt.ylabel('stations')
plt.xlim([DATASTART.year,DATAEND.year])
watermarkstring="/".join(os.getcwd().split('/')[4:])+'/'+os.path.basename( __file__ )+"   "+dt.datetime.strftime(dt.datetime.now(), "%d-%b-%Y %H:%M")
plt.figtext(0.01,0.01,watermarkstring,size=6)
leg=plt.legend(loc='lower center',ncol=5, bbox_to_anchor=(0.5,-0.1),frameon=False, title='',prop={'size':10})

plt.savefig('{}/HadISD_reporting_{}_{}.png'.format(IMAGE_LOCS, HADISD_VERSION, var))

plt.clf()
plt.fill_between(np.arange(nyears)+DATASTART.year,results[4,:],color='b',label='none')
plt.fill_between(np.arange(nyears)+DATASTART.year,results[4,:],results[3,:]+results[4,:],color='c',label='unknown')
plt.fill_between(np.arange(nyears)+DATASTART.year,results[3,:]+results[4,:],results[2,:]+results[3,:]+results[4,:],color='lime',label='6 hourly')
plt.fill_between(np.arange(nyears)+DATASTART.year,results[2,:]+results[3,:]+results[4,:],results[1,:]+results[2,:]+results[3,:]+results[4,:],color='orange',label='3 hourly')
plt.fill_between(np.arange(nyears)+DATASTART.year,results[1,:]+results[2,:]+results[3,:]+results[4,:],results[0,:]+results[1,:]+results[2,:]+results[3,:]+results[4,:],color='r',label='1 hourly')
plt.plot(np.arange(nyears)+DATASTART.year,results[4,:],'b-',label='none')
plt.plot(np.arange(nyears)+DATASTART.year,results[3,:]+results[4,:],'c-',label='unknown')
plt.plot(np.arange(nyears)+DATASTART.year,results[2,:]+results[3,:]+results[4,:],c='lime',ls='-',label='6 hourly')
plt.plot(np.arange(nyears)+DATASTART.year,results[1,:]+results[2,:]+results[3,:]+results[4,:],c='orange',ls='-',label='3 hourly')
plt.plot(np.arange(nyears)+DATASTART.year,results[0,:]+results[1,:]+results[2,:]+results[3,:]+results[4,:],'r-',label='1 hourly')
plt.plot(np.arange(nyears)+DATASTART.year,results[4,:],'k-')
plt.plot(np.arange(nyears)+DATASTART.year,results[3,:]+results[4,:],'k-')
plt.plot(np.arange(nyears)+DATASTART.year,results[2,:]+results[3,:]+results[4,:],'k-')
plt.plot(np.arange(nyears)+DATASTART.year,results[1,:]+results[2,:]+results[3,:]+results[4,:],'k-')
plt.plot(np.arange(nyears)+DATASTART.year,results[0,:]+results[1,:]+results[2,:]+results[3,:]+results[4,:],'k-')

plt.ylabel('stations')
plt.xlim([DATASTART.year,DATAEND.year])
watermarkstring="/".join(os.getcwd().split('/')[4:])+'/'+os.path.basename( __file__ )+"   "+dt.datetime.strftime(dt.datetime.now(), "%d-%b-%Y %H:%M")
plt.figtext(0.01,0.01,watermarkstring,size=6)
leg=plt.legend(loc='lower center',ncol=5, bbox_to_anchor=(0.5,-0.1),frameon=False, title='',prop={'size':10})

plt.savefig('{}/HadISD_reporting_{}_{}_filled.png'.format(IMAGE_LOCS, HADISD_VERSION, var))


total=np.sum(results[:4], axis=0)
results=results[:4]/total

plt.clf()
plt.fill_between(np.arange(nyears)+DATASTART.year,results[3,:],color='c',label='unknown')
plt.fill_between(np.arange(nyears)+DATASTART.year,results[3,:],results[2,:]+results[3,:],color='lime',label='6 hourly')
plt.fill_between(np.arange(nyears)+DATASTART.year,results[2,:]+results[3,:],results[1,:]+results[2,:]+results[3,:],color='orange',label='3 hourly')
plt.fill_between(np.arange(nyears)+DATASTART.year,results[1,:]+results[2,:]+results[3,:],results[0,:]+results[1,:]+results[2,:]+results[3,:],color='r',label='1 hourly')
plt.plot(np.arange(nyears)+DATASTART.year,results[3,:],'c-',label='unknown')
plt.plot(np.arange(nyears)+DATASTART.year,results[2,:]+results[3,:],c='lime',ls='-',label='6 hourly')
plt.plot(np.arange(nyears)+DATASTART.year,results[1,:]+results[2,:]+results[3,:],c='orange',ls='-',label='3 hourly')
plt.plot(np.arange(nyears)+DATASTART.year,results[0,:]+results[1,:]+results[2,:]+results[3,:],'r-',label='1 hourly')
plt.plot(np.arange(nyears)+DATASTART.year,results[3,:],'k-')
plt.plot(np.arange(nyears)+DATASTART.year,results[2,:]+results[3,:],'k-')
plt.plot(np.arange(nyears)+DATASTART.year,results[1,:]+results[2,:]+results[3,:],'k-')
plt.plot(np.arange(nyears)+DATASTART.year,results[0,:]+results[1,:]+results[2,:]+results[3,:],'k-')


plt.ylabel('proportion')
plt.xlim([DATASTART.year,DATAEND.year])
watermarkstring="/".join(os.getcwd().split('/')[4:])+'/'+os.path.basename( __file__ )+"   "+dt.datetime.strftime(dt.datetime.now(), "%d-%b-%Y %H:%M")
plt.figtext(0.01,0.01,watermarkstring,size=6)
leg=plt.legend(loc='lower center',ncol=5, bbox_to_anchor=(0.5,-0.1),frameon=False, title='',prop={'size':10})

plt.savefig('{}/HadISD_reporting_{}_{}_filled__proportion.png'.format(IMAGE_LOCS, HADISD_VERSION, var))

#************************************************************************
#                                 END
#************************************************************************
