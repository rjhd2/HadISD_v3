#!/usr/local/sci/bin/python
#************************************************************************
#                    SVN Info
#$Rev:: 61                                            $:  Revision of last commit
#$Author:: rdunn                                      $:  Author of last commit
#$Date:: 2015-04-28 18:11:08 +0100 (Tue, 28 Apr 2015) $:  Date of last commit
#************************************************************************

'''
Python script to read the ISD ASCII text format and output netcdf files

Runs with no inputs in current version. 

Compared to IDL output using compare.py on 30May2012 and found to match
except for total_cloud_flags - but on investigation with raw ISD files
the python extraction is the correct one. RJHD

Could change data types to match IDL, but adapting QC so that works with
floats and doubles as appropriate.RJHD

Logic change to match IDL so that overwrite only acted upon if writing
real data, rather than missing. RJHD

'''


import numpy as np
import datetime as dt
import glob
import gzip
import subprocess
import math
import netCDF4 as ncdf
import sys
import os
import argparse
import datetime, calendar

# Globals
INPUT_DATA_DIR =r'/project/hadobs2/hadisd/v200_2014/isd_files_v200_2014/'
OUTPUT_DATA_DIR=r'/project/hadobs2/hadisd/v200_2014/netcdf_files_v200_2014/'
OUTPUT_FILE_DIR=r'/project/hadobs2/hadisd/v200_2014/suppl_files_v200_2014/'

FILE_LOCS = r'input_files/'

station_list_filename= 'candidate_stations.txt'
merger_list_filename = 'final_mergers.txt'

INTMDI=-999
FLTMDI=-1.e30
hours = True # output time axis in hours.

NCDC_FLAGS={'A':10,'U':11,'P':12,'I':13,'M':14,'C':15,'R':16}

#---------------------------------------------------------------------


#************************************************************************
def ReadStations(filename):
    ''' 
    Read Station Information

    :param string filename: name and location of input file

    :returns: numpy array of file contents

    Use numpy genfromtxt reading to read in all station 
    data in ID,Lat,Lon,Elev list
    '''
    return np.genfromtxt(filename, dtype=(str)) #  ReadStations


#************************************************************************
def ReadComposites(filename):
    ''' 
    Read Composite Station Information

    :param string filename: name and location of input file

    :returns: list of lists containing composites
    '''
    composites=[]
    try:
        with open(filename) as infile:
            for line in infile:
                split_line=line.split()
                composites.append(split_line[:])                           
    except IOError:
        print "File not found: ",filename
        raise IOError
    return composites # ReadComposites

#************************************************************************
def RepresentsInt(s):
    '''
    Tests if string is an integer

    :param string s: string to test
    :returns: boolean if string is valid integer
    '''
    try: 
        int(s)
        return True
    except ValueError:
        return False # RepresentsInt

#************************************************************************
def TimeMatch(timearray,testtime, lower,upper):
    '''
    Do matching of np array to find time step

    :param array timearray: np array of timestamps
    :param float testtime: timestep to find
    :return: int of location
    '''
    return np.argwhere(timearray[lower:upper]==testtime)[0] # TimeMatch


#************************************************************************
def ExtractValues(missing,line,location,length,test,divisor=1.,flagoffset=0, doflag=True):
    '''
    Extract the appropriate values from the line string.

    Assumes that usually it is a value,flag pair, with no scaling factor and that
    the flag follows directly on from the value.  Can be adjusted.
    
    :param float/int missing: mdi
    :param string line: input line from file
    :param int location: location of start of string subset to read in
    :param int length: length of string to convert to value
    :param int test: value of string if data missing
    :param float divisor: scaling factor for data, default=1
    :param int flagoffset: shift of flag from end of data, default=0
    :param boolean doflag: to extract a flag value, default=True

    :returns: tuple of value, flag OR value if doflag=False
    '''

    temp_value=line[location:location+length]

    value=missing
    flag=INTMDI

    if temp_value != test:
        if missing==FLTMDI:
            value=float(temp_value)/divisor
        elif missing==INTMDI:
            value=int(temp_value)/divisor
        elif missing=='':
            value=temp_value

        if doflag:
            
            flag=line[location+length+flagoffset:location+length+flagoffset+1]

            if RepresentsInt(flag):
                flag=int(flag)

            else:
                try:
                    flag=NCDC_FLAGS[flag]
                except KeyError:
                    print 'ALPHA FLAG CONVERSION FAILED'
                    print 'input flag is: '+flag
                    print 'Line in raw ISD record reads:'
                    print line
                    flag=20


    if doflag:
        return value,flag
    else:
        return value # ExtractValues


#************************************************************************
def TestToExtract(data,missing,overwrite):
    '''
    Test if need to extract the data

    :param float/int data: data to test
    :param float/int missing: missing value
    :param boolean overwrite: to overwrite or not

    :returns: boolean if condition met
    '''

    if data==missing or overwrite:
        return True
    else:
        return False # TestToExtract


#************************************************************************
def ExtractionProcess(data, flags, time, missing, missingtest, line, location, length,divisor=1.,flagoffset=0, doflag=True):
    '''
    Run the extraction, and write the values if the extracted ones are not empty

    :param array data: the data array
    :param array flags: the flags array
    :param int time: the time stamp    
    :param float/int missing: mdi
    :param float/int missingtest: value of string if missing
    :param string line: input line from file
    :param int location: location of start of string subset to read in
    :param int length: length of string to convert to value
    :param int test: value of string if data missing
    :param float divisor: scaling factor for data, default=1
    :param int flagoffset: shift of flag from end of data, default=0
    :param boolean doflag: to extract a flag value, default=True
    '''

    if doflag:
        value,flag=ExtractValues(missing,line,location,length,missingtest,divisor=divisor,flagoffset=flagoffset,doflag=doflag)
        # no longer want to test that there was a value - as all data taken from 
        # single observation, regardless of who complete it is.
        # left old code extant in case changes to happen in future
        if value != missing:
            data[time],flags[time]=value,flag
        else:
            data[time],flags[time]=value,flag
    
    else:
        value=ExtractValues(missing,line,location,length,missingtest,divisor=divisor,flagoffset=flagoffset,doflag=doflag)
        if value != missing:
            data[time]=value
        else:
            data[time]=value
        
    return # ExtractionProcess

#************************************************************************
def WriteDubious(outfile,infile,code, station, time):
    '''
    Write note to dubious file list.
    
    :param string outfile: filename to be written to
    :param string infile: filename of dubious file
    :param string code: text identifier of variables being tested
    :param string station: station ID being processed
    :param string time: time of the dubious data

    :returns: int of flag status.
    '''
    flagged=0
    try:
        with open(outfile,'a') as of:
            of.write(station+' '+time+' '+code+' variables are first, but not nec. only problem '+infile+'\n')
            of.close()
            flagged=1
    except IOError:
        with open(outfile,'w') as of:
            of.write(station+' '+time+' '+code+' variables are first, but not nec. only problem '+infile+'\n')
            of.close()
            flagged=1
        
        
    return flagged # WriteDubious

#************************************************************************
def SortClouds(cloud_cover,cloud_flags, time, amounts, flags, clouds):
    '''
    Convert the raw cloud data into oktas for each level

    :param array cloud_cover: final cloud_cover array
    :param array cloud_flags: final cloud_flags array
    :param int time_loc: time stamp
    :param array amounts: raw cloud amounts - in oktas
    :param array flags: raw cloud flags
    :param array clouds: locations of where cloud heights match this level
    '''

    if len(clouds)>=1 and cloud_cover[time]==INTMDI:
        cloud_cover[time]=np.max(amounts[clouds])
        cloud_flags[time]=np.max(flags[clouds])

    return # SortClouds                              

#************************************************************************
def SortClouds2(cloud_cover,cloud_flags, time, amounts, amounts2, flags, clouds):
    '''
    Convert the raw cloud data into oktas and for each level

    :param array cloud_cover: final cloud_cover array
    :param array cloud_flags: final cloud_flags array
    :param int time_loc: time stamp
    :param array amounts: raw cloud amounts - in other units - see ISD documentation
    :param array amounts2: raw cloud amounts - in oktas
    :param array flags: raw cloud flags
    :param array clouds: locations of where cloud heights match this level
    '''

    inoktas=np.where(np.array(amounts2[clouds]) != INTMDI)[0]
    if len(inoktas)>=1 and cloud_cover[time]==INTMDI:
        cloud_cover[time]=np.max(amounts[clouds][inoktas])
        cloud_flags[time]=np.max(flags[clouds][inoktas])
    elif cloud_cover[time]==INTMDI:
        # convert to oktas
        cloud_cover[time]=np.max(amounts[clouds])*2.
        cloud_flags[time]=np.max(flags[clouds])

    return # SortClouds2                         

#************************************************************************
def WriteAttributes(variable,long_name,missing_value,units,axis,vmin,vmax,standard_name = ''):
    '''
    Write given attributes into ncdf variable
    
    :param object variable: netcdf Variable
    :param string long_name: long_name value for variable to be written
    :param float/int missing_value: missing_value value for variable to be written
    :param string units: units value for variable to be written
    :param string axis: axis value for variable to be written
    :param float/int vmin: valid_min value for variable to be written
    :param float/int vmax: valid_max value for variable to be written
    :param string standard_name: standard_name value for variable to be written

    '''
    variable.long_name=long_name
    variable.missing_value=missing_value
    variable.axis=axis
    variable.units=units
    variable.valid_min=vmin
    variable.valid_max=vmax
    
    if standard_name != '':
        variable.standard_name=standard_name

    return # WriteAttributes

#************************************************************************
def WriteFlagAttributes(variable,long_name,missing_value,axis):
    '''
    Write given attributes into ncdf variable
    
    :param object variable: netcdf Variable
    :param string long_name: long_name value for variable to be written
    :param float/int missing_value: missing_value value for variable to be written
    :param string axis: axis value for variable to be written
    '''
    variable.long_name=long_name
    variable.missing_value=missing_value
    variable.units="no_unit"
    variable.axis=axis

    return # WriteFlagAttributes


#************************************************************************
def MakeNetcdfFiles(STARTYEAR, ENDYEAR, restart_id="", end_id="", do_zip = True, Extra=False): 
    '''
    Parse the ASCII files and do the NetCDF file creation

    :param string restart_id: string for starting station, default=""
    :param string end_id: string for ending station, default=""
    :param boolean do_zip: make netCDF4 files with internal zipping
    :param boolean Extra: setting to extract extra variables
    '''

    StationInfo=ReadStations(FILE_LOCS+station_list_filename)

    StationIDs=np.array(StationInfo[:,0])

    # sorted in case of unordered lists
    sort_order=np.argsort(StationIDs)

    StationIDs=StationIDs[sort_order]

    StationLat=np.array([float(x) for x in StationInfo[sort_order,1]])
    StationLon=np.array([float(x) for x in StationInfo[sort_order,2]])
    StationElv=np.array([float(x) for x in StationInfo[sort_order,3]])

    print "Read in %i stations" % len(StationIDs)

    # reduce station list to start and end stations

    if restart_id != "":
        startindex=np.where(StationIDs==restart_id)[0][0]
        StationIDs=StationIDs[startindex:]
        StationLat=StationLat[startindex:]
        StationLon=StationLon[startindex:]
        StationElv=StationElv[startindex:]
    if end_id != "":
        endindex=np.where(StationIDs==end_id)[0][0]
        StationIDs=StationIDs[:endindex+1]
        StationLat=StationLat[:endindex+1]
        StationLon=StationLon[:endindex+1]
        StationElv=StationElv[:endindex+1]

    if restart_id !="" or end_id!="":
        print "Truncated run selected"
        print "Processing %i stations" % len(StationIDs)

    nstations=len(StationIDs)

    Composites=ReadComposites(FILE_LOCS+merger_list_filename)

    DaysBetween=dt.datetime(ENDYEAR+1,1,1,0,0)-dt.datetime(STARTYEAR,1,1,0,0)
    HoursBetween=int(DaysBetween.days*24.)

    TimeStamps=np.linspace(0,HoursBetween-1,HoursBetween) # keep in integer hours

    ValidYears=np.linspace(STARTYEAR,ENDYEAR,(ENDYEAR-STARTYEAR+1))
    dubiousfile=OUTPUT_FILE_DIR+'dubious_ISD_data_files.txt'

    # read in Canadian station list
    Canadian_stations_info = np.genfromtxt(FILE_LOCS + "Canada_time_ranges.dat", dtype=(str), delimiter = [12,20,20])
    Canadian_station_ids = Canadian_stations_info[:,0]
    Canadian_station_start = np.array([dt.datetime.strptime(d.strip(), "%Y-%m-%d %H:%M:%S") for d in Canadian_stations_info[:,1]])
    Canadian_station_end   = np.array([dt.datetime.strptime(d.strip(), "%Y-%m-%d %H:%M:%S") for d in Canadian_stations_info[:,2]])
  
    dbg_sttime=dt.datetime.now()

    for st,station in enumerate(StationIDs):

        print '%s, number %i of %i' %(station, st+1, nstations)

        temperatures=np.zeros(HoursBetween)
        temperature_flags=np.zeros(HoursBetween, dtype=np.int)
        dewpoints=np.zeros(HoursBetween)
        dewpoint_flags=np.zeros(HoursBetween, dtype=np.int)
        total_cloud_cover=np.zeros(HoursBetween, dtype=np.int)
        total_cloud_flags=np.zeros(HoursBetween, dtype=np.int)
        low_cloud_cover=np.zeros(HoursBetween, dtype=np.int)
        low_cloud_flags=np.zeros(HoursBetween, dtype=np.int)
        mid_cloud_cover=np.zeros(HoursBetween, dtype=np.int)
        mid_cloud_flags=np.zeros(HoursBetween, dtype=np.int)
        high_cloud_cover=np.zeros(HoursBetween, dtype=np.int)
        high_cloud_flags=np.zeros(HoursBetween, dtype=np.int)
        cloud_base=np.zeros(HoursBetween)
        cloud_base_flags=np.zeros(HoursBetween, dtype=np.int)
        windspeeds=np.zeros(HoursBetween)
        windspeeds_flags=np.zeros(HoursBetween, dtype=np.int)
        winddirs=np.zeros(HoursBetween, dtype=np.int)
        winddirs_flags=np.zeros(HoursBetween, dtype=np.int)
        past_sigwx1=np.zeros(HoursBetween, dtype=np.int)
        past_sigwx1_period=np.zeros(HoursBetween, dtype=np.int)
        past_sigwx1_flag=np.zeros(HoursBetween, dtype=np.int)
        precip1_period=np.zeros(HoursBetween, dtype=np.int)
        precip1_depth=np.zeros(HoursBetween)
        precip1_condition=['null' for i in range(HoursBetween)]
        precip1_flag=np.zeros(HoursBetween, dtype=np.int)
        slp=np.zeros(HoursBetween)
        slp_flag=np.zeros(HoursBetween, dtype=np.int)
        sun_duration=np.zeros(HoursBetween)
        sun_durationqc=np.zeros(HoursBetween, dtype=np.int)
        wind_gust_period=np.zeros(HoursBetween)
        wind_gust_value=np.zeros(HoursBetween)
        wind_gust_flags=np.zeros(HoursBetween, dtype=np.int)

        # Tells you what the true input station id was for the duplicate
        # using list as string array.
        input_station_id=['null' for i in range(HoursBetween)]


        temperatures.fill(FLTMDI)
        temperature_flags.fill(INTMDI)
        dewpoints.fill(FLTMDI)
        dewpoint_flags.fill(INTMDI)
        total_cloud_cover.fill(INTMDI)
        total_cloud_flags.fill(INTMDI)
        low_cloud_cover.fill(INTMDI)
        low_cloud_flags.fill(INTMDI)
        mid_cloud_cover.fill(INTMDI)
        mid_cloud_flags.fill(INTMDI)
        high_cloud_cover.fill(INTMDI)
        high_cloud_flags.fill(INTMDI)
        cloud_base.fill(INTMDI)
        cloud_base_flags.fill(INTMDI)
        windspeeds.fill(FLTMDI)
        windspeeds_flags.fill(INTMDI)
        winddirs.fill(INTMDI)
        winddirs_flags.fill(INTMDI)
        past_sigwx1.fill(INTMDI)
        past_sigwx1_period.fill(INTMDI)
        past_sigwx1_flag.fill(INTMDI)
        precip1_period.fill(INTMDI)
        precip1_depth.fill(FLTMDI)
        precip1_flag.fill(INTMDI)
        slp.fill(FLTMDI)
        slp_flag.fill(INTMDI)
        sun_duration.fill(INTMDI)
        sun_durationqc.fill(INTMDI)
        wind_gust_period.fill(FLTMDI)
        wind_gust_value.fill(FLTMDI)
        wind_gust_flags.fill(INTMDI)

        if Extra:
            windtypes=['null' for i in range(HoursBetween)]
            present_sigwx=np.zeros(HoursBetween, dtype=np.int)
            present_sigwx_flags=np.zeros(HoursBetween, dtype=np.int)
            past_sigwx2=np.zeros(HoursBetween, dtype=np.int)
            past_sigwx2_period=np.zeros(HoursBetween, dtype=np.int)
            past_sigwx2_flag=np.zeros(HoursBetween, dtype=np.int)
            precip2_period=np.zeros(HoursBetween, dtype=np.int)
            precip2_depth=np.zeros(HoursBetween)
            precip2_condition=['null' for i in range(HoursBetween)]
            precip2_flag=np.zeros(HoursBetween, dtype=np.int)
            precip3_period=np.zeros(HoursBetween, dtype=np.int)
            precip3_depth=np.zeros(HoursBetween)
            precip3_condition=['null' for i in range(HoursBetween)]
            precip3_flag=np.zeros(HoursBetween, dtype=np.int)
            precip4_period=np.zeros(HoursBetween, dtype=np.int)
            precip4_depth=np.zeros(HoursBetween)
            precip4_condition=['null' for i in range(HoursBetween)]
            precip4_flag=np.zeros(HoursBetween, dtype=np.int)
            maximum_temp_period=np.zeros(HoursBetween)
            maximum_temp_value=np.zeros(HoursBetween)
            maximum_temp_flags=np.zeros(HoursBetween, dtype=np.int)
            minimum_temp_period=np.zeros(HoursBetween)
            minimum_temp_value=np.zeros(HoursBetween)
            minimum_temp_flags=np.zeros(HoursBetween, dtype=np.int)

            present_sigwx.fill(INTMDI)
            present_sigwx_flags.fill(INTMDI)
            past_sigwx2.fill(INTMDI)
            past_sigwx2_period.fill(INTMDI)
            past_sigwx2_flag.fill(INTMDI)
            precip2_period.fill(INTMDI)
            precip2_depth.fill(FLTMDI)
            precip2_flag.fill(INTMDI)
            precip3_period.fill(INTMDI)
            precip3_depth.fill(FLTMDI)
            precip3_flag.fill(INTMDI)
            precip4_period.fill(INTMDI)
            precip4_depth.fill(FLTMDI)
            precip4_flag.fill(INTMDI)  
            maximum_temp_period.fill(FLTMDI)
            maximum_temp_value.fill(FLTMDI)
            maximum_temp_flags.fill(INTMDI)
            minimum_temp_period.fill(FLTMDI)
            minimum_temp_value.fill(FLTMDI)
            minimum_temp_flags.fill(INTMDI)




        # extract stations to process, including composites.
        is_composite=next((i for i, sublist in enumerate(Composites) if station in sublist), -1)
        if is_composite!=-1:
            consider_these=Composites[is_composite]
            print 'This is a duplicate station containing %s ' % ' '.join(consider_these)
        else:
            consider_these=[station]

        # get listing of all files to process
        raw_files=[]
        for cstn in consider_these:
            if cstn[0:3] >= '725' and cstn[0:3] <= '729':
                raw_files.extend(glob.glob(INPUT_DATA_DIR+'station725/'+cstn+'*'))
            else:
                raw_files.extend(glob.glob(INPUT_DATA_DIR+'station'+cstn[0:2]+'s/'+cstn+'*'))


        raw_files.sort()

        dbg_lasttime=dt.datetime.now()

        for rfile in raw_files:
            done_print = False # for output of Canadian station skipping

            a=dt.datetime.now()-dbg_lasttime
            print rfile, a

            dbg_lasttime=dt.datetime.now()

            raw_station=rfile.split('/')[-1][0:12]

            rfile_year=int(rfile.split('-')[-1].split('.')[0])

            rfile_days=dt.datetime(rfile_year,1,1,0,0)-dt.datetime(STARTYEAR,1,1,0,0)
            rfile_hours=rfile_days.days*24.

            rfile_ydays=dt.datetime(rfile_year+1,1,1,0,0)-dt.datetime(rfile_year,1,1,0,0)
            rfile_yhours=rfile_ydays.days*24

            if rfile_year in ValidYears:
                dubious_flagged=0

                if rfile[-2:]!='gz':
                    subprocess.call(['gzip','-f','-9',rfile])
                    rfile=rfile+'.gz'
                    # note - this amends the file identifier in the loop


                last_obs_time=0.
                
                try:
                    with gzip.open(rfile,'r') as infile:
                        for rawline in infile:

                            # main processing

                            # check for remarks
                            cleanline=rawline[0:rawline.find('REM')]

                            # find which timestamp we're working at
                            year=int(cleanline[15:19])
                            month=int(cleanline[19:21])
                            day=int(cleanline[21:23])
                            hour=int(cleanline[23:25])
                            minute=int(cleanline[25:27])

                            # found error in minute value in 030910-99999, 035623-99999
                            if minute < 0:
                                hour=hour-1
                                minute=60+minute
                            elif minute > 59:
                                hour=hour+1
                                minute=minute-60
                            if hour < 0:
                                day=day-1
                                hour=24+hour
                            elif hour > 23:
                                day=day+1
                                hour=hour-24
                            
                            dummy, ndays = calendar.monthrange(year, month)
                            if day <= 0:
                                month = month -1
                                dummy, ndays = calendar.monthrange(year, month)
                                day=ndays - day
                            elif day > ndays:
                                month=month+1
                                day=day - ndays

                            if month <= 0:
                                month = 12 - month
                                year = year - 1
                            elif month > 12:
                                month=month - 12
                                year = year + 1

                            dt_time = dt.datetime(year, month, day, hour, minute)

                            if raw_station in Canadian_station_ids:
                                # then test for restrictions on start/end times
                                loc, = np.where(Canadian_station_ids == raw_station)

                                if dt_time < Canadian_station_start[loc[0]]:
                                    if not done_print:
                                        print "skipping year {} of station {} as identified as undocumented move by Environment Canada".format(year, raw_station)
                                        done_print = True
                                    continue
                                if dt_time > Canadian_station_end[loc[0]]:
                                    if not done_print:
                                        print "skipping year {} of station {} as identified as undocumented move by Environment Canada".format(year, raw_station)
                                        done_print = True
                                    continue


                            # integer hours
                            obs_time=ncdf.date2num(dt_time, units='hours since '+str(STARTYEAR)+'-01-01 00:00:00', calendar='julian')

                            string_obs_time=dt.datetime.strftime(dt_time,"%d-%m-%Y, %H:%M")

                            # working in hours, so just round the hours since
                            # start date to get the time stamp 13/6/2012
                            time_loc=int(round(obs_time))

                            # test if this time_loc out of bounds:
                            #  e.g. if 2350 on 31/12/ENDYEAR then should not
                            #       takes the obs as it belongs to following year
                            if time_loc != HoursBetween:

                                # test if closer to timestamp than previous observation
                                # overwrite only acted upon if newer real data closer to full hour
                                currentT=temperatures[time_loc]
                                currentD=dewpoints[time_loc]

                                newT=ExtractValues(FLTMDI,cleanline,87,5,'+9999',divisor=10.,doflag=False)
                                newD=ExtractValues(FLTMDI,cleanline,93,5,'+9999',divisor=10.,doflag=False)

                                Extract=False

                                # no extract happened for this time stamp as yet - so extract
                                if input_station_id[time_loc] =='null':
                                    # this is not an overwrite, so nothing extra needs doing
                                    Extract=True

                                # if currently have no T or D data
                                elif currentT==FLTMDI and currentD==FLTMDI:
                                    # if overwriting, only do so with observation closer to time stamp
                                    if input_station_id[time_loc] == raw_station: # tests if already read into this time stamp
                                        if (newT != FLTMDI) or (newD != FLTMDI):
                                            # if updated data would have T or D, then take it, even if further from the time stamp
                                            Extract=True
                                        elif last_obs_time != 0.: # not really necessary as already have at least 1 obs, but still...
                                            # if time stamp closer than last one
                                            if np.abs(TimeStamps[time_loc]-obs_time) < np.abs(TimeStamps[time_loc]-last_obs_time):
                                                Extract=True
                                    # else just take the line - no observations read into this time stamp yet
                                    else:
                                        Extract=True

                                # if already have T but _no_ D OR D but _no_ T, but new one has T and D, take this line
                                #    this is an overwrite - so also check that overwriting with the same station
                                elif ((currentT!=FLTMDI and currentD==FLTMDI) or (currentT==FLTMDI and currentD!=FLTMDI)) \
                                        and (newT!=FLTMDI and newD!=FLTMDI):
                                    if input_station_id[time_loc] == raw_station:
                                        # preference to both values over just one
                                        Extract=True

                                # have D but no T, and new observation comes up with T, select if closer
                                elif (currentT==FLTMDI and currentD!=FLTMDI) and (newT!=FLTMDI):
                                    # if overwriting, only do so with observation closer to time stamp
                                    if input_station_id[time_loc] == raw_station: # tests if already read into this time stamp
                                        if last_obs_time != 0.: # not really necessary as already have at least 1 obs, but still...
                                            # if time stamp closer than last one
                                            if np.abs(TimeStamps[time_loc]-obs_time) < np.abs(TimeStamps[time_loc]-last_obs_time):
                                                Extract=True

                                # have T but no D, and new observation comes up with T, select if closer
                                elif (currentT!=FLTMDI and currentD==FLTMDI) and (newT!=FLTMDI):
                                    # if overwriting, only do so with observation closer to time stamp
                                    if input_station_id[time_loc] == raw_station: # tests if already read into this time stamp
                                        if last_obs_time != 0.: # not really necessary as already have at least 1 obs, but still...
                                            # if time stamp closer than last one
                                            if np.abs(TimeStamps[time_loc]-obs_time) < np.abs(TimeStamps[time_loc]-last_obs_time):
                                                Extract=True

                                # if already have T and D, and new one also has T and D, but at closer time stamp, take this line
                                #    this is an overwrite - so also check that overwriting with the same station
                                elif (currentT!=FLTMDI and currentD!=FLTMDI) and (newT!=FLTMDI and newD!=FLTMDI):
                                    if input_station_id[time_loc] == raw_station:
                                        if last_obs_time != 0.: # not really necessary as already have at least 1 obs, but still...
                                            # if time stamp closer than last one

                                            if np.abs(TimeStamps[time_loc]-obs_time) < np.abs(TimeStamps[time_loc]-last_obs_time):
                                                Extract=True

                                else:
                                    Extract=False # just in case

                                # sort last obs_time - 
                                last_obs_time=obs_time

                                if input_station_id[time_loc]=='null':
                                    input_station_id[time_loc]=raw_station                                  

                               # main variables
                                dummyflag=0
                                # if allowed to extract:
                                if Extract:
                                    ExtractionProcess(temperatures, temperature_flags,time_loc,FLTMDI,'+9999',cleanline,87,5, divisor=10.)

                                if Extract:
                                    ExtractionProcess(dewpoints, dewpoint_flags,time_loc,FLTMDI,'+9999',cleanline,93,5, divisor=10.)


                                if Extract:
                                    ExtractionProcess(slp, slp_flag,time_loc,FLTMDI,'99999',cleanline,99,5, divisor=10.)


                                if Extract:
                                    ExtractionProcess(winddirs, winddirs_flags,time_loc,INTMDI,'999',cleanline,60,3)
                                    if Extra:
                                        ExtractionProcess(windtypes, dummyflag, time_loc,'','-',cleanline,64,1,doflag=False)

                                if Extract:
                                    ExtractionProcess(windspeeds, windspeeds_flags,time_loc,FLTMDI,'9999',cleanline,65,4, divisor=10.)


                                if Extract:
                                    ExtractionProcess(cloud_base, cloud_base_flags,time_loc,INTMDI,'99999',cleanline,70,5)

                                # Optional Variables - need to hunt for start point

                                # CLOUDs
                                text_ident='GF1'
                                exists=cleanline.find(text_ident)
                                if exists!=-1:
                                    try:
                                        if RepresentsInt(cleanline[exists+3]):

                                            if Extract:
                                                ExtractionProcess(total_cloud_cover, total_cloud_flags,time_loc,INTMDI,'99',cleanline,
                                                                  exists+3,2,flagoffset=2)

                                            if Extract:
                                                ExtractionProcess(low_cloud_cover, low_cloud_flags,time_loc,INTMDI,'99',cleanline,
                                                                  exists+8,2)
                                    except IndexError:
                                        # string following data marker doesn't exist
                                        if dubious_flagged==0:
                                            dubious_flagged=WriteDubious(dubiousfile,rfile,text_ident, station, string_obs_time)
                                        

                                elif dubious_flagged==0:
                                    dubious_flagged=WriteDubious(dubiousfile,rfile,text_ident, station, string_obs_time)

                                text_ident='GA1'
                                exists_overall=cleanline.find(text_ident)

                                if exists_overall!=-1:
                                    cloud_amts=np.array([INTMDI for i in range(6)])
                                    cloud_hghts=np.array([INTMDI for i in range(6)])
                                    cloud_flags=np.array([INTMDI for i in range(6)])
                                    flagvals=['GA1','GA2','GA3','GA4','GA5','GA6']
                                    for cl,flg in enumerate(flagvals):
                                        exists=cleanline.find(flg)
                                        if exists!=-1:
                                            try:
                                                if RepresentsInt(cleanline[exists+3]):
                                                    if Extract:

                                                        ExtractionProcess(cloud_amts, cloud_flags,cl,INTMDI,'99',cleanline,
                                                                          exists+3,2)

                                                        ExtractionProcess(cloud_hghts, dummyflag,cl,INTMDI,'+99999',cleanline,
                                                                          exists+6,6,doflag=False)

                                                        # remove hard coded values?
                                                        if cloud_hghts[cl]!=INTMDI:
                                                            if cloud_hghts[cl]<=2000:
                                                                cloud_hghts[cl]=1
                                                            elif cloud_hghts[cl]>=6000:
                                                                cloud_hghts[cl]=3
                                                            elif cloud_hghts[cl]>=4:
                                                                cloud_hghts[cl]=2
                                            except IndexError:
                                                # string following data marker doesn't exist
                                                if dubious_flagged==0:
                                                    dubious_flagged=WriteDubious(dubiousfile,rfile,text_ident+'-'+flg, station, string_obs_time)
                                               
                                        elif dubious_flagged==0:
                                            dubious_flagged=WriteDubious(dubiousfile,rfile,text_ident+'-'+flg, station, string_obs_time)

                                    # end for loop
                                    SortClouds(total_cloud_cover, total_cloud_flags, time_loc, cloud_amts, cloud_flags, range(len(cloud_amts))) # select all using this slice

                                    lowclouds=np.where(np.array(cloud_hghts) == 1)[0]
                                    SortClouds(low_cloud_cover, low_cloud_flags, time_loc, cloud_amts, cloud_flags, lowclouds)
                                    medclouds=np.where(np.array(cloud_hghts) == 2)[0]
                                    SortClouds(mid_cloud_cover, mid_cloud_flags, time_loc, cloud_amts, cloud_flags, medclouds)
                                    hiclouds=np.where(np.array(cloud_hghts) == 3)[0]
                                    SortClouds(high_cloud_cover, high_cloud_flags, time_loc, cloud_amts, cloud_flags, hiclouds)                    


                                text_ident='GD1'
                                exists_overall=cleanline.find(text_ident)

                                if exists_overall!=-1:
                                    if (total_cloud_cover[time_loc] == INTMDI):
                                        cloud_amts=np.array([INTMDI for i in range(6)])
                                        cloud_amts2=np.array([INTMDI for i in range(6)])
                                        cloud_hghts=np.array([INTMDI for i in range(6)])
                                        cloud_flags=np.array([INTMDI for i in range(6)])
                                        flagvals=['GD1','GD2','GD3','GD4','GD5','GD6']
                                        for cl,flg in enumerate(flagvals):
                                            exists=cleanline.find(flg)
                                            if exists!=-1:
                                                try:
                                                    if RepresentsInt(cleanline[exists+3]):
                                                        if Extract:
                                                        # if TestToExtract(cloud_amts[cl],INTMDI,overwrite):


                                                            ExtractionProcess(cloud_amts, cloud_flags,cl,INTMDI,'9',cleanline,
                                                                              exists+3,1,flagoffset=1)

                                                            if cloud_amts[cl] >= 5 :
                                                                cloud_amts[cl]=INTMDI


                                                            ExtractionProcess(cloud_amts2, dummyflag,cl,INTMDI,'99',cleanline,
                                                                          exists+4,2,doflag=False)

                                                            ExtractionProcess(cloud_hghts, dummyflag,cl,INTMDI,'+99999',cleanline,
                                                                          exists+7,6,doflag=False)

                                                            # remove hard coded values?
                                                            if cloud_hghts[cl]!=INTMDI:
                                                                if cloud_hghts[cl]<=2000:
                                                                    cloud_hghts[cl]=1
                                                                elif cloud_hghts[cl]>=6000:
                                                                    cloud_hghts[cl]=3
                                                                elif cloud_hghts[cl]>=4:
                                                                    cloud_hghts[cl]=2
                                                except IndexError:
                                                    # string following data marker doesn't exist
                                                    if dubious_flagged==0:
                                                        dubious_flagged=WriteDubious(dubiousfile,rfile,text_ident+'-'+flg, station, string_obs_time)
                                                    
                                            elif dubious_flagged==0:
                                                dubious_flagged=WriteDubious(dubiousfile,rfile,text_ident+'-'+flg, station, string_obs_time)
                                      # end for loop


                                        SortClouds2(total_cloud_cover, total_cloud_flags, time_loc, cloud_amts, cloud_amts2, cloud_flags, range(len(cloud_amts))) # select whole list with slice

                                        lowclouds=np.where(np.array(cloud_hghts) == 1)[0]
                                        if len(lowclouds)>=1: 
                                            SortClouds2(low_cloud_cover, low_cloud_flags, time_loc, cloud_amts, cloud_amts2, cloud_flags, lowclouds)

                                        medclouds=np.where(np.array(cloud_hghts) == 2)[0]
                                        if len(medclouds)>=1:  
                                            SortClouds2(mid_cloud_cover, mid_cloud_flags, time_loc, cloud_amts, cloud_amts2, cloud_flags, medclouds)

                                        hiclouds=np.where(np.array(cloud_hghts) == 3)[0]
                                        if len(hiclouds)>=1:  
                                            SortClouds2(high_cloud_cover, high_cloud_flags, time_loc, cloud_amts, cloud_amts2, cloud_flags, hiclouds)


                                # PAST-SIGWX
                                text_ident='AY1'
                                exists=cleanline.find(text_ident)

                                if exists!=-1:
                                    try:
                                        if RepresentsInt(cleanline[exists+3]):
                                            if Extract:

                                                ExtractionProcess(past_sigwx1, past_sigwx1_flag,time_loc,INTMDI,'-',cleanline,
                                                                  exists+3,1)

                                                ExtractionProcess(past_sigwx1_period, dummyflag,time_loc,INTMDI,'99',cleanline,
                                                                  exists+5,2,doflag=False)

                                    except IndexError:
                                        if dubious_flagged==0:
                                            dubious_flagged=WriteDubious(dubiousfile,rfile,text_ident, station, string_obs_time)

                                elif dubious_flagged==0:
                                    dubious_flagged=WriteDubious(dubiousfile,rfile,text_ident, station, string_obs_time)

                                text_ident='AZ1'
                                exists=cleanline.find(text_ident)

                                if exists!=-1:
                                    try:
                                        if RepresentsInt(cleanline[exists+3]):
                                            if Extract:

                                                ExtractionProcess(past_sigwx1, past_sigwx1_flag,time_loc,INTMDI,'-',cleanline,
                                                                  exists+3,1)

                                                ExtractionProcess(past_sigwx1_period, dummyflag,time_loc,INTMDI,'99',cleanline,
                                                                  exists+5,2,doflag=False)

                                    except IndexError:
                                        if dubious_flagged==0:
                                            dubious_flagged=WriteDubious(dubiousfile,rfile,text_ident, station, string_obs_time)

                                elif dubious_flagged==0:
                                    dubious_flagged=WriteDubious(dubiousfile,rfile,text_ident, station, string_obs_time)

                                # PRECIP

                                text_ident='AA1'
                                exists=cleanline.find(text_ident)

                                if exists!=-1:
                                    try:
                                        if RepresentsInt(cleanline[exists+3]):
                                            if Extract:

                                                ExtractionProcess(precip1_period, precip1_flag,time_loc,INTMDI,'99',cleanline,
                                                                  exists+3,2, flagoffset=5, doflag=True)

                                                if precip1_period[time_loc] < 0:
                                                    precip1_period[time_loc]=INTMDI

                                                ExtractionProcess(precip1_depth, dummyflag,time_loc,FLTMDI,'9999',cleanline,
                                                                  exists+5,4,doflag=False,divisor=10.)


                                                # these two pass in empty strings as missing data test
                                                ExtractionProcess(precip1_condition, dummyflag,time_loc,'','9',cleanline,
                                                                  exists+9,1,doflag=False)

                                    except IndexError:
                                        if dubious_flagged==0:
                                            dubious_flagged=WriteDubious(dubiousfile,rfile,text_ident, station, string_obs_time)


                                elif dubious_flagged==0:
                                    dubious_flagged=WriteDubious(dubiousfile,rfile,text_ident, station, string_obs_time)


                                # SUN DURATION

                                text_ident='GJ1'
                                exists=cleanline.find(text_ident)

                                if exists!=-1:
                                    try:
                                        if RepresentsInt(cleanline[exists+3]):
                                            if Extract:

                                                ExtractionProcess(sun_duration, sun_durationqc,time_loc,INTMDI,'9999',cleanline,
                                                                  exists+3,4)
                                    except IndexError:
                                        if dubious_flagged==0:
                                            dubious_flagged=WriteDubious(dubiousfile,rfile,text_ident, station, string_obs_time)


                                elif dubious_flagged==0:
                                    dubious_flagged=WriteDubious(dubiousfile,rfile,text_ident, station, string_obs_time)


                                # WIND GUST

                                for text_ident in ['OA1','OA2','OA3','OA4']:
                                    # test all of the possible locations for wind gust
                                    exists=cleanline.find(text_ident)

                                    if exists!=-1:
                                        try:
                                            if RepresentsInt(cleanline[exists+3]):
                                                if cleanline[exists+3] == "4":
                                                    # then this is the maximum gust speed - see ish-format-document

                                                    if Extract:

                                                        ExtractionProcess(wind_gust_period, dummyflag,time_loc,INTMDI,'99',cleanline,
                                                                          exists+4,2,doflag=False)

                                                        if wind_gust_period[time_loc] < 0:
                                                            wind_gust_period[time_loc]=FLTMDI

                                                        ExtractionProcess(wind_gust_value, wind_gust_flags,time_loc,FLTMDI,'9999',cleanline,
                                                                          exists+6,4,divisor=10.)
                                        except IndexError:
                                            if dubious_flagged==0:
                                                dubious_flagged=WriteDubious(dubiousfile,rfile,text_ident, station, string_obs_time)

                                    elif dubious_flagged==0:
                                        dubious_flagged=WriteDubious(dubiousfile,rfile,text_ident, station, string_obs_time)


                                if Extra:
                                    # PRESENT SIGWX
                                    text_ident='AW1'
                                    exists=cleanline.find(text_ident)

                                    if exists!=-1:
                                        try:
                                            if RepresentsInt(cleanline[exists+3]):
                                                if Extract:

                                                    ExtractionProcess(present_sigwx, present_sigwx_flags,time_loc,INTMDI,'--',cleanline,
                                                                  exists+3,2)

                                        except IndexError:
                                            if dubious_flagged==0:
                                                dubious_flagged=WriteDubious(dubiousfile,rfile,text_ident, station, string_obs_time)

                                    elif dubious_flagged==0:
                                        dubious_flagged=WriteDubious(dubiousfile,rfile,text_ident, station, string_obs_time)

                                    # PAST-SIGWX2

                                    for text_ident in ['AY1','AY2','AZ1','AZ2']:

                                        exists=cleanline.find(text_ident)

                                        if exists!=-1:
                                            try:
                                                if RepresentsInt(cleanline[exists+3]):
                                                    if Extract:

                                                        ExtractionProcess(past_sigwx2, past_sigwx2_flag,time_loc,INTMDI,'-',cleanline,
                                                                      exists+3,1)
                                                        ExtractionProcess(past_sigwx2_period,dummyflag,time_loc,INTMDI,'99',cleanline,
                                                                      exists+5,2,doflag=False)


                                                        value=ExtractValues(INTMDI, cleanline,exists+5,2,'99',doflag=False)
                                                        if value!=INTMDI:
                                                            past_sigwx2_period[time_loc]=value

                                            except IndexError:
                                                if dubious_flagged==0:
                                                    dubious_flagged=WriteDubious(dubiousfile,rfile,text_ident, station, string_obs_time)


                                        elif dubious_flagged==0:
                                            dubious_flagged=WriteDubious(dubiousfile,rfile,text_ident, station, string_obs_time)


                                    # PRECIP

                                    text_ident='AA2'
                                    exists=cleanline.find(text_ident)

                                    if exists!=-1:
                                        try:
                                            if RepresentsInt(cleanline[exists+3]):
                                                if Extract:

                                                    ExtractionProcess(precip2_period, dummyflag,time_loc,INTMDI,'99',cleanline,
                                                                      exists+3,2,doflag=False)

                                                    if precip2_period[time_loc] < 0:
                                                        precip2_period[time_loc]=INTMDI

                                                    ExtractionProcess(precip2_depth,dummyflag,time_loc,FLTMDI,'9999',cleanline,
                                                                      exists+5,4,doflag=False,divisor=10.)

                                                    # leave this as is because of separate value and flag tests
                                                    value,flag=ExtractValues('', cleanline,exists+9,1,'9')
                                                    if value!='':
                                                        precip2_condition[time_loc]=value
                                                    if flag in range(20):
                                                        precip2_flag[time_loc]=flag

                                        except IndexError:
                                            if dubious_flagged==0:
                                                dubious_flagged=WriteDubious(dubiousfile,rfile,text_ident, station, string_obs_time)

                                    elif dubious_flagged==0:
                                        dubious_flagged=WriteDubious(dubiousfile,rfile,text_ident, station, string_obs_time)

                                    text_ident='AA3'
                                    exists=cleanline.find(text_ident)

                                    if exists!=-1:
                                        try:
                                            if RepresentsInt(cleanline[exists+3]):
                                                if Extract:

                                                    ExtractionProcess(precip3_period, dummyflag,time_loc,INTMDI,'99',cleanline,
                                                                      exists+3,2,doflag=False)

                                                    if precip3_period[time_loc] < 0:
                                                        precip3_period[time_loc]=INTMDI

                                                    ExtractionProcess(precip3_depth,dummyflag,time_loc,FLTMDI,'9999',cleanline,
                                                                      exists+5,4,doflag=False,divisor=10.)

                                                    # leave this as is because of separate value and flag tests
                                                    value,flag=ExtractValues('', cleanline,exists+9,1,'9')
                                                    if value!='':
                                                        precip3_condition[time_loc]=value
                                                    if flag in range(20):
                                                        precip3_flag[time_loc]=flag

                                        except IndexError:
                                            if dubious_flagged==0:
                                                dubious_flagged=WriteDubious(dubiousfile,rfile,text_ident, station, string_obs_time)


                                    elif dubious_flagged==0:
                                        dubious_flagged=WriteDubious(dubiousfile,rfile,text_ident, station, string_obs_time)

                                    text_ident='AA3'
                                    exists=cleanline.find(text_ident)

                                    if exists!=-1:
                                        try:
                                            if RepresentsInt(cleanline[exists+3]):
                                                if Extract:

                                                    ExtractionProcess(precip4_period, dummyflag,time_loc,INTMDI,'99',cleanline,
                                                                      exists+3,2,doflag=False)

                                                    if precip4_period[time_loc] < 0:
                                                        precip4_period[time_loc]=INTMDI

                                                    ExtractionProcess(precip4_depth,dummyflag,time_loc,FLTMDI,'9999',cleanline,
                                                                      exists+5,4,doflag=False,divisor=10.)

                                                    # leave this as is because of separate value and flag tests
                                                    value,flag=ExtractValues('', cleanline,exists+9,1,'9')
                                                    if value!='':
                                                        precip4_condition[time_loc]=value
                                                    if flag in range(20):
                                                        precip4_flag[time_loc]=flag

                                        except IndexError:
                                            if dubious_flagged==0:
                                                dubious_flagged=WriteDubious(dubiousfile,rfile,text_ident, station, string_obs_time)


                                    elif dubious_flagged==0:
                                        dubious_flagged=WriteDubious(dubiousfile,rfile,text_ident, station, string_obs_time)

                                    # EXTREME TEMPERATURES

                                    # RJHD - 11 March 2014 - these could be converted to Tmax and Tmin using the code information

                                    for text_ident in ['KA1','KA2']:
                                        exists=cleanline.find(text_ident)

                                        if exists!=-1:
                                            try:
                                                if RepresentsInt(cleanline[exists+3]):
                                                    if cleanline[exists+6] == "N":
                                                        # then this is the minimum temperature - see ish-format-document

                                                        if Extract:

                                                            ExtractionProcess(minimum_temp_period, dummyflag,time_loc,INTMDI,'999',cleanline,
                                                                              exists+3,3,doflag=False,divisor=10.)

                                                            if minimum_temp_period[time_loc] < 0:
                                                                minimum_temp_period[time_loc]=FLTMDI

                                                            ExtractionProcess(minimum_temp_value, minimum_temp_flags,time_loc,FLTMDI,'+9999',cleanline,
                                                                              exists+7,5,divisor=10.)

                                                    elif cleanline[exists+6] == "M":
                                                        # then this is the minimum temperature - see ish-format-document

                                                        if Extract:
 
                                                            ExtractionProcess(maximum_temp_period, dummyflag,time_loc,INTMDI,'999',cleanline,
                                                                              exists+3,3,doflag=False,divisor=10.)

                                                            if maximum_temp_period[time_loc] < 0:
                                                                maximum_temp_period[time_loc]=FLTMDI

                                                            ExtractionProcess(maximum_temp_value, maximum_temp_flags,time_loc,FLTMDI,'+9999',cleanline,
                                                                              exists+7,5,divisor=10.)

                                            except IndexError:
                                                if dubious_flagged==0:
                                                    dubious_flagged=WriteDubious(dubiousfile,rfile,text_ident, station, string_obs_time)

                                        elif dubious_flagged==0:
                                            dubious_flagged=WriteDubious(dubiousfile,rfile,text_ident, station, string_obs_time)


                                # end Extra variables    
                            # end if time_loc != HoursBetween
                        # end line in file loop
                except IOError:
                     print "Cannot find file: ", rfile

            # end file loop

        print "parsed all files for station ", station



        mask_vals=np.where(np.array(input_station_id) != 'null')[0]

        if hours:
            times_out=TimeStamps[mask_vals]
        else:
            times_out=TimeStamps[mask_vals]/24.

        # apply the mask
        input_station_id=np.array(input_station_id)[mask_vals]
        temperatures=temperatures[mask_vals]
        temperature_flags=temperature_flags[mask_vals]
        dewpoints=dewpoints[mask_vals]
        dewpoint_flags=dewpoint_flags[mask_vals]
        total_cloud_cover=total_cloud_cover[mask_vals]
        total_cloud_flags=total_cloud_flags[mask_vals]
        low_cloud_cover=low_cloud_cover[mask_vals]
        low_cloud_flags=low_cloud_flags[mask_vals]
        mid_cloud_cover=mid_cloud_cover[mask_vals]
        mid_cloud_flags=mid_cloud_flags[mask_vals]
        high_cloud_cover=high_cloud_cover[mask_vals]
        high_cloud_flags=high_cloud_flags[mask_vals]
        cloud_base=cloud_base[mask_vals]
        cloud_base_flags=cloud_base_flags[mask_vals]
        windspeeds=windspeeds[mask_vals]
        windspeeds_flags=windspeeds_flags[mask_vals]  
        winddirs=winddirs[mask_vals]
        winddirs_flags=winddirs_flags[mask_vals]
        past_sigwx1=past_sigwx1[mask_vals]
        past_sigwx1_period=past_sigwx1_period[mask_vals]
        past_sigwx1_flag=past_sigwx1_flag[mask_vals]
        precip1_period=precip1_period[mask_vals]
        precip1_depth=precip1_depth[mask_vals]
        precip1_condition=np.array(precip1_condition)[mask_vals]
        precip1_flag=precip1_flag[mask_vals]
        slp=slp[mask_vals]
        slp_flag=slp_flag[mask_vals]
        sun_duration=sun_duration[mask_vals]
        sun_durationqc=sun_durationqc[mask_vals]
        wind_gust_period=wind_gust_period[mask_vals]
        wind_gust_value=wind_gust_value[mask_vals]
        wind_gust_flags=wind_gust_flags[mask_vals]

        if Extra:
            windtypes=np.array(windtypes)[mask_vals]
            present_sigwx=present_sigwx[mask_vals]
            present_sigwx_flags=present_sigwx_flags[mask_vals]
            past_sigwx2=past_sigwx2[mask_vals]
            past_sigwx2_period=past_sigwx2_period[mask_vals]
            past_sigwx2_flag=past_sigwx2_flag[mask_vals]
            precip2_period=precip2_period[mask_vals]
            precip2_depth=precip2_depth[mask_vals]
            precip2_condition=np.array(precip2_condition)[mask_vals]
            precip2_flag=precip2_flag[mask_vals]
            precip3_period=precip3_period[mask_vals]
            precip3_depth=precip3_depth[mask_vals]
            precip3_condition=np.array(precip3_condition)[mask_vals]
            precip3_flag=precip3_flag[mask_vals]
            precip4_period=precip4_period[mask_vals]
            precip4_depth=precip4_depth[mask_vals]
            precip4_condition=np.array(precip4_condition)[mask_vals]
            precip4_flag=precip4_flag[mask_vals]
            maximum_temp_period=maximum_temp_period[mask_vals]
            maximum_temp_value=maximum_temp_value[mask_vals]
            maximum_temp_flags=maximum_temp_flags[mask_vals]
            minimum_temp_period=minimum_temp_period[mask_vals]
            minimum_temp_value=minimum_temp_value[mask_vals]
            minimum_temp_flags=minimum_temp_flags[mask_vals]

        netcdf_filename=OUTPUT_DATA_DIR+'/'+station+'.nc'

        if do_zip:
            netcdf_outfile = ncdf.Dataset(netcdf_filename,'w', format='NETCDF4')
        else:
            netcdf_outfile = ncdf.Dataset(netcdf_filename,'w', format='NETCDF3_CLASSIC')

        time=netcdf_outfile.createDimension('time',len(times_out))
        char_len=netcdf_outfile.createDimension('character_length',4)
        long_char_len=netcdf_outfile.createDimension('long_character_length',12)

        # create variables
        timesvar=netcdf_outfile.createVariable('time','f8',('time',), zlib = do_zip)
        stationsvar=netcdf_outfile.createVariable('input_station_id','S1',('time','long_character_length',), zlib = do_zip)
        tempsvar=netcdf_outfile.createVariable('temperatures','f8',('time',), zlib = do_zip)
        tempsflagsvar=netcdf_outfile.createVariable('temperature_flags','i4',('time',), zlib = do_zip)
        dewsvar=netcdf_outfile.createVariable('dewpoints','f8',('time',), zlib = do_zip)
        dewsflagsvar=netcdf_outfile.createVariable('dewpoint_flags','i4',('time',), zlib = do_zip)
        tcvar=netcdf_outfile.createVariable('total_cloud_cover','i4',('time',), zlib = do_zip)
        tcfvar=netcdf_outfile.createVariable('total_cloud_flags','i4',('time',), zlib = do_zip)
        lcvar=netcdf_outfile.createVariable('low_cloud_cover','i4',('time',), zlib = do_zip)
        lcfvar=netcdf_outfile.createVariable('low_cloud_flags','i4',('time',), zlib = do_zip)
        mcvar=netcdf_outfile.createVariable('mid_cloud_cover','i4',('time',), zlib = do_zip)
        mcfvar=netcdf_outfile.createVariable('mid_cloud_flags','i4',('time',), zlib = do_zip)
        hcvar=netcdf_outfile.createVariable('high_cloud_cover','i4',('time',), zlib = do_zip)
        hcfvar=netcdf_outfile.createVariable('high_cloud_flags','i4',('time',), zlib = do_zip)
        cbvar=netcdf_outfile.createVariable('cloud_base','f8',('time',), zlib = do_zip)
        cbfvar=netcdf_outfile.createVariable('cloud_base_flags','i4',('time',), zlib = do_zip)
        wsvar=netcdf_outfile.createVariable('windspeeds','f8',('time',), zlib = do_zip)
        wsfvar=netcdf_outfile.createVariable('windspeeds_flags','i4',('time',), zlib = do_zip)
        wdvar=netcdf_outfile.createVariable('winddirs','i4',('time',), zlib = do_zip)
        wdfvar=netcdf_outfile.createVariable('winddirs_flags','i4',('time',), zlib = do_zip)
        pswx1var=netcdf_outfile.createVariable('past_sigwx1','i4',('time',), zlib = do_zip)
        pswx1pvar=netcdf_outfile.createVariable('past_sigwx1_period','i4',('time',), zlib = do_zip)
        pswx1fvar=netcdf_outfile.createVariable('past_sigwx1_flag','i4',('time',), zlib = do_zip)
        ppt1pvar=netcdf_outfile.createVariable('precip1_period','i4',('time',), zlib = do_zip)
        ppt1dvar=netcdf_outfile.createVariable('precip1_depth','f8',('time',), zlib = do_zip)
        ppt1cvar=netcdf_outfile.createVariable('precip1_condition','S1',('time','character_length',), zlib = do_zip)
        ppt1fvar=netcdf_outfile.createVariable('precip1_flag','i4',('time',), zlib = do_zip)
        slpvar=netcdf_outfile.createVariable('slp','f8',('time',), zlib = do_zip)
        slpfvar=netcdf_outfile.createVariable('slp_flag','i4',('time',), zlib = do_zip)
        sdvar=netcdf_outfile.createVariable('sun_duration','f8',('time',), zlib = do_zip)
        sdfvar=netcdf_outfile.createVariable('sun_durationqc','i4',('time',), zlib = do_zip)
        wgstpvar=netcdf_outfile.createVariable('wind_gust_period','f8',('time',), zlib = do_zip)
        wgstvvar=netcdf_outfile.createVariable('wind_gust','f8',('time',), zlib = do_zip)
        wgstfvar=netcdf_outfile.createVariable('wind_gust_flag','i4',('time',), zlib = do_zip)

        if Extra:
            pswx2var=netcdf_outfile.createVariable('past_sigwx2','i4',('time',), zlib = do_zip)
            pswx2pvar=netcdf_outfile.createVariable('past_sigwx2_period','i4',('time',), zlib = do_zip)
            pswx2fvar=netcdf_outfile.createVariable('past_sigwx2_flag','i4',('time',), zlib = do_zip)
            wtvar=netcdf_outfile.createVariable('windtypes','S1',('time','character_length'), zlib = do_zip)
            swxvar=netcdf_outfile.createVariable('present_sigwx','i4',('time',), zlib = do_zip)
            swxfvar=netcdf_outfile.createVariable('present_sigwx_flags','i4',('time',), zlib = do_zip)
            ppt2pvar=netcdf_outfile.createVariable('precip2_period','i4',('time',), zlib = do_zip)
            ppt2dvar=netcdf_outfile.createVariable('precip2_depth','f8',('time',), zlib = do_zip)
            ppt2cvar=netcdf_outfile.createVariable('precip2_condition','S1',('time','character_length',), zlib = do_zip)
            ppt2fvar=netcdf_outfile.createVariable('precip2_flag','i4',('time',), zlib = do_zip)
            ppt3pvar=netcdf_outfile.createVariable('precip3_period','i4',('time',), zlib = do_zip)
            ppt3dvar=netcdf_outfile.createVariable('precip3_depth','f8',('time',), zlib = do_zip)
            ppt3cvar=netcdf_outfile.createVariable('precip3_condition','S1',('time','character_length',), zlib = do_zip)
            ppt3fvar=netcdf_outfile.createVariable('precip3_flag','i4',('time',), zlib = do_zip)
            ppt4pvar=netcdf_outfile.createVariable('precip4_period','i4',('time',), zlib = do_zip)
            ppt4dvar=netcdf_outfile.createVariable('precip4_depth','f8',('time',), zlib = do_zip)
            ppt4cvar=netcdf_outfile.createVariable('precip4_condition','S1',('time','character_length',), zlib = do_zip)
            ppt4fvar=netcdf_outfile.createVariable('precip4_flag','i4',('time',), zlib = do_zip)
            maxtpvar=netcdf_outfile.createVariable('maximum_temp_period','f8',('time',), zlib = do_zip)
            maxtvvar=netcdf_outfile.createVariable('maximum_temp_value','f8',('time',), zlib = do_zip)
            maxtfvar=netcdf_outfile.createVariable('maximum_temp_flag','i4',('time',), zlib = do_zip)
            mintpvar=netcdf_outfile.createVariable('minimum_temp_period','f8',('time',), zlib = do_zip)
            mintvvar=netcdf_outfile.createVariable('minimum_temp_value','f8',('time',), zlib = do_zip)
            mintfvar=netcdf_outfile.createVariable('minimum_temp_flag','i4',('time',), zlib = do_zip)

        # variables attributes
        print "Writing Attributes"

        timesvar.long_name='time'
        timesvar.standard_name='time'
        if hours:
            timesvar.units='hours since {}'.format(dt.datetime.strftime(dt.datetime(STARTYEAR,1,1,0,0), "%Y-%m-%d %H:%M"))
        else:
            timesvar.units='days since {}'.format(dt.datetime.strftime(dt.datetime(STARTYEAR,1,1,0,0), "%Y-%m-%d %H:%M"))
        timesvar.axis='T'
        timesvar.calendar='gregorian'
        timesvar.valid_min=0.
        stationsvar.long_name='Primary source for timestep (may be multiple sources for composite stations)'
        stationsvar.units='USAF - WBAN from ISD source'
        stationsvar.axis='T'
        stationsvar.missing_value='null'

        try:
            tmin,tmax=np.min(temperatures[np.where(temperatures != FLTMDI)[0]]),np.max(temperatures[np.where(temperatures != FLTMDI)[0]])
        except ValueError:
            tmin,tmax=FLTMDI,FLTMDI
        WriteAttributes(tempsvar,'Dry bulb temperatures nearest to reporting hour',FLTMDI,'Degrees C','T',tmin,tmax,standard_name = 'surface_temperature')
        WriteFlagAttributes(tempsflagsvar,'ISD flags for temperature - see ISD documentation',INTMDI,'T')

        try:
            dmin,dmax=np.min(dewpoints[np.where(dewpoints != FLTMDI)[0]]),np.max(dewpoints[np.where(dewpoints != FLTMDI)[0]])
        except ValueError:
            dmin,dmax=FLTMDI,FLTMDI
        WriteAttributes(dewsvar,'Dew point temperatures nearest to reporting hour',FLTMDI,'Degrees C','T',dmin,dmax, standard_name =  'dew_point_temperature')
        WriteFlagAttributes(dewsflagsvar,'ISD flags for dewpoint temperature - see ISD documentation',INTMDI,'T')

        WriteAttributes(tcvar,'Total cloud cover (oktas) derived in priority order GA, GF, GD - see ISD documentation', INTMDI, '1', 'T', 0,8, standard_name = "cloud_area_fraction")
        WriteFlagAttributes(tcfvar,'ISD flags for total cloud - see ISD documentation',INTMDI,'T')

        WriteAttributes(lcvar,'Low cloud cover (oktas) derived in priority order GA, GF, GD - see ISD documentation', INTMDI, '1', 'T', 0,8, standard_name = "low_type_cloud_area_fraction")
        WriteFlagAttributes(lcfvar,'ISD flags for low cloud - see ISD documentation',INTMDI,'T')

        WriteAttributes(mcvar,'Mid cloud cover (oktas) derived in priority order GA, GF, GD - see ISD documentation', INTMDI, '1', 'T', 0,8, standard_name = "medium_type_cloud_area_fraction")
        WriteFlagAttributes(mcfvar,'ISD flags for mid cloud - see ISD documentation',INTMDI,'T')

        WriteAttributes(hcvar,'High cloud cover (oktas) derived in priority order GA, GF, GD - see ISD documentation', INTMDI, '1', 'T', 0,8, standard_name = "high_type_cloud_area_fraction")
        WriteFlagAttributes(hcfvar,'ISD flags for high cloud - see ISD documentation',INTMDI,'T')

        try:
            cbmin,cbmax=np.min(cloud_base[np.where(cloud_base != INTMDI)[0]]),np.max(cloud_base[np.where(cloud_base != INTMDI)[0]])
        except ValueError:
            cbmin,cbmax=FLTMDI,FLTMDI
        WriteAttributes(cbvar,'Cloud base of lowest cloud layer nearest to reporting hour', INTMDI, 'meters', 'T', cbmin, cbmax, standard_name = 'cloud_base_altitude')
        WriteFlagAttributes(cbfvar,'ISD flags for cloud base - see ISD documentation',INTMDI,'T')

        try:
            wsmin,wsmax=np.min(windspeeds[np.where(windspeeds != FLTMDI)[0]]),np.max(windspeeds[np.where(windspeeds != FLTMDI)[0]])
        except ValueError:
            wsmin,wsmax=FLTMDI,FLTMDI
        WriteAttributes(wsvar,'Wind speed', FLTMDI, 'meters per second', 'T', wsmin, wsmax,standard_name = 'wind_speed')
        WriteFlagAttributes(wsfvar,'ISD flags for windspeed - see ISD documentation',INTMDI,'T')

        WriteAttributes(wdvar,'Wind Direction', INTMDI, 'Degrees', 'T', 0, 360, standard_name = 'wind_from_direction')
        WriteFlagAttributes(wdfvar,'ISD flags for wind direction - see ISD documentation',INTMDI,'T')

        WriteAttributes(pswx1var,'Station reports of past significant weather phenomena', INTMDI, '1', 'T', 0, 9)
        WriteFlagAttributes(pswx1fvar,'ISD flags for reported past significant weather - see ISD documentation',INTMDI,'T')
        WriteAttributes(pswx1pvar,'Period of significant weather report', INTMDI, 'Hours', 'T', 0, 24)

        WriteAttributes(ppt1pvar,'Period of Precipitation Report', INTMDI, 'Hours', 'T', 0, 98)
        WriteAttributes(ppt1dvar,'Depth of Precipitation Reported', FLTMDI, 'mm', 'T', 0, 999.8, standard_name = 'lwe_thickness_of_precipitation_amount')
        WriteFlagAttributes(ppt1cvar,'Precipitation Code (denotes if trace amount)', 'null','T')
        WriteFlagAttributes(ppt1fvar,'ISD flags for first precip field - see ISD documentation', INTMDI,'T')
        ppt1cvar.units='no_unit'

        try:
            smin,smax=np.min(slp[np.where(slp != FLTMDI)[0]]),np.max(slp[np.where(slp != FLTMDI)[0]])
        except ValueError:
            smin,smax=FLTMDI,FLTMDI
        WriteAttributes(slpvar,'Reported Sea Level Pressure',FLTMDI, 'hPa', 'T', smin, smax, standard_name = 'air_pressure_at_sea_level')
        WriteFlagAttributes(slpfvar,'ISD flags for slp field - see ISD documentation',INTMDI,'T')

        WriteAttributes(sdvar,'Reported Sunshine Duration', INTMDI, 'minutes', 'T', 0, 6000, standard_name = 'duration_of_sunshine')
        WriteFlagAttributes(sdfvar,'ISD flags sun duration field - see ISD documentation',INTMDI,'T')

        WriteAttributes(wgstpvar,'Period of Maximum Wind Gust Speed', INTMDI, 'Hours', 'T', 0, 48)
        WriteAttributes(wgstvvar,'Maximum Wind Gust Speed Reported',FLTMDI, 'meters per second', 'T', 0, 200.0, standard_name = 'wind_speed_of_gust')
        WriteFlagAttributes(wgstfvar,'ISD flags for wind gust field - see ISD documentation', INTMDI,'T')
        if Extra:
            WriteFlagAttributes(wtvar,'Wind observation type - see ISD documentation','null','T')

            WriteAttributes(pswx2var,'Station reports of past significant weather phenomena (2)', INTMDI, 'no_unit', 'T', 0, 9)
            WriteFlagAttributes(pswx2fvar,'ISD flags for reported past significant weather - see ISD documentation',INTMDI,'T')
            WriteAttributes(pswx2pvar,'Period of significant weather report', INTMDI, 'Hours', 'T', 0, 24)

            WriteAttributes(swxvar,'Station reports of present significant weather phenomena', INTMDI, 'no_unit', 'T', 0, 99)
            WriteFlagAttributes(swxfvar,'ISD flags for reported present significant weather - see ISD documentation',INTMDI,'T')

            WriteAttributes(ppt2pvar,'Period of Precipitation Report', INTMDI, 'Hours', 'T', 0, 98)
            WriteAttributes(ppt2dvar,'Depth of Precipitation Reported', FLTMDI, 'mm', 'T', 0, 999.8)
            WriteFlagAttributes(ppt2cvar,'Denotes if trace amount', 'null','T')
            WriteFlagAttributes(ppt2fvar,'ISD flags for second precip field - see ISD documentation', INTMDI,'T')
            ppt2cvar.units='no_unit'

            WriteAttributes(ppt3pvar,'Period of Precipitation Report', INTMDI, 'Hours', 'T', 0, 98)
            WriteAttributes(ppt3dvar,'Depth of Precipitation Reported', FLTMDI, 'mm', 'T', 0, 999.8)
            WriteFlagAttributes(ppt3cvar,'Denotes if trace amount', 'null','T')
            WriteFlagAttributes(ppt3fvar,'ISD flags for third precip field - see ISD documentation', INTMDI,'T')
            ppt3cvar.units='no_unit'

            WriteAttributes(ppt4pvar,'Period of Precipitation Report', INTMDI, 'Hours', 'T', 0, 98)
            WriteAttributes(ppt4dvar,'Depth of Precipitation Reported', FLTMDI, 'mm', 'T', 0, 999.8)
            WriteFlagAttributes(ppt4cvar,'Denotes if trace amount', 'null','T')
            WriteFlagAttributes(ppt4fvar,'ISD flags for fourth precip field - see ISD documentation', INTMDI,'T')
            ppt4cvar.units='no_unit'

            try:
                xtmin,xtmax=np.min(maximum_temp_value[np.where(maximum_temp_value != FLTMDI)[0]]),np.max(maximum_temp_value[np.where(maximum_temp_value != FLTMDI)[0]])
            except ValueError:
                xtmin,xtmax=FLTMDI,FLTMDI
            WriteAttributes(maxtpvar,'Period of maximum temperature report', FLTMDI, 'Hours', 'T', 0, 48)
            WriteAttributes(maxtvvar,'Reported dry bulb maximum temperature', FLTMDI,'Degrees C','T',xtmin,xtmax)
            WriteFlagAttributes(maxtfvar,'ISD flags for maximum temperature field - see ISD documentation', INTMDI, 'T')


            try:
                ntmin,ntmax=np.min(minimum_temp_value[np.where(minimum_temp_value != FLTMDI)[0]]),np.max(minimum_temp_value[np.where(minimum_temp_value != FLTMDI)[0]])
            except ValueError:
                ntmin,ntmax=FLTMDI,FLTMDI
            WriteAttributes(mintpvar,'Period of minimum temperature report', FLTMDI, 'Hours', 'T', 0, 48)
            WriteAttributes(mintvvar,'Reported dry bulb minimum temperature', FLTMDI,'Degrees C','T',ntmin,ntmax)
            WriteFlagAttributes(mintfvar,'ISD flags for minimum temperature field - see ISD documentation', INTMDI, 'T')



        # global attributes
        netcdf_outfile.station_information='Where station is a composite the station id refers to the primary source used in the timestep and does apply to all elements'

        netcdf_outfile.file_created=dt.datetime.strftime(dt.datetime.now(), "%a %b %d, %H:%M %Y")
        netcdf_outfile.station_id=station
        netcdf_outfile.latitude=StationLat[st]
        netcdf_outfile.longitude=StationLon[st]
        netcdf_outfile.elevation=StationElv[st]
        netcdf_outfile.date_created = dt.datetime.strftime(dt.datetime.now(), "%Y-%m-%d, %H:%M")
        netcdf_outfile.history = "Created by mk_netcdf_files.py \n"
 
        print "Writing data to netcdf file"

        # write data into file
        # ncdf.stringtochar changes np array of strings
        #    to times x character_length array of single
        #    character strings


        timesvar[:]=times_out
        stationsvar[:]=ncdf.stringtochar(input_station_id)
        tempsvar[:]=temperatures
        tempsflagsvar[:]=temperature_flags
        dewsvar[:]=dewpoints
        dewsflagsvar[:]=dewpoint_flags
        tcvar[:]=total_cloud_cover
        tcfvar[:]=total_cloud_flags
        lcvar[:]=low_cloud_cover
        lcfvar[:]=low_cloud_flags
        mcvar[:]=mid_cloud_cover
        mcfvar[:]=mid_cloud_flags
        hcvar[:]=high_cloud_cover
        hcfvar[:]=high_cloud_flags
        cbvar[:]=cloud_base
        cbfvar[:]=cloud_base_flags
        wsvar[:]=windspeeds
        wsfvar[:]=windspeeds_flags
        wdvar[:]=winddirs
        wdfvar[:]=winddirs_flags
        pswx1var[:]=past_sigwx1
        pswx1pvar[:]=past_sigwx1_period
        pswx1fvar[:]=past_sigwx1_flag
        ppt1pvar[:]=precip1_period
        ppt1dvar[:]=precip1_depth
        ppt1cvar[:]=ncdf.stringtochar(precip1_condition)
        ppt1fvar[:]=precip1_flag
        slpvar[:]=slp
        slpfvar[:]=slp_flag
        sdvar[:]=sun_duration
        sdfvar[:]=sun_durationqc
        wgstpvar[:]=wind_gust_period
        wgstfvar[:]=wind_gust_flags
        wgstvvar[:]=wind_gust_value


        if Extra:
            pswx2var[:]=past_sigwx2
            pswx2pvar[:]=past_sigwx2_period
            pswx2fvar[:]=past_sigwx2_flag
            wtvar[:]=ncdf.stringtochar(windtypes)
            swxvar[:]=present_sigwx
            swxfvar[:]=present_sigwx_flags
            ppt2pvar[:]=precip2_period
            ppt2dvar[:]=precip2_depth
            ppt2cvar[:]=ncdf.stringtochar(precip2_condition)
            ppt2fvar[:]=precip2_flag
            ppt3pvar[:]=precip3_period
            ppt3dvar[:]=precip3_depth
            ppt3cvar[:]=ncdf.stringtochar(precip3_condition)
            ppt3fvar[:]=precip3_flag
            ppt4pvar[:]=precip4_period
            ppt4dvar[:]=precip4_depth
            ppt4cvar[:]=ncdf.stringtochar(precip4_condition)
            ppt4fvar[:]=precip4_flag
            maxtpvar[:]=maximum_temp_period
            maxtvvar[:]=maximum_temp_value
            maxtfvar[:]=maximum_temp_flags
            mintpvar[:]=minimum_temp_period
            mintvvar[:]=minimum_temp_value
            mintfvar[:]=minimum_temp_flags
        #  extra

        netcdf_outfile.close()

        # gzip file
        cmd='gzip -f -9 '+netcdf_filename
        print cmd
        subprocess.call(cmd,shell=True)

        print "Done station "+station

        print dt.datetime.now()-dbg_sttime
        print dt.datetime.now()
        
    return # MakeNetcdfFiles
#--------------------------------

if __name__=="__main__":

    '''
    Calls creation of netCDF files.  

    Uses sys.argv to input parameters - positional

    restart_id, end_id, extra

    Use direct call to test, and hard code the station ID here
    '''

    parser = argparse.ArgumentParser()
    parser.add_argument('--start', dest='STARTYEAR', action='store', default = 1931,
                        help='Start year, default = 1931')
    parser.add_argument('--end', dest='ENDYEAR', action='store', default = datetime.datetime.now().year,
                        help='End year, default = current')
    parser.add_argument('--restart_id', dest='restart_id', action='store', default = "",
                        help='Restart ID for truncated run, default = ""')
    parser.add_argument('--end_id', dest='end_id', action='store', default = "",
                        help='End ID for truncated run, default = ""')
    parser.add_argument('--extra', dest='extra', action='store_true', default = False,
                        help='Include extra parameters, default = False')

    args = parser.parse_args()

    STARTYEAR=int(args.STARTYEAR)
    ENDYEAR=int(args.ENDYEAR)
    restart_id=args.restart_id
    end_id=args.end_id
    Extra=args.extra

    print "Reading data from %s" % INPUT_DATA_DIR
    print "Writing data to %s" % OUTPUT_DATA_DIR
    print "Start year %i, End year %i (inclusive)" % (STARTYEAR,ENDYEAR)
    print "Restart ID = {}, End ID = {}, Include Extra parameters = {}".format(restart_id, end_id, Extra)
  
    MakeNetcdfFiles(STARTYEAR, ENDYEAR, restart_id=restart_id,end_id=end_id,Extra=Extra)
























