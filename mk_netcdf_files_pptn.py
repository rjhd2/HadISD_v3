#!/usr/local/sci/bin/python
#************************************************************************
#                    SVN Info
#$Rev:: 219                                           $:  Revision of last commit
#$Author:: rdunn                                      $:  Author of last commit
#$Date:: 2019-05-20 16:56:47 +0100 (Mon, 20 May 2019) $:  Date of last commit
#************************************************************************

'''
mk_netcdf_files_pptn.py invoked by typing::

  python2.7 mk_netcdf_files_pptn.py --start 1931 --endyear 2018 --endmonth 10 --restart_id 000000-99999 --end_id 999999-99999 [--extra] [--doCanada]

Python script to read the ISD ASCII text format and output netcdf files

Runs with no inputs in current version. 

Input arguments:

--start             First year to include

--endyear           Last year to include

--endmonth          Last month of last year

--restart_id        First station to process

--end_id            Last station to process

--extra             [False] Include extra variables

--doCanada          [False] Include extra Canadian metadata
'''

#*************************************
# Compared to IDL output using compare.py on 30May2012 and found to match
# except for total_cloud_flags - but on investigation with raw ISD files
# the python extraction is the correct one. RJHD
#
# Could change data types to match IDL, but adapting QC so that works with
# floats and doubles as appropriate.RJHD
#
# Logic change to match IDL so that overwrite only acted upon if writing
# real data, rather than missing. RJHD
#*************************************

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

# RJHD utils
from set_paths_and_vars import *


# Globals
INTMDI=-999
FLTMDI=-1.e30
hours = True # output time axis in hours.

NCDC_FLAGS={'A':10,'U':11,'P':12,'I':13,'M':14,'C':15,'R':16, 'E':17, 'J':18}

#---------------------------------------------------------------------


#************************************************************************
def ReadStations(filename):
    """ 
    Read Station Information

    :param string filename: name and location of input file

    :returns: numpy array of file contents

    Use numpy genfromtxt reading to read in all station 
    data in ID,Lat,Lon,Elev list
    """
    return np.genfromtxt(filename, dtype=(str)) #  ReadStations


#************************************************************************
def ReadComposites(filename):
    """ 
    Read Composite Station Information

    :param string filename: name and location of input file

    :returns: list of lists containing composites
    """
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
    """
    Tests if string is an integer

    :param string s: string to test
    :returns: boolean if string is valid integer
    """
    try: 
        int(s)
        return True
    except ValueError:
        return False # RepresentsInt

#************************************************************************
def TimeMatch(timearray,testtime, lower,upper):
    """
    Do matching of np array to find time step

    :param array timearray: np array of timestamps
    :param float testtime: timestep to find
    :return: int of location
    """
    return np.argwhere(timearray[lower:upper]==testtime)[0] # TimeMatch


#************************************************************************
def ExtractValues(missing,line,location,length,test,divisor=1.,flagoffset=0, doflag=True):
    """
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
    """

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
    """
    Test if need to extract the data

    :param float/int data: data to test
    :param float/int missing: missing value
    :param boolean overwrite: to overwrite or not

    :returns: boolean if condition met
    """

    if data==missing or overwrite:
        return True
    else:
        return False # TestToExtract


#************************************************************************
def ExtractionProcess(data, flags, time, missing, missingtest, line, location, length,divisor=1.,flagoffset=0, doflag=True):
    """
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
    """

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
    """
    Write note to dubious file list.
    
    :param string outfile: filename to be written to
    :param string infile: filename of dubious file
    :param string code: text identifier of variables being tested
    :param string station: station ID being processed
    :param string time: time of the dubious data

    :returns: int of flag status.
    """
    flagged=0
    try:
        with open(outfile,'a') as of:
            of.write(station+' '+time+' '+code+' variables are first, but not nec. only problem '+infile+'\n')
            of.close()
            flagged=1
    except IOError:
        # file doesn't exist as yet, so make a new one
        with open(outfile,'w') as of:
            of.write(station+' '+time+' '+code+' variables are first, but not nec. only problem '+infile+'\n')
            of.close()
            flagged=1
        
        
    return flagged # WriteDubious

#************************************************************************
def SortClouds(cloud_cover,cloud_flags, time, amounts, flags, clouds):
    """
    Convert the raw cloud data into oktas for each level

    :param array cloud_cover: final cloud_cover array
    :param array cloud_flags: final cloud_flags array
    :param int time_loc: time stamp
    :param array amounts: raw cloud amounts - in oktas
    :param array flags: raw cloud flags
    :param array clouds: locations of where cloud heights match this level
    """

    if len(clouds)>=1 and cloud_cover[time]==INTMDI:
        cloud_cover[time]=np.max(amounts[clouds])
        cloud_flags[time]=np.max(flags[clouds])

    return # SortClouds                              

#************************************************************************
def SortClouds2(cloud_cover,cloud_flags, time, amounts, amounts2, flags, clouds):
    """
    Convert the raw cloud data into oktas and for each level

    :param array cloud_cover: final cloud_cover array
    :param array cloud_flags: final cloud_flags array
    :param int time_loc: time stamp
    :param array amounts: raw cloud amounts - in other units - see ISD documentation
    :param array amounts2: raw cloud amounts - in oktas
    :param array flags: raw cloud flags
    :param array clouds: locations of where cloud heights match this level
    """

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
def ProcessPrecip(depth, flag, condition, time, line, exists):
    """
    Extract the precipitation information - allows neater if clause.

    :param array depth: precipitation depth array
    :param array flag: precipitation flag array
    :param array condition: precipitation condition array
    :param int time: time stamp
    :param str line: the line to process
    :param int exists: the index of the entry in the line
    """

    ExtractionProcess(depth, flag, time, FLTMDI, '9999', line, 
                      exists+5, 4, doflag=True, flagoffset=1, divisor=10.)
    # these two pass in empty strings as missing data test
    dummyflag = 0
    ExtractionProcess(condition, dummyflag, time, '', '9', line,
                      exists+9, 1, doflag=False)

    return # ProcessPrecip                         

#************************************************************************
def WriteAttributes(variable,long_name,cell_methods,missing_value,units,axis,vmin,vmax,coordinates,standard_name = ''):
    """
    Write given attributes into ncdf variable
    
    :param object variable: netcdf Variable
    :param string long_name: long_name value for variable to be written
    :param string cell_methods: cell_methods value for variable to be written
    :param float/int missing_value: missing_value value for variable to be written
    :param string units: units value for variable to be written
    :param string axis: axis value for variable to be written
    :param float/int vmin: valid_min value for variable to be written
    :param float/int vmax: valid_max value for variable to be written
    :param string standard_name: standard_name value for variable to be written
    :param string coordinates: coordinates to associate to variable
    """
    variable.long_name=long_name
    variable.cell_methods=cell_methods
    variable.missing_value=missing_value
#    variable.axis=axis # 12/1/17 RJHD - not required for CF compliance.
    variable.units=units
    variable.valid_min=vmin
    variable.valid_max=vmax
    variable.coordinates=coordinates
    
    if standard_name != '':
        variable.standard_name=standard_name

    return # WriteAttributes

#************************************************************************
def WriteFlagAttributes(variable,long_name,missing_value,axis):
    """
    Write given attributes into ncdf variable
    
    :param object variable: netcdf Variable
    :param string long_name: long_name value for variable to be written
    :param float/int missing_value: missing_value value for variable to be written
    :param string axis: axis value for variable to be written
    """
    variable.long_name=long_name
    variable.missing_value=missing_value
    variable.units="1"
#    variable.axis=axis # 12/1/17 RJHD - not required for CF compliance.

    # for future [September 2015]
    # http://cfconventions.org/Data/cf-conventions/cf-conventions-1.6/build/cf-conventions.html#flags

    return # WriteFlagAttributes

#************************************************************************
def write_coordinates(outfile, short_name, standard_name, long_name, units, axis, data, coordinate_length = 1, do_zip = True):
    """
    Write coordinates as variables

    :param str outfile: output netcdf file
    :param str short_name: netcdf short_name
    :param str standard_name: netcdf standard_name
    :param str long_name: netcdf long_name
    :param str units: netcdf units
    :param str axis: netcdf axis
    :param flt data: coordinate 
    :param int coordinate_length: length of dimension
    :param bool do_zip: allow for zipping
    """


    if "coordinate_length" not in outfile.dimensions:
        coord_dim = outfile.createDimension('coordinate_length', coordinate_length)

    nc_var = outfile.createVariable(short_name, np.dtype('float'), ('coordinate_length',), zlib = do_zip)
    nc_var.standard_name = standard_name
    nc_var.long_name = long_name
    nc_var.units = units
    nc_var.axis = axis

    if short_name == "alt":
        nc_var.positive = "up"

    nc_var[:] = data


    return # write_coordinates

#************************************************************************
def MakeNetcdfFiles(STARTYEAR, ENDYEAR, ENDMONTH, restart_id="", end_id="", do_zip = True, Extra = False, doCanada = True): 
    """
    Parse the ASCII files and do the NetCDF file creation

    :param int STARTYEAR: first year of data
    :param int ENDYEAR: last year of data (inclusive)
    :param int ENDMONTH: last month of data (inclusive)
    :param string restart_id: string for starting station, default=""
    :param string end_id: string for ending station, default=""
    :param boolean do_zip: make netCDF4 files with internal zipping
    :param boolean Extra: setting to extract extra variables
    :param boolean doCanada: use extra information from ECCC
    """

    print "Note to Self on re-write (26/2/2018)"
    print "Use the netcdf_procs to write the netcdf file."
    print "Have text lookup file for variable attributes (allow for free text if necessary)"
    print "Standardise across the precip1-4 variables OR read all 4 and split into 1/3/6/12/24 hourly accumulations."
    print "Cope with updated ISD format."
 
    StationInfo=ReadStations(os.path.join(INPUT_FILE_LOCS, STATION_LIST))

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

    Composites=ReadComposites(os.path.join(INPUT_FILE_LOCS, MERGER_LIST))

    # adjust the end points for clarity and ease of use
    if ENDMONTH == 12:
        # easier to use the 1st January at 00:00
        ENDYEAR += 1
        ENDMONTH = 1
    else:
        # easier to use the 1st of the following month at 00:00
        ENDMONTH += 1


    DaysBetween=dt.datetime(ENDYEAR, ENDMONTH, 1, 0, 0) - dt.datetime(STARTYEAR, 1, 1, 0, 0)
    HoursBetween=int(DaysBetween.days*24.)

    TimeStamps=np.linspace(0,HoursBetween-1,HoursBetween) # keep in integer hours

    ValidYears=np.arange(STARTYEAR, ENDYEAR+1)
    dubiousfile=LOG_OUTFILE_LOCS+'dubious_ISD_data_files.txt'

    # read in Canadian station list
    if doCanada:
        Canadian_stations_info = np.genfromtxt(INPUT_FILE_LOCS + "Canada_time_ranges.dat", dtype=(str), delimiter = [12,20,20])
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
        precip1_depth=np.zeros(HoursBetween)
        precip1_condition=['null' for i in range(HoursBetween)]
        precip1_flag=np.zeros(HoursBetween, dtype=np.int)
        precip2_depth=np.zeros(HoursBetween)
        precip2_condition=['null' for i in range(HoursBetween)]
        precip2_flag=np.zeros(HoursBetween, dtype=np.int)
        precip3_depth=np.zeros(HoursBetween)
        precip3_condition=['null' for i in range(HoursBetween)]
        precip3_flag=np.zeros(HoursBetween, dtype=np.int)
        precip6_depth=np.zeros(HoursBetween)
        precip6_condition=['null' for i in range(HoursBetween)]
        precip6_flag=np.zeros(HoursBetween, dtype=np.int)
        precip9_depth=np.zeros(HoursBetween)
        precip9_condition=['null' for i in range(HoursBetween)]
        precip9_flag=np.zeros(HoursBetween, dtype=np.int)
        precip12_depth=np.zeros(HoursBetween)
        precip12_condition=['null' for i in range(HoursBetween)]
        precip12_flag=np.zeros(HoursBetween, dtype=np.int)
        precip15_depth=np.zeros(HoursBetween)
        precip15_condition=['null' for i in range(HoursBetween)]
        precip15_flag=np.zeros(HoursBetween, dtype=np.int)
        precip18_depth=np.zeros(HoursBetween)
        precip18_condition=['null' for i in range(HoursBetween)]
        precip18_flag=np.zeros(HoursBetween, dtype=np.int)
        precip24_depth=np.zeros(HoursBetween)
        precip24_condition=['null' for i in range(HoursBetween)]
        precip24_flag=np.zeros(HoursBetween, dtype=np.int)
        slp=np.zeros(HoursBetween)
        slp_flag=np.zeros(HoursBetween, dtype=np.int)
        stnlp=np.zeros(HoursBetween)
        stnlp_flag=np.zeros(HoursBetween, dtype=np.int)
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
        precip1_depth.fill(FLTMDI)
        precip1_flag.fill(INTMDI)
        precip2_depth.fill(FLTMDI)
        precip2_flag.fill(INTMDI)
        precip3_depth.fill(FLTMDI)
        precip3_flag.fill(INTMDI)
        precip6_depth.fill(FLTMDI)
        precip6_flag.fill(INTMDI)
        precip9_depth.fill(FLTMDI)
        precip9_flag.fill(INTMDI)
        precip12_depth.fill(FLTMDI)
        precip12_flag.fill(INTMDI)
        precip15_depth.fill(FLTMDI)
        precip15_flag.fill(INTMDI)
        precip18_depth.fill(FLTMDI)
        precip18_flag.fill(INTMDI)
        precip24_depth.fill(FLTMDI)
        precip24_flag.fill(INTMDI)
        slp.fill(FLTMDI)
        slp_flag.fill(INTMDI)
        stnlp.fill(FLTMDI)
        stnlp_flag.fill(INTMDI)
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
                raw_files.extend(glob.glob(ISD_DATA_LOCS+'station725/'+cstn+'*'))
            else:
                raw_files.extend(glob.glob(ISD_DATA_LOCS+'station'+cstn[0:2]+'s/'+cstn+'*'))


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
                            
                            # adjust the day (and month)
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

                            if dt_time > dt.datetime(ENDYEAR, ENDMONTH, 1, 0, 0):
                                # data beyond the end of the timeperiod.  Ignored
                                if not done_print:
                                    print "skipping rest of station {} as {} after end of selected period ({})".format(raw_station, dt.datetime.strftime(dt_time, "%Y-%m-%d %H:%M"), dt.datetime.strftime(dt.datetime(ENDYEAR, ENDMONTH, 1, 0, 0), "%Y-%m-%d %H:%M"))
                                    done_print = True
                                continue

                            if doCanada:
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

                                # Station Level Pressure
                                text_ident='MA1' # 3,5,1,5,1
                                exists=cleanline.find(text_ident)
                                if exists!=-1:
                                    try:
                                        if RepresentsInt(cleanline[exists+3]):

                                            if Extract:
                                                ExtractionProcess(stnlp, stnlp_flag,time_loc,FLTMDI,'99999',cleanline,
                                                                  exists+9, 5, divisor = 10.)
                                                
                                    except IndexError:
                                        # string following data marker doesn't exist
                                        if dubious_flagged==0:
                                            dubious_flagged=WriteDubious(dubiousfile,rfile,text_ident, station, string_obs_time)
                                        

                                elif dubious_flagged==0:
                                    dubious_flagged=WriteDubious(dubiousfile,rfile,text_ident, station, string_obs_time)

                                # CLOUDs
                                text_ident='GF1' # 3,2,2,1,2,1,2,1,5,1,2,1,2,1
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

                                text_ident='GA1' # 3,2,1,6,1,2,1
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


                                text_ident='GD1' # 3,1,2,1,6,1,1
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
                                text_ident='AY1' # 3,1,1,2,1
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

                                text_ident='AZ1' # 3,1,1,2,1
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

                                # PRECIP - split into different accumulation periods; run all four entries
                                for text_ident in ["AA1", "AA2", "AA3", "AA4"]: # 3,2,4,1,1
                                    exists=cleanline.find(text_ident)

                                    if exists!=-1:
                                        try:
                                            if RepresentsInt(cleanline[exists+3]):
                                                if Extract:
                                                    period = np.zeros(1) # need an array to use this routine
                                                    ExtractionProcess(period, dummyflag, 0, INTMDI, '99', cleanline,
                                                                      exists+3, 2, doflag = False)
                                                    # test if the accumulation period exists.
                                                    if period[0] == 1:
                                                        ProcessPrecip(precip1_depth, precip1_flag, precip1_condition, time_loc, cleanline, exists)
                                                    elif period[0] == 2:
                                                        ProcessPrecip(precip2_depth, precip2_flag, precip2_condition, time_loc, cleanline, exists)
                                                    elif period[0] == 3:
                                                        ProcessPrecip(precip3_depth, precip3_flag, precip3_condition, time_loc, cleanline, exists)
                                                    elif period[0] == 6:
                                                        ProcessPrecip(precip6_depth, precip6_flag, precip6_condition, time_loc, cleanline, exists)
                                                    elif period[0] == 9:
                                                        ProcessPrecip(precip9_depth, precip9_flag, precip9_condition, time_loc, cleanline, exists)
                                                    elif period[0] == 12:
                                                        ProcessPrecip(precip12_depth, precip12_flag, precip12_condition, time_loc, cleanline, exists)
                                                    elif period[0] == 15:
                                                        ProcessPrecip(precip15_depth, precip15_flag, precip15_condition, time_loc, cleanline, exists)
                                                    elif period[0] == 18:
                                                        ProcessPrecip(precip18_depth, precip18_flag, precip18_condition, time_loc, cleanline, exists)
                                                    elif period[0] == 24:
                                                        ProcessPrecip(precip24_depth, precip24_flag, precip24_condition, time_loc, cleanline, exists)


                                        except IndexError:
                                            if dubious_flagged==0:
                                                dubious_flagged=WriteDubious(dubiousfile,rfile,text_ident, station, string_obs_time)


                                    elif dubious_flagged==0:
                                        dubious_flagged=WriteDubious(dubiousfile,rfile,text_ident, station, string_obs_time)


                                # SUN DURATION

                                text_ident='GJ1' # 3,4,1
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

                                for text_ident in ['OA1','OA2','OA3','OA4']: # 3,1,2,4,1
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
                                    text_ident='AW1'  # 3,2,1
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

                                    for text_ident in ['AY1','AY2','AZ1','AZ2']: # 3,1,1,2,1

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


                                    # EXTREME TEMPERATURES

                                    # RJHD - 11 March 2014 - these could be converted to Tmax and Tmin using the code information

                                    for text_ident in ['KA1','KA2']:  #3,3,1,5,1
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
        precip1_depth=precip1_depth[mask_vals]
        precip1_condition=np.array(precip1_condition)[mask_vals]
        precip1_flag=precip1_flag[mask_vals]
        precip2_depth=precip2_depth[mask_vals]
        precip2_condition=np.array(precip2_condition)[mask_vals]
        precip2_flag=precip2_flag[mask_vals]
        precip3_depth=precip3_depth[mask_vals]
        precip3_condition=np.array(precip3_condition)[mask_vals]
        precip3_flag=precip3_flag[mask_vals]
        precip6_depth=precip6_depth[mask_vals]
        precip6_condition=np.array(precip6_condition)[mask_vals]
        precip6_flag=precip6_flag[mask_vals]
        precip9_depth=precip9_depth[mask_vals]
        precip9_condition=np.array(precip9_condition)[mask_vals]
        precip9_flag=precip9_flag[mask_vals]
        precip12_depth=precip12_depth[mask_vals]
        precip12_condition=np.array(precip12_condition)[mask_vals]
        precip12_flag=precip12_flag[mask_vals]
        precip15_depth=precip15_depth[mask_vals]
        precip15_condition=np.array(precip15_condition)[mask_vals]
        precip15_flag=precip15_flag[mask_vals]
        precip18_depth=precip18_depth[mask_vals]
        precip18_condition=np.array(precip18_condition)[mask_vals]
        precip18_flag=precip18_flag[mask_vals]
        precip24_depth=precip24_depth[mask_vals]
        precip24_condition=np.array(precip24_condition)[mask_vals]
        precip24_flag=precip24_flag[mask_vals]
        slp=slp[mask_vals]
        slp_flag=slp_flag[mask_vals]
        stnlp=stnlp[mask_vals]
        stnlp_flag=stnlp_flag[mask_vals]
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
            maximum_temp_period=maximum_temp_period[mask_vals]
            maximum_temp_value=maximum_temp_value[mask_vals]
            maximum_temp_flags=maximum_temp_flags[mask_vals]
            minimum_temp_period=minimum_temp_period[mask_vals]
            minimum_temp_value=minimum_temp_value[mask_vals]
            minimum_temp_flags=minimum_temp_flags[mask_vals]

        netcdf_filename = os.path.join(NETCDF_DATA_LOCS, "hadisd.{}_19310101-{}_{}_raw.nc".format(LONG_VERSION, END_TIME, station))

        if do_zip:
            netcdf_outfile = ncdf.Dataset(netcdf_filename,'w', format='NETCDF4')
        else:
            netcdf_outfile = ncdf.Dataset(netcdf_filename,'w', format='NETCDF3_CLASSIC')

        time=netcdf_outfile.createDimension('time',len(times_out))
        char_len=netcdf_outfile.createDimension('character_length',4)
        long_char_len=netcdf_outfile.createDimension('long_character_length',12)
        coords_len=netcdf_outfile.createDimension('coordinate_length',1)


        # write the coordinates
        write_coordinates(netcdf_outfile, "latitude", "latitude", "station_latitude", "degrees_north", "Y", StationLat[st])
        write_coordinates(netcdf_outfile, "longitude", "longitude", "station_longitude", "degrees_east", "X", StationLon[st])
        write_coordinates(netcdf_outfile, "elevation", "surface_altitude", "vertical distance above the surface", "meters", "Z", StationElv[st])

        # station ID as base variable
        nc_var = netcdf_outfile.createVariable("station_id", np.dtype('S1'), ('long_character_length',), zlib = do_zip)
#        nc_var.standard_name = "station_identification_code"
        nc_var.long_name = "Station ID number"
        nc_var[:] = ncdf.stringtochar(StationIDs[st])

        # create variables
        timesvar=netcdf_outfile.createVariable('time','f8',('time',), zlib = do_zip)
        # lonsvar=netcdf_outfile.createVariable('lon','f8',('coordinate_length',), zlib = do_zip)
        # latsvar=netcdf_outfile.createVariable('lat','f8',('coordinate_length',), zlib = do_zip)
        # altsvar=netcdf_outfile.createVariable('alt','f8',('coordinate_length',), zlib = do_zip)
        # idsvar=netcdf_outfile.createVariable('station_id','S1',('long_character_length',), zlib = do_zip)
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
        ppt1dvar=netcdf_outfile.createVariable('precip1_depth','f8',('time',), zlib = do_zip)
        ppt1cvar=netcdf_outfile.createVariable('precip1_condition','S1',('time','character_length',), zlib = do_zip)
        ppt1fvar=netcdf_outfile.createVariable('precip1_flag','i4',('time',), zlib = do_zip)
        ppt2dvar=netcdf_outfile.createVariable('precip2_depth','f8',('time',), zlib = do_zip)
        ppt2cvar=netcdf_outfile.createVariable('precip2_condition','S1',('time','character_length',), zlib = do_zip)
        ppt2fvar=netcdf_outfile.createVariable('precip2_flag','i4',('time',), zlib = do_zip)
        ppt3dvar=netcdf_outfile.createVariable('precip3_depth','f8',('time',), zlib = do_zip)
        ppt3cvar=netcdf_outfile.createVariable('precip3_condition','S1',('time','character_length',), zlib = do_zip)
        ppt3fvar=netcdf_outfile.createVariable('precip3_flag','i4',('time',), zlib = do_zip)
        ppt6dvar=netcdf_outfile.createVariable('precip6_depth','f8',('time',), zlib = do_zip)
        ppt6cvar=netcdf_outfile.createVariable('precip6_condition','S1',('time','character_length',), zlib = do_zip)
        ppt6fvar=netcdf_outfile.createVariable('precip6_flag','i4',('time',), zlib = do_zip)
        ppt9dvar=netcdf_outfile.createVariable('precip9_depth','f8',('time',), zlib = do_zip)
        ppt9cvar=netcdf_outfile.createVariable('precip9_condition','S1',('time','character_length',), zlib = do_zip)
        ppt9fvar=netcdf_outfile.createVariable('precip9_flag','i4',('time',), zlib = do_zip)
        ppt12dvar=netcdf_outfile.createVariable('precip12_depth','f8',('time',), zlib = do_zip)
        ppt12cvar=netcdf_outfile.createVariable('precip12_condition','S1',('time','character_length',), zlib = do_zip)
        ppt12fvar=netcdf_outfile.createVariable('precip12_flag','i4',('time',), zlib = do_zip)
        ppt15dvar=netcdf_outfile.createVariable('precip15_depth','f8',('time',), zlib = do_zip)
        ppt15cvar=netcdf_outfile.createVariable('precip15_condition','S1',('time','character_length',), zlib = do_zip)
        ppt15fvar=netcdf_outfile.createVariable('precip15_flag','i4',('time',), zlib = do_zip)
        ppt18dvar=netcdf_outfile.createVariable('precip18_depth','f8',('time',), zlib = do_zip)
        ppt18cvar=netcdf_outfile.createVariable('precip18_condition','S1',('time','character_length',), zlib = do_zip)
        ppt18fvar=netcdf_outfile.createVariable('precip18_flag','i4',('time',), zlib = do_zip)
        ppt24dvar=netcdf_outfile.createVariable('precip24_depth','f8',('time',), zlib = do_zip)
        ppt24cvar=netcdf_outfile.createVariable('precip24_condition','S1',('time','character_length',), zlib = do_zip)
        ppt24fvar=netcdf_outfile.createVariable('precip24_flag','i4',('time',), zlib = do_zip)
        slpvar=netcdf_outfile.createVariable('slp','f8',('time',), zlib = do_zip)
        slpfvar=netcdf_outfile.createVariable('slp_flag','i4',('time',), zlib = do_zip)
        stnlpvar=netcdf_outfile.createVariable('stnlp','f8',('time',), zlib = do_zip)
        stnlpfvar=netcdf_outfile.createVariable('stnlp_flag','i4',('time',), zlib = do_zip)
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
            maxtpvar=netcdf_outfile.createVariable('maximum_temp_period','f8',('time',), zlib = do_zip)
            maxtvvar=netcdf_outfile.createVariable('maximum_temp_value','f8',('time',), zlib = do_zip)
            maxtfvar=netcdf_outfile.createVariable('maximum_temp_flag','i4',('time',), zlib = do_zip)
            mintpvar=netcdf_outfile.createVariable('minimum_temp_period','f8',('time',), zlib = do_zip)
            mintvvar=netcdf_outfile.createVariable('minimum_temp_value','f8',('time',), zlib = do_zip)
            mintfvar=netcdf_outfile.createVariable('minimum_temp_flag','i4',('time',), zlib = do_zip)

        # variables attributes
        print "Writing Attributes"

        timesvar.long_name='time_of_measurement'
        timesvar.standard_name='time'
        if hours:
            timesvar.units='hours since {}'.format(dt.datetime.strftime(dt.datetime(STARTYEAR,1,1,0,0), "%Y-%m-%d %H:%M"))
        else:
            timesvar.units='days since {}'.format(dt.datetime.strftime(dt.datetime(STARTYEAR,1,1,0,0), "%Y-%m-%d %H:%M"))
        timesvar.axis='T'
        timesvar.calendar='gregorian'
        timesvar.valid_min=0.
        timesvar.start_year = "{}".format(STARTYEAR)
        timesvar.end_year = "{}".format(ENDYEAR)
        timesvar.start_month = "1"
        timesvar.end_month = "{}".format(ENDMONTH)
#        timesvar.coordinates = "time"

        # lonsvar.standard_name = "longitude"
        # lonsvar.long_name = "station_longitude"
        # lonsvar.units = "degrees_east"
        # lonsvar.axis = "X"

        # latsvar.standard_name = "latitude"
        # latsvar.long_name = "station_latitude"
        # latsvar.units = "degrees_north"
        # latsvar.axis = "Y"

        # altsvar.long_name = "vertical distance above the surface"
        # altsvar.standard_name = "height"
        # altsvar.units = "meters"
        # altsvar.positive = "up"
        # altsvar.axis = "Z"
        
        # idsvar.standard_name = "station_identification_code"
        # idsvar.long_name = "Station ID number"
        # idsvar.cf_role='timeseries_id'

#        stationsvar.standard_name='station_identification_code'
        stationsvar.long_name='Primary source for timestep (may be multiple sources for composite stations). USAF-WBAN from ISD source'
#        stationsvar.units='USAF - WBAN from ISD source'
        stationsvar.missing_value='null'

        try:
            tmin,tmax=np.min(temperatures[np.where(temperatures != FLTMDI)[0]]),np.max(temperatures[np.where(temperatures != FLTMDI)[0]])
        except ValueError:
            tmin,tmax=FLTMDI,FLTMDI
        WriteAttributes(tempsvar,'Dry bulb air temperature at screen height (~2m)','latitude: longitude: time: point (nearest to reporting hour)',FLTMDI,'degree_Celsius','T',tmin,tmax,'latitude longitude elevation',standard_name = 'surface_temperature')
        WriteFlagAttributes(tempsflagsvar,'ISD flags for temperature - see ISD documentation',INTMDI,'T')

        try:
            dmin,dmax=np.min(dewpoints[np.where(dewpoints != FLTMDI)[0]]),np.max(dewpoints[np.where(dewpoints != FLTMDI)[0]])
        except ValueError:
            dmin,dmax=FLTMDI,FLTMDI
        WriteAttributes(dewsvar,'Dew point temperature at screen height (~2m)','latitude: longitude: time: point (nearest to reporting hour)',FLTMDI,'degree_Celsius','T',dmin,dmax,'latitude longitude elevation',standard_name =  'dew_point_temperature')
        WriteFlagAttributes(dewsflagsvar,'ISD flags for dewpoint temperature - see ISD documentation',INTMDI,'T')

        WriteAttributes(tcvar,'Total cloud cover (oktas)','latitude: longitude: time: point (derived in priority order GA, GF, GD - see ISD documentation, nearest to reporting hour)', INTMDI, '1', 'T', 0,8, 'latitude longitude elevation',standard_name = "cloud_area_fraction")
        WriteFlagAttributes(tcfvar,'ISD flags for total cloud - see ISD documentation',INTMDI,'T')

        WriteAttributes(lcvar,'Low cloud cover (oktas)','latitude: longitude: time: point (derived in priority order GA, GF, GD - see ISD documentation, nearest to reporting hour)', INTMDI, '1', 'T', 0,8, 'latitude longitude elevation',standard_name = "low_type_cloud_area_fraction")
        WriteFlagAttributes(lcfvar,'ISD flags for low cloud - see ISD documentation',INTMDI,'T')

        WriteAttributes(mcvar,'Mid cloud cover (oktas)','latitude: longitude: time: point (derived in priority order GA, GF, GD - see ISD documentation, nearest to reporting hour)', INTMDI, '1', 'T', 0,8, 'latitude longitude elevation',standard_name = "medium_type_cloud_area_fraction")
        WriteFlagAttributes(mcfvar,'ISD flags for mid cloud - see ISD documentation',INTMDI,'T')

        WriteAttributes(hcvar,'High cloud cover (oktas)','latitude: longitude: time: point (derived in priority order GA, GF, GD - see ISD documentation, nearest to reporting hour)', INTMDI, '1', 'T', 0,8, 'latitude longitude elevation',standard_name = "high_type_cloud_area_fraction")
        WriteFlagAttributes(hcfvar,'ISD flags for high cloud - see ISD documentation',INTMDI,'T')

        try:
            cbmin,cbmax=np.min(cloud_base[np.where(cloud_base != INTMDI)[0]]),np.max(cloud_base[np.where(cloud_base != INTMDI)[0]])
        except ValueError:
            cbmin,cbmax=FLTMDI,FLTMDI
        WriteAttributes(cbvar,'Cloud base of lowest cloud layer','latitude: longitude: time: point (nearest to reporting hour)', INTMDI, 'meters', 'T', cbmin, cbmax, 'latitude longitude elevation',standard_name = 'cloud_base_altitude')
        WriteFlagAttributes(cbfvar,'ISD flags for cloud base - see ISD documentation',INTMDI,'T')

        try:
            wsmin,wsmax=np.min(windspeeds[np.where(windspeeds != FLTMDI)[0]]),np.max(windspeeds[np.where(windspeeds != FLTMDI)[0]])
        except ValueError:
            wsmin,wsmax=FLTMDI,FLTMDI
        WriteAttributes(wsvar,'Wind speed at mast height (~10m)','latitude: longitude: time: point (nearest to reporting hour)', FLTMDI, 'meters per second', 'T', wsmin, wsmax,'latitude longitude elevation',standard_name = 'wind_speed')
        WriteFlagAttributes(wsfvar,'ISD flags for windspeed - see ISD documentation',INTMDI,'T')

        WriteAttributes(wdvar,'Wind Direction at mast height (~10m)','latitude: longitude: time: point (nearest to reporting hour)', INTMDI, 'degree', 'T', 0, 360, 'latitude longitude elevation',standard_name = 'wind_from_direction')
        WriteFlagAttributes(wdfvar,'ISD flags for wind direction - see ISD documentation',INTMDI,'T')

        WriteAttributes(pswx1var,'Reported past significant weather phenomena','latitude: longitude: point (interval: 1 day)', INTMDI, '1', 'T', 0, 9,'latitude longitude elevation')
        WriteFlagAttributes(pswx1fvar,'ISD flags for reported past significant weather - see ISD documentation',INTMDI,'T')
        WriteAttributes(pswx1pvar,'Reported period over which significant weather report was recorded','latitude: longitude: point (interval: 1 day)', INTMDI, 'Hours', 'T', 0, 24,'latitude longitude elevation')

        WriteAttributes(ppt1dvar,'Depth of Precipitation Reported in 1 hour (from all four ISD fields)','latitude: longitude: time: sum ', FLTMDI, 'mm', 'T', 0, 999.8,'latitude longitude elevation', standard_name = 'lwe_thickness_of_precipitation_amount')
        WriteFlagAttributes(ppt1cvar,'Precipitation Code for 1 hour accumulation (denotes if trace amount)', 'null','T')
        WriteFlagAttributes(ppt1fvar,'ISD flags for 1 hour precipitation accumulation (from all four ISD fields - see ISD documentation)', INTMDI,'T')

        WriteAttributes(ppt2dvar,'Depth of Precipitation Reported in 2 hour (from all four ISD fields)','latitude: longitude: time: sum ', FLTMDI, 'mm', 'T', 0, 999.8,'latitude longitude elevation', standard_name = 'lwe_thickness_of_precipitation_amount')
        WriteFlagAttributes(ppt2cvar,'Precipitation Code for 2 hour accumulation (denotes if trace amount)', 'null','T')
        WriteFlagAttributes(ppt2fvar,'ISD flags for 2 hour precipitation accumulation (from all four ISD fields - see ISD documentation)', INTMDI,'T')

        WriteAttributes(ppt3dvar,'Depth of Precipitation Reported in 3 hours (from all four ISD fields)','latitude: longitude: time: sum ', FLTMDI, 'mm', 'T', 0, 999.8,'latitude longitude elevation', standard_name = 'lwe_thickness_of_precipitation_amount')
        WriteFlagAttributes(ppt3cvar,'Precipitation Code for 3 hour accumulation (denotes if trace amount)', 'null','T')
        WriteFlagAttributes(ppt3fvar,'ISD flags for 3 hour precipitation accumulation (from all four ISD fields - see ISD documentation)', INTMDI,'T')

        WriteAttributes(ppt6dvar,'Depth of Precipitation Reported in 6 hour (from all four ISD fields)','latitude: longitude: time: sum ', FLTMDI, 'mm', 'T', 0, 999.8,'latitude longitude elevation', standard_name = 'lwe_thickness_of_precipitation_amount')
        WriteFlagAttributes(ppt6cvar,'Precipitation Code for 6 hour accumulation (denotes if trace amount)', 'null','T')
        WriteFlagAttributes(ppt6fvar,'ISD flags for 6 hour precipitation accumulation (from all four ISD fields - see ISD documentation)', INTMDI,'T')

        WriteAttributes(ppt9dvar,'Depth of Precipitation Reported in 9 hour (from all four ISD fields)','latitude: longitude: time: sum ', FLTMDI, 'mm', 'T', 0, 999.8,'latitude longitude elevation', standard_name = 'lwe_thickness_of_precipitation_amount')
        WriteFlagAttributes(ppt9cvar,'Precipitation Code for 9 hour accumulation (denotes if trace amount)', 'null','T')
        WriteFlagAttributes(ppt9fvar,'ISD flags for 9 hour precipitation accumulation (from all four ISD fields - see ISD documentation)', INTMDI,'T')

        WriteAttributes(ppt12dvar,'Depth of Precipitation Reported in 12 hour (from all four ISD fields)','latitude: longitude: time: sum ', FLTMDI, 'mm', 'T', 0, 999.8,'latitude longitude elevation', standard_name = 'lwe_thickness_of_precipitation_amount')
        WriteFlagAttributes(ppt12cvar,'Precipitation Code for 12 hour accumulation (denotes if trace amount)', 'null','T')
        WriteFlagAttributes(ppt12fvar,'ISD flags for 12 hour precipitation accumulation (from all four ISD fields - see ISD documentation)', INTMDI,'T')

        WriteAttributes(ppt15dvar,'Depth of Precipitation Reported in 15 hour (from all four ISD fields)','latitude: longitude: time: sum ', FLTMDI, 'mm', 'T', 0, 999.8,'latitude longitude elevation', standard_name = 'lwe_thickness_of_precipitation_amount')
        WriteFlagAttributes(ppt15cvar,'Precipitation Code for 15 hour accumulation (denotes if trace amount)', 'null','T')
        WriteFlagAttributes(ppt15fvar,'ISD flags for 15 hour precipitation accumulation (from all four ISD fields - see ISD documentation)', INTMDI,'T')

        WriteAttributes(ppt18dvar,'Depth of Precipitation Reported in 18 hour (from all four ISD fields)','latitude: longitude: time: sum ', FLTMDI, 'mm', 'T', 0, 999.8,'latitude longitude elevation', standard_name = 'lwe_thickness_of_precipitation_amount')
        WriteFlagAttributes(ppt18cvar,'Precipitation Code for 18 hour accumulation (denotes if trace amount)', 'null','T')
        WriteFlagAttributes(ppt18fvar,'ISD flags for 18 hour precipitation accumulation (from all four ISD fields - see ISD documentation)', INTMDI,'T')

        WriteAttributes(ppt24dvar,'Depth of Precipitation Reported in 24 hour (from all four ISD fields)','latitude: longitude: time: sum ', FLTMDI, 'mm', 'T', 0, 999.8,'latitude longitude elevation', standard_name = 'lwe_thickness_of_precipitation_amount')
        WriteFlagAttributes(ppt24cvar,'Precipitation Code for 24 hour accumulation (denotes if trace amount)', 'null','T')
        WriteFlagAttributes(ppt24fvar,'ISD flags for 24 hour precipitation accumulation (from all four ISD fields - see ISD documentation)', INTMDI,'T')

        try:
            smin,smax=np.min(slp[np.where(slp != FLTMDI)[0]]),np.max(slp[np.where(slp != FLTMDI)[0]])
        except ValueError:
            smin,smax=FLTMDI,FLTMDI

        WriteAttributes(slpvar,'Reported Sea Level Pressure at screen height (~2m)','latitude: longitude: time: point (nearest to reporting hour)',FLTMDI, 'hPa', 'T', smin, smax, 'latitude longitude elevation',standard_name = 'air_pressure_at_sea_level')
        WriteFlagAttributes(slpfvar,'ISD flags for slp field - see ISD documentation',INTMDI,'T')

        WriteAttributes(stnlpvar,'Reported Station Level Pressure at screen height (~2m)','latitude: longitude: time: point (nearest to reporting hour)',FLTMDI, 'hPa', 'T', smin, smax, 'latitude longitude elevation',standard_name = 'surface_air_pressure')
        WriteFlagAttributes(stnlpfvar,'ISD flags for stnlp field - see ISD documentation',INTMDI,'T')

        WriteAttributes(sdvar,'Reported Sunshine Duration','latitude: longitude: time: point (nearest to reporting hour)', INTMDI, 'minutes', 'T', 0, 6000, 'latitude longitude elevation', standard_name = 'duration_of_sunshine')
        WriteFlagAttributes(sdfvar,'ISD flags sun duration field - see ISD documentation',INTMDI,'T')

        WriteAttributes(wgstpvar,'Period of Maximum Wind Gust Speed', 'latitude: longitude: time: point (nearest to reporting hour)', INTMDI, 'Hours', 'T', 0, 48, 'latitude longitude elevation') #,standard_name = 'period_of_wind_gust')
        WriteAttributes(wgstvvar,'Wind Gust Speed at mast height (~10m)','latitude: longitude: time: point (nearest to reporting hour)',FLTMDI, 'meters per second', 'T', 0, 200.0, 'latitude longitude elevation',standard_name = 'wind_speed_of_gust')
        WriteFlagAttributes(wgstfvar,'ISD flags for wind gust field - see ISD documentation', INTMDI,'T')

        if Extra:
            WriteFlagAttributes(wtvar,'Wind observation type - see ISD documentation','null','T')

            WriteAttributes(pswx2var,'Station reports of past significant weather phenomena (2)','latitude: longitude: point (interval: 1 day)', INTMDI, '1', 'T', 0, 9,'latitude longitude elevation')
            WriteFlagAttributes(pswx2fvar,'ISD flags for reported past significant weather - see ISD documentation',INTMDI,'T')
            WriteAttributes(pswx2pvar,'Period of significant weather report','latitude: longitude: point (interval: 1 day)', INTMDI, 'hour', 'T', 0, 24,'latitude longitude elevation')

            WriteAttributes(swxvar,'Station reports of present significant weather phenomena','latitude: longitude: point (interval: 1 day)', INTMDI, '1', 'T', 0, 99,'latitude longitude elevation')
            WriteFlagAttributes(swxfvar,'ISD flags for reported present significant weather - see ISD documentation',INTMDI,'T')

            try:
                xtmin,xtmax=np.min(maximum_temp_value[np.where(maximum_temp_value != FLTMDI)[0]]),np.max(maximum_temp_value[np.where(maximum_temp_value != FLTMDI)[0]])
            except ValueError:
                xtmin,xtmax=FLTMDI,FLTMDI
            WriteAttributes(maxtpvar,'Reported period over which maximum temperature was recorded','latitude: longitude: point (interval: 1 day)', FLTMDI, 'hour', 'T', 0, 48,'latitude longitude elevation')
            WriteAttributes(maxtvvar,'Dry bulb maximum temperature reported over time period','latitude: longitude: time: point (interval: 1 day)', FLTMDI,'degrees_Celsius','T',xtmin,xtmax)
            WriteFlagAttributes(maxtfvar,'ISD flags for maximum temperature field - see ISD documentation', INTMDI, 'T')


            try:
                ntmin,ntmax=np.min(minimum_temp_value[np.where(minimum_temp_value != FLTMDI)[0]]),np.max(minimum_temp_value[np.where(minimum_temp_value != FLTMDI)[0]])
            except ValueError:
                ntmin,ntmax=FLTMDI,FLTMDI
            WriteAttributes(mintpvar,'Reported period over which minimum temperature was recorded','latitude: longitude: point (interval: 1 day)', FLTMDI, 'hour', 'T', 0, 48,'latitude longitude elevation')
            WriteAttributes(mintvvar,'Dry bulb minimum temperature reported over time period','latitude: longitude: time: point (interval: 1 day)', FLTMDI,'degree_Celsius','T',ntmin,ntmax)
            WriteFlagAttributes(mintfvar,'ISD flags for minimum temperature field - see ISD documentation', INTMDI, 'T')



        # global attributes
        netcdf_outfile.station_information='Where station is a composite the station id refers to the primary source used in the timestep and does apply to all elements'

        netcdf_outfile.file_created=dt.datetime.strftime(dt.datetime.now(), "%a %b %d, %H:%M %Y")
        netcdf_outfile.station_id=station
        netcdf_outfile.latitude=StationLat[st]
        netcdf_outfile.longitude=StationLon[st]
        netcdf_outfile.elevation=StationElv[st]
        netcdf_outfile.Conventions="CF-1.6"
        netcdf_outfile.date_created = dt.datetime.strftime(dt.datetime.now(), "%Y-%m-%d, %H:%M")
        netcdf_outfile.history = "Created by mk_netcdf_files.py \n"
 
        print "Writing data to netcdf file"

        # write data into file
        # ncdf.stringtochar changes np array of strings
        #    to times x character_length array of single
        #    character strings


        timesvar[:]=times_out
#        lonsvar[:]=StationLon[st]
#        latsvar[:]=StationLat[st]
#        altsvar[:]=StationElv[st]
#        idsvar[:]=StationIDs[st]
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
        ppt1dvar[:]=precip1_depth
        ppt1cvar[:]=ncdf.stringtochar(precip1_condition)
        ppt1fvar[:]=precip1_flag
        ppt2dvar[:]=precip2_depth
        ppt2cvar[:]=ncdf.stringtochar(precip2_condition)
        ppt2fvar[:]=precip2_flag
        ppt3dvar[:]=precip3_depth
        ppt3cvar[:]=ncdf.stringtochar(precip3_condition)
        ppt3fvar[:]=precip3_flag
        ppt6dvar[:]=precip6_depth
        ppt6cvar[:]=ncdf.stringtochar(precip6_condition)
        ppt6fvar[:]=precip6_flag
        ppt9dvar[:]=precip9_depth
        ppt9cvar[:]=ncdf.stringtochar(precip9_condition)
        ppt9fvar[:]=precip9_flag
        ppt12dvar[:]=precip12_depth
        ppt12cvar[:]=ncdf.stringtochar(precip12_condition)
        ppt12fvar[:]=precip12_flag
        ppt15dvar[:]=precip15_depth
        ppt15cvar[:]=ncdf.stringtochar(precip15_condition)
        ppt15fvar[:]=precip15_flag
        ppt18dvar[:]=precip18_depth
        ppt18cvar[:]=ncdf.stringtochar(precip18_condition)
        ppt18fvar[:]=precip18_flag
        ppt24dvar[:]=precip24_depth
        ppt24cvar[:]=ncdf.stringtochar(precip24_condition)
        ppt24fvar[:]=precip24_flag
        slpvar[:]=slp
        slpfvar[:]=slp_flag
        stnlpvar[:]=stnlp
        stnlpfvar[:]=stnlp_flag
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
        print "\n END"
        
    return # MakeNetcdfFiles
#--------------------------------

if __name__=="__main__":

    """
    Calls creation of netCDF files.  

    Uses sys.argv to input parameters - positional

    restart_id, end_id, extra

    Use direct call to test, and hard code the station ID here
    """

    parser = argparse.ArgumentParser()
    parser.add_argument('--start', dest='STARTYEAR', action='store', default = 1931,
                        help='Start year, default = 1931', type = int)
    parser.add_argument('--endyear', dest='ENDYEAR', action='store', default = datetime.datetime.now().year,
                        help='End year (inclusive), default = current', type = int)
    parser.add_argument('--endmonth', dest='ENDMONTH', action='store', default = datetime.datetime.now().month,
                        help='End month (inclusive), default = current', type = int)
    parser.add_argument('--restart_id', dest='restart_id', action='store', default = "",
                        help='Restart ID for truncated run, default = ""')
    parser.add_argument('--end_id', dest='end_id', action='store', default = "",
                        help='End ID for truncated run, default = ""')
    parser.add_argument('--extra', dest='extra', action='store_true', default = False,
                        help='Include extra parameters, default = False')
    parser.add_argument('--doCanada', dest='doCanada', action='store_false', default = True,
                        help='Include adjustments for Canadian stations, default = True')

    args = parser.parse_args()

    print "\n Making NetCDF files from ISD ASCII files \n"

    print "Reading data from %s" % ISD_DATA_LOCS
    print "Writing data to %s" % NETCDF_DATA_LOCS
    print "Start year %i, End year %i (inclusive)" % (args.STARTYEAR,args.ENDYEAR)
    print "Restart ID = {}, End ID = {}, Include Extra parameters = {}".format(args.restart_id, args.end_id, args.extra)
  
    MakeNetcdfFiles(args.STARTYEAR, args.ENDYEAR, args.ENDMONTH, restart_id = args.restart_id, end_id = args.end_id, Extra = args.extra, doCanada = args.doCanada)
























