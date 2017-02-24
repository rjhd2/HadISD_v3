#!/usr/local/sci/bin/python
#*****************************
#
# netCDF utilities for Python QC.
#
#
#************************************************************************
#                    SVN Info
#$Rev:: 112                                           $:  Revision of last commit
#$Author:: rdunn                                      $:  Author of last commit
#$Date:: 2017-01-13 14:47:17 +0000 (Fri, 13 Jan 2017) $:  Date of last commit
#************************************************************************


# assume class already set up?

import numpy as np
import scipy as sp
import datetime as dt
import netCDF4 as ncdf
import os


# RJHD 
import qc_utils as utils

#*******************************************************
def read(filename, station, var_list, opt_var_list = [], diagnostics = False, read_input_station_id = True, read_qc_flags = True, read_flagged_obs = True):
    '''
    Reads the netcdf file and appends attributes to station object
    
    variables to be read passed in as list
    '''

    try:
        ncfile = ncdf.Dataset(filename,'r')

    except RuntimeError:
        print "File not available {}".format(filename)
        raise RuntimeError
    
#     try:
#         # old method
#         if ncfile.latitude != station.lat:
#             print "Station longitudes do not match"
#             raise RuntimeError
#         if ncfile.longitude != station.lon:
#             print "Station longitudes do not match"
#             raise RuntimeError
#         if ncfile.elevation != station.elev:
#             print "Station elevations do not match"
#             raise RuntimeError
#         if ncfile.id != station.id:
#             print "Station IDs do not match"
#             raise RuntimeError
#     except AttributeError:
#         # new method
#         # change in how these are stored
#         if not utils.nearly_equal(ncfile.variables['lat'][:][0], station.lat):
#             print "Station latitudes do not match"
#             raise RuntimeError
#         if not utils.nearly_equal(ncfile.variables['lon'][:][0], station.lon):
#             print "Station longitudes do not match"
#             raise RuntimeError
#         if not utils.nearly_equal(ncfile.variables['alt'][:][0], station.elev):
#             print "Station elevations do not match"
#             raise RuntimeError
# #        if ncfile.variables['station_id'][:][0] != station.id:
# #            print "Station IDs do not match"
# #            raise RuntimeError
      
    # print ncfile.variables

    # if optional/carry through variables given, then set to extract these too
    if opt_var_list != []:
        full_var_list = np.append(var_list, opt_var_list)
    else:
        full_var_list = var_list


    if read_input_station_id:
        final_var_list = np.append(full_var_list, ["time", "input_station_id"])
    else:
        final_var_list = np.append(full_var_list, ["time"])


    for variable in final_var_list:
        if diagnostics: print "reading {}".format(variable)

        var = ncfile.variables[variable] # this is a masked array
        
        this_var = utils.MetVar(variable, var.long_name)

        try:
            this_var.units = var.units
        except AttributeError:
            pass

        this_var.dtype = var.dtype

        if variable in ["input_station_id","precip1_condition","windtypes","precip2_condition","precip3_condition","precip4_condition"]:
            # this is a slow step for input_station_id
            this_var.data = np.ma.array(["".join(i) for i in var[:]])
        else:
            this_var.data = np.ma.array(var[:]) # keep as masked array           
            
        # section for non-default netcdf attributes
        if variable == "time":
            this_var.calendar = var.calendar
        
        try:
            this_var.mdi = var.missing_value
        except AttributeError:
            if variable in ["temperatures"]:
                this_var.mdi=-1.e30
            elif variable in ["total_cloud_cover"]:
                this_var.mdi=-999
            elif variable in ["input_station_id"]:
                this_var.mdi="null"
            pass

        try:
            this_var.valid_max = var.valid_max
        except AttributeError:
            pass
        try:
            this_var.valid_min = var.valid_min
        except AttributeError:
            pass
        try:
            this_var.standard_name = var.standard_name
        except AttributeError:
            pass
        try:
            this_var.coordinates = var.coordinates
        except AttributeError:
            pass
        try:
            this_var.cell_methods = var.cell_methods
        except AttributeError:
            pass
        try:
            this_var.fdi = var.flagged_value
        except AttributeError:
            if variable in ["temperatures","dewpoints","slp","windspeeds"]:
                this_var.fdi=-2.e30
            elif variable in ["total_cloud_cover","low_cloud_cover","mid_cloud_cover","high_cloud_cover", "winddirs"]:
                this_var.fdi=-888
            pass

        try:
            this_var.flags = var.flags
        except AttributeError:
            this_var.flags = np.zeros(len(this_var.data))
            pass

        # now read in all info for this variable
        # exec("station."+this_var.name+" = this_var") # replaced with setattr
        setattr(station, variable, this_var)

    # read in the qc_flags array
    if read_qc_flags == True:
        try:
            qc_flags = ncfile.variables["quality_control_flags"]
            setattr(station, "qc_flags", qc_flags[:])
            
        except KeyError:
            if diagnostics:
                print "no QC flags available"

    # read in the flagged_obs array   
    if read_flagged_obs == True:
        try:
            flagged_obs = ncfile.variables["flagged_value"]
        except KeyError:
            print "no flagged obs available in netcdf file"
            # if doesn't exist, make an empty array
            flagged_obs = np.zeros([len(station.time.data),len(var_list)])

            # and set empty array to be correct MDI for variable
            for v,var in enumerate(var_list):
                st_var = getattr(station, var)
                flagged_obs[:,v] = st_var.mdi

        # push array into relevant attribute
        for v,var in enumerate(var_list):
            st_var = getattr(station, var)

            try:
                dummy = flagged_obs[:,v].mask
                st_var.flagged_obs = flagged_obs[:,v].data
            except AttributeError:
                # no mask so just put in the data
                st_var.flagged_obs = flagged_obs[:,v]

    # read in reporting statistics - just to carry through
    try:
        reporting_stats = ncfile.variables["reporting_stats"]

        # push array into relevant attributes
        for v,var in enumerate(var_list):
            st_var = getattr(station, var)
        
            st_var.reporting_stats = reporting_stats[v]

    except KeyError:
        print "no reporting stats available in netcdf file"
        
    # other station attributes:
    try:
        station.history = ncfile.history
    except AttributeError:
        station.history = ""

    ncfile.close()

    return # read


#************************************************************************
def read_global_attributes(attr_file):
    '''
    Reads attributes file and returns a dictionary of key/values
    '''
        
    try:
        with open(attr_file,'r') as infile:        
            lines = infile.readlines()
        
    except IOError:
        print "Attributes file not found at " + attr_file
    
    
    attributes = {}
    
    for line in lines:
        split_line = line.split()
        
        attributes[split_line[0]] = " ".join(split_line[1:])    
        
    return attributes

#************************************************************************
def write_coordinates(outfile, short_name, standard_name, long_name, units, axis, data, coordinate_length = 1, do_zip = True, least_significant_digit = None):
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
    :param int least_significant_digit: smallest reliable decimal place
    """

    if "coordinate_length" not in outfile.dimensions:
        coord_dim = outfile.createDimension('coordinate_length', coordinate_length)

    nc_var = outfile.createVariable(short_name, np.dtype('float'), ('coordinate_length',), zlib = do_zip, least_significant_digit = least_significant_digit)
    try:
        nc_var.standard_name = standard_name
    except AttributeError:
        pass
    nc_var.long_name = long_name
    nc_var.units = units
    nc_var.axis = axis

    if short_name == "elevation":
        nc_var.positive = "up"

    nc_var[:] = data


    return

#************************************************************************
def write(filename, station, var_list, attr_file, processing_date = '', qc_code_version = '', opt_var_list = [], compressed = [], do_zip = True, write_QC_flags = True, write_flagged_obs = True, least_significant_digit = 0):
    '''
    Writes the netcdf file.  

    compressed - compress the time axis
    '''

    # remove file
    if os.path.exists(filename):
        os.remove(filename)

    # decide on format
    if do_zip:
        outfile = ncdf.Dataset(filename,'w', format='NETCDF4')
    else:
        outfile = ncdf.Dataset(filename,'w', format='NETCDF3_CLASSIC')
        
    # sort length of time axis if compressed time
    if compressed != []:
        time_dim = outfile.createDimension('time',len(np.where(compressed == True)[0]))
    else:
        time_dim = outfile.createDimension('time',len(station.time.data))

    # sort character dimensions
    long_character_length = 12
    character_length = 4
    
    long_char_dim = outfile.createDimension('long_character_length',long_character_length)
    char_dim = outfile.createDimension('character_length',character_length)
 
    # sort other dimensions if required
    if write_QC_flags:
        qc_test_length = station.qc_flags.shape[1]
        test_dim = outfile.createDimension('test',qc_test_length)
        
    if write_flagged_obs:
        flagged_dim = outfile.createDimension('flagged',len(var_list))

    # set up reporting stats if attribute available
    try:
        st_var = getattr(station,var_list[0])
        reportingT_dim = outfile.createDimension('reporting_t',st_var.reporting_stats.shape[0]) # N months
        reportingV_dim = outfile.createDimension('reporting_v',len(var_list))
        reporting2_dim = outfile.createDimension('reporting_2', 2) # accuracy and frequency

    except AttributeError:
        print "no reporting information - cannot set up dimensions"


    # write station coordinates
    write_coordinates(outfile, "longitude", "longitude", "station longitude", "degrees_east", "X", station.lon)
    write_coordinates(outfile, "latitude", "latitude", "station latitude", "degrees_north", "Y", station.lat)
    write_coordinates(outfile, "elevation", "surface_altitude", "vertical distance above the surface", "meters", "Z", station.elev)
   

    # station ID as base variable
    nc_var = outfile.createVariable("station_id", np.dtype('S1'), ('long_character_length',), zlib = do_zip)
#    nc_var.standard_name = "station_identification_code"
    nc_var.long_name = "Station ID number"
    nc_var[:] = station.id



    # if optional/carry through variables given, then set to extract these too
    if opt_var_list != []:
        full_var_list = np.append(var_list, opt_var_list)
    else:
        full_var_list = var_list

    # spin through all attributes and put them into the file
    for var in np.append(full_var_list, ["time", "input_station_id"]):
        
        st_var = getattr(station,var)
        
        if var == "input_station_id":
            nc_var = outfile.createVariable(st_var.name, st_var.dtype, ('time','long_character_length',), zlib = do_zip)
        elif var in ["precip1_condition","windtypes","precip2_condition","precip3_condition","precip4_condition",""]:
            nc_var = outfile.createVariable(st_var.name, st_var.dtype, ('time','character_length',), zlib = do_zip, fill_value = st_var.mdi)           
        else:

            if var == "time":
                nc_var = outfile.createVariable(st_var.name, st_var.dtype, ('time',), zlib = do_zip, least_significant_digit = least_significant_digit)
            elif least_significant_digit != 0:
                nc_var = outfile.createVariable(st_var.name, st_var.dtype, ('time',), zlib = do_zip, least_significant_digit = least_significant_digit, fill_value = st_var.mdi)
            else:
                nc_var = outfile.createVariable(st_var.name, st_var.dtype, ('time',), zlib = do_zip, fill_value = st_var.mdi)

        nc_var.long_name = st_var.long_name
        #nc_var.axis = st_var.axis #  apparently not needed according to CF Checker
        try:
            nc_var.units = st_var.units
        except AttributeError:
            pass
        try:
            nc_var.missing_value = st_var.mdi
        except AttributeError:
            pass
        try:
            nc_var.flagged_value = st_var.fdi
        except AttributeError:
            pass
        try:
            nc_var.valid_min = st_var.valid_min
        except AttributeError:
            pass
        try:
            nc_var.valid_max = st_var.valid_max
        except AttributeError:
            pass
        try:
            nc_var.standard_name = st_var.standard_name
        except AttributeError:
            pass
        try:
            nc_var.coordinates = st_var.coordinates
        except AttributeError:
            pass
        try:
            nc_var.cell_methods = st_var.cell_methods
        except AttributeError:
            pass

        # extra attributes

        # have to expand string array out into individual characters if appropriate
        if compressed != []:
            if var == "input_station_id":
                reformatted_string_data = st_var.data[compressed]
                reformatted_string_data = np.ma.array([list("{:12s}".format(i)) for i in reformatted_string_data])
                nc_var[:] = reformatted_string_data
            else:
                nc_var[:] = st_var.data[compressed]
        
        else:
            if var == "input_station_id":
                reformatted_string_data = np.ma.array([list("{:12s}".format(i)) for i in st_var.data.filled()])
                nc_var[:] = reformatted_string_data[:]
            else:
                nc_var[:] = st_var.data
        
        if var == "time":
            nc_var.calendar = st_var.calendar
            
    # write QC flag information if available
    if write_QC_flags:
        # for future - http://cfconventions.org/Data/cf-conventions/cf-conventions-1.6/build/cf-conventions.html#flags
        # use a binary list in base 2 to build up unique flag single values that allow unpicking of each test status

        try:
            nc_var = outfile.createVariable("quality_control_flags", np.dtype('float'), ('time','test',), zlib = do_zip, fill_value = -999)
            nc_var.units = '1'
            nc_var.missing_value = -999
            nc_var.long_name = "Quality Control status for individual obs"
#            nc_var.cell_methods = "lat: lon: time: point"
#            nc_var.coordinates = "time test"
            if compressed != []:
                nc_var[:] = station.qc_flags[compressed,:]
            else:
                nc_var[:] = station.qc_flags
        except AttributeError:
            print "qc_flags attribute doesn't exist"

    # combine all flagged observations together to output as single array in netcdf file - if available
    # only done on masking?
    if write_flagged_obs:
        try:
            flagged_obs = np.zeros([len(station.time.data),len(var_list)])
            for v,var in enumerate(var_list):
                st_var = getattr(station, var)

                flagged_obs[:,v] = st_var.flagged_obs       
            
            if least_significant_digit != 0:
                nc_var = outfile.createVariable("flagged_value", np.dtype('float'), ('time','flagged',), zlib = do_zip, least_significant_digit = least_significant_digit, fill_value = -1.e30)
            else:
                nc_var = outfile.createVariable("flagged_value", np.dtype('float'), ('time','flagged',), zlib = do_zip, fill_value = -1.e30)
            nc_var.units = '1'
            nc_var.missing_value = -1.e30
            nc_var.long_name = "Observation Values removed by QC flags "+" ".join(var_list)
            nc_var.cell_methods = "latitude: longitude: time: point"
            nc_var.coordinates = "latitude longitude elevation"
            if compressed != []:
                nc_var[:] = flagged_obs[compressed,:]
            else:
                nc_var[:] = flagged_obs
        except AttributeError:
            print "no flagged observations."

    # combine all reporting accuracies together to output as single array in netcdf file - if available
    try:
        reporting_stats = []
        for v,var in enumerate(var_list):
            st_var = getattr(station, var)
            
            reporting_stats += [st_var.reporting_stats]

        reporting_stats = np.array(reporting_stats)

        nc_var = outfile.createVariable("reporting_stats", np.dtype('float'), ('reporting_v','reporting_t','reporting_2',), zlib = do_zip)
        nc_var.units = '1'
        nc_var.missing_value = -999
        nc_var.long_name = "Reporting frequency and accuracy for each month and variable in following order: "+" ".join(var_list)
        nc_var[:] = reporting_stats
        
    except AttributeError:
        print "no reporting accuracy information"   

    # Global Attributes
    # from file
    attribs = read_global_attributes(attr_file)
    
    for attr in attribs:
        
        outfile.__setattr__(attr, attribs[attr])
    
    # from code

    # removed these as added as coordinates
#    outfile.station_id = station.id
#    outfile.latitude = station.lat
#    outfile.longitude = station.lon
#    outfile.elevation = station.elev
    outfile.date_created = dt.datetime.strftime(dt.datetime.now(), "%Y-%m-%d, %H:%M")
    outfile.qc_code_version = qc_code_version
    outfile.station_information = 'Where station is a composite the station id refers to the primary source used in the timestep and does apply to all elements'
    outfile.Conventions = 'CF-1.6' 
    outfile.Metadata_Conventions = 'Unidata Dataset Discovery v1.0,CF Discrete Sampling Geometries Conventions'
    outfile.featureType = 'timeSeries'
    outfile.processing_date = processing_date
    outfile.history = station.history
    
    outfile.close()
       
    return # write

