#!/usr/local/sci/bin/python
#*****************************
#
# All humidity variables
#
#
#************************************************************************
#                    SVN Info
#$Rev:: 107                                           $:  Revision of last commit
#$Author:: rdunn                                      $:  Author of last commit
#$Date:: 2016-07-29 15:38:24 +0100 (Fri, 29 Jul 2016) $:  Date of last commit
#************************************************************************

import numpy as np
import datetime as dt

# RJHD utilities
import qc_utils as utils

# note - using e_v for vapour pressure to avoid confusion with exponent

#************************************************************************
def calculate_e_v_wrt_water(t, P):
    '''
    Calculate vapour pressure wrt water

    Buck, A. L.: New equations for computing vapor pressure and enhancement factor, J. Appl. Meteorol., 20, 1527?1532, 1981.

    :param array t: temperature (or dewpoint temperature for saturation e_v) (deg C)
    :param array P: station level pressure (hPa)

    :returns array e_v: vapour pressure (or saturation vapour pressure if dewpoint temperature used) (hPa)
    '''

    f = 1 + (7.e-4) + ((3.46e-6) * P)

    e_v = 6.1121 * f * np.ma.exp(((18.729 - (t / 227.3)) * t) / (257.87 + t))

    e_v.mask = [any(tup) for tup in zip(t.mask, P.mask)]

    return e_v # calculate_e_v_wrt_water

#************************************************************************
def calculate_e_v_wrt_ice(t, P):
    '''
    Calculate vapour pressure wrt ice

    Buck, A. L.: New equations for computing vapor pressure and enhancement factor, J. Appl. Meteorol., 20, 1527?1532, 1981.

    :param array t: temperature (or dewpoint temperature for saturation e_v) (deg C)
    :param array P: station level pressure (hPa)

    :returns array e_v: vapour pressure (or saturation vapour pressure if dewpoint temperature used) (hPa)
    '''

    f = 1 + (3.e-4) + ((4.18e-6) * P)

    e_v = 6.1115 * f * np.ma.exp(((23.036 - (t / 333.7)) * t) / (279.82 + t))

    e_v.mask = [any(tup) for tup in zip(t.mask, P.mask)]

    return e_v # calculate_e_v_wrt_ice

#************************************************************************
def calculate_Td_wrt_water(e_v, P):
    '''
    Calculate dewpoint temperature wrt water

    Buck, A. L.: New equations for computing vapor pressure and enhancement factor, J. Appl. Meteorol., 20, 1527?1532, 1981.

    :param array e_v: vapour pressure (or saturation vapour pressure for temperature) (hPa)
    :param array P: station level pressure (hPa)

    :returns array Td: dewpoint temperature (or temperature if saturation e_v used) (deg C)
    '''
  
    f = 1 + (7.e-4) + ((3.46e-6) * P)

    a = 1
    b = (227.3 * np.ma.log(e_v / (6.1121 * f))) - (18.729 * 227.3)
    c = (257.87 * 227.3 * np.ma.log(e_v / (6.1121 * f)))
    
    Td = (-b - np.ma.sqrt(b**2 - (4 * a * c))) / (2 * a)

    Td.mask = [any(tup) for tup in zip(e_v.mask, P.mask)]

    return Td # calculate_Td_wrt_water

#************************************************************************
def calculate_Td_wrt_ice(e_v, P):
    '''
    Calculate dewpoint temperature wrt ice

    Buck, A. L.: New equations for computing vapor pressure and enhancement factor, J. Appl. Meteorol., 20, 1527?1532, 1981.

    :param array e_v: vapour pressure (or saturation vapour pressure for temperature) (hPa)
    :param array P: station level pressure (hPa)

    :returns array Td: dewpoint temperature (or temperature if saturation e_v used) (deg C)
    '''
  
    f = 1 + (3.e-4) + ((4.18e-6) * P)

    a = 1
    b = (333.7 * np.ma.log(e_v / (6.1115 * f))) - (23.036 * 333.7)
    c = (279.82 * 333.7 * np.ma.log(e_v / (6.1115 * f)))
    
    Td = (-b - np.ma.sqrt(b**2 - (4 * a * c))) / (2 * a)

    Td.mask = [any(tup) for tup in zip(e_v.mask, P.mask)]

    return Td # calculate_Td_wrt_ice

#************************************************************************
def calculate_Tw(e_v, P, Td, t):
    '''
    Calculate the webulb temperature
    
    Jensen, M. E., Burman, R. D., and Allen, R. G. (Eds.): Evapotranspiration and Irrigation Water Requirements: ASCE Manuals and Reports on Engineering Practices No. 70, American Society of Civil Engineers, New York, 360 pp., 1990.

    :param array e_v: vapour pressure (hPa)
    :param array P: station level pressure (hPa)
    :param array Td: dewpoint temperature (deg C)
    :param array t: dry-bulb temperature (deg C)
    
    :returns array Tw: wetbulb temperature (deg C)
    '''

    a = 0.000066 * P
    b = ((409.8 * e_v) / ((Td + 237.3)**2))

    Tw = (((a * t) + (b * Td)) / (a + b))

    Tw.mask = [any(tup) for tup in zip(e_v.mask, P.mask, Td.mask, t.mask)]

    return Tw # calculate_Tw

#************************************************************************
def calculate_q(e_v, P):
    '''
    Calculate specific humidity

    Peixoto, J. P. and Oort, A. H.: The climatology of relative humidity in the atmosphere, J. Climate, 9, 3443?3463, 1996.

    :param array e_v: vapour pressure (hPa)
    :param array P: station level pressure (hPa)
    
    :returns array q: specific humidity (g/kg)
    '''

    q = 1000. * ((0.622 * e_v) / (P - ((1 - 0.622) * e_v)))

    q.mask = [any(tup) for tup in zip(e_v.mask, P.mask)]

    return q # calculate_q

#************************************************************************
def calculate_e_v_from_q(q, P):
    '''
    Calculate vapour pressure from specific humidity

    :param array q: specific humidity (g/kg)
    :param array P: station level pressure (hPa)

    :returns array e_v: vapour pressure (hPa)
    '''

    e_v = ((q / 1000.) * P) / (0.622 + (0.388 * (q / 1000.)))

    e_v.mask = [any(tup) for tup in zip(q.mask, P.mask)]

    return e_v # calculate_e_v_from_q

#************************************************************************
def calculate_rh(e_v, es):
    '''
    Calculate relative humidity

    :param array e_v: vapour pressure (hPa)
    :param array es : saturation vapour pressure (hPa)

    :returns array rh: relative humidity (%rh)
    '''

    rh = (e_v / es) * 100.

    rh.mask = [any(tup) for tup in zip(e_v.mask, es.mask)]

    return rh # calculate_rh

#************************************************************************
def calculate_e_v_from_rh_es(es, rh):
    '''
    Calculate vapour pressure from rh and saturation vapour pressure

    :param array es: saturation vapour pressure (hPa)
    :param array rh: relative humidity (%rh)

    :returns array e_v: vapour pressure (hPa)
    '''

    e_v = (rh * es) / 100.

    e_v.mask = [any(tup) for tup in zip(es.mask, rh.mask)]

    return e_v # calculate_e_v_from_rh_es

#************************************************************************
def calculate_es_from_rh_e_v(e_v, rh):
    '''
    Calculate saturation vapour pressure from rh and vapour pressure

    :param array e_vs: vapour pressure (hPa)
    :param array rh: relative humidity (%rh)

    :returns array es: saturation vapour pressure (hPa)
    '''

    es = (e_v / rh) * 100.

    es.mask = [any(tup) for tup in zip(e_v.mask, rh.mask)]

    return es # calculate_es_from_rh_e_v

#************************************************************************
def fix_wrt_ice_or_water(temperatures, dewpoints, station_pressure):
    '''
    Calculate the vapour pressures and wet-bulb temperatures, adjusting
    for an ice- or water-bulb as appropriate from the Tw

    :param array temperatures: temperature array
    :param array dewpoints: dewpoint temperature array
    :param array station_pressure: station pressure array

    :returns: e_v, e_s, Tw - vapour pressure, saturation vapour pressure and wet-bulb temperature
    '''

    # get vapour pressures
    e_v = calculate_e_v_wrt_water(dewpoints, station_pressure)
    e_v_ice = calculate_e_v_wrt_ice(dewpoints, station_pressure)
    
    # saturation
    e_s = calculate_e_v_wrt_water(temperatures, station_pressure)
    e_s_ice = calculate_e_v_wrt_ice(temperatures, station_pressure)  
    
    # get wet-bulb temperatures
    Tw = calculate_Tw(e_v, station_pressure, dewpoints, temperatures)
    Tw_ice = calculate_Tw(e_v_ice, station_pressure, dewpoints, temperatures)
    
    # adjust for ice-bulbs
    e_v[Tw <= 0] = e_v_ice[Tw <= 0]
    e_s[Tw <= 0] = e_s_ice[Tw <= 0]
    Tw[Tw <= 0] = Tw_ice[Tw <= 0]
    
    return e_v, e_s, Tw # fix_wrt_ice_or_water


#************************************************************************
def get_station_level_pressure(station):
    '''
    Convert from sea level pressure back to station level
    List, R. J.: Smithsonian Meteorological Tables, Vol. 114, 6th Edn.,
    Smithsonian Institution, Washington DC, 268 pp., 1963.

    :param object station: station object
    
    :returns: station level pressure (hPa)
    '''

    slp = utils.apply_flags_to_mask(station, "slp")
    temperatures = utils.apply_flags_to_mask(station, "temperatures")

    t_in_kelvin = temperatures.data + 273.15

    elevation = station.elev

    pressure = slp.data * np.ma.power(t_in_kelvin / (t_in_kelvin + 0.0065 * elevation), 5.625)

    return pressure # get_station_level_pressure


#************************************************************************
def run_calcs(station, logfile, plots = False, diagnostics = False):
    '''
    Run the humidity calculations and add the attributes to the station file

    :param object station: station object
    :param file logfile: logfile to store outputs
    :param boolean diagnostics: output diagnostic information
    :param boolean plots: make a plot

    :returns: station - updated with humidity variables
    '''

    temperatures = utils.apply_flags_to_mask(station, "temperatures")
    dewpoints = utils.apply_flags_to_mask(station, "dewpoints")
   
    # adjust from sea-level to station-level
    station_pressure = get_station_level_pressure(station)
    

    e_v = utils.set_MetVar_attributes("vapor_pressure", "Vapor pressure calculated w.r.t water", "water_vapor_pressure", "hPa", temperatures.mdi, np.dtype('float64'))
    e_s = utils.set_MetVar_attributes("saturation_vapor_pressure", "Saturation vapor pressure calculated w.r.t. water", "water_vapor_pressure", "hPa", temperatures.mdi, np.dtype('float64'))
    Tw = utils.set_MetVar_attributes("wet_bulb_temperature", "Wet bulb temperatures nearest to reporting hour", "wet_bulb_temperature", "C", temperatures.mdi, np.dtype('float64'))
    q = utils.set_MetVar_attributes("specific_humidity", "Specific humidity", "specific_humidity", "g/kg", temperatures.mdi, np.dtype('float64'))
    rh = utils.set_MetVar_attributes("relative_humidity", "Relative humidity", "relative_humidity", "%rh", temperatures.mdi, np.dtype('float64'))


    # sort the vapour pressures and wet-bulb --> ice or water?
    e_v.data, e_s.data, Tw.data = fix_wrt_ice_or_water(temperatures.data, dewpoints.data, station_pressure)

    
    # get relative and specific humidity
    q.data = calculate_q(e_v.data, station_pressure)
    rh.data = calculate_rh(e_v.data, e_s.data)
    
    if plots or diagnostics:
        print "Humidity variables calculated, setting attributes\n"
    else:
        logfile.write("Humidity variables calculated, setting attributes\n")

    setattr(station, "vapour_pressure", e_v)
    setattr(station, "saturation_vapour_pressure", e_s)
    setattr(station, "wetbulb_temperature", Tw)
    setattr(station, "specific_humidity", q)
    setattr(station, "relative_humidity", rh)

    station = utils.append_history(station, "Humidity Calculations")

    return station # run_calcs
