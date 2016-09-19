#!/usr/local/sci/bin/python
#*****************************
#
# All heat stress indices
#
#
#************************************************************************
#                    SVN Info
#$Rev:: 107                                           $:  Revision of last commit
#$Author:: rdunn                                      $:  Author of last commit
#$Date:: 2016-07-29 15:38:24 +0100 (Fri, 29 Jul 2016) $:  Date of last commit
#************************************************************************

import numpy as np

# RJHD utilities
import qc_utils as utils

# note - using e_v for vapour pressure to avoid confusion with exponent

#************************************************************************
def calculate_thi(t, rh):
    '''
    Calculates the Temperature-Humidity Index as used for cattle heat stress

    NRC, 1971 - see Dikmen & Hansen 2009

    :param array t: dry bulb temperature (deg C)
    :param array rh: %rh 

    :returns: thi (0-100)
    '''

    thi = ((1.8 * t) + 32.) - ((0.55 - (0.0055 * rh)) * ((1.8 * t) - 26.))

    return thi # calculate_thi


#************************************************************************
def calculate_wbgt(t, e_v):
    '''
    Calculate the pseudo wet-bulb globe temperature (WGBT)

    ACSM. 1984. Prevention of thermal injuries during distance running.  American College of Sports Medicine. Medical Journal of Australia, 141: 876?879.

    :param array t: drybulb temperature (deg C)
    :param array e: vapour pressure (hPa)

    :returns array wbgt: pseudo wet-bulb globe temperature (deg C)
    '''


    wbgt = (0.567 * t) + (0.393 * e_v) + 3.94

    return wbgt # calculate_wbgt


#************************************************************************
def calculate_e_from_wbgt(t, wbgt):
    '''
    Calculate vapour pressure (e) from pseudo WBGT

    ACSM. 1984. Prevention of thermal injuries during distance running.  American College of Sports Medicine. Medical Journal of Australia, 141: 876?879.

    :param array t: drybulb temperature (deg C)
    :param array wbgt: wet bulb globe temperature (deg C)

    :returns array e_v: vapour pressure (hPa)
    '''

    e_v = (wbgt - 3.94 - (0.567 * t)) / 0.393	

    return e_v # calculate_e_from_wbgt

#************************************************************************
def calculate_humidex(t, e_v):
    '''
    Calculate Humidex

    Masterson J, Richardson FA. 1979. Humidex, A method of Quantifying Human Discomfort Due to Excessive Heat and Humidity. Environment Canada: Downsview, Ontario, pp 45.

    :param array t: drybulb temperature (deg C)
    :param array e_v: vapour pressure (hPa)

    :returns array h: humidex (deg C?)
    '''

    h = t + (0.5555 * (e_v - 10.))

    return h # calculate_humidex

#************************************************************************
def calculate_apparent_t(t, e_v, w):
    '''
    Calculate apparent temperature

    Steadman N. 1994. Norms of apparent temperature in Australia. Australian Meteorology Magazine 43: 1?16.

    :param array t: drybulb temperature (deg C)
    :param array e_v: vapour pressure (hPa)
    :param array w: wind speed (m/s)

    :returns array ta: apparent temperature (deg C)
    '''
    

    apparent_t = t + (0.33 * e_v) - (0.7 * w) - 4.

    return apparent_t # calculate_apparent_t

#************************************************************************
def calculate_heat_index(t, rh):
    '''
    Calculate heat index

    Rothfusz LP. 1990. The Heat Index Equation, SR Technical Attachment, 94-19, pp 6.
    http://www.hpc.ncep.noaa.gov/html/heatindex_equation.shtml

    :param array t: drybulb temperature (deg C)
    :param array rh: relative humidity (%)

    :returns array heat_index: heat index (deg C)

    NOTE: This gives meaning less results when t < 26.6 oC and rh < 40%
    '''

    t_fahrenheit = t * (9./5.) + 32

    heat_index_fahrenheit = -42.379 + (2.04901523 * t_fahrenheit) + (10.14333127 * rh) + \
        (-0.22475541 * t_fahrenheit * rh) + (-0.006837837 * t_fahrenheit * t_fahrenheit) + \
        (-0.05481717 * rh * rh) + (0.001228747 * t_fahrenheit * t_fahrenheit * rh) + \
        (0.00085282 * t_fahrenheit * rh * rh) + (-0.00000199 * t_fahrenheit * t_fahrenheit * rh * rh)

    # apply adjustment factors
    locs = np.ma.where(np.ma.logical_and((rh < 13), (t_fahrenheit > 80), (t_fahrenheit < 112)))
    if len(locs[0]) > 0:
        heat_index_fahrenheit[locs] = heat_index_fahrenheit[locs] - (((13.- rh[locs]) / 4.) * np.ma.sqrt((17. - np.ma.abs(t_fahrenheit[locs] - 95.)) / 17.)) 


    locs = np.ma.where(np.ma.logical_and((rh > 85), (t_fahrenheit > 80), (t_fahrenheit < 87)))
    if len(locs[0]) > 0:
        heat_index_fahrenheit[locs] = heat_index_fahrenheit[locs] - (((rh[locs ] - 85) / 10.) * ((87. - t_fahrenheit[locs]) / 5.))
    
    locs = np.ma.where(heat_index_fahrenheit < 80)
    if len(locs[0]) > 0:
        heat_index_fahrenheit[locs] = 0.5 * (t_fahrenheit[locs] + 61. + ((t_fahrenheit[locs] - 68.) * 1.2) + (rh[locs] * 0.094))
 
    heat_index = (heat_index_fahrenheit - 32) / (9./5.)

    locs = np.ma.where(t < 26.6667) # 80F
    if len(locs[0]) > 0:
        heat_index[locs] = -99
    locs = np.ma.where(rh < 40.0)
    if len(locs[0]) > 0:
        heat_index[locs] = -99
 

    '''
    FURTHER NOTES:

The computation of the heat index is a refinement of a result obtained by multiple regression analysis carried out by Lans P. Rothfusz and described in a 1990 National Weather Service (NWS) Technical Attachment (SR 90-23).  The regression equation of Rothfusz is

HI = -42.379 + 2.04901523*T + 10.14333127*RH - .22475541*T*RH - .00683783*T*T - .05481717*RH*RH + .00122874*T*T*RH + .00085282*T*RH*RH - .00000199*T*T*RH*RH 

where T is temperature in degrees F and RH is relative humidity in percent.  HI is the heat index expressed as an apparent temperature in degrees F.  If the RH is less than 13% and the temperature is between 80 and 112 degrees F, then the following adjustment is subtracted from HI:

ADJUSTMENT = [(13-RH)/4]*SQRT{[17-ABS(T-95.)]/17} 

where ABS and SQRT are the absolute value and square root functions, respectively.  On the other hand, if the RH is greater than 85% and the temperature is between 80 and 87 degrees F, then the following adjustment is added to HI:

ADJUSTMENT = [(RH-85)/10] * [(87-T)/5] 

The Rothfusz regression is not appropriate when conditions of temperature and humidity warrant a heat index value below about 80 degrees F. In those cases, a simpler formula is applied to calculate values consistent with Steadman's results:

HI = 0.5 * {T + 61.0 + [(T-68.0)*1.2] + (RH*0.094)} 

In practice, the simple formula is computed first and the result averaged with the temperature. If this heat index value is 80 degrees F or higher, the full regression equation along with any adjustment as described above is applied.

The Rothfusz regression is not valid for extreme temperature and relative humidity conditions beyond the range of data considered by Steadman. '''

    return heat_index # calculate_heat_index

#************************************************************************
def run_calcs(station, logfile, plots = False, diagnostics = False):
    '''
    Run the heat stress calculations and add the new attributes to the station

    :param obj station: station object
    :param file logfile: logfile to store outputs
    :param boolean diagnostics: output diagnostic information
    :param boolean plots: make a plot

    :returns: updated station object with heat stress values.
    '''

    temperatures = utils.apply_flags_to_mask(station, "temperatures")
    rh = getattr(station, "relative_humidity") # no separate flags using fdi
    e_v = getattr(station, "vapour_pressure") # no separate flags using fdi
    windspeeds = utils.apply_flags_to_mask(station, "windspeeds")
   
    thi = utils.set_MetVar_attributes("temperature_humidity_index", "Temperature Humidity Index (THI)", "temperature_humidity_index", "1", temperatures.mdi, np.dtype('float64'))
    wbgt = utils.set_MetVar_attributes("wet_bulb_globe_temperature", "Wet Bulb Globe Temperature (WBGT)", "wet_bulb_globe_temperature", "C", temperatures.mdi, np.dtype('float64'))
    humidex = utils.set_MetVar_attributes("humidex", "Humidex", "humidex", "1", temperatures.mdi, np.dtype('float64'))
    apparent_t = utils.set_MetVar_attributes("apparent_temperature", "Apparent Temperature", "apparent_temperature", "C", temperatures.mdi, np.dtype('float64'))
    heat_index = utils.set_MetVar_attributes("heat_index", "Heat Index", "heat_index", "C", temperatures.mdi, np.dtype('float64'))



    thi.data = calculate_thi(temperatures.data, rh.data)
    wbgt.data = calculate_wbgt(temperatures.data, e_v.data)
    humidex.data = calculate_humidex(temperatures.data, e_v.data)
    apparent_t.data = calculate_apparent_t(temperatures.data, e_v.data, windspeeds.data)
    heat_index.data = calculate_heat_index(temperatures.data, rh.data)

    if plots or diagnostics:
        print "Heat stress variables calculated, setting attributes\n"
    else:
        logfile.write("Heat stress variables calculated, setting attributes\n")

    setattr(station, "THI", thi)
    setattr(station, "WBGT", wbgt)
    setattr(station, "humidex", humidex)
    setattr(station, "apparent_t", apparent_t)
    setattr(station, "heat_index", heat_index)

    station = utils.append_history(station, "Heat Stress Calculations")

    return station # run_calcs


def do_plots():
    '''
    Run plots of T-q or T-ev for the different heat-stress measures

    '''

    import matplotlib.pyplot as plt

    plt.clf()

    Ts = np.arange(-20,60)
    RHs = np.arange(0,100)

    meshT, meshRH = np.meshgrid(Ts, RHs)

    THI = calculate_thi(meshT, meshRH)
    levels = np.arange(0,200,10)
    CS = plt.contour(Ts, RHs, THI, levels, colors='k', linestyles = 'solid')
    plt.clabel(CS, levels[1::2], fmt='%1.0f', inline=1, fontsize=10)

    HI = calculate_heat_index(meshT, meshRH)
    HI = np.ma.masked_where(HI < 0, HI)
    levels = np.arange(0,400,20)
    CS = plt.contour(Ts, RHs, HI, levels, colors='r', linestyles = 'solid')
    plt.clabel(CS, levels[1::2], fmt='%1.0f', inline=1, fontsize=10)

    plt.xlabel("Temperature")
    plt.ylabel("RH")

    plt.show()

    plt.clf()

    Ts = np.arange(-20,60)
    evs = np.arange(0,100)

    meshT, meshev = np.meshgrid(Ts, evs)

    WBGT = calculate_wbgt(meshT, meshev)
    levels = np.arange(0,200,10)
    CS = plt.contour(Ts, evs, WBGT, levels, colors='k', linestyles = 'solid')
    plt.clabel(CS, levels[1::2], fmt='%1.0f', inline=1, fontsize=10)

    humidex = calculate_humidex(meshT, meshev)
    levels = np.arange(0,200,10)
    CS = plt.contour(Ts, evs, humidex, levels, colors='r', linestyles = 'dashed')
    plt.clabel(CS, levels[1::2], fmt='%1.0f', inline=1, fontsize=10)

    Ta = calculate_apparent_t(meshT, meshev, 0)
    levels = np.arange(0,200,10)
    CS = plt.contour(Ts, evs, Ta, levels, colors='b', linestyles = 'dotted')
    plt.clabel(CS, levels[1::2], fmt='%1.0f', inline=1, fontsize=10)

    plt.xlabel("Temperature")
    plt.ylabel("ev")
    plt.show()

if __name__ == "__main__":
    do_plots()

