#!/usr/local/sci/bin/python
#*****************************
#
# controller for internal QC checks.
#
#
#************************************************************************
#                    SVN Info
#$Rev:: 219                                           $:  Revision of last commit
#$Author:: rdunn                                      $:  Author of last commit
#$Date:: 2019-05-20 16:56:47 +0100 (Mon, 20 May 2019) $:  Date of last commit
#************************************************************************

'''
internal_checks.py invoked by typing::

  python2.7 internal_checks.py --restart_id 000000-99999 --end_id 999999-99999 --all [--doMonth] [--plots] [--diagnostics]

Input arguments:

--restart_id        First station to process

--end_id            Last station to process

--all               Run all tests

--doMonth           [False] Adjust QC tests to cope with appended months

--plots             [False] Create plots into the image directory

--diagnostics       [False] Verbose output
'''

#
import numpy as np
import scipy as sp
import os
import sys
import datetime as dt
import subprocess
import time
import gc

# RJHD utilities
import netcdf_procs as ncdfp
import qc_utils as utils
import qc_tests
from set_paths_and_vars import *

#*********************************************
def internal_checks(restart_id = "", end_id = "",
                    all_checks = True,
                    duplicate = False,
                    odd = False,
                    frequent = False,
                    diurnal = False,
                    gap = False,
                    records = False,
                    streaks = False,
                    climatological = False,
                    spike = False,
                    humidity = False,
                    cloud = False,
                    variance = False, 
                    winds = False, 
                    pressure = False,
                    precipitation = False,
                    diagnostics = False,
                    plots = False,
                    doMonth = False):
    '''
    Run through internal checks on list of stations passed
    
    :param str restart_id: which station to start on
    :param str end_id: which station to end on

    :param bool all_checks: run all the checks

    :param bool duplicate/odd/frequent/diurnal/gap/records/streaks/climatological/spike/humidity/cloud/variance/winds/pressure/precipitation: run each test separately
    :param bool diagnostics: print extra material to screen
    :param bool plots: create plots from each test [many files if all stations/all tests]
    :param bool doMonth: a monthly append process

    '''

    if all_checks:
        duplicate = True
        odd = True
        frequent = True
        diurnal = True
        gap = True
        records = True
        streaks = True
        climatological = True
        spike = True
        humidity = True
        cloud = True
        variance = True
        winds = True
        pressure = True
        precipitation = True
    else:
        print "single tests selected"
        
#    qc_code_version = subprocess.check_output(['svnversion']).strip()
    qc_code_version = subprocess.check_output(['svn', 'info', 'file:///home/h05/rdunn/svn/hadisd_py_qc/branches/monthly/'])
    for line in qc_code_version.split("\n"):
        if line.split(":")[0] == "Revision":
            qc_code_version = line.split(":")[1]
            break

        
    # get station information
    try:
        station_info = np.genfromtxt(os.path.join(INPUT_FILE_LOCS, STATION_LIST), dtype=(str))
    except IOError:
        print "station list not found"
        sys.exit()

    # sort truncated run
    startindex = [0]
    if restart_id != "":
        startindex, = np.where(station_info[:,0] == restart_id)


    if end_id != "":
        endindex, = np.where(station_info[:,0] == end_id)
        if endindex != len(station_info) -1:
            station_info = station_info[startindex[0]: endindex[0]+1]
        else:
            station_info = station_info[startindex[0]:]
    else:
        station_info = station_info[startindex[0]:]
        

    for st,stat in enumerate(station_info):     

        # if st%100 != 0: continue # do every nth station
  
        print dt.datetime.strftime(dt.datetime.now(), "%A, %d %B %Y, %H:%M:%S")
        print "{:35s} {:d}/{:d}".format("Station Number : ", st + 1, len(station_info))
        print "{:35s} {}".format("Station Identifier :", stat[0])
        if doMonth: print "Running with incomplete final year"

        # set up the log file
        logfile = file(LOG_OUTFILE_LOCS+stat[0]+'.log','w')
        logfile.write(dt.datetime.strftime(dt.datetime.now(), "%A, %d %B %Y, %H:%M:%S\n"))
        logfile.write("Internal Checks\n")
        logfile.write("{:35s} {}\n".format("Station Identifier :", stat[0]))

        process_start_time = time.time()

        station = utils.Station(stat[0], float(stat[1]), float(stat[2]), float(stat[3]))

        # latitude and longitude check
        if np.abs(station.lat) > 90.:
            if plots or diagnostics:
                print "{} {} {} {} {} {} {}\n".format(\
                        station.id,"Latitude Check",DATASTART.year, DATAEND.year,"All", "Unphysical latitude {}".format(station.lat))
            else:
                logfile.write("{} {} {} {} {} {} {}\n".format(\
                        station.id,"Latitude Check",DATASTART.year, DATAEND.year,"All", "Unphysical latitude {}".format(station.lat)))
                logfile.close()

            continue

        # check if station longitude outside of bounds
        if np.abs(station.lon) > 180.:       
            if plots or diagnostics:
                print "{} {} {} {} {} {} {}\n".format(\
                    station.id,"Longitude Check",DATASTART.year, DATAEND.year,"All", "Unphysical longitude {}".format(station.lon))
            else:
                logfile.write("{} {} {} {} {} {} {}\n".format(\
                        station.id,"Longitude Check",DATASTART.year, DATAEND.year,"All", "Unphysical longitude {}".format(station.lon)))
                logfile.close()
            continue

        # check if file is zipped
        if os.path.exists(os.path.join(NETCDF_DATA_LOCS, "hadisd.{}_19310101-{}_{}_raw.nc.gz".format(LONG_VERSION, END_TIME, station.id))):
            # if gzip file, unzip here
            subprocess.call(["gunzip",os.path.join(NETCDF_DATA_LOCS, "hadisd.{}_19310101-{}_{}_raw.nc.gz".format(LONG_VERSION, END_TIME, station.id))])
            time.sleep(5) # make sure it is unzipped before proceeding

        # read in the data
        ncdfp.read(os.path.join(NETCDF_DATA_LOCS, "hadisd.{}_19310101-{}_{}_raw.nc".format(LONG_VERSION, END_TIME, station.id)), station, process_vars, opt_var_list = carry_thru_vars, diagnostics = diagnostics)

        if plots or diagnostics:
            print "{:35s}  {}\n".format("Total station record size :",len(station.time.data))
        else:
            logfile.write("{:35s}  {}\n".format("Total station record size :",len(station.time.data)))

        match_to_compress = utils.create_fulltimes(station, process_vars, DATASTART, DATAEND, carry_thru_vars)

        station.qc_flags = np.zeros([len(station.time.data),71]) # changed to include updated wind tests, station level pressure & precipitation

        # get reporting accuracies and frequencies.

        for var in process_vars:

            st_var = getattr(station, var)
            st_var.reporting_stats = utils.monthly_reporting_statistics(st_var, DATASTART, DATAEND)


        # Add history text to netcdf file
        # Reporting Changes - TODO

        # Duplicate months - check on temperature ONLY
        if duplicate:
            # no change as result of incomplete year
            qc_tests.duplicate_months.dmc(station, ['temperatures'], process_vars, [0], DATASTART, DATAEND, logfile, diagnostics = diagnostics, plots = plots)

        # Odd Clusters
        if odd:
            # no change as result of incomplete year
            qc_tests.odd_cluster.occ(station,['temperatures','dewpoints','windspeeds','slp'], [54,55,56,57], DATASTART, logfile, diagnostics = diagnostics, plots = plots)
            utils.apply_flags_from_A_to_B(station, "slp", "stnlp", diagnostics = diagnostics)
            utils.apply_flags_from_A_to_B(station, "windspeeds", "winddirs", diagnostics = diagnostics)

        # Frequent Values
        if frequent:
            qc_tests.frequent_values.fvc(station, ['temperatures', 'dewpoints','slp'], [1,2,3], DATASTART, DATAEND, logfile, diagnostics = diagnostics, plots = plots, doMonth = doMonth)
            utils.apply_flags_from_A_to_B(station, "slp", "stnlp", diagnostics = diagnostics)

        # Diurnal Cycle 
        if diurnal:
            if np.abs(station.lat) <= 60.:
                qc_tests.diurnal_cycle.dcc(station, ['temperatures'], process_vars, [4], DATASTART, DATAEND, logfile, diagnostics = diagnostics, plots = plots, doMonth = doMonth)
                
            else:
                if plots or diagnostics:
                    print "Diurnal Cycle Check not run as station latitude ({}) > 60\n".format(station.lat)
                else:
                    logfile.write("Diurnal Cycle Check not run as station latitude ({}) > 60\n".format(station.lat))

        # Distributional Gap
        if gap:
            qc_tests.distributional_gap.dgc(station, ['temperatures','dewpoints','slp'], [5,6,7], DATASTART, DATAEND, logfile, diagnostics = diagnostics, plots = plots, GH = True, doMonth = doMonth)
            utils.apply_flags_from_A_to_B(station, "slp", "stnlp", diagnostics = diagnostics)

        # Records 
        if records:
            qc_tests.records.krc(station, ['temperatures','dewpoints','windspeeds','slp'], [8,9,10,11], logfile, diagnostics = diagnostics, plots = plots)
            utils.apply_flags_from_A_to_B(station, "windspeeds", "winddirs", diagnostics = diagnostics)
            utils.apply_flags_from_A_to_B(station, "slp", "stnlp", diagnostics = diagnostics)

        # Streaks and Repetitions 
        if streaks:
            qc_tests.streaks.rsc(station, ['temperatures','dewpoints','windspeeds','slp','winddirs'], [[12,16,20],[13,17,21],[14,18,22],[15,19,23],[66,67,68]], DATASTART, DATAEND, logfile, diagnostics = diagnostics, plots = plots, doMonth = doMonth)
            utils.apply_flags_from_A_to_B(station, "windspeeds", "winddirs", diagnostics = diagnostics)
            utils.apply_flags_from_A_to_B(station, "slp", "stnlp", diagnostics = diagnostics)

        # Climatological Outlier
        if climatological:
            qc_tests.climatological.coc(station, ['temperatures','dewpoints'], [24,25], DATASTART, DATAEND, logfile, diagnostics = diagnostics, plots = plots, doMonth = doMonth)
            # column 26 kept spare for slp

        # Spike
        if spike:
            qc_tests.spike.sc(station, ['temperatures','dewpoints','slp','windspeeds'], [27,28,29,65], DATASTART, DATAEND, logfile, diagnostics = diagnostics, plots = plots, doMonth = doMonth)
            utils.apply_flags_from_A_to_B(station, "windspeeds", "winddirs", diagnostics = diagnostics)
            utils.apply_flags_from_A_to_B(station, "slp", "stnlp", diagnostics = diagnostics)

        # Humidity cross checks
        if humidity:
            qc_tests.humidity.hcc(station, [30,31,32], DATASTART, DATAEND, logfile, diagnostics = diagnostics, plots = plots)

        # Cloud cross check
        if cloud:
            qc_tests.clouds.ccc(station, [33,34,35,36,37,38,39,40], logfile, diagnostics = diagnostics, plots = plots)

        # Variance
        if variance:
            qc_tests.variance.evc(station, ['temperatures','dewpoints','slp','windspeeds'], [58,59,60,61], DATASTART, DATAEND, logfile, diagnostics = diagnostics, plots = plots, doMonth = doMonth) 
            utils.apply_flags_from_A_to_B(station, "windspeeds", "winddirs", diagnostics = diagnostics)
            utils.apply_flags_from_A_to_B(station, "slp", "stnlp", diagnostics = diagnostics)

        # Winds
        if winds:
            qc_tests.winds.wdc(station, [62,63,64], DATASTART, DATAEND, logfile, diagnostics = diagnostics, plots = plots, doMonth = doMonth)

        # Pressure
        if pressure:
            qc_tests.pressure.spc(station, [69], DATASTART, DATAEND, logfile, diagnostics = diagnostics, plots = plots, doMonth = doMonth)

        # Precipitation
        if precipitation:
            qc_tests.precipitation.pcc(station, [70], DATASTART, DATAEND, logfile, diagnostics = diagnostics, plots = plots)


        # are flags actually applied?
        sys.stdout.flush()
        if diagnostics or plots: raw_input("stop")

        # write to file
        ncdfp.write(os.path.join(NETCDF_DATA_LOCS, "hadisd.{}_19310101-{}_{}_internal.nc".format(LONG_VERSION, END_TIME, station.id)), station, process_vars, os.path.join(INPUT_FILE_LOCS,'attributes.dat'), opt_var_list = carry_thru_vars, compressed = match_to_compress, processing_date = dt.datetime.strftime(dt.datetime.now(), "%d-%b-%Y"), qc_code_version = qc_code_version)
        # gzip the raw file
        subprocess.call(["gzip",os.path.join(NETCDF_DATA_LOCS, "hadisd.{}_19310101-{}_{}_raw.nc".format(LONG_VERSION, END_TIME, station.id))])


        logfile.write(dt.datetime.strftime(dt.datetime.now(), "%A, %d %B %Y, %H:%M:%S\n"))
        logfile.write("processing took {:4.0f}s\n\n".format(time.time() - process_start_time))
        logfile.close()

        # clean up
        gc.collect()

    print "Internal Checks completed\n"

    return # internal_checks

#************************************************************************
if __name__=="__main__":

    import argparse

    # set up keyword arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--restart_id', dest='restart_id', action='store', default = "",
                        help='Restart ID for truncated run, default = ""')
    parser.add_argument('--end_id', dest='end_id', action='store', default = "",
                        help='End ID for truncated run, default = ""')
    parser.add_argument('--diagnostics', dest='diagnostics', action='store_true', default = False,
                        help='Run diagnostics (will not write out file)')
    parser.add_argument('--plots', dest='plots', action='store_true', default = False,
                        help='Run plots (will not write out file)')
    parser.add_argument('--all', dest='all', action='store_true', default = False,
                        help='Run all checks')
    parser.add_argument('--duplicate', dest='duplicate', action='store_true', default = False,
                        help='Run duplicate months check')
    parser.add_argument('--odd', dest='odd', action='store_true', default = False,
                        help='Run odd cluster check')
    parser.add_argument('--frequent', dest='frequent', action='store_true', default = False,
                        help='Run frequent value check')
    parser.add_argument('--diurnal', dest='diurnal', action='store_true', default = False,
                        help='Run diurnal cycle check')
    parser.add_argument('--gap', dest='gap', action='store_true', default = False,
                        help='Run distributional gap check')
    parser.add_argument('--records', dest='records', action='store_true', default = False,
                        help='Run world records check')
    parser.add_argument('--streaks', dest='streaks', action='store_true', default = False,
                        help='Run streak check')
    parser.add_argument('--climatological', dest='climatological', action='store_true', default = False,
                        help='Run climatological outlier check')
    parser.add_argument('--spike', dest='spike', action='store_true', default = False,
                        help='Run spike check')
    parser.add_argument('--humidity', dest='humidity', action='store_true', default = False,
                        help='Run humidity cross checks')
    parser.add_argument('--cloud', dest='cloud', action='store_true', default = False,
                        help='Run cloud cross checks')
    parser.add_argument('--variance', dest='variance', action='store_true', default = False,
                        help='Run variance check')
    parser.add_argument('--winds', dest='winds', action='store_true', default = False,
                        help='Run winds checks')
    parser.add_argument('--pressure', dest='pressure', action='store_true', default = False,
                        help='Run pressure check')
    parser.add_argument('--precipitation', dest='precipitation', action='store_true', default = False,
                        help='Run precipitation cross check')
    parser.add_argument('--doMonth', dest='doMonth', action='store_true', default = False,
                        help='Run monthly update rather than annual one')
    args = parser.parse_args()

    if args.all:
        # check that no other test is set on top
        if args.duplicate: sys.exit("all tests and single test set - what did you want to do?")
        if args.odd: sys.exit("all tests and single test set - what did you want to do?")
        if args.frequent: sys.exit("all tests and single test set - what did you want to do?")
        if args.diurnal: sys.exit("all tests and single test set - what did you want to do?")
        if args.gap: sys.exit("all tests and single test set - what did you want to do?")
        if args.records: sys.exit("all tests and single test set - what did you want to do?")
        if args.streaks: sys.exit("all tests and single test set - what did you want to do?")
        if args.climatological: sys.exit("all tests and single test set - what did you want to do?")
        if args.spike: sys.exit("all tests and single test set - what did you want to do?")
        if args.humidity: sys.exit("all tests and single test set - what did you want to do?")
        if args.cloud: sys.exit("all tests and single test set - what did you want to do?")
        if args.variance: sys.exit("all tests and single test set - what did you want to do?")
        if args.winds: sys.exit("all tests and single test set - what did you want to do?")
        if args.pressure: sys.exit("all tests and single test set - what did you want to do?")
        if args.precipitation: sys.exit("all tests and single test set - what did you want to do?")

    internal_checks(restart_id = args.restart_id, 
                    end_id = args.end_id, 
                    all_checks = args.all,
                    duplicate = args.duplicate,
                    odd = args.odd,
                    frequent = args.frequent,
                    diurnal = args.diurnal,
                    gap = args.gap,
                    records = args.records,
                    streaks = args.streaks,
                    climatological = args.climatological,
                    spike = args.spike,
                    humidity = args.humidity,
                    cloud = args.cloud,
                    variance = args.variance,
                    winds = args.winds, 
                    pressure = args.pressure,
                    precipitation = args.precipitation,
                    diagnostics = args.diagnostics,
                    plots = args.plots,
                    doMonth = args.doMonth)


#************************************************************************
