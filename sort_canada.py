#!/usr/local/sci/bin/python2.7

#-------------------------------------------------
#
# Code to try and split out the different ways the 
#  Canadian stations present themselves in the list
#  provided by EC.  Not all possible options are 
#  gone through - leaving around 10 stations with 
#  undetermined issues.
#
# 30 May 2013, RJHD, Exeter
#
#-------------------------------------------------
#                START
#-------------------------------------------------

import numpy as np
import datetime as dt
import struct
import os
import re

# RJHD utilities
from station_selection import *

FILE_LOCS="input_files/"

#*********************************************
def writefile(filename, lines):
    '''
    Write the lines string to the filename file
    '''

    try:
        with open(filename,"a") as outfile:
            for line in lines:
                outfile.write(line)
            
    except IOError:
        print "cannot find %s" % filename

    return # writefile


#*********************************************
def name_test(name1,name2):
    '''
    Test if the strings in names1 and names2 are plausibly the same.
    '''

    delimiters=[" ",". ",".",'\\']
    regexPattern = '|'.join(map(re.escape, delimiters))

    name1_split=re.split(regexPattern,name1)
    name2_split=re.split(regexPattern,name2)

    suffixes=['A','UA','AWOS','AUT','(AUT)']
    
    if name1_split[-1] in suffixes:
        del name1_split[-1]
    if name2_split[-1] in suffixes:
        del name2_split[-1]

    result=False
    # test if first word the same
    if name1_split[0] == name2_split[0]:
        if len(name1_split)==1 or len(name2_split)==1:
            # if only one, then has to be correct
            result=True
        else:
            if name1_split[1] == name2_split[1]:
                # if second also the same then can be confident
                result=True
                
    return result # name_test

#****************************************************          
def process_EC_file():

    # set up the columns to read in 

    fieldwidths = (9,35,5,11,11,11,20)
    fmtstring = ''.join('%ds' % f for f in fieldwidths)
    parse = struct.Struct(fmtstring).unpack_from

    #print 'fmtstring: {!r}, record size: {} chars'.format(fmtstring,
    #                                             struct.calcsize(fmtstring))

    # storage files

    SINGLE=FILE_LOCS+"Canada_single.dat"
    ONOFF=FILE_LOCS+"Canada_onoff.dat"
    REMAINING=FILE_LOCS+"Canada_rem.dat"
    HOMOGENISATION=FILE_LOCS+"Canada_homogenisation.dat"
    DATES=FILE_LOCS+"Canada_dates.dat"
    GOODMOVE=FILE_LOCS+"Canada_goodmove.dat"
    QUESTIONABLEMOVE=FILE_LOCS+"Canada_questionablemove.dat"
    OVERLAP=FILE_LOCS+"Canada_overlap.dat"

    station_count = np.zeros(8)

    for outfile in [SINGLE, ONOFF, REMAINING, HOMOGENISATION, DATES, GOODMOVE, QUESTIONABLEMOVE, OVERLAP]:
        try:
            # remove files in each run
            os.remove(outfile)
        except OSError:
            print "file does not exist ", outfile


    stn_id=[]
    location=[]
    state=[]
    lat=[]
    lon=[]
    active=[]
    date=[]
    datestr=[]
    lines=[]

    # allow option of no EC file available - as may not be able to share this
    try:

        with open(FILE_LOCS+"WMOHistory2014.txt",'r') as infile:

            for l,line in enumerate(infile):

                try:

                    if line.strip() == "":
                        continue

                    # extract the details
                    fields = parse(line)

                    if len(stn_id) > 0:

                        if int(fields[0].strip()) != stn_id[-1]:
                            # have change in station ID - process and reset

                            if len(stn_id)==1:
                                # nothing wrong with this one - so print to good listing
                                writefile(SINGLE,lines)
                                station_count[0] += 1

                            else:

                                sort_order = np.array(date).argsort()

                                stn_id = np.array(stn_id)[sort_order]
                                location = np.array(location)[sort_order]
                                state = np.array(state)[sort_order]
                                lat = np.array(lat)[sort_order]
                                lon = np.array(lon)[sort_order]
                                active = np.array(active)[sort_order]
                                date = np.array(date)[sort_order]
                                lines = np.array(lines)[sort_order]

                                names=False
                                # test the names for consistency
                                if len(stn_id)>=2:
                                    for l,loc in enumerate(location[:-1]):
                                        names=name_test(location[l],location[l+1])


                                if len(set(date)) != len(date):
                                    # duplicated dates

                                    writefile(DATES,lines)
                                    station_count[4] += 1

                                elif len(stn_id) == 2:

                                    if (active[0]=="Active" and active[1]=="Inactive") or (active[1]=="Active" and active[0]=="Inactive") and names:
                                        # 2 entries, names match, but have active and inactive settings
                                        writefile(ONOFF, lines)
                                        station_count[1] += 1

                                    elif active[0]==active[1] and names:
                                        # 2 entries, names match, but two identical statuses
                                        writefile(HOMOGENISATION, lines)
                                        station_count[3] += 1

                                    elif not names:
                                        # 2 entries, no match - could just be a name change
                                        writefile(QUESTIONABLEMOVE, lines)
                                        station_count[6] += 1
                                    else:
                                        raw_input("missing a 2 entry solution")

                                elif (len(stn_id) == 3) and \
                                        (active[0]==active[2]) and (active[1]=='Inactive') and name_test(location[0],location[1]):
                                    # 3 entries, with 2 locations and correct active/inactive order

                                    writefile(GOODMOVE, lines)
                                    station_count[5] += 1

                                elif (len(stn_id) == 3) and \
                                        (active[0]==active[2]) and (active[1]=='Active') and name_test(location[2],location[1]):
                                    # 3 entries, with 2 locations and correct active/inactive order

                                    writefile(GOODMOVE, lines)
                                    station_count[5] += 1

                                elif (len(stn_id) == 3) and names and len(set(lat))==1 and len(set(lon))==1:
                                    # 3 entries, all the same location - as per lat and long as well as names

                                    writefile(HOMOGENISATION, lines)
                                    station_count[6] += 1

                                elif (len(stn_id) == 3) and \
                                        (active[0]!=active[2]) and name_test(location[0],location[2]) and \
                                        (active[0]==active[1]) and not name_test(location[1],location[2]) :
                                    # 3 entries, new location/name active before old one inactive

                                    writefile(OVERLAP, lines)
                                    station_count[7] += 1

                                elif (len(stn_id) == 4) and name_test(location[0],location[1]) and name_test(location[2],location[3]) and \
                                       not name_test(location[1],location[2]) and \
                                        (active[0]==active[2]) and (active[1]==active[3]):
                                    # 4 entries, with 2 locations and correct active/inactive order

                                    writefile(GOODMOVE, lines)
                                    station_count[5] += 1

                                elif (len(stn_id) == 4) and \
                                        (active[0]==active[1]) and (active[2]==active[3]) and (active[0]!=active[2]) and names:
                                    # 4 entries, active/inactive order incorrect

                                    writefile(OVERLAP, lines)
                                    station_count[7] += 1

                                elif (len(stn_id) == 4) and name_test(location[0],location[1]) and name_test(location[2],location[3]) and \
                                       not name_test(location[1],location[2]) and \
                                        (active[0]!=active[1]) and (active[2]==active[3]) :
                                    # 4 entries, with 2 locations and correct active/inactive order for the move, and homogenisation for first location

                                    writefile(GOODMOVE, lines)
                                    station_count[5] += 1

                                elif (len(stn_id) == 4) and name_test(location[0],location[1]) and name_test(location[1],location[2]) and \
                                       not name_test(location[2],location[3]) and \
                                        (active[0]==active[1]) and (active[2]!=active[3]) :
                                    # 4 entries, with 2 locations and correct active/inactive order for the move, and homogenisation for second location

                                    writefile(GOODMOVE, lines)
                                    station_count[5] += 1

                                elif (len(stn_id) == 4) and names and len(set(lat))==1 and len(set(lon))==1:
                                    # 4 entries, all in same location - as per lat and long as well as names

                                    writefile(HOMOGENISATION, lines)
                                    station_count[3] += 1

                                elif (len(stn_id) == 4) and names :
                                    # 4 entries, names correct, but some station moves 

                                    writefile(QUESTIONABLEMOVE, lines)
                                    station_count[6] += 1

                                elif (len(stn_id) == 5) and \
                                        name_test(location[0],location[1]) and \
                                        name_test(location[2],location[3]) and \
                                        not name_test(location[0],location[2]) and \
                                        not name_test(location[4],location[2]) and \
                                        (active[0]==active[2]==active[4]) and (active[1]==active[3]=='Inactive'):
                                    # 5 entries, with 3 locations and correct active/inactive order

                                    writefile(GOODMOVE, lines)
                                    station_count[5] += 1

                                elif (len(stn_id) == 5) and \
                                        name_test(location[0],location[1]) and \
                                        name_test(location[2],location[3]) and \
                                        name_test(location[0],location[2]) and \
                                        not name_test(location[4],location[2]) and \
                                        (active[0]==active[2]==active[4]) and (active[1]==active[3]=='Inactive'):
                                    # 5 entries, with 2 locations and correct active/inactive order

                                    writefile(GOODMOVE, lines)
                                    station_count[5] += 1

                                elif (len(stn_id) == 5) and \
                                        (active[0]==active[1]) and (active[2]==active[3]) and (active[0]!=active[2]):
                                    # 5 entries, active/inactive order incorrect

                                    writefile(OVERLAP, lines)
                                    station_count[7] += 1

                                elif (len(stn_id) == 5) and names and len(set(lat))==1 and len(set(lon))==1:
                                    # 5 entries, names, latitudes and longitudes all correct - assume multiple statuses

                                    writefile(HOMOGENISATION, lines)
                                    station_count[3] += 1

                                elif (len(stn_id) == 5) and names :
                                    # 5 entries, names correct, but some station moves 

                                    writefile(QUESTIONABLEMOVE, lines)
                                    station_count[3] += 1

                                elif (len(stn_id) == 6) and \
                                        name_test(location[0],location[1]) and \
                                        name_test(location[2],location[3]) and \
                                        name_test(location[4],location[5]) and \
                                        not name_test(location[0],location[2]) and \
                                        not name_test(location[2],location[4]) and \
                                        (active[0]==active[2]==active[4]) and (active[1]==active[3]==active[5]):
                                    # 6 entries, with 3 locations and correct active/inactive order

                                    writefile(GOODMOVE, lines)
                                    station_count[5] += 1

                                elif (len(stn_id) == 6) and names and len(set(lat))==1 and len(set(lon))==1:
                                    # 6 entries, names, latitudes and longitudes all correct - assume multiple statuses

                                    writefile(HOMOGENISATION, lines)
                                    station_count[3] += 1

                                elif (len(stn_id) == 6) and names :
                                    # 6 entries, names correct, but some station moves

                                    writefile(QUESTIONABLEMOVE, lines)
                                    station_count[6] += 1

                                else:
                                    writefile(REMAINING,lines)
                                    station_count[2] += 1



                            # and reset the lists
                            stn_id=[]
                            location=[]
                            state=[]
                            lat=[]
                            lon=[]
                            active=[]
                            date=[]
                            datestr=[]
                            lines=[]
                        # change in ID loop

                    # at least one station read in loop

                    # read in the data into the list
                    stn_id+=[int(fields[0].strip())]
                    location+=[fields[1].strip()]
                    state+=[fields[2].strip()]
                    lat+=[float(fields[3].strip())]
                    lon+=[float(fields[4].strip())]
                    active+=[fields[5].strip()]
                    date+=[dt.datetime.strptime(fields[6].strip(),"%Y-%m-%d %H:%M:%S")]
                    datestr+=[fields[6].strip()]
                    lines+=[line]

                    #print stn_id, location,state, lat,lon, active, date

                except ValueError:
                    # skip non-numeric lines - as these are not WMO ids
                    continue

    except IOError:
        print "No EC file available"

    print "single, onoff, remaining, homog, dates, good move, quest. move, overlap"
    print station_count
    return # process_EC_file

#*********************************************
def read_canada_info(filename, DATA_START):
    '''
    Read station information from processed Canadian files

    :param str filename: filename to read
    :returns:
       ids, names, lats, lons, active and date arrays
    '''

    sorted_station_info = np.genfromtxt(FILE_LOCS + filename,delimiter=[9,35,5,11,11,11,20],dtype=(str))

    ids = sorted_station_info[:,0].astype(int)
    names = sorted_station_info[:,1].astype(str)
    lats = sorted_station_info[:,3].astype(float)
    lons = sorted_station_info[:,4].astype(float)
    active = sorted_station_info[:,5].astype(str)
    date = np.array([dt.datetime.strptime(d.strip(), "%Y-%m-%d %H:%M:%S") for d in sorted_station_info[:,6].astype(str)])
    
    # fix 1777 start date which can't be handled by datetime
    newdates = []
    for d in date:
        if d.year == 1777:
            newdates += [dt.datetime(DATA_START,1,1,0,0)]
        else:
            newdates += [d]
        

    return ids, names, lats, lons, active, np.array(newdates) # read_canada_info

#*********************************************
def process_canadian_stations(all_stations, DATA_START):

    '''
    if in Canada_single - then single occurrance - USE
    if in Canada_onoff - then single occurrance with time limits - USE
    if in Canada_homogenisation - then multiple occurances in single location - USE
    if in Canada_goodmove - then well documented station move - USE - but only section which matches the dates
    else DON'T USE (_overlap, _questionablemove, _rem, _dates)
    '''

    station_ids = np.array([stn.id for stn in all_stations])
    # can only work on 71???0-99999 station numbers
    canadian_ids = np.array([s for s, stn in enumerate(station_ids) if stn[:2] == "71" and stn[-7:] == "0-99999"])

    use = np.array([0 for i in range(len(canadian_ids))])# 1 - use, 0 - not tested, -1 - don't use
    # only reject those we are sure about, and keep those that we can't test as we don't know

    # single, onoff & homogenisation
    for category_files in ["Canada_single.dat", "Canada_onoff.dat", "Canada_homogenisation.dat"]:

        test_ids, test_names, test_lats, test_lons, test_active, test_date = read_canada_info(category_files, DATA_START)

        for c,cid in enumerate(canadian_ids):
            id_to_compare = int(station_ids[cid][:5])

            if id_to_compare in test_ids:
                # ID present

                loc = np.where(test_ids == id_to_compare)[0][0]

                test_station = Station(test_ids[loc], test_lats[loc], test_lons[loc], all_stations[cid].elev) # fake the elev
                test_station.name = test_names[loc].strip()
                test_station.call = ""
                probs = do_match(test_station, all_stations[cid])

                if np.product(probs) > PROB_THRESHOLD:
                    use[c] = 1
                else:
                    print  all_stations[cid].name, all_stations[cid]
                    print  test_station.name.strip(), test_station
                    print probs, "\n"


    # overlap, questionablemove, rem, dates
    for category_files in ["Canada_rem.dat", "Canada_dates.dat", "Canada_questionablemove.dat", "Canada_overlap.dat"]:

        test_ids, test_names, test_lats, test_lons, test_active, test_date = read_canada_info(category_files, DATA_START)

        for c,cid in enumerate(canadian_ids):
            id_to_compare = int(station_ids[cid][:5])

            if id_to_compare in test_ids:
                # ID present - don't use this station
                use[c] = -1

                # restrict start and end times so won't be selected.
                all_stations[cid].start = dt.datetime.today()
                all_stations[cid].end = dt.datetime.today()

    # goodmove
    for category_files in ["Canada_goodmove.dat"]:

        outfilename = FILE_LOCS+"Canada_time_ranges.dat"
        try:
            os.remove(outfilename)
        except OSError:
            print "file does not exist ", outfilename
        outfile = file(outfilename, "w")
            

        test_ids, test_names, test_lats, test_lons, test_active, test_date = read_canada_info(category_files, DATA_START)

        for c,cid in enumerate(canadian_ids):
            id_to_compare = int(station_ids[cid][:5])

            if id_to_compare in test_ids:
                # ID present
                locs, = np.where(test_ids == id_to_compare)
                
                # test which locations match
                all_probs=[]

                for loc in locs:
                    test_station = Station(test_ids[loc], test_lats[loc], test_lons[loc], all_stations[cid].elev) # fake the elev
                    test_station.name = test_names[loc]
                    test_station.call = ""
                    probs = do_match(test_station, all_stations[cid])

                    all_probs += [np.product(probs)]
                
                good_locs, = np.where(np.array(all_probs) > PROB_THRESHOLD)

                # need to test which range of dates can be taken.
                start = dt.datetime(DATA_START,1,1,0,0)
                end   = dt.datetime(dt.datetime.now().year+1,1,1,0,0)

 
                # single "Active" and last entry - then use as start date
                if (len(good_locs) == 1) and (good_locs[0] + 1 == len(all_probs)) and (test_active[locs][good_locs[0]].strip() == "Active"):
                    use[c] = 1
                    start = test_date[locs][good_locs[0]]
                    
                    # adjust start date and write out useful range
                    all_stations[cid].start = start
                    outfile.write("{} {} {}\n".format(station_ids[cid], dt.datetime.strftime(start ,"%Y-%m-%d %H:%M:%S"), dt.datetime.strftime(end, "%Y-%m-%d %H:%M:%S")))

                elif len(good_locs) == 2:
                    
                    # "Active" followed by "Inactive" - use as range 71100
                    if (test_active[locs][good_locs[0]].strip() == "Active") and (test_active[locs][good_locs[1]].strip() == "Inactive"):
                        use[c] = 1
                        start = test_date[locs][good_locs[0]]
                        end   = test_date[locs][good_locs[1]]
                        
                        # adjust start and end date and write out useful range
                        all_stations[cid].start = start
                        all_stations[cid].end = end
                        outfile.write("{} {} {}\n".format(station_ids[cid], dt.datetime.strftime(start ,"%Y-%m-%d %H:%M:%S"), dt.datetime.strftime(end, "%Y-%m-%d %H:%M:%S")))


                    # 2 "Active"s, then if in final place, use first as start date 71038
                    elif (test_active[locs][good_locs[0]].strip() == "Active") and (test_active[locs][good_locs[1]].strip() == "Active") and (good_locs[0] + 1 == len(all_probs) - 1) and (good_locs[1] + 1 == len(all_probs)) :                 
                        use[c] = 1
                        start = test_date[locs][good_locs[0]]

                        # adjust start date and write out useful range
                        all_stations[cid].start = start
                        outfile.write("{} {} {}\n".format(station_ids[cid], dt.datetime.strftime(start ,"%Y-%m-%d %H:%M:%S"), dt.datetime.strftime(end, "%Y-%m-%d %H:%M:%S")))


                # 3 "Active"s or combination of "Active" and "Inactive" - if no other IDs present, then use all
                elif len(locs) == len(good_locs):
                    use[c] = 1
                    # no change to start/end times


        outfile.close()
                    

    print "{} Canadian stations processed - {} kept, {} not tested, {} rejected".format(len(use), len(use[use == 1]), len(use[use == 0]), len(use[use == -1]))
 
    return all_stations # process_canadian_stations

#*********************************************
if __name__ == "__main__":

    process_EC_file()

#-------------------------------------------------
#                END
#-------------------------------------------------
