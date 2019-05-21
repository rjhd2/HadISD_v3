#!/usr/local/sci/bin/python2.7
#------------------------------------------------------------
#                    SVN Info
#$Rev:: 219                                           $:  Revision of last commit
#$Author:: rdunn                                      $:  Author of last commit
#$Date:: 2019-05-20 16:56:47 +0100 (Mon, 20 May 2019) $:  Date of last commit
#------------------------------------------------------------
# History:
#  16 - Created and tested
#       27 April 2012 RJHD
#  32 - Sorted copying functionality - not working.
#       Added a failsafe doftp to NCDC if the python download fails (filesize issue)
#       7 January 2013 RJHD
#  33 - Forgot to add in the station725 directory special bit.  Stops after 10 attempts
#       in case remote files have changed on server!
#       14 January 2014 RJHD
#  
#   
#
#------------------------------------------------------------
# START
'''
get_isd_files.py invoked by typing::

  python2.7 get_isd_file.py --start 1931 --end 2018 --restart_id 000000-99999 --end_id 999999-99999 --get_dir_list

Re-run the FTP scrape, checking if files need updating, 
and only download those that do.

Input Arguments:

--start		[1931] first year downloaded

--end		[current] last year downloaded

--restart_id	first station ID to download

--end_id	last station ID to download

--get_dir_list	[True] get file listings of server.  This takes time so
		if run for a single/subset of stations cached
		versions may be useful

Rsync style of capability - downloads new, or copies old if no change, or leaves alone if
new files still most recent.  Deletes files locally if not present on server
'''


# Had issues with urlretrieve and file sizes and so have a doftp call
# using os.system as a failsafe backup (Jan2013)


import urllib2
import datetime
import os
import shutil
import sys
import subprocess
import argparse
import calendar

# RJHD utils
from set_paths_and_vars import *


# Globals
HOST='ftp://ftp.ncdc.noaa.gov'
ISD_LOC='/pub/data/noaa/'

STATIONLISTFILES=[STATION_LIST, MERGER_LIST]

        
#---------------------------------------------------------------------
class DataFile(object):
    """
    Class to hold filename, size and modification date

    creates target_dir attribute from filen value
    """
    def __init__(self,filen,date,size):
        self.filen=filen
        self.date=date
        self.size=size

        # set target_dir using filen
        if (self.filen[0:3] == '725') or (self.filen[0:3] == '726') or (self.filen[0:3] == '727') or (self.filen[0:3] == '728') or (self.filen[0:3] == '729'):
            self.target_dir=r'station725/'
        else:
            self.target_dir=r'station'+self.filen[0:2]+r's/'


    def __str__(self):
        if self.size/(1024.**2)>1:
            hsize=str(int(self.size/(1024.**2)))+"MB"
        elif self.size/(1024.)>1:
            hsize=str(int(self.size/(1024.)))+"kB"
        else:
            hsize=str(int(self.size))+"B"
      
        return "station %s, created on %s, with size %s" % \
            (self.filen,self.date.strftime("%Y %b %d, %H:%M"), hsize)

    __repr__=__str__


#---------------------------------------------------------------------

#*********************************************
def GetFileListing(host,location,year = ""):
    """Log in to Server and return file listing"""

    try:
        ftpsite=urllib2.urlopen(host+location+str(year))

        html_dump=ftpsite.readlines()
    except KeyboardInterrupt:
        raise UserWarning("Aborted listing")
    except:
        raise UserWarning("FTP file listing failure")

    # write out to file to save doing this for each station.
    if year != "":
        outfile = file("{}/isd_ftp_list_{}.txt".format(INPUT_FILE_LOCS, year), "w")
    else:
        outfile = file("{}/isd_ftp_list_parent.txt".format(INPUT_FILE_LOCS), "w")
        
    for line in html_dump:
        outfile.write(line)
    outfile.close()

    return # GetFileListing

#*********************************************
def ExtractFileData(file_list, directory = False):
    """
    Extract all the file/directory names, creation dates (and sizes)

    Convert the sizes to integers and the file modified 
    date to datetime objects

    Just parses HTML - possibly could be improved

    directory = sets if keeping only directory info
    """

    data_files=[]

    nlines = len(file_list)

    for table_start,line in enumerate(file_list):
        if line[:6] == "</pre>":
            break

    for row in range(table_start+2, nlines):

        # first line of footer
        if file_list[row][:6] == "</pre>":
            break
        
        # first line of data mixed into header
        if file_list[row][:5] == "<pre>":
            line = file_list[row].split("</h4>")[1]
        else:
            line = file_list[row]

        try:
            # two creation formats
            file_creation_date = datetime.datetime.strptime(line.split("<")[0][-13:-1] ,"%b %d %Y")    
        except:
            # if date in last 6 months then linux gives no year - so need to take care on year boundary
            month_abbr = line.split("<")[0][-13:-1].split(" ")[0]
            month = list(calendar.month_abbr).index(month_abbr)

            if month > dt.datetime.now().month:
                file_creation_date = datetime.datetime.strptime("{} {}".format(line.split("<")[0][-13:-1], dt.datetime.now().year-1) ,"%b %d %H:%M %Y")    
            else:
                file_creation_date = datetime.datetime.strptime("{} {}".format(line.split("<")[0][-13:-1], dt.datetime.now().year) ,"%b %d %H:%M %Y")    


        if directory:
            size = line[33:41]
        else:
            size = float(line[33:41])

        file_name = line.split(">")[1][:-3]       

        if directory:
            if size == "Directory":
                data_files.append(DataFile(file_name[:-1],file_creation_date,0))
        else:
            data_files.append(DataFile(file_name,file_creation_date,size))

    ########## OLD 2
    # for table_start,line in enumerate(file_list):
    #     if line[:7] == "<table>":
    #         break

    # for row in range(table_start+1,nlines,3):

    #     if file_list[row][:8] == "</tbody>":
    #         break

    #     file_creation_date = datetime.datetime.strptime(file_list[row][8:-6] ,"%a %b %d %H:%M:%S %Y")    

    #     if directory:
    #         size = file_list[row+1][4:-6]
    #     else:
    #         size = float(file_list[row+1][4:-6])

    #     file_name = file_list[row+2][file_list[row+2].find('"')+1:file_list[row+2].rfind('"')]       

    #     if directory:
    #         if size == "Directory":
    #             data_files.append(DataFile(file_name[:-1],file_creation_date,0))
    #     else:
    #         data_files.append(DataFile(file_name,file_creation_date,size))

    ########## OLD 1
    # for line in file_list:
    #     if line[0]=="<":
    #         # ignore all HTML tag lines
    #         continue
    #     else:
    #         time_string=line[0:17]
    #         size=int(line[18:30].strip())

    #         file_name=line[line.find('"')+1:line.rfind('"')]

    #         file_creation_date=datetime.datetime.strptime(time_string, '%b %d %Y %H:%M')

    #         data_files.append(DataFile(file_name,file_creation_date,size))
        
    return data_files # ExtractFileData


#*********************************************
def CombineHadISDStations(filelist, outfile, uk=False, restart_id="", end_id=""):
    """
    Combine all the input station lists into one
    master list
    """
    if restart_id == "":
        add_to_list = True
    else:
        add_to_list = False

    stationlist=[]
    # candidate_stations.txt
    try:
        with open(INPUT_FILE_LOCS+filelist[0],'r') as ifile:
            for line in ifile:
                station = line.split()[0].strip()

                if restart_id == station:
                    # if start set, then add, and all subsequent
                    add_to_list = True

                if add_to_list:
                    # add all intervening
                    stationlist+=[station]

                # separate test to terminate addition to list.
                if end_id == station:
                    # if end set, then add, but no subsequent
                    #    all previous set before try loop
                    add_to_list = False

    except IOError:
        print "Could not find station list file ", filelist[0]

    # final_mergers.txt
    try:
        with open(INPUT_FILE_LOCS+filelist[1],'r') as ifile:
            for line in ifile:
                station = line.split()[0].strip()
                line=line.split()

                if station in stationlist:
                    # if the primary source station present, then add, else not
                    stationlist.extend(line[1:])
                    # removes need for restart and end IDs as only add to list
                    #  if the primary station is present in the first place
                    #  which has been set by the loop above.                    

    except IOError:
        print "Could not find station list file ", filelist[1]
    except IndexError:
        print "Single file only", filelist

    stationlist=sorted(set(stationlist))

    try:
        with open(INPUT_FILE_LOCS+outfile,'w') as ofile:
            for station in stationlist:
                if uk:
                    if station[:2] == "03":
                        ofile.write("%s\n" % station)
                else:
                    ofile.write("%s\n" % station)

    except IOError:
        print "Could not create outfile, ", outfile

    return # CombineHadISDStations

#*********************************************
def HadISDStationExtract(mlistfile):
    """
    Extract the list of stations used in HadISD
    """

    try:
        with open(INPUT_FILE_LOCS + mlistfile,'r') as ifile:
            MasterStationList=[]
            for line in ifile:
                if line[0]=="#":
                    # escape clause for testing
                    break
                MasterStationList.append(line.strip().split()[0])

    except IOError:
        print "Could not find master list ", mlistfile
    except EOFError:
        print "End of file ", mlistfile

    return MasterStationList # HadISDStationExtract

#*********************************************
def FilterStations(MasterStationList,data_files):
    """
    Filter the files to download by those stations required

    If station of remote file exists in master list, want to download, 
    else not.    
    """

    DownloadStationList=[]
    for dfile in data_files:
        station=dfile.filen.split('.')[0][:-5]
        if station in MasterStationList:
            # We want this file
            DownloadStationList.append(dfile)

    return DownloadStationList # FilterStations

#*********************************************
def CheckToDelete(MasterStationList,data_files):
    """
    Extract file names on server, and check if any which are not in
    master station list

    These already may not exist on the local system
    """

    RemoteStations=[]

    for dfile in data_files:
        RemoteStations+=[dfile.filen.split('.')[0][:-5]]

    StationDeleteList=[]
    for station in MasterStationList:
        if station not in RemoteStations:
            StationDeleteList+=[station]
 
    return StationDeleteList # CheckToDelete 

#*********************************************
def DownloadChecks(file_path, datafile, lastrun, CheckSize = False):
    """
    Check whether this file should be re-downloaded

    Compares file date and size (where appropriate) to set
    download flag, which is returned.

    NOTE - remote files are gzip -9, local files are gzip -6 
    (the default) so remote files are smaller than local
    files that have been gunzip'd and gzip'd.  Therefore
    there is the CheckSize setting.
    """

    Status=False

    if os.path.exists(file_path):

        size=os.path.getsize(file_path)

        if lastrun < datafile.date:
            print "Different dates {} in {}".format(datafile.filen, "/".join(file_path.split("/")[:-1]))
            print "Local: {},  Server: {}".format(lastrun, datafile.date)
            # local date is before remote date
            Status=True
        elif CheckSize and size != datafile.size:
            # if checking file sizes & are different
            print "Different sizes {} in {}".format(datafile.filen, "/".join(file_path.split("/")[:-1]))
            print "Local: {},  Server: {}".format(size, datafile.size)
            Status=True

        else:
            # to ensure that all has been set
            print "Up to date local {} exists in {}".format(datafile.filen, "/".join(file_path.split("/")[:-1])) 
            Status=False
    else:
        print "No local version {} in {}".format(datafile.filen, "/".join(file_path.split("/")[:-1]))
        # no local version, download
        Status=True

    return Status # DownloadChecks

#*********************************************
def CheckToDownload(DownloadStationList, lastrun, copy_only):
    """
    Check if remote file is different from local file

    Check using file ages and sizes if remote file
    is different (newer/better) than local file.  

    NOTE - creation time not easily accessed in Unix and
    modification time will change each time the files are
    processed, so using UseFileAge setting can select whether
    to use a hard-coded date.
    NOTE - remote files are gzip -9, local files are gzip -6 
    (the default) so remote files are smaller than local
    files that have been gunzip'd a and gzip'd.  Therefore
    there is the CheckSize setting.
    """

    ToDownload=[False for dsl in range(len(DownloadStationList))]

    # if wanting to keep a previous version of the ISD files, then return all as Download = False
    if copy_only:
        return ToDownload

    else:

        for df, dfile in enumerate(DownloadStationList):

            old_file = OLD_ISD_DATA_LOCS + dfile.target_dir + dfile.filen

            ToDownload[df] = DownloadChecks(old_file, dfile, lastrun)
            # also check if file recently downloaded to same directory   
            #   and if this needs overwriting
            #   allows for scrape to be performed into same directory
            new_file = ISD_DATA_LOCS + dfile.target_dir + dfile.filen

            # only check if the file exists there.
            if os.path.exists(new_file):
                ToDownload[df] = DownloadChecks(new_file, dfile, lastrun,  CheckSize = False)

        return ToDownload # CheckToDownload

#*********************************************
def write_last_download(station, year):
    """
    (Over)Writes file containing last station download was attempted for a given year
    """

    outfile=os.path.join(ISD_DATA_LOCS, 'failed_download_{}.txt'.format(year))
    with open(outfile,'w') as lf:
        lf.write("{}".format(station))

    return
    

#*********************************************
def DownloadStations(DownloadStationList,ToDownload,year, verbose = False):
    """
    Download the required files if necessary for this year

    In filtered file list (DownloadStationList), if the files 
    have changed (ToDownload), then download them.  Check size
    after completed download to ensure that it completed 
    successfully.

    If downloading not necessary, files are copied from the
    old location to the new location.

    :param bool verbose: allows for verbose setting
    """
   
    # test if file left over from previous download
    fail_start_index = 0
    if os.path.exists(os.path.join(ISD_DATA_LOCS, 'failed_download_{}.txt'.format(year))):
        
        logfile=os.path.join(ISD_DATA_LOCS, 'failed_download_{}.txt'.format(year))
        with open(logfile,'r') as lf:
            for l, line in enumerate(lf):
                if l == 0:
                    fail_start_id = line.strip()
        
        # find location of station in DownloadStationList  
        station_id_list = [d.filen for d in DownloadStationList]
        try:
            fail_start_index = station_id_list.index(fail_start_id)
        except ValueError:
            # not in the list
            fail_start_index = 0

    # download the files
    for df,dfile in enumerate(DownloadStationList):

        # if fail_start_index set, then skip all until that station
        if df < fail_start_index:
            continue


        new_loc=ISD_DATA_LOCS + dfile.target_dir

        if not os.path.exists(new_loc):
            # if storage directory doesn't exist, make it
            os.mkdir(new_loc)

        if ToDownload[df]:

            if verbose:
                print "Downloading ", dfile.filen

            downloaded_size=0

            attempts=0
            while downloaded_size != dfile.size:

                # write logfile for each downloaded file
                write_last_download(dfile.filen, year)
                
                try:
                    print 'doftp -host '+HOST[6:]+' -cwd /pub/data/noaa/'+str(year)+'/ -get '+dfile.filen+'='+new_loc + dfile.filen
                    subprocess.check_call(["doftp",'-host', HOST[6:], '-cwd', '/pub/data/noaa/'+str(year)+'/', '-get', dfile.filen+'='+new_loc + dfile.filen])
                    
                # handle the error
                except subprocess.CalledProcessError:
                    print subprocess.CalledProcessError.message
                    print "waiting 10 sec and trying again"
                    import time
                    time.sleep(10) # wait 10 seconds and onto next while loop

 
                except OSError:
                    # executable not found
                    print "Issue with doftp"
                    raise OSError("doftp not found")

                try:
                    downloaded_size=os.path.getsize(new_loc + dfile.filen)
                except OSError:
                    # file doesn't exist
                    downloaded_size = 0

                attempts += 1

                if attempts > 10:
                    print "Cancelling download after 10 attempts as filesizes don't match"
                    raise Exception("Cancelling download after 10 attempts")

            if verbose:
                print "  Download Complete"
        else:
            if verbose:
                print "File %s not downloaded" % (dfile.filen)

            if os.path.exists(new_loc+dfile.filen):
                # file already exists
                # do not copy
                if verbose:
                    print "  Up-to-date file already exists in ", new_loc
            else:
                old_loc = OLD_ISD_DATA_LOCS + dfile.target_dir
                if verbose:
                    print "  Copying file from {} to {}".format(old_loc, new_loc)
            
                shutil.copy(old_loc + dfile.filen, new_loc + dfile.filen)

    return # DownloadStations

#*********************************************
def DeleteFiles(ToDelete, year, verbose = False):
    """
    Delete local files if they are not on the server

    If local files marked as not existing on server
    exist, then delete them.  They may already not exist
    locally - and no output printed in this case    
    """

    for st, station in enumerate(ToDelete):

        target_dir=r'station'+station[0:2]+r's/'

        new_loc=ISD_DATA_LOCS+target_dir

        st_filen=station+'-'+str(year)+'.gz'

        if os.path.exists(new_loc+st_filen):
            # if file exists
            os.remove(new_loc+st_filen)

            if verbose:
                print "Removed %s " % st_filen
                
    return # DeleteFiles

#*********************************************
def LastRunDate(infile):

    try:
        with open(infile,'r') as lf:
            for line in lf:
                date = datetime.datetime.strptime(line.split("at")[-1].strip(), '%d-%m-%Y %H:%M') # the last time this code was run prior to planned run
    except:
        print "Could not read log file {}".format(infile)
        date = DATASTART

    return date # LastRunDate


#*********************************************
def main(STARTYEAR, ENDYEAR, restart_id="", end_id="", GetDirList = False, copy_only = False):
    """
    Call all other routines

    Set start and end years to download
    Allow for range of stations to be processed
    """
    if restart_id != "":
        print "Restart ID: {}".format(restart_id)
    if end_id != "":
        print "End ID: {}".format(end_id)

    # Allow for check before running and overwriting files
    print "Old Data Location {}".format(OLD_ISD_DATA_LOCS)
    print "New Data Location {}".format(ISD_DATA_LOCS)

    print "Downloading or copying data from {} to {}".format(OLD_ISD_DATA_LOCS,ISD_DATA_LOCS)

    # extra diagnostic output
    IsVerbose=True  

    # run main programme
    
    # test all required files exist
    for slfile in STATIONLISTFILES:
        if not os.path.exists(os.path.join(INPUT_FILE_LOCS, slfile)):
            print "Missing file required for run {}".format(os.path.join(INPUT_FILE_LOCS, slfile))
            sys.exit("Please amend and re-run script")
        
    # data stored in years

    print datetime.datetime.now().strftime('%d-%m-%Y %H:%M')

    # make overall list - required if any of the station input files change
    CombineHadISDStations(STATIONLISTFILES,'MasterList.txt', restart_id = restart_id, end_id = end_id)
    
    # get list of folders on server (parent directory to all of the years)
    GetFileListing(HOST, ISD_LOC)
    folder_list_file = file("{}/isd_ftp_list_parent.txt".format(INPUT_FILE_LOCS), "r")
    folder_list = folder_list_file.readlines()
    folder_list = ExtractFileData(folder_list, directory = True)
                              
    # do backwards so that get most recent years (i.e. updated ones) first
    for year in range(ENDYEAR, STARTYEAR-1, -1):
        YearStartTime = dt.datetime.now()
        
        # When was last download time?
        #    test both old and new location (in case it's an updating download)
        old_lastrun = LastRunDate(os.path.join(OLD_ISD_DATA_LOCS, 'last_download_{}.txt'.format(year)))
        new_lastrun = LastRunDate(os.path.join(ISD_DATA_LOCS, 'last_download_{}.txt'.format(year)))
        
        lastrun = max(old_lastrun, new_lastrun)
    
        print "Script last successfully run on {}".format(datetime.datetime.strftime(lastrun, '%d-%b-%Y %H:%M'))

        print year
        last_folder_change = min(old_lastrun, new_lastrun)
        for folder in folder_list:
            if folder.filen == str(year):
                last_folder_change = folder.date
                break

        # have there been any updates to this folder since the last run
        if last_folder_change > lastrun or GetDirList:
            print "getting FTP list"
            # if so, get list of files on server
            GetFileListing(HOST,ISD_LOC,year = year)

        # now read in the list.
        dir_list_file = file("{}/isd_ftp_list_{}.txt".format(INPUT_FILE_LOCS, year), "r")
        dir_list = dir_list_file.readlines()

        # get list of files & properties on server using DataFile Class
        isd_files=ExtractFileData(dir_list)

        hadisd_stations=HadISDStationExtract(r'MasterList.txt')

        FileGetList=FilterStations(hadisd_stations,isd_files)

        FileDeleteList=CheckToDelete(hadisd_stations,isd_files)

        FileFilter=CheckToDownload(FileGetList, lastrun, copy_only = copy_only)
        print "\n*****Downloading*****\n"
        DownloadStations(FileGetList,FileFilter,year, verbose=IsVerbose)

        print "\n*****Deleting*****\n"
        DeleteFiles(FileDeleteList, year, verbose=IsVerbose)

        # output record of when last download occurred - using start time for this year
        #   in case of interruption or no successful completion, then this isn't overwritten
        try:
            logfile=os.path.join(ISD_DATA_LOCS, 'last_download_{}.txt'.format(year))
            with open(logfile,'w') as lf:
                lf.write("Last download commenced at %s \n" %  YearStartTime.strftime('%d-%m-%Y %H:%M'))
        except IOError:
            print "Could not create log file {}".format(os.path.join(ISD_DATA_LOCS, 'last_download_{}.txt'.format(year)))
            raise IOError

        # successful download, so remove the log
        try:
            os.remove(os.path.join(ISD_DATA_LOCS, 'failed_download_{}.txt'.format(year)))
        except OSError:
            # if doesn't exist, then just continue
            pass

        print "\nCompleted year %i\n %s" % (year, datetime.datetime.now().strftime('%d-%m-%Y %H:%M'))
        # move to next year

    print "finished"
    print datetime.datetime.now().strftime('%d-%m-%Y %H:%M')

#------------------------------------------------------------


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--start', dest='STARTYEAR', action='store', default = 1931,
                        help='Start year, default = 1931')
    parser.add_argument('--end', dest='ENDYEAR', action='store', default = datetime.datetime.now().year,
                        help='End year, default = current')
    parser.add_argument('--restart_id', dest='restart_id', action='store', default = "",
                        help='Restart ID for truncated run, default = ""')
    parser.add_argument('--end_id', dest='end_id', action='store', default = "",
                        help='End ID for truncated run, default = ""')
    parser.add_argument('--get_dir_list', dest='get_dir_list', action='store_true', default = False,
                        help='update the directory listings, default = False')
    parser.add_argument('--copy_only', dest='copy_only', action='store_true', default = False,
                        help='Only copy files from archive, no download, default = False')

    args = parser.parse_args()

    STARTYEAR=int(args.STARTYEAR)
    ENDYEAR=int(args.ENDYEAR)
    restart_id=args.restart_id
    end_id=args.end_id
    get_dir_list=args.get_dir_list
    copy_only=args.copy_only

    main(STARTYEAR, ENDYEAR, restart_id = restart_id, end_id = end_id, GetDirList = get_dir_list, copy_only = copy_only)

    print "***************************************************"
    print " NOW CHECK THAT THE FILES DOWNLOADED SUCCESSFULLY "
    print " SEE IF FILES ON FTP WERE UPDATED DURING TRANSFER "
    print "***************************************************"

#------------------------------------------------------------
# END
