#########################################################################
# HadISD_v3

[May 2019]
Code to create HadISDv3 from start to finish.  Manual intervention still
currently required to this code base as internal procedures are slightly
different.

This is research code and not supported, however, if you have issues,
suggestions or think you've found a bug please get in touch with the
maintainers.

For details on what HadISD is please see the following references (all OA).

HadISD version 3: Monthly Updates
RJH Dunn
Hadley Centre Technical Notes #103
https://www.metoffice.gov.uk/research/library-and-archive/publications/science/climate-science-technical-notes

Expanding HadISD: quality-controlled, sub-daily station data from 1931
RJH Dunn, KM Willett, DE Parker, L Mitchell,
Geosci. Instrum. Method. Data Syst., 5, 473-491, 2016
http://www.geosci-instrum-method-data-syst.net/5/473/2016/

HadISD: a quality controlled global synoptic report database for selected
variables at long-term stations from 1973-2010 
RJH Dunn, KM Willett, PW Thorne, EV Woolley, I Durre, A Dai, DE Parker, RS Vose
Climate of the Past Discussions 8, 1763-1833 (2012)
http://www.clim-past.net/8/1649/2012/cp-8-1649-2012.html

Pairwise homogeneity assessment of HadISD
RJH Dunn, KM Willett, CP Morice, DE Parker
Climate of the Past 10 (4), 1501-1522, (2014)
http://www.clim-past.net/10/1501/2014/cp-10-1501-2014.html

#########################################################################
# STRUCTURE

As an improvement over version 2, there are configuration and parameter files
which set most of the paths and directories.  However, these are stored in a way
that is separate from the code base for operational updates at the Met Office. 
A copy of these files is in the input_files directory, but the
set_paths_and_vars.py script will need adapting to point to these files.

In the main code directory two other folders are required:

-> trunk                        Contains main scripts
   |
    -> qc_tests			Contains the QC test scripts
   |
    -> input_files		Contains the input files, station lists etc
   |
    -> images			Could contain any output images

# INPUT FILES

The input files are extra text files required to run the code.  Some are
downloaded by the station_selection.py from the ISD repository for further
processing.  

At a minimum, an attributes.dat file is required which will be
translated into global attributes of the final netCDF files.  A
configuration.txt file sets versions, paths, settings and dates. And a
parameters.json contains the lists of which variables to process.


#########################################################################
# RUNNING HADISDv2

1) Station Selection

Using the ISD FTP site, process the raw ISD station listings to obtain the set
of stations matching HadISD selection criteria and update the files.  

$> python2.7 station_selection.py --plots --update_files


2) Raw Data Download & NetCDF conversion

Uses the output of the station selection to download the files.  

$> python2.7 get_isd_files.py

Using the raw ISD files creates netCDF files and outputs in locations specified
in configuration.txt.  Last year/month in data can be specified with --endyear and
--endmonth keywords

$> python2.7 mk_netcdf_files_pptn.py --endyear YYYY --endmonth MM


3) Quality Control

Run the internal checks.  Individual tests can be
specified with keywords to the internal_checks, and the flags applied with the
--masking keyword for the neighbour checks.  The internal_checks keyword
--doMonth runs a monthly update where distribution/population parameters are
only calculated from data up to the end of the last complete year

$> python2.7 internal_checks.py --all

$> python2.7 neighbour_checks.py --masking

4) Analysis & Plotting

Print summary information, maps of the fail rates and also calculate the
humidity and heat stress indices.

$> python2.7 qc_summary.py

$> python2.7 plot_fail_rate_map.py

$> python2.7 calculate_humidity_and_indices.py


#########################################################################
# RUNNING HADISDv2 with Rose/Cylc

Internally, the Cylc scheduling engine combined with the Rose suite management
system are used to run the updated in unattended mode.  If you are interested in
seeing how this works please get in touch with the maintainers.


# END
#########################################################################
