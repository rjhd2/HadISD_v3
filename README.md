#########################################################################
# HadISD_v2

[September 2015]
Code to create HadISDv2 from start to finish.  Manual intervention still
currently required.

This is research code and not supported, however, if you have issues,
suggestions or think you've found a bug please get in touch with the
maintainers.

For details on what HadISD is please see the following references (all OA).

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

Unfortunately, as it stands, the file paths to the data, supplementary file and
image directories are hard-coded into the Python scripts.  I intend to simplify
this in the future.  These will need to be adjusted for individual runs. 
Similarly, the start and end times of the data processing.

In the main code directory two other folders are required:

-> trunk                        Contains main scripts
   |
    -> qc_tests			Contains the QC test scripts
   |
    -> input_files		Contains the input files, station lists etc
   |
    -> images			Will contain any output images

# INPUT FILES

The input files are extra text files required to run the code.  Some are
downloaded by the station_selection.py from the ISD repository for further
processing.  At a minimum, an attributes.dat file is required which will be
translated into global attributes of the final netCDF files.


#########################################################################
# RUNNING HADISDv2

1) Station Selection

Using the ISD FTP site, process the raw ISD station listings to obtain the set
of stations matching HadISD selection criteria and update the files.  

$> python2.7 station_selection.py --plots --update_files


2) Raw Data Download & NetCDF conversion

Uses the output of the station selection to download the files.  Will need to
update the location of the files in lines 47-48 of this script

$> python2.7 get_isd_files.py

Using the raw ISD files creates netCDF files and outputs in locations specified
in lines 40-42.  Last year in data can be specified with --end keyword

$> python2.7 mk_netcdf_files.py --end 2014


3) Quality Control

Run the internal checks - check file locations in lines 33-35 (37-39 in
neighbour_checks) and also the start and end times.  Individual tests can be
specified with keywords to the internal_checks, and the flags applied with the
--masking keyword for the neighbour checks.

$> python2.7 internal_checks.py --all

$> python2.7 neighbour_checks.py --masking

4) Analysis & Plotting

Print summary information, maps of the fail rates and also calculate the
humidity and heat stress indices.

$> python2.7 qc_summary.py

$> python2.7 plot_fail_rate_map.py

$> python2.7 calculate_humidity_and_indices.py


#########################################################################
# ALL FILES

calculate_humidity_and_indices.py  	Make the humidity and heat stress files
convert_to_esgf.py                 	Convert standard output to ESGF format
diurnal_example.py                  	Plot how diurnal check works
get_isd_files.py                   	Download the ISD files
heat_stress.py                     	Heat stress utilities
humidity_vars.py                   	Humidity utilities
internal_checks.py                 	Run the internal QC checks 
mk_netcdf_files.py                 	Make the raw NetCDF files
neighbour_checks.py			Run the neighbour checks
neighbour_utils.py			Utilities for the neighbour checks
netcdf_procs.py				NetCDF routines
plot_fail_rate_map.py			Plot percentage of removed data
plot_yearly_fails.py			Plot a year at at time showing flags
qc_summary.py				Create summary file of flagged obs
qc_utils.py				Utilities for QC suite
sort_canada.py				Process Canadian station list
spike_example.py			Plot how spike check works
station_selection.py			Select stations from ISD

qc_tests/
	clouds.py			Cloud Check
        distributional_gap.py		Distributional Gap Check
        diurnal_cycle.py		Diurnal Cycle Check
        streaks.py			Streaks and Repeated Value Check
        variance.py			Unusual Variance Check
        winds.py			Wind Specific Checks
        clean_up.py			Clean Up Months with high flagging rates
        humidity.py			Humidity Cross Checks
        climatological.py		Climatological Outlier Check
        duplicate_months.py		Duplicate Months Check
        frequent_values.py		Frequent Values Check
        odd_cluster.py			Odd Cluster Check
        records.py			World Record Check
        spike.py			Spike Check
        __init__.py			system file

# END
#########################################################################
