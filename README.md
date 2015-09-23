# HadISD_v2

[September 2015]
Code to create HadISDv2 from start to finish.  Manual intervention still
currently required.

This is research code and not supported, however, if you have issues,
suggestions or think you've found a bug please get in touch with the
maintainers.

For details on what HadISD is please see the following references (all OA).

HadISD: a quality controlled global synoptic report database for selected variables at long-term stations from 1973-2010
RJH Dunn, KM Willett, PW Thorne, EV Woolley, I Durre, A Dai, DE Parker, RS Vose
Climate of the Past Discussions 8, 1763-1833 (2012)
http://www.clim-past.net/8/1649/2012/cp-8-1649-2012.html

Pairwise homogeneity assessment of HadISD
RJH Dunn, KM Willett, CP Morice, DE Parker
Climate of the Past 10 (4), 1501-1522, (2014)
http://www.clim-past.net/10/1501/2014/cp-10-1501-2014.html


RUNNING HADISDv2

-- Station Selection

Using the ISD FTP site, process the raw ISD station listings to obtain the set
of stations matching HadISD selection criteria and update the files.  

>>> python2.7 station_selection.py --plots --update_files


-- Raw Data Download

Uses the output of the station selection to download the files.  Will need to
update the location of the files in lines 47-48 of this script

>>> python2.7 get_isd_files.py

Using the raw ISD files creates netCDF files and outputs in locations specified
in lines 40-42

>>> python2.7 mk_netcdf_files.py --end 2014


-- Quality Control

Run the internal checks - check file locations in lines 33-35 (37-39 in
neighbour_checks) and also the start
and end times

>>> python2.7 internal_checks.py --all

>>> python2.7 neighbour_checks.py --masking

-- Analysis & Plotting

Print summary information, maps of the fail rates and also calculate the
humidity and heat stress indices.

>>> python2.7 print_flags_summary.py

>>> python2.7 plot_fail_rate_map.py

>>> python2.7 calculate_humidity_and_indices.py
