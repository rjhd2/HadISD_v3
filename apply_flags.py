#!/usr/local/sci/bin/python
#*****************************
#
# applies the flags to the station data
#
#************************************************************************
#                    SVN Info
#$Rev:: 21                                            $:  Revision of last commit
#$Author:: rdunn                                      $:  Author of last commit
#$Date:: 2014-09-26 11:55:42 +0100 (Fri, 26 Sep 2014) $:  Date of last commit
#************************************************************************

import numpy as np
import scipy as sp
import os
import sys
import datetime as dt

# RJHD utilities
import netcdf_procs as ncdf
import qc_utils as utils

