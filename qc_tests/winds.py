#!/usr/local/sci/bin/python
#*****************************
#
# Winds Check (WDC)
#
#   Winds checks - cross checks on speed and direction and also wind-rose shape
#
#
#************************************************************************
#                    SVN Info
#$Rev:: 104                                           $:  Revision of last commit
#$Author:: rdunn                                      $:  Author of last commit
#$Date:: 2014-12-01 16:20:58 +0000 (Mon, 01 Dec 2014) $:  Date of last commit
#************************************************************************
import numpy as np
import scipy as sp
import datetime as dt
import os
import scipy.stats as stats
from scipy.optimize import leastsq

# RJHD routines
import qc_utils as utils

DEGREEBINS = 45
PROB_THRESHOLD = 0.01

#************************************************************************
def logical_checks(station, flag_col, logfile, plots = False, diagnostics = False):
    """
    Select occurrences of wind speed and direction which are 
    logically inconsistent with measuring practices.

    From Table 2 - DeGaetano, JOAT, 14, 308-317, 1997

    :param Station station: station object
    :param array flag_col: which columns to use in QC flag array
    :param file logfile: logfile to output to
    :param bool plots: do plots?
    :param bool diagnostics: do diagnostics?

    """

    speed = getattr(station, "windspeeds")
    direction = getattr(station, "winddirs")

    # recover direction information where the speed is Zero
    fix_zero_direction = np.ma.where(np.logical_and(speed.data == 0, direction.data.mask == True))
    direction.data[fix_zero_direction] = 0
    direction.data.mask[fix_zero_direction] = False
    station.qc_flags[fix_zero_direction, flag_col[1]] = -1 # to make a note of these

    # negative speeds
    negative_speed = np.ma.where(speed.data < 0)
    station.qc_flags[negative_speed, flag_col[0]] = 1

    # negative directions (don't try to adjust)
    negative_direction = np.ma.where(direction.data < 0)
    station.qc_flags[negative_direction, flag_col[1]] = 1

    # wrapped directions (don't try to adjust)
    wrapped_direction = np.ma.where(direction.data > 360)
    station.qc_flags[wrapped_direction, flag_col[1]] = 1
    
    # no direction possible if speed == 0
    bad_direction = np.ma.where(np.logical_and(speed.data == 0, direction.data != 0))
    station.qc_flags[bad_direction, flag_col[1]] = 1

    # northerlies given as 360, not 0 --> calm
    bad_speed = np.ma.where(np.logical_and(direction.data == 0, speed.data != 0))
    station.qc_flags[bad_speed, flag_col[0]] = 1

    # and output to file/screen
    flag_locs0, = np.where(station.qc_flags[:, flag_col[0]] > 0)    # in case of direction fixes
    flag_locs1, = np.where(station.qc_flags[:, flag_col[1]] > 0)    # in case of direction fixes

    if plots or diagnostics:
        utils.print_flagged_obs_number(logfile, "Wind Logical Checks", "windspeeds", len(flag_locs0), noWrite=True)
        utils.print_flagged_obs_number(logfile, "Wind Logical Checks", "winddirs", len(flag_locs1), noWrite=True)
    else:
        utils.print_flagged_obs_number(logfile, "Wind Logical Checks", "windspeeds", len(flag_locs0))
        utils.print_flagged_obs_number(logfile, "Wind Logical Checks", "winddirs", len(flag_locs1))

    # copy flags into attribute
    station.windspeeds.flags[flag_locs0] = 1
    station.winddirs.flags[flag_locs1] = 1
    
    return # logical_checks


#************************************************************************
def plot_wind_rose(speed, direction, title, label = ""):
    """
    Plot a wind rose for each year

    :param array speed: wind speed
    :param array direction : wind direction
    :param str title: plot title
    :param str label: optional label
    """
    import matplotlib.pyplot as plt

    binspeeds=np.array([0.,2.5,5.,7.5,10.,15.,20.,30.,100.])
    semi_separation=(2.5/180.)*np.pi

    # 30 degree bins
    
    angles=np.arange(0.0,360.+DEGREEBINS,DEGREEBINS)
    angles=[((i+(DEGREEBINS/2.))/180.)*np.pi for i in angles]

    # 2-D histogram (x, y, bins=[xbins,ybins])

    good=np.where((direction.mask == False) & (speed.mask == False))
    

    hist2d, xedge, yedge=np.histogram2d(np.deg2rad(direction[good]), speed[good], bins = [angles,binspeeds], normed = True)


    # normalisation is by bin area => (np.pi/6.)*5  angle*speed
    

    plot_angles=[i+semi_separation for i in angles[:-1]]

    LegendBars=[]
    LegendLines=[]

    # YlGnBu from ColorBrewer with 8 steps
    colours = ['#FFFFD9',\
                '#EDF8B1',\
                '#C7E9B4',\
                '#7FCDBB',\
                '#41B6C4',\
                '#1D91C0',\
                '#225EA8',\
                '#0C2C84']

    fig=plt.figure(figsize=(8,8.5))
    plt.clf()
    ax=fig.add_axes([0.1,0.15,0.8,0.8], polar=True)
    for yb,ybins in enumerate(yedge[:-1]):
        if yb==0:
            bars=ax.bar(plot_angles,hist2d[:,yb], \
                        width=(DEGREEBINS/180.)*np.pi-2.*semi_separation, \
                        bottom=0,color=colours[yb])
        else:
            # left,height, width, bottom
            bars=ax.bar(plot_angles,hist2d[:,yb], \
                        width=(DEGREEBINS/180.)*np.pi-2.*semi_separation, \
                        bottom=np.sum(hist2d[:,0:yb],axis=1),color=colours[yb])
        
        LegendBars+=[bars[0]]
        if yb != len(yedge[:-1])-1:
            LegendLines+=["%i-%i" %(binspeeds[yb],binspeeds[yb+1])]
        else:
            LegendLines+=[">%i" %(binspeeds[yb])]

    ax.set_theta_direction(-1)
    ax.set_theta_zero_location("N")
    plt.figlegend(LegendBars,LegendLines,loc=8,frameon=False,ncol=4, title="Wind Speed (m/s)")

    watermarkstring="/".join(os.getcwd().split('/')[4:])+'/'+os.path.basename( __file__ )+"   "+dt.datetime.strftime(dt.datetime.now(), "%d-%b-%Y %H:%M")
    plt.figtext(0.01,0.01,watermarkstring,size=6)
    plt.figtext(0.7,0.95,title)

    if label != "":
        plt.figtext(0.7,0.9, label)
    
    

    plt.show()

    return

#************************************************************************
def wind_create_bins(indata):
    ''' create bins and bin centres from data 
    given bin width covers entire range'''

    # set up the bins
    bmins = np.ma.min(indata)
    bmaxs = np.ma.max(indata)
    
    binwidth = (bmaxs - bmins)/20.

    # get power of ten
    decimal_places = len(str(int(1/binwidth))) 
    binwidth = np.round(binwidth, decimal_places)

    bins = np.arange(0, bmaxs + (3. * binwidth), binwidth)
    bincenters = 0.5 * (bins[1:] + bins[:-1])

    return bins, bincenters # wind_create_bins

#************************************************************************
def get_histogram_norm(indata, bins):
    '''
    Adjust the normalisation of the fitted PDF to the histogram in natural units

    :param array indata: input data
    :param array bins: bins to use for histogram
    
    :returns: norm
    '''

    natural_hist, dummy = np.histogram(indata,  bins = bins)
    normed_hist, dummy = np.histogram(indata,  bins = bins, density = True)
    
    maxloc = np.argmax(natural_hist)

    return natural_hist[maxloc]/normed_hist[maxloc] # get_histogram_norm


#************************************************************************
def wind_rose_check(station, flag_col, start, end, logfile, plots = False, diagnostics = False):
    '''
    Checks for large differences in the year-to-year wind-rose shape.  
    Uses RMSE and fits Gaussian.  Finds gap in distribution to flag beyond

    :param MetStation station: station object
    :param int flag_col: which column to store the flags in
    :param datetime start: start of data
    :param datetime end: end of data
    :param bool plots: run the plots
    :param bool diagnostics: run the diagnostics

    '''

    direction = station.winddirs.data
    speed = station.windspeeds.data
    flags = station.qc_flags[:,flag_col]

    month_ranges = utils.month_starts_in_pairs(start, end)
    month_ranges_years = month_ranges.reshape(-1,12,2)

    # histogram of wind directions ( ~ unravelled wind-rose)
    bw=20
    bins = range(0,360+bw,bw)
    full_hist, binEdges = np.histogram(direction, bins = bins, normed = True)

    # use rmse as this is known (Chi-sq remains just in case)
    rmse, chisq = -np.ma.ones([month_ranges_years.shape[0]]), -np.ma.ones([month_ranges_years.shape[0]])

    # run through each year to extract RMSE's
    for y,year in enumerate(month_ranges_years):

        if len(direction[year[0][0]:year[-1][0]].compressed()) > 0:

            hist, binEdges = np.histogram(direction[year[0][0]:year[-1][0]],  bins = bins, normed = True)

            chisq[y] = np.sum((full_hist-hist)**2/(full_hist+hist))/2.
            rmse[y] = np.sqrt(np.mean((full_hist-hist)**2))

        else:
            rmse.mask[y] = True

    # now to bin up the differences and see what the fit is.
    # need to have values spread so can bin!
    if len(np.unique(rmse.compressed())) > 1:
        binEdges, bincenters = wind_create_bins(rmse)
        hist, binEdges = np.histogram(rmse,  bins = binEdges)#, density=True)

        norm = get_histogram_norm(rmse, binEdges)

        # inputs for fit
        mu = np.mean(rmse)
        std = np.std(rmse)

        # try to get decent fit to bulk of obs.
    #    initial_values = [np.max(hist), np.mean(rmse), np.std(rmse), stats.skew(rmse), stats.kurtosis(rmse)] # norm, mean, std, sk#ew, kurtosis
    #    fit = leastsq(utils.residualsGH, initial_values, [bincenters, hist, np.ones(len(hist))])
    #    res = utils.hermite2gauss(fit[0])
    #    plot_gaussian = utils.funcGH(fit[0], bincenters)

        fit = stats.rice.fit(rmse.compressed(), loc = 0, scale = np.ma.std(rmse))
        dist_pdf = stats.rice.pdf(bincenters, fit[:-2], loc=fit[-2], scale=fit[-1]) * norm

        gaussian = utils.fit_gaussian(bincenters, hist, max(hist), mu = mu, sig = std)

        # invert Gaussian to find initial threshold, then hunt for first gap beyond
        # threshold = utils.invert_gaussian(PROB_THRESHOLD, gaussian)

        # invert Rician to find initial threshold, then hunt for first gap beyond
        if dist_pdf[-1] < PROB_THRESHOLD:
            # then curve has dropped below the threshold, so can find some updated ones.
            threshold = -np.where(dist_pdf[::-1] > PROB_THRESHOLD)[0][0]
        else:
            threshold = bincenters[-1]

        n = 0
        center = np.argmax(hist)
        gap = bincenters[-1] # nothing should be beyond this

        while True:
            if center + n + 1 == len(bincenters): 
                # gone beyond edge - nothing to flag, so just break
                break

            if bincenters[center + n] < threshold:
                n += 1
                # continue moving outwards
                continue

            if hist[center + n] == 0:
                # found one
                if center + n + 1 == len(bincenters):
                    # gone beyond edge - nothing to flag - escape
                    break
                elif hist[center + n + 1] == 0:
                    # has to be two bins wide?
                    gap = bincenters[center + n]
                    break
            n += 1

        # run through each year to extract RMSE's
        for y,year in enumerate(month_ranges_years):

                if rmse[y] > gap:

                    # only flag where there are observations
                    good, = np.where(np.logical_or(direction.mask[year[0][0]:year[-1][0]] == False, speed.mask[year[0][0]:year[-1][0]] == False))

                    flags[year[0][0]:year[-1][0]][good] = 1

                    if diagnostics or plots:
                        print "Flagging {}  RMSE {} > {}".format(y+start.year, rmse[y], gap)
                elif rmse.mask[y] == False: 
                    if diagnostics or plots:
                        print "{}".format(y+start.year)



        if plots:
            import matplotlib.pyplot as plt
            # plot underlying histogram
            plt.clf()
            plot_hist = np.array([float(x) if x != 0 else 1e-1 for x in hist])
            plt.step(binEdges[1:], plot_hist, color = 'k')

            # plot the Rician distribution on top
            plt.plot(bincenters, dist_pdf, "r-", label = "Rician") 

            # plot the gaussian on top
            plt.plot(binEdges[1:], utils.gaussian(bincenters, gaussian), color = 'b', ls = ":", label = "Gaussian")
            plt.yscale("log")
            plt.ylim([0.001, 2*max(plot_hist)])

            # plot the thresholds
            plt.axvline(threshold, color = 'g')
            plt.axvline(gap, color = 'r')

            # plot flagged values in different colour
            if len(rmse[rmse > gap]) > 0:
                plt.step(binEdges[1:][bincenters >= gap], plot_hist[bincenters >= gap], color = 'r')

            # prettify
            plt.xlabel("RMSE between complete record and each year")
            plt.ylabel("Frequency")
            plt.title(station.id + " annual wind rose differences")
            plt.xlim([0, 1.1*np.ma.max(rmse)])
            plt.legend(loc = "lower right", frameon = False)

            plt.show()


            # plot all the annual wind roses, flattened out.

            plt.clf()

            hist, binEdges = np.histogram(direction,  bins = np.arange(0.0,360.+DEGREEBINS,DEGREEBINS), normed = True) 
            bincenters = (binEdges[:-1] + binEdges[1:])/2.
            plt.plot(bincenters, hist, "k-", lw = 2)

            for y,year in enumerate(month_ranges_years):
                if len(speed[year[0][0]:year[-1][0]].compressed() > 0):
                    hist, binEdges = np.histogram(direction[year[0][0]:year[-1][0]],  bins = binEdges, normed = True) 
                    plt.plot(bincenters, hist)

            plt.xlabel("Direction (degrees)")
            plt.show()

            # plot wind roses as wind roses

            plot_wind_rose(speed, direction, "{} - {}".format(station.id, "all years"))

            for y,year in enumerate(month_ranges_years):
                if len(speed[year[0][0]:year[-1][0]].compressed() > 0):
                    plot_wind_rose(speed[year[0][0]:year[-1][0]], direction[year[0][0]:year[-1][0]], "{} - {}".format(station.id, start.year + y), label = "RMSE {:6.4f}\nThreshold {:6.4f}".format(rmse[y], gap))

    # and apply the flags and output text

    flag_locs, = np.where(flags != 0)
    if plots or diagnostics:
        utils.print_flagged_obs_number(logfile, "Wind Rose Check", "windspeeds/dirs", len(flag_locs), noWrite=True)
    else:
        utils.print_flagged_obs_number(logfile, "Wind Rose Check", "windspeeds/dirs", len(flag_locs))
    
    station.qc_flags[:,flag_col] = flags

    # and flag the variables
    station.windspeeds.flags[flag_locs] = 1
    station.winddirs.flags[flag_locs] = 1

    return # wind_rose_check


#************************************************************************
def wdc(station, flag_col, start, end, logfile, plots = False, diagnostics = False):
    """
    Specific wind speed and direction checks.
    """

    # what to do about synergistic flagging - can there be more speed obs than dir obs or vv?

    logical_checks(station, flag_col[:2], logfile, plots = plots, diagnostics = diagnostics)

    wind_rose_check(station, flag_col[2], start, end, logfile, plots = plots, diagnostics = diagnostics)

    return

#************************************************************************

if __name__ == "__main__":

    print "running wind checks"
