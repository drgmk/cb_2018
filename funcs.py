import subprocess
import numpy as np
import scipy.stats
import rebound


class CBSystem(object):
    '''Class to hold a circumbinary planet setup.'''

    def __init__(self,m1 = None, f1 = None, m2 = None, f2 = None,
                      ab = None, r1 = None, r2 = None, eb = None,
                      ib = np.pi/2., wb = None, fb = None, Wb = 0.0,
                      ld11 = 0.85, ld21 = -0.09, ld12 = 0.85, ld22 = -0.09,
                      mp = 0.0 , ap = None, rp = None, ep = 0.0,
                      ip = np.pi/2., wp = None, fp = None, Wp = 0.0, t0 = 0.0):
        '''Initialisation.
            
        Parameters
        ----------
        m1, m2, mp : float
            Mass of binary components and planet, in Solar masses.
        ab, eb, ib, fb, wb, Wb : float
            Orbital elements of binary, ab in au and radians for angles.
        ap, ep, ip, fp, wp, Wp : float
            Orbital elements of planet, ap in au and radians for angles.
        fl1, fl2 : float
            Fractional flux of binary components, planet is non-luminous.
        r1, r2, rp : float
            Radius of binary components and planet, in au.
        ld11, ld21, etc. : float
            Quadratic limb darkening coefficients, ld12 is first coeff
            for second body.
        t0 : float
            Time at which orbital elements apply, in days.
        '''
        self.m1 = m1
        self.f1 = f1
        self.m2 = m2
        self.f2 = f2
        self.ab = ab
        self.r1 = r1
        self.r2 = r2
        self.eb = eb
        self.ib = ib
        self.wb = wb
        self.fb = fb
        self.Wb = Wb
        self.ld11 = ld11
        self.ld21 = ld21
        self.ld12 = ld12
        self.ld22 = ld22
        self.mp = mp
        self.ap = ap
        self.rp = rp
        self.ep = ep
        self.ip = ip
        self.wp = wp
        self.fp = fp
        self.Wp = Wp
        self.t0 = t0


def reb_cb_dvm(cb, baseBody, transitingBody, tmin, tmax, timing_precision,
               close=0.5):
    '''Return transit times for a circumbinary system.

    This method is typically faster than functions that just integrate
    forward in time.
    
    Function returns arrays of transit times and durations, in days. 
    
    Parameters
    ----------
    cb : CBSys object
        Class holding system configuration.
    baseBody, transitingBody : int
    	Which transits to extract. Use 0, 1, 2 for primary star, secondary star, planet.
    tmin, tmax : float
        Start and end times of simulation, in days.
    timing_precision : float
    	Desired precision on transit times, in years/2pi (yes the units need updating).
    close : float
        Close encounter distance in Hill units of secondary.
    '''
    # Load some initial stuff
    transittimes = []
    transitdurations = []

    # create the rebound simulation and set the units
    sim = rebound.Simulation()

    # set close encouter distance
    r_hill = cb.ab * (1-cb.eb) * ( cb.m2 / (cb.m1+cb.m2) )**(1/3.)
    sim.exit_min_distance = close * r_hill

    # put the radii into an array, to be used later for transit calculations
    R = [cb.r1,cb.r2,cb.rp]
	
    # add the three bodies to rebound
    sim.t = 0.0
    sim.add(m = cb.m1)
    comp = rebound.Particle(simulation=sim, m=cb.m2, a=cb.ab, e=cb.eb,
                            inc=cb.ib, omega=cb.wb, f=cb.fb, Omega=cb.Wb)
    sim.add(comp)
    sim.move_to_com()
    planet = rebound.Particle(simulation=sim, m=cb.mp, a=cb.ap, e=cb.ep,
                              inc=cb.ip, omega=cb.wp, f=cb.fp, Omega=cb.Wp)
    sim.add(planet)

    # below you can set the number of active particles
    # active particle means it has mass and you therefore have to calculate its gravitational force
    # on other bodies. 3=all, 2=massless planet
    sim.N_active = 3 

    # put everything to the centre of mass	
    sim.move_to_com()

    # integrate to start time
    sim.integrate(tmax = (tmin - cb.t0)/365.25 * 2 * np.pi )
	
    # p is a reference to the positions and velocities of all three bodies
    p = sim.particles
	
    # Get a reference to the period of the transiting body, used for integration timing
    if (transitingBody == 0 or transitingBody == 1):
        P_bin = (cb.ab**3/(cb.m1+cb.m2))**(1./2.)/(2*np.pi) #in years/2pi
        P_transiter = P_bin
    elif (transitingBody == 2):
        P_p = p_p0 = (cb.ap**3/(cb.m1+cb.m2))**(1./2.)/(2*np.pi) #in years/2pi
        P_transiter = P_p
		
    # Integrate the system from time 0 to tMax
    while sim.t<(tmax - cb.t0)/365.25 * 2 * np.pi :
        # The old x position of the transiting body with respect to the base body
        x_old = p[transitingBody].x - p[baseBody].x 
        # and the corresponding time
        t_old = sim.t
        # Integrate over a quarter of the planet's period
        sim.integrate(sim.t + P_transiter/4.) 
        # Calculate a new position and time
        x_new = p[transitingBody].x - p[baseBody].x # old x position of the transiting body with respect to the base body
        t_new = sim.t
        # Check if the sign has changed on the x axis (a crossing) and the planet is in front of the star (z<0)
        # Remember that as an observer we are looking down the positive z-axis
        if x_old*x_new<0. and (p[transitingBody].z - p[baseBody].z) > 0: 
            while t_new-t_old>timing_precision:	# do a bisection to a set precision
                if x_old*x_new<0.:
                    t_new = sim.t
                else:
                    t_old = sim.t
                sim.integrate((t_new+t_old)/2.)
                x_new = p[transitingBody].x - p[baseBody].x 
            # finished the bisection because the t_old and t_new are now within the precision
            # now we can store the time
            transittimes.append(sim.t) 
            # Now want to calculate the transit duration
            # Calculate impact parameter over the star being transited
            # Calculated just by looking at the y parameter of the planet with respect to the star
            bi = np.abs(p[baseBody].y - p[transitingBody].y)/R[baseBody] 
            if bi < 1: # a transit occurred because the impact parameter is less than 1
                # the transit duration is just a simple distance/speed calculation
                transitdurations.append(2*((R[baseBody]+R[transitingBody])**2. - (bi*R[baseBody])**2.)**(1./2.)/((p[transitingBody].vx-p[baseBody].vx)**2. + (p[transitingBody].vy-p[baseBody].vy)**2.)**(1./2.))
            else: # no transit occured
                transitdurations.append(0)
            # Note that the transit time is always stored, regardless of whether or not a transit actually occurs (impact parameter < 1)
            # The transit time is purely due to a crossing of the x-axis
            # A transit duration of 0 indicates that the planet missed the star.
            sim.integrate(sim.t + P_transiter/10.) # add a 10th of the planet's orbital period just to push past the transit

    # get transit times and convert back to days
    tts = np.array(transittimes)
    tds = np.array(transitdurations)

    tts /= 2.0 * np.pi / 365.25
    tts += cb.t0
    if np.sum(tds>0)>0:
        tds /= 2.0 * np.pi / 365.25
    return tts, tds


def convfm(f_in,ecc):
    '''Convert true anomaly to mean anomaly.
    
    Parameters
    ----------
    f_in : float
        True anomaly, in radians.
    ecc : float
        Eccentricity.
    '''
    tf2 = np.tan(0.5*f_in)
    fact = np.sqrt( (1.0+ecc)/(1.0-ecc) )
    bige = 2.0*np.arctan2(tf2,fact)
    bigm = bige - ecc*np.sin(bige)
    return bigm


def pd_cb(cb, times=None, filed=None, filet=None, cleanup=True,
          run='~/code/photodynam/photodynam'):
    '''Return flux from photodynam for a circumbinary system.
    
    Parameters
    ----------
    cb : CBSys object
        Class holding system configuration.
    times : list, tuple, or array
        Times to compute flux, in days.
    filed, filet : str
        File names for running photodynam.
    cleanup : bool, optional
        Remove files after running.
    run : str
        Path to photodynam.
    '''

    mcon = (365.25/2./np.pi)**2 # mass conversion factor, time is in days

    rstr = str(np.random.randint(1e5))
    if filed is None and filet is None:
        #filed = '/Users/davidarmstrong/Software/git-repos/cb/pd_input.txt'
        #filet = '/Users/davidarmstrong/Software/git-repos/cb/pd_times.txt'
        filed = '/tmp/pd1-'+rstr+'-.txt'
        filet = '/tmp/pd2-'+rstr+'-.txt'
        fileout = '/tmp/pd3-'+rstr+'-.txt'

    # write the dynamics file
    with open(filed,'w') as f:
        f.write('3 {}\n'.format(cb.t0))
        f.write('0.01 1e-16\n')                                 
        f.write('{} {} {}\n'.format(cb.m1/mcon, cb.m2/mcon, cb.mp/mcon))
        f.write('{} {} {}\n'.format(cb.r1, cb.r2, cb.rp))
        f.write('{} {} {}\n'.format(cb.f1, cb.f2, 0.0))
        f.write('{} {} {}\n'.format(cb.ld11, cb.ld12, 1.0))
        f.write('{} {} {}\n'.format(cb.ld21, cb.ld22, 0.0))
        # orbital elements, a, e, i, w, O=0.0, M
        f.write('{} {} {} {} 0.0 {}\n'.format(cb.ab, cb.eb, cb.ib, cb.wb,
                                              np.mod(convfm( cb.fb, cb.eb ),2*np.pi)))
        f.write('{} {} {} {} {} {}\n'.format(cb.ap, cb.ep, cb.ip, cb.wp, cb.Wp,
                                              np.mod(convfm( cb.fp, cb.ep ),2*np.pi)))
                 
    # write the times file, we want just flux and RV as we know the times
    with open(filet,'w') as f:
        f.write('F v\n')
        f.write(' '.join( map(str,times) ))
    
    with open(fileout,'w') as outf:
        # run the code and clean up (1.17s per call, +0.05s on file writing)
        x = subprocess.run([run,filed,filet],stdout=outf,check=True)

    res = np.genfromtxt(fileout,comments='#')
    
    flux = res[:,0]
    vel = -res[:,3]
    
    if cleanup:
        subprocess.run(['rm',filed,filet,fileout])
    
    return flux, vel


def running_mean(x, N):
    '''Efficient running mean (needs gap free data)
    
    Parameters
    ----------
    x : ndarray
    	Data array to calculate running mean on (cannot have gaps)
    N : number of consecutive datapoints to average
    '''
    cumsum = numpy.cumsum(numpy.insert(x, 0, 0)) 
    return (cumsum[N:] - cumsum[:-N]) / float(N)


def running_mean_gaps(x, y, window, minpoints=3, blur=0):
    ''' Slow-ish implementation of a running mean - but deals with data gaps
    
    Parameters
    ----------
    x : 		ndarray
    			x-axis values (e.g. time)
    y : 		ndarray
    			y-axis values
    window : 	float
    			range in x to average over
    minpoints : int
    			minimum points to consider for an average. Fills in with 1000 if below
    blur :		float
    			window over which to blur the running mean. Creates two arrays,
    			one using maximum running mean values obver this window, the other using
    			minimum
    	
    Returns
    ----------
    stat : 				running mean
    blurredstat : 		if blur is not 0, the blurred stat using maximum values
    sourcetimes : 		if blur is not 0, the time from which the maximum came from 
    					for each point in blurredstat
    antiblurredstat: 	if blur is not 0, the blurred stat using minimum values
    '''
    stat = np.zeros(len(x))
    blurredstat = np.zeros(len(x))
    antiblurredstat = np.zeros(len(x))
    sourcetimes = np.zeros(len(x))
    scan = window / 2.
    blurscan = blur / 2.
    for point in np.arange(len(x)):
        start = np.searchsorted(x,x[point]-scan)
        end = np.searchsorted(x,x[point]+scan)
        if end-start >= minpoints:
            stat[point] = np.mean(y[start:end])
        else:
            stat[point] = 1000.
    
    if blur:
        for point in np.arange(len(x)):
            start = np.searchsorted(x,x[point]-blurscan)
            end = np.searchsorted(x,x[point]+blurscan)
            if end>start:
                segidx = np.argmin(stat[start:end])
                antisegidx = np.argmax(stat[start:end])
                blurredstat[point] = stat[start:end][segidx]
                sourcetimes[point] = x[start:end][segidx]
                antiblurredstat[point] = stat[start:end][antisegidx]
            else:
                blurredstat[point] = 1000.
                sourcetimes[point] = 1000.
                antiblurredstat[point] = 1000.
    return stat, blurredstat, sourcetimes, antiblurredstat


def extract_transit_window_sourcetimes(tt,dur,time,flux,windows,sourcetimes,statdict,ndurs=11.):
    '''
    Redoes scan over a window to find what event was selected for a transit.
    Uses source of statistic times calculated when statistic was made.
    
    Parameters
    ----------
    tt : 			float
    				transit time
    dur : 			float
    				transit duration
    time : 			ndarray
    				time array
    flux : 			ndarray
    				flux data
    windows :		list
    				list of durations used to construct statdict
    sourcetimes : 	dict
    				each key is a value of windows, containing the sourcetimes output of
    				a call to running_mean_gaps
    statdict :		dict
    				each key is a value of windows, containing the statistic being used
    				to analyse the lightcurve
    ndurs :			int
    				number of transit durations to normalise output to
    	
    Returns
    ----------
    time_window	:	segment of time array around transit
    flux_window	:	segment of flux array around transit
    timescale	:	normalied time centred on transit time, normalised by transit
    				duration * ndurs
    '''
    window = windows[np.argmin(np.abs(windows-dur))]

    timeindex = np.searchsorted(time,tt)
    sourcetime = sourcetimes[window][timeindex] 
    
    if timeindex < len(statdict[window]):    
        if time[timeindex] - tt < window:
                       
            #extract that flux window (and a bracket)
            start_w = np.searchsorted(time,sourcetime - window*(ndurs/2.))
            end_w = np.searchsorted(time,sourcetime + window*(ndurs/2.))

            time_window = time[start_w:end_w]
            flux_window = flux[start_w:end_w]
            timescale = time_window - sourcetime
            timescale /= window*(ndurs/2.)
            return time_window,flux_window,timescale
        else:
            return [],[],[]
    else:
        return [],[],[] 

                 
def normalise_stat(statdict,normstatdict,window=20):
    '''
    Normalises the statistic in statdict by the local std of another statistic (same keys)
    Goes point-by-point, so gaps not accounted for.
    
    Parameters
    ----------
    statdict		:	dict
    					each key a transit duration, containing the lightcurve
    					statistic for that duration
    normstatdict	:	dict
    					each key a transit duration, containing a lightcurve statistic
    					for that duration. Used to normalise statdict.
    window			:	int
    					number of data points to consider when taking the std of
    					normstatdict
    '''
    output = {}
    for key in statdict.keys():
        output[key] = np.zeros(len(statdict[key]))
        for point in range(len(statdict[key])):
            if point > window/2:
                start = int(point - window/2)
            else:
                start = 0
            end = int(point+window/2)
            #norm = np.std(normstatdict[key][start:end])
            dev = normstatdict[key][start:end] - np.median(normstatdict[key][start:end])
            MAD = np.median(np.abs(dev))
            output[key][point] = (statdict[key][point]) / MAD
            #if statdict[key][point] < -0.001:
            #    print(point)
            #    print(MAD)
            #    print(statdict[key][point])
            #    print((statdict[key][point]) / MAD)
    return output


def find_nearest(array,value):
    '''Find the nearest array entry(ies) to value(s).
    
    Both array and value must be numpy arrays.
    '''
    idx = np.searchsorted(array, value, side="left")
    idx[idx==len(array)] -= 1
    idx = idx - (np.abs(value - array[idx-1]) < np.abs(value - array[idx]))
    return array[idx]
    
def make_periodogram(tts_all,tds_all,time,windows,statdict):
    '''
    Turns a set of transit times and lightcurve stats into a periodogram
    This is for 2d only - just period and true anomaly
    
    Parameters
    ----------
    tts_all 	: 	dict
    				transit times from Nbody run
    tds_all 	:	dict
    				matching durations from Nbody run
    time		:	array
    				timestamps
    windows		:	list
					Lightcurve windows used
    statdict	:	dict
    				Calculated lightcurve statistics (blurred or otherwise)
    '''
    ppset = tts_all.keys()
    fpset = tts_all[list(tts_all.keys())[0]].keys()

    periodogram = np.zeros([len(ppset),len(fpset)])
    for ipp,pp in enumerate(ppset):
        for ifp,fp in enumerate(fpset):
            tts = tts_all[pp][fp]
            tds = tds_all[pp][fp]
            stat = 0
            timeindexes = np.searchsorted(time,tts) 
            #window = windows[np.argmin(np.abs(windows-tds[t]))]
            window = find_nearest(windows,tds)
            for t in range(len(tts)):
                if timeindexes[t] < len(statdict[window[t]]):
                    if time[timeindexes[t]] - tts[t] < window[t]:
                        #if np.isfinite(statdict[window[t]][timeindexes[t]]):
                        stat += statdict[window[t]][timeindexes[t]]

            periodogram[ipp,ifp] = stat
    return periodogram


def make_periodogram_pertransit(tts_all,tds_all,time,ppset,fpset,windows,statdict):
    '''
    Turns a set of transit times and lightcurve stats into a periodogram
    This is for 2d only - just period and true anomaly
    
    Parameters
    ----------
    tts_all 	: 	dict
    				transit times from Nbody run
    tds_all 	:	dict
    				matching durations from Nbody run
    time		:	array
    				timestamps
    ppset		:	list or array
    				planet periods trialled
    fpset		:	list or array
    				true anomalies trialled
    windows		:	list
					Lightcurve windows used
    statdict	:	dict
    				Calculated lightcurve statistics (blurred or otherwise)
    '''
    periodogram = np.zeros([len(ppset),len(fpset)])
    for ipp,pp in enumerate(ppset):
        for ifp,fp in enumerate(fpset):
            tts = tts_all[str(pp)[:6]][str(fp)[:6]]
            tds = tds_all[str(pp)[:6]][str(fp)[:6]]
            stat = 0
            tcount = 0
            timeindexes = np.searchsorted(time,tts) 
            #window = windows[np.argmin(np.abs(windows-tds[t]))]
            window = find_nearest(windows,tds)
            for t in range(len(tts)):
                if timeindexes[t] < len(statdict[window[t]]):
                    if time[timeindexes[t]] - tts[t] < window[t]:
                        #then we're not in a gap
                        #if np.isfinite(statdict[window[t]][timeindexes[t]]):
                        stat += statdict[window[t]][timeindexes[t]]
                    tcount += 1 #gaps count for per transit normalisation
            if tcount > 0:
                periodogram[ipp,ifp] = stat / tcount 
    return periodogram


def make_periodogram_prisec(tts_all,tts_all_2,tds_all,tds_all_2,
							time,ppset,fpset,windows,statdict):
    '''
    Turns a set of transit times and lightcurve stats into a periodogram
    This is for 2d only - just period and true anomaly. 
    Combines times of eclipses on primary and secondary.
    
    Parameters
    ----------
    tts_all 	: 	dict
    				transit times from Nbody run
    tts_all_2 	: 	dict
    				transit times from Nbody run
    tds_all 	:	dict
    				matching durations from Nbody run
    tds_all_2 	:	dict
    				matching durations from Nbody run
    time		:	array
    				timestamps
    ppset		:	list or array
    				planet periods trialled
    fpset		:	list or array
    				true anomalies trialled
    windows		:	list
					Lightcurve windows used
    statdict	:	dict
    				Calculated lightcurve statistics (blurred or otherwise)
    '''
    periodogram = np.zeros([len(ppset),len(fpset)])
    for ipp,pp in enumerate(ppset):
        for ifp,fp in enumerate(fpset):
            tts = tts_all[str(pp)[:6]][str(fp)[:6]]
            tts = np.hstack((tts,tts_all_2[str(pp)[:6]][str(fp)[:6]]))
            tds = tds_all[str(pp)[:6]][str(fp)[:6]]
            tds = np.hstack((tds,tds_all_2[str(pp)[:6]][str(fp)[:6]]))
            stat = 0
            timeindexes = np.searchsorted(time,tts) 
            #window = windows[np.argmin(np.abs(windows-tds[t]))]
            window = find_nearest(windows,tds)
            for t in range(len(tts)):
                if timeindexes[t] < len(statdict[window[t]]):
                    if time[timeindexes[t]] - tts[t] < window[t]:
                        #if np.isfinite(statdict[window[t]][timeindexes[t]]):
                        stat += statdict[window[t]][timeindexes[t]]
            periodogram[ipp,ifp] = stat
    return periodogram


def make_periodogram_pertransit_prisec(tts_all,tts_all_2,tds_all,tds_all_2,
										time,ppset,fpset,windows,statdict):
    '''
    Turns a set of transit times and lightcurve stats into a periodogram, per transit.
    This is for 2d only - just period and true anomaly. 
    Combines times of eclipses on primary and secondary.
    
    Parameters
    ----------
    tts_all 	: 	dict
    				transit times from Nbody run
    tts_all_2 	: 	dict
    				transit times from Nbody run
    tds_all 	:	dict
    				matching durations from Nbody run
    tds_all_2 	:	dict
    				matching durations from Nbody run
    time		:	array
    				timestamps
    ppset		:	list or array
    				planet periods trialled
    fpset		:	list or array
    				true anomalies trialled
    windows		:	list
					Lightcurve windows used
    statdict	:	dict
    				Calculated lightcurve statistics (blurred or otherwise)
    '''
    periodogram = np.zeros([len(ppset),len(fpset)])
    for ipp,pp in enumerate(ppset):
        for ifp,fp in enumerate(fpset):
            tts = tts_all[str(pp)[:6]][str(fp)[:6]]
            tts = np.hstack((tts,tts_all_2[str(pp)[:6]][str(fp)[:6]]))
            tds = tds_all[str(pp)[:6]][str(fp)[:6]]
            tds = np.hstack((tds,tds_all_2[str(pp)[:6]][str(fp)[:6]]))
            stat = 0
            tcount = 0
            timeindexes = np.searchsorted(time,tts) 
            #window = windows[np.argmin(np.abs(windows-tds[t]))]
            window = find_nearest(windows,tds)
            for t in range(len(tts)):
                if timeindexes[t] < len(statdict[window[t]]):
                    if time[timeindexes[t]] - tts[t] < window[t]:
                        #if np.isfinite(statdict[window[t]][timeindexes[t]]):
                        stat += statdict[window[t]][timeindexes[t]]
                    tcount += 1 #gaps count for per transit normalisation
            if tcount > 0:
                periodogram[ipp,ifp] = stat / tcount 
    return periodogram
