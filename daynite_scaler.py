# from mpl_toolkits.basemap.solar import daynight_terminator
import datetime as dt
import numpy as np
from coordinate_structure import transform_coords
from coordinate_structure import coordinate_structure
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
# -------------------------------------------------------
# daynite_scaler.py
# 
#   Generates a scaling factor along longitude, based on 
#   whether the local point is day or night.
#
#   Completely copied from mpl_toolkits.basemap.solar, 
#   with the addition of geomagnetic coordinate rotation.
# 
#   6.15.2016 APS
#--------------------------------------------------------



class daynite_scaler():
    def __init__(self, in_time):
#         self.in_time = in_time
#         self.grid_lats = grid_lats
#         self.grid_lons = grid_lons
        
        self.daynite_map = self.compute_at(in_time)
  
    def JulianDayFromDate(self, date,calendar='standard'):
        """
        creates a Julian Day from a 'datetime-like' object.  Returns the fractional
        Julian Day (resolution 1 second).

        if calendar='standard' or 'gregorian' (default), Julian day follows Julian 
        Calendar on and before 1582-10-5, Gregorian calendar after 1582-10-15.

        if calendar='proleptic_gregorian', Julian Day follows gregorian calendar.

        if calendar='julian', Julian Day follows julian calendar.

        Algorithm:

        Meeus, Jean (1998) Astronomical Algorithms (2nd Edition). Willmann-Bell,
        Virginia. p. 63
        """
        # based on redate.py by David Finlayson.
        year=date.year; month=date.month; day=date.day
        hour=date.hour; minute=date.minute; second=date.second
        # Convert time to fractions of a day
        day = day + hour/24.0 + minute/1440.0 + second/86400.0
        # Start Meeus algorithm (variables are in his notation)
        if (month < 3):
            month = month + 12
            year = year - 1
        A = int(year/100)
        jd = int(365.25 * (year + 4716)) + int(30.6001 * (month + 1)) + \
             day - 1524.5
        # optionally adjust the jd for the switch from 
        # the Julian to Gregorian Calendar
        # here assumed to have occurred the day after 1582 October 4
        if calendar in ['standard','gregorian']:
            if jd >= 2299170.5:
                # 1582 October 15 (Gregorian Calendar)
                B = 2 - A + int(A/4)
            elif jd < 2299160.5:
                # 1582 October 5 (Julian Calendar)
                B = 0
            else:
                raise ValueError('impossible date (falls in gap between end of Julian calendar and beginning of Gregorian calendar')
        elif calendar == 'proleptic_gregorian':
            B = 2 - A + int(A/4)
        elif calendar == 'julian':
            B = 0
        else:
            raise ValueError('unknown calendar, must be one of julian,standard,gregorian,proleptic_gregorian, got %s' % calendar)
        # adjust for Julian calendar if necessary
        jd = jd + B
        return jd 

        
    def epem(self, date):
        """
        input: date - datetime object (assumed UTC)
        ouput: gha - Greenwich hour angle, the angle between the Greenwich
               meridian and the meridian containing the subsolar point.
               dec - solar declination.
        """
        dg2rad = np.pi/180.
        rad2dg = 1./dg2rad
        # compute julian day from UTC datetime object.
        # datetime objects use proleptic gregorian calendar.
        jday = self.JulianDayFromDate(date,calendar='proleptic_gregorian')
        jd = np.floor(jday) # truncate to integer.
        # utc hour.
        ut = date.hour + date.minute/60. + date.second/3600.
        # calculate number of centuries from J2000
        t = (jd + (ut/24.) - 2451545.0) / 36525.
        # mean longitude corrected for aberration
        l = (280.460 + 36000.770 * t) % 360
        # mean anomaly
        g = 357.528 + 35999.050 * t
        # ecliptic longitude
        lm = l + 1.915 * np.sin(g*dg2rad) + 0.020 * np.sin(2*g*dg2rad)
        # obliquity of the ecliptic
        ep = 23.4393 - 0.01300 * t
        # equation of time
        eqtime = -1.915*np.sin(g*dg2rad) - 0.020*np.sin(2*g*dg2rad) \
                + 2.466*np.sin(2*lm*dg2rad) - 0.053*np.sin(4*lm*dg2rad)
        # Greenwich hour angle
        gha = 15*ut - 180 + eqtime
        # declination of sun
        dec = np.arcsin(np.sin(ep*dg2rad) * np.sin(lm*dg2rad)) * rad2dg
        return gha, dec

    def daynight_terminator(self, date, delta, lonmin, lonmax):
        """
        date is datetime object (assumed UTC).
        nlons is # of longitudes used to compute terminator."""
        dg2rad = np.pi/180.
        lons = np.arange(lonmin,lonmax+0.5*delta,delta,dtype=np.float32)
        # compute greenwich hour angle and solar declination
        # from datetime object (assumed UTC).
        tau, dec = self.epem(date)
        # compute day/night terminator from hour angle, declination.
        longitude = lons + tau
        lats = np.arctan(-np.cos(longitude*dg2rad)/np.tan(dec*dg2rad))/dg2rad
        return lons, lats, tau, dec

    def compute_at(self, in_time):
        # Get terminator
        delta = 1
        lons, lats, tau, dec = self.daynight_terminator(in_time, delta, -180, 179)
        
        # # rotate to geomagnetic
        cs = coordinate_structure(lats, lons, np.zeros_like(lats),'geographic')
        # dnt_coords_geom = transform_coords(lats, lons, np.zeros_like(lats),'geographic','geomagnetic')
        cs.transform_to('geomagnetic')
        # print np.min(lons), np.max(lons)
        # print np.min(cs.lon()), np.max(cs.lon())

        # This is a pretty hacky thing -- interp1d won't extrapolate past the bounds of
        # cs.lon(), which is about 178 instead of 180. SO, we fill with a value that's
        # totally ridiculous, such that nearest_index will never select it. Without it,
        # we get little slivers outside the fill region. Go back and write an extrapolater.
        interpolator = interp1d(cs.lon(), cs.lat(), bounds_error=False, fill_value=-100000)
        # interpolator = splrep(cs.lon(), cs.lat(), k=1)
        lats = interpolator(lons)
        # print lats
        # #lats = cs.lat()
        #lons = cs.lon()




        # #plt.plot(cs.lon(), cs.lat())
        # plt.plot(lons, lats)
        # plt.figure()
        # plt.plot(lons)

        lats2 = np.arange(-90,90, delta,dtype=np.float32)
        self.grid_lats = lats2
        nlons = len(lons); nlats = len(lats2)
        lons2, lats2 = np.meshgrid(lons,lats2)
        lats = lats[np.newaxis,:]*np.ones((nlats,nlons),dtype=np.float32)
        daynight = np.ones(lons2.shape, np.int8)
        if dec > 0: # NH summer
            daynight = np.where(lats2>lats,0, daynight)
        else: # NH winter
            daynight = np.where(lats2<lats,0, daynight)
    
        # self.grid_lats = lats2[:,0]
        # self.grid_lons = lons2[0,:]
        self.grid_lons = lons
        self.in_time   = in_time
        self.daynight = daynight

        # plt.figure()
        # plt.pcolor(self.daynight)        # latz, lonz = np.meshgrid(self.grid_lats, self.grid_lons)

        # cs = coordinate_structure(latz.ravel(), lonz.ravel(), np.zeros_like(latz.ravel()),'geographic')
        # cs.transform_to('geomagnetic')
        # latz_geom = cs.lat().reshape(np.shape(latz))
        # lonz_geom = cs.lon().reshape(np.shape(lonz))
        # print np.shape(latz_geom)
        # # intrp = interp2d(latz_geom, lonz_geom, self.daynight_cartesian)



    def nearest_index(self, grid, values):
        # Find closest index of a value in an array (i.e., quick quantize to grid value)
        idx = np.searchsorted(grid, values, side="left")
        idx = np.clip(idx, 0, len(grid) - 1)
        idx_l = np.clip(idx - 1, 0, len(grid) - 1)

        idx[abs(values - grid[idx_l]) < abs(values - grid[idx])] -= 1
        return idx

    def nearest_index_wrap(self, grid, values):
        # Find closest index of a value in an array (i.e., quick quantize to grid value)
        idx = np.searchsorted(grid, values, side="left")
        idx[idx >= len(grid)] -= len(grid) - 1
        idx[idx < 0] += len(grid) - 1

        idx_l = idx - 1
        idx_l[idx_l==-1] = len(grid) -1
        #idx = np.clip(idx, 0, len(grid) - 1)
        #idx_l = np.clip(idx - 1, 0, len(grid) - 1)

        idx[abs(values - grid[idx_l]) < abs(values - grid[idx])] -= 1
        return idx
    
    def scaling_vector_at(self, lat):

        day_atten_factor = 0.05

        # # print lat
        # geom_coords = transform_coords([lat], self.grid_lons, [0],'geomagnetic','geographic')
        # # print np.shape(geom_coords)
        # # print geom_coords
        # # print self.grid_lats
        # # print self.grid_lons

        # lat_inds = self.nearest_index(self.grid_lats, geom_coords[:,0])
        # # lat_ind = self.nearest_index(self.grid_lats, [lat])
        # # print np.shape(lat_inds)
        # print np.shape(self.daynight)
        # lon_inds = self.nearest_index(self.grid_lons, geom_coords[:,1])
        # # print self.grid_lats[lat_inds]
        # # print self.grid_lons[lon_inds]
        # # print lat_inds
        # # print lon_inds
        # is_night = self.daynight[lat_inds,lon_inds]

        # plt.figure()
        # plt.imshow(self.daynight,origin='lower')
        # # plt.plot(self.grid_lons, is_night)
        # return is_night

        # Completely without transforming:

        #print "lat orig",lat
        # print self.grid_lons
        # cs = coordinate_structure([lat], self.grid_lons, [0], 'geomagnetic')
        # cs.transform_to('geographic')
        # cs.transform_to('ecef')
        # cs.transform_to('geomagnetic')
        # print "lat tran,",cs.lat()

        # plt.scatter(cs.lon(), cs.lat())
        lat_ind = self.nearest_index(self.grid_lats, [lat])
        # lat_inds = self.nearest_index_wrap(self.grid_lats, cs.lat())
        # lon_inds = self.nearest_index_wrap(self.grid_lons, cs.lon())

        # geom_coords = transform_coords([lat], [lon], [0],'geomagnetic','geographic')
        #lat_ind  = np.argmin(abs(lat - self.grid_lats))
        # is_night = [self.daynight[x,y] for x, y in zip(lat_inds, lon_inds)]
        is_night = self.daynight[lat_ind, :].squeeze()
        return (1 - day_atten_factor)*is_night + day_atten_factor
