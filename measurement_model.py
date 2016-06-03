import numpy as np
import pickle
#from build_database import flux_obj
from GLD_file_tools import GLD_file_tools
from satellite import Satellite
import datetime
import ephem
from coordinate_structure import coordinate_structure
# from coordinate_structure import transform_coords
#from longitude_scaling import longitude_scaling
import os
from precip_model import precip_model

# def longitude_scaling(flash_coords, out_coords):
#     ''' Returns dB attenuation of a wave, based on the WIPP ionosphere input model '''

#     H_IONO = 1e5
#     # H_E = 5000.0
#     # Z0  = 377.0
#     # A = 5e3
#     # B = 1e5
#     path_atten = 12 # db per 1000 km (Completely handwaved, but I can't get the trig from the C code right)

#     R_earth = geopy.distance.EARTH_RADIUS
#     D2R = np.pi/180.0

#     #f = 4000 #Hz  
#     # Ionospheric absorption at output points
#     #iono_atten = ionoAbsorp(flash_coords.lat(),f)

#     # Separation in longitude, kilometers
#     dist_lon = (R_earth)*np.abs(flash_coords.lon() - out_coords.lon())*D2R*np.cos(D2R*flash_coords.lat())*1e-3
#     #return path_atten*dist_lon/1000.0# - iono_atten

#     return dist_lon*path_atten

#     # Approx attenuation factor (db per 1000 km):
#     # Note: Still ignoring vertical propagation losses

class measurement_model(object):
    '''Instantiate this bad boy to make precipitation measurements'''
    
    def __init__(self,
                 database='database.pkl',
                 GLD_root = 'alex/array/home/Vaisala/feed_data/GLD',
                 multiple_bands=False):
    
        self.m = precip_model(database, multiple_bands)
        
        self.RES_DT = self.m.db[self.m.db.keys()[0]]['RES_DT']
        self.RES_FINT = self.m.db[self.m.db.keys()[0]]['RES_FINT']
        
        
        self.td = datetime.timedelta(seconds = 5 + self.RES_FINT) # Maximum previous time to examine flashes in. 
                                                                  # Ideally 2*self.RES_INT to assure we don't miss anything.
        
        # Lightning database
        self.gld = GLD_file_tools(GLD_root, prefix='GLD')
        # Column indices
        self.lat_ind = 7;
        self.lon_ind = 8;
        self.mag_ind = 9;
        
    def get_measurement(self, in_time, coordinates, mode='continuous', bands=None):
        ''' Take a flux measurement at a given time and location, with a given sensor setting'''
        # Get flashes within timeframe:
        flashes, flash_times = self.gld.load_flashes(in_time, self.td)
        if flashes is None:
            print "No flashes found at ", in_time
            return 0.0

        flashes = flashes[:,(self.lat_ind, self.lon_ind, self.mag_ind, self.mag_ind)]
        flash_coords = transform_coords(flashes[:,0], flashes[:,1], np.zeros_like(flashes[:,0]), 'geographic', 'geomagnetic')
        flashes[:,:2] = flash_coords[:,:2]
        flashes[:,3] = [(in_time - s).microseconds*1e-6 + (in_time - s).seconds for s in flash_times]

        print "total flashes: ", np.shape(flashes)[0]
        # So! No we have an array of relevant flashes -- lat, lon, mag, time offset.
        # Let's model the flux at the satellite.
        flux = 0

        time_sampling_vector = np.linspace(-self.RES_FINT,0,np.round(self.RES_FINT/self.RES_DT))
        if mode=='continuous':
            for f in flashes:
                #print td.seconds - f[3]   
                flux += np.sum( self.m.get_precip_at(f[0], coordinates.lat(), time_sampling_vector + f[3]) *
                          self.m.get_longitude_scaling(f[0], f[1], coordinates.lon(), I0=f[2]) * self.RES_DT )

        if mode=='banded':
            for f in flashes:
                flux += np.sum(( np.array([self.m.get_multiband_precip_at(f[0],
                    coordinates.lat(), energy,
                    time_sampling_vector + f[3]) for energy in bands]) *
                    self.m.get_longitude_scaling(f[0], f[1], coordinates.lon(), I0=f[2]) * self.RES_DT ))
# #            
        return flux

if __name__== "__main__":
# -------------- Here's how to create a satellite and take some flux measurements: -------------
    #GLD_root  = 'alex/array/home/Vaisala/feed_data/GLD'
    #NLDN_root = 'alex/array/home/Vaisala/feed_data/NLDN'
    GLD_root = 'GLD'
    sat_TLE  = ["1 40378U 15003C   15293.75287141  .00010129  00000-0  48835-3 0  9990",
                "2 40378  99.1043 350.5299 0153633 201.4233 158.0516 15.09095095 39471"]

    # Satellite object:
    sat = Satellite(sat_TLE[0], sat_TLE[1],'Firebird 4')

    # Measurement object:
    f = measurement_model(database = "database_dicts.pkl", GLD_root = GLD_root, multiple_bands = True)

    # ---- Do The Thing:
    inTime = "2015-11-01T00:45:00"
    plottime = datetime.datetime.strptime(inTime,  "%Y-%m-%dT%H:%M:%S")

    sat.compute(plottime)
    sat.coords.transform_to('geomagnetic')

    # bands is a list of energy bands to sample at (depending on database, 1 thru 8)
    print "From banded measurement (all on):"
    print f.get_measurement(plottime, sat.coords, mode='banded',bands=f.m.E_bands)
    print "From single measurement:"
    print f.get_measurement(plottime, sat.coords, mode='continuous',bands=f.m.E_bands)


