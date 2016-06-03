
import numpy as np
import pickle
#from build_database import flux_obj
from scipy import interpolate
# from sklearn.svm import SVR
# from sklearn.svm import NuSVR
from matplotlib import pyplot as plt
from coordinate_structure import coordinate_structure
import itertools
import geopy.distance



class precip_model(object):
    def __init__(self,database="database.pkl", multiple_bands=False):

        self.R_earth = geopy.distance.EARTH_RADIUS
        self.d2r = np.pi/180.0
        self.path_atten = -12 # db per 1000 km attenuation (approximation of decay for earth-ionsphere waveguide)

        with open(database,'rb') as file:
            self.db = pickle.load(file)

        # in_lats = sorted(self.db.keys())
        # self.multiple_bands = multiple_bands

        self.in_lats = self.db['in_lats']
        self.out_lats= self.db['out_lats']
        self.t = self.db['t']

        N_el = self.db['N_el']
        S_el = self.db['S_el']

        self.N_interp = interpolate.RegularGridInterpolator((self.in_lats, self.out_lats, self.t), N_el, fill_value=0, bounds_error=False)
        self.S_interp = interpolate.RegularGridInterpolator((self.in_lats, self.out_lats, self.t), S_el, fill_value=0, bounds_error=False)




    def get_multiple_precip_at(self, in_lats, out_lats, t):
        '''A vectorized version of get_precip_at(). This should be faster!
            in_lats, out_lats, t are all vectors.
            Returns: A numpy array of dimension [in_lats x out_lats x t]
        '''

        tx, ty, tz = np.meshgrid(in_lats, out_lats, t)

        keys =  np.array([np.abs(tx.ravel()),np.abs(ty.ravel()), tz.ravel()]).T

        # Model is symmetric around northern / southern hemispheres (mag. dipole coordinates):
        # If in = N, out = N  --> Northern hemisphere
        #    in = N, out = S  --> Southern hemisphere
        #    in = S, out = N  --> Southern hemisphere
        #    in = S, out = S  --> Northern hemisphere
        use_southern_hemi = np.array(((tx > 0) ^ (ty > 0)).ravel())

        keys_N = keys[~use_southern_hemi,:]
        keys_S = keys[ use_southern_hemi,:]

        out_data = np.zeros([len(in_lats)*len(out_lats)*len(t)])

        out_data[ use_southern_hemi] = self.S_interp(keys_S)
        out_data[~use_southern_hemi] = self.N_interp(keys_N)

        # Indexes need to be swapped because of how meshgrid arranges things. Works!
        return out_data.reshape(len(out_lats),len(in_lats),len(t)).swapaxes(0,1)




    # ------ Previous interpolation call -- works fine! but doesn't support multiple vector input
    def get_precip_at(self, in_lat, out_lat, t):
        ''' in_lat:  Flash latitude (degrees)
            out_lat: Satellite latitude (degrees)
            t:       Time elapsed from flash (seconds)
            '''
        # Model is symmetric around northern / southern hemispheres (mag. dipole coordinates):
        # If in = N, out = N  --> Northern hemisphere
        #    in = N, out = S  --> Southern hemisphere
        #    in = S, out = N  --> Southern hemisphere
        #    in = S, out = S  --> Northern hemisphere
        use_southern_hemi = (in_lat > 0) ^ (out_lat > 0)

        inps = self.tile_keys(np.array([np.abs(in_lat), np.abs(out_lat)]), t)

        if use_southern_hemi:
            return self.S_interp(inps)
        else:
            return self.N_interp(inps)


    # def tile_keys(self, key1, key2):
    #     return np.vstack([np.outer(key1, np.ones(np.size(key2))), key2]).T

    def get_longitude_scaling(self, inp_lat, inp_lon, out_lon, I0=None, db_scaling = False):
        ''' inp_lat: Scalar latitude
            inp_lon: Scalar longitude 
            out_lon: vector of longitudes to compute at
        '''

        # ----------- Old version: In all your computed measurements, rats (12.2.2015) ----
        #dist_lon = (self.R_earth)*np.abs(inp_lon - out_lon)*self.d2r*np.cos(self.d2r*inp_lat)*1e-3
        #vals = dist_lon*self.path_atten

        # ----------- New version: Actually works, does wraparound properly -------
        b = np.cos(self.d2r*inp_lat)*np.sin(self.d2r*(inp_lon - out_lon)/2.0)
        dist_lon = self.R_earth*2*np.arcsin(np.abs(b))

        vals = dist_lon*self.path_atten/1000.0

        if db_scaling:
                return vals
        else:
            if not I0:
                return np.power(10,vals/10.0)
            else:
                return np.power(10,vals/10.0)*(I0**2)/(self.I0**2)
                

if __name__ =="__main__":

    m = precip_model("database_dicts.pkl",multiple_bands=True)

    t = np.linspace(0,30,300)
    tmp = m.get_multiband_precip_at(30,45,5,t)

    plt.plot(t,tmp)
    plt.show()
    #tmp = [m.get_precip_at(30,45,x,t,"N") for x in m.E_bands]
    #print tmp

    # t = np.linspace(0,30,100)
    # out_lats = np.linspace(30,70,60)
    # in_lat = 19
    # res = []

    # for o in out_lats:
    #     points = m.tile_keys((in_lat, o), t)
    #     res.append(m.N_interp(points))

    # res = np.log10(np.array(res))

    # plt.pcolor(t, out_lats, res)
    # plt.clim([-4,4])
    # plt.show()











