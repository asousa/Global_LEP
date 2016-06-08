
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
    def __init__(self,database="database.pkl", multiple_bands=False, cumsum = False):

        self.R_earth = geopy.distance.EARTH_RADIUS
        self.D2R = np.pi/180.0
        self.path_atten = -12 # db per 1000 km attenuation (approximation of decay for earth-ionsphere waveguide)

        with open(database,'r') as file:
            self.db = pickle.load(file)

        # in_lats = sorted(self.db.keys())
        # self.multiple_bands = multiple_bands

        self.in_lats = self.db['in_lats']
        self.out_lats= self.db['out_lats']
        self.t = self.db['t']
        self.sc = self.db['consts']

        # print dir(self.sc)

        if cumsum:
            print "Modeling cumulative sums"
            N_el = np.cumsum(self.db['N_el'],axis=2)*self.sc.T_STEP
            S_el = np.cumsum(self.db['S_el'],axis=2)*self.sc.T_STEP
        else:
            N_el = self.db['N_el']
            S_el = self.db['S_el']


        self.N_interp = interpolate.RegularGridInterpolator((self.in_lats, self.out_lats, self.t), N_el, fill_value=0, bounds_error=False)
        self.S_interp = interpolate.RegularGridInterpolator((self.in_lats, self.out_lats, self.t), S_el, fill_value=0, bounds_error=False)

        self.precalculated = None
        self.pc_in_lats = None
        self.pc_out_lats = None
        self.pc_t = None

    def precalculate_gridded_values(self, in_lats, out_lats, t):
        print "Precalculating..."
        
        self.precalculated = self.get_multiple_precip_at(in_lats, out_lats, t)
        self.pc_in_lats = in_lats
        self.pc_out_lats = out_lats
        self.pc_t = t



    def get_multiple_precip_at(self, in_lats, out_lats, t):
        '''A vectorized version of get_precip_at(). This should be faster!
            in_lats, out_lats, t are all vectors.
            Returns: A numpy array of dimension [in_lats x out_lats x t]
        '''

        tx, ty, tz = np.meshgrid(in_lats, out_lats, t)

        keys =  np.array([np.abs(tx.ravel()),np.abs(ty.ravel()), tz.ravel()]).T
        # print np.shape(tx), np.shape(ty), np.shape(tz)

        # keys = cartesian([in_lats, out_lats, t])

        
        # Model is symmetric around northern / southern hemispheres (mag. dipole coordinates):
        # If in = N, out = N  --> Northern hemisphere
        #    in = N, out = S  --> Southern hemisphere
        #    in = S, out = N  --> Southern hemisphere
        #    in = S, out = S  --> Northern hemisphere
        use_southern_hemi = np.array(((tx > 0) ^ (ty > 0)).ravel())
        keys = np.abs(keys)
        # use_southern_hemi = np.array(((keys[:,0] > 0) ^ (keys[:,1] > 0)).ravel())
        # keys = np.abs(keys)

        keys_N = keys[~use_southern_hemi,:]
        keys_S = keys[ use_southern_hemi,:]

        out_data = np.zeros([len(in_lats)*len(out_lats)*len(t)])

        out_data[ use_southern_hemi] = self.S_interp(keys_S)
        out_data[~use_southern_hemi] = self.N_interp(keys_N)

        # out_data = self.get_precip_at(tx, ty, tz)


        # Indexes need to be swapped because of how meshgrid arranges things. Works tho!
        # return out_data.reshape(len(in_lats),len(out_lats),len(t))
        return out_data.reshape(len(out_lats),len(in_lats),len(t)).swapaxes(0,1)





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

        assert len(in_lat) == len(out_lat) == len(t), "Length mismatch!"

        use_southern_hemi = np.array(((in_lat > 0) ^ (out_lat > 0)).ravel())
        
        keys = np.array([np.abs(in_lat.ravel()),np.abs(out_lat.ravel()), t.ravel()]).T

        keys_N = keys[~use_southern_hemi,:]
        keys_S = keys[ use_southern_hemi,:]

        out_data = np.zeros([len(in_lat)*len(out_lat)*len(t)])

        out_data[ use_southern_hemi] = self.S_interp(keys_S)
        out_data[~use_southern_hemi] = self.N_interp(keys_N)

        return out_data



        # use_southern_hemi = (in_lat > 0) ^ (out_lat > 0)

        # inps = self.tile_keys(np.array([np.abs(in_lat), np.abs(out_lat)]), t)

        # if use_southern_hemi:
        #     return self.S_interp(inps)
        # else:
        #     return self.N_interp(inps)


    # # ------ Previous interpolation call -- works fine! but doesn't support multiple vector input
    # def get_precip_at(self, in_lat, out_lat, t):
    #     ''' in_lat:  Flash latitude (degrees)
    #         out_lat: Satellite latitude (degrees)
    #         t:       Time elapsed from flash (seconds)
    #         '''
    #     # Model is symmetric around northern / southern hemispheres (mag. dipole coordinates):
    #     # If in = N, out = N  --> Northern hemisphere
    #     #    in = N, out = S  --> Southern hemisphere
    #     #    in = S, out = N  --> Southern hemisphere
    #     #    in = S, out = S  --> Northern hemisphere


    #     use_southern_hemi = (in_lat > 0) ^ (out_lat > 0)

    #     inps = self.tile_keys(np.array([np.abs(in_lat), np.abs(out_lat)]), t)

    #     if use_southern_hemi:
    #         return self.S_interp(inps)
    #     else:
    #         return self.N_interp(inps)


    # def tile_keys(self, key1, key2):
    #     return np.vstack([np.outer(key1, np.ones(np.size(key2))), key2]).T

    def get_longitude_scaling(self, inp_lat, inp_lon, out_lat, out_lon, I0=None, db_scaling = False, dlong=0.7):
        ''' inp_lat: Scalar latitude
            inp_lon: Scalar longitude 
            out_lon: vector of longitudes to compute at
        '''
        dlong_sim = 0.7   # extra sidestep added to avoid null
        # ------------- More-realistic attempt, using same scaling factor as latitude.
        # Ratio of (2.1) --> equation (5.5) in Jacob's thesis
        dlat  = self.D2R*(out_lat - inp_lat)
        dlong = self.D2R*(out_lon - inp_lon)
        clat1 = np.cos(self.D2R*inp_lat)
        clat2 = np.cos(self.D2R*out_lat)
        slat1 = np.sin(self.D2R*inp_lat)
        slat2 = np.sin(self.D2R*out_lat)



        # Compute previous (latitude-dependent) weighting:
        dist_lat_0  = (self.sc.R_E + self.sc.H_IONO/2.0)*dlat  
        dist_long_0 = (self.sc.R_E + self.sc.H_IONO/2.0)*dlong_sim*self.D2R
        dist_iono_0 = np.hypot(dist_lat_0, dist_long_0)

        R_0  = np.hypot(  dist_iono_0, self.sc.H_IONO)
        xi_0 = np.arctan2(dist_iono_0, self.sc.H_IONO)

        # Compute current (latitude and longitude dependent) weighting:
        # (Use Haversine formula since we're moving further out than reasonable for Cartesian)
        # ----> Does not wrap around the north / south poles! What gives?
        a = pow(np.sin(dlat/2.0),2)
        b = (clat1*np.outer(clat2, pow(np.sin(dlong/2.0),2)))
        dist_iono_1 = 2.0*self.sc.R_E * np.arcsin(np.sqrt(a[:,np.newaxis] + b))
        R_1  = np.hypot(  dist_iono_1, self.sc.H_IONO)
        xi_1 = np.arctan2(dist_iono_1, self.sc.H_IONO)


        # # print np.shape(R_1)

        new_weight = np.sin(xi_1)/R_1
        old_weight = np.sin(xi_0)#/R_0

        ratio = new_weight/(old_weight[:,np.newaxis])

        # ratio = R_1*self.path_atten/1000.0

        # ratio = np.power(10, ratio/10.0)*abs((I0)/(self.sc.I0))
        if I0:
            ratio = ratio * np.abs(I0/self.sc.I0)

        # return ratio, new_weight, old_weight
        return ratio

        # # ----------- New version: Actually works, does wraparound properly -------
        # b = np.cos(self.D2R*inp_lat)*np.sin(self.D2R*(inp_lon - out_lon)/2.0)
        # dist_lon = self.R_earth*2*np.arcsin(np.abs(b))

        # vals = dist_lon*self.path_atten/1000.0

        # if db_scaling:
        #         return vals
        # else:
        #     if not I0:
        #         return np.power(10,vals/10.0)
        #     else:
        #         return np.power(10,vals/10.0)*(I0**2)/(self.sc.I0**2)
                

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







def cartesian(arrays, out=None):
    """
    Generate a cartesian product of input arrays.

    Parameters
    ----------
    arrays : list of array-like
        1-D arrays to form the cartesian product of.
    out : ndarray
        Array to place the cartesian product in.

    Returns
    -------
    out : ndarray
        2-D array of shape (M, len(arrays)) containing cartesian products
        formed of input arrays.

    Examples
    --------
    >>> cartesian(([1, 2, 3], [4, 5], [6, 7]))
    array([[1, 4, 6],
           [1, 4, 7],
           [1, 5, 6],
           [1, 5, 7],
           [2, 4, 6],
           [2, 4, 7],
           [2, 5, 6],
           [2, 5, 7],
           [3, 4, 6],
           [3, 4, 7],   
           [3, 5, 6],
           [3, 5, 7]])

    """

    arrays = [np.asarray(x) for x in arrays]
    dtype = arrays[0].dtype

    n = np.prod([x.size for x in arrays])
    if out is None:
        out = np.zeros([n, len(arrays)], dtype=dtype)

    m = n / arrays[0].size
    out[:,0] = np.repeat(arrays[0], m)
    if arrays[1:]:
        cartesian(arrays[1:], out=out[0:m,1:])
        for j in xrange(1, arrays[0].size):
            out[j*m:(j+1)*m,1:] = out[0:m,1:]
    return out








