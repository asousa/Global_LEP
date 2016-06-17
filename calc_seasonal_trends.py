import numpy as np
import pickle
import os
import datetime as dt



pklpath = '/shared/users/asousa/WIPP/global_precip/outputs/2015/saves/'
outpath = '/shared/users/asousa/WIPP/global_precip/outputs/2015/'

allfiles = os.listdir(pklpath)
filenames = sorted([f for f in allfiles if f.endswith('.pkl')])

print "we have %d files"%len(filenames)

R_e = 6378 # km^2
km2cm = 1e5 
millierg2joule=1e-10
D2R = np.pi/180.0
dLon = 1
dLat = 1

grid_lats = np.arange(-90,90, step=dLat)
grid_lons = np.arange(-180,180, step=dLon)



# Index where latitude grid is above zero
eq_ind = np.where(grid_lats > 0)[0][0]

# Area in km^2 of each latitude-longitude cell
A_lat = (R_e**2)*dLon*D2R*(np.sin(grid_lats*D2R) - np.sin((grid_lats - dLat)*D2R))
# plt.plot(A_lat)

print (360.0/dLon)*np.sum(A_lat)*1e-8 # total area of earth ~5e8 km^2


pwr_times     = np.zeros(len(filenames),dtype=dt.datetime)
NH_pwr_vector = np.zeros(len(filenames),dtype=np.float32)
SH_pwr_vector = np.zeros(len(filenames),dtype=np.float32)


for ind, filename in enumerate(allfiles):
    with open(os.path.join(pklpath,filename),'rb') as file:
        print "%d: %s"%(ind, filename)
        r = pickle.load(file)
        timestr = r[0][0]
#         print timestr
        data = r[1]
        NH_pwr = np.sum(np.dot(A_lat[eq_ind:], data[eq_ind:,:]))  # mErg/sec/(km^2/cm^2)
        NH_pwr_watts = NH_pwr*millierg2joule*(km2cm**2)

        SH_pwr = np.sum(np.dot(A_lat[0:eq_ind], data[0:eq_ind,:]))  # mErg/sec/(km^2/cm^2)
        SH_pwr_watts = SH_pwr*millierg2joule*(km2cm**2)

        time_obj = dt.datetime.strptime(timestr, '%Y-%m-%dT%H:%M:%S')
        
        pwr_times[ind] = time_obj
        NH_pwr_vector[ind] = NH_pwr_watts
        SH_pwr_vector[ind] = SH_pwr_watts
        
        
out_file = dict()

out_file['times'] = pwr_times
out_file['NH_pwr'] = NH_pwr_vector
out_file['SH_pwr'] = SH_pwr_vector

with open(os.path.join(outpath,'trends.pkl'),'wb') as file:
    print "Saving..."
    pickle.dump(out_file, file)


#         print NH_pwr_watts
#         print SH_pwr_watts

# print np.shape(mglats)
# # dLon = mglons - np.roll(mglons, -1, axis=1)
# dLon = 1
# lat2 = grid_lats
# lat1 = np.roll(grid_lats, -1)

# A = (R_e**2)*np.outer((dLon*D2R)*np.ones_like(grid_lons), (np.sin(lat2*D2R) - np.sin(lat1*D2R)))
# print np.max(A)
# print np.min(A)
# print np.sum(A[:,0:-1])
# plt.figure()
# # plt.imshow(A, origin='lower')
# plt.pcolormesh(mglats, mglons, A.T)
# # print '%e'%A