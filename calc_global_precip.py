import numpy as np
from GLD_file_tools import GLD_file_tools
from precip_model import precip_model
import datetime as dt
from coordinate_structure import transform_coords



def calc_global_precip(p, gld, in_time, window_time, grid_lats, grid_lons):

    # Instantiate objects:
    # p = precip_model(database="db_test.pkl", cumsum=True)
    # gld = GLD_file_tools('GLD_mount',prefix='GLD')

    # How far back to load flashes from:
    lookback_time = dt.timedelta(seconds=p.t[-1] + window_time)

    in_lat_grid  = np.arange(-70,70, step=0.1)


    print "Precalculating..."
    p.precalculate_gridded_values(in_lat_grid, grid_lats, p.t)


    lat_ind = 7
    lon_ind = 8
    mag_ind = 9

    print "Loading flashes..."
    flashes, flash_times = gld.load_flashes(in_time, lookback_time)
    if flashes is None:
        print "No flashes found at ", in_time
    else:
        print np.shape(flashes)

    flashes = flashes[:,(lat_ind, lon_ind, mag_ind, mag_ind)]
    flash_coords = transform_coords(flashes[:,0], flashes[:,1], np.zeros_like(flashes[:,0]), 'geographic', 'geomagnetic')
    flashes[:,:2] = flash_coords[:,:2]
    flashes[:,3] = [(in_time - s).microseconds*1e-6 + (in_time - s).seconds for s in flash_times]


    # Mask out flashes outside the range of the interpolator:
    mask = (  (np.abs(flashes[:,0]) > 10) 
            & (np.abs(flashes[:,0]) < 60)) 
            # (flashes[:,2] > 50))

    print "%g flashes (post-filter)" % np.sum(mask)

    masked_flashes = flashes[mask, :]


    flux = np.zeros([len(grid_lats),len(grid_lons)])

    in_lat_inds  = nearest_index(p.pc_in_lats,  masked_flashes[:,0])
    t_end_inds   = nearest_index(p.pc_t,masked_flashes[:,3])
    t_start_inds = nearest_index(p.pc_t, masked_flashes[:,3] - 1.0)

    # (num_flashes x num_outlats)
    lv = (p.precalculated[in_lat_inds, : , t_end_inds] - p.precalculated[in_lat_inds, : , t_start_inds])

    flux = np.zeros([len(out_lat_grid), len(out_lon_grid)])

    for ind, f in enumerate(masked_flashes):
        scalefactor, _, _ = p.get_longitude_scaling(f[0], f[1], out_lat_grid, out_lon_grid, I0 = f[2])
        lv_single = lv[ind,:].squeeze()
        scalefactor*lv_single[:,np.newaxis]







    # for f in masked_flashes:
    #     t_end   = np.min([p.t[-1], f[3]])   # Min to account for the interpolator returning 0 when outside of range
    #     t_start = np.max([0, t_end - window_time]) # np.max([0,f[3] - p.t[-1]])
        
    #     lv_1 = p.get_multiple_precip_at([f[0]], grid_lats, [t_end])
    #     lv_0 = p.get_multiple_precip_at([f[0]], grid_lats, [t_start])
    #     lv = (lv_1 - lv_0).squeeze()
    #     scalefactor, _,_ = p.get_longitude_scaling(f[0], f[1], grid_lats, grid_lons, I0=f[2])

    #     flux += scalefactor*lv[:,np.newaxis]    

    return flux

def nearest_index(grid, values):
    # Find closest index of a value in an array (i.e., quick quantize to grid value)
    idx = np.searchsorted(grid, values, side="left")
    idx = np.clip(idx, 0, len(grid) - 1)
    idx_l = np.clip(idx - 1, 0, len(grid) - 1)
    idx[abs(values - grid[idx_l]) < abs(values - grid[idx])] -= 1
    return idx