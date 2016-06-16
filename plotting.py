from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt
import sys
from coordinate_structure import transform_coords
from mpl_toolkits.axes_grid1 import make_axes_locatable
import datetime



def plot_flux_basemap(flux, in_lats, in_lons, flashes=None, plottime=None,
                      logscale=False, clims=[0,1], num_contours=10, mode='counts'):

    lons, lats = np.meshgrid(in_lons, in_lats)

    # print np.shape(lons)
    # print np.shape(lats)

    new_coords = transform_coords(lats.ravel(), lons.ravel(),
                 100*np.ones_like(lons.ravel()),'geomagnetic','geographic')

    mag_lons = new_coords[:,1].T.reshape(np.shape(lons))
    mag_lats = new_coords[:,0].T.reshape(np.shape(lats))

    # print np.shape(mag_lons)
    # print np.shape(mag_lats)


    fig = plt.figure()
    ax1 = plt.subplot(111)
    m = Basemap(ax1,resolution='c',projection='robin', lon_0 = 0)
    x, y = m(mag_lons, mag_lats)

    # CS1 = m.contourf(x,y,flux,cmap=plt.cm.jet,extend='both')


    m.drawcoastlines(color='white')
    m.drawmapboundary()
    m.fillcontinents(color='grey',alpha=0.3)
    if plottime is not None:
        if isinstance(plottime,str):
            plottime = datetime.datetime.strptime(plottime,'%Y-%m-%dT%H:%M:%S')

        m.nightshade(plottime, alpha=0.25)


    contours = np.linspace(clims[0],clims[1],num_contours+1)

    # log scale?
    if logscale:
        pd = np.log10(flux).ravel()
    else:
        pd = flux.ravel()
    
    # Clip flux to min and max contours
    pd = np.clip(pd, contours[0], contours[-1])

    # Plot flux
    CS1 = plt.tricontourf(x.ravel(),y.ravel(),pd, contours,cmap=plt.cm.jet)
    cbar = m.colorbar(CS1)

    if logscale:
        logstr = 'log${_10}$'
    else:
        logstr=''

    if mode == 'energy':
        cbar.set_label('Energy flux %s[mErg/cm$^2$ sec]'%logstr)
    elif mode == 'counts':
        cbar.set_label('Particle flux %s[el/cm$^2$ sec]'%logstr)


    if flashes is not None:
        new_flash_coords= transform_coords(flashes[:,0], flashes[:,1], 
                           100*np.ones_like(flashes[:,0]), 'geomagnetic','geographic')



        xf, yf = m(new_flash_coords[:,1], new_flash_coords[:,0])
        msize = abs(flashes[:,2])
        mcolor= flashes[:,3]
        p2  = m.scatter(xf, yf, marker='o', c=mcolor,s=msize,alpha=0.8,edgecolor='None', cmap=plt.cm.Reds_r)
        # divider = make_axes_locatable(ax1)
        # cax2 = divider.append_axes("bottom",size="2%",pad=0.5)
        # plt.colorbar(p2, cax=cax2)

    ax1.set_title(plottime)
    plt.tight_layout()
    #plt.subplots_adjust(0,0,1,1,0,0)
    return fig


