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
    cbar = m.colorbar(CS1, size="2%", pad=0.5)

    if logscale:
        ctix = np.arange(clims[0], clims[1]+1)
#         print ctix
        cbar.set_ticks(ctix)
        cbar.set_ticklabels(['$10^{%d}$'%f for f in ctix])


    # if logscale:
    #     logstr = 'log$_{10}$'
    # else:
    #     logstr=''

    if mode == 'energy':
        cbar.set_label('Energy flux [mErg/cm$^2$ sec]')
    elif mode == 'counts':
        cbar.set_label('Particle flux [el/cm$^2$ sec]')


    if flashes is not None:
        new_flash_coords= transform_coords(flashes[:,0], flashes[:,1], 
                           100*np.ones_like(flashes[:,0]), 'geomagnetic','geographic')



        xf, yf = m(new_flash_coords[:,1], new_flash_coords[:,0])
        msize = abs(flashes[:,2])*2
        mcolor= flashes[:,3]
        p2  = m.scatter(xf, yf, marker='o', c=mcolor,s=msize,alpha=0.8,edgecolor='None', cmap=plt.cm.Reds_r)
        # divider = make_axes_locatable(ax1)
        # cax2 = divider.append_axes("bottom",size="2%",pad=0.5)
        # plt.colorbar(p2, cax=cax2)

    ax1.set_title(plottime)
    plt.tight_layout()
    #plt.subplots_adjust(0,0,1,1,0,0)
    return fig




def plot_flux_polar(flux, in_lats, in_lons, logscale=True, clims=[-7,-3], num_contours=20):
    ''' two-up polar plots, Northern / Southern hemispehre. No terminator.
    '''
    lons, lats = np.meshgrid(in_lons, in_lats)

    # print np.shape(lons)
    # print np.shape(lats)

    new_coords = transform_coords(lats.ravel(), lons.ravel(),
                 100*np.ones_like(lons.ravel()),'geomagnetic','geographic')

    mag_lons = new_coords[:,1].T.reshape(np.shape(lons))
    mag_lats = new_coords[:,0].T.reshape(np.shape(lats))

    fig = plt.figure()
    m_list = []
    ax = []
    
    for plt_ind, plt_hem in enumerate(['N','S']):
#         print plt_hem
        ax.append(plt.subplot(1,2,plt_ind+1))
        if plt_hem=='N':
            m = Basemap(projection='npstere',boundinglat=40, lon_0=270,resolution='l')
        else:
            m = Basemap(projection='spstere',boundinglat=-40,lon_0=270,resolution='l')
        m_list.append(m)
        x, y = m(mag_lons, mag_lats)

        # CS1 = m.contourf(x,y,flux,cmap=plt.cm.jet,extend='both')

        m.fillcontinents(color='white',lake_color='aqua',alpha=0.5)
#         m.drawcoastlines(color='white',alpha=0.5)


        # draw parallels and meridians.
        m.drawparallels(np.arange(-80.,81.,20.),color='white',labels=[0,0,0,0])
        m.drawmeridians(np.arange(-180.,181.,20.),color='white', labels=[0,0,0,0])
    #     m1.drawcoastlines(color='white')
    #     m1.drawmapboundary()
    #     m1.fillcontinents(color='grey',alpha=0.3)

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
#         cbar = m.colorbar(CS1)
#         if plt_ind == 0:
#             cbar.set_visible(False)


#     divider = make_axes_locatable(ax[1])
#     cax = divider.append_axes("right", size='5%',pad=0.05)
    plt.subplots_adjust(wspace=0.05, hspace=0)
    cb = fig.colorbar(CS1, ax=ax, shrink=0.335, pad=0.02)

    if logscale:
        ctix = np.arange(clims[0], clims[1]+1)
#         print ctix
        cb.set_ticks(ctix)
        cb.set_ticklabels(['$10^{%d}$'%f for f in ctix])
    
    cb.set_label('Energy flux[mErg/cm$^2$ sec]',fontsize=11)
    
#     print dir(cb.config_axis)
    return fig





def plot_flux_polar_4up(flux1, flux2, in_lats, in_lons, logscale=True, clims=[-7,-3], num_contours=20):
    ''' four-up polar plots, Northern / Southern hemispehre. No terminator.
    '''
    lons, lats = np.meshgrid(in_lons, in_lats)

    # print np.shape(lons)
    # print np.shape(lats)

    new_coords = transform_coords(lats.ravel(), lons.ravel(),
                 100*np.ones_like(lons.ravel()),'geomagnetic','geographic')

    mag_lons = new_coords[:,1].T.reshape(np.shape(lons))
    mag_lats = new_coords[:,0].T.reshape(np.shape(lats))

    fig = plt.figure()
    m_list = []
    ax = []
    
    for plt_ind, plt_hem in enumerate(['N1','S1','N2','S2']):
#         print plt_hem
        ax.append(plt.subplot(1,4,plt_ind+1))

        if plt_ind in [0,1]:
            m = Basemap(projection='npstere',boundinglat=40, lon_0=270,resolution='l')
        else:
            m = Basemap(projection='spstere',boundinglat=-40,lon_0=270,resolution='l')
        m_list.append(m)
        x, y = m(mag_lons, mag_lats)

        # CS1 = m.contourf(x,y,flux,cmap=plt.cm.jet,extend='both')

        m.fillcontinents(color='white',lake_color='aqua',alpha=0.5)
#         m.drawcoastlines(color='white',alpha=0.5)


        # draw parallels and meridians.
        m.drawparallels(np.arange(-80.,81.,20.),color='white',labels=[0,0,0,0])
        m.drawmeridians(np.arange(-180.,181.,20.),color='white', labels=[0,0,0,0])
    #     m1.drawcoastlines(color='white')
    #     m1.drawmapboundary()
    #     m1.fillcontinents(color='grey',alpha=0.3)

        contours = np.linspace(clims[0],clims[1],num_contours+1)

        if plt_ind in [0,2]:
            flux = flux1
        else:
            flux = flux2

        # log scale?
        if logscale:
            pd = np.log10(flux).ravel()
        else:
            pd = flux.ravel()

        # Clip flux to min and max contours
        pd = np.clip(pd, contours[0], contours[-1])

        # Plot flux
        CS1 = plt.tricontourf(x.ravel(),y.ravel(),pd, contours,cmap=plt.cm.jet)
#         cbar = m.colorbar(CS1)
#         if plt_ind == 0:
#             cbar.set_visible(False)


#     divider = make_axes_locatable(ax[1])
#     cax = divider.append_axes("right", size='5%',pad=0.05)
    # plt.subplots_adjust(wspace=0, hspace=0)
    plt.tight_layout(pad=1, w_pad=0.5, h_pad=0.5)

    cb = fig.colorbar(CS1, ax=ax, shrink=0.72, pad=0.02)

    if logscale:
        ctix = np.arange(clims[0], clims[1]+1)
#         print ctix
        cb.set_ticks(ctix)
        cb.set_ticklabels(['$10^{%d}$'%f for f in ctix])
    
    cb.set_label('Energy flux[mErg/cm$^2$ sec]')
    
#     print dir(cb.config_axis)
    return fig





