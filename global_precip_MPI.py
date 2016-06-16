# Calculate global precipitation due to GLD data
# (in parallel, on nansen)

import os
import datetime
import random
import string
import matplotlib
matplotlib.use("Agg")

import time
import numpy as np
import itertools
import commands
import subprocess # For calling shell commands with returned output
from partition import partition

# import pandas as pd
import pickle

from mpi4py import MPI  # parallel toolbox
from calc_global_precip import calc_global_precip
from precip_model import precip_model
from GLD_file_tools import GLD_file_tools

from plotting import plot_flux_basemap

import matplotlib.pyplot as plt


# Initialize MPI:
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
host = commands.getoutput("hostname")

in_lat_grid = np.arange(-60, 60, step=0.1)
out_lat_grid = np.arange(-90,90,step=1)
out_lon_grid = np.arange(-180,180, step=1)
window_time = 60  # seconds




task_spacing = datetime.timedelta(minutes=30) # Spacing between tasks
save_spacing = datetime.timedelta(days=2)     # Interval between file saves

out_dir = '/shared/users/asousa/WIPP/global_precip/outputs/2016'

fig_dir = os.path.join(out_dir, 'figures')

db_file = '/shared/users/asousa/WIPP/global_precip/db8_downsampled_oldrun.pkl'
# --------------------- Prep jobs to scatter ------------------------
if rank==0:


    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    if not os.path.exists(fig_dir):
        os.mkdir(fig_dir)



    print "available nodes: ",comm.Get_size()
    print "Setting up parallel runs..."
    
    sim_start = datetime.datetime(2016, 01, 1, 0, 0, 0)
    sim_stop  = datetime.datetime(2016, 06, 1, 0, 0, 0)

    run_starttime = time.time()
    print "start time: %s"%run_starttime
    sss = sim_start

    # tasklist = []

    # while sss<= sim_stop:
    #     tasklist.append(sss)
    #     sss += task_spacing

    tasklist_dict = dict()
    daylist = []

    # Days to do:
    while sss <= sim_stop:
        daylist.append(sss)
        sss += save_spacing

    print "Days:"
    print daylist

    # Tasks to do, organized into dict per day:
    for day_ind, day in enumerate(daylist):
        sss = day
        tasklist = []
        while sss < min(day + save_spacing, sim_stop):
            tasklist.append(sss)
            sss += task_spacing
        tasklist_dict[day.strftime('%Y-%m-%d_%H-%M-%S')] = tasklist


    for k in sorted(tasklist_dict.keys()):
        print k
        print len(tasklist_dict[k]), tasklist_dict[k][0], tasklist_dict[k][-1]

    # tasklist = pd.date_range(start = sim_start, end = sim_stop, freq = '%ss'%task_spacing).tolist()

    print "%d total days"%(len(tasklist_dict.keys()))
    # print "%d total tasks"%(len(tasklist))

    # # np.random.shuffle(tasklist)
    # nTasks = 1.0*len(tasklist)
    # nProcs = 1.0*comm.Get_size()
    # nSteps = np.ceil(nTasks/nProcs).astype(int)

    # chunks = partition(tasklist, nProcs)
    # # chunks = [tasklist[i:i+nSteps] for i in range(0, len(tasklist), nSteps)]

    # print "%d total chunks"%(len(chunks))
else:
    # chunks = None
    tasklist_dict = None

# -------------------- Prep model bits -----------------------------
if rank==0:
    print "Setting up model stuff..."
    p = precip_model(database=db_file, cumsum=True, mode='energy')
    p.precalculate_gridded_values(in_lat_grid, out_lat_grid, p.t)

else:
    p = None


gld = GLD_file_tools('GLD_mount', prefix='GLD')


# Broadcast from root to nodes:
p = comm.bcast(p, root=0)

# Broadcast tasklist from root to nodes:
tasklist_dict = comm.bcast(tasklist_dict, root=0)


for day in tasklist_dict.keys():
    if rank==0:
        print "Queueing jobs for %s"%day

        tasklist = tasklist_dict[day]
        nTasks = 1.0*len(tasklist)
        nProcs = 1.0*comm.Get_size()
        nSteps = np.ceil(nTasks/nProcs).astype(int)

        chunks = partition(tasklist, nProcs)
    else:
        tasklist = None
        chunks=None

    # chunks = comm.bcast(chunks, root=0)
    chunks = comm.bcast(chunks, root=0)


    comm.Barrier()

    chunk = comm.scatter(chunks)

    print "Process %d on host %s, doing %g jobs"%(rank, host, len(chunk))


    # vv----------------- Do stuff here --------------------------------

    sendbuf = []

    for in_time in chunk:
        print "starting at %s"%in_time
        try:
            flux, flashes = calc_global_precip(p, gld, in_time, window_time, out_lat_grid, out_lon_grid)
            # sendbuf.append([[in_time.isoformat()], flux, flashes])
            sendbuf.append([[in_time.isoformat()], flux])

            print "flux range:",np.min(flux),np.max(flux)
            fig = plot_flux_basemap(flux, out_lat_grid, out_lon_grid, flashes,
                                    plottime=in_time, logscale=True, clims=[-8,-2],
                                    num_contours=20, mode='energy')

            plt.savefig('%s/%s.png'%(fig_dir, in_time.strftime('%Y-%m-%d_%H-%M-%S')),bb_inches='tight')
            # sendbuf.append([in_time.isoformat()])
        except:
            print "Something went bad at %s"%in_time

    recvbuf = comm.gather(sendbuf, root=0)

    if rank==0:
        print "recvbuf is: ", np.shape(recvbuf)

        print "Concatenating results..."
        
        results = []
        [results.extend(r) for r in recvbuf]



        print "Saving results..."
        with open(os.path.join(out_dir,'%s.pkl'%day),'wb') as file:
            pickle.dump(results, file)

        run_stoptime = time.time()

        print "Total runtime: %d seconds"%(run_stoptime - run_starttime)


# ^^----------------- Do stuff here --------------------------------

# for in_time in chunk:
#     print "Process %d on host %s, doing time %s"%(rank, host, in_time)



# data = [rank, host, chunk]

# data = comm.gather(data, root=0)


# if rank==0:
#     print "received:"
#     for d in data:
#         print d
#         print " ---- "






