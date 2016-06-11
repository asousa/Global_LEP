# Calculate global precipitation due to GLD data
# (in parallel, on nansen)

import os
import datetime
import random
import string
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



# Initialize MPI:
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
host = commands.getoutput("hostname")

in_lat_grid = np.arange(-70, 70, step=0.1)
out_lat_grid = np.arange(-90,90,step=1)
out_lon_grid = np.arange(-180,180, step=1)
window_time = 10  # seconds




task_spacing = datetime.timedelta(minutes=15)

out_dir = '/shared/users/asousa/WIPP/global_precip/outputs/15min'
db_file = '/shared/users/asousa/WIPP/global_precip/db3.pkl'
# --------------------- Prep jobs to scatter ------------------------
if rank==0:


    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    print "available nodes: ",comm.Get_size()
    print "Setting up parallel runs..."
    
    sim_start = datetime.datetime(2015,11,1,0,0,0)
    sim_stop  = datetime.datetime(2015,11,2,0,0,0)

    run_starttime = time.time()
    print "start time: %s"%run_starttime
    sss = sim_start
    tasklist = []
    while sss <= sim_stop:
        tasklist.append(sss)
        sss += task_spacing



    # tasklist = pd.date_range(start = sim_start, end = sim_stop, freq = '%ss'%task_spacing).tolist()

    print "%d total tasks"%(len(tasklist))

    # np.random.shuffle(tasklist)
    nTasks = 1.0*len(tasklist)
    nProcs = 1.0*comm.Get_size()
    nSteps = np.ceil(nTasks/nProcs).astype(int)

    chunks = partition(tasklist, nProcs)
    # chunks = [tasklist[i:i+nSteps] for i in range(0, len(tasklist), nSteps)]

    print "%d total chunks"%(len(chunks))
else:
    chunks = None

# -------------------- Prep model bits -----------------------------
if rank==0:
    print "Setting up model stuff..."
    p = precip_model(database=db_file, cumsum=True)
    p.precalculate_gridded_values(in_lat_grid, out_lat_grid, p.t)

else:
    p = None


gld = GLD_file_tools('GLD_mount', prefix='GLD')


# Broadcast from root to nodes:
p = comm.bcast(p, root=0)
chunks = comm.bcast(chunks, root=0)


comm.Barrier()

chunk = comm.scatter(chunks)

print "Process %d on host %s, doing %g jobs"%(rank, host, len(chunk))


# vv----------------- Do stuff here --------------------------------

sendbuf = []

for in_time in chunk:
    print "starting at %s"%in_time
    flux, flashes = calc_global_precip(p, gld, in_time, window_time, out_lat_grid, out_lon_grid)
    sendbuf.append([[in_time.isoformat()], flux, flashes])

    # flux = []
    # flashes = []
recvbuf = comm.gather(sendbuf, root=0)

if rank==0:
    print "recvbuf is: ", np.shape(recvbuf)

    print "Concatenating results..."
    
    results = []
    [results.extend(r) for r in recvbuf]



    print "Saving results..."
    with open(os.path.join(out_dir,'dump.pkl'),'wb') as file:
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






