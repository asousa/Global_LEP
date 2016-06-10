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

import pandas as pd
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
window_time = 60  # seconds




task_spacing = 1800 # seconds

out_dir = 'outputs/60sec'
# --------------------- Prep jobs to scatter ------------------------
if rank==0:
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    print "available nodes: ",comm.Get_size()
    print "Setting up parallel runs..."
    sim_start = datetime.datetime(2015,11,1,0,0,0)
    sim_stop  = datetime.datetime(2015,11,2,0,0,0)



    tasklist = pd.date_range(start = sim_start, end = sim_stop, freq = '%ss'%task_spacing).tolist()

    print "%d total tasks"%(len(tasklist))

    # np.random.shuffle(tasklist)
    nTasks = 1.0*len(tasklist)
    nProcs = 1.0*comm.Get_size()
    nSteps = np.ceil(nTasks/nProcs).astype(int)

    chunks = [tasklist[i:i+nSteps] for i in range(0, len(tasklist), nSteps)]

    print "%d total chunks"%(len(chunks))
else:
    chunks = None

# -------------------- Prep model bits -----------------------------
if rank==0:
    print "Setting up model stuff..."
    p = precip_model(database="db3.pkl", cumsum=True)
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

# if rank==0:
#     # Root node prepares to receive results
#     recvbuf = np.empty([nProcs, len(out_lat_grid), len(out_lon_grid)],dtype='d')
#     flux = np.zeros(n)
# else:
#     # Worker nodes get workin'
#     recvbuf = None

#     flux = calc_global_precip(p, gld, chunk[0], window_time, out_lat_grid, out_lon_grid)

# comm.Gather(flux, recvbuf, root=0)


# sendbuf = np.zeros([len(chunk), len(out_lat_grid), len(out_lon_grid)], dtype='i') + rank
# recvbuf = None
# if rank == 0:
#     recvbuf = np.empty([len(tasklist), len(out_lat_grid), len(out_lon_grid)], dtype='i')
# comm.Gather(sendbuf, recvbuf, root=0)
# if rank == 0:
#     for i in range(comm.Get_size()):
#         assert np.allclose(recvbuf[i,:], i)

# sendbuf = np.zeros([len(chunk), len(out_lat_grid), len(out_lon_grid)], dtype='i') + rank

# recvbuf = comm.gather(sendbuf, root=0)

# print "recvbuf is:", np.shape(recvbuf)



# flux, flashes = calc_global_precip(p, gld, chunk[0], window_time, out_lat_grid, out_lon_grid)

# # flux = np.ones([len(out_lat_grid), len(out_lon_grid)])

# print "flux size: ", np.shape(flux)
# sendbuf = flux
# # sendbuf = np.zeros([len(out_lat_grid), len(out_lon_grid)], dtype='d') + rank
# recvbuf = None


# if rank == 0:
#     recvbuf = np.empty([comm.Get_size(), len(out_lat_grid), len(out_lon_grid)], dtype='d')


# comm.Gather(sendbuf, recvbuf, root=0)
# if rank == 0:
#     for i in range(comm.Get_size()):
#         assert np.allclose(recvbuf[i,:], i), "node %d fucked up"%i

#     print "Receive buffer is:"
#     print recvbuf


sendbuf = []
for in_time in chunk:
    flux, flashes = calc_global_precip(p, gld, in_time.to_datetime(), window_time, out_lat_grid, out_lon_grid)

    sendbuf.append([[in_time.isoformat()], flux, flashes])
# for k in chunk:
#     sendbuf[k] = flux

recvbuf = comm.gather(sendbuf, root=0)

if rank==0:
    print "recvbuf is: ", np.shape(recvbuf)

    print "Concatenating results..."
    
    results = []
    [results.extend(r) for r in recvbuf]


    print "Saving results..."
    with open(os.path.join(out_dir,'dump.pkl'),'wb') as file:
        pickle.dump(results, file)



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
