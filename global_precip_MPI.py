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

# --------------------- Prep jobs to scatter ------------------------
if rank==0:
    print "available nodes: ",comm.Get_size()
    print "Setting up parallel runs..."
    sim_start = datetime.datetime(2015,11,1,0,0,0)
    sim_stop  = datetime.datetime(2015,11,2,0,0,0)

    window_time = 3600  # seconds


    tasklist = pd.date_range(start = sim_start, end = sim_stop, freq = '%ds'%window_time).tolist()

    # np.random.shuffle(tasklist)
    nTasks = 1.0*len(tasklist)
    nProcs = 1.0*comm.Get_size()
    nSteps = np.ceil(nTasks/nProcs).astype(int)

    chunks = [tasklist[i:i+nSteps] for i in range(0, len(tasklist), nSteps)]
else:
    chunks = None

# -------------------- Prep model bits -----------------------------
if rank==0:
    print "Setting up model stuff..."
    p = precip_model(database="db3.pkl", cumsum=True)
    # p.precalculate_gridded_values(in_lat_grid, out_lat_grid, p.t)

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





# ^^----------------- Do stuff here --------------------------------

# for in_time in chunk:
#     print "Process %d on host %s, doing time %s"%(rank, host, in_time)



data = [rank, host, chunk]

data = comm.gather(data, root=0)


if rank==0:
    print "received:"
    for d in data:
        print d
        print " ---- "
