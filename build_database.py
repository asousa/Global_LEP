import os
import glob
import re
import numpy as np
from coordinate_structure import coordinate_structure
from matplotlib import pyplot as plt
from load_phi_files import load_phi_files
from load_sim_constants import load_sim_constants
import subprocess
import re
import sys
import pickle

# Build a database (okay, just a pickle file) of flux simulations.
# This version added to global_precip project, starting edits ~6.1.2016

# class flux_obj(object):
#     def __init__(self):
#         self.t = None
#         self.coords = coordinate_structure()
#         self.data = None
#         self.N = None
#         self.N_E = None
#         self.S_E = None
#         self.S = None
#         self.L = None
#         self.consts_list = None
#         self.I0 = None
#         self.T_STEP = None
#         self.T_MAX = None
#         self.NUM_STEPS = None


def build_database(input_dir ='outputs', output_filename='database.pkl', t_new_step = None, num_L = 33, old_format = False):

    ev2joule = (1.60217657)*1e-19 # Joules / ev
    joule2millierg = 10*1e10 

    print ""
    # rootDir = os.getcwd() + '/' + input_dir_name + '/'
    d = os.listdir(input_dir)

    database = dict()

    # d = os.listdir(rootDir + "/" + r)
    ins = sorted([f for f in d if 'in_' in f])

    # Parse constants file
    consts_file = os.path.join(input_dir, "codesrc/consts.h")
    sc = load_sim_constants(consts_file, old_format = old_format)

    in_lats = np.array(sorted([int(i[3:]) for i in ins]))
    print "in lats:", in_lats


    # Should we downsample?
    if t_new_step is not None:
        t = np.arange(0,sc.T_MAX, step=t_new_step) # New time vector
        intervals = np.round(t*sc.NUM_STEPS/sc.T_MAX).astype(int) # Intervals to sum between
    else:
        t = np.arange(0,sc.T_MAX, step=sc.T_STEP)

    # Electron flux arrays -- [in_lats x out_lats x t]
    N_arr = np.zeros([len(in_lats), num_L, len(t)])
    S_arr = np.zeros([len(in_lats), num_L, len(t)])

    # Energy flux arrays  --  [in_lats x out_lats x t]
    N_e_arr = np.zeros([len(in_lats), num_L, len(t)])
    S_e_arr = np.zeros([len(in_lats), num_L, len(t)])

    # Energy at center of bins:
    E_centers = pow(10, sc.E_EXP_BOT + (sc.DE_EXP/2.0) + sc.DE_EXP*np.arange(0,sc.NUM_E))
    # bin widths in keV (matches np.diff(sc.E_tot_arr))
    dE = (1e-3)*(pow(10, sc.E_EXP_BOT + (sc.DE_EXP/2.0)) *
                np.exp(sc.DE_EXP*np.arange(0,sc.NUM_E)/np.log10(np.e))  *
                sc.DE_EXP / np.log10(np.e) )

    for in_lat_ind, in_lat in enumerate(in_lats):
    # Load some actual data!

        cur_dir = os.path.join(input_dir, "in_%g"%in_lat)
        N, S, L = load_phi_files(dataDir=cur_dir, sc=sc)            
        
        NUM_L = len(L)

        # New version (6.2016) returns N, S ~ [n_E x n_T x n_L]
        # electron flux totals:
        N_el_totals = np.inner(N.swapaxes(0,2), dE)
        S_el_totals = np.inner(S.swapaxes(0,2), dE)
        # N_el_totals = np.sum(N,axis=0).swapaxes(1,0)
        # S_el_totals = np.sum(S,axis=0).swapaxes(1,0)

        N_energy_totals = np.inner(N.swapaxes(0,2), dE*E_centers*ev2joule*joule2millierg)
        S_energy_totals = np.inner(S.swapaxes(0,2), dE*E_centers*ev2joule*joule2millierg)




        # Should we downsample?
        if t_new_step is not None:
            for t_ind, Tp in enumerate(zip(intervals[0:], intervals[1:])):
                
                # HEY:
                # I think the model units are el/(cm^2 kev sec)...
                # multipy by T_STEP to get total electrons per bin, then divide by new timestep
                N_arr[in_lat_ind,:,t_ind] = np.sum(N_el_totals[:,Tp[0]:Tp[1]],axis=1)*sc.T_STEP/t_new_step
                S_arr[in_lat_ind,:,t_ind] = np.sum(S_el_totals[:,Tp[0]:Tp[1]],axis=1)*sc.T_STEP/t_new_step
                
                N_e_arr[in_lat_ind,:,t_ind] = np.sum(N_energy_totals[:,Tp[0]:Tp[1]],axis=1)*sc.T_STEP/t_new_step
                S_e_arr[in_lat_ind,:,t_ind] = np.sum(S_energy_totals[:,Tp[0]:Tp[1]],axis=1)*sc.T_STEP/t_new_step

        else:
            N_arr[in_lat_ind,:,:] = N_el_totals
            S_arr[in_lat_ind,:,:] = S_el_totals

            N_e_arr[in_lat_ind,:,:] = N_energy_totals
            S_e_arr[in_lat_ind,:,:] = S_energy_totals


        # N_energy_totals = np.dot(N, E_centers)
        # S_energy_totals = np.dot(S, E_centers)

        # Load each (lat, Lk, I0, L-shell) vector separately
        coords = coordinate_structure(L,[0],[100],'L_dipole')
        coords.transform_to('geomagnetic')

    database['N_el'] = N_arr        # el/(cm^2 sec)
    database['S_el'] = S_arr
    database['N_energy'] = N_e_arr  # mErg/(cm^2 sec)
    database['S_energy'] = S_e_arr
    database['t'] = t
    database['L'] = L
    database['in_lats'] = in_lats
    database['out_lats']= coords.lat()
    database['consts']  = sc


    print "Saving database"
    with open(output_filename,'wb') as f:
        pickle.dump(database,f,pickle.HIGHEST_PROTOCOL)


if __name__ == "__main__":
    build_database(input_dir_name="outputs/probably/",output_filename="database_dicts.pkl")

