"""
Created on Feburary 11, 2020

@author: Yi Wang
"""

import matplotlib.pyplot as plt
import numpy as np
import sys

sys.path.append('/Dedicated/jwang-data/ywang/opt/anaconda3/lib/\
python3.7/site-packages')
from mylib.colormap.gbcwpry_map import gbcwpry_map
from mylib.grid_utility import get_center_index
from mylib.io import read_nc
from mylib.modified_taylor_diagrams import ModTaylorDiagram

#######################
# Start user parameters
#

# OMI GC AK
gc_omi_file = '../data/omi_gc_so2_201310.nc'
gc_omi_varname_list = [
        'Latitude', 'Latitude_e', 'Longitude', 'Longitude_e',
        'prior_gc_v_so2', 'prior_omi_v_gc_so2', 
        'post_gc_v_so2', 'post_omi_v_gc_so2',
        'joint_s200_gc_v_so2', 'joint_s200_omi_v_gc_so2'
        ]

# OMI GC no AK
no_AK_dir = '../data/sub_OMI_SO2_no_AK/'


verbose = True

figdir = '../figure/'

#
# End user parameters
#####################


# read GC, OMI data
gc_omi_vars = read_nc(gc_omi_file, gc_omi_varname_list, verbose=verbose)
print('Latitude_e:' )
print(gc_omi_vars['Latitude_e'][:,0])
print('Latitude:' )
print(gc_omi_vars['Latitude'][:,0])
print('Longitude_e:' )
print(gc_omi_vars['Longitude_e'][0,:])
print('Longitude:' )
print(gc_omi_vars['Longitude'][0,:])
prior_gc_v_so2       = gc_omi_vars['prior_gc_v_so2']
post_gc_v_so2        = gc_omi_vars['post_gc_v_so2']
joint_s200_gc_v_so2  = gc_omi_vars['joint_s200_gc_v_so2']
prior_omi_v_gc_so2       = gc_omi_vars['prior_omi_v_gc_so2']
post_omi_v_gc_so2        = gc_omi_vars['post_omi_v_gc_so2']
joint_s200_omi_v_gc_so2  = gc_omi_vars['joint_s200_omi_v_gc_so2']

# faltten data
prior_gc_v_so2_1D          = prior_gc_v_so2.flatten()
post_gc_v_so2_1D           = post_gc_v_so2.flatten()
joint_s200_gc_v_so2_1D     = joint_s200_gc_v_so2.flatten()
prior_omi_v_gc_so2_1D      = prior_omi_v_gc_so2.flatten()
post_omi_v_gc_so2_1D       = post_omi_v_gc_so2.flatten()
joint_s200_omi_v_gc_so2_1D = joint_s200_omi_v_gc_so2.flatten()

# no AK data
no_AK_sub_OMI_SO2          = np.load(no_AK_dir + '/sub_OMI_SO2.npy'         ).flatten()
no_AK_sub_prior_GC_SO2_OMI = np.load(no_AK_dir + '/sub_prior_GC_SO2_OMI.npy').flatten()
no_AK_sub_post_GC_SO2_OMI  = np.load(no_AK_dir + '/sub_post_GC_SO2_OMI.npy' ).flatten()
no_AK_sub_joint_s200_GC_SO2_OMI  = \
        np.load(no_AK_dir + '/sub_joint_s200_GC_SO2_OMI.npy' ).flatten()
no_AK_sub_OMI_num          = np.load(no_AK_dir + '/sub_OMI_num.npy'         ).flatten()

no_AK_sub_prior_GC_SO2_OMPS = np.load(no_AK_dir + '/sub_prior_GC_SO2_OMPS.npy').flatten()

flag = no_AK_sub_prior_GC_SO2_OMPS > 0.1

# plot taylor diaggram
fig = plt.figure(figsize=(7,7))
plt.rcParams.update({'font.size': 12})

mtd = ModTaylorDiagram(fig=fig,max_normed_std=2.5, \
        bias_vmin=-300, bias_vmax=300, \
        std_ratios=[1,2], cmap=gbcwpry_map, 
        title_expected='Observation')

# point (no AK)
mtd.add_prediction(no_AK_sub_OMI_SO2[flag], no_AK_sub_prior_GC_SO2_OMI[flag], r'', '1', 's')
mtd.add_prediction(no_AK_sub_OMI_SO2[flag], no_AK_sub_post_GC_SO2_OMI[flag],  r'', '1', '^')
mtd.add_prediction(no_AK_sub_OMI_SO2[flag], 
        no_AK_sub_joint_s200_GC_SO2_OMI[flag],  r'', '1', 'd')


# point (considering AK)
mtd.add_prediction(prior_omi_v_gc_so2_1D[flag], 
        prior_gc_v_so2_1D[flag], r'', '2', 's')
mtd.add_prediction(post_omi_v_gc_so2_1D[flag],
        post_gc_v_so2_1D[flag],  r'', '2', '^')
mtd.add_prediction(joint_s200_omi_v_gc_so2_1D[flag],
        joint_s200_gc_v_so2_1D[flag],  r'', '2', 'd')


mtd.plot()

contours = mtd.add_contours(levels=[0.5, 1.0, 1.5, 2.0, 2.5], colors='0.5')
plt.clabel(contours, inline=1, fontsize=12, fmt='%1.1f')

marker_labels = ( ('s', 'Prior GC'), \
                  ('^', r'E-SO$_2$ GC'), \
                  ('^', r'E-joint GC') )
mtd.marker_legend( marker_labels )

stringID_labels = ( (r'$1$', r'OMI L3 SO$_2$'), \
                    (r'$2$', r'OMI L3 SO$_2$ (AK)') )
mtd.stringID_legend( stringID_labels, loc='upper left' )


cbar = mtd.bias_colorbar(orientation='horizontal', shrink=0.77, pad=0.08, extend='neither')

figname = figdir + 'taylor_SO2_AK.png'
plt.savefig(figname, format='png', dpi=300)
