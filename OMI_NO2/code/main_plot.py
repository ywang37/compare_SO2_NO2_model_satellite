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
gc_omi_file = '../data/monthly/model_satellite_2013-10-01_2013-10-31.nc'
gc_omi_varname_list = [
        'Latitude', 'Latitude_e', 'Longitude', 'Longitude_e',
        'mod_NO2Trop_tp_sat_prior', 'mod_NO2Trop_tp_sat_post',
        'mod_NO2Trop_tp_sat_joint_s200',
        'sat_ColumnAmountNO2Trop_tp_sat_prior',
        'sat_ColumnAmountNO2Trop_tp_sat_post',
        'sat_ColumnAmountNO2Trop_tp_sat_joint_s200',
        ]

# OMI GC no AK
no_AK_dir = '../data/sub_OMI_NO2_no_AK/'


region_limit = [26.0, 102.5, 52.0, 120.0]

verbose = True

figdir = '../figure/'

#
# End user parameters
#####################


# read GC, OMI data
gc_omi_vars = read_nc(gc_omi_file, gc_omi_varname_list, verbose=verbose)
# subset data
lat_e_1D = gc_omi_vars['Latitude_e'][:,0]
lon_e_1D = gc_omi_vars['Longitude_e'][0,:]
i1 = get_center_index(lat_e_1D, region_limit[0])
i2 = get_center_index(lat_e_1D, region_limit[2])
j1 = get_center_index(lon_e_1D, region_limit[1])
j2 = get_center_index(lon_e_1D, region_limit[3])
for varname in gc_omi_vars:
    if varname in ['Latitude_e', 'Longitude_e']:
        gc_omi_vars[varname] = gc_omi_vars[varname][i1:i2+2,j1:j2+2]
    else:
        gc_omi_vars[varname] = gc_omi_vars[varname][i1:i2+1,j1:j2+1]
print('Latitude_e:' )
print(gc_omi_vars['Latitude_e'][:,0])
print('Latitude:' )
print(gc_omi_vars['Latitude'][:,0])
print('Longitude_e:' )
print(gc_omi_vars['Longitude_e'][0,:])
print('Longitude:' )
print(gc_omi_vars['Longitude'][0,:])
mod_NO2Trop_tp_sat_prior = gc_omi_vars['mod_NO2Trop_tp_sat_prior']
mod_NO2Trop_tp_sat_post  = gc_omi_vars['mod_NO2Trop_tp_sat_post']
mod_NO2Trop_tp_sat_joint_s200 = gc_omi_vars['mod_NO2Trop_tp_sat_joint_s200']
sat_ColumnAmountNO2Trop_tp_sat_prior = \
        gc_omi_vars['sat_ColumnAmountNO2Trop_tp_sat_prior']
sat_ColumnAmountNO2Trop_tp_sat_post  = \
        gc_omi_vars['sat_ColumnAmountNO2Trop_tp_sat_post']
sat_ColumnAmountNO2Trop_tp_sat_joint_s200  = \
        gc_omi_vars['sat_ColumnAmountNO2Trop_tp_sat_joint_s200']

# faltten data
mod_NO2Trop_tp_sat_prior_1D = mod_NO2Trop_tp_sat_prior.flatten()
mod_NO2Trop_tp_sat_post_1D  = mod_NO2Trop_tp_sat_post.flatten()
mod_NO2Trop_tp_sat_joint_s200_1D  = mod_NO2Trop_tp_sat_joint_s200.flatten()
sat_ColumnAmountNO2Trop_tp_sat_prior_1D = \
        sat_ColumnAmountNO2Trop_tp_sat_prior.flatten()
sat_ColumnAmountNO2Trop_tp_sat_post_1D  = \
        sat_ColumnAmountNO2Trop_tp_sat_post.flatten()
sat_ColumnAmountNO2Trop_tp_sat_joint_s200_1D  = \
        sat_ColumnAmountNO2Trop_tp_sat_joint_s200.flatten()

# no AK data
no_AK_sub_OMI_NO2          = np.load(no_AK_dir + '/sub_OMI_NO2.npy'         ).flatten()
no_AK_sub_prior_GC_NO2_OMI = np.load(no_AK_dir + '/sub_prior_GC_NO2_OMI.npy').flatten()
no_AK_sub_post_GC_NO2_OMI  = np.load(no_AK_dir + '/sub_post_GC_NO2_OMI.npy' ).flatten()
no_AK_sub_joint_s200_GC_NO2_OMI  = \
        np.load(no_AK_dir + '/sub_joint_s200_GC_NO2_OMI.npy' ).flatten()
no_AK_sub_OMI_num          = np.load(no_AK_dir + '/sub_OMI_num.npy'         ).flatten()

# plot taylor diaggram
fig = plt.figure(figsize=(7,7))
plt.rcParams.update({'font.size': 12})

mtd = ModTaylorDiagram(fig=fig,max_normed_std=1.2, \
        bias_vmin=-30, bias_vmax=30, \
        std_ratios=[1], cmap=gbcwpry_map, 
        title_expected='Observation')

# point (no AK)
mtd.add_prediction(no_AK_sub_OMI_NO2, no_AK_sub_prior_GC_NO2_OMI, r'', '1', 's')
mtd.add_prediction(no_AK_sub_OMI_NO2, no_AK_sub_post_GC_NO2_OMI,  r'', '1', '^')
mtd.add_prediction(no_AK_sub_OMI_NO2, 
        no_AK_sub_joint_s200_GC_NO2_OMI,  r'', '1', 'd')


# point (considering AK)
mtd.add_prediction(sat_ColumnAmountNO2Trop_tp_sat_prior_1D, 
        mod_NO2Trop_tp_sat_prior_1D, r'', '2', 's')
mtd.add_prediction(sat_ColumnAmountNO2Trop_tp_sat_post_1D,
        mod_NO2Trop_tp_sat_post_1D,  r'', '', '^')
mtd.add_prediction(sat_ColumnAmountNO2Trop_tp_sat_joint_s200_1D,
        mod_NO2Trop_tp_sat_joint_s200_1D,  r'', '', 'd')


mtd.plot()

contours = mtd.add_contours(levels=[0.3, 0.6, 0.9, 1.2, 1.5], colors='0.5')
plt.clabel(contours, inline=1, fontsize=12, fmt='%1.1f')

marker_labels = ( ('s', 'Prior GC'), \
                  ('^', r'E-NO$_2$ GC'), \
                  ('d', r'E-joint GC') )
mtd.marker_legend( marker_labels )

stringID_labels = ( (r'$1$', r'OMI L3 NO$_2$'), \
                    (r'$2$', r'OMI L2 NO$_2$ (AK)') )
mtd.stringID_legend( stringID_labels, loc='upper left' )

cbar = mtd.bias_colorbar(orientation='horizontal', shrink=0.77, pad=0.08, extend='neither')

figname = figdir + 'taylor_NO2_AK.png'
plt.savefig(figname, format='png', dpi=300)
