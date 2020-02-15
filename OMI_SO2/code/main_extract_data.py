"""
Created on Feburary 11, 2020

@author: Yi Wang
"""

import numpy as np
import sys

from mylib.gc_io.read_nd49 import read_bpch
from mylib.grid_utility import get_center_index

sys.path.append("/Dedicated/jwang-data/ywang/soil_NOx/shared_code/")
from sn_io import save_ave

def get_region_omi_gc_so2(filename, region_limit):
    """
    """

    var_dict = read_bpch(filename, ['IJ-AVG-$_SO2'], squeeze=True)

    latitude_e  = var_dict['latitude_e']
    longitude_e = var_dict['longitude_e']


    i1 = get_center_index(latitude_e,  region_limit[0])
    i2 = get_center_index(latitude_e,  region_limit[2])
    j1 = get_center_index(longitude_e, region_limit[1])
    j2 = get_center_index(longitude_e, region_limit[3])

    out_dict = {}

    latitude_e  = latitude_e[i1:i2+2]
    longitude_e = longitude_e[j1:j2+2]
    longitude_e, latitude_e = np.meshgrid(longitude_e, latitude_e)
    out_dict['Latitude_e']  = latitude_e
    out_dict['Longitude_e'] = longitude_e

    latitude  = var_dict['latitude'][i1:i2+1]
    longitude = var_dict['longitude'][j1:j2+1]
    longitude, latitude = np.meshgrid(longitude, latitude)
    out_dict['Latitude']  = latitude
    out_dict['Longitude'] = longitude

    var_dict['IJ-AVG-$_SO2'] = var_dict['IJ-AVG-$_SO2'] / 1e9

    gc_v_ind = 0 # GEOS-Chem VCD
    gc_v_so2 = var_dict['IJ-AVG-$_SO2'][gc_v_ind,i1:i2+1,j1:j2+1]
    out_dict['gc_v_so2'] = gc_v_so2

    omi_v_gc_ind = 5 # OMI VCD corrected by GC profile
    omi_v_gc_so2 = var_dict['IJ-AVG-$_SO2'][omi_v_gc_ind,i1:i2+1,j1:j2+1]
    out_dict['omi_v_gc_so2'] = omi_v_gc_so2

    return out_dict

#######################
# Start user parameters
#

prior_file = '/Dedicated/jwang-data/ywang/compare_SO2_NO2_model_satellite/\
OMI_SO2/GC_adj/gcadj_std/runs/v8-02-01/geos5_201310_pressure_prior/\
diagadj/gctm.omi.so2.ave.20131001.01'

post_file = '/Dedicated/jwang-data/ywang/compare_SO2_NO2_model_satellite/\
OMI_SO2/GC_adj/gcadj_std/runs/v8-02-01/geos5_201310_pressure_OMPS_SO2_iter05/\
diagadj/gctm.omi.so2.ave.20131001.01'

joint_s200_file = '/Dedicated/jwang-data/ywang/\
compare_SO2_NO2_model_satellite/OMI_SO2/GC_adj/gcadj_std/runs/v8-02-01/\
geos5_201310_pressure_joint_S200_iter07/\
diagadj/gctm.omi.so2.ave.20131001.01'

out_file = '../data/omi_gc_so2_201310.nc'

varname_list = ['IJ-AVG-$_SO2']

region_limit = [26.0, 102.5, 52.0, 120.0]

#
# End user parameters
#####################

out_dict = {}

# prior data
prior_dict = get_region_omi_gc_so2(prior_file, region_limit)
out_dict['Latitude']           = prior_dict['Latitude']
out_dict['Longitude']          = prior_dict['Longitude']
out_dict['Latitude_e']         = prior_dict['Latitude_e']
out_dict['Longitude_e']        = prior_dict['Longitude_e']
out_dict['prior_gc_v_so2']     = prior_dict['gc_v_so2']
out_dict['prior_omi_v_gc_so2'] = prior_dict['omi_v_gc_so2']

# posterior data
post_dict = get_region_omi_gc_so2(post_file, region_limit)
out_dict['post_gc_v_so2']     = post_dict['gc_v_so2']
out_dict['post_omi_v_gc_so2'] = post_dict['omi_v_gc_so2']

# joint_s200 data
joint_s200_dict = get_region_omi_gc_so2(joint_s200_file, region_limit)
out_dict['joint_s200_gc_v_so2']     = joint_s200_dict['gc_v_so2']
out_dict['joint_s200_omi_v_gc_so2'] = joint_s200_dict['omi_v_gc_so2']


# output data
save_ave(out_file, out_dict)
