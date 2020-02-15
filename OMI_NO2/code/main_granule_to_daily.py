"""
Created on January 13, 2020

@author: Yi Wang
"""

import datetime
import copy
import glob
import numpy as np
import sys

from mylib.io import read_nc

sys.path.append("/Dedicated/jwang-data/ywang/soil_NOx/shared_code/")
from sn_io import save_ave

#######################
# Start user parameters
#

scene_list = ['prior', 'post', 'joint_s200']

tp_list = ['tp_sat']

granule_dir = '../data/granule/'

daily_dir = '../data/daily/'

startDate = '2013-10-01'
endDate   = '2013-10-31'

# model species name
mod_spe_name = 'mod_NO2Trop'

# model AMF name
mod_amf_name = 'mod_AmfTrop'

# number var name
num_name = 'count'

# satellite species
sat_spe_name = 'sat_ColumnAmountNO2Trop'

# some variables
varname_list = [ \
        num_name, \
        'sat_AmfTrop', \
        sat_spe_name, \
        'sat_TropopausePressure', \
        ]

coord_name_list = [ 
        'Latitude',
        'Latitude_e',
        'Longitude',
        'Longitude_e'
        ]

count_thre = 3

#
# End user parameters
#####################

# add model and satellite variable names of difference scenarios.
for scene in scene_list:
    for tp in tp_list:
        
        tmp_name = '_' + tp + '_' + scene
        # model species
        varname_list.append( mod_spe_name + tmp_name )
        varname_list.append( mod_spe_name + '_AK' + tmp_name )
        # amf
        varname_list.append( mod_amf_name + tmp_name )
        # satellite species
        varname_list.append( sat_spe_name + tmp_name )

# Date
currDate   = startDate
currDate_D = datetime.datetime.strptime(currDate, '%Y-%m-%d')
endDate_D  = datetime.datetime.strptime(endDate,  '%Y-%m-%d')

out_dict = {}
while currDate_D <= endDate_D:

    # current date
    currDate = str(currDate_D)[0:10]
    print(''.join(np.full((79,), '-')))
    print('processing ' + currDate)

    # find all files
    wildcard = granule_dir + 'model_satellite_' + currDate[0:4] \
            + 'm' + currDate[5:7] + currDate[8:10] + 't*.nc'
    all_files = glob.glob(wildcard)
    all_files.sort()

    # read  data
    in_data_list = []
    in_count_list = []
    latlon_flag = True
    for i in range(len(all_files)):

        # read data
        in_filename = all_files[i]
        if latlon_flag:
            latlon_flag = False
            in_data = read_nc(in_filename, varname_list + coord_name_list, 
                    verbose=True)
            for varname in coord_name_list:
                out_dict[varname] = in_data.pop(varname)
        else:
            in_data = read_nc(in_filename, varname_list, verbose=True)
        in_data_list.append(in_data)
        in_count_list.append(in_data[num_name])

    # process granule data to daily data
    for varname in varname_list:

        print('  processing ' + varname)
        in_var_list = []

        if (varname != num_name):

            for i in range(len(all_files)):

                in_var = in_data_list[i][varname]
                in_count = in_count_list[i]

                # only use data when there are sufficient
                # pixels in grid
                in_var[in_count<count_thre] = np.nan
                in_var_list.append(in_var)

            # calculate daily
            in_var_ave = np.array(in_var_list)
            in_var_ave = np.nanmean(in_var_ave, axis=0)
            out_dict[varname] = in_var_ave 

        else:

            in_count_sum = np.array(in_count_list, dtype=int)
            in_count_sum = (in_count_sum >= count_thre)
            in_count_sum = np.sum(in_count_sum, axis=0)
            out_dict[varname] = in_count_sum

    # output data
    out_filename = daily_dir + 'model_satellite_' + \
            currDate + '.nc'
    save_ave(out_filename, out_dict)



    # go to next day
    currDate_D = currDate_D + datetime.timedelta(days=1)
