"""
Created on January 27, 2020

@author: Yi Wang
"""

import datetime
import glob
import numpy as np

from mylib.conversion import vmr_to_DU
from mylib.gc_io.read_nd49 import read_nd49_resample
from mylib.pro_omps_no2_l2.io_omps_no2_l2 import read_OMPS_NO2_L2


#######################
# Start user parameters
#

gc_root_dir = '/Dedicated/jwang-data/ywang/OMPS_NO2_inverse/GC_adj/\
gcadj_std/runs/v8-02-01/geos5_201310_pressure/'

sat_dir = '/Dedicated/jwang-data/shared_satData/OMPS/NMNO2-L2/AK201310/'

gc_trop_file = gc_root_dir + 'ctm.bpch'

gc_nd49_dir = gc_root_dir + 'ND49.1/'

startDate = '2013-10-08'
endDate   = '2013-10-08'

# species names
mod_spename = 'TIME-SER_NO2'

# all model variable names
mod_varname_list = [mod_spename, 'TIME-SER_AIRDEN', \
        'BXHGHT-$_BXHEIGHT', 'PEDGE-$_NOx']

verbose = True

#
# End user parameters
#####################

# Date
currDate   = startDate
currDate_D = datetime.datetime.strptime(currDate, '%Y-%m-%d')
endDate_D  = datetime.datetime.strptime(endDate,  '%Y-%m-%d')


while currDate_D <= endDate_D:

    # current date
    currDate = str(currDate_D)[0:10]
    print(''.join(np.full((79,), '-')))
    print('processing ' + currDate)

    # A satelite file may span two days. Thus, we need to
    # find model data from previous day and next day
    preDate_D  = currDate_D + datetime.timedelta(days=-1)
    nextDate_D = currDate_D + datetime.timedelta(days=1)
    preDate  = str(preDate_D)
    nextDate = str(nextDate_D)
    # all dates
    p_date = preDate[0:4]  + preDate[5:7]  + preDate[8:10]
    c_date = currDate[0:4] + currDate[5:7] + currDate[8:10]
    n_date = nextDate[0:4] + nextDate[5:7] + nextDate[8:10]
    # filenames for previous, current, and next day.
    pre_files  = gc_nd49_dir + 'ts' + p_date + '.bpch'
    curr_files = gc_nd49_dir + 'ts' + c_date + '.bpch'
    next_files = gc_nd49_dir + 'ts' + n_date + '.bpch'

    # read model data
    model_files = [pre_files, curr_files, next_files]
    model_data = \
            read_nd49_resample(model_files, mod_varname_list)

    # unit conversion (ppbv => molec/cm2)
    # at every layer
    # and prepare model data for resmapling
    mod_var_dict = {}
    mod_NO2 = vmr_to_DU(model_data[mod_spename],
            model_data['TIME-SER_AIRDEN'], model_data['BXHGHT-$_BXHEIGHT'])
    mod_var_dict['NO2'] = mod_NO2

    # preprae additional data for resampling
    mod_TAI93 = model_data['TAI93']
    mod_var_dict['PEdge_Bot'] = model_data['PEDGE-$_NOx']

    # find satellite files
    sat_wildcard = sat_dir + 'OMPS-NPP_NMNO2-L2_' + currDate[0:4] \
            + 'm' + currDate[5:7] + currDate[8:10] + '_o*.h5'
    all_sat_files = glob.glob(sat_wildcard)
    all_sat_files.sort()
    print('All satellite files on ' + currDate)
    for i in range(len(all_sat_files)):
        print('  ' + all_sat_files[i])

    # process satellite files
    print('Process satellites on ' + currDate)
    #for i in range(len(all_sat_files)):
    for i in [3]:

        # read satellite file
        sat_file = all_sat_files[i]
        print('  reading ' + sat_file)
        sat_data = read_OMPS_NO2_L2(sat_file, verbose=verbose)

        # quality control



    # go to next day
    currDate_D = currDate_D + datetime.timedelta(days=1)
