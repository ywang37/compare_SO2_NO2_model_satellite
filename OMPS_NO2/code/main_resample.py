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
from mylib.pro_omps_no2_l2.pro_omps_no2_l2 import QC_OMPS_NO2_L2
from mylib.pro_satellite.sat_model_sample import sat_model_sample
from mylib.pro_satellite.sat_model_sample import save_sat_model_sample

#######################
# Start user parameters
#

gc_root_dir = '/Dedicated/jwang-data/ywang/OMPS_NO2_inverse/GC_adj/\
gcadj_std/runs/v8-02-01/geos5_201310_pressure/'

sat_dir = '/Dedicated/jwang-data/shared_satData/OMPS/NMNO2-L2/AK201310/'

out_dir = '../data/granule/'

gc_trop_file = gc_root_dir + 'ctm.bpch'

gc_nd49_dir = gc_root_dir + 'ND49.1/'

startDate = '2013-10-08'
endDate   = '2013-10-08'

# species names
mod_spename = 'TIME-SER_NO2'

# all model variable names
mod_varname_list = [mod_spename, 'TIME-SER_AIRDEN', \
        'BXHGHT-$_BXHEIGHT', 'PEDGE-$_NOx']

# mod_coord_dict
mod_coord_dict = {}
mod_coord_dict['coord_format']  = 'ses'
mod_coord_dict['mod_lat_start'] =    9.0
mod_coord_dict['mod_lat_end']   =   61.0
mod_coord_dict['mod_lat_step']  =    2.0
mod_coord_dict['mod_lon_start'] =   88.75
mod_coord_dict['mod_lon_end']   =  151.25
mod_coord_dict['mod_lon_step']  =    2.5

# Satellite variables that need to be vertically reversed
sat_rev_varname_list = ['AveragingKernel',
        'PressureLevel', 'TROPO_NO2_ShapeFactor']

sat_aprior_varname_list = [
        '/aPriori/AveragingKernel',
        '/aPriori/PressureLevel',
        '/aPriori/TROPO_NO2_ShapeFactor',
        '/aPriori/nlyrs',
        ]


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
    #for i in [3]:
    for i in [2,3]:

        # read satellite file
        sat_file = all_sat_files[i]
        print('  reading ' + sat_file)
        sat_data = read_OMPS_NO2_L2(sat_file, varnames=sat_aprior_varname_list,
                verbose=verbose)

        # In the OMPS dataset, layer order is from top to bottom
        # We reverse the layer order, hence it becomes from bottom to top
        for sat_rev_varname in sat_rev_varname_list:
            # reverse vertically
            sat_data[sat_rev_varname] = sat_data[sat_rev_varname][:,:,::-1]
            # assign unavailable valule as np.nan
            sat_data[sat_rev_varname][sat_data[sat_rev_varname]<-99999.] = \
                    np.nan


        # quality control
        sat_NO2_trop = sat_data['ColumnAmountNO2tropo']
        sza          = sat_data['SolarZenithAngle']
        vza          = sat_data['ViewingZenithAngle']
        CF           = sat_data['CloudFraction']
        QF_pixel     = sat_data['PixelQualityFlags']
        sat_flag = QC_OMPS_NO2_L2(sat_NO2_trop, sza, vza, CF, QF_pixel)

        # prepare satellite data for resmapling
        sat_lat = sat_data['Latitude']
        sat_lon = sat_data['Longitude']
        sat_TAI93 = np.tile(sat_data['ImageMidpoint_TAI93'], sat_lat.shape[1])
        sat_TAI93 = sat_TAI93.reshape(sat_lat.shape[::-1])
        sat_TAI93 = sat_TAI93.T
        sat_obs_dict = {}
        sat_obs_dict['ColumnAmountNO2tropo'] = sat_data['ColumnAmountNO2tropo']
        sat_obs_dict['AveragingKernel']      = sat_data['AveragingKernel']
        sat_obs_dict['PressureLevel']        = sat_data['PressureLevel']
        sat_obs_dict['TROPO_NO2_ShapeFactor'] = \
                sat_data['TROPO_NO2_ShapeFactor']
        sat_obs_dict['nlyrs']                = sat_data['nlyrs']

        # Sample model results according satellite observations
        # and regrid satellite observations to model grids.
        sat_mod_dict = \
                sat_model_sample(mod_coord_dict, mod_TAI93, mod_var_dict,
                sat_lat, sat_lon, sat_TAI93, sat_obs_dict,
                sat_flag=sat_flag)

        # save data
        lon, lat = np.meshgrid(model_data['longitude'], model_data['latitude'])
        sat_mod_dict['Latitude']  = lat
        sat_mod_dict['Longitude'] = lon
        lon_e, lat_e = \
                np.meshgrid(model_data['longitude_e'],
                        model_data['latitude_e'])
        sat_mod_dict['Latitude_e']   = lat_e
        sat_mod_dict['Longitude_e']  = lon_e
        out_file = out_dir + 'model_satellite_' + \
                sat_file.split('/')[-1][18:34] + '.nc'
        save_sat_model_sample(out_file, sat_mod_dict)


    # go to next day
    currDate_D = currDate_D + datetime.timedelta(days=1)
