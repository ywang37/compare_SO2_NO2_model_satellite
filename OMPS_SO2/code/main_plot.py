"""
Created on Feburary 3, 2020

@author: Yi Wang
"""

import datetime
import copy
import glob
import matplotlib.pyplot as plt
import numpy as np
import os
import sys

from mylib.io import read_nc
sys.path.append("/Dedicated/jwang-data/ywang/\
compare_SO2_NO2_model_satellite/shared_code/")
from com_plot import compare_two_varilables


#######################
# Start user parameters
#


monthly_dir = '../data/monthly/'

startDate = '2018-04-01'
endDate   = '2018-04-30'

varname1 = 'sat_ColumnAmountSO2_PBL'
varname2 = varname1 + '_correction'

region_limit = [24.0, 92.5, 54.0, 130.0]

diff_min = -0.15
diff_max = -diff_min

vmax=0.6

sc_xlabel=r'OMPS SO$_2$ [DU]'
sc_ylabel=r'OMPS SO$_2$ (AK) [DU]'

filename = monthly_dir + 'model_satellite_' + startDate \
        + '_' + endDate + '.nc'

figname = '../figure/SO2_VCD_' + startDate + '_' + endDate + '.png'

#
# End user parameters
#####################



compare_two_varilables(filename, varname1, varname2,
        diff_min=diff_min, diff_max=diff_max,
        vmax=vmax,
        sc_xlabel=sc_xlabel, sc_ylabel=sc_ylabel,
        region_limit=region_limit)

plt.savefig(figname, format='png', dpi=300)
