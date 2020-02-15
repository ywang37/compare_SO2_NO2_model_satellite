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

startDate = '2013-10-01'
endDate   = '2013-10-31'

varname1 = 'sat_ColumnAmountNO2tropo'
varname2 = varname1 + '_correction'

region_limit = [24.0, 92.5, 54.0, 130.0]

diff_min = -0.12
diff_max = -diff_min

sc_xlabel=r'OMPS NO$_2$ [DU]'
sc_ylabel=r'OMPS NO$_2$ (AK) [DU]'

filename = monthly_dir + 'model_satellite_' + startDate \
        + '_' + endDate + '.nc'

figname = '../figure/NO2_VCD_' + startDate + '_' + endDate + '.png'

#
# End user parameters
#####################



compare_two_varilables(filename, varname1, varname2,
        diff_min=diff_min, diff_max=diff_max,
        sc_xlabel=sc_xlabel, sc_ylabel=sc_ylabel,
        region_limit=region_limit)

plt.savefig(figname, format='png', dpi=300)
