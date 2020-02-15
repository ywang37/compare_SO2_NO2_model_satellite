"""
Created on Feburary 3, 2020

@author: Yi Wang
"""

from copy import deepcopy
import matplotlib.pyplot as plt
import numpy as np

from mylib.borders.China.china import add_China_province
from mylib.cartopy_plot import add_geoaxes
from mylib.cartopy_plot import cartopy_plot
from mylib.colormap.WhGrYlRd_map import WhGrYlRd_map
from mylib.grid_utility import get_center_index
from mylib.io import read_nc
from mylib.layout import h_1_ax, h_2_ax, right_center_label
from mylib.scatter_plot import scatter

def compare_two_varilables(filename, varname1, varname2,
        vmin=0.0, vmax=0.8,
        diff_min=None, diff_max=None,
        seperate_cbar=False,
        region_limit=None,
        sc_xlabel='',
        sc_ylabel='',
        sc_max=1.0,
        verbose=True):
    """
    """

    # read data
    coord_varns = ['Latitude_e', 'Longitude_e']
    some_varns = [varname1, varname2]
    varnames = coord_varns + some_varns
    data_dict = read_nc(filename, varnames, verbose=verbose)


    # region_limit
    if region_limit is not None:
        
        lat_e_1D = data_dict['Latitude_e'][:,0]
        lon_e_1D = data_dict['Longitude_e'][0,:]

        # index 
        i1 = get_center_index(lat_e_1D, region_limit[0])
        i2 = get_center_index(lat_e_1D, region_limit[2])
        j1 = get_center_index(lon_e_1D, region_limit[1])
        j2 = get_center_index(lon_e_1D, region_limit[3])

        for varn in varnames:
            if varn in coord_varns:
                data_dict[varn] = data_dict[varn][i1:i2+2,j1:j2+2]
            if varn in some_varns:
                data_dict[varn] = data_dict[varn][i1:i2+1,j1:j2+1]
        
    # plot
    fig = plt.figure(figsize=(8,7))
    plt.subplots_adjust(top=0.98)
    ax_list = []
    xtick = [100, 110, 120]
    ytick = [30, 40, 50]
    xtick_list = [xtick, xtick, xtick]
    ytick_list = [ytick, [], ytick]
    for i in range(len(xtick_list)):
        ax = add_geoaxes(fig, int('22'+str(i+1)),
                xtick=xtick_list[i], ytick=ytick_list[i])
        ax_list.append(ax)

    lat_e = data_dict['Latitude_e']
    lon_e = data_dict['Longitude_e']

    # colorbar height
    h = 0.02

    # variables
    for i in range(len(some_varns)):

        varn = some_varns[i]
        var = data_dict[varn]
        ax = ax_list[i]
        var_out = cartopy_plot(lon_e, lat_e, var, ax=ax,
                vmin=vmin, vmax=vmax,
                cbar=seperate_cbar, cmap=deepcopy(WhGrYlRd_map))
        add_China_province(ax)

    # colorbar for variables
    ax1 = ax_list[0]
    ax2 = ax_list[1]
    cax_v = h_2_ax(fig, ax1, ax2, y_off=-0.06, height=h)
    plt.colorbar(var_out['mesh'], cax=cax_v, orientation='horizontal')
    right_center_label(cax_v, '[DU]')


    # difference
    ax = ax_list[2]
    var_diff = data_dict[varname2] - data_dict[varname1]
    diff_out = cartopy_plot(lon_e, lat_e, var_diff, ax=ax,
            vmin=diff_min, vmax=diff_max,
            cbar=seperate_cbar, cmap=plt.get_cmap('seismic'))
    add_China_province(ax)

    # colorbar of difference
    ax = ax_list[2]
    cax_diff = h_1_ax(fig, ax, y_off=-0.06, height=h)
    plt.colorbar(diff_out['mesh'], cax=cax_diff, orientation='horizontal')
    right_center_label(cax_diff, '[DU]')

    # scatter plot
    ax_sc = fig.add_subplot(224)
    ax_sc.set_aspect('equal')
    ax_sc.set_xlabel(sc_xlabel)
    ax_sc.set_ylabel(sc_ylabel)
    ax_sc.set_xlim([0,1.0])
    ax_sc.set_ylim([0,1.0])
    label_ul = ['R', 'linear_eq', 'rmse', 'nmb', 'mb', 'N']
    scatter(ax_sc, data_dict[varname1], data_dict[varname2], s=3,
            label_ul=label_ul)


    # region
    if region_limit is not None:
        for ax in ax_list:
            ax.set_xlim(region_limit[1], region_limit[3])
            ax.set_ylim(region_limit[0], region_limit[2])





