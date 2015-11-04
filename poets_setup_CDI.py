# Copyright (c) 2014, Vienna University of Technology (TU Wien), Department
# of Geodesy and Geoinformation (GEO).
# All rights reserved.

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
#
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
#
# * Neither the name of the Vienna University of Technology - Department of
#   Geodesy and Geoinformation nor the names of its contributors may be used to
#   endorse or promote products derived from this software without specific
#   prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL VIENNA UNIVERSITY OF TECHNOLOGY,
# DEPARTMENT OF GEODESY AND GEOINFORMATION BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
# OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
# EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# Author: Thomas Mistelbauer thomas.mistelbauer@geo.tuwien.ac.at
# Creation date: 2014-11-05

"""
Setup of poets for the SATIDA project.
"""

import os
import sys
import getopt
from datetime import datetime
from cdi_poets import CDIPoet

if __name__ == "__main__":

    # POETS SETTINGS
    spatial_resolution = 0.25
    temporal_resolution = 'dekad'
    rootpath = '/home/tuwien/poets'
    nan_value = -99
    start_date = datetime(1992, 1, 1)
    regions = ['ET', 'CT', 'SG']
    region_names = ['Ethiopia', 'Central African Republic', 'Senegal']
    delete_rawdata = True
    valid_range = (0, 1)
    ip = [9, 18]
    refparam = 'Vegetation_Status_dataset_DI'

    staticsources = ['ECDI_9', 'ECDI_18', 'Vegetation_Status_DI_9',
                     'Vegetation_Status_DI_18',
                     'ECV_DI_9', 'ECV_DI_18', 'MODIS_LST_DI_9',
                     'MODIS_LST_DI_18', 'TAMSAT_DI_9', 'TAMSAT_DI_18',
                     'WARNING_LEVELS_ECDI_9',
                     'WARNING_LEVELS_ECDI_18']

    p = CDIPoet(ip, refparam, staticsources, rootpath, regions,
                spatial_resolution, temporal_resolution, start_date, nan_value,
                delete_rawdata=delete_rawdata, region_names=region_names)

    # SOURCE SETTINGS TAMSAT
    name = 'TAMSAT'
    filename = "rfe{YYYY}_{MM}-dk{P}.nc"
    filedate = {'YYYY': (3, 7), 'MM': (8, 10), 'P': (13, 14)}
    temp_res = 'dekad'
    host = "http://www.met.reading.ac.uk"
    protocol = 'HTTP'
    directory = '~tamsat/public_data'
    dirstruct = ['YYYY', 'MM']
    begin_date = datetime(1983, 1, 11)
    colorbar = 'jet_r'

    p.add_source(name, filename, filedate, temp_res, host, protocol,
                 directory=directory, dirstruct=dirstruct,
                 begin_date=begin_date, colorbar=colorbar)

    # TAMSAT_DI
    for ipe in ip:
        src = 'TAMSAT'
        name = src + '_DI_' + str(ipe)
        filename = 'nix'
        filedate = {}
        temp_res = 'dekad'
        host = 'nix'
        protocol = 'FTP'
        colorbar = 'jet_r'
        begin_date = datetime(1992, 1, 1)
        variables = []
        for var in p.sources[src].get_variables():
            variables.append(var + '_DI_' + str(ipe))
        valid_range = (0, 1)
        src_file = {}
        labels = ['     deficit', 'surplus      .']
        xticks = [0, 1]
        unit = 'Rainfall'
        for reg in regions:
            src_file[reg] = os.path.join(p.di_path, reg + '_' + src +
                                         '_DI_' + str(ipe) + '.nc')

        p.add_source(name, filename, filedate, temp_res, host, protocol,
                     colorbar=colorbar, variables=variables,
                     begin_date=begin_date, valid_range=valid_range,
                     src_file=src_file, xticks=xticks, labels=labels,
                     unit=unit)

    # SOURCE SETTINGS BOKU NDVI
    name = 'Vegetation_Status'
    filename = "HOA_MCD13Q1.A{{YYYY}}{{RP}}.10_days.NDVIgd.tif"
    filedate = {'YYYY': (13, 17), 'RP': (17, 19)}
    temp_res = 'dekad'
    host = "ivfl-rio.boku.ac.at"
    protocol = 'FTP'
    directory = "/SATIDA/MODIS_NRT/HOA/"
    begin_date = datetime(2002, 1, 1)
    nan_value = 255
    data_range = (0, 250)
    valid_range = (-0.2, 1)
    unit = "NDVI"
    colorbar = 'BrBG'
    regs = ['ET', 'CT']

    p.add_source(name, filename, filedate, temp_res, host, protocol,
                 directory=directory, begin_date=begin_date,
                 nan_value=nan_value, valid_range=valid_range,
                 data_range=data_range, unit=unit, colorbar=colorbar,
                 regions=regs)

    # NDVI_DI
    for ipe in ip:
        src = 'Vegetation_Status'
        name = src + '_DI_' + str(ipe)
        filename = 'nix'
        filedate = {}
        temp_res = 'dekad'
        host = 'nix'
        protocol = 'FTP'
        colorbar = 'jet_r'
        begin_date = datetime(1992, 1, 1)
        variables = []
        for var in p.sources[src].get_variables():
            variables.append(var + '_DI_' + str(ipe))
        valid_range = (0, 1)
        src_file = {}
        labels = ['     bad', 'good      .']
        xticks = [0, 1]
        unit = 'Vegetation health'
        for reg in regs:
            src_file[reg] = os.path.join(p.di_path, reg + '_' + src +
                                         '_DI_' + str(ipe) + '.nc')

        p.add_source(name, filename, filedate, temp_res, host, protocol,
                     colorbar=colorbar, variables=variables,
                     begin_date=begin_date, valid_range=valid_range,
                     src_file=src_file, xticks=xticks, labels=labels,
                     unit=unit, regions=regs)

    # SOURCE SETTINGS MODIS LST
    name = 'MODIS_LST'
    filename = "MOD11C1_D_LSTDA_{YYYY}_{MM}-{DD}.png"
    filedate = {'YYYY': (16, 20), 'MM': (21, 23), 'DD': (24, 26)}
    temp_res = 'daily'
    host = "neoftp.sci.gsfc.nasa.gov"
    protocol = 'FTP'
    directory = "/gs/MOD11C1_D_LSTDA/"
    begin_date = datetime(2000, 2, 24)
    nan_value = 255
    data_range = (0, 254)
    valid_range = (-25, 45)
    unit = "degree Celsius"

    p.add_source(name, filename, filedate, temp_res, host, protocol,
                 directory=directory, begin_date=begin_date,
                 nan_value=nan_value, valid_range=valid_range,
                 data_range=data_range, unit=unit)

    # MODIS_LST_DI
    for ipe in ip:
        src = 'MODIS_LST'
        name = src + '_DI_' + str(ipe)
        filename = 'nix'
        filedate = {}
        temp_res = 'dekad'
        host = 'nix'
        protocol = 'FTP'
        colorbar = 'jet_r'
        begin_date = datetime(1992, 1, 1)
        variables = []
        for var in p.sources[src].get_variables():
            variables.append(var + '_DI_' + str(ipe))
        valid_range = (0, 1)
        src_file = {}
        labels = ['     high', 'low    .']
        xticks = [0, 1]
        unit = 'Temperature anomaly'
        for reg in regions:
            src_file[reg] = os.path.join(p.di_path, reg + '_' + src +
                                         '_DI_' + str(ipe) + '.nc')

        p.add_source(name, filename, filedate, temp_res, host, protocol,
                     colorbar=colorbar, variables=variables,
                     begin_date=begin_date, valid_range=valid_range,
                     src_file=src_file, xticks=xticks, labels=labels,
                     unit=unit)

    # SOURCE SETTINGS ECV
    name = 'ECV'
    temp_res = 'daily'
    host = 'ftp.ipf.tuwien.ac.at'
    # =========================================================================
    # #Data Archive:
    #directory = "_down/daily_files/COMBINED/"
    #dirstruct = ['YYYY']
    #protocol = 'SFTP'
    #username = 'esacci_sm_v022'
    #password = '5ZVJuN4rxsWD'
    #filename = ("ESACCI-SOILMOISTURE-L3S-SSMV-COMBINED-{YYYY}{MM}{TT}{hh}"
    #            "{mm}{ss}-fv02.1.nc")
    #filedate = {'YYYY': (38, 42), 'MM': (42, 44), 'DD': (44, 46),
    #            'hh': (46, 48), 'mm': (48, 50), 'ss': (50, 52)}
    # =========================================================================
    # NRT data:
    #==========================================================================
    directory = "_down/CCI_ASCAT_AMSR2_merged"
    protocol = 'SFTP'
    username = 'crts'
    password = 'Ibonaguma343'
    filename = ("merged_NRT_{YYYY}{MM}{TT}.nc")
    filedate = {'YYYY': (11, 15), 'MM': (15, 17), 'DD': (17, 19)}
    #==========================================================================
    nan_value = -9999
    port = 22
    begin_date = datetime(1992, 1, 01)
    variables = ['sm']
    valid_range = (0, 0.6)
    colorbar = 'jet_r'

    p.add_source(name, filename, filedate, temp_res, host, protocol,
                 username=username, password=password, port=port,
                 directory=directory, begin_date=begin_date,
                 variables=variables, colorbar=colorbar,
                 valid_range=valid_range, nan_value=nan_value)#,)
                 #dirstruct=dirstruct)

    staticsources = []

    # ECV_DI
    for ipe in ip:
        src = 'ECV'
        name = src + '_DI_' + str(ipe)
        filename = 'nix'
        filedate = {}
        temp_res = 'dekad'
        host = 'nix'
        protocol = 'FTP'
        colorbar = 'jet_r'
        begin_date = datetime(1992, 1, 1)
        variables = []
        for var in p.sources[src].get_variables():
            variables.append(var + '_DI_' + str(ipe))
        valid_range = (0, 1)
        src_file = {}
        labels = ['     deficit', 'surplus      .']
        xticks = [0, 1]
        unit = 'Soil Moisture'
        for reg in regions:
            src_file[reg] = os.path.join(p.di_path, reg + '_' + src +
                                         '_DI_' + str(ipe) + '.nc')

        p.add_source(name, filename, filedate, temp_res, host, protocol,
                     colorbar=colorbar, variables=variables,
                     begin_date=begin_date, valid_range=valid_range,
                     src_file=src_file, xticks=xticks, labels=labels,
                     unit=unit)

    for ipe in ip:
        # CDI parameter
        name = 'ECDI_' + str(ipe)
        filename = 'nix'
        filedate = {}
        temp_res = 'dekad'
        host = 'nix'
        protocol = 'FTP'
        colorbar = 'jet_r'
        begin_date = datetime(1992, 1, 1)
        variables = ['ECDI_' + str(ipe)]
        valid_range = (0, 1)
        src_file = {}
        xticks = [0, 1]
        labels = ['     high', 'low     .']
        unit = 'Drought risk'

        for reg in regs:
            src_file[reg] = os.path.join(p.cdi_path, reg + '_ECDI_' + str(ipe)
                                         + '.nc')

        p.add_source(name, filename, filedate, temp_res, host, protocol,
                     colorbar=colorbar, variables=variables,
                     begin_date=begin_date, valid_range=valid_range,
                     src_file=src_file, labels=labels, xticks=xticks,
                     unit=unit, regions=regs)

        # WARNING LEVELS
        name = 'WARNING_LEVELS_ECDI_' + str(ipe)
        filename = 'nix'
        filedate = {}
        temp_res = 'dekad'
        host = 'nix'
        protocol = 'FTP'
        colorbar = 'RdYlGn_r'
        begin_date = datetime(1992, 1, 1)
        valid_range = (0, 3)
        src_file = {}
        xticks = [0, 1, 2, 3]
        # blanks in labels are for formatting reasons
        labels = ['     normal', 'increased', 'severe', 'extreme      .']

        for reg in regs:
            src_file[reg] = os.path.join(p.cdi_path, reg + '_WARNING_LEVELS_' +
                                         'ECDI_' + str(ipe) + '.nc')

        p.add_source(name, filename, filedate, temp_res, host, protocol,
                     colorbar=colorbar, labels=labels, xticks=xticks,
                     begin_date=begin_date, valid_range=valid_range,
                     src_file=src_file, regions=regs)

    argv = sys.argv[1:]

    try:
        opts, args = getopt.getopt(argv, 'ho:', ['help', 'action='])
    except:
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            print('-o [start_app, fetch_data, calc_weights, calc_di, calc_cdi,'
                  'fill_gaps]')
        elif opt == '-o':
            if arg == 'start_app':
                p.start_app('192.168.122.81', port=80, r_host='148.251.42.101',
                            r_port=8580, debug=True)
            elif arg == 'fetch_data':
                p.fetch_data()
            elif arg == 'calc_weights':
                p.WeightsToNetCDF(refparam, region=regs,
                                  exclude='Vegetation_Status')
            elif arg == 'calc_di':
                p.DItoNetCDF()
            elif arg == 'calc_cdi':
                p.CDItoNetCDF(region=regs, exclude='Vegetation_Status')
            elif arg == 'fill_gaps':
                p.fill_gaps()
            elif arg == 'custom':
                p.sources['ECV'].download_and_resample(begin=datetime(2014, 3, 21), 
                                                       end=datetime(2014, 12, 31),
                                                       delete_rawdata=True)
            else:
                print 'unknown argument'
        else:
            print 'please give an argument'

