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

import matplotlib
matplotlib.use('Agg')
from datetime import datetime
from poets.poet import Poet

if __name__ == "__main__":

    # POETS SETTINGS
    spatial_resolution = 0.25
    temporal_resolution = 'dekad'
    rootpath = '/home/tuwien/poets'
    nan_value = -99
    start_date = datetime(1992, 1, 1)
    regions = ['CT']
    delete_rawdata = True
    valid_range = (0, 1)

    p = Poet(rootpath, regions, spatial_resolution, temporal_resolution,
             start_date, nan_value, delete_rawdata=delete_rawdata)

    # SOURCE SETTINGS BOKU NDVI
    name = 'Vegetation_Status'
    filename = "CF_MCD13Q1.A{{YYYY}}{{RP}}.10_days.NDVIgd.tif"
    filedate = {'YYYY': (12, 16), 'RP': (16, 18)}
    temp_res = 'dekad'
    host = "ivfl-rio.boku.ac.at"
    protocol = 'FTP'
    directory = "/SATIDA/MODIS_NRT/CF/"
    begin_date = datetime(2015, 1, 1)
    nan_value = 255
    data_range = (0, 254)
    valid_range = (0, 1)
    unit = "NDVI"
    colorbar = 'Greens'

    p.add_source(name, filename, filedate, temp_res, host, protocol,
                 directory=directory, begin_date=begin_date,
                 nan_value=nan_value, valid_range=valid_range,
                 data_range=data_range, unit=unit, colorbar=colorbar)

    p.sources['Vegetation_Status'].download_and_resample(begin=begin_date, 
                                                 delete_rawdata=delete_rawdata)
