# Copyright (c) 2015, Vienna University of Technology (TU Wien), Department
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

# Author: Thomas Mistelbauer Thomas.Mistelbauer@geo.tuwien.ac.at
# Creation date: 2015-01-16

"""
POETS extension for Drought Monitoring.
"""

import os
import cdi
import numpy as np
import pandas as pd
from netCDF4 import Dataset, date2num, num2date
from poets.poet import Poet
from poets.grid import grids
from poets.timedate.dateindex import get_dtindex
from pytesmo.grid.netcdf import save_grid


class CDIPoet(Poet):
    """
    Extends Poet base class with CDI specific parameters and methods.

    Parameters
    ----------
    ip : int, list of int
        Interest periods for calculating the CDI.
    refparam : str
        Reference parameter for calculating the CDI.
    lags : list of int, optional
        Time periods to shift parameter against refparam, defaults to [0, 10].

    Attributes
    ----------
    ip : list of int
        Interest periods for calculating the CDI.
    refparam : str
        Reference parameter for calculating the CDI.
    lags : list of int
        Time periods to shift parameter against refparam.
    cdi_sources : dict of poets.io.BasicSource objects
        Sources for CDI parameters used by poets given as BasicSource class.
    """

    def __init__(self, ip, refparam, staticsources, *args, **kwargs):
        super(CDIPoet, self).__init__(*args, **kwargs)

        if isinstance(ip, int):
            self.ip = [ip]
        else:
            self.ip = ip

        self.lags = [0, 10]
        self.refparam = refparam
        self.di_path = os.path.join(self.data_path, 'DI')
        self.weights_path = os.path.join(self.data_path, 'Weights')
        self.cdi_path = os.path.join(self.data_path, 'CDI')
        self.cdi_sources = {}
        self.staticsources = staticsources

    def fetch_data(self, begin=None, end=None, delete_rawdata=None):
        """Starts download and resampling of input data for sources as added
        to `Poets.sources`.

        Parameters
        ----------
        begin : datetime, optional
            Start date of data to download, defaults to start date as defined
            in poets class.
        end : datetime, optional
            End date of data to download, defaults to current datetime.
        delete_rawdata : bool, optional
            Original files will be deleted from rawdata_path if set True.
            Defaults to value of delete_rawdata attribute as set in Poet class.
        """

        if delete_rawdata is None:
            delete_rawdata = self.delete_rawdata

        for source in self.sources.keys():
            if source in self.staticsources:
                continue
            src = self.sources[source]
            print '[INFO] Download data for source ' + source
            src.download_and_resample(begin=begin, end=end,
                                      shapefile=self.shapefile,
                                      delete_rawdata=delete_rawdata)

        print '[SUCCESS] Download and resampling complete!'

    def fill_gaps(self):
        """
        Detects gaps in data and tries to fill them by downloading and
        resampling the data within these periods.
        """

        for source in self.sources.keys():
            if source in self.staticsources:
                continue
            src = self.sources[source]
            print '[INFO] Scanning ' + source + ' for gaps'
            src.fill_gaps()

    def DItoNetCDF(self, region=None, source=None, ip=None):
        """
        Calculates the Drought Index (DI) for a given source over one or more
        regions and stores them as NetCDF files.

        Parameters
        ----------
        region : str, list of str, optional
            Region(s) of of interest; must be one of the regions as set in the
            CDIPoet instance; Defaults the regions attribute value of the
            CDIPoet instance.
        source : str, list of str, optional
            Source parameter(s) for which to calculate the DI; must be one of
            the sources set in the CDIPoet instance; Defaults to all
            sources set in the CDIPoet instance.
        ip : int, list of int, optional
            Interest period for calculating the DI; must be one of the ip
            as set in the CDIPoet instance; Defaults to the ip  attribute value
            in the CDIPoet instance.
        """

        if region is None:
            region = self.regions
        else:
            if isinstance(region, str):
                region = [region]

        if source is None:
            source = self.sources.keys()
        else:
            if isinstance(source, str):
                source = [source]

        if ip is None:
            ip = self.ip
        else:
            if isinstance(ip, int):
                ip = [ip]

        if not os.path.exists(self.di_path):
            os.mkdir(self.di_path)

        for reg in region:

            grid = grids.ShapeGrid(reg, self.spatial_resolution)
            gps = grid.get_gridpoints().index

            for ipe in ip:
                for src in source:
                    if src in self.staticsources:
                        continue
                    if reg not in self.sources[src].valid_regions:
                        continue
                    print ('[INFO] calc DI ' + reg + ' IP' + str(ipe) + ' ' +
                           src),
                    self.__writeDI(reg, src, gps, grid, ipe)
                    print ' done!'

    def __writeDI(self, region, src, gridpoints, grid, ip, suffix='',
                  scaled=True, modf_all=True, start=None):

        if start is not None:
            dt = get_dtindex('dekad', start)
        else:
            dt = get_dtindex('dekad', self.start_date)

        dest_file = os.path.join(self.di_path, region + '_' + src + '_DI'
                                 + '_' + str(ip) + '.nc')

        if not os.path.isfile(dest_file):
            save_grid(dest_file, grid)

        for i, gp in enumerate(gridpoints):

            if i % 100 == 0:
                print '.',

            ts = self.read_timeseries(src, gp, region)
            if start is not None:
                sel = (ts.index >= start)
                ts = ts[sel]

            inverse = False
            if src == 'MODIS_LST':
                inverse = True

            ts_di = cdi.calc_DI(ts.copy(), inverse, [ip], scale_zero=False,
                                scaled=scaled, modf_all=modf_all)

            with Dataset(dest_file, 'r+', format='NETCDF4') as nc:
                if 'time' not in nc.dimensions.keys():
                    nc.createDimension("time", None)

                    times = nc.createVariable('time', 'uint16', ('time',))

                    times.units = 'days since ' + str(self.start_date)
                    times.calendar = 'standard'
                    times[:] = date2num(dt.tolist(), units=times.units,
                                        calendar=times.calendar)

                else:
                    times = nc.variables['time']

                dim = ('time', 'lat', 'lon')

                position = np.where(nc.variables['gpi'][:] == gp)
                lat_pos = position[0][0]
                lon_pos = position[1][0]

                # extend times variable in NetCDF
                tsdates = date2num(ts_di.index.tolist(), units=times.units,
                                   calendar=times.calendar).astype(int)
                begin = np.where(times == tsdates[0])[0][0]
                times[begin:] = tsdates

                for dataset in ts_di.keys():

                    if dataset not in nc.variables.keys():
                        var = nc.createVariable(dataset,
                                                ts_di[dataset].dtype.char,
                                                dim, fill_value=self.nan_value)
                    else:
                        var = nc.variables[dataset]

                    var[begin:, lat_pos, lon_pos] = ts_di[dataset].values

    def WeightsToNetCDF(self, refparam, region=None, ip=None, exclude=None):
        """
        Parameters
        ----------
        exclude : string, optional
            Variable which should not be used for calculation of the weights.
        """

        if region is None:
            region = self.regions
        else:
            if isinstance(region, str):
                region = [region]

        if ip is None:
            ip = self.ip
        else:
            if isinstance(ip, int):
                ip = [ip]

        if not os.path.exists(self.weights_path):
            os.mkdir(self.weights_path)

        for reg in region:

            grid = grids.ShapeGrid(reg, self.spatial_resolution)
            gps = grid.get_gridpoints().index

            for ipe in ip:
                print '[INFO] calc weights ' + reg + ' IP' + str(ipe),
                for i, gp in enumerate(gps):
                    if i % 100 == 0:
                        print '.',
                    self.__writeWeight(gp, reg, refparam, ipe, exclude)
                print ' done!'

    def __writeWeight(self, gp, region, refparam, ip, exclude=None):
        """
        Parameters
        ----------
        exclude : string, optional
            Variable which should not be used for calculation of the weights.
        """

        refparam += '_' + str(ip)

        df = pd.DataFrame()

        for param in self.sources.keys():

            difile = os.path.join(self.di_path,
                                  region + '_' + param + '_DI_' + str(ip) +
                                  '.nc')

            if not os.path.exists(difile):
                continue

            with Dataset(difile, 'r', format='NETCDF4') as nc:
                if len(df.index.values) == 0:
                    time = nc.variables['time']
                    dates = num2date(time[:], units=time.units,
                                     calendar=time.calendar)
                    df = pd.DataFrame(index=pd.DatetimeIndex(dates))

                ncvar = None
                for var in nc.variables.keys():
                    if param in var:
                        ncvar = var
                        continue

                position = np.where(nc.variables['gpi'][:] == gp)
                lat_pos = position[0][0]
                lon_pos = position[1][0]

                df[ncvar] = np.NAN
                for i in range(0, nc.variables[ncvar].shape[0] - 1):
                    df[ncvar][i] = nc.variables[ncvar][i, lat_pos, lon_pos]

                    if 'scaling_factor' in nc.variables[ncvar].ncattrs():
                        vvar = nc.variables[ncvar]
                        if vvar.getncattr('scaling_factor') < 0:
                            df[ncvar] = (df[ncvar] *
                                         float(vvar.getncattr('scaling_factor')))
                        else:
                            df[ncvar] = (df[ncvar] /
                                         float(vvar.getncattr('scaling_factor')))

        weights = cdi.calc_weights(df, refparam, lags=self.lags,
                                   exclude=exclude)

        dest_file = os.path.join(self.weights_path,
                                 region + '_weights_' + str(ip) + '.nc')

        if not os.path.isfile(dest_file):
            grid = grids.ShapeGrid(region, self.spatial_resolution)
            save_grid(dest_file, grid)

        with Dataset(dest_file, 'r+', format='NETCDF4') as nc:
            dim = ('lat', 'lon')

            position = np.where(nc.variables['gpi'][:] == gp)
            lat_pos = position[0][0]
            lon_pos = position[1][0]

            keys = []
            if exclude is not None:
                for par in df.keys():
                    if exclude in par:
                        continue
                    keys.append(par)
            else:
                keys = df.keys()

            for i, dataset in enumerate(keys):

                if dataset not in nc.variables.keys():
                    var = nc.createVariable(dataset, 'd', dim,
                                            fill_value=self.nan_value)
                else:
                    var = nc.variables[dataset]

                var[lat_pos, lon_pos] = weights[i]

    def CDItoNetCDF(self, region=None, ip=None, separatefile=True,
                    exclude=None):
        """
        Creates NetCDF that contains CDI for all timestamps.

        Parameters
        ----------
        region : str, list of str, optional
            Region(s) of of interest; must be one of the regions as set in the
            CDIPoet instance; Defaults the regions attribute value of the
            CDIPoet instance.
        ip : int, list of int, optional
            Interest period for calculating the DI; must be one of the ip
            as set in the CDIPoet instance; Defaults to the ip  attribute value
            in the CDIPoet instance.
        separatefile : bool
            If True, writes weights to separate file; If False, writes weights
            to NetCDF database file.
        exclude : string, optional
            Variable which should not be used for calculation of CDI.
        """

        if region is None:
            region = self.regions
        else:
            if isinstance(region, str):
                region = [region]

        if ip is None:
            ip = self.ip
        else:
            if isinstance(ip, int):
                ip = [ip]

        if not os.path.exists(self.cdi_path):
            os.mkdir(self.cdi_path)

        for reg in region:
            grid = grids.ShapeGrid(reg, self.spatial_resolution)
            gps = grid.get_gridpoints().index

            for ipe in ip:
                key = 'ECDI_' + str(ipe)

                print ('[INFO] calc ECDI ' + reg + ' IP' + str(ipe))

                if separatefile:
                    dest_file = os.path.join(self.cdi_path,
                                             reg + '_' + key + '.nc')
                else:
                    dest_file = os.path.join(self.data_path, reg + '_' +
                                             str(self.spatial_resolution) +
                                             '_' + self.temporal_resolution +
                                             '.nc')

                wfile = os.path.join(self.weights_path, reg + '_weights_'
                                     + str(ipe) + '.nc')

                if not os.path.isfile(dest_file):
                    grid = grids.ShapeGrid(reg, self.spatial_resolution)
                    save_grid(dest_file, grid)

                with Dataset(dest_file, 'r+', format='NETCDF4') as cdifile:

                    if 'time' not in cdifile.dimensions.keys():
                        dt = get_dtindex(self.temporal_resolution,
                                         self.start_date)
                        cdifile.createDimension("time", None)

                        times = cdifile.createVariable('time', 'uint16',
                                                       ('time',))

                        times.units = 'days since ' + str(self.start_date)
                        times.calendar = 'standard'
                        times[:] = date2num(dt.tolist(), units=times.units,
                                            calendar=times.calendar)

                    else:
                        times = cdifile.variables['time']

                    if key not in cdifile.variables.keys():
                        dim = ('time', 'lat', 'lon')
                        cdi = cdifile.createVariable(key, 'f8',
                                                     dim, fill_value=-99)
                    else:
                        cdi = cdifile.variables[key]

                    for k, gp in enumerate(gps):

                        if k % 100 == 0:
                            print '.',

                        position = np.where(cdifile.variables['gpi'][:] == gp)
                        lat_pos = position[0][0]
                        lon_pos = position[1][0]

                        weights = {}

                        parnum = (len(self.sources.keys()) -
                                  len(self.staticsources))

                        if exclude is not None:
                            parnum = parnum - 1

                        dat = np.zeros((parnum, cdi.shape[0]), dtype=np.float)

                        # dat = np.zeros((len(self.sources.keys()), cdi.shape[0]),
                        #               dtype=np.float)
                        dat[dat == 0] = self.nan_value
                        dat = np.ma.masked_values(dat, self.nan_value)

                        # extract data from DI files and calc weights
                        i = 0

                        for param in self.sources.keys():
                            if param in self.staticsources:
                                continue
                            if param == exclude:
                                continue

                            difile = os.path.join(self.di_path,
                                                  reg + '_' + param
                                                  + '_DI_' + str(ipe) + '.nc')

                            with Dataset(difile, 'r', format='NETCDF4') as nc:
                                for var in nc.variables.keys():
                                    if param in var:
                                        for j in range(0,
                                                       nc.variables[var].shape[0]):
                                            dat[i, j] = (nc.variables[var]
                                                         [j, lat_pos, lon_pos])

                            with Dataset(wfile, 'r', format='NETCDF4') as nc:
                                for var in nc.variables.keys():
                                    if param in var:
                                        weights[param] = (nc.variables[var]
                                                          [lat_pos, lon_pos])
                            i += 1

                        dat = np.ma.masked_where(dat == self.nan_value, dat)
                        dat = np.nan_to_num(dat)
                        dat = np.ma.masked_where(dat == 0., dat)

                        avg = np.ma.average(dat, axis=0,
                                            weights=weights.values())

                        cdi[:, lat_pos, lon_pos] = avg

                    print 'Done!'

        print 'Done!'


