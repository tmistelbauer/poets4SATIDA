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

# Author: Thomas Mistelbauer
# Creation date: 2014-07-29

"""
This module includes the poets base class `Poet`.
"""

from datetime import datetime
from netCDF4 import Dataset
from poets.grid.grids import ShapeGrid
from poets.io.source_base import BasicSource
import numpy as np
import os
import pandas as pd
import poets.web.app as app

valid_temp_res = ['dekad', 'dekade', 'month', 'monthly', 'week', 'weekly',
                  'day', 'daily']


class Poet(object):
    """POETS base class.

    Provides methods to download and resample data using parameters as defined
    in this class. Resampled outputfiles will be saved as NetCDF4 files.

    Parameters
    ----------
    rootpath : str
        path to the directory where data should be stored
    regions : list of str, str, optional
        Identifier of the region in the shapefile. If the default shapefile is
        used, this would be the FIPS country code. Defaults to global.
    region_names : list of str, optional
        Full name of the regions. If set, must have same size as regions
        parameter. Defaults to regions parameter.
    sub_regions : list of str, optional
        Areas or administrative units within the region.
    spatial_resolution : float, optional
        spatial resolution in degree, defaults to 0.25
    temporal_resolution : str, optional
        temporal resolution of the data, possible values: day, week, month,
        dekad; defaults to dekad.
    start_date : datetime.datetime, optional
        first date of the dataset, defaults to 2000-01-01
    nan_value : int
        NaN value to use, defaults to -99
    shapefile : str, optional
        Path to shape file, uses "world country admin boundary shapefile" by
        default. Custom shapefile must use WGS84 as reference system.
    delete_rawdata : bool, optional
        Original files will be deleted from rawdata_path if set True. Defaults
        to False.

    Attributes
    ----------
    rootpath : str
        path to the directory where data should be stored
    regions : list of str
        Identifier of the region in the shapefile.
    region_names : list of str
        Full name of the regions.
    sub_regions : list of str
        Areas or administrative units within the region.
    spatial_resolution : float
        Spatial resolution in degree.
    temporal_resolution : str
        Temporal resolution of the data.
    data_path : str
        Path where resampled NetCDF file is stored.
    rawdata_path : str
        Path where original files are stored and downloaded.
    tmp_path : str
        Path where temporary files are stored.
    nan_value : int
        NaN value to use, defaults to -99.
    start_date : datetime.datetime
        First date of the dataset.
    shapefile : str
        Path to shape file.
    sources : dict of poets.io.BasicSource objects
        Sources used by poets given as BasicSource class.
    delete_rawdata : bool
        Original files will be deleted from rawdata_path if True.
    """

    def __init__(self, rootpath, regions=['global'],
                 spatial_resolution=0.25, temporal_resolution='dekad',
                 start_date=datetime(2000, 1, 1), nan_value=-99, url=None,
                 shapefile=None, delete_rawdata=False, region_names=None,
                 sub_regions=None):

        self.rootpath = rootpath
        if isinstance(regions, str):
            self.regions = [regions]
        else:
            self.regions = regions
        if region_names is not None:
            if isinstance(region_names, str):
                self.region_names = [region_names]
            else:
                self.region_names = region_names
        else:
            self.region_names = self.regions

        self.sub_regions = sub_regions
        self.spatial_resolution = spatial_resolution

        if temporal_resolution not in valid_temp_res:
            raise ValueError("Temporal resulution must be one of " +
                             str(valid_temp_res))

        self.temporal_resolution = temporal_resolution
        self.rawdata_path = os.path.join(rootpath, 'RAWDATA')
        self.data_path = os.path.join(rootpath, 'DATA')
        self.tmp_path = os.path.join(rootpath, 'TMP')
        self.nan_value = nan_value
        self.start_date = start_date
        self.shapefile = shapefile
        self.delete_rawdata = delete_rawdata
        self.sources = {}
        self.url = url

        if not os.path.exists(self.rawdata_path):
            os.mkdir(self.rawdata_path)
        if not os.path.exists(self.tmp_path):
            os.mkdir(self.tmp_path)
        if not os.path.exists(self.data_path):
            os.mkdir(self.data_path)

    def add_source(self, name, filename, filedate, temp_res, host, protocol,
                   username=None, password=None, port=22, directory=None,
                   dirstruct=None, begin_date=None, regions=None,
                   variables=None, nan_value=None,
                   valid_range=None, unit=None, ffilter=None, data_range=None,
                   colorbar=None, src_file=None, labels=None, xticks=None):
        """Creates BasicSource class and adds it to `Poet.sources`.

        Parameters
        ----------
        name : str
            Name of the data source.
        filename : str
            Structure/convention of the file name.
        filedate : dict
            Position of date fields in filename, given as tuple.
        temp_res : str
            Temporal resolution of the source.
        host : str
            Link to data host.
        protocol : str
            Protocol for data transfer.
        username : str, optional
            Username for data access.
        password : str, optional
            Password for data access.
        port : int, optional
            Port to data host, defaults to 22.
        directory : str, optional
            Path to data on host.
        dirstruct : list of strings
            Structure of source directory, each list item represents a
            subdirectory.
        regions : list of str, optional
            List of regions where data from source is available. Uses all
            regions as given in Poet.regions attribute if None.
        begin_date : datetime.date, optional
            Date from which on data is available.
        variables : string or list of strings, optional
            Variables used from data source.
        nan_value : int, float, optional
            Nan value of the original data as given by the data provider.
        valid_range : tuple of int of float, optional
            Valid range of data, given as (minimum, maximum).
        data_range : tuple of int of float, optional
            Range of the values as data given in rawdata (minimum, maximum).
            Will be scaled to valid_range.
        ffilter : str, optional
            Pattern that apperas in filename. Can be used to select out not
            needed files if multiple files per date are provided.
        colorbar : str, optional
            Colorbar to use, use one from
            http://matplotlib.org/examples/color/colormaps_reference.html;
            defaults to jet.
        labels : list, optional
            Custom tick-labels for the legend in the web-app; must have same
            dimension as xticks and only works if xticks is set; Defaults
            to None.
        xticks : list of int or float, optional
            Custom tick locations for the legend in the web-app; must have same
            dimension as labels and only works if labels is set; Defaults
            to None.
        unit : str, optional
            Unit of dataset for displaying in legend. Does not have to be set
            if unit is specified in input file metadata. Defaults to None.
        src_file : dict of str, optional
            Path to file that contains source. Uses default NetCDF file if
            None.
        """

        source = BasicSource(name, filename, filedate, temp_res, self.rootpath,
                             host, protocol, username=username,
                             password=password, port=port, ffilter=ffilter,
                             directory=directory, dirstruct=dirstruct,
                             begin_date=begin_date, variables=variables,
                             nan_value=nan_value, valid_range=valid_range,
                             unit=unit, regions=regions,
                             data_range=data_range, colorbar=colorbar,
                             dest_nan_value=self.nan_value,
                             dest_regions=self.regions,
                             dest_sp_res=self.spatial_resolution,
                             dest_temp_res=self.temporal_resolution,
                             dest_start_date=self.start_date,
                             src_file=src_file, labels=labels, xticks=xticks)

        self.sources[name] = source

    def fetch_data(self, begin=None, end=None, delete_rawdata=None,
                   source=None):
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
        source = str, optional
            Name of the source to fetch_data. If None, all sources will be
            fetched.
        """

        if delete_rawdata is None:
            delete_rawdata = self.delete_rawdata

        for src in self.sources.keys():
            if source is not None:
                if source != src:
                    continue
            s = self.sources[src]
            print '[INFO] Download data for source ' + src
            s.download_and_resample(begin=begin, end=end,
                                    shapefile=self.shapefile,
                                    delete_rawdata=delete_rawdata)

        print '[SUCCESS] Download and resampling complete!'

    def download(self, begin=None, end=None, source=None):
        """Starts download of input data for sources as added
        to `Poets.sources`.

        Parameters
        ----------
        begin : datetime, optional
            Start date of data to download, defaults to start date as defined
            in poets class.
        end : datetime, optional
            End date of data to download, defaults to current datetime.
        source = str, optional
            Name of the source download. If None, all sources will be
            downloaded.
        """

        for src in self.sources.keys():
            if source is not None:
                if source != src:
                    continue
            s = self.sources[src]
            print '[INFO] Download data for source ' + src
            s.download(begin=begin, end=end)

        print '[SUCCESS] Download complete!'

    def resample(self, begin=None, end=None, delete_rawdata=None, source=None):
        """Starts download of input data for sources as added
        to `Poets.sources`.

        Parameters
        ----------
        begin : datetime, optional
            Start date of data to download, defaults to start date as defined
            in poets class.
        end : datetime, optional
            End date of data to download, defaults to current datetime.
        source = str, optional
            Name of the source to resample. If None, all sources will be
            resampled.
        """

        for src in self.sources.keys():
            if source is not None:
                if source is not src:
                    continue
            s = self.sources[src]
            print '[INFO] Resampling data for source ' + src
            s.resample(begin=begin, end=end,
                       shapefile=self.shapefile,
                       delete_rawdata=delete_rawdata)

        print '[SUCCESS] Resampling complete!'

    def fill_gaps(self, source=None):
        """
        Detects gaps in data and tries to fill them by downloading and
        resampling the data within these periods.
        source = str, optional
            Name of the source to fill gaps. If None, all sources will be
            fixed.
        """

        for src in self.sources.keys():
            if source is not None:
                if source is not src:
                    continue
            s = self.sources[src]
            print '[INFO] Scanning ' + src + ' for gaps'
            s.fill_gaps()

    def get_gridpoints(self):
        """Returns gridpoints from NetCDF file.

        Parameters
        ----------
        region : str
            Identifier of the region in the NetCDF file.

        Returns
        -------
        gridpoints : dict of pandas.DataFrame
            Dict containing Dataframes with gridpoint index as index,
            longitutes and latitudes as columns for each region.
        """

        gridpoints = {}

        if self.regions == ['global']:
            filename = (self.regions[0] + '_' + str(self.spatial_resolution)
                        + '_' + str(self.temporal_resolution) + '.nc')
            ncfile = os.path.join(self.data_path, filename)

            with Dataset(ncfile, 'r+', format='NETCDF4') as nc:
                gpis = nc.variables['gpi'][:]
                lons = nc.variables['lon'][:]
                lats = nc.variables['lat'][:]
                gpis = gpis.flatten()
                lons, lats = np.meshgrid(lons, lats)
                lons = lons.flatten()
                lats = lats.flatten()

            points = pd.DataFrame(index=gpis)
            points['lon'] = lons
            points['lat'] = lats
            gridpoints['global'] = points
        else:
            for region in self.regions:
                grid = ShapeGrid(region, self.spatial_resolution,
                                 self.shapefile)
                points = grid.get_gridpoints()
                gridpoints[region] = points

        return gridpoints

    def read_image(self, source, date, region=None, variable=None):
        """Gets images from netCDF file for certain date

        Parameters
        ----------
        date : datetime
            Date of the image.
        source : str
            Data source from which image should be read.
        region : str, optional
            Region of interest, set to first defined region if None.
        variable : str, optional
            Variable to display, set to first variable of source if None.

        Returns
        -------
        img : numpy.ndarray
            Image of selected date.
        lon : numpy.array
            Array with longitudes.
        lat : numpy.array
            Array with latitudes.
        metadata : dict
            Dictionary containing metadata of the variable.
        """

        img, lon, lat, metadata = self.sources[source].read_img(date, region,
                                                                variable)

        return img, lon, lat, metadata

    def read_timeseries(self, source, location, region=None, variable=None):
        """
        Gets timeseries from netCDF file for a gridpoint.

        Parameters
        ----------
        source : str
            Data source from which time series should be read.
        location : int or tuple of floats
            Either Grid point index as integer value or Longitude/Latitude
            given as tuple.
        region : str, optional
            Region of interest, set to first defined region if None.
        variable : str, optional
            Variable to display, set to first variable of source if None.

        Returns
        -------
        ts : pd.DataFrame
            Timeseries for the selected data.
        """

        ts = self.sources[source].read_ts(location, region, variable,
                                          shapefile=self.shapefile)

        return ts

    def bulkread_timeseries(self, source, locations, region=None,
                            variable=None, grid=None):
        """
        Gets timeseries from netCDF file for a gridpoint.

        Parameters
        ----------
        source : str
            Data source from which time series should be read.
        location : int or tuple of floats
            Either Grid point index as integer value or Longitude/Latitude
            given as tuple.
        region : str, optional
            Region of interest, set to first defined region if None.
        variable : str, optional
            Variable to display, set to first variable of source if None.

        Returns
        -------
        ts : pd.DataFrame
            Timeseries for the selected data.
        """

        ts = self.sources[source].bulkread_ts(locations, region, variable,
                                              shapefile=self.shapefile,
                                              grid=grid)

        return ts

    def average_timeseries(self, source, region, variable=None):
        """
        Calculates mean of all time series in a region.

        Parameters
        ----------
        source : str
            Data source from which time series should be read.
        region : str, optional
            Region of interest, set to first defined region if None.
        variable : str, optional
            Variable to display, set to first variable of source if None.

        Returns
        -------
        ts : pd.DataFrame
            Timeseries for the selected data.
        """

        grid = ShapeGrid(region, self.spatial_resolution, self.shapefile)
        points = grid.get_gridpoints()

        df = pd.DataFrame()

        locations = []

        if region in self.regions:
            region = region
            locations = list(points.index)
        else:
            for i, sr in enumerate(self.sub_regions):
                if region in sr:
                    idx = i
                    break
            region = self.regions[idx]
            lat = points['lat'].tolist()
            lon = points['lon'].tolist()
            for i in range(0, points.shape[0]):
                locations.append((lon[i], lat[i]))

        if len(points) < 1:
            return 'ERROR: No points available in the selected region.'

        df, gpis = self.bulkread_timeseries(source, locations, region,
                                            variable=variable)

        df_mean = pd.DataFrame()

        for var in df.keys():
            df_mean[var + '_mean'] = df[var].mean(axis=1)

        return df_mean, len(gpis)

    def get_variables(self, region=None):
        """
        Returns all variables available.

        Parameters
        ----------
        region : str, optional
            Region to check for variables.

        Returns
        -------
        variables : list of str
            Sorted list of all variables.
        """

        variables = []

        for src in self.sources.keys():
            if region is None:
                for var in self.sources[src].get_variables():
                    variables.append(var)
            else:
                if self.sources[src].regions is not None:
                    if region in self.sources[src].regions:
                        for var in self.sources[src].get_variables():
                            variables.append(var)
                else:
                    for var in self.sources[src].get_variables():
                        variables.append(var)

        variables.sort()

        return variables

    def start_app(self, host='127.0.0.1', port=None, r_host=None, r_port=None,
                  url=None, debug=False):
        """Starts web interface.

        Parameters
        ----------
        host : str, optional
            Host that is used by the app, defaults to 127.0.0.1.
        port : int, optional
            Port where app runs on, defaults to None.
        r_host : str, optional
            IP of router that is between host and internet.
        r_port : int, optional
            Port of router that is between host and internet.
        url : str, optional
            URL where poets web interface is hosted. If not set, poets will run
            on the IP of the hosting machine.
        debug : bool, optional
            Starts app in debug mode if set True, defaults to False.
        """

        app.start(self, host, port, r_host, r_port, url, debug)
