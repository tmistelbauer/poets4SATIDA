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
# Creation date: 2014-08-04

"""
Description of module.
"""

import pandas as pd
import numpy as np
from poets.timedate.dekad import get_dekad_period


def calc_CDI(data, refparam=None, lags=[0, 10]):
    """
    Calculates a weighted average over all columns of a pandas DataFrame.

    Parameters
    ----------
    data : pandas.DataFrame
        Pandas DataFrame containing data to be averaged.
    refparam : str, optional
        Reference parameter. If not set, parameters will be weighted
        equally.
    lags : list of int, optional
        Time periods to shift parameter against refparam, defaults to [0, 10].

    Returns
    -------
    df : pandas DataFrame
        Return the average of data
    """

    cols = data.keys()
    dat = np.array(data[cols])
    dat = np.ma.masked_invalid(dat)

    weights = calc_weights(data, refparam, lags)

    if refparam is None:
        avg = np.ma.average(dat, axis=1)
    else:
        avg = np.ma.average(dat, axis=1, weights=weights)
    df = pd.DataFrame(avg, columns=['CDI'], index=data.index)

    return df


def calc_weights(data, refparam, lags=[0, 10], exclude=None):
    """
    Calculates the weights of parameters for weighted averaging. Weights
    are calculated using correlation and time shift of each parameter
    against the reference parameter. Parameters must be direct proportional
    to reference parameter!

    Parameters
    ----------
    data : pandas.DataFrame
        DataFrame containing data in columns.
    refparam : str
        Reference parameter.
    lags : list of int, optional
        Time periods to shift parameter against refparam,
        defaults to [0, 10].
    exclude : string, optional
        Variable which should not be used for calculation of the weights.

    Returns
    -------
    sorted_weights : list of int
        Weights associated with the parameters in data.
    """

    params = data.keys()

    maxlag = {}
    maxcorr = {}
    weights = {}
    sorted_weights = []

    correlations = calc_correlation(data, refparam, lags, exclude)

    for param in params:
        if exclude is not None and exclude in param:
            continue
        maxlag[param] = correlations[param]['lag']
        maxcorr[param] = correlations[param]['corr']

    print maxlag.keys()
    print maxcorr.keys()

    for key in maxlag.keys():
        weights[key] = (float(maxlag[key])) / sum(maxlag.values()) * 100

    for key in maxcorr.keys():
        weights[key] = ((weights[key] +
                         (float(maxcorr[key]) / sum(maxcorr.values())) * 100)
                        / 2)

    for param in params:
        if exclude is not None and exclude in param:
            continue
        sorted_weights.append(weights[param])

    return sorted_weights


def calc_correlation(data, refparam, lags=[0, 10], exclude=None):
    """
    Calculates the correlations between parameters and a reference
    parameter given as columns in a DataFrame.

    Parameters
    ----------
    data : pandas.DataFrame
        DataFrame containing data in columns.
    refparam : str
        Reference parameter.
    lags : list of int, optional
        Time periods to shift parameter against refparam,
        defaults to [0, 10].
    exclude : string, optional
        Variable which should not be used for calculation of the correlation.

    Returns
    -------
    correlation : dict
        Dictionary containing correlations and max time lags.
    """

    correlation = {}

    for param in data.keys():
        if exclude is not None and exclude in param:
            continue
        correlation[param] = {'corr': None, 'lag': None}
        for i in range(lags[0], lags[1]):
            i += abs(lags[0]) + 1
            corr = data[param].corr(data[refparam].shift(periods=i),
                                    method='pearson')
            if correlation[param]['corr'] is None:
                correlation[param]['corr'] = abs(corr)
                correlation[param]['lag'] = i
            if abs(corr) > abs(correlation[param]['corr']):
                correlation[param]['corr'] = abs(corr)
                correlation[param]['lag'] = i
            if abs(corr) == abs(correlation[param]['corr']):
                if abs(i) < abs(correlation[param]['lag']):
                    correlation[param]['corr'] = abs(corr)
                    correlation[param]['lag'] = i

    return correlation


def calc_DI(data, inverse=False, interest_period=[6, 12, 24], scaled=False,
            scale_zero=False, modf_all=False):
    """
    Calculates a Drought Index based on an algorithm developed by
    FAO SWALIM.

    Parameters
    ----------
    data : pandas.DataFrame
        Input data as Pandas DataFrame, must come with column names.
    inverse : bool
        Inverts the input time series; set True if time series is indirect
        proportional to the expected  output, e.g. Temperature with output
        Temperature Drought Index.
    interest_period : list of int, optional
        interest periods used to calculate drought index,
        defaults to [6, 12, 24]
    scaled : boolean, optional
        If True values will be scaled between 0 and 1.
    scale_zero : boolean, optional
        If True values will be shifted around zero, defaults to False.
    modf_all : boolean, optional
        If True values will be modified, independent of their min.
    """

    ts_date = data.index
    variables = data.keys()
    data['period'] = get_dekad_period(ts_date)

    for var in variables:

        if inverse is True:
            data[var] = ((data[var].max() + 1) - data[var])

        if modf_all is True:
            data['modf'] = data[var] + 1
            del data[var]
        elif data[var].min() == 0:
            data['modf'] = data[var] + 1
            del data[var]
        else:
            data['modf'] = data[var]
            del data[var]

        data['modf_avg'] = (data.groupby('period').modf
                            .transform(lambda x: x.mean()))

        # Excess
        # Dekads below long term average. If the statement is true the
        # program return 1
        data['exc'] = np.choose((data['modf_avg'] / data['modf']) >= 1,
                                [0, 1])

        # Run length
        # Maximum number of successive dekads below long term average
        for ip in interest_period:
            data['rlen'] = pd.rolling_apply(data['exc'], ip,
                                            (lambda x:
                                             len(max((''.join(str(j)
                                                              for j in map(int,
                                                                           x)))
                                                     .split('0')))),
                                            ip)

            # get modified run length
            max_rlen = data['rlen'].max()
            data['rlen'] = (max_rlen + 1) - data['rlen']

            # average run lenghts
            rlen_avg = (data.groupby('period').modf
                        .transform(lambda x: x.mean()))
            data['form'] = data['rlen'] / rlen_avg

            # sumip matrix
            # calculates sum of the values for each interest period
            data['sumip'] = pd.rolling_apply(data['modf'], ip,
                                             lambda x: sum(x), ip)

            # average values for each interest period over all years
            sumip_avg = (data.groupby('period')['sumip']
                         .transform(lambda x: x.mean()))
            data['nrl'] = data['sumip'] / sumip_avg

            # calculating PDI/TDI
            data['val'] = data['nrl'] * np.sqrt(data['form'])

            # scaled index
            dkey = var + '_DI_' + str(ip)
            if scaled:
                data[dkey] = ((data['val'] - data['val'].min()) /
                              (data['val'].max() - data['val'].min()))
            else:
                data[dkey] = data['val']

            if scale_zero:
                data[dkey] = data[dkey] - data[dkey].mean()

            del (data['val'], data['nrl'], data['sumip'], data['rlen'],
                 data['form'])

        # deletes not further relevant columns
        del data['modf'], data['modf_avg'], data['exc']

    del data['period']

    return data


if __name__ == "__main__":
    pass
