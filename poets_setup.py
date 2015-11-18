import os
import sys
import getopt
from datetime import datetime
from poets.poet import Poet
import ConfigParser


def curpath():
    pth, _ = os.path.split(os.path.abspath(__file__))
    return pth

if __name__ == "__main__":

    # os.environ['LD_PRELOAD'] = '/usr/local/lib/libgdal.so.1'

    cfg_parser = ConfigParser.RawConfigParser()
    cfg_path = os.path.join(curpath(), 'config.ini')
    cfg_parser.read(cfg_path)
    cfg_rootpath = cfg_parser.get('settings', 'rootpath')
    cfg_shapefile = cfg_parser.get('settings', 'shapefile')
    cfg_host = cfg_parser.get('settings', 'host')
    cfg_port = cfg_parser.get('settings', 'port')
    cfg_url = cfg_parser.get('settings', 'url')
    cfg_rhost = cfg_parser.get('settings', 'r_host')
    cfg_rport = cfg_parser.get('settings', 'r_port')

    if cfg_url == 'None':
        cfg_url = None

    if cfg_port == 'None':
        cfg_port = None
    else:
        cfg_port = int(cfg_port)

    if cfg_rport == 'None':
        cfg_rport = None
    else:
        cfg_rport = int(cfg_rport)

    if cfg_host == 'None':
        cfg_host = None

    if cfg_rhost == 'None':
        cfg_rhost = None

    # POETS SETTINGS
    spatial_resolution = 0.25
    temporal_resolution = 'dekad'
    rootpath = cfg_rootpath
    shapefile = cfg_shapefile
    nan_value = -99
    start_date = datetime(1992, 1, 1)
    regions = ['CE', 'IN']
    region_names = ['Sri Lanka', 'India']
    sub_regions = [['Ampara', 'Anuradhapura', 'Badulla', 'Batticaloa', 'Colombo',
                    'Vavuniya', 'Galle', 'Gampaha', 'Hambantota', 'Jaffna',
                    'Kalutara', 'Kandy', 'Kegalle', 'Kilinochchi', 'Kurunegala',
                    'Mannar', 'Matale', 'Matara', 'Moneragala', 'Mullaitivu',
                    'Nuwara Eliya', 'Polonnaruwa', 'Puttalam', 'Ratnapura',
                    'Trincomalee', 'Vavuniya'],
                   ['Andhra Pradesh', 'Arunachal Pradesh', 'Assam', 'Bihar',
                    'Chhattisgarh', 'Goa', 'Gujarat', 'Haryana',
                    'Himachal Pradesh', 'Jammu and Kashmir', 'Jharkhand',
                    'Karnataka', 'Kerala', 'Madhya Pradesh', 'Maharashtra',
                    'Manipur', 'Meghalaya', 'Mizoram', 'Nagaland', 'Odisha',
                    'Punjab', 'Rajasthan', 'Sikkim', 'Tamil', 'Telangana',
                    'Tripura', 'Uttar Pradesh', 'Uttarakhand', 'West Bengal']]
    delete_rawdata = True

    p = Poet(rootpath, regions, spatial_resolution, temporal_resolution,
             start_date, nan_value, region_names=region_names,
             shapefile=shapefile, sub_regions=sub_regions,
             delete_rawdata=delete_rawdata)

    # SOURCE SETTINGS ECV
    name = 'satellite-derived surface soil moisture'
    temp_res = 'daily'
    host = "ftp.ipf.tuwien.ac.at/"
    directory = "_down/CCI_ASCAT_AMSR2_merged"
    port = 22
    protocol = 'SFTP'
    username = 'crts'
    password = 'Ibonaguma343'
    filename = ("merged_NRT_{YYYY}{MM}{TT}.nc")
    filedate = {'YYYY': (11, 15), 'MM': (15, 17), 'DD': (17, 19)}
    begin_date = datetime(1978, 11, 01)
    variables = ['sm']
    valid_range = (0, 0.6)
    colorbar = 'jet_r'
    nan_value = -9999

    p.add_source(name, filename, filedate, temp_res, host, protocol,
                 username=username, password=password, port=port,
                 directory=directory, begin_date=begin_date,
                 colorbar=colorbar, variables=variables,
                 valid_range=valid_range, nan_value=nan_value)

    # SOURCE SETTINGS ECV SWI
    name = 'satellite-derived soil water index'
    temp_res = 'daily'
    host = '/home/eodc/Datapool/SWI_092_dailyimages/'
    dirstruct = ['YYYY']
    protocol = 'LOCAL'
    filename = 'ESACCI-SOILMOISTURE-L3S-SWIV-COMBINED-{YYYY}{MM}{DD}000000-fv02.2.nc'
    filedate = {'YYYY': (38, 42), 'MM': (42, 44), 'DD': (44, 46)}
    begin_date = datetime(1978, 11, 01)
    variables = ['swi']
    valid_range = (0, 0.6)
    colorbar = 'jet_r'
    nan_value = -9999
    src_file = {}
    for reg in regions:
        src_file[reg] = os.path.join(p.data_path, 'DATA', reg + '_' +
                                     str(p.spatial_resolution) + '_' +
                                     p.temporal_resolution + '_SWI.nc')

    p.add_source(name, filename, filedate, temp_res, host, protocol,
                 dirstruct=dirstruct, begin_date=begin_date,
                 variables=variables, valid_range=valid_range,
                 colorbar=colorbar, nan_value=nan_value, src_file=src_file)

    argv = sys.argv[1:]

    try:
        opts, args = getopt.getopt(argv, 'ho:', ['help', 'action='])
    except:
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            print '-o [start_app, fetch_data]'
        elif opt == '-o':
            if arg == 'start_app':
                p.start_app(host=cfg_host, port=cfg_port, r_host=cfg_rhost,
                            r_port=cfg_rport, url=cfg_url)
            elif arg == 'fetch_data':
                p.fetch_data()
            elif arg == 'fill_gaps':
                p.fill_gaps()
            else:
                print 'unknown argument'
        else:
            print 'please give an argument'

    if len(opts) == 0:
        p.start_app()

