#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from math import radians, cos, sin, asin, sqrt
import ftplib
import pycountry_convert as pc
from geopy.geocoders import Nominatim
import urllib
import xmltodict
from statsmodels.robust.robust_linear_model import RLM
from datetime import datetime as dt
import pandas as pd
import os
import numpy as np


def get_spotgins_dir():
    spotgins_dir = os.getenv('SPOTGINS_DIR')
    if not spotgins_dir:
        print("WARN: environment variable $SPOTGINS_DIR is not defined")
        print("      It contains the path to the spotgins root folder")
    return spotgins_dir


def haversine(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points on the earth (specified in decimal degrees)
    :return: distance in km
    """

    # convert decimal degrees to radians
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
    # haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a))
    km = 6367 * c
    return km

def read_ngl(acro4,start=None,end=None):

    try:
        data = pd.read_csv('http://geodesy.unr.edu/gps_timeseries/tenv3/IGS14/' + acro4 + '.tenv3',delim_whitespace=True,skiprows=[0],usecols=[1,2,8,10,12,14,15,16],names=['datestr','decy','e','n','u','se','sn','su'],index_col='date_dt',parse_dates={'date_dt':['datestr']},date_format='%y%b%d')
    except:
        return None

    if start is None:
        start = data.index[0]
    if end is None:
        end = data.index[-1]

    return data[start:end]

def read_uga(filename):
    data = pd.read_csv(filename,delim_whitespace=True,skiprows=50,usecols=[0,15,16,17],header=None,names=['datestr','dN','dE','dU'],parse_dates={'date_dt':['datestr']},date_format='%Y%m%d',index_col='date_dt')

    return data

def get_codomes(local_file = None):
    """
    Get the list of official domes numbers from:
        - ftp://igs-rf.ign.fr/pub/DOMES/codomes_gps_coord.snx if local_file is None
        - local_file if local_file is not None
    :return: dictionnary (keys: 4char acronyms) of dictionarry (keys: long, lat, domes)
    """

    if local_file is None:
        #getting the file from ftp
        try:
            ftp = ftplib.FTP('igs-rf.ign.fr')
            ftp.login('anonymous', 'spotgins')

            ftp.cwd('/pub/DOMES')

            try:
                ftp.retrbinary('RETR codomes_gps_coord.snx', open('codomes_gps_coord.snx', 'wb').write)
                print ('codomes_gps_coord.snx download ok')
                ftp.quit()
            except:
                print ('FATAL: Unable to download codomes_gps_coord.snx')
                ftp.quit()
                return

        except:
            print ('FATAL: Unable to connect to igs-rf.ign.fr')
            return

        data = open('codomes_gps_coord.snx')

    else:
        data = open(local_file)

    #initialising output dictionnary
    codomes = {}

    for line in data:

        #4char acronym
        acro = line[0:4]
        #domes number
        domes = line[8:17]
        #extracting coordinates
        #detecting the index of the first coordinate (i.e: first float in the rest of the line)
        f = line[18:].split()
        found = False
        for (ind, el) in enumerate(f):
            if not found:
                if '.' in el:
                    try:
                        float(el)
                        found = True
                    except:
                        continue
            else:
                break
        long = f[ind - 1]
        lat = f[ind]

        if acro in codomes.keys():
            codomes[acro].append({'domes':domes,'long':float(long),'lat':float(lat)})
        else:
            #incrementation fo the output dictionnary
            codomes[acro] = [{'domes':domes,'long':float(long),'lat':float(lat)}]

    data.close()

    return codomes

def get_country_code(long,lat):
    """
    Get the 3char ISO code of the country from longitude and latitude
    :return: 3char ISO country code
    """

    geolocator = Nominatim(user_agent="spotgins")

    try:
        #2char code
        alpha2 = geolocator.reverse(str(lat) + ', ' + str(long),timeout=10).raw['address'] ['country_code'].upper()
        #country name
        country_name = pc.country_alpha2_to_country_name(alpha2)
        #3char code
        alpha3 = pc.country_name_to_country_alpha3(country_name)
    except:
        alpha3 = None

    return alpha3

def get_existing_longnames(acro):

    existing = {}

    m3glist = pd.read_html('https://gnss-metadata.eu/data/station/xml/')[0]
    m3glist = m3glist.drop(0)
    m3glist = m3glist[m3glist['File Name ↓'].str[0:4] == acro]

    for xmlfile in m3glist['File Name ↓'].values:
        acro9 = xmlfile[0:9]
        (mlat,mlon,mh) = get_M3G_pos(acro9)
        existing[acro9]=[mlon,mlat]

    igs = pd.read_csv('https://network.igs.org/api/public/stations.csv?draw=1&length=514&current=true&ordering=name&fields%5B%5D=name&fields%5B%5D=country&fields%5B%5D=latitude&fields%5B%5D=longitude&start=0',sep=',',header=None,skiprows=[0],usecols=[0,2,3],names=['acro9','lat','long'])

    igs = igs[igs['acro9'].str[0:4] == acro]

    for (acro9,ilat,ilong) in igs.values:
        if acro9 not in existing.keys():
            existing[acro9] = [ilong,ilat]

    return existing

def get_country_code_longname(acro,long,lat):
    """
    Get the 3char ISO code of the country from longitude and latitude for editing station longname
    :return: customized 3char ISO country code
    """

    #already exist?
    existing = get_existing_longnames(acro)

    # m3glist = pd.read_html('https://gnss-metadata.eu/data/station/xml/')[0]
    # igs = pd.read_csv('https://network.igs.org/api/public/stations.csv?draw=1&length=514&current=true&ordering=name&fields%5B%5D=name&fields%5B%5D=country&fields%5B%5D=latitude&fields%5B%5D=longitude&start=0', sep=',')
    #
    # m3glist = m3glist[m3glist['File Name ↓'].str[0:4] == acro]
    # igs = igs[igs['Site Name'].str[0:4] == acro]
    found_sta = False


    for acro9 in existing.keys():
        (lon_tmp,lat_tmp) = existing[acro9]

        if haversine(long, lat, lon_tmp, lat_tmp) < 0.5:
            found_sta=True
            country_code = acro9[6:9]
            break

    if not found_sta:
        geolocator = Nominatim(user_agent="spotgins")

        try:
            loc = geolocator.reverse(str(lat) + ', ' + str(long),timeout=10)
            # 2char code
            alpha2 = loc.raw['address'] ['country_code'].upper()
            # country name
            country_name = pc.country_alpha2_to_country_name(alpha2)
            # 3char code
            alpha3 = pc.country_name_to_country_alpha3(country_name)

            #particular cases following the IGS existing names
            if 'state' in loc.raw['address']:
                if loc.raw['address']['state'] == 'Polynésie Française':
                    alpha3 = 'PYF'
                if loc.raw['address']['state'] == 'Guadeloupe':
                    alpha3 = 'GLP'
                if loc.raw['address']['state'] == 'Martinique':
                    alpha3 = 'MTQ'
                if loc.raw['address']['state'] == 'Saint-Pierre-et-Miquelon':
                    alpha3 = 'SPM'
                if loc.raw['address']['state'] == 'La Réunion':
                    alpha3 = 'REU'
                if loc.raw['address']['state'] == 'Guyane':
                    alpha3 = 'GUF'
                if loc.raw['address']['state'] == 'Mayotte':
                    alpha3 = 'MYT'
                if loc.raw['address']['state'] == 'United States Virgin Islands':
                    alpha3 = 'VIR'
                if loc.raw['address']['state'] == 'Guam':
                    alpha3 = 'GUM'
                if loc.raw['address']['state'] == 'Terres australes et antarctiques françaises':
                    alpha3 = 'ATF'
                if loc.raw['address']['state'] == 'Saint Helena':
                    alpha3 = 'GBR'
                if loc.raw['address']['state'] == 'Wallis-et-Futuna':
                    alpha3 = 'WLF'
                if 'Hong Kong' in loc.raw['address']['state']:
                    alpha3 = 'HKG'
            if 'archipelago' in loc.raw['address'].keys():
                if loc.raw['address']['archipelago'] == 'Nouvelle-Calédonie':
                    alpha3 = 'NCL'

            if 'country' in loc.raw['address'].keys():
                if loc.raw['address']['country'] == 'British Indian Ocean Territory':
                    alpha3 = 'GBR'
                if loc.raw['address']['country'] == 'Palestinian Territory':
                    alpha3 = 'ISR'

        except:
            alpha3 = 'XXX'

        #Antartica
        if lat < -60:
            alpha3 = 'ATA'

        country_code = alpha3

    return country_code

def read_rcvant():

    data = open('D:/Projets/SPOTGINS/rcvant.dat')
    rcvant_data= {'receivers':[],'antennas':[]}
    read_rcv = False
    read_ant = False
    for line in data:

        # print (line)
        # print (read_rcv)
        # print (read_ant)
        if line[0:1] == ' ' and line.strip() != '':
            if line[0:10] == ' RECEIVERS':
                read_rcv = True
                continue
            if line[0:9] == ' ANTENNAS':
                read_ant = True
                continue
            if line[0:4] == 'END':
                read_rcv = False
                read_ant = False
                continue

            if read_rcv:
                rec = line[15:36].strip()
                if rec != '':
                    rcvant_data['receivers'].append(line[15:36].strip())
            if read_ant:
                ant = line[15:36].strip()
                if ant != '':
                    rcvant_data['antennas'].append(line[15:36].strip())

    return rcvant_data

def extract_from_stafile(sta,date_dt,acro9):
    """
    Extract receiver, antenna and eccentricities from a station file at a given date

    Also implemented as a method `extract` for StationFile objects
    """

    found = False
    for nch in sta.data[acro9].keys():
        staline = sta.data[acro9][nch]

        if staline['start'] <= date_dt < staline['end']:
            rec = staline['rec']
            ant = staline['ant']
            ant[0] = ant[0].upper()
            ant[1] = ant[1].upper().replace('----','NONE').replace('    ','NONE').replace('UNKN','NONE')
            ecc = {'u':staline['ecc_u'],'e':staline['ecc_e'],'n':staline['ecc_n']}
            found = True


            return (rec.upper(),ant,ecc)
    if not found:
        #print ('Unable to extract information at ' + str(date_dt) + ' for ' + acro9)
        return (None,None,None)

def list_cmtfile():
    #ndk files at https://www.ldeo.columbia.edu/~gcmt/projects/CMT/catalog/NEW_MONTHLY/
    eqfiles = []
    for eqfile in os.listdir('D:/Donnees/CMT'):
        if eqfile.endswith('.ndk'):
            eqfiles.append('D:/Donnees/CMT/' + eqfile)
    return eqfiles

def get_continent_code(country_code):
    """
    Get the 2char continent code from 3char ISO country code
    :return: 2char continent code
    """
    #2char country code
    alpha2 = pc.country_alpha3_to_country_alpha2(country_code)
    #continent code
    continent_code = pc.country_alpha2_to_continent_code(alpha2)

    return continent_code

def get_domes_region_code(country_code,continent_code,lat):
    """
    Get the region code (1=>9) for editing the domes number followinf this convention:
        01 = Europe (West of Russia and without Turkey).
        02 = Asie (with Turkey, all East of the Red Sea, include Archipels in South-East Asia)
        03 = Africa
        04 = North America (USA, Mexico, Canada and Alaska)
        05 = Oceania (Australia, New-Zealand and near)
        06 = Antarctica
        07 = South America (starts in the south of Mexico)
        08 = Isolated territories / Islands (all that are not included in any other regions)
        09 = Nothern Lands (Svalbard, Greenland, Iceland)
    :return: region code
    """

    if continent_code is not None:
        codes = {'EU':1,'AS':2,'AF':3,'NA':4,'OC':5,'SA':7}
        region_code = codes[continent_code]
    else:
        region_code = 8

    #Greenland, Iceland
    if country_code in ['ISL','GRL']:
        region_code = 9

    #Svalbard
    if country_code == 'NOR' and lat > 75:
        region_code = 9

    #Antartica
    if lat < -60:
        region_code = 6

    return str(region_code)

def is_M3G(long,lat,m3glist):


    found = False
    for xmlfile in m3glist['File Name ↓'].values:
        acro9 = xmlfile[0:9]
        (mlat,mlon,mh) = get_M3G_pos(acro9)

        if haversine(long,lat,mlon,mlat) < 0.5:
            #same station
            found = True
            return acro9
    if not found:
        return None

def get_M3G_pos(acro9):


    file = urllib.request.urlopen("https://gnss-metadata.eu/data/station/xml/" + acro9 + ".xml")
    data = file.read()
    file.close()
    data = xmltodict.parse(data)

    (lat, long, h) = data['geo:GeodesyML']['geo:siteLog']['geo:siteLocation']['geo:SiteLocation']['geo:approximatePositionITRF'][
        'geo:geodeticPosition']['gml:Point']['gml:pos']['#text'].split()

    return (float(lat),float(long),float(h))

def build_domes(long,lat,existing_domes):
    '''
    Constructs a domes numbers from longitude, latitude, and a list of existing domes, following these rules:
    1st digit: 0 (=> homemade domes)
    2nd digit: region code (from get_domes_region_code function)
    3rd:5th digits: number between 001 and 999. Follow the chronological order of homemade domes in the corresponding region
    6th digit: X
    7th:9th digits: occupation number. For now: '001' hard-written.
    :return: domes number
    '''
    try:
        #country 3char code
        cntry_code = get_country_code(long,lat)
        #continent 2char code
        cont_code = get_continent_code(cntry_code)

    except:
        cntry_code = None
        cont_code = None

    #region code
    domes_region_code =  get_domes_region_code(cntry_code,cont_code,lat)

    #existing homemade domes in that region
    order_in_region = 1
    for dome in existing_domes:
        if dome.startswith('0' + domes_region_code):
            n = float(dome[2:5])

            if n >= order_in_region:
                order_in_region = n+1

    final_domes = '0' + domes_region_code + str(int(order_in_region)).zfill(3) + 'X001'

    return final_domes

def build_id(long,lat,existing_ids):
    '''
    Constructs an id from longitude, latitude, and a list of existing ids, following these rules:
    1st digit: region code (from get_domes_region_code function)
    2nd:5th digits: number between 0001 and 9999. Follow the chronological order of homemade ids in the corresponding region
    :return: id
    '''
    try:
        #country 3char code
        cntry_code = get_country_code(long,lat)
        #continent 2char code
        cont_code = get_continent_code(cntry_code)

    except:
        cntry_code = None
        cont_code = None

    #region code
    region_code =  get_domes_region_code(cntry_code,cont_code,lat)

    #existing homemade ids in that region
    order_in_region = 1
    for id in existing_ids:

        id = str(id)
        if id.startswith(region_code):
            n = float(id[1:5])

            if n >= order_in_region:
                order_in_region = n+1

    final_id = region_code + str(int(order_in_region)).zfill(4)

    return int(float(final_id))

def build_longname(acro,long,lat):
    """
    Constructs the 9char longname of GNSS station
    For the country code, it tries to follow the conventions of the existing IGS stations. See get_country_code_longname function.
    :return: long name
    """

    cntry_code = get_country_code_longname(acro,long,lat)

    print(cntry_code)
    if acro == 'GUUG':
        cntry_code = 'USA'
    long_name = acro + '00' + cntry_code

    return long_name

#-------------------------------------------------------------------------------
# Routine : geo2cart
# Purpose : Transform geographical coordinates into cartesian coordinates
# Author  : P. Rebischung
# Created : 25-May-2011
#
# Changes :
#
# Input   : - phi : Latitudes  (rad)
#           - lam : Longitudes (rad)
#           - h   : Ellipsoidal heights (m)
# Output  : - X   : [X, Y, Z] (cartesian coordinates in m)
#-------------------------------------------------------------------------------
def geo2cart(phi, lam, h):
    '''
    Converts geographic coordinates (in radians) to cartesian coordinates (in meters)
    From P. Rebischung, 25-May-2011
    :return: [X, Y, Z] (cartesian coordinates in m)
    '''
    # GRS80 parameters
    ae = 6378137.
    fe = 0.00335281068118
    ee = sqrt(2 * fe - fe ** 2)

    N = ae / np.sqrt(1 - (ee * np.sin(phi)) ** 2)

    if (type(phi) == type(float())):
        X = np.zeros(3, )
    else:
        X = np.zeros(phi.shape + (3,))

    X[..., 0] = (N + h) * np.cos(phi) * np.cos(lam)
    X[..., 1] = (N + h) * np.cos(phi) * np.sin(lam)
    X[..., 2] = (N * (1 - ee ** 2) + h) * np.sin(phi)

    return X

def linear_trend(dates, h):
    '''

    Calcul d'une tendance lineaire robuste
    Attention: bien verifier que la librairie statsmodels est a jour.
    Derniere version: 0.12.2

    Entree:
        - dates: list de dates en annees decimales, sans nan
        - h: liste des observations, sans nan
    Sortie:
        - trend: tendance lineaire
        - unc: incertitude sur la tendance lineaire
        - fit: liste contenant le modele lineaire estime
    '''

    # mise en forme des vecteurs
    try:
        float(dates[0])
    except:
        #dates as datetime
        dates = dt2decyear(dates)

    x1 = np.linspace(1, 1, len(h))
    X = np.column_stack((dates, x1))

    # regression lineaire robuste
    reg = RLM(h, X).fit()

    # tendance
    trend = reg.params[0]
    # incertitude sur la tendance
    unc = reg.bse[0]
    # valeurs estimées
    fit = reg.fittedvalues

    return (trend, unc, fit)


def dt2decyear(dates_dt):

    dates_dy = []
    for date_dt in dates_dt:
        year_part = date_dt - dt(year=date_dt.year, month=1, day=1)
        year_length = (
            dt(year=date_dt.year + 1, month=1, day=1)
            - dt(year=date_dt.year, month=1, day=1)
        )

        dates_dy.append(date_dt.year + year_part / year_length)
    return dates_dy

def get_spotgins_files(master_path=None, stafile_path=None):
    """
    Get the path of the masterfile and stationfile from the given path or from the $SPOTGINS_DIR
    """
    if not master_path or not stafile_path:
        spotgins_dir = get_spotgins_dir()
        if not spotgins_dir:
            print("ERROR: no $SPOTGINS_DIR in env & no explicit masterfile/stationfile path given. Abort.")
            return None
        if not stafile_path:
            stafile_path = os.path.join(spotgins_dir, 'metadata/stations/station_file.dat')
        if not master_path:
            master_path = os.path.join(spotgins_dir, 'metadata/stations/station_master_file.dat')

    if not os.path.isfile(master_path) or not os.path.isfile(stafile_path):
        print("ERROR: no masterfile/stationfile exist in given path. Abort.")
        print(master_path)
        print(stafile_path)
        return None, None

    return master_path, stafile_path