# -*- coding: utf-8 -*-
"""
Config file for globals used in the future weather project.

@author: paul
"""


#=========================================================
#            MACAV2METDATA FILE PARAMETERS
#=========================================================
# VARNAME, VARLONGNAME, MODELNAME, SCENNAME, YEAR_START, YEAR_END, RUN_NUM, DOMAIN

VARNAME =('tasmax','tasmin','rhsmax','rhsmin','pr','rsds','uas','vas','huss')
VARLONGNAME=('air_temperature','air_temperature', 'relative_humidity','relative_humidity',\
                'precipitation','surface_downwelling_shortwave_flux_in_air',\
                'eastward_wind','northward_wind','specific_humidity')
MODELNAME=('bcc-csm1-1','bcc-csm1-1-m','BNU-ESM','CanESM2','CCSM4','CNRM-CM5','CSIRO-Mk3-6-0',\
            'GFDL-ESM2G','GFDL-ESM2M','HadGEM2-CC365','HadGEM2-ES365',\
                'inmcm4','IPSL-CM5A-MR','IPSL-CM5A-LR','IPSL-CM5B-LR',\
                'MIROC5','MIROC-ESM','MIROC-ESM-CHEM','MRI-CGCM3','NorESM1-M')
SCENNAME =('historical','rcp45','rcp85') 
YEAR_START=('1950','2006','2006')
YEAR_END =('2005','2099','2099')
RUN_NUM = [1] * 20
RUN_NUM[4] = 6 #setting CCSM4 with run 6 # this is just what MACA used

DOMAIN='CONUS'

#=========================================================
#               Time Constants
#=========================================================
HRS_IN_MONTH = [744,672,744,720,744,720,744,744,720,744,720,744];

HRS_IN_YEAR = 8760;

DAYS_IN_MONTH = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];

DAYS_IN_YEAR = 365;

MONTHS_IN_YEAR = 12;


#=========================================================
#               Other Variables
#=========================================================

# Range of the latitude and longitude that are imported from OPeNDAP
LAT_RANGE = [250, 584];
LON_RANGE = [0, 550];
    
    
METHODNAMES = ('NPCC', 'Belcher', 'NPCCAdjusted');

  
# List of cities for NPCC project, with names of tmy3 files in .tmy2 format
CITY = [ "WASeattle3", "ORCorvallis3", "IDBoise3", "ORRedmond3", "NVElko3",
               "IDBurley3", "IDSodaSprings3", "MTHavre3", "MTMilesCity3",
               "ORPortland3", "WASpokane3", "MTKalispell3"];

TMY3NUMBER = ["727930", "726945", "726810", "726835", "725825", "725867", 
              "725868", "727770", "742300",
              "726980", "727850", "727790"] 

TMY2EXT = '.tm2'
TMY3EXT = 'TYA.CSV'



