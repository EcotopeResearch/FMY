# -*- coding: utf-8 -*-
"""
Config file for globals used in the future weather project.

@author: paul
"""\


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

METHODNAMES = ('NPCC', 'Belcher', 'NPCC_Adjusted')
