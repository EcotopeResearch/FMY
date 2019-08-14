#This script was run using Python 2.7.13, Juypter Notebook from Anaconda Installation
#Author: Katherine Hegewisch (khegewisch@uidaho.edu) 
#Updated: 03/13/17
#Description: This script extracts a time series of data from the MACAv2-METDATA dataset from 1950-2099 for a single model, 
#but is easily generalizable to all models, all variables by looping over the models/variables in VARNAME, MODELNAME and 
#over a set of point locations by looping over LAT_TARGETS
#=========================================================
#            MACAV2METDATA FILE PARAMETERS
#=========================================================
dirPath='http://thredds.northwestknowledge.net:8080/thredds/dodsC/'  
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
RUN_NUM=[1] * 20
RUN_NUM[4] = 6 #setting CCSM4 with run 6
DOMAIN='CONUS'

#lat/lon pairs for point locations
LAT_TARGETS=[46.7317]
LON_TARGETS=[116.9972]
#=========================================================
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
#--------------------------------------------------------
#   MODEL/VAR/SCEN CHOICES
#--------------------------------------------------------
model=1
var=0
scen=1
#--------------------------------------------------------
#   FORM FILENAME AND GET HANDLE TO FILE AND LAT/LON 
#--------------------------------------------------------
Time=YEAR_START[scen]+'_'+YEAR_END[scen]
fileName = ('agg_macav2metdata_'+str(VARNAME[var])+'_'+str(MODELNAME[model])+'_r'+str(RUN_NUM[model])+'i1p1_'+str(SCENNAME[scen])+'_'+Time+'_'+DOMAIN+'_daily.nc')  
fullfilename= dirPath+fileName
print(fullfilename)
filehandle=Dataset(fullfilename,'r')
lat=filehandle.variables['lat']
lon=filehandle.variables['lon']

#--------------------------------------------------------
#   LAT/LON CHOICES
#--------------------------------------------------------
#this is how you can loop over the point locations
#for station in enumerate(LAT_TARGETS):
station=0
lat_targ = LAT_TARGETS[station]   
lon_targ=360-LON_TARGETS[station] #should be in [0 360 ]E

#--------------------------------------------------------
#   LAT/LON DATA
#--------------------------------------------------------
lathandle=filehandle.variables['lat']
lonhandle=filehandle.variables['lon']
lat = lathandle[:]
lon = lonhandle[:]

#--------------------------------------------------------
#  GET LAT/LON NEAREST SELECTED LAT/LON
#--------------------------------------------------------
lat_index =   (np.abs(lat-lat_targ)).argmin() 
lon_index =   (np.abs(lon-lon_targ)).argmin()

if(lat[lat_index]>lat_targ):
    if(lat_index!=0):
        lat_index = lat_index - 1
if(lon[lon_index]>lon_targ):
    if(lon_index!=0):
        lon_index = lon_index - 1 
print(lat[lat_index] )
print(lon[lon_index] -360)



#this is how you can loop over the MODELNAMEs, VARNAMES
#for var,varName in enumerate(VARNAME):
var=0  #tasmax
#for model,modelName in enumerate(MODELNAME):
model=0 #bcc-csm-1-1
#--------------------------------------------------------
#   GET HISTORIC DATA
#--------------------------------------------------------
scen=0
Time=YEAR_START[scen]+'_'+YEAR_END[scen]
fileName = ('agg_macav2metdata_'+str(VARNAME[var])+'_'+str(MODELNAME[model])+'_r'+str(RUN_NUM[model])+'i1p1_'+str(SCENNAME[scen])+'_'+Time+'_'+DOMAIN+'_daily.nc')   
fullfilename= dirPath+fileName 
filehandle=Dataset(fullfilename,'r')
datahandle=filehandle.variables[VARLONGNAME[var]]
histdata = datahandle[:,lat_index,lon_index]
#--------------------------------------------------------
#   GET FUTURE DATA
#--------------------------------------------------------  
scen=1
Time=YEAR_START[scen]+'_'+YEAR_END[scen]
fileName = ('agg_macav2metdata_'+str(VARNAME[var])+'_'+str(MODELNAME[model])+'_r'+str(RUN_NUM[model])+'i1p1_'+str(SCENNAME[scen])+'_'+Time+'_'+DOMAIN+'_daily.nc')   
fullfilename= dirPath+fileName 
filehandle=Dataset(fullfilename,'r')
datahandle=filehandle.variables[VARLONGNAME[var]]
futdata = datahandle[:,lat_index,lon_index] 