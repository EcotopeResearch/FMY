# -*- coding: utf-8 -*-
"""
Created on Tue Jul  2 10:08:49 2019

@author: paul
"""

import pandas as pd;
from datetime import datetime, date, timedelta
from netCDF4 import Dataset
import numpy as np
import time
from cfg import VARNAME, VARLONGNAME, MODELNAME, SCENNAME, YEAR_START, YEAR_END, RUN_NUM, DOMAIN


#################################################################################
def get_data(df, scen, var, models, LAT_TARGETS, LON_TARGETS, stations, daily = 0, interp = False):
    
    # Iterate over the combination of stations scenarios and models chosen.
    
    for ss in stations:
        
        lat_targ = LAT_TARGETS[ss];   
        lon_targ = 360-LON_TARGETS[ss]; #should be in [0 360 ]E
    
        for sc in scen:
            
            for mm in models:
                print('model #: %d' % mm) 
                  
                ddfftemp = pd.DataFrame();
            
                for vv in var:
                    
                    dat = load_data(sc, vv, mm, lat_targ, lon_targ, daily, interp);
                    
                    if vv == 0 or vv == 1: #These are temperatures. Convert them to C
                        dat[VARNAME[vv]] = dat[VARNAME[vv]] - 273.15;
                        
                    ddfftemp[VARNAME[vv]] = dat[VARNAME[vv]];
                
                dl = len(dat['timehandle']);   
                ddff = pd.DataFrame( {'station':[ss]*dl, 
                                       'model':[mm]*dl,
                                       'scenario':[sc]*dl})
                ddff =  pd.concat([ddff, ddfftemp], axis=1);
                
                #Manage time formats   
                temp = from_1900_to_doy(dat['timehandle'], sc, daily);
                ddff['year']  = temp[0];
                ddff['month'] = temp[1];
                ddff['doy']   = temp[2];

                
                df = df.append(ddff);


    return(df);

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Load data
def load_data(scen, var, model, lat_target, lon_target, daily, interp):

    start = time.time()

    ll_range = 10;
    #=========================================================
    #               Load data at station
    #=========================================================
    
    dirPath='http://thredds.northwestknowledge.net:8080/thredds/dodsC/'  

    # Let's start with mean temperature from  min and max daily
    #--------------------------------------------------------
    #   FORM FILENAME AND GET HANDLE TO FILE AND LAT/LON 
    #--------------------------------------------------------
    Time=YEAR_START[scen]+'_'+YEAR_END[scen]
    if daily == 1:
        fileName = ('agg_macav2metdata_'+str(VARNAME[var])+'_'+str(MODELNAME[model])+'_r'+str(RUN_NUM[model])+'i1p1_'+str(SCENNAME[scen])+'_'+Time+'_'+DOMAIN+'_daily.nc')  
    elif daily == 0:
        fileName = ('agg_macav2metdata_'+str(VARNAME[var])+'_'+str(MODELNAME[model])+'_r'+str(RUN_NUM[model])+'i1p1_'+str(SCENNAME[scen])+'_'+Time+'_'+DOMAIN+'_monthly.nc')  

    fullfilename= dirPath+fileName

    # Sort data
    filehandle=Dataset(fullfilename,'r',format="NETCDF4")
    lathandle=filehandle.variables['lat']
    lonhandle=filehandle.variables['lon']
    timehandle=filehandle.variables['time']
    datahandle=filehandle.variables[VARLONGNAME[var]]
    
    print(fullfilename)
     
    #=========================================================
    #get data
    lat = lathandle[:]
    lon = lonhandle[:]
    #=========================================================
    #find indices of target lat/lon/day
    lat_index = (np.abs(lat-lat_target)).argmin()
    lon_index = (np.abs(lon-lon_target)).argmin()
    #check final is in right bounds
    if(lat[lat_index]>lat_target):
    	if(lat_index!=0):
    		lat_index = lat_index - 1
    if(lat[lat_index]<lat_target):
    	if(lat_index!=len(lat)):
    		lat_index =lat_index +1
    if(lon[lon_index]>lon_target):
    	if(lon_index!=0):
    		lon_index = lon_index - 1
    if(lon[lon_index]<lon_target):
    	if(lon_index!=len(lon)):
    		lon_index = lon_index + 1
    
    #=========================================================
    # If interpolating to the station or not.
    if interp:
        
        #adjust range for interpolation so that its not out of the index bounds (i.e. not a negative index)
        if lat_index - ll_range < 0:
            ll_range = lat_index;
        elif lat_index + ll_range >= len(datahandle[0,:,0]):
            ll_range = len(datahandle[0,:,0]) - lat_index;
        if lon_index - ll_range < 0:
            ll_range = lon_index
        elif lon_index + ll_range >= len(datahandle[0,0,:]):
            ll_range =  len(datahandle[0,0,:]) - lon_index
            
        #set lat and lon index such that they are a range for gridded interpolation
        lat_index = list(range(lat_index - ll_range, lat_index + ll_range))
        lon_index = list(range(lon_index - ll_range, lon_index + ll_range))
        
        lat = lat[lat_index]
        lon = lon[lon_index] 
        
        #=====================
        
        station_data = interp2station(lat, lon, lat_target, lon_target, datahandle[:,lat_index,lon_index]);
    
        #=====================
        
        ls = {'timehandle': timehandle,
              str(VARNAME[var]): station_data}
    
    else:
        lat = lat[lat_index]
        lon = lon[lon_index]
        
        ls = {'timehandle': timehandle,
              str(VARNAME[var]): datahandle[:,lat_index,lon_index]}
    
    print('Time to import: %.2f s' % (time.time()-start))
    
    return(ls);
    
    # import matplotlib.pyplot as plt;
   # fig, ax0 = plt.subplots(1)
   # im = ax0.pcolormesh(lon.data,lat.data,datahandle[2,lat_index,lon_index])
   # fig.colorbar(im, ax=ax0)

   # plt.xlim(235, 250)
    #plt.ylim(40, 49)
   # plt.show()
   
#=========================================================
#=========================================================
#=========================================================
    
#Convert months since 1900-01-01 00:00:00
def from_1900_to_doy(time1, scen, daily):
    leap_years = np.array([1904, 1908, 1912, 1916, 1920, 1952, 1956, 1960, 1964, 1968, 1972, 1976, 1980, 1984, 1988, 1992, 1996, 
                           2000, 2004, 2008, 2012, 2016, 2020, 2024, 2028, 2032, 2036, 2040, 2044, 2048, 
                           2052, 2056, 2060, 2064, 2068, 2072, 2076, 2080, 2084, 2088, 2092, 2096]);
    
    start   = date(1900,1,1);     # This is the "days since" part
    tempdoy = [0]*(len(time1));
    tempy   = [0]*(len(time1));
    tempm   = [0]*(len(time1));
        
    day = time1[:].data;
    
    # Account for leap years before the time array starts counting days.
    # scenario == 0 starting at 1950, scenario > 0 == 2006
    if scen == 0: 
        day = day + sum(leap_years < 1950); 
    elif scen > 0:   
        day = day + sum(leap_years < 2006);
    
    ind = 0;
    for ii in range(0,len(time1)):
        delta = timedelta(int(day[ii]));     # Create a time delta object from the number of days

        offset = start + delta;     # Add the specified number of days to 1990
        
        tempdoy[ind] = offset.timetuple().tm_yday;
        tempm[ind]   = offset.timetuple().tm_mon;
        tempy[ind]   = offset.timetuple().tm_year;    
        
        ind = ind + 1;
        
        if daily == 0 and sum(offset.timetuple().tm_year == leap_years) and offset.timetuple().tm_mon == 1:
            day = day + 1; # Add a day to account for leap years.
        
       # if daily == 1 and sum(offset.timetuple().tm_year == leap_years) and offset.timetuple().tm_mon == 2 and offset.timetuple().tm_mday == 29:
       #     ind = ind - 1; # If it was a leap day let's just rewrite it next step.
        
    return(tempy, tempm, tempdoy)
    
###############################################################################  
    
def interp2station(lat, lon, target_lat, target_lon, data):
    from scipy import interpolate;
    
    lent = len(data[:,0,0].data);
    data_t = np.zeros(lent);
    
    for tt in range(0,lent):
        f = interpolate.interp2d(lon, lat, data[tt,:,:].data, kind = "linear");
        
        data_t[tt] =  f(target_lon, target_lat)[0];
        
    return( data_t );
    
    
    