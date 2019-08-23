# -*- coding: utf-8 -*-
"""
Methods to get data from MACA catalogs using the OPenDAP frame works.

Will retrieve data for a list of stations, a list of scenarios, a list of 
variables, a list of models. And will either get daily data or monthly data
from MACA outputs.

interp controls the interpolation to the lat/lon cooridinates of the stations
vs just using the closest grid point.

@author: paul
"""
import os;
import pandas as pd;
from datetime import datetime, date, timedelta
from netCDF4 import Dataset
import numpy as np
import time
import json
from cfg import VARNAME, VARLONGNAME, MODELNAME, SCENNAME, YEAR_START, YEAR_END, RUN_NUM, DOMAIN, LAT_RANGE, LON_RANGE, CITY


#################################################################################
def get_data(df, scen, var, models, LAT_TARGETS, LON_TARGETS, stations, datapath, daily = 0, interp = True, download_data = False):
    station_inds = list(range(0,len(stations)))
       
    # Convert to numpy arrays so we can do maths
    LAT_TARGETS = np.array(LAT_TARGETS);   
    LON_TARGETS = 360-abs(np.array(LON_TARGETS)); #should be in [0 360 ]E 
    
    
    # Iterate over the combination of stations scenarios and models chosen.
    for sc in scen:
        
        for mm in models:
            dfstack = pd.DataFrame();
            
            for vv in var:
                ddffvariables = pd.DataFrame();
    
                #Get data for all the stations off MACA with OPenDAP
                dat = load_data(sc, vv, mm, LAT_TARGETS, LON_TARGETS, stations, daily, interp, datapath, download_data);
                  
                
                for ss in station_inds:
                    ddffstations  = pd.DataFrame();
                    
                    if vv == 0 or vv == 1: #These are temperatures. Convert them to C
                         dat[ss][VARNAME[vv]] = dat[ss][VARNAME[vv]] - 273.15;
                    
                    ddffstations[VARNAME[vv]] = dat[ss][VARNAME[vv]];
            
                    # If first time through the variables manage the data frame setup
                    if vv == var[0]:
                        
                        
                        dl = len(dat[ss]['timehandle']);   
                        ddffinputs = pd.DataFrame( {'station':[stations[ss]]*dl, 
                                                    'model':[mm]*dl,
                                                    'scenario':[sc]*dl})
                        
                        # Manage time formats if first time through 
                        temp = from_1900_to_doy(dat[ss]['timehandle'], sc, daily);
                        ddffinputs['year']        = temp[0];
                        ddffinputs['month']       = temp[1];
                        ddffinputs['doy']         = temp[2];
                        
                        ddffstations = pd.concat([ddffinputs, ddffstations], axis=1)
                        
                        # Data frame with all the stations and just one variable.      
                        ddffvariables = pd.concat([ddffvariables, ddffstations], axis=0).reset_index(drop = True);
                        
                    else:
                        #Append the variable as a column
                        ddffvariables = pd.concat([ddffvariables, ddffstations], axis=0).reset_index(drop = True)
                    
                    # End of for stations
                    
                
                dfstack = pd.concat([dfstack, ddffvariables], axis=1);
                #End of for variables
                
            df = df.append(dfstack);   
            #End of for models
        #End of for scenarios
    


    return(df);

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    

#Load data
def load_data(scen, var, model, lat_targets, lon_targets, stations, daily, interp, datapath, dl_dat ):

    start = time.time()

    ll_range = 15;
    #=========================================================
    #               Load data at station
    #=========================================================
    
    dirPath='http://thredds.northwestknowledge.net:8080/thredds/dodsC/'  

    # Let's start with mean temperature from  min and max daily
    #--------------------------------------------------------
    #   FORM FILENAME AND GET HANDLE TO FILE AND LAT/LON 
    #--------------------------------------------------------
    Time=YEAR_START[scen]+'_'+YEAR_END[scen]
    
    # fileName = "agg_macav2metdata_rsds_MIROC5_r1i1p1_rcp85_2006_2099_CONUS_monthly.nc
    # ?lat[0:1:584],lon[0:1:1385],crs[0:1:0],time[0:1:1127],surface_downwelling_shortwave_flux_in_air[0:1:1127][0:1:584][0:1:1385]"
    # LAT_RANGE = [250, 584],
    # LON_RANGE = [0, 550],
    
    fileName = ('agg_macav2metdata_' + str(VARNAME[var]) + '_' + 
                str(MODELNAME[model]) + '_r' + str(RUN_NUM[model]) + 'i1p1_' +
                str(SCENNAME[scen]) + '_' + Time + '_' + DOMAIN ) ;
                 
    if daily == 1: # Load the daily data
        time_index = ['20453', '34332', '34332'];
        fileName = fileName + '_daily.nc';
        part2    = ('?lat[' +  str(LAT_RANGE[0]) + ':1:' + str(LAT_RANGE[1]) +
                    '],lon[' + str(LON_RANGE[0]) + ':1:' + str(LON_RANGE[1]) + 
                    '],crs[0:1:0],time[0:1:' + time_index[scen] + '],'
                    + VARLONGNAME[var] + '[0:1:' + time_index[scen] + '][' +
                    str(LAT_RANGE[0]) + ':1:' + str(LAT_RANGE[1]) +
                    '][' + str(LON_RANGE[0]) + ':1:' + str(LON_RANGE[1]) + ']' );
                    
    elif daily == 0: # Load the monthly data
        time_index = ['671', '1127', '1127'];
        fileName = fileName + '_monthly.nc';
        part2    = ('?lat[' + str(LAT_RANGE[0]) + ':1:' + str(LAT_RANGE[1]) +
                    '],lon[' + str(LON_RANGE[0]) + ':1:' + str(LON_RANGE[1]) + 
                    '],crs[0:1:0],time[0:1:' + time_index[scen] + '],'
                    + VARLONGNAME[var] + '[0:1:' + time_index[scen] + '][' +
                    str(LAT_RANGE[0]) + ':1:' + str(LAT_RANGE[1]) +
                    '][' + str(LON_RANGE[0]) + ':1:' + str(LON_RANGE[1]) + ']' );
                    
    else: print('ERROR: incorrect value for daily used must be 0, or 1');
    
    fullfilename = dirPath + fileName + part2

    # Check locally for the GCM data else go look on the MACA site.
    sst = False;
    for ss in range(0,len(stations)): 
        if not os.path.exists(datapath + CITY[stations[ss]] + fileName):
            sst = True;


    # If there's one station missing data load it
    lst = [];
    if sst:   
        print('Loading... ' + fullfilename)

        # Get data from the MACA address           
        filehandle = Dataset(fullfilename,'r',format="NETCDF4")     
            
        # Sort data
        lathandle   = filehandle.variables['lat']
        lonhandle   = filehandle.variables['lon']
        timehandle  = filehandle.variables['time']
        datahandle  = filehandle.variables[VARLONGNAME[var]]
        
        #=========================================================
        # get lat and lon in a usable form
        lat = lathandle[:]
        lon = lonhandle[:]
        
        #=========================================================
        # get the data for the different stations
        for ss in range(0,len(stations)): 
            lat_target  = lat_targets[ss] 
            lon_target  = lon_targets[ss]
            
            #=========================================================
            # check the target is inbounds
            if not (min(lat) <= lat_target <= max(lat)):
                print('ERROR: The target station, number '+ str(ss) +' is out of the lat bounds, consider changing the station or the index bounds in cfg.py')
            if not (min(lon) <= lon_target <= max(lon)):
                print('ERROR: The target station, number '+ str(ss) +' is out of the lon bounds, consider changing the station or the index bounds in cfg.py')
    
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
                    lat_range = 0;
                elif lat_index + ll_range >= len(datahandle[0,:,0]):
                    lat_range = len(datahandle[0,:,0]) - lat_index - 1;
                else:
                    lat_range = ll_range;
                if lon_index - ll_range < 0:
                    lon_range = 0
                elif lon_index + ll_range >= len(datahandle[0,0,:]):
                    lon_range =  len(datahandle[0,0,:]) - lon_index - 1
                else:
                    lon_range = ll_range;
                #set lat and lon index such that they are a range for gridded interpolation
                lat_index = list(range(lat_index - lat_range, lat_index + lat_range))
                lon_index = list(range(lon_index - lon_range, lon_index + lon_range))
                
                lat_interp = lat[lat_index]
                lon_interp = lon[lon_index] 
                
                #=====================
                
                station_data = interp2station(lat_interp, lon_interp, lat_target, lon_target, 
                                              datahandle[:,lat_index,lon_index]);
            
                #=====================
                
                lst.append({'timehandle': timehandle[:].data,
                      str(VARNAME[var]): station_data})
            
                if dl_dat: # Write to file if downloading data.
                    with open(datapath+CITY[stations[ss]]+fileName, 'w') as out_file:
                            json.dump({'timehandle': timehandle[:].data.tolist(),
                                       str(VARNAME[var]): station_data.tolist()}, out_file)

            else: #on if interp:
                lst.append({'timehandle': timehandle[:].data,
                      str(VARNAME[var]): datahandle[:,lat_index,lon_index]})
        
                if dl_dat: # Write to file if downloading data.
                    with open(datapath+CITY[stations[ss]]+fileName, 'w') as out_file:
                            json.dump({'timehandle': timehandle[:].data.tolist(),
                                       str(VARNAME[var]): station_data.tolist()}, out_file) 
    
    # All the data is on the computer
    else: 
        for ss in range(0,len(stations)): 
            fullfilename = datapath+CITY[stations[ss]]+fileName;
            
            print('Loading... ' + fullfilename)

            with open(fullfilename, 'r') as in_file:
                temp = json.load(in_file);
                lst.append({'timehandle': np.array(temp['timehandle']),
                      str(VARNAME[var]): np.array(temp[VARNAME[var]])})
        
    print('Time to import: %.2f s\n' % (time.time()-start))
        
    return(lst);
    
    
#    import matplotlib.pyplot as plt;
#    fig, ax0 = plt.subplots(1)
#    im = ax0.pcolormesh(lon.data,lat.data,datahandle[2,:,:])
#    im = ax0.pcolormesh(datahandle[2,:,:])
#    fig.colorbar(im, ax=ax0)
#
#    plt.xlim(235, 250)
#    plt.ylim(40, 49)
#    plt.show()
   
#=========================================================
#=========================================================
#=========================================================
 
def write_file(file, src):
    
    #output file
    dst = Dataset(file, "w", format="NETCDF4")
    
    # copy global attributes all at once via dictionary
    dst.setncatts(src.__dict__)
    # copy dimensions
    for name, dimension in src.dimensions.items():
        dst.createDimension(
            name, (len(dimension) if not dimension.isunlimited() else None))
    # copy all file data 
    for name, variable in src.variables.items():
        x = dst.createVariable(name, variable.datatype, variable.dimensions)
        dst[name][:] = src[name][:]
        # copy variable attributes all at once via dictionary
        #dst[name].setncatts(src[name].__dict__)
        
    # close the output file
    dst.close()
    
###############################################################################    

###############################################################################
                            
#Convert months since 1900-01-01 00:00:00
def from_1900_to_doy(time1, scen, daily):
    leap_years = np.array([1904, 1908, 1912, 1916, 1920, 1952, 1956, 1960, 1964, 1968, 1972, 1976, 1980, 1984, 1988, 1992, 1996, 
                           2000, 2004, 2008, 2012, 2016, 2020, 2024, 2028, 2032, 2036, 2040, 2044, 2048, 
                           2052, 2056, 2060, 2064, 2068, 2072, 2076, 2080, 2084, 2088, 2092, 2096]);
    
    start   = date(1900,1,1);     # This is the "days since" part
    tempdoy = [0]*(len(time1));
    tempy   = [0]*(len(time1));
    tempm   = [0]*(len(time1));
        
    day = time1;
    
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
    
    
    