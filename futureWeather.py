"""
#Author:	Paul Kintner
#Updated: 	8/14/2019
#Description: 	This script uses OPeNDAP to download the specified subset of 
the MACAv2-METDATA data then applies multiple methods for temporally downscalling 
data to hourly features.

MACA advocates for using at least 10 models in analyses. It may be worth the effort 
to evaluate the relevant variables against observations, just to be cognizant of 
model biases, but recognize that most studies have found little or no difference 
in culling or weighting model outputs.

The important outputs changed in the new future TMY file are:
1.	Total horizontal solar Btu/h ft^2
2.	Direct normal solar: Btu/h ft^2
3.	Diffuse Horizontal Solar: Btu/h ft^2
4.	Tdrybulb
5.	Relative Humidity

and maybe
6.	Windspeed
7.	Wind direction azimuth
8.	Could cover fraction

And prints into TMY3  and TMY2 formatting

Takes Cities: 
    0. Seattle (WA)
    1. Corvallis (OR)
    2. Boise (ID)
    3. Redmond(OR)
    4. Elko (NV)
    5. Burley (ID)
    6. Soda_Springs (ID)
    7. Havre (MT)
    8. Miles City (MT)
 
city        = [ "WASeattle3", "ORCorvallis3", "IDBoise3", "ORRedmond3", "NVElko3",
               "IDBurley3", "IDSodaSprings3", "MTHavre3", "MTMilesCity3" ];
   
workflow:
    Set specifics (i.e. variables, cities, models, scenarios, and methods)
    Load tmy2 file
    Load GCM historical and future
    Determine future peroid of interest and cull the data
    Adjust the variables asked to be adjusted
    Write output file in specified formats
"""

#=========================================================
import os
import numpy as np;
import pandas as pd;
import matplotlib.pyplot as plt;
import seaborn as sns;

workingdir = 'R:/NPCC/2019_Future_Meteorological_Years (19-028-BL)/Python'

os.chdir(workingdir)

from cfg import METHODNAMES, VARNAME, MODELNAME, SCENNAME,HRS_IN_YEAR;

from OPeNDAP_load import get_data;
from variable_conversions import get_monthly_mmms;
from weather import weather; 
from bias_correction import bias_correct;
from tmy_to_fmy import adjust_to_fmy;
#
#import metpy.calc as mpcalc
#from metpy.units import units

#=========================================================

# Paths
weatherpath = "E:/TMY2DATA/"; # Path to local TMY data in .tm2 format
graphpath = "./Graphs/";      # Output path for graphs
outputpath = "./output_FMY/"  # Output path for future .tm2 and .tm3 files

# Output formats
outformats = ['csv','tmy2']; # List containing ['csv','tmy2', 'tmy3'], and any combination of them or none.

# Takes Cities: 
#    0. Seattle (WA)
#    1. Corvallis (OR)
#    2. Boise (ID)
#    3. Redmond(OR)
#    4. Elko (NV)
#    5. Burley (ID)
#    6. Soda_Springs (ID)
#    7. Havre (MT)
#    8. Miles City (MT)
stations    = [ 0, 1, 2, 3, 4, 5, 6, 7, 8 ];
stations    = [ 0, 6 ]; 

# Affects plots and loads
method      = 2; # 1 is NPCC, 
                 # 2 is Belcher  
                 # 3 is adjusted NPCC with monthly means of daily max and min
                 
# Set to 'tmy' or 'gcm' to set which current climate is used to calculate the ajustments to the mean for            
which_current_climate = 'gcm';

interpolate_to_station = True

bias_correction_method = 0; #0 no correction / 1 is EDCDFm / 2 should be CDFm / 3 is linear Hawkins

hourly_plots = [0, 1, 2, 3, 5, 8];

suppress_all_plots = 1;
plt.ioff()

#--------------------------------------------------------
#   MODEL/VAR/SCEN CHOICES
#--------------------------------------------------------

# bcc-csm1-1 (0) did not have daily data available for 12/31/2099 for the RCP8.5 scenario only.
# CCSM4 (4) and NorESM1-M (19) did not have relative humidity available at daily timescales.
models      = [1, 2, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18] # All complete sets
models      = [3, 4, 5, 6, 8, 9, 10, 11, 12, 15] #The RMJOC-II “10”. Note issues with model 4 in scenario=1
models      = [3] # for testing

scen        = [ 1 ] # 1 is RCP4.5, 2 is RCP8.5

var         = [0, 1, 2, 3, 5, 8]
#var         = [ 0, 1 ]

# Years to average over
tmy3_years   = [ 1976, 2005 ]; # Years that the tmy3 weather files are taken from. 
future_years = [ 2020, 2049 ]; # Future years from 2006 - 2099.


#--------------------------------------------------------
#   Initialize other things.
#--------------------------------------------------------

# List of cities
city        = [ "WASeattle3", "ORCorvallis3", "IDBoise3", "ORRedmond3", "NVElko3",
               "IDBurley3", "IDSodaSprings3", "MTHavre3", "MTMilesCity3" ];
               
# Get the indexs for talking to the station list
station_inds = list(range(0,len(stations)))

#lat/lon pairs for point locations
LAT_TARGETS=[0] * len(station_inds) #wd.lat
LON_TARGETS=[0] * len(station_inds)# wd.lon

# Data frame for tmy data.
dft     = pd.DataFrame();
# Data frame for historic model data.
df      = pd.DataFrame();
# Data frame for predicted model data.
df1     = pd.DataFrame();
# Data frame for FMY data.
df_fmy  = pd.DataFrame([]);

# Switch Data Collection based on method
if method == 1:
    daily       = 1; # Gets daily data (1) or monthly data (0).
    daily_max   = 0; # For stats (0) gets the max/min for the whole month, 1 gets the monthly mean of the daily max/min
elif method == 2:
    daily       = 0; # Gets daily data (1) or monthly data (0).
    daily_max   = 1; # For stats (0) gets the max/min for the whole month, 1 gets the monthly mean of the daily max/min
elif method == 3:
    daily       = 0; # Gets daily data (1) or monthly data (0).
    daily_max   = 1; # For stats (0) gets the max/min for the whole month, 1 gets the monthly mean of the daily max/min
    
#--------------------------------------------------------
#   Check Directories.
#--------------------------------------------------------

os.makedirs(graphpath, exist_ok=True)
os.makedirs(outputpath, exist_ok=True)

#===============================================================================
#===============================================================================
#                               Load Data!
#===============================================================================
#===============================================================================
print('Loading TMY data...')
# Get tmy weather file
tmy = [0] * len(station_inds);
for ss in station_inds:
    tmy[ss] = weather(weatherpath, city[stations[ss]]);
    tmy[ss].get_weather();
    
    dft = pd.concat( [ dft, tmy[ss].write_to_monthly_df(stations[ss], daily_max) ],  axis=0); #Stack stations on top of each other.
    
    LAT_TARGETS[ss] = tmy[ss].lat;
    LON_TARGETS[ss] = tmy[ss].lon;
    
    print('Loaded TMY dataset of ' + city[ss] + '.\n')

print('Loading GCM data...')
# Get historical data, scenario = 0
dfhist = get_data(df, [0], var, models, LAT_TARGETS, LON_TARGETS, stations, daily, interpolate_to_station);

# and get future data for a scenario specified above
dffut = get_data(df1, scen ,var, models, LAT_TARGETS, LON_TARGETS, stations, daily, interpolate_to_station);

#scatter_matrix(df)
#plt.show()

#=========================================================
#               Cut data to study "area"
#=========================================================
"""
Cut to look at specific years for models
"""
#tmy3_years   = [ 1990, 2005 ]; # Years that the tmy3 weather files are taken from.
dfhist = dfhist[ dfhist.year >= tmy3_years[0]];
dfhist = dfhist[ dfhist.year <= tmy3_years[1]];

#future_years = [ 2060, 2090 ]; # Future years from 2006 - 2099.
dffut = dffut[ dffut.year >= future_years[0]];
dffut = dffut[ dffut.year <= future_years[1]];


#=========================================================
#               Bias correct station data
#=========================================================
"""
The MACA dataset is bias correct but they do recommend 
doing another bias correction to a specific station to avoid large geographic 
features.
"""
# Select bias selection method
if bias_correction_method > 0 and bias_correction_method <= 3:
    tempdf1 = pd.DataFrame()
    for ss in stations:
        #Get hourly tmy data for bias correction
        dft_temp = tmy[ss].write_to_daily_df(ss);
        
        for vv in var:
           # if vv == 0 or vv == 1:
                tempdf1[str(VARNAME[vv]+'_c')] = bias_correct( 
                        dfhist[ (dfhist['station']==ss) ][VARNAME[vv]],
                        dffut[ (dffut['station']==ss)] [VARNAME[vv]], 
                        dft_temp[ (dft_temp['station']==ss) ][VARNAME[vv]], 
                        vv, bias_correction_method)
           # if vv == 2 or vv == 3:
           #     df1[str(VARNAME[vv]+'_c')] = df1[VARNAME[vv]];

    dffut = pd.concat([dffut, tempdf1], axis=1);
        
elif bias_correction_method == 0: #The corrected variables are just the variables
        for vv in var:
            dffut[str(VARNAME[vv]+'_c')] = dffut[VARNAME[vv]];

else: print("\nWARNING: Invalid bias_correction_method chose can be 0,1,2,3, continuing without a correction.\n")

# Give the historical data frame that same name for their values so we can use one fuction for both
for vv in var:
    dfhist[str(VARNAME[vv]+'_c')] = dfhist[VARNAME[vv]];
            
#=========================================================
#                Check aggregation Method
#=========================================================
dfhist    = get_monthly_mmms(dfhist, var, daily_max);
dffut     = get_monthly_mmms(dffut, var, daily_max);

#=========================================================
#                Mean Temperature
#=========================================================
if 0 in var or 1 in var:
    if 0 in var and 1 in var:
        dfhist['tasmean'] = ( dfhist['tasmax_c'] + dfhist['tasmin_c'] ) / 2.;

        dffut['tasmean'] = ( dffut['tasmax_c'] + dffut['tasmin_c'] ) / 2.;
    else:
        print("\nWARNING: Selecting the max/min temperature without the other add 0/1 to var")  
        
#=========================================================
#                Mean Relative Humidity
#=========================================================
if 2 in var or 3 in var:
    if 2 in var and 3 in var:
        dfhist['rhsmean'] = ( dfhist['rhsmax_c'] + dfhist['rhsmin_c'] ) / 2.;

        dffut['rhsmean'] = ( dffut['rhsmax_c'] + dffut['rhsmin_c'] ) / 2.;
    else:
        print("\nWARNING: Selecting the max/min relative humidity without the other add 2/3 to var")    

#=========================================================
#                Find Month Closest to Average
#=========================================================
df        = dfhist.copy(); # New data frames to transform
df1       = dffut.copy();  # New data frames to transform

# For NPCC method first find the most average month for all variables
# and then find the constants for transformation.
df1['err']   = 0; #Initialize "error" column 

for vv in var:             
    # Grab the monthly means of a variable 
    monthly_mean = df1.groupby(['station', 'month', 'model', 'scenario'])[VARNAME[vv]+'_c'].agg(['mean'])
    
    # Rearange the data frame to merge into tempdf
    monthly_mean = monthly_mean.reset_index();
    monthly_mean = monthly_mean.rename(index=str, columns={"mean" : VARNAME[vv]+"_mean"})
    
    # Have to merge monthly_mmm with tempdf
    df1=pd.merge(df1, monthly_mean)
    
    #Check the normalized vector distance for every row in tempdf
    df1['err'] = df1['err'] + ( (df1[VARNAME[vv]+'_c'] - df1[VARNAME[vv]+'_mean']) / df1[VARNAME[vv]+'_mean'] ) **2;

# Final transmoration for the errors        
df1['err'] = np.sqrt(df1['err']);              
    
df1['avg_month'] = df1.groupby(['station', 'month', 'model', 'scenario'])['err'].transform(min) == df1['err']

dfam = (df1[ df1['avg_month'] == True ].sort_values(by='month')).reset_index(drop=True);

# Quick Check on the unique number of months in dfam
if not (sum(df1.groupby('month').nunique().station == [len(stations)]*12) == 12):
    print('\nWARNING: Each station was not given an average month for each month of the year')

###############################################################################
#           Historical Plot
###############################################################################     
if not suppress_all_plots:
    for mm in models:
        for ss in stations:
            for sc in scen:
                for vv in var:
                    # Plot temperatures
                    if vv == 1:
                        fig = plt.figure()
                        fig.suptitle('Historical Monthly Max and Min Temperatures for \n'+ city[ss] + ' for ' + 
                             MODELNAME[mm] + ' from ' + str(tmy3_years[0]) +' - ' + 
                             str(tmy3_years[1]) + '\n running scenario: ' + SCENNAME[sc])
                        ax = fig.add_subplot(111)
                        sns.set_style("ticks")                     
                        sns.stripplot(x='month', y='tasmax_c', data=df[(df['station']==ss) & (df['model']==mm)], jitter=0.25, color = 'royalblue', alpha=.60, marker = 's', edgecolor="gray"  )
                        sns.stripplot(x='month', y='tasmin_c', data=df[(df['station']==ss) & (df['model']==mm)], jitter=0.25, color = 'royalblue',alpha=.50, marker = 'o')
                        sns.stripplot(x='month', y='tasmax', data=dft[dft['station']==ss], jitter=False,alpha=.90, color = 'black',marker="D")
                        sns.stripplot(x='month', y='tasmin', data=dft[dft['station']==ss], jitter=False, alpha=.90,color = 'black',marker="D")
                        ax.grid()
                        ax.set_xlabel(u'Month')
                        ax.set_ylabel(u'Max / Min Temperature (degC)')

                        fig.savefig(graphpath + city[ss] + '/' + city[ss]+'_'+METHODNAMES[method-1]+'_Baseline_'+which_current_climate+'_Historical_Monthly_Mean_of_Temperatures')
                        plt.close()
                        
                    # Plot relative humidity
                    if vv == 3:
                        fig = plt.figure()
                        fig.suptitle('Historical Monthly Max and Min Relative Humidity for\n'+ city[ss]+ ' for ' + 
                             MODELNAME[mm] + ' from ' + str(tmy3_years[0]) +' - ' + 
                             str(tmy3_years[1]) + '\n running scenario: ' + SCENNAME[sc])
                        ax = fig.add_subplot(111)
                        sns.set_style("ticks")                     
                        sns.stripplot(x='month', y='rhsmax_c', data=df[(df['station']==ss) & (df['model']==mm)], jitter=0.25, color = 'royalblue', alpha=.60, marker = 's', edgecolor="gray"  )
                        sns.stripplot(x='month', y='rhsmin_c', data=df[(df['station']==ss) & (df['model']==mm)], jitter=0.25,  color = 'royalblue', alpha=.60, marker = 'o')
                        sns.stripplot(x='month', y='rhsmax', data=dft[dft['station']==ss], jitter=False, alpha=.90, color = 'black',marker="D")
                        sns.stripplot(x='month', y='rhsmin', data=dft[dft['station']==ss], jitter=False,  alpha=.90,color = 'black',marker="D")
                        ax.set_ylim(0,110)
                        ax.grid() 
                        ax.set_xlabel(u'Month')
                        ax.set_ylabel(u'Max / Min Relative Humidity (%)')

                        fig.savefig(graphpath + city[ss] + '/' + city[ss]+'_'+METHODNAMES[method-1]+'_Baseline_'+which_current_climate+'_Historical_Monthly_Mean_of_Relative_Humidity')
                        plt.close()
                    # Plot specific humidity
                    if vv == 8:
                        fig = plt.figure()
                        fig.suptitle('Historical Monthly Specific Humidity for\n'+ city[ss]+ ' for ' + 
                             MODELNAME[mm] + ' from ' + str(tmy3_years[0]) +' - ' + 
                             str(tmy3_years[1]) + '\n running scenario: ' + SCENNAME[sc])
                        ax = fig.add_subplot(111)
                        sns.set_style("ticks")                     
                        sns.stripplot(x='month', y='huss_c', data=df[(df['station']==ss) & (df['model']==mm)], jitter=0.25, color = 'royalblue', alpha=.60, marker = 's', edgecolor="gray"  )
                        sns.stripplot(x='month', y='huss', data=dft[dft['station']==ss], jitter=False,  alpha=.90,color = 'black',marker="D")
                        ax.grid() 
                        ax.set_xlabel(u'Month')
                        ax.set_ylabel(u'Specific Humidity (%)')
                        
                        fig.savefig(graphpath + city[ss] + '/' + city[ss]+'_'+METHODNAMES[method-1]+'_Baseline_'+which_current_climate+'_Historical_Monthly_Mean_of_Specific_Humidity')
                        plt.close() 
                    # Plot Solar Terms
                    if vv == 5:
                        fig = plt.figure()
                        fig.suptitle('Historical Monthly Total Horizontal Radiation for\n'+ city[ss]+ ' for ' + 
                             MODELNAME[mm] + ' from ' + str(tmy3_years[0]) +' - ' + 
                             str(tmy3_years[1]) + '\n running scenario: ' + SCENNAME[sc])
                        ax = fig.add_subplot(111)
                        sns.set_style("ticks")                     
                        sns.stripplot(x='month', y='rsds_c', data=df[(df['station']==ss) & (df['model']==mm)], jitter=0.25, color = 'royalblue', alpha=.60, marker = 's', edgecolor="gray"  )
                        sns.stripplot(x='month', y='rsds', data=dft[dft['station']==ss], jitter=False,  alpha=.90,color = 'black',marker="D")
                        ax.grid() 
                        ax.set_xlabel(u'Month')
                        ax.set_ylabel(u'Total Horizontal Radiation (Wh/m2)')
                        fig.savefig(graphpath + city[ss] + '/' + city[ss]+'_'+METHODNAMES[method-1]+'_Baseline_'+which_current_climate+'_Historical_Monthly_Mean_of_Global_Radiation')
                        plt.close() 

############################################################################
#               Future plots
############################################################################
if not suppress_all_plots:
    for mm in models:
        for ss in stations:
            for sc in scen:
                dfplt = df1[(df1['station']==ss) & (df1['model']==mm) & (df1['scenario']==sc)]
                dfhistplot =  df[(df['station']==ss) & (df['model']==mm)]   ;              
                for vv in var:
                    #Plot Temperatures
                    if vv == 1:
                        fig = plt.figure()
                        fig.suptitle('Future Monthly Max and Min Temperatures for \n'+ city[ss]+ ' for ' + 
                                     MODELNAME[mm] + ' from ' + str(future_years[0]) +' - ' + 
                                     str(future_years[1]) + '\n running scenario: ' + SCENNAME[sc])
                        ax = fig.add_subplot(111)
                        sns.set_style("ticks")
                        
                        sns.stripplot(x='month', y='tasmax_c', data=dfhistplot, jitter=0.33, color = 'purple', alpha=.45, marker = '^' )
                        sns.stripplot(x='month', y='tasmin_c', data=dfhistplot, jitter=0.33, color = 'purple', alpha=.45, marker = 'v' )

                        sns.stripplot(x='month', y='tasmax_c', data=dfplt, jitter=.33, hue = 'avg_month', alpha=.60, marker = 's' )
                        sns.stripplot(x='month', y='tasmin_c', data=dfplt, jitter=.33, hue = 'avg_month', alpha=.60, marker = 'o')
                        
                        sns.stripplot(x='month', y='tasmax', data=dft[(dft['station']==ss)], jitter=False,  alpha=.90, color = 'black',marker="D")
                        sns.stripplot(x='month', y='tasmin', data=dft[(dft['station']==ss)], jitter=False, alpha=.90,color = 'black',marker="D")
                        ax.grid()
                        ax.set_xlabel(u'Month')
                        ax.set_ylabel(u'Max / Min Temperature (degC)')
                        plt.show()
                        
                        fig.savefig(graphpath + city[ss] + '/' + city[ss]+'_'+METHODNAMES[method-1]+'_Baseline_'+which_current_climate+'_Future_Monthly_Mean_of_Temperature')
                        plt.close()
                    # Plot Relative humidities
                    if vv == 3:
                        fig = plt.figure()
                        ax = fig.add_subplot(111)
                        fig.suptitle('Future Monthly Max and Min Relative Humidity for\n'+ city[ss]+ ' for ' + 
                                     MODELNAME[mm] + ' from ' + str(future_years[0]) +' - ' + 
                                     str(future_years[1]) + '\n running scenario: ' + SCENNAME[sc])
                        sns.set(style="ticks")
                
                        sns.stripplot(x='month', y='rhsmax_c', data=dfplt,jitter=.25, hue = 'avg_month', alpha=.60, marker = 's', edgecolor="gray"  )
                        sns.stripplot(x='month', y='rhsmin_c', data=dfplt, jitter=.25, hue = 'avg_month', alpha=.60, marker = 'o')
                        sns.stripplot(x='month', y='rhsmax', data=dft[(dft['station']==ss)], jitter=False,  alpha=.90, color = 'black', marker="D")
                        sns.stripplot(x='month', y='rhsmin', data=dft[(dft['station']==ss)], jitter=False,  alpha=.90, color = 'black', marker="D")
                        ax.set_ylim(0,110)
                        ax.grid()       
                        ax.set_xlabel(u'Month')
                        ax.set_ylabel(u'Max / Min Relative Humidity (%)')
                        
                        fig.savefig(graphpath + city[ss] + '/' + city[ss]+'_'+METHODNAMES[method-1]+'_Baseline_'+which_current_climate+'_Future_Monthly_Mean_of_Relative_Humidity')
                        plt.close()
                    # Plot Solar Terms
                    if vv == 8:
                        fig = plt.figure()
                        ax = fig.add_subplot(111)
                        fig.suptitle('Future Monthly Specific Humidity for\n'+ city[ss]+ ' for ' + 
                                     MODELNAME[mm] + ' from ' + str(future_years[0]) +' - ' + 
                                     str(future_years[1]) + '\n running scenario: ' + SCENNAME[sc])
                        sns.set(style="ticks")
                
                        sns.stripplot(x='month', y='huss_c', data=dfplt, jitter=.25, hue = 'avg_month', alpha=.60, marker = 'o')
                        sns.stripplot(x='month', y='huss', data=dft[(dft['station']==ss)], jitter=False,  alpha=.90, color = 'black', marker="D")
                        ax.grid()
                        ax.set_xlabel(u'Month')
                        ax.set_ylabel(u'Specific Humidity (%)')
                        
                        fig.savefig(graphpath + city[ss] + '/' + city[ss]+'_'+METHODNAMES[method-1]+'_Baseline_'+which_current_climate+'_Future_Monthly_Mean_of_Global_Radiation')
                        plt.close()
                    # Plot Solar Terms
                    if vv == 5:
                        fig = plt.figure()
                        fig.suptitle('Future Monthly Total Horizontal Radiation for\n'+ city[ss]+ ' for ' + 
                             MODELNAME[mm] + ' from ' + str(future_years[0]) +' - ' + 
                             str(future_years[1]) + '\n running scenario: ' + SCENNAME[sc])
                        ax = fig.add_subplot(111)
                        sns.set_style("ticks")                     
                        sns.stripplot(x='month', y='rsds_c', data=dfplt, jitter=0.25, hue = 'avg_month', alpha=.60, marker = 's', edgecolor="gray"  )
                        sns.stripplot(x='month', y='rsds', data=dft[dft['station']==ss], jitter=False,  alpha=.90,color = 'black',marker="D")
                        ax.grid()    
                        ax.set_xlabel(u'Month')
                        ax.set_ylabel(u'Total Horizontal Radiation (Wh/m2)')
           
                        
                        fig.savefig(graphpath + city[ss] + '/' + city[ss]+'_'+METHODNAMES[method-1]+'_Baseline_'+which_current_climate+'_Future_Monthly_Mean_of_Specific_Humidity')
                        plt.close()  


#=========================================================
#               Adjust Variables
#=========================================================
# What to use as the current climate monthly means to adjust to FMY with
if which_current_climate == 'tmy' or which_current_climate == 'TMY': 
    dfc0 = dft; #Already aggregated by month, station, model and scenario.
elif which_current_climate == 'gcm' or which_current_climate == 'GCM': 
    # Need to get monthly means for each station, model and scenario.
    # Start with the mean variables
    dfc0 = df.groupby(['station', 'month', 'model', 'scenario'])['tasmean'].agg(['mean']).reset_index()
    dfc0 = dfc0.rename(index=str, columns={"mean": 'tasmean'}) #Rename column from mean to variable name
    
    dfc0['rhsmean'] = list( df.groupby(['station', 'month', 'model', 'scenario'])['rhsmean'].agg(['mean']).reset_index()['mean'] )
    
    # Do rest of variables
    for vv in var:             
        # Grab the monthly means of a variable and put it into dataframe 
        dfc0[VARNAME[vv]] = list( df.groupby(['station', 'month', 'model', 'scenario'])[VARNAME[vv]].agg(['mean']).reset_index()['mean'] )
else: 
    print('\nWARNING: which_current_climate is ,' + which_current_climate + ' should be set to ''tmy'' or ''gcm''. Defaulting to tmy.\n')
    dfc0 = dft;

df_fmy  = pd.DataFrame([]);   
fmy = [0] * len(station_inds);  

# CONVERT TO FMY!

for mm in models:    
    for ss in stations:
        for sc in scen:

            # Initialize new data frame for hourly
            tempdf = pd.DataFrame({'station'    : [city[ss]] * HRS_IN_YEAR, 
                                   'model'      : [MODELNAME[mm]] * HRS_IN_YEAR,
                                   'scenario'   : [SCENNAME[sc]] * HRS_IN_YEAR,
                                   'hoy'        : tmy[stations.index(ss)].hoy,
                                   'doy'        : tmy[stations.index(ss)].doy,
                                   'month'      : tmy[stations.index(ss)].month});
    
            #Add the new adjust weather variables data from tmy to fmy
            tempdf = pd.concat([tempdf, 
                                adjust_to_fmy(var, 
                                              method, 
                                              tmy[stations.index(ss)], 
                                              dfc0[dfc0.station == ss].reset_index(drop = True), 
                                              dfam[(dfam.station == ss) & (dfam.model == mm) & (dfam.scenario == sc)].reset_index(drop = True), 
                                              df1[(df1.station == ss) & (df1.model == mm) & (df1.scenario == sc)].reset_index(drop = True),
                                              hourly_plots, 
                                              future_years, 
                                              city[ss],
                                              mm, sc,
                                              which_current_climate,
                                              graphpath, 
                                              suppress_all_plots) ],
                                axis = 1);
            
            #Concat the new dataframe onto the bottom.
            df_fmy = pd.concat([df_fmy, tempdf], axis = 0);
            
            if 'csv' in outformats:
                tempdf.to_csv(outputpath + 'FMY_' + city[ss] + '_' + MODELNAME[mm] +
                          '_' +  SCENNAME[sc] + '_' + METHODNAMES[method-1] +
                          '_' + which_current_climate + '.csv')

            if 'tmy2' in outformats or 'tmy3' in outformats:
                #Write to fmy
                fmy = weather(weatherpath, city[stations[ss]]);
                fmy.get_weather();
                fmy.tdry    = tempdf.T_fmy;
                fmy.tdew    = tempdf.Tdew_fmy;
                fmy.rhs     = tempdf.RHS_fmy;
                fmy.tothor  = tempdf.Rg_fmy;
                fmy.dirnorm = tempdf.Rdir_fmy;
                fmy.difhor  = tempdf.Rdiff_fmy;
                
                if 'tmy2' in outformats:
                    fmy.write_tmy2( outputpath, MODELNAME[mm], SCENNAME[sc] );
                if 'tmy3' in outformats:
                    fmy.write_tmy3( outputpath, MODELNAME[mm], SCENNAME[sc] );
