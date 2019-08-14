"""
#Author:	Paul Kintner
#Updated: 	7/01/2019
#Description: 	This script uses OPeNDAP to download the specified subset of 
the MACAv2-METDATA data then applies multiple methods for temporally downscalling 
data to hourly features.

MACA advocates for using at least 10 models in analyses.It may be worth the effort 
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
    1. Seattle (WA)
    2. Corvallis (OR)
    3. Boise (ID)
    4. Redmond(OR)
    5. Elko (NV)
    6. Burley (ID)
    7. Soda_Springs (ID)
    8. Havre (MT)
    9. Miles CIty (MT)
    
workflow:
    load tmy2 file
    load models historical and future
    determin future peroid of interest and cull the data
    bias correct historical data to tmy2 file
"""

#=========================================================
import os
os.chdir('R:/NPCC/2019_Future_Meteorological_Years (19-028-BL)/Python')
import numpy as np;
import pandas as pd;
from pandas.plotting import scatter_matrix;
import matplotlib.pyplot as plt;
import seaborn as sns;

from cfg import VARNAME, MODELNAME, SCENNAME, DAYS_IN_MONTH, MONTHS_IN_YEAR;

from OPeNDAP_load import get_data;
from variable_conversions import get_monthly_mmms, months2hrs;
#import Belcher; 
from weather import weather; 
from bias_correction import bias_correct;


import metpy.calc as mpcalc
from metpy.units import units

#=========================================================

weatherpath = "E:/TMY2DATA/";

write_tmy2 = True;


interpolate_to_station = True

bias_correction_method = 0; #0 no correction / 1 is EDCDFm / 2 should be CDFm / 3 is linear Hawkins


tmy3_years   = [ 1990, 2005 ]; # Years that the tmy3 weather files are taken from.
future_years = [ 2020, 2049 ]; # Future years from 2006 - 2099.


city        = ["WASeattle3","ORCorvallis3","IDBoise3"];
stations    = [ 0 ]

#Affects plots and loads
method      = 3; # 1 is NPCC, 
                 # 2 is Belcher  
                 # 3 is adjusted NPCC with monthly means of daily max and min

if method == 1:
    daily       = 1; # Gets daily data (1) or monthly data (0).
    daily_max   = 0; # For stats (0) gets the max/min for the whole month, 1 gets the monthly mean of the daily max/min
elif method == 2:
    daily       = 1; # Gets daily data (1) or monthly data (0).
    daily_max   = 1; # For stats (0) gets the max/min for the whole month, 1 gets the monthly mean of the daily max/min
elif method == 3:
    daily       = 1; # Gets daily data (1) or monthly data (0).
    daily_max   = 1; # For stats (0) gets the max/min for the whole month, 1 gets the monthly mean of the daily max/min
    
#--------------------------------------------------------
#   MODEL/VAR/SCEN CHOICES
#--------------------------------------------------------

# bcc-csm1-1 (0) did not have daily data available for 12/31/2099 for the RCP8.5 scenario only.
# CCSM4 (4) and NorESM1-M (19) did not have relative humidity available at daily timescales.
models      = [1, 2, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18] # All complete sets
models      = [3, 4, 5, 6, 8, 9, 10, 11, 12, 15] #The RMJOC-II “10”. Note issues with model 4 in scenario=1
models      = [3] # for testing

var         = [0, 1, 2, 3, 5, 8]
var         = [0, 1, 2, 3]
#var         = [ 2, 3, 8 ]
scen        = [ 1 ]

hourly_plots = [8];

#--------------------------------------------------------
#   Initialize other things.
#--------------------------------------------------------

#lat/lon pairs for point locations
LAT_TARGETS=[0] * (max(stations)+1) #wd.lat
LON_TARGETS=[0] * (max(stations)+1) # wd.lon

# Data frame for tmy data.
dft = pd.DataFrame();
# Data frame for historic model data.
df = pd.DataFrame();
# Data frame for predicted model data.
df1 = pd.DataFrame();

#=========================================================
#                    Load Data!
#=========================================================

# Get tmy weather file
tmy = [0]*len(stations);
for ss in stations:
    tmy[ss] = weather(weatherpath, city[ss]);
    tmy[ss].get_weather();
    dft = pd.concat( [ dft, tmy[ss].write_to_monthly_df(ss, daily_max) ],  axis=0); #Stack stations on top of each other.
    
    LAT_TARGETS[ss] = tmy[ss].lat;
    LON_TARGETS[ss] = tmy[ss].lon;

# Check that the mean daily temperature makes sense for the tmy
#dftd = dft.groupby(['station','doy']).tasmean.agg(['max', 'min','mean'])
#dftd['tmean2'] = (dftd['max'] + dftd['min'])/2
#plt.figure();
#dftd.plot();

# Get historical data, scenario = 0
df = get_data(df, [0], var, models, LAT_TARGETS, LON_TARGETS, stations, daily, interpolate_to_station);

# and get future data for a scenario specified above
df1 = get_data(df1, scen,var,models,LAT_TARGETS, LON_TARGETS, stations, daily, interpolate_to_station);

#scatter_matrix(df)
#plt.show()

#scatter_matrix(df1)
#plt.show()

#plotmonthlytemps.py

#=========================================================
#               Cut data to study "area"
#=========================================================
"""
Cut to look at specific years for models
"""
#tmy3_years   = [ 1990, 2005 ]; # Years that the tmy3 weather files are taken from.
df = df[ df.year >= tmy3_years[0]];
df = df[ df.year <= tmy3_years[1]];

#future_years = [ 2060, 2090 ]; # Future years from 2006 - 2099.
df1 = df1[ df1.year >= future_years[0]];
df1 = df1[ df1.year <= future_years[1]];


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
        dft_temp = tmy[ss].write_to_daily_df(ss, daily_max);
        
        for vv in var:
           # if vv == 0 or vv == 1:
                tempdf1[str(VARNAME[vv]+'_c')] = bias_correct( 
                        df[(df['station']==ss)][VARNAME[vv]],
                        df1[(df1['station']==ss)][VARNAME[vv]], 
                        dft_temp[(dft_temp['station']==ss)][VARNAME[vv]], 
                        vv, bias_correction_method)
           # if vv == 2 or vv == 3:
           #     df1[str(VARNAME[vv]+'_c')] = df1[VARNAME[vv]];

    df1 = pd.concat([df1, tempdf1], axis=1);
        
elif bias_correction_method == 0: #The corrected variables are just the variables
        for vv in var:
            df1[str(VARNAME[vv]+'_c')] = df1[VARNAME[vv]];

else: print("Invalid bias_correction_method chose can be 0,1,2,3")
    
# Give the historical data frame that same name for their values
for vv in var:
    df[str(VARNAME[vv]+'_c')] = df[VARNAME[vv]];
            
 
#=========================================================
#                Check aggregation Method
#=========================================================
df      = get_monthly_mmms(df, var, daily_max);
df1     = get_monthly_mmms(df1, var, daily_max);

#=========================================================
#                Mean Temperature
#=========================================================
if 0 in var or 1 in var:
    if 0 in var and 1 in var:
        df['tasmean'] = ( df['tasmax'] + df['tasmin'] ) / 2.;

        df1['tasmean'] = ( df1['tasmax_c'] + df1['tasmin_c'] ) / 2.;
    else:
        print("\nWARNING: Selecting the max/min temperature without the other add 0/1 to var")  
        
#=========================================================
#                Mean Relative Humidity
#=========================================================
if 2 in var or 3 in var:
    if 2 in var and 3 in var:
        df['rhsmean'] = ( df['rhsmax'] + df['rhsmin'] ) / 2.;

        df1['rhsmean'] = ( df1['rhsmax_c'] + df1['rhsmin_c'] ) / 2.;
    else:
        print("\nWARNING: Selecting the max/min relative humidity without the other add 2/3 to var")    

############################################################################
#           Historical Plot
############################################################################
for mm in models:
    for ss in stations:
        for sc in scen:
            for vv in var:
                # Plot temperatures
                if vv == 1:
                    fig = plt.figure()
                    fig.suptitle('Historical Monthly Max and Min Temperatures for \n'+ city[ss]+ ' for ' + 
                         MODELNAME[mm] + ' from ' + str(tmy3_years[0]) +' - ' + 
                         str(tmy3_years[1]) + '\n running scenario: ' + SCENNAME[sc])
                    ax = fig.add_subplot(111)
                    sns.set_style("ticks")                     
                    sns.stripplot(x='month', y='tasmax_c', data=df[(df['station']==ss) & (df['model']==mm)], jitter=0.25, color = 'royalblue', alpha=.60, marker = 's', edgecolor="gray"  )
                    sns.stripplot(x='month', y='tasmin_c', data=df[(df['station']==ss) & (df['model']==mm)], jitter=0.25, color = 'royalblue',alpha=.50, marker = 'o')
                    sns.stripplot(x='month', y='tasmax', data=dft[dft['station']==ss], jitter=False,alpha=.90, color = 'black',marker="D")
                    sns.stripplot(x='month', y='tasmin', data=dft[dft['station']==ss], jitter=False, alpha=.90,color = 'black',marker="D")
                    
                    ax.set_ylim(-5,35)
                    ax.grid()
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
                
                # Plot specific humidity
                if vv == 8:
                    fig = plt.figure()
                    fig.suptitle('Historical Monthly Specific Humidity for\n'+ city[ss]+ ' for ' + 
                         MODELNAME[mm] + ' from ' + str(tmy3_years[0]) +' - ' + 
                         str(tmy3_years[1]) + '\n running scenario: ' + SCENNAME[sc])
                    ax = fig.add_subplot(111)
                    sns.set_style("ticks")                     
                    sns.stripplot(x='month', y='huss', data=df[(df['station']==ss) & (df['model']==mm)], jitter=0.25, color = 'royalblue', alpha=.60, marker = 's', edgecolor="gray"  )
                    sns.stripplot(x='month', y='huss', data=dft[dft['station']==ss], jitter=False,  alpha=.90,color = 'black',marker="D")
                    ax.grid()    

#=========================================================
#                Find Month Closest to Average
#=========================================================
                    
# For NPCC method first find the most average month and then find the constants for transformation.
# convert temperatures how they describe with the linear transformations for temperature and Relative Humidity.
tempdf = df1;
tempdf['err'] = 0; #Initialize "error" column 

for vv in var:             
    # Grab the monthly means of a variable and store in repeating dataframe
    monthly_mmm = df1.groupby(['station', 'month', 'model', 'scenario'])[VARNAME[vv]+'_c'].agg(['max', 'min', 'mean'])
    monthly_mmm = pd.concat([monthly_mmm]*int(len(df1)/MONTHS_IN_YEAR))
    monthly_mmm = monthly_mmm.reset_index();
    
    #Check the normalized vector distance for everymonth
    tempdf['err'] = tempdf['err'] + ( (tempdf[VARNAME[vv]+'_c'] - monthly_mmm['mean']) / monthly_mmm['mean'] ) **2;
         # (monthly_mmm['max'] - monthly_mmm['min']) ) ** 2
        
tempdf['err'] = np.sqrt(tempdf['err']);              
    
df1['avg_month'] = tempdf.groupby(['station', 'month', 'model', 'scenario'])['err'].transform(min) == tempdf['err']

#df of average months
dfam = (df1[ df1['avg_month'] == True ].sort_values(by='month')).reset_index();

############################################################################
#               Future plots
############################################################################
for mm in models:
    for ss in stations:
        for sc in scen:
            for vv in var:
                if vv == 1:
                    fig = plt.figure()
                    fig.suptitle('Future Monthly Max and Min Temperatures for \n'+ city[ss]+ ' for ' + 
                                 MODELNAME[mm] + ' from ' + str(future_years[0]) +' - ' + 
                                 str(future_years[1]) + '\n running scenario: ' + SCENNAME[sc])
                    ax = fig.add_subplot(111)
                    sns.set_style("ticks")
                    sns.stripplot(x='month', y='tasmax_c', data=df1[(df1['station']==ss) & (df1['model']==mm)], jitter=.25, hue = 'avg_month', alpha=.60, marker = 's', edgecolor="gray"  )
                    sns.stripplot(x='month', y='tasmin_c', data=df1[(df1['station']==ss) & (df1['model']==mm)], jitter=.25, hue = 'avg_month', alpha=.60, marker = 'o')
                    sns.stripplot(x='month', y='tasmax', data=dft[(dft['station']==ss)], jitter=False,  alpha=.90, color = 'black',marker="D")
                    sns.stripplot(x='month', y='tasmin', data=dft[(dft['station']==ss)], jitter=False, alpha=.90,color = 'black',marker="D")
                    ax.set_ylim(-5,35)
                    ax.grid()
                if vv == 3:
                    fig = plt.figure()
                    ax = fig.add_subplot(111)
                    fig.suptitle('Future Monthly Max and Min Relative Humidity for\n'+ city[ss]+ ' for ' + 
                                 MODELNAME[mm] + ' from ' + str(future_years[0]) +' - ' + 
                                 str(future_years[1]) + '\n running scenario: ' + SCENNAME[sc])
                    sns.set(style="ticks")
            
                    sns.stripplot(x='month', y='rhsmax_c', data=df1[(df1['station']==ss) & (df1['model']==mm)],jitter=.25, hue = 'avg_month', alpha=.60, marker = 's', edgecolor="gray"  )
                    sns.stripplot(x='month', y='rhsmin_c', data=df1[(df1['station']==ss) & (df1['model']==mm)], jitter=.25, hue = 'avg_month', alpha=.60, marker = 'o')
                    sns.stripplot(x='month', y='rhsmax', data=dft[(dft['station']==ss)], jitter=False,  alpha=.90, color = 'black', marker="D")
                    sns.stripplot(x='month', y='rhsmin', data=dft[(dft['station']==ss)], jitter=False,  alpha=.90, color = 'black', marker="D")
                    ax.set_ylim(0,110)
                    ax.grid() 
                if vv == 8:
                    fig = plt.figure()
                    ax = fig.add_subplot(111)
                    fig.suptitle('Future Monthly Specific Humidity for\n'+ city[ss]+ ' for ' + 
                                 MODELNAME[mm] + ' from ' + str(future_years[0]) +' - ' + 
                                 str(future_years[1]) + '\n running scenario: ' + SCENNAME[sc])
                    sns.set(style="ticks")
            
                    sns.stripplot(x='month', y='huss', data=df1[(df1['station']==ss) & (df1['model']==mm)], jitter=.25, hue = 'avg_month', alpha=.60, marker = 'o')
                    sns.stripplot(x='month', y='huss', data=dft[(dft['station']==ss)], jitter=False,  alpha=.90, color = 'black', marker="D")
                    ax.grid()
       
os.system("pause")
 
#=========================================================
#                Initialize FMY
#=========================================================
tfmy        = pd.DataFrame({'T_tmy':tmy[0].tdry, 'RHS_tmy':100*tmy[0].rhs});
T_fmy       = np.array([]);
RHS_fmy     = np.array([]);
Tdew_fmy    = np.array([]);

#=========================================================
#                Find A, B, C, D for each month and plot!
#=========================================================                   
if method == 1 or method == 3:
    # Monthly Temperature Constants for Transformation
    if 0 in var and 1 in var:
        dfam['A_tas'] = ( dfam.tasmax_c - dfam.tasmin_c ) / ( dft.tasmax - dft.tasmin)               
        dfam['B_tas'] = dfam.tasmax_c - dfam.A_tas * dft.tasmax;
    # Monthly Relative Humidity Constants for Transformation
    if 2 in var and 3 in var:
        dfam['C_tas'] = ( dfam.rhsmax_c - dfam.rhsmin_c ) / ( dft.rhsmax - dft.rhsmin)               
        dfam['D_tas'] = dfam.rhsmax_c/100 - dfam.C_tas * dft.rhsmax/100; #Divide by 100 because to go from %
    
    for mm in models:
        for ss in stations:
            for sc in scen:
                for ii in range(1,MONTHS_IN_YEAR+1):
                    
                    if 0 in var and 1 in var:
                        T_fmy = np.append(T_fmy, dfam.A_tas[ii-1]*tmy[ss].tdry[tmy[ss].month==ii] + dfam.B_tas[ii-1])
                        
                        if 0 in var and 1 in hourly_plots:
                        #Temperature Plots
                            fig = plt.figure(figsize=(11, 7))
                            fig.suptitle('Transformed Hourly TMY Temperatures for Month ' + str(ii) + '\n'+ city[ss]+ ' for ' + 
                                         MODELNAME[mm] + ' from ' + str(future_years[0]) +' - ' + 
                                         str(future_years[1]) + '\n running scenario: ' + SCENNAME[sc])
                        
                            ax = fig.add_subplot(111)
                            ax.set_xlabel(u'Hour of Year')
                            ax.set_ylabel(u'Temperature (C)')
                            ax.ticklabel_format(style='plain')
                            
                            ax.plot(tmy[ss].hoy[tmy[ss].month==ii], tmy[ss].tdry[tmy[ss].month==ii])
                            ax.plot(tmy[ss].hoy[tmy[ss].month==ii], dfam.A_tas[ii-1]*tmy[ss].tdry[tmy[ss].month==ii] + dfam.B_tas[ii-1])
                            
                            ax.plot(tmy[ss].hoy[tmy[ss].month==ii], dfam.A_tas[ii-1]*tmy[ss].tdry[tmy[ss].month==ii] + dfam.B_tas[ii-1] - tmy[ss].tdry[tmy[ss].month==ii], '--' )
                            
                            ax.grid()
                            plt.show()
                            
                            if method == 1:
                                fig.savefig('Graphs/'+city[ss]+'_NPCC_Transformed_Hourly_TMY_Temperatures_for_Month_'+ str(ii))
                            elif method == 3:
                                fig.savefig('Graphs/'+city[ss]+'_Adjusted_NPCC_Transformed_Hourly_TMY_Temperatures_for_Month_'+ str(ii))

                    if 2 in var and 3 in var:
                        
                        RHS_fmy = np.append(RHS_fmy, 100*(dfam.C_tas[ii-1]*tmy[ss].rhs[tmy[ss].month==ii] + dfam.D_tas[ii-1]))
                        
                        # Force high values down.
                        RHS_fmy[RHS_fmy>100] = 100;

                        if 2 in var and 3 in hourly_plots:

                                #Relative Humidity Plots
                                fig = plt.figure(figsize=(11, 7))
                                fig.suptitle('Transformed Hourly TMY Relative Humidity for Month ' + str(ii) + '\n'+ city[ss]+ ' for ' + 
                                             MODELNAME[mm] + ' from ' + str(future_years[0]) +' - ' + 
                                             str(future_years[1]) + '\n running scenario: ' + SCENNAME[sc])
                            
                                ax = fig.add_subplot(111)
                                ax.set_xlabel(u'Hour of Year')
                                ax.set_ylabel(u'Relative Humidity (%)')
                                ax.ticklabel_format(style='plain')
                                
                                ax.plot(tmy[ss].hoy[tmy[ss].month==ii], 100*tmy[ss].rhs[tmy[ss].month==ii])
                                ax.plot(tmy[ss].hoy[tmy[ss].month==ii], 100*(dfam.C_tas[ii-1]*tmy[ss].rhs[tmy[ss].month==ii] + dfam.D_tas[ii-1]))
                                
                                ax.plot(tmy[ss].hoy[tmy[ss].month==ii], 100*(dfam.C_tas[ii-1]*tmy[ss].rhs[tmy[ss].month==ii] + dfam.D_tas[ii-1]) - 100*tmy[ss].rhs[tmy[ss].month==ii], '--' )
                                
                                ax.grid()
                                plt.show()
                                
                                if method == 1:
                                    fig.savefig('Graphs/'+city[ss]+'_NPCC_Transformed_Hourly_TMY_Relative_Humidity_for_Month_'+ str(ii))
                                elif method == 3:
                                    fig.savefig('Graphs/'+city[ss]+'_Adjusted_NPCC_Transformed_Hourly_TMY_Relative_Humidity_for_Month_'+ str(ii))

#=========================================================
#               Adjust Belecher Variables
#========================================================= 
elif method == 2:
    monthlyf1 = df1.groupby(['station', 'month', 'scenario']).mean().reset_index()
    monthlytmy = dft.groupby(['station', 'month']).mean().reset_index()
    
    for ss in stations:
        for sc in scen:
            
            temp1 = monthlyf1[ ((monthlyf1.station == ss) | (monthlyf1.scenario == sc))] # Get relavent MACA data
            tempt = monthlytmy[ (monthlytmy.station == ss)] #Get relevant station data
    
            # Adjust Temperatures
            if 0 in var and 1 in var:
                #Scaling
                adbt =  months2hrs(((temp1.tasmax_c - tempt.tasmax) - (temp1.tasmin_c - tempt.tasmin)) / ( tempt.tasmax - tempt.tasmin));
                #fmy dry bulb ttem
                T_fmy = tmy[ss].tdry + months2hrs((temp1.tasmean - tempt.tasmean)) + adbt*(tmy[ss].tdry - months2hrs(tempt.tasmean))
                                     
                if 0 in var and 1 in hourly_plots:
    
                    # Plot fmy and tmy temperature
                    for ii in range(1,MONTHS_IN_YEAR+1):
        
                        fig = plt.figure(figsize=(11, 7))
                        fig.suptitle('Belcher Hourly TMY Temperatures for Month ' + str(ii) + '\n'+ city[ss]+ ' for ' + 
                                     MODELNAME[mm] + ' from ' + str(future_years[0]) +' - ' + 
                                     str(future_years[1]) + '\n running scenario: ' + SCENNAME[sc])
                    
                        ax = fig.add_subplot(111)
                        ax.set_xlabel(u'Hour of Year')
                        ax.set_ylabel(u'Temperature (C)')
                        ax.ticklabel_format(style='plain')
                        
                        ax.plot(tmy[ss].hoy[tmy[ss].month==ii], tmy[ss].tdry[tmy[ss].month==ii])
                        
                        ax.plot(tmy[ss].hoy[tmy[ss].month==ii], T_fmy[tmy[ss].month==ii])
                        ax.plot(tmy[ss].hoy[tmy[ss].month==ii], T_fmy[tmy[ss].month==ii] -  tmy[ss].tdry[tmy[ss].month==ii], '--' )
                        ax.grid()
                        plt.show()
                        
                        fig.savefig('Graphs/'+city[ss]+'_Belcher_Transformed_Hourly_TMY_Temperatures_for_Month_'+ str(ii))


            # Adjust Relative Humidity
            if 2 in var and 3 in var:
                #Scaling
                arhs = months2hrs( ((temp1.rhsmax_c - tempt.rhsmax) - (temp1.rhsmin_c - tempt.rhsmin) )/ ( tempt.rhsmax - tempt.rhsmin) );
                arhs = months2hrs( 1 + ((temp1.rhsmean) - (tempt.rhsmean) ) / 100. );
                arhs = months2hrs( 1 + ((temp1.rhsmean) - (tempt.rhsmean) ) / tempt.rhsmean )
                
                
                #fmy relative humdity
                RHS_fmy = 100.0 * arhs * tmy[ss].rhs;
                
                #Force high values down
                RHS_fmy[RHS_fmy>100] = 100;
                
                if 2 in hourly_plots or 3 in hourly_plots:
                    #Plot fmy and tmy RHS
                    for ii in range(1,MONTHS_IN_YEAR+1):
        
                        fig = plt.figure(figsize=(11, 7))
                        fig.suptitle('Belcher Hourly TMY Relative Humidity for Month ' + str(ii) + '\n'+ city[ss]+ ' for ' + 
                                     MODELNAME[mm] + ' from ' + str(future_years[0]) +' - ' + 
                                     str(future_years[1]) + '\n running scenario: ' + SCENNAME[sc])
                    
                        ax = fig.add_subplot(111)
                        ax.set_xlabel(u'Hour of Year')
                        ax.set_ylabel(u'Temperature (C)')
                        ax.ticklabel_format(style='plain')
                        
                        ax.plot(tmy[ss].hoy[tmy[ss].month==ii], 100*tmy[ss].rhs[tmy[ss].month==ii])
                        
                        ax.plot(tmy[ss].hoy[tmy[ss].month==ii], RHS_fmy[tmy[ss].month==ii])
                        ax.plot(tmy[ss].hoy[tmy[ss].month==ii], RHS_fmy[tmy[ss].month==ii] -  100*tmy[ss].rhs[tmy[ss].month==ii], '--' )
                        ax.grid()
                        plt.show()
                        
                        fig.savefig('Graphs/'+city[ss]+'_Belcher_Transformed_Hourly_TMY_Temperatures_for_Month_'+ str(ii))
                

#=========================================================
#       Adjust Dew Point or Dew Temperature
#=========================================================
if 8 in var:
    if 0 not in var:  print('\nWARNING: Can not find fmy dew point without drybulb, include 0 and 1 in var selection.\n')

    monthlyf1 = df1.groupby(['station', 'month', 'scenario']).mean().reset_index()
    monthlytmy = dft.groupby(['station', 'month']).mean().reset_index()
    
    for ss in stations:
        for sc in scen:
            
            temp1 = monthlyf1[ ((monthlyf1.station == ss) | (monthlyf1.scenario == sc))] # Get relavent MACA data
            tempt = monthlytmy[ (monthlytmy.station == ss)] #Get relevant station data
            
            ahuss = months2hrs( 1 + (temp1.huss - tempt.huss) )
                            
            huss_fmy = ahuss * tmy[ss].huss;
            Tdew_fmy = mpcalc.dewpoint_from_specific_humidity(huss_fmy, T_fmy * units.degC, tmy[ss].press * units.mbar).magnitude
            
            #Force high values down
            Tdew_fmy[Tdew_fmy > T_fmy] = T_fmy[Tdew_fmy > T_fmy];
            
            if 8 in hourly_plots:
    
                #Plot fmy and tmy tdew
                for ii in range(1,MONTHS_IN_YEAR+1):
            
                    fig = plt.figure(figsize=(11, 7))
                    fig.suptitle('Belcher Hourly TMY Dew Point for Month ' + str(ii) + '\n'+ city[ss]+ ' for ' + 
                                 MODELNAME[mm] + ' from ' + str(future_years[0]) +' - ' + 
                                 str(future_years[1]) + '\n running scenario: ' + SCENNAME[sc])
                
                    ax = fig.add_subplot(111)
                    ax.set_xlabel(u'Hour of Year')
                    ax.set_ylabel(u'Dew Point (C)')
                    ax.ticklabel_format(style='plain')
                    
                    ax.plot(tmy[ss].hoy[tmy[ss].month==ii], tmy[ss].tdew[tmy[ss].month==ii])
                    
                    ax.plot(tmy[ss].hoy[tmy[ss].month==ii], Tdew_fmy[tmy[ss].month==ii])
                    ax.plot(tmy[ss].hoy[tmy[ss].month==ii], Tdew_fmy[tmy[ss].month==ii] - tmy[ss].tdew[tmy[ss].month==ii], '--' )
                    ax.grid()
                    plt.show()
                    
                    fig.savefig('Graphs/'+city[ss]+'_Belcher_Transformed_Hourly_TMY_DewPoint_for_Month_'+ str(ii))

else:
    Tdew_fmy = tmy[ss].tdew;
                    
#=========================================================
#               Check and Finilize FMY
#=========================================================
    
if sum(RHS_fmy < 0) > 0 or sum(RHS_fmy > 100) > 0:
    print('\nWARNING: The FMY of Relative Humidity is less than 0% ('+str(sum(RHS_fmy < 0))+ 'values) and/or greater than 100% ('+str(sum(RHS_fmy > 100))+ 'values)\n')
                
if sum(Tdew_fmy > T_fmy) > 0:  
    print('\nWARNING: The FMY of Dew Point is greater than the Dry Bulb ('+str(sum(Tdew_fmy > T_fmy))+ 'values)\n')

if 0 in var and 1 in var: tfmy['T_fmy'] = T_fmy;
if 2 in var and 3 in var: tfmy['RHS_fmy'] = RHS_fmy;

#
## Create an instance of the PairGrid class.
#grid = sns.PairGrid(data= tfmy)
##grid = grid.map_upper(plt.scatter, color = 'darkred')
#grid = grid.map_upper(sns.regplot,  line_kws={'color': 'red'}, scatter_kws={'s':2})
#grid = grid.map_upper(reg_stats)
#
## Map a histogram to the diagonal
#grid = grid.map_diag(plt.hist, bins = 20,  edgecolor = 'k')
## Map a density plot to the lower triangle
#grid = grid.map_lower(sns.kdeplot, cmap = "Blues_d")
#grid = grid.map_lower(corr)
#
#
## Function to calculate correlation coefficient between two arrays
#def corr(x, y, **kwargs):
#    
#    # Calculate the value
#    coef = np.corrcoef(x, y)[0][1]
#    # Make the label
#    label = r'$\rho$ = ' + str(round(coef, 2))
#    
#    # Add the label to the plot
#    ax = plt.gca()
#    ax.annotate(label, xy = (0.2, 0.95), size = 20, xycoords = ax.transAxes)
#    
## Function to calculate linear regression between variables
#def reg_stats(x, y, **kwargs):
#    from scipy import stats
#
#    # Calculate the values
#    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
#    
#    # Make the label
#    label = 'slope: ' + str(round(slope, 2)) + '\nintercept: '+ str(round(intercept, 2)) + '\nR^2: '+str(round(r_value**2, 2)) 
#    
#    # Add the label to the plot
#    ax = plt.gca()
#    ax.annotate(label, xy = (0.2, 0.85), size = 12, xycoords = ax.transAxes)
#    