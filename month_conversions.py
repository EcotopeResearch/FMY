# -*- coding: utf-8 -*-
"""
Monthly variable conversions!

"""
import numpy as np;
import pandas as pd;
from cfg import VARNAME;

###############################################################################   

# Group into monthly means for the different methods, 
# Either daily = 1 and we get the monthly absoulte max and min or daily = 0 and we get 
# the monthly average of the maxs and mins
def get_monthly_mmms(df1, var, daily):
    
    holddf1 = df1;
    dfT = None; dfR = None; tempH = None; tempRg = None;
    
    #If looking at absolute max/min from month
    if daily == 0:
        if 0 in var and 1 in var:
            # Temperature
            tempT1 = holddf1.groupby(['station','model','scenario','year','month']).tasmax_c.agg(['max'])
            tempT2 = holddf1.groupby(['station','model','scenario','year','month']).tasmin_c.agg(['min'])
        if 2 in var and 3 in var:
            # Relative Humidity
            tempR1 = holddf1.groupby(['station','model','scenario','year','month']).rhsmax_c.agg(['max'])
            tempR2 = holddf1.groupby(['station','model','scenario','year','month']).rhsmin_c.agg(['min'])
        
    #If looking at monthly mean of max/min from days
    else:
        if 0 in var and 1 in var:
            # Temperature
            tempT1 = holddf1.groupby(['station','model','scenario','year','month','doy']).tasmax_c.agg(['max']) #Get daily max
            tempT1 = tempT1.groupby(level=[0,1,2,3,4]).mean() # Returns the monthly averages
            tempT2 = holddf1.groupby(['station','model','scenario','year','month','doy']).tasmin_c.agg(['min']) #Get daily min
            tempT2 = tempT2.groupby(level=[0,1,2,3,4]).mean() # Returns the monthly averages
        if 2 in var and 3 in var:
            # Relative Humidity
            tempR1 = holddf1.groupby(['station','model','scenario','year','month','doy']).rhsmax_c.agg(['max']) #Get daily max
            tempR1 = tempR1.groupby(level=[0,1,2,3,4]).mean() # Returns the monthly averages
            tempR2 = holddf1.groupby(['station','model','scenario','year','month','doy']).rhsmin_c.agg(['min']) #Get daily min
            tempR2 = tempR2.groupby(level=[0,1,2,3,4]).mean() # Returns the monthly averages
      
     # Look at Mean Global Radiation because that's all we got!         
    if 5 in var:
        tempH = holddf1.groupby(['station','model','scenario','year','month']).rsds_c.agg(['mean']);
        tempH = tempH.rename(index=str, columns={"mean": "rsds"})
               
    # Look at Mean Specific Humidity, that's all we got!         
    if 8 in var:
        tempRg = holddf1.groupby(['station','model','scenario','year','month']).huss_c.agg(['mean']);
        tempRg = tempRg.rename(index=str, columns={"mean": "huss"})
            
    #Combine the data frames, temperature, relative humditiy
    if 0 in var and 1 in var:
        dfT = pd.concat([tempT1,tempT2], axis = 1)
        dfT = dfT.rename(index=str, columns={"max": "tasmax", "min": "tasmin"})
    if 2 in var and 3 in var:       
        dfR = pd.concat([tempR1, tempR2], axis = 1);
        dfR = dfR.rename(index=str, columns={"max": "rhsmax", "min": "rhsmin"})
       
    df1 = pd.concat([ dfT, dfR, tempH, tempRg ], axis = 1)
  
    # Add important columns back to dataframe and make sure they're integers for some reason have to mape through floats...
    df1 = df1.reset_index();
    df1['station']  = list(map(int,map(float,df1['station'])));
    df1['model']    = list(map(int,map(float,df1['model'])));
    df1['scenario'] = list(map(int,map(float,df1['scenario'])));
    df1['year']     = list(map(int,map(float,df1['year'])));
    df1['month']    = list(map(int,map(float,df1['month'])));

    for vv in var:
            df1[str(VARNAME[vv]+'_c')] = df1[VARNAME[vv]];

    return df1;



###############################################################################   

# Months to Hours!
# Takes a variable with monthly data (i.e. a array of len 12) and converts it 
# to hourly data (i.e. len of 8760) by repeating the constants through a month.
def months2hrs(month_data):
    from cfg import HRS_IN_MONTH, HRS_IN_YEAR, MONTHS_IN_YEAR
    
    temp            = np.array([0.0]*HRS_IN_YEAR)
    hrs_tot         = 0;
    
    for ii in range(0,MONTHS_IN_YEAR):
        
        temp[ hrs_tot : (hrs_tot + HRS_IN_MONTH[ii]) ] = month_data[ii];
        
        hrs_tot                             = hrs_tot + HRS_IN_MONTH[ii]

    return(temp);