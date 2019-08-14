# -*- coding: utf-8 -*-
"""
Created on Mon Jul 29 15:55:39 2019

@author: paul
"""
import pandas as pd;
import numpy as np;
import matplotlib.pyplot as plt;
import seaborn as sns;

# This is just temperature by hour. So can look at monthly averages, check cdfs 


# Read csv for cities.
obs = pd.read_csv('R:/NPCC/2019_Future_Meteorological_Years (19-028-BL)/Data/temperature_4_cities_1948_2017_1.csv');


#tmy3_years   = [ 1990, 2005 ]; # Years that the tmy3 weather files are taken from.
obs = obs[ obs.Year >= tmy3_years[0]];
obs = obs[ obs.Year <= tmy3_years[1]];

obs = obs.drop(columns = ['Boise','Spokane','Portland'])

#To Celcius
obs.Seattle = (obs.Seattle-32)*5/9;

###################

tempT = obs.groupby(['Year','Month','Day']).Seattle.agg(['max', 'min', 'mean'])
tempT = tempT.groupby(level=[0,1]).mean() # Returns the monthly averages
#tempT = tempT.groupby(level=[0]).agg(['mean','std']) # Returns the monthly averages

tempT = tempT.reset_index()

###################
var         = [0, 1]

df = pd.DataFrame()
df = get_data(df, [0], var, models, LAT_TARGETS, LON_TARGETS, stations, daily, interpolate_to_station);

df = df[ df.year >= tmy3_years[0]];
df = df[ df.year <= tmy3_years[1]];

dft = tmy[ss].write_to_daily_df(ss, daily_max)
###################
col = ['#ffbb78', '#aec7e8', '#ff7f0e', '#c5b0d5'];
fig = plt.figure()
ax = fig.add_subplot(111) 
ax.set_xlabel(u'Max Temperatures (C)')
ax.set_ylabel(u'CDF')
ax.hist(tempT['max'], bins = 200, color = col[2], density = True, cumulative = True, histtype = 'step', label = 'Observations' )
ax.hist(dft.tasmax, bins = 200, color = col[0], density = True, cumulative = True, histtype = 'step', label = 'TMY' )
ax.hist((df.tasmax), bins = 200,color = col[1], density = True, cumulative = True,histtype = 'step', label = 'Model Historical' )
ax.legend(loc='upper left');
plt.show()

fig = plt.figure()
ax = fig.add_subplot(111) 
ax.set_xlabel(u'Min Temperatures (C)')
ax.set_ylabel(u'CDF')
ax.hist(tempT['min'], bins = 200, color = col[2], density = True, cumulative = True, histtype = 'step', label = 'Observations' )
ax.hist(dft.tasmin,bins = 200, color = col[0], density = True, cumulative = True, histtype = 'step', label = 'TMY' )
ax.hist((df.tasmin), bins = 200,color = col[1], density = True, cumulative = True,histtype = 'step', label = 'Model Historical' )
ax.legend(loc='upper left');
plt.show()


################
fig = plt.figure()
fig.suptitle('Historical Monthly Max and Min Temperatures for \n'+ city[ss]+ ' for ' + 
             MODELNAME[mm] + ' from ' + str(tmy3_years[0]) +' - ' + 
             str(tmy3_years[1]) + '\n running scenario: ' + SCENNAME[sc])
ax = fig.add_subplot(111)
sns.set_style("ticks")                     
sns.stripplot(x='month', y='tasmax_c', data=df[(df['station']==ss) & (df['model']==mm)], jitter=0.25, color = 'royalblue', alpha=.50, marker = 's', edgecolor="gray"  )
sns.stripplot(x='month', y='tasmin_c', data=df[(df['station']==ss) & (df['model']==mm)], jitter=0.25, color = 'royalblue',alpha=.50, marker = 'o')
sns.stripplot(x='month', y='tasmax', data=dft[dft['station']==ss], jitter=False,alpha=.90, color = 'black',marker="D")
sns.stripplot(x='month', y='tasmin', data=dft[dft['station']==ss], jitter=False, alpha=.90,color = 'black',marker="D")

sns.stripplot(x='Month', y='max', data=tempT, jitter=0.25, alpha=.50, color = 'red', marker="D")
sns.stripplot(x='Month', y='min', data=tempT, jitter=0.25, alpha=.50, color = 'red', marker="D")
ax.set_ylim(-5,35)
ax.grid()
plt.show()




