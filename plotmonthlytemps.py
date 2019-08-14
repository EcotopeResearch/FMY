# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 16:31:15 2019

@author: paul
"""

#=========================================================
#               Historical Mean temp all years mean
#=========================================================
df['tasmean'] = ( df['tasmax'] + df['tasmin'] ) / 2.;

data =  df.groupby(['month','model']).mean()
days = np.arange(0,12);

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlabel(u'Day of Year')
ax.set_ylabel(u'Mean Monthly Tempreatures (K)')
ax.ticklabel_format(style='plain')
for ii in models:
    data = df[ df['model'] == ii ]
    data = data.groupby(['month']).mean()
    ax.plot(days, data['tasmean'])

plt.show()


#=========================================================
#               Future Mean Temp year mean
#=========================================================
scen        = [ 1 ]

df1 = pd.DataFrame();
df1 = get_data(df1, scen,var,models,LAT_TARGETS, LON_TARGETS, stations) ;

df1['tasmean'] = ( df1['tasmax'] + df1['tasmin'] ) / 2.;

days = np.arange(0,12);

###################
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlabel(u'Day of Year')
ax.set_ylabel(u'Mean Monthly Tempreatures (K)')
ax.ticklabel_format(style='plain')
for ii in models:
    print(ii)
    data1 = df1[ df1['model'] == ii ]
    data1 = data1.groupby(['month']).mean()
    ax.plot(days, data1['tasmean'])
plt.show()
#################

fig = plt.figure()
ax = fig.add_subplot(111) 
ax.set_xlabel(u'Mean Monthly Tempreatures (C)')
ax.set_ylabel(u'PDF')
data = df[ df['model'] == 2 ]
ax.hist(data['tasmean']-273.15, bins = 100, density = True, histtype = 'step' )
#ax.hist(data['tasmax']-273.15, bins = 400, density = True, cumulative = True, histtype = 'step' )
for ii in models:
    data1 = df1[ df1['model'] == ii ]

   # ax.hist(data1['tasmax']-273.15, bins = 400, density = True, cumulative = True, histtype = 'step' )
    ax.hist(data1['tasmax']-273.15, bins = 100, density = True, histtype = 'step' )
plt.show()

#=========================================================
#               Change in Mean temp
#=========================================================
dfmean_models   =  df.groupby(['month']).mean()['tasmean']
dfstd_models    =  df.groupby(['month']).std()['tasmean']

days = np.arange(0,12);

fig, ax = plt.subplots(1)
ax.set_xlabel(u'Day of Year')
ax.set_ylabel(u'Change in Mean Daily Tempreatures (C)')
#ax.plot(days, dfmean_models, lw=2, label='baseline average', color='blue')
#ax.fill_between(days, dfmean_models+dfstd_models, dfmean_models-dfstd_models, facecolor='blue', alpha=0.5)

for ii in models:
    data        =  df[ df['model'] == ii ]
    data1       =  df1[ df1['model'] == ii ]
    
    mu          =  data.groupby(['month']).mean()['tasmean']
    mu1         =  data1.groupby(['month']).mean()['tasmean']
    
    ax.plot(days, mu1-mu, lw=2, label=ii)
   # ax.fill_between(days, mu1+sigma1, mu1-sigma1, alpha=0.5)
mu1         =  df1.groupby(['month']).mean()['tasmean']
ax.plot(days, mu1-mu, lw=2, label="Model Mean")

ax.legend(loc='lower left',ncol=6)
plt.show()
##########################################################
# Calculate statistics of the average changes on the set of 
# models (multi-model mean gives signal, multi-model
# standard deviation gives measure of uncertainty)
