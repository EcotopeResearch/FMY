# -*- coding: utf-8 -*-
"""
Created on Mon Jul  8 10:44:37 2019

Correct biases to station data from gridcell.

The MACA dataset is gridded data (either on 4-km or ~6-km grid cells). The 
incongruence between gridded data and station data can be rather pronounced in
complex terrain, necessitating that the downscaled projections be adjusted for
the station locations. 

This process should be done by first identifying the grid-cell which is 
co-located with each station and then performing a secondary bias correction 
procedure of the grid-cell data to that station data. For variables such as 
radiation and winds where station data is not available, interpolation of the 
4-km value at the center of the cell to the station location is acceptable. 

@author: paul
"""

VARNAME =('tasmax','tasmin','rhsmax','rhsmin','pr','rsds','uas','vas','huss');
VARLONGNAME=('air_temperature','air_temperature', 'relative_humidity','relative_humidity',\
                'precipitation','surface_downwelling_shortwave_flux_in_air',\
                'eastward_wind','northward_wind','specific_humidity')

###############################################################################
def bias_correct( mod_c, mod_p, obs_data, var, method ):
    """ 
    mod_c is the modeled climate current climate data
    mod_p is the modeled projected climate data
    obs_data is the observed data tmy in this case
    var is the variable being looked at
    method is the method to be used should be 1 or 2
    """
    station_data = None; 
    
    
    if var < 5 or var == 8 :     # If radiation or wind just interp to station as MACA says
        
        if method == 1:
             # Bias correct to the station using EDCDFm
             station_data = bias2station_EDCDFm( mod_c, mod_p, obs_data, var);
            
        elif method == 2:
             # Bias correct to the station using CDFm
             station_data =bias2station_CDFm( mod_c, mod_p, obs_data, var );
             
        elif method == 3:
             station_data = bias2station_linear( mod_c, mod_p, obs_data, var);

    #else 
    return(station_data);

    
###############################################################################    
    #        bias_correct( df['tasmean'], df1['tasmean'], dft['tasmean'], vv)


###############################################################################
def bias2station_linear( mod_c, mod_p, obs_data, var):
    # Bias correctoction method following Hawkins et al., 2013, and following 
    # the describtion from http://www.ccafs-climate.org/bias_correction/
    
    import numpy as np;
    import matplotlib.pyplot as plt;
    from cfg import VARNAME;

    data_adj = np.average(obs_data) + np.std(obs_data) / np.std(mod_c) * (mod_p - np.average(mod_c));
    
    
    ###########################################################################
    col = ['#ffbb78', '#aec7e8', '#ff7f0e', '#c5b0d5'];
    fig = plt.figure()
    fig.suptitle(VARNAME[var])
    ax = fig.add_subplot(111) 
    ax.set_xlabel(u'Temperatures (C)')
    ax.set_ylabel(u'CDF')
    ax.hist(obs_data, color = col[0], density = True, histtype = 'step', label = 'tmy' )
    ax.hist(mod_c, color = col[1], density = True, histtype = 'step', label = 'model historical' )
    ax.hist(mod_p, color = col[2], density = True, histtype = 'step', label = 'model predicted' )
    ax.hist(data_adj, color = col[3], density = True, histtype = 'step', label = 'model corrected' )
    ax.legend();
    plt.show()
    ###########################################################################
    
    return(data_adj);

###############################################################################
def bias2station_CDFm( mod_c, data, obs_data, var ):
    import numpy as np;
    from scipy import stats;
    from scipy import interpolate;
    import matplotlib.pyplot as plt;
    from cfg import VARNAME
    
    minx = min( [min(mod_c), min(data), min(obs_data)] ) - abs(min(mod_c));
    maxx = max( [max(mod_c), max(data), max(obs_data)] ) + abs(max(mod_c));
    lnspc = np.linspace(minx, maxx, len(data))
    
     # Fits the pdf of the data using a 2-parameter beta distribution adjusting 
    # the beta distribution with variables for α, β, and loc (lower limit), 
    # scale (upper limit - lower limit). Documentation also states that the four
    # parameter distribution scales the range of the beta distribuition. This
    # is similar to what Li et al. (2010) and Watterson (2008) describe where
    # their four parameters are  α, β, and the minimum and maximum value of x. 
    #For purposes here, these two methods are considered to be equivalent. 
    
    abmc,bbmc,cbmc,dbmc = stats.beta.fit(mod_c);    # Model current climate
    abo,bbo,cbo,dbo     = stats.beta.fit(obs_data); # Observed current climate 

    # The cumulative distrubtion functions! 
    cdf_o = stats.beta.cdf(lnspc, abo,bbo,cbo,dbo);      # Observed current climate
    cdf_m = stats.beta.cdf(lnspc, abmc,bbmc,cbmc,dbmc);      # Observed current climate

    # Create functions to act as the inverse of the cdfs
    Finv_o = interpolate.interp1d( cdf_o, lnspc );
    
    # The data adjustment!
    data_adj = Finv_o(stats.beta.cdf(data, abmc, bbmc, cbmc, dbmc))
    
# =============================================================================
    col = ['#ffbb78', '#aec7e8', '#ff7f0e', '#c5b0d5'];
    fig = plt.figure()
    fig.suptitle(VARNAME[var])
    ax = fig.add_subplot(111) 
    ax.set_ylabel(u'CDF')
    ax.hist(obs_data, color = col[0], bins = 300, density = True, cumulative = True, histtype = 'step', label = 'tmy' )
    ax.plot(lnspc,  cdf_o, '--', color = col[0])
    ax.hist(mod_c, color = col[1], bins = 300, density = True, cumulative = True, histtype = 'step', label = 'model historical' )
    ax.hist(data, color = col[2], bins = 300, density = True, cumulative = True, histtype = 'step', label = 'model predicted' )
    ax.plot(lnspc,  cdf_m, '--', color = col[1])
    ax.hist(data_adj, color = col[3], bins = 300, density = True, cumulative = True, histtype = 'step', label = 'model corrected' )
    ax.legend();
    plt.show()
#=============================================================================
    
    return(data_adj);

   
    
############################################################################### 
def bias2station_EDCDFm( mod_c, data, obs_data, var ):
    import numpy as np;
    from scipy import stats;
    from scipy import interpolate;
    import matplotlib.pyplot as plt;
    from cfg import VARNAME
    
    #If RHS just force between 0 and 100%
    if var == 2 or var == 3:
            lnspc = np.linspace(0.0, 1.0, len(data))
            mod_c = mod_c/100.0
            data = data/100.0
            obs_data = obs_data/100.0
    else:
        minx = min( [min(mod_c), min(data), min(obs_data)] ) - abs(min(mod_c));
        maxx = max( [max(mod_c), max(data), max(obs_data)] ) + abs(max(mod_c));
        lnspc = np.linspace(minx, maxx, len(data))

    # Fits the pdf of the data using a 2-parameter beta distribution adjusting 
    # the beta distribution with variables for α, β, and loc (lower limit), 
    # scale (upper limit - lower limit). Documentation also states that the four
    # parameter distribution scales the range of the beta distribuition. This
    # is similar to what Li et al. (2010) and Watterson (2008) describe where
    # their four parameters are  α, β, and the minimum and maximum value of x. 
    #For purposes here, these two methods are considered to be equivalent. 
    
    abmc,bbmc,cbmc,dbmc = stats.beta.fit(mod_c);    # Model current climate
    abmp,bbmp,cbmp,dbmp = stats.beta.fit(data);     # Model predicted 
    abo,bbo,cbo,dbo     = stats.beta.fit(obs_data); # Observed current climate 

    # The cumulative distrubtion functions! 
    cdf_c = stats.beta.cdf(lnspc, abmc,bbmc,cbmc,dbmc);  # Model current climate
    cdf_p = stats.beta.cdf(lnspc, abmp,bbmp,cbmp,dbmp);  # Model predicted 
    cdf_o = stats.beta.cdf(lnspc, abo,bbo,cbo,dbo);      # Observed current climate

    # Prevent rounding issues and force the top and bottom points to the edges 
    if var == 2 or var == 3:
        if cdf_c[0] < 0.001: cdf_c[0] = 0.0; 
        if cdf_c[-1] > 0.999: cdf_c[-1] = 1.0;
        if cdf_o[0] < 0.001: cdf_o[0] = 0.0; 
        if cdf_o[-1] > 0.999: cdf_o[-1] = 1.0;
        
    # Create functions to act as the inverse of the cdfs
    Finv_c = interpolate.interp1d( cdf_c, lnspc );
    Finv_o = interpolate.interp1d( cdf_o, lnspc );
    
    #subtract and interpolate values to get the correct x adjst via Li et al. 2010
    OCMP = Finv_o(stats.beta.cdf(data, abmp, bbmp, cbmp, dbmp));
    MCMP = Finv_c(stats.beta.cdf(data, abmp, bbmp, cbmp, dbmp));
    
    # The data adjustment!
    data_adj = data + OCMP - MCMP;
 
# =============================================================================
    
    col = ['#ffbb78', '#aec7e8', '#ff7f0e', '#c5b0d5'];
    fig = plt.figure()
    fig.suptitle(VARNAME[var])
    ax = fig.add_subplot(111) 
    ax.set_xlabel(u'Temperatures (C)')
    ax.set_ylabel(u'CDF')
    ax.hist(obs_data, color = col[0], bins = 300, density = True, cumulative = True, histtype = 'step', label = 'tmy' )
    ax.plot(lnspc,  cdf_o, '--', color = col[0])
    ax.hist(mod_c, color = col[1], bins = 300, density = True, cumulative = True, histtype = 'step', label = 'model historical' )
    ax.plot(lnspc, cdf_c, '--', color = col[1])
    ax.hist(data, color = col[2], bins = 300, density = True, cumulative = True, histtype = 'step', label = 'model predicted' )
    ax.plot(lnspc, cdf_p, '--', color = col[2])
    ax.hist(data_adj, color = col[3], bins = 300, density = True, cumulative = True, histtype = 'step', label = 'model corrected' )
    ax.legend();
    plt.show()
#=============================================================================
    if var == 2 or var == 3: return( 100.0 * data_adj)

    else: return( data_adj );


    
    
    