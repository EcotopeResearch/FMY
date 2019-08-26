# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 09:20:27 2019

@author: paul
"""

from cfg import MODELNAME, SCENNAME, METHODNAMES, MONTHS_IN_YEAR;
from month_conversions import months2hrs;

import numpy as np;
import pandas as pd;
import matplotlib.pyplot as plt;
import os;

import metpy.calc as mpcalc
from metpy.units import units

###############################################################################

def adjust_to_fmy(var, methods, tmy, dft, dfam, df1, 
                  hr_plots, ftr_yr, mm, sc, 
                  cc0, graphdir, city_name, suppress_plts):
    """
    Function for adjusting variables from tmy to fmy, which uses various methods 
    for the different variables specified in methods. 
    # var -  a list of variables
    # methods - method to use for conversion of temperatures and relative humidity. 
    # tmy - a single weather class tmy year for a station
    # dft - a data from for the tmy year
    # dfam - the monthly averages for the NPCC method
    # df1 - a data frame of future GCM projections 
    
    # Variables for plotting:
    # hr_plots - a list of variables to be plotted hourly for each month
    # ftr_yr - a list of starting and ending years of interest for the GCM
    # mm - the model number
    # sc - the scenario number
    # cc0 - which current climate is used
    """
    
    # Turn interactive plotting off. Will slow program way down if on.
    plt.ioff()

    #=========================================================
    #                Initialize FMY
    #=========================================================
    tfmy        = pd.DataFrame({'T_tmy'     : tmy.tdry, 
                                'RHS_tmy'   : 100*tmy.rhs,
                                'Tdew_tmy'  : tmy.tdew,
                                'Rg_tmy'    : tmy.tothor,
                                'Rdir_tmy'  : tmy.dirnorm,
                                'Rdiff_tmy' : tmy.difhor
                                });
    
    
    # Check for temperature adjustments
    if 0 in var and 1 in var:
        if methods == 1 or methods == 3:
            T_fmy = adjust_temperature_NPCC(tmy, dft, dfam);
        elif methods == 2:
            T_fmy = adjust_temperature_Belcher(tmy, dft, df1);
    else:
        T_fmy = tmy.tdry;

    # Check for relative humidity adjustments
    if 2 in var and 3 in var:
        if methods == 1 or methods == 3:
            RHS_fmy = adjust_relative_humidity_NPCC(tmy, dft, dfam);
        elif methods == 2:
            RHS_fmy = adjust_relative_humidity_Belcher(tmy, dft, df1);
    else:
        RHS_fmy = 100.0 * tmy.rhs;
        
    # Check for solar radiation adjustments
    if 5 in var:
        Rg_fmy, Rdir_fmy, Rdiff_fmy =  adjust_solar_Belcher( tmy, dft, df1 );
    else:
        Rg_fmy      = tmy.tothor;
        Rdir_fmy    = tmy.dirnorm;
        Rdiff_fmy   = tmy.difhor;

    # Check for dew point adjustments
    if 8 in var:
        if 0 not in var:  
            raise Exception('\nERROR:: Can not find fmy dew point without drybulb, include 0 and 1 in var selection.\n');
        else:
            Tdew_fmy = adjust_dew_point_Belcher(tmy, dft, df1, T_fmy);
    else:
        Tdew_fmy = tmy.tdew;


    #=========================================================
    #                Plots
    #=========================================================
    os.makedirs(graphdir + city_name+ '/', exist_ok=True)
    
    #Temperature Plots
    if 0 in var and 1 in hr_plots and not suppress_plts:
        for ii in range(1,MONTHS_IN_YEAR+1):

            fig = plt.figure(figsize=(11, 7))
            fig.suptitle(METHODNAMES[methods-1]+' Transformed Hourly TMY Temperatures for Month ' + str(ii) + '\n'+ city_name+ ' for ' + 
                         MODELNAME[mm] + ' from ' + str(ftr_yr[0]) +' - ' + 
                         str(ftr_yr[1]) + '\n running scenario: ' + SCENNAME[sc])
            ax = fig.add_subplot(111)
            ax.set_xlabel(u'Hour of Year')
            ax.set_ylabel(u'Temperature (C)')
            ax.ticklabel_format(style='plain')
            ax.plot(tmy.hoy[tmy.month==ii], tmy.tdry[tmy.month==ii])
            ax.plot(tmy.hoy[tmy.month==ii], T_fmy[tmy.month==ii])
            ax.plot(tmy.hoy[tmy.month==ii], T_fmy[tmy.month==ii] - tmy.tdry[tmy.month==ii], '--' )
            ax.grid()
            #plt.show()
            
            fig.savefig(graphdir + city_name+ '/' + city_name+'_'+MODELNAME[mm]+'_'+METHODNAMES[methods-1]+'_Baseline_'+cc0+'_Transformed_Hourly_TMY_Temperatures_for_Month_'+ str(ii))
            plt.close()  

    ########################################################################### 
    #Relative Humidity Plots
    if 2 in hr_plots or 3 in hr_plots and not suppress_plts:
        for ii in range(1,MONTHS_IN_YEAR+1):

            fig = plt.figure(figsize=(11, 7))
            fig.suptitle(METHODNAMES[methods-1]+' Transformed Hourly TMY Relative Humidity for Month ' + str(ii) + '\n'+ city_name+  ' for ' + 
                         MODELNAME[mm] + ' from ' + str(ftr_yr[0]) +' - ' + 
                         str(ftr_yr[1]) + '\n running scenario: ' + SCENNAME[sc])
            ax = fig.add_subplot(111)
            ax.set_xlabel(u'Hour of Year')
            ax.set_ylabel(u'Relative Humidity (%)')
            ax.ticklabel_format(style='plain')
            ax.plot(tmy.hoy[tmy.month==ii], 100*tmy.rhs[tmy.month==ii])
            ax.plot(tmy.hoy[tmy.month==ii], RHS_fmy[tmy.month==ii])
            ax.plot(tmy.hoy[tmy.month==ii], RHS_fmy[tmy.month==ii] -  100*tmy.rhs[tmy.month==ii], '--' )
            ax.grid()
            #plt.show()
            
            fig.savefig(graphdir + city_name+ '/' +city_name+'_'+MODELNAME[mm]+'_'+METHODNAMES[methods-1]+'_Baseline_'+cc0+'_Transformed_Hourly_TMY_Relative_Humidity_for_Month_'+ str(ii))
            plt.close() 
    ###########################################################################
    #Dew Point Plots
    if 8 in hr_plots and not suppress_plts:
        for ii in range(1,MONTHS_IN_YEAR+1):
    
            fig = plt.figure(figsize=(11, 7))
            fig.suptitle('Belcher Hourly TMY Dew Point for Month ' + str(ii) + '\n'+ city_name + ' for ' + 
                         MODELNAME[mm] + ' from ' + str(ftr_yr[0]) +' - ' + 
                         str(ftr_yr[1]) + '\n running scenario: ' + SCENNAME[sc])
            ax = fig.add_subplot(111)
            ax.set_xlabel(u'Hour of Year')
            ax.set_ylabel(u'Dew Point (C)')
            ax.ticklabel_format(style='plain')
            ax.plot(tmy.hoy[tmy.month==ii], tmy.tdew[tmy.month==ii])
            ax.plot(tmy.hoy[tmy.month==ii], Tdew_fmy[tmy.month==ii])
            ax.plot(tmy.hoy[tmy.month==ii], Tdew_fmy[tmy.month==ii] - tmy.tdew[tmy.month==ii], '--' )
            ax.grid()
            #plt.show()
            
            fig.savefig(graphdir + city_name+ '/' +city_name+'_'+MODELNAME[mm]+'_'+METHODNAMES[methods-1]+'_Baseline_'+cc0+'_Transformed_Hourly_TMY_DewPoint_for_Month_'+ str(ii))
            plt.close() 
    ###########################################################################
    # Total Global Solar Plots
    if 5 in hr_plots and not suppress_plts:
        for ii in range(1,MONTHS_IN_YEAR+1):
            # Total horizontal Solar
            fig = plt.figure(figsize=(11, 7))
            fig.suptitle('Belcher Hourly TMY Global Horizontal for Month ' + str(ii) + '\n'+ city_name + ' for ' + 
                         MODELNAME[mm] + ' from ' + str(ftr_yr[0]) +' - ' + 
                         str(ftr_yr[1]) + '\n running scenario: ' + SCENNAME[sc])
            ax = fig.add_subplot(111)
            ax.set_xlabel(u'Hour of Year')
            ax.set_ylabel(u'Total Horizontal Irradiance (Wh/m2)')
            ax.ticklabel_format(style='plain')
            ax.plot(tmy.hoy[tmy.month==ii], tmy.tothor[tmy.month==ii])
            ax.plot(tmy.hoy[tmy.month==ii], Rg_fmy[tmy.month==ii])
            ax.plot(tmy.hoy[tmy.month==ii], Rg_fmy[tmy.month==ii] - tmy.tothor[tmy.month==ii], '--' )
            ax.grid()
            #plt.show()
            
            fig.savefig(graphdir + city_name+ '/' +city_name+'_'+MODELNAME[mm]+'_'+METHODNAMES[methods-1]+'_Baseline_'+cc0+'_Transformed_Hourly_TMY_Total_Horizontal_for_Month_'+ str(ii))
            plt.close() 

            # Direct Normal, dirnorm
            fig = plt.figure(figsize=(11, 7))
            fig.suptitle('Belcher Hourly TMY Direct Normal for Month ' + str(ii) + '\n'+ city_name + ' for ' + 
                         MODELNAME[mm] + ' from ' + str(ftr_yr[0]) +' - ' + 
                         str(ftr_yr[1]) + '\n running scenario: ' + SCENNAME[sc])
            ax = fig.add_subplot(111)
            ax.set_xlabel(u'Hour of Year')
            ax.set_ylabel(u'Direct Normal Irradiance (Wh/m2)')
            ax.ticklabel_format(style='plain')
            ax.plot(tmy.hoy[tmy.month==ii], tmy.dirnorm[tmy.month==ii])
            ax.plot(tmy.hoy[tmy.month==ii], Rdir_fmy[tmy.month==ii])
            ax.plot(tmy.hoy[tmy.month==ii], Rdir_fmy[tmy.month==ii] - tmy.dirnorm[tmy.month==ii], '--' )
            ax.grid()
            #plt.show()
            
            fig.savefig(graphdir + city_name+ '/' +city_name+'_'+MODELNAME[mm]+'_'+METHODNAMES[methods-1]+'_Baseline_'+cc0+'_Transformed_Hourly_TMY_Direct_Normal_for_Month_'+ str(ii))
            plt.close() 

            #Diffuse Horizontal, difhor
            fig = plt.figure(figsize=(11, 7))
            fig.suptitle('Belcher Hourly TMY Diffuse Horizontal for Month ' + str(ii) + '\n'+ city_name + ' for ' + 
                         MODELNAME[mm] + ' from ' + str(ftr_yr[0]) +' - ' + 
                         str(ftr_yr[1]) + '\n running scenario: ' + SCENNAME[sc])
            ax = fig.add_subplot(111)
            ax.set_xlabel(u'Hour of Year')
            ax.set_ylabel(u'Diffuse Horizontal Irradiance (Wh/m2)')
            ax.ticklabel_format(style='plain')
            ax.plot(tmy.hoy[tmy.month==ii], tmy.difhor[tmy.month==ii])
            ax.plot(tmy.hoy[tmy.month==ii], Rdiff_fmy[tmy.month==ii])
            ax.plot(tmy.hoy[tmy.month==ii], Rdiff_fmy[tmy.month==ii] - tmy.difhor[tmy.month==ii], '--' )
            ax.grid()
            #plt.show()
            
            fig.savefig(graphdir + city_name+ '/' +city_name+'_'+MODELNAME[mm]+'_'+METHODNAMES[methods-1]+'_Baseline_'+cc0+'_Transformed_Hourly_TMY_Diffuse_Horizontal_for_Month_'+ str(ii))
            plt.close() 
            
    #=========================================================
    #               Check and Finilize FMY
    #=========================================================
        
    if sum(np.isnan(T_fmy)) > 0:     raise Exception('\nERROR: T_fmy contains NAN for ' + city_name + '\n')
    if sum(np.isnan(RHS_fmy)) > 0:   raise Exception('\nERROR: RHS_fmy contains NAN for ' + city_name + '\n') 
    if sum(np.isnan(Tdew_fmy)) > 0:  raise Exception('\nERROR: Tdew_fmy contains NAN for ' + city_name + '\n') 
    if sum(np.isnan(Rg_fmy)) > 0:    raise Exception('\nERROR: Rg_fmy contains NAN for ' + city_name + '\n') 
    if sum(np.isnan(Rdir_fmy)) > 0:  raise Exception('\nERROR: Rdir_fmy contains NAN for ' + city_name + '\n')   
    if sum(np.isnan(Rdiff_fmy)) > 0: raise Exception('\nERROR: Rdiff_fmy contains NAN for ' + city_name + '\n') 
    
    if sum(RHS_fmy < 0) > 0 or sum(RHS_fmy > 100) > 0:
        print('\nWARNING: The FMY of Relative Humidity is less than 0% ('+str(sum(RHS_fmy < 0))+ 
                ' values) and/or greater than 100% ('+str(sum(RHS_fmy > 100))+
                ' values) for ' + city_name + '\n')
                    
    if sum(Tdew_fmy > T_fmy) > 0:  
        print('\nWARNING: The FMY of Dew Point is greater than the Dry Bulb ('+
               str(sum(Tdew_fmy > T_fmy))+ ' values) for ' + city_name + '')
        if sum(tfmy.Tdew_tmy > tfmy.T_tmy) > 0:
            print('\nWARNING: This is a result of the TMY data where the Dew Point is greater than Dry Bulb (' + 
                  str(sum(tfmy.Tdew_tmy > tfmy.T_tmy))+ ' values) for ' + city_name + '\n')
            
    
    if sum(Rdiff_fmy > Rg_fmy) > 0:
        print('\nWARNING: The FMY of the diffuse horizontal radiation is greater than the total global horizontal radiation (' + 
               str(sum(Rdiff_fmy > Rg_fmy))+ ' values) for ' + city_name + '')
        if sum(tfmy.Rdiff_tmy > tfmy.Rg_tmy) > 0:
            print('\nWARNING: This is a result of the TMY data where the diffuse horizontal radiation is greater than the total global horizontal radiation (' + 
                  str(sum(tfmy.Rdiff_tmy > tfmy.Rg_tmy))+ ' values) for ' + city_name + '\n')
            
    # Store the new variables in a data frame to export which will become a weather object...
    tfmy['T_fmy']       = T_fmy;
    tfmy['RHS_fmy']     = RHS_fmy;
    tfmy['Tdew_fmy']    = Tdew_fmy;
    tfmy['Rg_fmy']      = Rg_fmy;
    tfmy['Rdir_fmy']    = Rdir_fmy;
    tfmy['Rdiff_fmy']   = Rdiff_fmy;
    
    return(tfmy);

#==============================================================================               
#==============================================================================     
            
def adjust_temperature_NPCC(tmy, df_c0, df_avgmth):    
    T_fmy       = np.array([]);

    # Find A, B,  for each month and plot!
    # Monthly Temperature Constants for Transformation
    df_avgmth['A_tas'] = ( df_avgmth.tasmax_c - df_avgmth.tasmin_c ) / ( df_c0.tasmax - df_c0.tasmin )               
    df_avgmth['B_tas'] = df_avgmth.tasmax_c - df_avgmth.A_tas * df_c0.tasmax;
    
    for ii in range(1,MONTHS_IN_YEAR+1):
        T_fmy = np.append(T_fmy, 
                          df_avgmth.A_tas[ii-1]*tmy.tdry[tmy.month==ii] + df_avgmth.B_tas[ii-1])
    
    return( T_fmy );
    
#==============================================================================               
#==============================================================================   
        
def adjust_relative_humidity_NPCC(tmy, df_c0, df_avgmth):
    RHS_fmy     = np.array([]);

    # Monthly Relative Humidity Constants for Transformation
    # Find C, D,  for each month and plot!
    df_avgmth['C_tas'] = ( df_avgmth.rhsmax_c - df_avgmth.rhsmin_c ) / ( df_c0.rhsmax - df_c0.rhsmin)               
    df_avgmth['D_tas'] = df_avgmth.rhsmax_c/100 - df_avgmth.C_tas * df_c0.rhsmax/100; #Divide by 100 because to go from %
        
    for ii in range(1,MONTHS_IN_YEAR+1):
        RHS_fmy = np.append(RHS_fmy, 100*(df_avgmth.C_tas[ii-1]*tmy.rhs[tmy.month==ii] + df_avgmth.D_tas[ii-1]))
    
    # Force high values down.
    RHS_fmy[RHS_fmy>100] = 100;

    return( RHS_fmy );
    
#==============================================================================               
#==============================================================================  
                            
#=========================================================
#               Adjust Belecher Variables
#========================================================= 
    
def adjust_temperature_Belcher(tmy, df_c0, df_cf):
    T_fmy       = np.array([]);
 
    # Grab monthly means
    temp1 = df_cf.groupby(['month']).mean().reset_index()
    tempt = df_c0.groupby(['month']).mean().reset_index()

    # Adjust Temperatures
    #Scaling
    adbt =  months2hrs(((temp1.tasmax_c - tempt.tasmax) - (temp1.tasmin_c - tempt.tasmin)) / ( tempt.tasmax - tempt.tasmin));
    #fmy dry bulb ttem
    T_fmy = tmy.tdry + months2hrs((temp1.tasmean - tempt.tasmean)) + adbt*(tmy.tdry - months2hrs(tempt.tasmean))
                         
    return( T_fmy );
   
#==============================================================================               
#==============================================================================    
    
def adjust_relative_humidity_Belcher(tmy, df_c0, df_cf):
    RHS_fmy     = np.array([]);

    # Grab monthly means
    temp1 = df_cf.groupby(['month']).mean().reset_index()
    tempt = df_c0.groupby(['month']).mean().reset_index()
             
    # Adjust Relative Humidity
    # Scaling
    #arhs = months2hrs( ((temp1.rhsmax_c - tempt.rhsmax) - (temp1.rhsmin_c - tempt.rhsmin) )/ ( tempt.rhsmax - tempt.rhsmin) );
    arhs = months2hrs( 1 + ((temp1.rhsmean) - (tempt.rhsmean) ) / 100. );
    
    #fmy relative humdity
    RHS_fmy = 100.0 * arhs * tmy.rhs;
    
    #Force high values down
    RHS_fmy[RHS_fmy>100] = 100.0;
    #Force low values up
    RHS_fmy[RHS_fmy<0] = 0.0;
    
    return( RHS_fmy );

#==============================================================================               
#==============================================================================
    
def adjust_dew_point_Belcher( tmy, df_c0, df_cf, T_fmy):
    Tdew_fmy    = np.array([]);
    
    # Grab monthly means
    temp1 = df_cf.groupby(['month']).mean().reset_index()
    tempt = df_c0.groupby(['month']).mean().reset_index()
    
    # Scaling
    ahuss = months2hrs( (temp1.huss_c / tempt.huss) )
    #ahuss = months2hrs( 1 + (temp1.huss - tempt.huss) )
    
    # FMY Specific Humidity
    huss_fmy = ahuss * tmy.huss;
    
    #Force high values down by finding specific humidity for 100% relative humidity
#    saturation_huss = mpcalc.specific_humidity_from_mixing_ratio( 
#            mpcalc.mixing_ratio_from_relative_humidity(1.00, T_fmy * units.degC, tmy.press * units.mbar)).magnitude;
#    huss_fmy[huss_fmy > saturation_huss] = saturation_huss[huss_fmy > saturation_huss];
#
#    #Force low values up by finding specific humidity for 0% relative humidity
#    dry_huss = mpcalc.specific_humidity_from_mixing_ratio( 
#            mpcalc.mixing_ratio_from_relative_humidity(0.00, T_fmy * units.degC, tmy.press * units.mbar)).magnitude;
#    huss_fmy[huss_fmy < dry_huss] = dry_huss[huss_fmy < dry_huss];

    # FMY dew point          
    Tdew_fmy = mpcalc.dewpoint_from_specific_humidity(huss_fmy, T_fmy * units.degC, tmy.press * units.mbar).magnitude
    
    #Force high values down
    Tdew_fmy[Tdew_fmy > T_fmy] = T_fmy[Tdew_fmy > T_fmy];
        
    return( Tdew_fmy );
   
#==============================================================================               
#==============================================================================
    
def adjust_solar_Belcher( tmy, df_c0, df_cf ):
    Rg_fmy      = np.array([]);
    Rdir_fmy    = np.array([]);
    Rdiff_fmy   = np.array([]);
    
    # Grab monthly means
    temp1 = df_cf.groupby(['month']).mean().reset_index()
    tempt = df_c0.groupby(['month']).mean().reset_index()
    
    # Scaling
    aRg = months2hrs( 1 + (temp1.rsds_c - tempt.rsds) / tempt.rsds)
           
    # FMY Solar Variables         
    Rg_fmy      = aRg * tmy.tothor;
    Rdir_fmy    = aRg * tmy.dirnorm;
    Rdiff_fmy   = aRg * tmy.difhor;
        
    return( Rg_fmy, Rdir_fmy, Rdiff_fmy );
   