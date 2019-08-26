# -*- coding: utf-8 -*-
"""
Weather Class for reading and writing TMY2 and TMY3 formatted files.

Note you will need the TMY2/3 file that you want to convert to FMY, to read and
to write to FMY. 
"""
from cfg import HRS_IN_MONTH, HRS_IN_YEAR, DAYS_IN_MONTH, DAYS_IN_YEAR, MONTHS_IN_YEAR, TMY3NUMBER, TMY2EXT, TMY3EXT,CITY
import numpy as np;
import pandas as pd;
import metpy.calc as mpcalc
from metpy.units import units

#Should I bother to initialize everything to length 8760? it'll run faster but...
class weather:
    
    def __init__(self, weatherpath, city, tmy23):
        self.city = city;
        self.weatherpath = weatherpath;
        
        if tmy23 == 2:
            self.file_ext = TMY2EXT;
        elif tmy23 == 3:
            self.file_ext = TMY3EXT;
        else:
            print('\nERROR: Invalid TMY Data type chosen, should be 2 or 3, you choose '
                  + str(tmy23) +'.');
                  
        self.lon        = 0.0;
        self.lat        = 0.0;
        
        self.hoy        = np.array(list(range(0, HRS_IN_YEAR)));
        self.doy        = np.floor(np.array(range(0,HRS_IN_YEAR))/24) % DAYS_IN_YEAR;
        
        self.dom        = np.array([1] * 24 * DAYS_IN_MONTH[0] + [2] * 24 * DAYS_IN_MONTH[1] + \
                            [3] * 24 * DAYS_IN_MONTH[2] + [4] * 24 * DAYS_IN_MONTH[3] + \
                            [5] * 24 * DAYS_IN_MONTH[4] + [6] * 24 * DAYS_IN_MONTH[5] + \
                            [7] * 24 * DAYS_IN_MONTH[6] + [8] * 24 * DAYS_IN_MONTH[7] + \
                            [9] * 24 * DAYS_IN_MONTH[8] + [10] * 24 * DAYS_IN_MONTH[9] + \
                            [11] * 24 * DAYS_IN_MONTH[10] + [12] * 24 * DAYS_IN_MONTH[11]);
        
        self.month      = np.array([1] * HRS_IN_MONTH[0] + [2] * HRS_IN_MONTH[1] + \
                            [3] * HRS_IN_MONTH[2] + [4] * HRS_IN_MONTH[3] + \
                            [5] * HRS_IN_MONTH[4] + [6] * HRS_IN_MONTH[5] + \
                            [7] * HRS_IN_MONTH[6] + [8] * HRS_IN_MONTH[7] + \
                            [9] * HRS_IN_MONTH[8] + [10] * HRS_IN_MONTH[9] + \
                            [11] * HRS_IN_MONTH[10] + [12] * HRS_IN_MONTH[11]);
                            
        self.tdry       = np.array([0.0] * HRS_IN_YEAR); #tdrybulb
        self.rhs        = np.array([0.0] * HRS_IN_YEAR); #relative humidity
        self.tdew       = np.array([0.0] * HRS_IN_YEAR); #dew point
        
        self.huss       = np.array([0.0] * HRS_IN_YEAR); #specfific humidity to be caluculated from tdew

        self.press      = np.array([0.0] * HRS_IN_YEAR); #Atmospheric Pressure

        self.tothor     = np.array([0.0] * HRS_IN_YEAR); #Total horizontal solar
        self.dirnorm    = np.array([0.0] * HRS_IN_YEAR); #Direct normal solar
        self.difhor     = np.array([0.0] * HRS_IN_YEAR); #Diffuse Horizontal Solar
        
       # self.wind_speed = np.array([0.0] * HRS_IN_YEAR); #windspeed m/s
       # self.wind_dir   = np.array([0.0] * HRS_IN_YEAR); #wind direction azimuth
            
       # self.cloud      = np.array([0.0] * HRS_IN_YEAR); #Could cover fraction
        
        
###############################################################################

    def get_weather(self):
        
        #------------------------------------------------------
        #------------------------------------------------------
        # Load TMY2
        if self.file_ext == TMY2EXT:
            f = open(self.weatherpath+self.city+self.file_ext,'r');

            # Header read for lat and lon
            head = f.readline();
            self.lat = int(head[39:41]) + int(head[42:44])/60.0;
            self.lon = int(head[47:50]) + int(head[51:53])/60.0;
            
            line = f.readline()
            ind = 0;
            while line:
                
                # Process the line
                self.tothor[ind]       = float(line[17:21])             #Total horizontal solar Wh/m2
                self.dirnorm[ind]      = float(line[23:27])             #Direct normal solar Wh/m2
                self.difhor[ind]       = float(line[29:33])	            #Diffuse Horizontal Solar Wh/m2
                
                self.tdry[ind]         = float(line[67:71]) * 0.1;	    #tdrybulb (deg C)
                self.rhs[ind]          = float(line[79:82]) * 0.01;		#relative humidity (%)
                self.tdew[ind]         = float(line[73:77]) * 0.1;		#tdew (deg C) to conform with TB code
    
                self.press[ind]        = float(line[84:88]);			#atmospheric pressure (mbar) mb = 100 Pascals
                #self.wind_speed[ind]   = float(line[95:98]) * 0.1;		#windspeed m/s
                #self.wind_dir[ind]     = float(line[90:93]);   			#wind direction azimuth
                
                #self.cloud[ind]        = float(line[59:61])/10.0;		    #Could cover fraction
                #wd.ocloud       = getint(line,63,2)/10.0;		        #Opaque cloud cover fraction
                #wd.ceilht       = getint(line,106,5);		            #Cloud ceiling height m
                
                # Calculate specfic humidity from dry bulb, dew point, and atm pressure using MetPy
                self.huss[ind] = mpcalc.specific_humidity_from_mixing_ratio(
                                    mpcalc.mixing_ratio_from_relative_humidity(
                                            mpcalc.relative_humidity_from_dewpoint(
                                                    self.tdry[ind] * units.degC, 
                                                    self.tdew[ind] * units.degC), 
                                             self.tdry[ind]  * units.degC,
                                             self.press[ind]  * units.mbar)
                                            )
                
                #Next line
                line = f.readline();
                ind = ind + 1;
    
            f.close();
            
        #------------------------------------------------------
        #------------------------------------------------------
        # Load TMY3
        elif self.file_ext == TMY3EXT:
            f = open(self.weatherpath+TMY3NUMBER[CITY.index(self.city)]+self.file_ext,'r');
            
            # Header read for lat and lon
            head = f.readline().split(',');
            self.lat = float(head[4]) ;
            self.lon = float(head[5]) ;
            
            #Burn a line for the second part of the header.
            line = f.readline() 
            
            line = f.readline(); 
            ind = 0;
            while line:
                
                line = line.split(',');
                
                if len(line) < 20:
                    print('line is short!')
                # Process the line
                self.tothor[ind]       = float(line[4])             #Total horizontal solar Wh/m2
                self.dirnorm[ind]      = float(line[7])             #Direct normal solar Wh/m2
                self.difhor[ind]       = float(line[10])	            #Diffuse Horizontal Solar Wh/m2
                
                self.tdry[ind]         = float(line[31]);	    #tdrybulb (deg C)
                self.rhs[ind]          = float(line[37]) * 0.01;		#relative humidity (%)
                self.tdew[ind]         = float(line[34]);		#tdew (deg C) to conform with TB code
    
                self.press[ind]        = float(line[40]);			#atmospheric pressure (mbar) mb = 100 Pascals
                #self.wind_speed[ind]   = float(line[46]);		#windspeed m/s
                #self.wind_dir[ind]     = float(line[43]);   			#wind direction azimuth
                
              
                # Calculate specfic humidity from dry bulb, dew point, and atm pressure using MetPy
                self.huss[ind] = mpcalc.specific_humidity_from_mixing_ratio(
                                    mpcalc.mixing_ratio_from_relative_humidity(
                                            mpcalc.relative_humidity_from_dewpoint(
                                                    self.tdry[ind] * units.degC, 
                                                    self.tdew[ind] * units.degC), 
                                             self.tdry[ind]  * units.degC,
                                             self.press[ind]  * units.mbar)
                                            )
                
                #Next line
                line = f.readline();
                ind = ind + 1;
    
            f.close();

###############################################################################

    def write_to_monthly_df(self,station, daily_maxs):
        #writes the weather class to a data frame.
        #Gets one station at the time
        # daily_maxs asks are we looking at the monthly means of the daily max/
        # mins OR are we looking at the max/min throughout every month. 
        
        ddff = pd.DataFrame(
                            {'station': [station]*HRS_IN_YEAR, 
                             'month'  : self.month,
                             'doy'    : self.doy,
                             'hoy'    : self.hoy,
                             
                             'tdry'   : self.tdry, #
                             
                             'tdew'   : self.tdew, 
                             'rhs'    : self.rhs,  #
                             'press'  : self.press,
                             
                             'huss'   : self.huss,
                             
                             'tothor' : self.tothor, #
                             'dirnorm': self.dirnorm, #
                             'difhor' : self.difhor #
                             
                             #'cloud'  : self.cloud
                             })
        
        # Step this down to monthly data so have tasmean/max/min
        #  ('tasmax', 'tasmin', 'rhsmax', 'rhsmin', 'pr', 'rsds', 'uas', 'vas', 'huss')
        
        # Get monthly maxs and mins for each variable
        if daily_maxs == 0: 
            # Temperatures
            tempT = ddff.groupby(['station','month']).tdry.agg(['max', 'min', 'mean'])
       
            # Relative humidities
            tempR = ddff.groupby(['station','month']).rhs.agg(['max', 'min', 'mean'])
            
            #Specific humidity
            tempH = ddff.groupby(['station','month']).huss.agg(['max', 'min', 'mean'])
            
            # Solar Variables
            tempTotHor  = ddff.groupby(['station','month']).tothor.agg(['mean'])
            tempDir     = ddff.groupby(['station','month']).dirnorm.agg(['mean'])
            tempDiff    = ddff.groupby(['station','month']).difhor.agg(['mean'])

        # Grab average monthly maxs and mins 
        else:
            # Temperatures
            tempT = ddff.groupby(['station','doy','month']).tdry.agg(['max', 'min', 'mean'])
            tempT = tempT.groupby(level=[0,2]).mean() # Returns the monthly averages

            # Relative humidities
            tempR = ddff.groupby(['station','doy','month']).rhs.agg(['max', 'min', 'mean'])
            tempR = tempR.groupby(level=[0,2]).mean()
            
            #Specific humidity
            tempH = ddff.groupby(['station','doy','month']).huss.agg(['max', 'min', 'mean'])
            tempH = tempH.groupby(level=[0,2]).mean()

            # Solar Variables
            tempTotHor  = ddff.groupby(['station','doy','month']).tothor.agg(['mean'])
            tempTotHor  = tempTotHor.groupby(level=[0,2]).mean()
            tempDir     = ddff.groupby(['station','doy','month']).dirnorm.agg(['mean'])
            tempDir     = tempDir.groupby(level=[0,2]).mean()
            tempDiff    = ddff.groupby(['station','doy','month']).difhor.agg(['mean'])
            tempDiff    = tempDiff.groupby(level=[0,2]).mean()

     # write to a dataframe
        df = pd.DataFrame({'station'  : [station]*MONTHS_IN_YEAR, 
                           'month'    : np.array(range(1, MONTHS_IN_YEAR+1)), 
                           
                           'tasmax'   : tempT['max'],
                           'tasmin'   : tempT['min'],
                           'tasmean'  : tempT['mean'],
                           
                           'rhsmax'   : tempR['max']*100,
                           'rhsmin'   : tempR['min']*100,
                           'rhsmean'  : tempR['mean']*100,
                           
                           'huss'     : tempH['mean'],
                           
                           'rsds'     : tempTotHor['mean'], 
                           'dirnorm'  : tempDir['mean'], 
                           'difhor'   : tempDiff['mean'] 
                           } )
                             
        df = df.reset_index(drop=True);
        
        return(df);
      
        
###############################################################################

    def write_to_daily_df(self, station):
        #writes the weather class to a data frame.
        #Gets one station at the time
        
        ddff = pd.DataFrame(
                            {'station': [station]*HRS_IN_YEAR, 
                             'month'  : self.month,
                             'doy'    : self.doy,
                             'hoy'    : self.hoy,
                             'tdry'   : self.tdry, #
                             
                             'tdew'   : self.tdew, 
                             'rhs'    : self.rhs,  #
                             'press'  : self.press,
                             
                             'huss'   : self.huss,
                             
                             'tothor' : self.tothor, #
                             'dirnorm': self.dirnorm, #
                             'difhor' : self.difhor #
                             # 'cloud'  : self.cloud 
                             })
 
        # Temperatures
        tempT = ddff.groupby(['station','doy']).tdry.agg(['max', 'min', 'mean'])
   
        # Relative humidities
        tempR = ddff.groupby(['station','doy']).rhs.agg(['max', 'min', 'mean'])

        # Specific humidity 
        tempH = ddff.groupby(['station','doy']).huss.agg(['mean'])

        # Solar Variables
        tempTotHor  = ddff.groupby(['station','doy']).tothor.agg(['mean'])
        tempDir     = ddff.groupby(['station','doy']).dirnorm.agg(['mean'])
        tempDiff    = ddff.groupby(['station','doy']).difhor.agg(['mean'])


        df = pd.DataFrame({'station'  : [station]*DAYS_IN_YEAR, 
                           'doy'      : np.array(range(1, DAYS_IN_YEAR+1)), 
                           'tasmax'   : tempT['max'],
                           'tasmin'   : tempT['min'],
                           'tasmean'  : tempT['mean'],
                           
                           'rhsmax'   : tempR['max']*100,
                           'rhsmin'   : tempR['min']*100,
                           'rhsmean'  : tempR['mean']*100,
                           
                           'huss'     : tempH['mean'],
                           
                           'rsds'     : tempTotHor['mean'], 
                           'dirnorm'  : tempDir['mean'], 
                           'difhor'   : tempDiff['mean'] 
                           } )
        
        df = df.reset_index(drop=True);
                 
        return(df);
      
        
 ###############################################################################
 
    def plot_hourly_variable(self,var):
        import matplotlib.pyplot as plt;
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_xlabel(u'Month')
        ax.ticklabel_format(style='plain')
        ax.plot(np.array(self.hoy)/8760*12+1, var,'-')
        plt.xticks(list(range(0,13)))
        plt.show()
###############################################################################
        
    def write_tmy_file(self, outdir, modelname, scenname):
        
        if self.file_ext == TMY2EXT:
            self.write_tmy2( outdir, modelname, scenname)
        elif self.file_ext == TMY3EXT:
            self.write_tmy3( outdir, modelname, scenname)
            
###############################################################################

    def check_tmy2_limits(self):
    # checks the limits that are imposed by TMY2

        # Limits from tmy2 file [min/max]
        tmy2_dbt    = [-50.0, 50.0];
        tmy2_dew    = [-60.0, 30.0];
        tmy2_rhs    = [0.0, 100.0];
        tmy2_ghr    = [0.0, 1200.0]; #Total horizontal solar
        tmy2_dir    = [0.0, 1100.0]; #Direct normal solar
        tmy2_dif    = [0.0, 700.0];  #Diffuse Horizontal Solar
        
        # Switch values over/under the limits to the limits.
        self.tdry = [max(min(x, tmy2_dbt[1]), tmy2_dbt[0]) for x in self.tdry]        
        self.tdew = [max(min(x, tmy2_dew[1]), tmy2_dew[0]) for x in self.tdew]        
        self.rhs  = [max(min(x, tmy2_rhs[1]), tmy2_rhs[0]) for x in self.rhs]        
        
        self.tothor  = [max(min(x, tmy2_ghr[1]), tmy2_ghr[0]) for x in self.tothor] 
        self.dirnorm = [max(min(x, tmy2_dir[1]), tmy2_dir[0]) for x in self.dirnorm] 
        self.difhor  = [max(min(x, tmy2_dif[1]), tmy2_dif[0]) for x in self.difhor] 
        

###############################################################################
    def write_fmy(self, outdir, mod, scen, fyrs):
        #Method to write output of to FMY file.
        
        if ( self.file_ext == TMY2EXT):
            self.write_tmy2(outdir, mod, scen, fyrs);
            
        elif ( self.file_ext == TMY3EXT):
            self.write_tmy3(outdir, mod, scen, fyrs);

###############################################################################
        
        
    def write_tmy2(self, outdir, modelname, scenname, future_yrs):
        
        self.check_tmy2_limits();
        
        f = open(self.weatherpath+self.city+self.file_ext,'r');
        g = open(outdir + self.city + "_future_" + 
                 str(future_yrs[0]) + "_" + str(future_yrs[0]) + "_" +
                 modelname + "_" +  scenname +  self.file_ext , 'w+');

        #####################
        # Write the header verbatum
        head = f.readline();
        g.write(head);
        
        ######################
        line = f.readline()
        ind = 0;
        
        while line:
            
            # Transfrom the string line to a list
            ls = list(line);

            # tdrybulb (deg C)
            ls[67:71]               = str(int(self.tdry[ind]*10)).zfill(4); 
            # relative humidity (%)
            ls[79:82]               = str(int(self.rhs[ind])).zfill(3); 
            # tdew (deg C) to conform with TB code
            ls[73:77]               = str(int(self.tdew[ind]*10)).zfill(4); 
            
            # Total horizontal solar Wh/m2
            ls[17:21]               = str(int(self.tothor[ind])).zfill(4); 
            # Direct normal solar Wh/m2
            ls[23:27]               = str(int(self.dirnorm[ind])).zfill(4); 
            # Diffuse Horizontal Solar Wh/m2
            ls[29:33]               = str(int(self.difhor[ind])).zfill(4); 

            # Write line
            g.write(''.join(ls));
            
            # Next line
            line = f.readline();
            ind = ind + 1;

        f.close();
        g.close();
        
###############################################################################

    def write_tmy3(self, outdir, modelname, scenname, future_yrs):
                
        f = open(self.weatherpath+TMY3NUMBER[CITY.index(self.city)]+self.file_ext,'r');
        g = open(outdir + self.city + "_future_" + 
                 str(future_yrs[0]) + "_" + str(future_yrs[0]) + "_" +
                 modelname + "_" +  scenname +  self.file_ext , 'w+');
                 
        #####################
        # Write the header lines verbatum
        head = f.readline();
        g.write(head);
        
        head = f.readline();
        g.write(head);
        
        ######################
        line = f.readline();
        ind = 0;
        while line:
            line = line.split(',');
            
            # tdrybulb (deg C)
            line[31]               = str(round(self.tdry[ind],1));
            # relative humidity (%)
            line[37]               = str(round(self.rhs[ind])); 
            # tdew (deg C) to conform with TB code
            line[34]               = str(round(self.tdew[ind],1)); 
            
            # Total horizontal solar Wh/m2
            line[4]               = str(round(self.tothor[ind]));
            # Direct normal solar Wh/m2
            line[7]               = str(round(self.dirnorm[ind]));
            # Diffuse Horizontal Solar Wh/m2
            line[10]              = str(round(self.difhor[ind])); 

            # Write line
            g.write(','.join(line));
            
            # Next line
            line = f.readline();
            ind = ind + 1;

        f.close();
        g.close();
        

