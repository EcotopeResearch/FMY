"""
Created on Wed Aug 14 16:42:45 2019

@author: paul
"""

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
import os
import subprocess
import importlib

#==============================================================================
#                           SETTINGS
#==============================================================================

#------------------------------------------------------------------------------
#   PATHS
#------------------------------------------------------------------------------

# Working directory, where the code lives.
workingdir = 'R:/NPCC/2019_Future_Meteorological_Years (19-028-BL)/Python/FMY'

# Path to local TMY data in .tm2 format
weatherpath = "E:/TMY2DATA/"; 

# Output path for graphs
graphpath = "../Graphs_FMY_all_models/";      

# Output path for future .tm2 and .tm3 files
outputpath = "../output_FMY_all_models/"  


#------------------------------------------------------------------------------
# Set Cities: 
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
#stations    = [ 0, 6 ]; 

# Years to average over
tmy3_years   = [ 1976, 2005 ]; # Years that the tmy3 weather files are taken from. This is the baseline time frame too.
future_years = [ 2020, 2049 ]; # Future years from 2006 - 2099.

#------------------------------------------------------------------------------
#   Output formats
#------------------------------------------------------------------------------

# FMY output formats, will print to csv, tmy2, and tmy3. Is a list of strings,
#  containing ['csv','tmy2', 'tmy3'], and any combination of them or none.
outformats = ['csv'];

# Which Variables hourly plots
hourly_plots = [0, 1, 2, 3, 5, 8];

suppress_all_plots = 0; #0 plots print to pdf which impeeds speed, 1 no plots are made at all.

#------------------------------------------------------------------------------
#   MODEL/VAR/SCEN CHOICES
#------------------------------------------------------------------------------
# bcc-csm1-1 (0) did not have daily data available for 12/31/2099 for the RCP8.5 scenario only.
# CCSM4 (4) and NorESM1-M (19) did not have relative humidity available at daily timescales.
models      = [1, 2, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18] # All complete sets
#models      = [3, 4, 5, 6, 8, 9, 10, 11, 12, 15] #The RMJOC-II “10”. Note issues with model 4 in scenario=1
#models      = [3] # for testing

scenarios        = [ 2 ] # 1 is RCP4.5, 2 is RCP8.5

variables         = [0, 1, 2, 3, 5, 8]
#var         = [ 0, 1 ]

#------------------------------------------------------------------------------
#   Experimental Method Options (Leave as is if you aren't sure)
#------------------------------------------------------------------------------
# Affects plots and loads
method      = 2; # 1 is NPCC, 
                 # 2 is Belcher  
                 # 3 is adjusted NPCC with monthly means of daily max and min
                 
# Set to 'tmy' or 'gcm' to set which current climate is used to calculate the ajustments to the mean for            
which_current_climate = 'gcm';

interpolate_to_station = True

# 0 no correction / 1 is EDCDFm / 2 should be CDFm / 3 is linear Hawkins
bias_correction_method = 0;

#==============================================================================
#                           MAIN FUNCTION
#==============================================================================

def main( wpth,  gpth, opth, outfrmt, 
                  stats, mods, scen, var, base_years, future_years,
                  meth, which_cc, interp_to, 
                  bias_correct, hrly_plts, supp_plts ):
    
    from futureWeather import futureWeather;
    
    pkgs = ['numpy',
            'pandas',
            'matplotlib',
            'seaborn',
            'datetime',
            'time',
            'metpy',
            'netCDF4']
    
    # Check for packages
    if importlib.util.find_spec('pip') is None: 
       print("installing pip")
       cmd = "sudo easy_install pip"
       os.system(cmd)
            
    for package in pkgs:
        if importlib.util.find_spec(package) is None: 
            print("installing " + package)
            subprocess.check_call(['python', '-m', 'pip', 'install', package])   
            
    # Run program    
    futureWeather( wpth,  gpth, opth, outfrmt, 
                  stats, mods, scen, var, base_years, future_years,
                  meth, which_cc, interp_to, 
                  bias_correct, hrly_plts, supp_plts );
      
                  
                  


#==============================================================================
#                           Code Below
#==============================================================================


#------------------------------------------------------------------------------
#       Check Directories.
#------------------------------------------------------------------------------

os.chdir(workingdir)

os.makedirs(graphpath, exist_ok=True)
os.makedirs(outputpath, exist_ok=True)
    
#------------------------------------------------------------------------------
#       Run scripts
#------------------------------------------------------------------------------
if __name__ == "__main__":
    main( weatherpath,  graphpath, outputpath, outformats, 
                  stations, models, scenarios, variables, 
                  tmy3_years, future_years,
                  method, which_current_climate, interpolate_to_station, 
                  bias_correction_method,
                  hourly_plots, suppress_all_plots );
    
    
    