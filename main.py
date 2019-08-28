
"""
#Author:	Paul Kintner
#Updated: 	8/28/2019
#Description: 	This script uses OPeNDAP to download the specified subset of 
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
# File names are state_abbreviation_city_2or3.format, for example WASeattle3.tm2
weatherpath = './TMY/'; #"E:/TMY2DATA/"; 
# Input files are in (2) TMY2 fromat or (3) TMY3 format
load_tmy23 = 3;

# Data path for GCM data if download_data = True
datapath = "./GCM_Data/"; 

# Output path for graphs
graphpath = "./Graphs_FMY/"; 
  
# Output path for future .tm2 and .csv files
outputpath = "./output_FMY/";


#------------------------------------------------------------------------------
#   Input Choices
#------------------------------------------------------------------------------
 
"""
 Set Cities: 
    0. Seattle (WA)
    1. Corvallis (OR)
    2. Boise (ID)
    3. Redmond(OR)
    4. Elko (NV)
    5. Burley (ID)
    6. Soda_Springs (ID)
    7. Havre (MT)
    8. Miles City (MT)
    9. Portland (OR)
    10.Spokane (WA)
    11.Kalispell (MT)
"""
stations    = [ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 ];


#------------------------------------------------------------------------------
#   MODEL/VAR/SCEN CHOICES
#------------------------------------------------------------------------------
"""
 The GCMs available in MACA are:
 0.   bcc-csm1-1
 1.   bcc-csm1-1-m
 2.   BNU-ESM
 3.   CanESM2
 4.   CCSM4
 5.   CNRM-CM5
 6.   CSIRO-Mk3-6-0
 7.   GFDL-ESM2G
 8.   GFDL-ESM2M
 9.   HadGEM2-CC365
 10.  HadGEM2-ES365
 11.  inmcm4
 12.  IPSL-CM5A-MR
 13.  IPSL-CM5A-LR
 14.  IPSL-CM5B-LR
 15.  MIROC5
 16.  MIROC-ESM
 17.  MIROC-ESM-CHEM
 18.  MRI-CGCM3
 19.  NorESM1-M
 
NOTE:
 bcc-csm1-1 (0) did not have daily data available for 12/31/2099 for the RCP8.5 scenario only.
 CCSM4 (4) and NorESM1-M (19) did not have relative humidity available at daily timescales.
"""
#models      = [1, 2, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18] # All complete sets
#models      = [3, 4, 5, 6, 8, 9, 10, 11, 12, 15] #The RMJOC-II “10”. Note issues with model 4 in scenario=0 (historical)
models      = [1, 2, 3, 5, 6, 7, 8, 9, 10, 11, 12, 15] # Selction of 12, RMJOC-II with extras
models      = [5] # for testing

# Years to average over
tmy3_years   = [ 1976, 2005 ]; # Years that the tmy3 weather files are taken from. This is the baseline time frame too.
future_years = [ 2020, 2049 ]; # Future years from 2006 - 2099.

# 1 is RCP4.5, 2 is RCP8.5
scenarios        = [ 1, 2 ] 
scenarios = [2]

# Variables:
#   0. Max Temperature
#   1. Min Temperature
#   2. Max Relative Humidity
#   3. Min Relative Humidity
#   4. Precipitation (Not FMY Supported)
#   5. Surface Downwelling Shortwave Flux in Air 
#   6. Eastward Wind Component (Not FMY Supported)
#   7. Northward Wind Component (Not FMY Supported)
#   8. Specific Humidity
variables         = [0, 1, 2, 3, 5, 8]

#------------------------------------------------------------------------------
#   Output formats
#------------------------------------------------------------------------------

# FMY output formats, will print to csv, tmy2, and tmy3. Is a list of strings,
# containing ['csv','tmy'], and any combination of them or none. 
# csv will output a csv file containing just the variables that are changed
# before and after the change. 'csv' DOES NOT write a TMY3 formated file
# tmy will output a tmy file corresponding to the input file. 
outformats = ['csv', 'tmy'];

# Which Variables hourly plots
hourly_plots =  [0, 1, 2, 3, 5, 8];

suppress_all_plots = True; #0 plots print to pdf which impeeds speed, 1 no plots are made at all.

# If download the GCM dataset, can be useful if you are worried about getting 
# booted from someone's server
download_data = False; 

#------------------------------------------------------------------------------
#   Experimental Method Options (Leave as is if you aren't sure)
#------------------------------------------------------------------------------
# Affects plots and loads
method      = 2; # 1 is NPCC, 
                 # 2 is Belcher  
                 # 3 is adjusted NPCC with monthly means of daily max and min
                 
# Set to 'tmy' or 'gcm' to set which current climate is used to calculate the ajustments to the mean for            
which_current_climate = 'gcm';

interpolate_to_station = True;

# 0 no correction / 1 is EDCDFm / 2 should be CDFm 
bias_correction_method = 0;

#==============================================================================
#                           MAIN FUNCTION
#==============================================================================

def main( wpth, dtpth, gpth, opth, outfrmt,tmy23,
                stats, mods, scen, var, base_years, future_years,
                meth, which_cc, interp_to, 
                bias_correct, hrly_plts, supp_plts, 
                dl_data):
    
    from futureWeather import futureWeather;
    
    pkgs = ['numpy',  'pandas',  'matplotlib',  'seaborn',  'json',
             'datetime',  'time',  'metpy',  'netCDF4'];
    
    # Check for packages
    if importlib.util.find_spec('pip') is None: 
       print("Installing pip... \n")
       cmd = "sudo easy_install pip"
       os.system(cmd)
            
    for package in pkgs:
        if importlib.util.find_spec(package) is None: 
            print("Installing " + package + "... \n")
            subprocess.check_call(['python', '-m', 'pip', 'install', package])   
            
    # Run program    
    futureWeather( wpth, dtpth, gpth, opth, outfrmt, tmy23,
                  stats, mods, scen, var, base_years, future_years,
                  meth, which_cc, interp_to, 
                  bias_correct, hrly_plts, supp_plts,
                  dl_data );
                       
                 

#==============================================================================
#                           Code Below
#==============================================================================


#------------------------------------------------------------------------------
#       Check Directories.
#------------------------------------------------------------------------------

os.chdir(workingdir)

os.makedirs(graphpath, exist_ok=True)
os.makedirs(outputpath, exist_ok=True)
os.makedirs(datapath, exist_ok=True)


#------------------------------------------------------------------------------
#       Check Inputs.
#------------------------------------------------------------------------------
# Check inputs that are lists are lists.
if  type(variables) != list:
    raise Exception("\nERROR: variables must be a list, it should have square brackets around the number(s), i.e. [3]")
if  type(models) != list:
    raise Exception("\nERROR: models must be a list, it should have square brackets around the number(s), i.e. [3]")
if  type(scenarios) != list:
    raise Exception("\nERROR: scenarios must be a list, it should have square brackets around the number(s), i.e. [2]")
if  type(future_years) != list:
    raise Exception("\nERROR: future_years must be a list, it should have square brackets around the number(s), i.e. [2010, 2019]")
if  type(tmy3_years) != list:
    raise Exception("\nERROR: tmy3_years must be a list, it should have square brackets around the number(s), i.e. [1990, 2005]")
if  type(outformats) != list:
    raise Exception("\nERROR: outformats must be a list, it should have square brackets around the number(s), i.e. ['tmy']")
if  type(hourly_plots) != list:
    raise Exception("\nERROR: hourly_plots must be a list, it should have square brackets around the number(s), i.e. [0,1]")

# Check variables 
if 0 in variables or 1 in variables:
    if not (0 in variables and 1 in variables):
        raise Exception("\nERROR: Selecting the max/min temperature without the other add 0/1 to var")  
if 2 in variables or 3 in variables:
    if not (2 in variables and 3 in variables):
        raise Exception("\nERROR:: Selecting the max/min relative humidity without the other add 2/3 to var")   

# Check years
if len(future_years) != 2 or max(future_years) > 2099 or min(future_years) < 2006:
    raise Exception("\nERROR: future_years must be a list of two years bounded by 2006, and 2099, i.e. [2020, 2049]")

if len(tmy3_years) !=  2or max(tmy3_years) > 2005 or min(tmy3_years) < 1950:     
    raise Exception("\nERROR: tmy3_years must be a list of two years bounded by 1950 and 2005, i.e. [1976, 2005]")

#Check things are in bounds 
if min(scenarios) < 1 or max(scenarios) > 2:
    raise Exception('\nERROR: Invalid scenarios chosen, can only choose 1 or 2, but user input was: '+ str(scenarios)+'.\n')
    
    
if not (0 < method <= 3):
    raise Exception('\nERROR: Invalid method chosen, can only choose 1, 2, or 3, but user input was: '+ str(method)+'.\n')
if not (0 <=  bias_correction_method <=3):  
    raise Exception("\nERROR: Invalid bias_correction_method chose can be 0, 1, 2, 3, user input was: "+ str(bias_correction_method)+'.\n')
    
if not (which_current_climate == 'tmy' or which_current_climate == 'TMY' or which_current_climate == 'gcm' or which_current_climate == 'GCM'): 
    raise Exception('\nERROR: which_current_climate is ' + which_current_climate + ' should be set to ''tmy'' or ''gcm''.\n')
    
if not not outformats: # List is not empty
    if not ('tmy' in (name.lower() for name in outformats) or 'csv' in (name.lower() for name in outformats)) :
        raise Exception('\nERROR: Invalid output format chosen.')
else:
    print('\nWarning no output formats chosen, will not write out data to a file.\n');
        
#------------------------------------------------------------------------------
#       Run scripts
#------------------------------------------------------------------------------
if __name__ == "__main__":
    main( weatherpath, datapath, graphpath, outputpath, outformats, load_tmy23, 
                  stations, models, scenarios, variables, 
                  tmy3_years, future_years,
                  method, which_current_climate, interpolate_to_station, 
                  bias_correction_method,
                  hourly_plots, suppress_all_plots,
                  download_data);
    
    
    