
## Future Meteorological Year (FMY)

Writen in Python 3.7

A script to make future hourly weather from the typical meteorological year (TMY) and global climate models (GCM)

FMY was developed for usage with housing simulation in mind, it is intended to read TMY2 and TMY3 files and to write FMY files in the TMY2 and TMY3 format for SEEM. 

## Set up
Clone the git repository locally or download and extract the code.

## Before Running
1. Open main.py in a python or text editor like Spyder.
2. Set the working directory to the location of the FMY folder. Optionally set the paths to the data directory, and the output paths for graphs and FMY files. The script will create the output directories if they don't exist.
3. Set which stations to look at, numbers correspond to the cities in the city list.
4. Set suppress_all_plots to 0 makes no plots, 1 creates graphs in the graph folder.
5. Set which models to look at, need to look at at least 1.
6. Set the baseline years and the future years to look at (recommend to use at least 30 years).
7. Set which scenario to look at (1 is RCP4.5, 2 is RCP8.5).
8. Set which variables to transform, a typical usage will use numbers [0,1,2,3,5,8].
9. Set output formats, 'csv' or 'tmy', FMY does not convert between TMY2 and TMY3.
Optional: Set download_data to True, which will download GCM data to the working directory, which is useful for running in batch since MACA will boot you if you're hogging their servers.

## Runing FMY 
Open main.py in a python development evironment (i.e. Spyder) and run.

or

Use the command line to move to the FMY directory, type 'python main.py' in the command line to run.

## Description

This script uses OPeNDAP to download the specified subset of 
the MACAv2-METDATA data then applies multiple methods for temporally downscalling 
data to hourly features.

The program reads weather variables from MACA:
	0.  Max Temperature
	1.  Min Temperature
	2.  Max Relative Humidity
	3.  Min Relative Humidity
	4.  Precipitation (Not FMY Supported)
	5.  Surface Downwelling Shortwave Flux in Air 
	6.  Eastward Wind Component (Not FMY Supported)
	7.  Northward Wind Component (Not FMY Supported)
	8.  Specific Humidity

And prints into TMY3  and TMY2 formatting. 

The script takes cities: 
    0.  Seattle (WA)
    1.  Corvallis (OR)
    2.  Boise (ID)
    3.  Redmond(OR)
    4.  Elko (NV)
    5.  Burley (ID)
    6.  Soda_Springs (ID)
    7.  Havre (MT)
    8.  Miles City (MT)
    9.  Portland (OR)
    10. Spokane (WA)
    11. Kalispell (MT)
   
Workflow:
    Set specifics (i.e. variables, cities, models, scenarios, and methods)
    Load TMY file
    Load GCM historical and future
    Determine future peroid of interest and cull the data
    Adjust the variables asked to be adjusted
    Write output file in specified formats
