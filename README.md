
# Future Meteorological Year (FMY)

Writen in Python 3.7

A script to make future hourly weather from the typical meteorological year (TMY) and global climate models (GCM)

FMY was developed for usage with simulation in mind, it is intended to read TMY2 and TMY3 files and to write FMY files in the TMY2 and TMY3 format for SEEM. 

## Before Running
1. Open main.py in a python or text editor.
2. Set the paths to the data directory, and the output paths for graphs and fmy.tm2 files. The script will create the output directories if they don't exist.
3. Set output formats, 'csv' or 'tmy'
4. Set which stations to look at, numbers correspond to the cities in the city list.
5. Set suppress_all_plots do it. 0 makes no plots, 1 creates graphs in the graph folder.
6. Set which models to look at, need to look at at least 1.
7. Set which scenario to look at (1 is RCP4.5, 2 is RCP8.5)
8. Set which variables to transform, a typical usage will use numbers [0,1,2,3,5,8].
9. Set the baseline years and the future years to look at (recommend to use at least 30 years).


## Runing FMY 

1. Clone the git repository or download and extract the source code.
2. Open a console in the directory or open in a python editor (i.e. Spyder).
3. Type 'python main.py' in the command line, or run in a python editor.