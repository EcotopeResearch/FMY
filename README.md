
# Future Meteorological Year (FMY)

Writen in Python 3.7

A script to make future hourly weather from the typical meteorological year (TMY) and global climate models (GCM)

FMY was developed for usage with whole house simulation in mind, it is intended to write FMY files in the .tm2 and .tm3 format. 

## Before Running
1. Open main.py in a python or text editor.
2. Set the paths to the tmy2 directory, and the output paths for graphs and fmy.tm2 files. The script will create the output directories if they don't exist.
3. Set output formats, 'csv', 'tmy2', 'tmy3'
4. Set which stations to look at, numbers correspond to the cities in the city list.
5. Decide on a method (NPCC (1) or Belcher (2)).
6. Decide on a baseline climate ('gcm' or 'tmy')
7. Set a bias correction method.
8. suppress_all_plots do it. 0 makes no plots, 1 creates graphs in the graph folder.
9. Set which models to look at, need to look at at least 1.
10. Set which scenario to look at (1 is RCP4.5, 2 is RCP8.5)
11. Set which variables to transform.
12. Set the baseline years and the future years to look at (recommend to use at least 30 years).


## Runing FMY 

1. Clone the git repository or download and extract the source code.
2. Open a console in the directory or open in a python editor (i.e. Spyder).
3. Type 'python main.py' in the command line, or run in a python editor.