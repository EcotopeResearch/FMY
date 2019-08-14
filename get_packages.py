# -*- coding: utf-8 -*-
"""
Created on Tue Aug  6 14:02:36 2019

@author: paul
"""

import os

try:
   import pip
except ImportError:
   print("installing pip")
   cmd = "sudo easy_install pip"
   os.system(cmd)


pkgs = ['numpy',
        'pandas',
        'matplotlib',
        'seaborn',
        'datetime',
        'time',
        'metpy',
        'netCDF4']

for package in pkgs:
    try: 
        import package;
    except ImportError:
        print("installing " + package)
        pip.main(['install', package])
