#loadconc.py - possibly these classes will be added to ajustador/loader.py when ready
# -*- coding:utf-8 -*-

from __future__ import print_function, division
import numpy as np

class trace(object):
    def __init__(self, molname, x, y):
        self.molname=molname
        wave=np.rec.fromarrays((x, y), names='x,y')
        self.wave=wave
        #may need features associated with each wave

class CSV_conc(object):
    """Load a series of concentration measurements from a CSV file
    Each CSV file contains data for one or more molecules:
      Time time_units, mol_name1 (nM), [mol_name2]
      read time_units (sec,msec,min allowed) and convert to msec
    need to update (add super object?) to allow multiple files, differing in stim protocol
    """
    def __init__(self, fname,features=[]):
        import pandas as pd
        self.name=fname
        self.features=[]
        csv = pd.read_csv(self.name, index_col=0)
        x_head=csv.index.name.split()
        if len(x_head)>1:
            time_units=x_head[-1]
            if time_units.startswith('sec'):
                time_factor=1e3
            elif time_units.startswith('min'):
                time_factor=1e3*60
            else:
                time_factor=1
            print('x column header: {}, time_units: {}, conversion factor: {}'.format(x_head,time_units,time_factor))
        else:
            time_factor=1
        x = csv.index.values*time_factor #time values
        self.waves = {col.split()[0]:trace(col.split()[0], x, csv[col].values) for col in csv.columns}

