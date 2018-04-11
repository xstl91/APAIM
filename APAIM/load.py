#!/usr/bin/env python
'''
Functions for loading energies.
'''

import os
from .constant import *
from .config import *
from .count_inter import *

def load_single(file_load = single_E_file):
    '''
    Loading isolated binding energies for each adsorbate.
    The file constructs as a binding type and the corresponding energy per line, e.g. N_brg  0.000.
        Argument:
            file_load: string for the route and file name
                       default as single_E_file defined in constant.py
        Return:
            dict of binding types and their energies, e.g. {'N_brg':0.000, ...}
    '''
    single_E = {}
    with open(file_load) as f:
        for line in f:
            k,v = line.split()
            single_E[k] = float(v)
    return single_E


def load_inter(dir_load = inter_directory):
    '''
    Loading interaction energies.
    The directory contains interaction files named after interaction pair type, e.g. N-N.
    Each file constructs as a interaction type and the corresponding energy per line, e.g. hh1  0.000.
        Argument:
            dir_load: string for the route to the directory containing interaction files
                      default as inter_directory defined in constant.py
        Return:
            dict of interaction, e.g. { 'N-N':{'hh1':0.000, 'hh1f':0.000, ...}, 'N-O':{...}, ...}
    '''
    inter_E = { x:{} for x in species_pairs}
    for pair in species_pairs:
        with open(dir_load + '/' + pair) as f:
            for line in f:
                k,v = line.split()
                inter_E[pair][k] = float(v)
    return inter_E


def load_testfile(file_load, modify = None):
    '''
    Loading testing dataset.
    The file should be named after species_str, e.g. 2O-1NO, with optional suffix after '_' used for remark.
    The file constructs as a configuration and the corresponding binding energy per line, e.g. 22-24-55 0.000.
        Argument:
            file_load: string for the route and file name
            modify   : string for modification type to be used, refer to inter_modify
        Return:
            list of instances of class Config
    '''
    config_all = []
    species_str = file_load.split('/')[-1].split('_')[0]
    with open(file_load) as f:
        for line in f:
            config,energy = line.split()
            config_all.append(Config(species_str, config, float(energy), modify))
    return config_all


def load_testdir(dir_load, modify = None):
    '''
    Loading testing dataset.
    The directory contains energy files named after species_str, e.g. 2O-1NO, with optional suffix after '_'.
    All the files under this directory will be loaded.
    Each file constructs as a configuration and the corresponding binding energy per line, e.g. 22-24-55 0.000.
        Argument:
            dir_load: string for the route to the directory containing testing dataset files
            modify  : string for modification type to be used, refer to inter_modify
        Return:
            list of instances of class Config
    '''
    config_all = []
    for species_str in ( x.rstrip() for x in os.popen('ls ' + dir_load) if x[0].isdigit() ):
        with open(dir_load + '/' + species_str) as f:
            for line in f:
                config,energy = line.split()
                config_all.append(Config(species_str.split('_')[0], config, float(energy), modify))
    return config_all
