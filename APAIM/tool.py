#!/usr/bin/env python
'''
Tool Functions.
'''

from .constant import *


def get_species(species):
    '''
    Get adsorbate type and number from initial string.
        Argument:
            species: str, e.g. '2N-1O-1NO'
        Return:
            list of adsorbate type and number in initial order, e.g. [['N',2],['O',1],['NO',1]]
    '''
    species_list = []
    for num_mol in species.split('-'):
        split_num = [ x.isdigit() for x in num_mol ].index(False)
        species_list.append([num_mol[split_num:],int(num_mol[:split_num])])
    return species_list


def sort_site(species, config):
    '''
    Sort the sites to standard adsorbate order defined in species_all.
        Argument:
            species: str, e.g. '1O-2N-1NO'
            config : list, e.g. ['36','44','24','76']
        Return:
            list of sorted configuration, e.g. ['24','44','36','76']
    '''
    species_list = get_species(species)
    if sum(x[1] for x in species_list) != len(config):
        raise Exception("Species and configuration don't match.")
    config_tmp = {}
    count = 0
    for mol,num in species_list:
        config_tmp[mol] = sorted(config[count:count+num])
        count += num
    config_sorted = []
    for mol in species_all:
        config_sorted.extend(config_tmp.get(mol,[]))
    return config_sorted


def site_convert(site):
    '''
    Changing site notation between str and list.
        Argument:
            site: str or list, e.g. '43' or [4,3]
        Return:
            correspond notation of another type
    '''
    correspond = {'a':10, 'b':11, 'c':12, 10:'a', 11:'b', 12:'c'}
    if isinstance(site,str):
        site_list = []
        for site_char in site:
            site_list.append( correspond[site_char] if site_char in correspond else int(site_char))
        return site_list
    elif isinstance(site,list):
        site_str = ''
        for site_int in site:
            site_str += ( correspond[site_int] if site_int in correspond else str(site_int))
        return site_str
    else:
        return site


def site_type(site):
    '''
    Detecting type of site.
        Argument:
            site: list, e.g. [2,4]
        Return:
            String of type name, e.g. 'top'
    '''
    if surface == '100' :
        if (site[0] + site[1]) % 2 == 1 : site_name = 'brg'
        elif site[0] % 2 == 0 : site_name = 'hol'
        else : site_name = 'top'
    elif surface == '111' :       
        if (site[0] - 2) % 3 == 0 : site_name = 'hcp'
        elif site[0] % 3 == 0 : site_name = 'fcc'
        else : site_name = 'top'
    return site_name


def site_to_offset(site1, site2):
    '''
    Detecting the shortest offset of two sites.
        Argument:
            site1: list, e.g. [2,2]
            site2: list, e.g. [6,4]
        Return:
            list of offset, e.g. [4,2]
    '''
    offset = []
    for i in [0,1]:
        if surface == '100' : period = 8
        elif i == 0 : period = 12
        else : period = 4
        offset_tmp = site2[i]-site1[i]
        if offset_tmp > period/2  : offset_tmp -= period
        elif offset_tmp <= -period/2 : offset_tmp += period
        offset.append(offset_tmp)
    return offset


def offset_to_site(site, offset):
    '''
    Get new site from offset and center site
        Argument:
            site  : list, e.g. [6,4]
            offset: list, e.g. [2,2]
        Return:
            list of new site, e.g. [8,6]
    '''
    site_new = []
    for i in [0,1]:
        if surface == '100' : period = 8
        elif i == 0 : period = 12
        else : period = 4
        new_tmp = site[i] + offset[i]
        if new_tmp > period : new_tmp -= period
        elif new_tmp <= 0 : new_tmp += period
        site_new.append(new_tmp)
    return site_new
