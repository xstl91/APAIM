#!/usr/bin/env python
'''
Functions for analyzing the model.
'''

from .constant import *

def contain_inter(config_list, inter_con):
    '''
    Function to get collection of configurations contain certain interaction pairs and types.
        Argument:
            config_list: list of instances of class Config
            inter_con  : dict of interaction pairs and types to collect instances, e.g. { 'N-N':['hh3','bb6'], 'N-O':['bb6'] }
        Return:
            list of instances of class Config
    '''
    configs = []
    for config in config_list:
        for k,v in inter_con.items():
            for i in v:
                if (k not in config.inter_pair or i not in config.inter_pair[k]) and (k not in config.inter_modify or i not in config.inter_modify[k]):
                    break
            else: continue
            break 
        else:
            configs.append(config)
    return configs


def error(config_list, err_type, E_type = 'T'):
    '''
    Function to calculate average error value for predicted binding energy.
        Argument:
            config_list: list of instances of class Config
            err_type   : string to define type of error
                         'A'    : average error value
                         'R'    : RMSE
                         'MAX n': maximum error value(s) with the corresponding instance(s), 'n' for number, e.g. 'MAX 3'
                         'MIN n': minimum error value(s) with the corresponding instance(s), 'n' for number
                         'D'    : distribution of errors in 20 equal bins with the median at 0
            E_type     : energy type to be used for analysis
                         'T': total adsorption energy
                         'A': average adsorption energy, i.e. the energy divided by number of adsorbates
                         'R': relative energy, i.e. the interaction energy percentage
        Return:
            float of error value for err_type 'A' and 'E'
            list of error value(s) and the corresponding instance(s) for err_type 'MAX n' and 'MIN n'
            list of numbers of instances and the overall error range for err_type 'D'
    '''
    err = []
    for config in config_list:
        if E_type == 'T':
            err.append(config.E_err)
        elif E_type == 'A':
            err.append(config.E_err/sum( config.species.values() ))
        elif E_type == 'R':
            err.append(config.E_err / config.E_bind_a)
        else:
            raise Exception('Unsupport choice for E_type')
    if err_type == 'A':
        return sum(err) / len(err)
    elif err_type == 'R':
        return ( sum(x**2 for x in err) / len(err) ) ** 0.5
    else:
        err_con = sorted( zip(err, config_list) )
        if err_type[:3] == 'MAX':
            return err_con[len(err_con)-1 : len(err_con)-int(err_type.split()[-1])-1 : -1]
        elif err_type[:3] == 'MIN':
            return err_con[:int(err_type.split()[-1])]
        elif err_type == 'D':
            bound = max(abs(err_con[0][0]),abs(err_con[-1][0]))
            count = [0] * 20
            for k in err_con:
                count[ int((k[0] + bound)*10.0/bound) if k[0] < bound else 19 ] += 1
            return [count,[-bound, bound]]
        else:
            raise Exception('Unsupport choice for E_type')
