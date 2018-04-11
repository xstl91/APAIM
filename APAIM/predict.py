#!/usr/bin/env python

from .constant import *

def init_predict(single_energy, inter_energy, inter_ex = {}):
    '''
    Function to initiate prediction function with certain isolated binding energies and interaction energies.
        Argument:
            single_E: dict of binding types and their energies, e.g. {'N_brg':0.000, ...}
                      can be get through function load_single
            inter_E : dict of interaction, e.g. { 'N-N':{'hh1':0.000, 'hh1f':0.000, ...}, 'N-O':{...}, ...}
                      can be get through function load_inter
            inter_ex: dict of interaction to be excluded, e.g. { 'N-N':['hh3','bb6'], 'N-O':['bb6'] }
                      default to blank dict
        Return:
            function to predict energy for single configuration.
    '''
    inter_E = inter_energy.copy()
    for pair in inter_ex:
        for inter_type in inter_ex[pair]:
            inter_E[pair][inter_type] = None
            for inter in inter_E[pair]:
                if inter[:-1] == inter_type:
                    inter_E[pair][inter] = None

    pred_ex_str = ''
    for pair in species_pairs:
        if pair in inter_ex:
            pred_ex_str += '{:>5}: {}\n'.format(pair,' '.join(inter_ex[pair]))
    
    def pred_energy(config):
        '''
        Function to get predicted energy for a single configuration.
            Argument:
                config: instance of class Config
            Return:
                dict of energy calculated:
                   { 'E_inter_a'  : float, actual interaction energy ,
                     'E_bind_p'   : float, predicted binding energy ,
                     'E_inter_p'  : float, predicted interaction energy ,
                     'E_inter_m'  : float, interaction energy change due to modification
                     'E_err'      : float, error between actual energy and predicted energy ,
                     'pred_ex_str': string, interaction pairs and types excluded for predicting energy }
        '''
        energy_all = { 'pred_ex_str': pred_ex_str }
        E_bind0 = sum( config.binding[x] * single_energy[x] for x in binding_types if config.binding[x] )
        energy_all['E_inter_a'] = ( config.E_bind_a - E_bind0 if config.E_bind_a != None else None )
        energy_all['E_inter_m'] = 0
        for pair_name in config.inter_modify:
            for modify_type in config.inter_modify[pair_name]:
                if inter_E[pair_name][modify_type]:
                    energy_all['E_inter_m'] += config.inter_modify[pair_name][modify_type] * inter_E[pair_name][modify_type]
        energy_all['E_inter_p'] = energy_all['E_inter_m']
        for pair_name in config.inter_pair:
            for inter_type in config.inter_pair[pair_name]:
                if inter_E[pair_name][inter_type]:
                    energy_all['E_inter_p'] += config.inter_pair[pair_name][inter_type][0] * inter_E[pair_name][inter_type]
        energy_all['E_bind_p'] = energy_all['E_inter_p'] + E_bind0
        energy_all['E_err'] = ( energy_all['E_bind_p'] - config.E_bind_a if config.E_bind_a != None else None )
        return energy_all

    pred_energy.pred_ex_str = pred_ex_str

    return pred_energy
