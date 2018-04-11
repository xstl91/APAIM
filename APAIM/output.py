#!/usr/bin/env python

from .constant import *

def print_energy(config_list, output_file = 'summary_E'):
    '''
    Function to print energy informations to a file.
        Argument:
            config_list: list of instances of class Config
            output_file: string for the route and file name
    '''
    with open(output_file,'w') as f:
        len_sp_str = max( max( len(x.species_str) for x in config_list ), 7)
        len_co_str = max( max( len(x.config_str) for x in config_list ), 13)
        f.write('{0:>{1}} {2:>{3}}  E_bind_a  E_bind_p E_inter_p E_inter_m     E_err E_err_ave\n'.format('species',len_sp_str,'configuration',len_co_str))
        for config in config_list:
            f.write('{2.species_str:>{0}} {2.config_str:>{1}} {2.E_bind_a:>9.3f} {2.E_bind_p:>9.3f} {2.E_inter_p:>9.3f} {2.E_inter_m:>9.3f} {2.E_err:>9.3f}'.format(len_sp_str,len_co_str,config))
            f.write(' {:>9.3f}\n'.format(config.E_err/sum( config.species.values() )))


def print_error(error_list, name_1c, output_file = 'summary_e'):
    '''
    Function to print error informations to a file.
        Argument:
            error_list : list of error information list, which contains name, average error, RMSE and distribution of error.
            name_1c    : name for first colomn
            output_file: string for the route and file name
    '''
    name_len = max(len(name_1c), max(len(x[0]) for x in error_list))
    with open(output_file,'w') as f:
        f.write('{0:{4}}  {1:9}  {2:8}  {3}\n'.format(name_1c,'error_ave','RMSE','distribution',name_len))
        for item in error_list:
            f.write('{0[0]:{1}}{0[1]:>11.2e}{0[2]:>10.2e}'.format(item,name_len))
            f.write('  [{:>9.2e},{:>9.2e}] {}\n'.format(item[3][1][0], item[3][1][1], str(item[3][0])))


def inter_example(inter = 'all'):
    '''
    Function to print example configuration for each interaction type.
        Argument:
            inter: string for interaction type name
                   default to 'all', which stands for all interaction types.
    ''' 
    if inter == 'all':
        for k in inter_types:
            print '{:<4} {}'.format(k,inter_examp[k][0])
    else:
        if inter in inter_examp:
            print '{:<4} {}'.format(inter,inter_examp[inter][0])
