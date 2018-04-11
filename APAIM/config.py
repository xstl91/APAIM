#!/usr/bin/env python
'''
Class for single configuration.
'''

from .constant import *
from .count_inter import *
from .tool import *
from .draw import draw

class Config(object):
    '''
    Data:
        self.species_str  : string, sorted in standard ortder, e.g. '2N-1O-1NO
        self.config_srt   : string, corresponding to self.species_str, e.g. '24-36-44-76'
        self.species      : dict, e.g. {'N':2, 'O':1, 'NO':1, 'CO':0}
        self.config       : list, e.g. [['N',[2,4]],['N',[3,6]],['O',[4,4]],['NO',[7,6]]]
        self.binding      : dict, e.g. {'N_brg':1, 'N_hol':1, 'O_brg':0, ...}
        self.inter_pair   : dict, e.g. { 'N-N':{'hh1':[2,[[[2,4],[2,6]],[[2,4],[4,4]]], 'bb2h':[1,...], ...}, 'N-O':{...}, ...}
        self.inter_modify : dict, e.g. { 'N-N':{'hh1f':0.5, 'hh1c':1, ...}, 'N-O':{...}, ...}
        self.modify_str   : string, modification type used
        self.E_bind_a     : float, actual binding energy

    Data set by self.set_energy(pred_fun)
        self.E_inter_a    : float, actual interaction energy
        self.E_bind_p     : float, predicted binding energy
        self.E_inter_p    : float, predicted interaction energy
        self.E_inter_m    : float, interaction energy change due to modification
        self.E_err        : float, error between actual energy and predicted energy
        self.pred_ex_str  : string, interaction pairs and types excluded for predicting energy
    '''


    def __init__(self, species_str, config_str, E_actual = None, modify = None):
        '''
        __init__(self, species_str, config_str, E_actual = None, )
            species_str: str, number and adsorbates separated by '-', e.g. '2N-1O-1NO'
            config_str : str, sites separated by '-' correspond to species_str, e.g. '24-36-44-76'
            E_actual   : float, actual binding energy, e.g. -7.325
            modify     : string for modification type to be used, refer to inter_modify
        '''
        self.E_bind_a = E_actual
        self._init_species(species_str)
        self._init_config(species_str, config_str)
        self._init_binding()
        self._init_pair()
        self._init_modify(modify)
        self.modify_str = modify


    def _init_species(self, species_str):
        self.species = { x:0 for x in species_all }
        for mol,num in get_species(species_str):
            self.species[mol] = num

        self.species_str = ''
        for mol in species_all:
            if self.species[mol]:
                self.species_str += (str(self.species[mol]) + mol + '-')
        self.species_str = self.species_str[:-1]


    def _init_config(self, species_str, config_str):
        sites = sort_site(species_str, config_str.split('-'))
        self.config_str = '-'.join(sites)
        self.config = []
        for mol in species_all:
            for i in range(self.species[mol]):
                self.config.append([mol, site_convert(sites.pop(0))])


    def _init_binding(self):
        self.binding = { x:0 for x in binding_types }
        for site in self.config:
            self.binding[ site[0] + '_' + site_type(site[1]) ] += 1


    def show_binding(self):
        '''Function for printing binding information in a readable way.'''
        for x in binding_types:
            if self.binding[x]:
                print '{:>6} {:>2d}'.format(x, self.binding[x])


    def _init_pair(self):
        self.inter_pair = {}
        for i in range(len(self.config)-1):
            for j in range(i+1,len(self.config)):
                pair_name = '-'.join(( self.config[x][0] for x in (i,j) ))
                inter_type,multi = inter_pair(self.config[i][1], self.config[j][1])
                if inter_type:
                    sites = [ self.config[x][1] for x in (i,j) ]
                    if pair_name not in self.inter_pair:
                        self.inter_pair[pair_name] = {}
                    if self.config[i][0] == self.config[j][0] and inter_type in same_exclude:
                        if inter_type == 'bb4t':
                            inter_type = 'bb4h'
                        elif inter_type == 'ff2h':
                            inter_type = 'ff2t'
                        elif inter_type == 'hh2t':
                            inter_type = 'hh2f'
                        elif inter_type == 'tt2f':
                            inter_type = 'tt2h'
                        else:
                            inter_type = ''.join((inter_type[1],inter_type[0],inter_type[2]))
                            sites = [sites[1], sites[0]]
                    if inter_type not in self.inter_pair[pair_name]:
                        self.inter_pair[pair_name][inter_type] = [0,[]]
                    self.inter_pair[pair_name][inter_type][0] += multi
                    self.inter_pair[pair_name][inter_type][1].append(sites)


    def show_inter_pair(self):
        '''Function for printing pairwise interaction types in a readable way.'''
        for x in species_pairs:
            if x in self.inter_pair:
                for y in inter_types:
                    if y in self.inter_pair[x]:
                        print '{:>5} {:>5} {:>5}'.format(x,y,self.inter_pair[x][y][0])


    def _init_modify(self, modify):
        self.inter_modify = {}
        for pair_name in self.inter_pair:
            for inter_type in self.inter_pair[pair_name]:
                for two_site in self.inter_pair[pair_name][inter_type][1]:
                    for modify_type, multi in inter_modify(pair_name, inter_type, two_site, self.config, modify):
                        if pair_name not in self.inter_modify:
                            self.inter_modify[pair_name] = {}
                        if modify_type not in self.inter_modify[pair_name]:
                            self.inter_modify[pair_name][modify_type] = 0
                        self.inter_modify[pair_name][modify_type] += multi


    def show_inter_modify(self):
        '''Function for printing modification interaction types in a readable way.'''
        for x in species_pairs:
            if x in self.inter_modify:
                for y in self.inter_modify[x]:
                    print '{:>5} {:>5} {:>5}'.format(x,y,self.inter_modify[x][y])


    def set_energy(self, pred_func):
        '''
        Function to set E_inter_a, E_bind_p, E_inter_p, E_inter_m, E_err, pred_ex_str for instance itself.
            Argument:
                pred_func: function to predict energy, can be get using function init_predict
                           can be get using function init_predict
        '''
        self.__dict__.update(pred_func(self))


    def get_equal(self):
        '''
        Function to get equal configurations.
            Return:
                list of list of sites, e.g. [['12','44'],['14','46'],...]
        '''
        equal_config = [self.config_str.split('-')]
        for sym in sym_operator:
            equal_site = []
            for site in equal_config[-1]:
                equal_site.append(symmetry_site[site][sym])
            equal_config.append(equal_site)
        for config in equal_config[:]:
            for x in trans_operator[0]:
                for y in trans_operator[1]:
                    if (x,y) != (0,0):
                        config_trans = []
                        for site in config:
                            config_trans.append(site_convert(offset_to_site(site_convert(site),[x,y])))
                        equal_config.append(config_trans)
        unieq_config = sorted(set( tuple(sort_site(self.species_str,x)) for x in equal_config ))
        return [ list(x) for x in unieq_config ]


    def draw(self, repetitions, name):
        '''
        Function to get the figure of the configuration.
            Argument:
                repetitions: list of two int, repetitions along x axis and y axis
                name: string for figure name
        '''
        try:
            draw(self, repetitions, name)
        except:
            print "Configuration visualization unsupported."


    def __eq__(self,another):
        '''Return True if two configurations are the same considering the symmetry.'''
        if not isinstance(another, Config): raise(TypeError('An instance of class Config is needed.'))
        if self.species_str == another.species_str and another.config_str.split('-') in self.get_equal():
            return True
        else:
            return False


    def __ne__(self,another):
        '''Return True if two configurations are not the same considering the symmetry.'''
        return not self == another

