#!/usr/bin/env python

import random
import os

# Set parameters used
num_ad = {'N': input('Number for N atom(s): '),
          'O': input('Number for O atom(s): '),
          'NO':input('Number for NO molecule(s): '),
          'CO':input('Number for CO molecule(s): ')}
if not 0 < sum(num_ad.values()) <= 16:
    raise Exception('Too many adsorbates.')
from APAIM import *
single = load_single()
inter = load_inter()
modify = 'all'
pred_func = init_predict(single, inter)


from APAIM.tool import offset_to_site

# Get species_str
species_str = ''
for mol in ['N','O','NO','CO']:
    if num_ad[mol] != 0: species_str += str(num_ad[mol]) + mol + '-'
species_str = species_str[:-1]
print
print species_str

# Define blocked offsets
def block_offsets(site):
    if surface == '100':
        return [[-1,-1],[-1,0],[-1,1],[0,-1],[0,1],[1,-1],[1,0],[1,1]]
    elif site[0] % 3 == 0 :
        return [[1,0],[2,0],[1,1],[-1,0],[-2,0],[-1,-1]]
    elif site[0] % 3 == 1 :
        return [[1,-1],[2,0],[1,0],[-1,0],[-2,-1],[-1,-1]]
    elif site[0] % 3 == 2 :
        return [[1,0],[2,1],[1,1],[-1,1],[-2,0],[-1,0]]
block_num = 8 if surface == '100' else 6

# Define offsets for nearest same sites
def same_offsets(site):
    if surface == '100':
        return [[-2,0],[0,-2],[2,0],[0,2]]
    else:
        return [[0,-1],[0,1],[-3,-1],[-3,0],[3,0],[3,1]]

# Define nearest sites
def block_sites(site):
    return [ offset_to_site(site, x) for x in block_offsets(site) ]
def same_sites(site):
    return [ offset_to_site(site, x) for x in same_offsets(site) ]

# Get config_str
def config_str(sites):
    config_str = ''
    for site_number in sites:
        site_str = str(site_number[0]*10 + site_number[1])
        if len(site_str) == 3:
            if site_str[1] == '0': site_str = 'a' + site_str[2]
            elif site_str[1] == '1': site_str = 'b' + site_str[2]
            elif site_str[1] == '2': site_str = 'c' + site_str[2]
        config_str += site_str + '-'
    config_str = config_str[:-1]
    return config_str

# Get energy for certain configuration:
def energy(sites):
    config = Config(species_str, config_str(sites), modify = modify)
    config.set_energy(pred_func)
    return config.E_bind_p

# Set initial configuration
if surface == '100':
    sites_origin = [ [x,y] for x in range(1,9) for y in range(1,9) ]
else:
    sites_origin = [ [x,y] for x in range(1,13) for y in range(1,5) ]
initial = True
while initial:
    initial = False
    sites_choose = []
    sites_block = { tuple(x):0 for x in sites_origin }
    sites_left = sites_origin[:]
    for mol in ['N','O','NO','CO']:
        if initial: break
        for i in range(num_ad[mol]):
            if surface == '100' and mol in ['N','O']:
                sites_from = [ x for x in sites_left if x[0] % 2 != 1 or x[1] % 2 != 1 ]
            elif surface == '111' and mol in ['N','O']:
                sites_from = [ x for x in sites_left if x[0] % 3 != 1 ]
            else:
                sites_from = sites_left
            if not sites_from:
                initial = True
                break
            site_choose = random.choice( sites_from )
            sites_left.remove(site_choose)
            for site_block in block_sites(site_choose):
                sites_block[tuple(site_block)] += 1
                if site_block in sites_left:
                    sites_left.remove(site_block)
            sites_choose.append(site_choose)

# Structure optimization
next_round = True
current_energy = energy(sites_choose)
print '{:.3f}  '.format(current_energy) + config_str(sites_choose)
while next_round:
    next_round = False
    for i,site in enumerate(sites_choose):
        site_next = None
        for j,site_to in enumerate(block_sites(site) + same_sites(site)):
            if ( i < num_ad['N'] + num_ad['O']
                 and ( site_to[0]%2 + site_to[1]%2 == 2 if surface == '100' else site_to[0]%3 == 1 )):
                continue
            if (( j < block_num and sites_block[tuple(site_to)] == 1 )
                 or ( j >= block_num and sites_block[tuple(site_to)] == 0 and site_to not in sites_choose )) :
                sites_to = sites_choose[:]
                sites_to[i] = site_to
                if energy(sites_to) < current_energy:
                    site_next = site_to
                    sites_choose = sites_to
                    current_energy = energy(sites_choose)
        if site_next:
            for site_block in block_sites(site):
                sites_block[tuple(site_block)] -= 1
            for site_block in block_sites(site_next):
                sites_block[tuple(site_block)] += 1
            print '{:.3f}  '.format(current_energy) + config_str(sites_choose)
            next_round = True

# Finalization
config = Config(species_str, config_str(sites_choose), modify = modify)
print '-'*20
print '{:.3f}  '.format(current_energy) + '-'.join(config.get_equal()[0])
