#!/usr/bin/env python
from math import sqrt

from APAIM import *

dataset = raw_input('Directory for dataset: ')
modify = 'all'

if surface == '100':
    cutoffs = (([], 3.0/2),
               (['ht3','th3','bb6'], sqrt(34)/4),
               (['bt5','tb5','bh5','hb5'], sqrt(2)),
               (['tt3','bb5h','bb5t','hh3'], sqrt(26)/4),
               (['bt4','tb4','bh4','hb4'], sqrt(5)/2),
               (['ht2','th2','bb4h','bb4t'], sqrt(18)/4),
               (['bt3','tb3','bh3','hb3'], 1),
               (['tt2','bb3','hh2'], sqrt(10)/4),
               (['bt2','tb2','bh2','hb2'], sqrt(2)/2),
               (['tt1','bb2h','bb2t','hh1'], 0))
elif surface == '111':
    cutoffs = (([], sqrt(24)/3),
               (['ft5','tf5','ht5','th5','fh5','hf5'], sqrt(78)/6),
               (['ft4','tf4','ht4','th4','fh4','hf4'], sqrt(2)),
               (['tt3','ff3','hh3'], sqrt(6)/2),
               (['tt2f','tt2h','ff2t','ff2h','hh2t','hh2f'], sqrt(42)/6),
               (['ft3','tf3','ht3','th3','fh3','hf3'], sqrt(6)/3),
               (['ft2','tf2','ht2','th2','fh2','hf2'], sqrt(2)/2),
               (['tt1','ff1','hh1'], 0))

single = load_single()
inter = load_inter()
print ('{:>8}'*5).format('distance','RMSE','Average','Max','Min')
configs = load_testdir(dataset, modify)
ex_list = []
for ex_a,dis in cutoffs:
    ex_list.extend(ex_a)
    inter_ex = { x:ex_list for x in species_pairs }
    pred_func = init_predict(single, inter, inter_ex)
    for config in configs:
        config.set_energy(pred_func)
    print ('{:>8.4f}'*5).format(dis,error(configs,'R','A'),error(configs,'A','A'),error(configs,'MAX 1','A')[0][0],error(configs,'MIN 1','A')[0][0])
print '* Unit for distance is the lattice constant of bulk crystal.'
