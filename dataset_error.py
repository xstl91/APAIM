#!/usr/bin/env python

from APAIM import *

dataset = raw_input('Directory for dataset: ')
single = load_single()
inter = load_inter()
pred_func = init_predict(single, inter)

print 'type    RMSE Average     Max     Min'
for modify, suffix in ((None,'_none'),('all','_all')):
    configs = load_testdir(dataset, modify)
    for config in configs:
        config.set_energy(pred_func)
    print_energy(configs, 'error' + suffix)
    print ('{:>4}'+'{:>8.4f}'*4).format(suffix[1:],error(configs,'R','A'),error(configs,'A','A'),error(configs,'MAX 1','A')[0][0],error(configs,'MIN 1','A')[0][0])
