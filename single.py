#!/usr/bin/env python

from APAIM import *
species_str = raw_input("Species (e.g. 2N-2O): ")
config_str = raw_input("Configuration (e.g. 22-26-62-68): ")
single = load_single()
inter = load_inter()
modify = 'all'
pred_func = init_predict(single, inter)
config = Config(species_str, config_str, modify = modify)
config.set_energy(pred_func)
print
print config.E_bind_p
