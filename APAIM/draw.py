#!/usr/bin/env python
'''
Function to get the figure of a configuration.
'''

try:
    import numpy as np
    import matplotlib.pyplot as plt
except:
    print "Configuration visualization unsupported."

from .constant import surface, species_all

def draw(config, repetitions, name):
    '''
    Function to get the figure of a configuration.
        Argument:
            config: an instance of class Config
            repetitions: list of two int, repetitions along x axis and y axis
            name: string for figure name
    '''
    color_def = dict( zip( ['sub'] + list(species_all), ['y','b','r','m','g','c','k']) )
    x_rep,y_rep = repetitions
    if surface == '100':
        fig = plt.figure(figsize=(10,10.0/x_rep*y_rep))
        ax = fig.add_subplot(1,1,1)
        ax.set_xlim([0,x_rep])
        ax.set_ylim([0,y_rep])
        substrate = [[(x,y) for x in sites for y in sites] for sites in ((2,4,6,8),(1,3,5,7))]
        alpha_def = {0:0.3, 1:0.9}
    elif surface == '111':
        fig = plt.figure(figsize=(10,10.0/(x_rep+y_rep*0.5)*y_rep*(3**0.5/2)))
        ax = fig.add_subplot(1,1,1)
        ax.set_xlim([-y_rep*0.5, x_rep])
        ax.set_ylim([0,y_rep*(3**0.5/2)])
        substrate = [[(x,y) for x in sites for y in (1,2,3,4)] for sites in ((3,6,9,12),(2,5,8,11),(1,4,7,10))]
        alpha_def = {0:0.2, 1:0.3, 2:0.9}
    plt.axis('off')
    for num,layer in enumerate(substrate):
        for site in layer:
            for coord in _site_to_coord(site, x_rep, y_rep):
                ax.add_patch(plt.Circle(coord, 0.11, color=color_def['sub'], alpha=alpha_def[num]))
    for ad,site in config.config:
        for coord in _site_to_coord(site, x_rep, y_rep):
            ax.add_patch(plt.Circle(coord, 0.05, color=color_def[ad]))
    if surface == '111':
        ax.add_patch(plt.Polygon([[-y_rep*0.5,0],[0,0],[-y_rep*0.5,y_rep*(3**0.5/2)]],color='w'))
        ax.add_patch(plt.Polygon([[x_rep,0],[x_rep-y_rep*0.5,y_rep*(3**0.5/2)],[x_rep,y_rep*(3**0.5/2)]],color='w'))
    fig.savefig(name+'.png',bbox_inches='tight')


def _site_to_coord(site, x_rep, y_rep):
    coords = []
    for x in range(-1,x_rep+1):
        for y in range(-1,y_rep+1):
            if surface == '100':
                coords.append(((site[1]-1)/8.0+x, (site[0]-1)/8.0+y))
            elif surface == '111':
                def x_off(site_y):
                    if site_y % 3 == 1: return 0
                    elif site_y % 3 == 2: return 1/6.0
                    elif site_y % 3 == 0: return 1/12.0
                coords.append(((site[1]-1)/4.0-((site[0]-1)/12.0+y)/2+x_off(site[0])+x, ((site[0]-1)/12.0+y)*3**0.5/2))
    return coords
