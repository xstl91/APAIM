#!/usr/bin/env python
'''
Function to count interaction type between two adsorbates.
'''

from .constant import *
from .tool import *


def inter_pair(site1,site2):
    '''
    Function to decide the interaction type between two adsorbates.
        Argument:
            site1: list of name and site for molecule 1, e.g. ['N',[2,2]]
            site2: list of name and site for molecule 2, e.g. ['N',[6,4]]
        Retrun:
            list of interaction type and multiplicity, e.g. ['hh1',1]
    '''
    inter_name = ''.join(( site_type(x)[0] for x in (site1, site2) ))  # e.g. inter_name = 'hh'
    offset = site_to_offset(site1,site2)

    if surface == '100' :
        offset_abs = [ abs(x) for x in offset ]
        for x in offset_abs_type:
            if tuple(offset_abs) in offset_abs_type[x][inter_name[0]]:
                inter_name += x[-1]
                break
        if inter_name in ['bb2','bb4','bb5']:
            if site1[0] % 2 == 0:
                inter_name += ( 'h' if offset_abs in [[0,2],[1,3],[0,4]] else 't' )
            else:
                inter_name += ( 't' if offset_abs in [[0,2],[1,3],[0,4]] else 'h' )
            
    elif surface == '111':
        for x in offset_type:
            if tuple(offset) in offset_type[x][inter_name[0]]:
                inter_name += x[-1]
                break
        if inter_name in ['ff2','hh2','tt2']:
            inter_name += offset_same2[tuple(offset)][inter_name[0]]
            

    if inter_name not in inter_types:
        return [None,None]
    else:
        multi = 1
        for x in multiplicity:
            if inter_name in multiplicity[x]:
                multi = x
        return [inter_name,multi]


if surface == '100':


    def inter_modify(pair_name, inter_type, two_site, config, modify):
        '''
        Function to find interaction types to modify interaction energy.
            Argument:
                pair_name : string for molecule types, e.g. 'O-O'
                inter_type: string for interaction type, e.g. 'hh1'
                two_site  : list for sites, e.g. [[2,2],[2,4]]
                config    : list for configuration, e.g. [['N',[4,6]],['O',[2,2]],['O',[2,4]]]
                modify    : string for modify types to be used, {'all', 'fix', 'recon'}
            Return:
                list of interaction types for modification and corresponding multiplicity,  e.g. [['hh1f',1],['hh1c',1]]
        '''
        modify_list = []
        if modify in ('all','fix'):
            fix = _modify_fix(inter_type, two_site, config)
            if fix:
                modify_list.append(fix)
        if modify in ('all','recon'):
            recon = _modify_recon(pair_name, inter_type, two_site, config)
            if recon:
                modify_list.append(recon)
        return modify_list


    def _modify_fix(inter_type, two_site, config):
        res = [[],[]]
        offset = site_to_offset(two_site[0],two_site[1])
        d = (lambda x: -1 if x == 0 else 1)
        if inter_type in ('bb2h','bb2t','hh1','tt1'):    # e.g. 44-46, res[0]=[42], res[1]=[48]
            for k in (0,1):
                res[k].append( offset_to_site(two_site[k],[ d(k)*x for x in offset ]) )
        elif inter_type in ('bh2','bt2','hb2','tb2'):    # e.g. 43-64, res[0]=[22,31,41,23], res[1]=[85,76,66,84]
            for k in (0,1):
                res[k].append( offset_to_site(two_site[k],[ d(k)*x for x in offset ]) )
                res[k].append( offset_to_site(two_site[k],[ d(k)*2*x if abs(x) == 1 else d(k)*x/2 for x in offset ]) )
                res[k].append( offset_to_site(two_site[k],[ d(k)*2*x if abs(x) == 1 else 0 for x in offset ]) )
                res[k].append( offset_to_site(two_site[k],[ 0 if abs(x) == 1 else d(k)*x for x in offset ]) )
        elif inter_type == 'hh2':                        # e.g. 44-66, res[0]=[22,23,24,32,42], res[1]=[88,87,86,78,68]
            for k in (0,1):
                res[k].append( offset_to_site(two_site[k],[ d(k)*x for x in offset ]) )
                res[k].append( offset_to_site(two_site[k],[ d(k)*offset[0]/2,d(k)*offset[1] ]) )
                res[k].append( offset_to_site(two_site[k],[ 0,d(k)*offset[1] ]) )
                res[k].append( offset_to_site(two_site[k],[ d(k)*offset[0],d(k)*offset[1]/2 ]) )
                res[k].append( offset_to_site(two_site[k],[ d(k)*offset[0],0 ]) )
        sites = [ x[1] for x in config ]
        multi = len([ x for x in res[0] + res[1] if x in sites ])
        if multi:
            return [inter_type+'f', multi]
        else:
            return None


    def _modify_recon(pair_name, inter_type, two_site, config):
        if (pair_name in recon_pair_hh1 and inter_type == 'hh1') or (pair_name in recon_pair_hh2 and inter_type == 'hh2'):
            hh1_offset = ((2,0),(-2,0),(0,2),(0,-2))
            hh2_offset = ((2,2),(-2,2),(2,-2),(-2,-2))
            hb2_offset = ((-2,-1),(-2,1),(-1,-2),(-1,2),(1,-2),(1,2),(2,-1),(2,1))
            hh1n, hh2n, hb2n = [], [], []
            for site in two_site:
                for site_l,off_l in ((hh1n, hh1_offset), (hh2n, hh2_offset), (hb2n, hb2_offset)):
                    for off in off_l:
                        site_a = offset_to_site(site, off)
                        if site_a not in site_l and site_a not in two_site:
                            site_l.append(site_a)
            h_sur = [ x[1] for x in config if x[1] in hh1n + hh2n ]
            b_sur = [ x[1] for x in config if x[1] in hb2n ]
            if inter_type == 'hh1':
                offset = site_to_offset(two_site[0],two_site[1])
                res = [ offset_to_site(two_site[0],[ -x for x in offset ]) ]
                res.append( offset_to_site(two_site[1], offset) )
                if len( [x[1] for x in config if x[1] in res] ) == 2:
                    return [inter_type+'c', -1]
                if len( [x[1] for x in config if x[1] in h_sur] ) >= 4:
                    return [inter_type+'c', -1]
            elif inter_type == 'hh2':
                if len( [x[1] for x in config if x[1] in hh1n] ) >= 2:
                    return [inter_type+'c', -1]
            if len(b_sur) != 0:
                return None
            for o_mol,o_site in zip(pair_name.split('-'), two_site):
                for mol,site in config:
                    if site not in two_site:
                        for pairs,offsets in ((recon_pair_hh1, hh1_offset),(recon_pair_hh2, hh2_offset)):
                            if ('-'.join((o_mol,mol)) in pairs or '-'.join((mol,o_mol)) in pairs) and tuple(site_to_offset(o_site,site)) in offsets:
                                return [inter_type+'c', 1]
            return None


elif surface == '111':


    def inter_modify(pair_name, inter_type, two_site, config, modify):
        '''
        Function to find interaction types to modify interaction energy.
            Argument:
                pair_name : string for molecule types, e.g. 'O-O'
                inter_type: string for interaction type, e.g. 'hh1'
                two_site  : list for sites, e.g. [[2,1],[2,2]]
                config    : list for configuration, e.g. [['N',[6,1]],['O',[2,1]],['O',[2,2]]]
                modify    : string for modify types to be used, {'all', 'restr', 'push'}
            Return:
                list of interaction types for modification and corresponding multiplicity,  e.g. [['tt1r',0.5],['hh1p',1]]
        '''
        modify_list = []
        if modify in ('all','restr'):
            restr = _modify_restr(inter_type, two_site, config)
            if restr:
                modify_list.append(restr)
        if modify in ('all','push'):
            push = _modify_push(inter_type, two_site, config)
            if push:
                modify_list.append(push)
        return modify_list


    def _modify_restr(inter_type, two_site, config):
        res = [[],[]]
        offset = site_to_offset(two_site[0],two_site[1])
        d = (lambda x: -1 if x == 0 else 1)
        if inter_type in ('ff1','hh1'):       # e.g. 62-63, res[0]=[61,41,31], res[1]=[64,44,33]
            s_type = inter_type[0]
            for k in (0,1):
                res[k].append(offset_to_site(two_site[k],[ d(k)*x for x in offset ]))
                tmp1 = set(tuple(offset_to_site(two_site[k],x)) for x in offset_type['dif2'][s_type])
                tmp2 = set(tuple(offset_to_site(res[k][0],x)) for x in offset_type['block'][s_type])
                res[k].extend([ list(x) for x in tmp1 & tmp2 if site_type(x)[0] == 't' ])
                tmp3 = set(tuple(offset_to_site(two_site[k],x)) for x in offset_type['same1'][s_type])
                tmp4 = set(tuple(offset_to_site(res[k][1],x)) for x in offset_type['block']['t'])
                res[k].extend([ list(x) for x in tmp3 & tmp4 if list(x) not in res[k] ])
        elif inter_type == 'tt1':             # e.g. 42-43, res[0]=[41,61,24], res[1]=[44,64,23]
            for k in (0,1):
                res[k].append( offset_to_site(two_site[k],[ d(k)*x for x in offset ]) )
                tmp1 = set(tuple(offset_to_site(two_site[k],x)) for x in offset_type['dif2']['t'])
                tmp2 = set(tuple(offset_to_site(res[k][0],x)) for x in offset_type['block']['t'])
                res[k].extend([ list(x) for x in tmp1 & tmp2 ])
        sites = [ x[1] for x in config ]
        multi = len([ x for x in res[0] + res[1] if x in sites ])
        if multi:
            return [inter_type+'r', multi]
        else:
            return None


    def _modify_push(inter_type, two_site, config):
        if inter_type in ('ff1','hh1'):
            s_type = inter_type[0]
            tmp1 = set(tuple(offset_to_site(two_site[0],x)) for x in offset_type['block'][s_type])
            tmp2 = set(tuple(offset_to_site(two_site[1],x)) for x in offset_type['block'][s_type])
            tmp_site = [ list(x) for x in tmp1 & tmp2 if site_type(x) == 'top' ][0]
            push_tmp = [ offset_to_site(tmp_site,x) for x in offset_type['block']['t']]
            push = [ x for x in push_tmp if site_type(x)[0] == s_type and x not in two_site ][0]
            sites = [ x[1] for x in config ]
            if push in sites:
                return [inter_type+'p',1]
            else:
                return None
        else:
            return None
