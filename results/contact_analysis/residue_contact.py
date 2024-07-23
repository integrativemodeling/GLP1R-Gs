'''
    Author: Chenxi Wang (chenxi.wang@salilab.org)
    Date: 2022-11-20

    The residue contacts in the GLP-1R-Gs interface is calculated by Pymol.
    
'''

#!/usr/bin/env python
# coding: utf-8

import __main__
import pymol
from pymol import cmd
from pymol import stored
from os import listdir
from os.path import isfile, join


def list_saltbridge(selection,cutoff=4.0,angle=180,mode=1,hb_list_name='saltbridges'):

    cutoff = float(cutoff)
    angle = float(angle)
    mode = float(mode)

    posresList = ["LYS","ARG","HIS","HIP"]
    negresList = ["GLU","ASP","CYM"]
    posatomList = ["NE","NH1","NH2","NZ","ND1","NE2"]
    negatomList = ["SG","OE1","OE2","OD1","OD2"]
    
    selection = selection + " & e. n+o+s"
    
    sb = cmd.find_pairs(selection,selection,mode = mode,cutoff = cutoff,angle = angle)
    sb.sort(key=lambda x:x[0][1])
    
    result = []
    for pairs in sb:
         stored.x = []
         stored.y = []
         cmd.iterate("%s and index %s" % (pairs[0][0],pairs[0][1]), 'stored.x += [chain,resi,resn,name]')
         cmd.iterate("%s and index %s" % (pairs[1][0],pairs[1][1]), 'stored.y += [chain,resi,resn,name]')
         if (stored.x[2] in posresList and stored.x[3] in posatomList and stored.y[2] in negresList and stored.y[3] in negatomList):
              d = cmd.distance(hb_list_name,"%s and index %s" % (pairs[0][0],pairs[0][1]),"%s and index %s" % (pairs[1][0],pairs[1][1]))
              result.append([stored.x,stored.y,float(d)])
         else:
              continue

    return result



def list_hbonds(selection, selection2=None, cutoff=3.5, angle=80, mode=1, hb_list_name='hbonds'):

    cutoff = float(cutoff)
    angle = float(angle)
    mode = float(mode)
    
    selection = selection + " & e. n+o+s & ! resn MET"
    
    if not selection2:
        hb = cmd.find_pairs(selection, selection, mode=mode, cutoff=cutoff, angle=angle)
    else:
        selection2 = selection2 + " & e. n+o & ! resn MET"
        hb = cmd.find_pairs(selection, selection2, mode=mode, cutoff=cutoff, angle=angle)
        
    # Sort the list for easier reading
    hb.sort(key=lambda x: x[0][1])
    
    result = []
    for pairs in hb:
        stored.x = []
        stored.y = []
        cmd.iterate("%s and index %s" % (pairs[0][0], pairs[0][1]), 'stored.x += [chain, resi, resn, name]')
        cmd.iterate("%s and index %s" % (pairs[1][0], pairs[1][1]), 'stored.y += [chain, resi, resn, name]')
        d = cmd.distance(hb_list_name, "%s and index %s" % (pairs[0][0], pairs[0][1]), "%s and index %s" % (pairs[1][0], pairs[1][1]))
        result.append([stored.x, stored.y, float(d)])

    return result



def list_hydrophobic(selection,cutoff=4.0,angle=180,mode=1,hb_list_name='hydrophobic'):

    cutoff = float(cutoff)
    angle = float(angle)
    mode = float(mode)
    hydrophobiclist = ["ALA","VAL","LEU","ILE","PHE","PRO","MET"]
    selection = selection + " & e. c+n+o+s"
    hp = cmd.find_pairs(selection,selection,mode = mode,cutoff = cutoff,angle = angle)
    hp.sort(key=lambda x:x[0][1])
    result = []
    for pairs in hp:
         stored.x = []
         stored.y = []
         cmd.iterate("%s and index %s" % (pairs[0][0],pairs[0][1]), 'stored.x += [chain,resi,resn,name]')
         cmd.iterate("%s and index %s" % (pairs[1][0],pairs[1][1]), 'stored.y += [chain,resi,resn,name]')
         if (stored.x[2] in hydrophobiclist and stored.y[2] in hydrophobiclist):
              d = cmd.distance(hb_list_name,"%s and index %s" % (pairs[0][0],pairs[0][1]),"%s and index %s" % (pairs[1][0],pairs[1][1]))
              result.append([stored.x,stored.y,float(d)])
         else:
              continue

    return result
    


def write_pair(bonds, chain_id1, chain_id2, filepath):

    interface = [item for item in bonds if (item[0][0],item[1][0]) in [(chain_id1,chain_id2),(chain_id2,chain_id1)]] 

    f = open(filepath, 'w')

    if interface != []:
        for item in interface:
            print(f'{item[0][0]}:{item[0][1]}_{item[0][2]}_{item[0][3]} {item[1][0]}:{item[1][1]}_{item[1][2]}_{item[1][3]} {item[2]}', file=f)
    else:
        print('No contact pair found.', file=f)

    f.close()



if __name__ == '__main__':
    
    input_pdb_path = '/public/home/wangchx/glp1r/out_pdbs/'
    inputfiles = [f for f in listdir(input_pdb_path) if isfile(join(input_pdb_path, f)) if '.pdb' in f]
    print(inputfiles)
    output_dir = './'

    chain_id1, chain_id2 = 'A', 'C'

    for pdbname in inputfiles:

        pdbname = pdbname[:-4]
        print(inputfiles.index(pdbname+'.pdb'))

        cmd.load(input_pdb_path + pdbname + '.pdb', pdbname)        
        cmd.h_add()
        cmd.extend("list_hbonds",list_hbonds)
        hb = list_hbonds(pdbname)
        write_pair(hb, chain_id1, chain_id2, output_dir+'hbonds_'+pdbname)
        
        cmd.extend("list_saltbridge",list_saltbridge)
        sb = list_saltbridge(pdbname)
        write_pair(sb, chain_id1, chain_id2, output_dir+'saltbridges_'+pdbname)

        cmd.extend("list_hydrophobic",list_hydrophobic)
        hp = list_hydrophobic(pdbname)
        write_pair(hp, chain_id1, chain_id2, output_dir+'hydrophobic_'+pdbname)

        cmd.delete(pdbname)
        