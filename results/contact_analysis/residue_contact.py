'''
    Author: Chenxi Wang (chenxi.wang@salilab.org)
    Date: 2022-11-20

    The residue contacts in the GLP-1R-Gs interface is calculated by Pymol.
    
'''

#!/usr/bin/env python
# coding: utf-8

from pymol import cmd
from pymol import stored
from os import listdir
from os.path import isfile, join


def list_hbonds(selection,selection2=None,cutoff=3.5,angle=80,mode=1,hb_list_name='hbonds'):

    cutoff = float(cutoff)
    angle = float(angle)
    mode = float(mode)
    
    selection = selection + " & e. n+o+s"
    
    if not selection2:
        hb = cmd.find_pairs(selection,selection,mode = mode,cutoff = cutoff,angle = angle)
    else:
        selection2 = selection2 + " & e. n+o"
        hb = cmd.find_pairs(selection,selection2,mode = mode,cutoff = cutoff,angle = angle)
        
    # sort the list for easier reading
    hb.sort(key=lambda x:x[0][1])
    
    result = []
    for pairs in hb:
        stored.x = []
        stored.y = []
        cmd.iterate("%s and index %s" % (pairs[0][0],pairs[0][1]), 'stored.x += [chain,resi,resn,name]')
        cmd.iterate("%s and index %s" % (pairs[1][0],pairs[1][1]), 'stored.y += [chain,resi,resn,name]')
        d = cmd.distance(hb_list_name,"%s and index %s" % (pairs[0][0],pairs[0][1]),"%s and index %s" % (pairs[1][0],pairs[1][1]))
        result.append([stored.x,stored.y,float(d)])

    return result
    
    


def list_saltbridge(selection,selection2=None,cutoff=4.0,angle=180,mode=1,hb_list_name='saltbridges'):

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




def list_hydrophobic(selection,selection2=None,cutoff=4.0,angle=180,mode=1,hb_list_name='hydrophobic'):

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
    


def select_interface_pair(bonds, filepath):

    R_A, R_B = [], []
    
    for item in bonds:
        pairs = (item[0][0],item[1][0])
        # change the index to chain index
        if pairs in [('A','C'),('C','A')]: R_A += [item] 
        if pairs in [('A','D'),('D','A')]: R_B += [item] 

    file = open(filepath, 'w')
    if R_A != []:
        for item in R_A:
            print('R_A', item[0][0], item[0][1], item[0][2], item[0][3], item[1][0], item[1][1], item[1][2], item[1][3], item[2], sep=' ', file=file)
    if R_B != []:
        for item in R_B:
            print('R_B', item[0][0], item[0][1], item[0][2], item[0][3], item[1][0], item[1][1], item[1][2], item[1][3], item[2], sep=' ', file=file)
    file.close()



########################### main ###########################

input_pdb_path = '/public/home/wangchx/glp1r/out_pdbs/'
inputfiles = [f for f in listdir(input_pdb_path) if isfile(join(input_pdb_path, f)) if '.pdb' in f]
output_dir = './'

for pdbname in inputfiles:
    if 'cl0' in pdbname:
        pdbname = pdbname[:-4]
        print(inputfiles.index(pdbname+'.pdb'))

        cmd.load(input_pdb_path + pdbname + '.pdb', pdbname)        
        cmd.h_add()
        cmd.extend("list_hbonds",list_hbonds)
        hb = list_hbonds(pdbname)
        select_interface_pair(hb, output_dir+'hbonds/hbonds_'+pdbname, 'hbond')
        
        cmd.extend("list_saltbridge",list_saltbridge)
        sb = list_saltbridge(pdbname)
        select_interface_pair(sb, output_dir+'saltbridges/saltbridges_'+pdbname, 'saltbridge')

        cmd.extend("list_hydrophobic",list_hydrophobic)
        hp = list_hydrophobic(pdbname)
        select_interface_pair(hp, output_dir+'hydrophobic/distance/hydrophobic_'+pdbname, 'hydrophobic')

        cmd.delete(pdbname)
        
