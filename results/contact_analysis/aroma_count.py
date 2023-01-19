'''
    Author: Chenxi Wang (chenxi.wang@salilab.org)
    Date: 2022-12-05

    This script is used to identify aromatic hydrophobic residue contact.
    The distance between the atoms in one residue and the center of the benzene ring in another aromatic residue is smaller than or equal to 4.0 Å, 
    or the distance between the centers of two benzene ring is smaller than or equal to 5.0 Å.
    
'''

from pymol import math
from os import listdir
from os.path import isfile, join


def get_dict_coord_from_pdb(f):

    ''' return a dict of minimal distances of heavy atoms for each pair of residues in a standard pdb file'''
    
    atom_lines = [l for l in f if l[0:6] == 'ATOM  ']
    dict_coord = {} # dict to store coordinates. dict_coord[res][atom] = (x,y,z,occ)
    atomnum_2_name = {} # map atom number to atom name, in order to find N, CA, C, O
    for line in atom_lines:
        # retrive info from each atom line
        atom_num = int(line[6:11].strip())
        atom_name = line[12:16].replace(' ', '_')
        res_name = line[17:20]
        res_num = int(line[22:26].strip())
        chain_id = line[21:22]
        x = float(line[30:38].strip())
        y = float(line[38:46].strip())
        z = float(line[46:54].strip())
        occ = float(line[54:60].strip())
        res = chain_id + ':' + str(res_num) + '_' + res_name 
        atomnum_2_name[atom_num] = atom_name
        if res_num <= 0:
            continue
        if res not in dict_coord:
            dict_coord[res] = {}
        dict_coord[res][atom_name] = (x, y, z, occ)
                
    return dict_coord



def get_center_of_atoms(atoms):

    x=0
    y=0
    z=0
    
    for atm in atoms:
        x+=atm[0]
        y+=atm[1]
        z+=atm[2]
    
    numatm=float(len(atoms))
    x=x/numatm
    y=y/numatm
    z=z/numatm

    return [x,y,z]



def get_dist(cen1,cen2):

    dx=cen1[0]-cen2[0]
    dy=cen1[1]-cen2[1]
    dz=cen1[2]-cen2[2]
    dist=dx*dx+dy*dy+dz*dz

    return math.sqrt(dist)


aromaticlist = ["PHE", "TYR", "TRP"]
aromatic_dict = {'PHE': ('_CG_', '_CD1', '_CD2', '_CE1', '_CE2', '_CZ_'), 
                 'TYR': ('_CG_', '_CD1', '_CD2', '_CE1', '_CE2', '_CZ_'),
                 'TRP': ('_CD2','_CE2', '_CE3', '_CZ2', '_CZ3', '_CH2')}


input_pdb_path = '/public/home/wangchx/glp1r/out_pdbs/'
inputfiles = [f for f in listdir(input_pdb_path) if isfile(join(input_pdb_path, f))]

for file in inputfiles:
    if 'cl0' in file:
        fout = open('./hydrophobic/aromatic/'+file[:-4], 'w')
        finput = open(input_pdb_path+file).read().splitlines()
        coord_dict = get_dict_coord_from_pdb(finput)

        for ires in list(coord_dict.keys()):
            if 'A:' in ires:
                if ires.split('_')[1] in aromaticlist:
                    for jres in list(coord_dict.keys()):
                        if 'C:' in jres or 'D:' in jres:
                            atom_in_ires = coord_dict[ires].keys()
                            atom_in_jres = coord_dict[jres].keys()
                            aroma_atom_coord = []
                            for iatom in atom_in_ires:
                                if iatom in aromatic_dict[ires.split('_')[1]]:
                                    (ix, iy, iz, iocc) = coord_dict[ires][iatom]
                                    aroma_atom_coord += [(ix, iy, iz, iocc)]
                            if len(aroma_atom_coord) == 0:
                                # print(ires, atom_in_ires)
                                pass
                            else:
                                aroma_center = get_center_of_atoms(aroma_atom_coord)
                                for jatom in atom_in_jres:                  
                                    (jx, jy, jz, jocc) = coord_dict[jres][jatom]
                                    d = get_dist(aroma_center, [jx, jy, jz])
                                    if 0 < d <= 4:
                                        print(*[ires, jres, jatom, d], file=fout)

                            ####### additional pi-pi count #####
                            if jres.split('_')[1] in aromaticlist:
                                iaroma_atom_coord = []
                                for iatom in atom_in_ires:
                                    if iatom in aromatic_dict[ires.split('_')[1]]:
                                        (ix, iy, iz, iocc) = coord_dict[ires][iatom]
                                        iaroma_atom_coord += [(ix, iy, iz, iocc)]
                                jaroma_atom_coord = []
                                for jatom in atom_in_jres:
                                    if jatom in aromatic_dict[jres.split('_')[1]]:
                                        (jx, jy, jz, jocc) = coord_dict[jres][jatom]
                                        jaroma_atom_coord += [(jx, jy, jz, jocc)]
                                if len(iaroma_atom_coord) == 0 or len(jaroma_atom_coord) == 0 :
                                    # print(ires, atom_in_ires)
                                    pass
                                else:
                                    iaroma_center = get_center_of_atoms(iaroma_atom_coord)
                                    jaroma_center = get_center_of_atoms(jaroma_atom_coord)
                                    d = get_dist(iaroma_center, jaroma_center)
                                    if 0 < d <= 5:
                                        print(*[ires, jres, 'pi-pi', d], file=fout)


            elif 'C:' in ires or 'D:' in ires:
                if ires.split('_')[1] in aromaticlist:
                    for jres in list(coord_dict.keys()):
                        if 'A:' in jres:
                            atom_in_ires = coord_dict[ires].keys()
                            atom_in_jres = coord_dict[jres].keys()
                            aroma_atom_coord = []
                            for iatom in atom_in_ires:
                                if iatom in aromatic_dict[ires.split('_')[1]]:
                                    (ix, iy, iz, iocc) = coord_dict[ires][iatom]
                                    aroma_atom_coord += [(ix, iy, iz, iocc)]
                            if len(aroma_atom_coord) == 0:
                                # print(ires, atom_in_ires)
                                pass
                            else:
                                aroma_center = get_center_of_atoms(aroma_atom_coord)
                                for jatom in atom_in_jres:                  
                                    (jx, jy, jz, jocc) = coord_dict[jres][jatom]
                                    d = get_dist(aroma_center, [jx, jy, jz])
                                    if 0 < d <= 4:
                                        print(*[ires, jres, jatom, d], file=fout)
                                        
                            ####### additional pi-pi count #####
                            if jres.split('_')[1] in aromaticlist:
                                iaroma_atom_coord = []
                                for iatom in atom_in_ires:
                                    if iatom in aromatic_dict[ires.split('_')[1]]:
                                        (ix, iy, iz, iocc) = coord_dict[ires][iatom]
                                        iaroma_atom_coord += [(ix, iy, iz, iocc)]
                                jaroma_atom_coord = []
                                for jatom in atom_in_jres:
                                    if jatom in aromatic_dict[jres.split('_')[1]]:
                                        (jx, jy, jz, jocc) = coord_dict[jres][jatom]
                                        jaroma_atom_coord += [(jx, jy, jz, jocc)]
                                if len(iaroma_atom_coord) == 0 or len(jaroma_atom_coord) == 0 :
                                    # print(ires, atom_in_ires)
                                    pass
                                else:
                                    iaroma_center = get_center_of_atoms(iaroma_atom_coord)
                                    jaroma_center = get_center_of_atoms(jaroma_atom_coord)
                                    d = get_dist(iaroma_center, jaroma_center)
                                    if 0 < d <= 5:
                                        print(*[ires, jres, 'pi-pi', d], file=fout)

        fout.close()
