'''
    Author: Chenxi Wang (chenxi.wang@salilab.org)
    Date: 2022-11-20

    The interface area was calculated by the program FreeSASA 2.0, 
    using the Sharke-Rupley algorithm with a probe radius of 1.4 Å, 
    and in accordance with the definition of buried surface area in PDBePISA (https://www.ebi.ac.uk/pdbe/pisa/)  
    
'''

import freesasa
from os import listdir
from os.path import isfile, join


def compute_area(f1, f2, f3):

    param = freesasa.Parameters()
    param.setProbeRadius(1.4)

    complex = freesasa.Structure(f1)
    result = freesasa.calc(complex, param)
    complex_sasa = result.totalArea()
    # print('complex', complex_sasa)

    Gab = freesasa.Structure(f2)
    result = freesasa.calc(Gab, param)
    Gab_SASA = result.totalArea()
    # print(Gab_SASA)

    Gr = freesasa.Structure(f3)
    result = freesasa.calc(Gr, param)
    Gr_SASA = result.totalArea()
    # print(Gr_SASA)

    interface_SASA = (Gab_SASA + Gr_SASA - complex_sasa)/2
    
    return interface_SASA



def generate_file(pdblines, ICL_idx):

    f_ICL_name = './temp_ICL.pdb'
    f_ICL_AB_name = './temp_ICL_AB.pdb'

    f_ICL= open(f_ICL_name, 'w')
    for line in pdblines:
        res_num = int(line[22:26].strip())
        if line[21:22] == receptor_domain_index and ICL_idx[0] <= res_num <= ICL_idx[1]:
            print(line, file=f_ICL)
    f_ICL.close()

    f_ICL_AB= open(f_ICL_AB_name, 'w')
    for line in pdblines:
        res_num = int(line[22:26].strip())
        if line[21:22] == receptor_domain_index and ICL_idx[0] <= res_num <= ICL_idx[1]:
            print(line, file=f_ICL_AB)
        elif line[21:22] in (alpha_domain_index, beta_domain_index):
            print(line, file=f_ICL_AB)
    f_ICL_AB.close()

    return f_ICL_name, f_ICL_AB_name
    



pdbpath = '/public/home/wangchx/glp1r/out_pdbs/'
pdbfiles = [f for f in listdir(pdbpath) if isfile(join(pdbpath, f)) and '.pdb' in f]

receptor_domain_index = 'A'
alpha_domain_index = 'C'
beta_domain_index = 'D'   

fout = open('./interface_area.csv','w')

for pdbfile in pdbfiles:
    if 'cl0' in pdbfile:
        pdb = open(pdbpath+pdbfile).read().splitlines()
        atom_lines = [l for l in pdb if l[0:6] == 'ATOM  ']

        ########### total buried interface area ###########
        fRAB_name = './temp_RAB.pdb'
        fRAB = open(fRAB_name,'w')
        for line in atom_lines:
            if line[21:22] in (receptor_domain_index, alpha_domain_index, beta_domain_index):
                print(line, file=fRAB)
        fRAB.close()

        fAB_name = './temp_AB.pdb'
        fAB = open(fAB_name,'w')
        for line in atom_lines:
            if line[21:22] in (alpha_domain_index, beta_domain_index):
                print(line, file=fAB)
        fAB.close()

        fR_name = './temp_R.pdb'
        fR = open(fR_name,'w')
        for line in atom_lines:
            if line[21:22] == receptor_domain_index:
                print(line, file=fR)
        fR.close()

        total = compute_area(fRAB_name, fAB_name, fR_name)

        ########### ICL1 interface area ###########
        
        fICL1_name, fICL1_AB_name = generate_file(atom_lines, ICL_idx=[146, 151])
        ICL1 = compute_area(fICL1_AB_name, fAB_name, fICL1_name)

        ########### ICL2 interface area ###########
        fICL2_name, fICL2_AB_name = generate_file(atom_lines, ICL_idx=[234, 239])
        ICL2 = compute_area(fICL2_AB_name, fAB_name, fICL2_name)

        ########### ICL3 interface area ###########
        fICL3_name, fICL3_AB_name = generate_file(atom_lines, ICL_idx=[316, 321])
        ICL3 = compute_area(fICL3_AB_name, fAB_name, fICL3_name)

        print(pdbfile[:-4], total, ICL1, ICL2, ICL3, sep=',', file=fout)

fout.close()
