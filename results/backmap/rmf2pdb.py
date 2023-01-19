import IMP
import IMP.atom
import IMP.core
import IMP.display
import IMP.rmf
import IMP.pmi
import IMP.pmi.analysis
import IMP.pmi.output
import RMF
import glob
import re
import numpy as np


color_dict = {'GLP1R':(0,0,1),
              'GLP1':(1,0,0),
              'GA':(0.627,0.125,0.941),
              'GB':(1,0.498,0),
              'GG':(0,1,1),
              'NB':(1,1,0)}


def color_rmf(rmf_file):
    rmf_file = rmf_file.split('.rmf3')[0]+'.rmf3'
    m = IMP.Model()
    prots = IMP.pmi.analysis.get_hiers_from_rmf(m, 0, rmf_file)
    print(rmf_file, prots)
    prots = prots[0]
    for p,col in color_dict.items():
        s=IMP.atom.Selection(prots,
                            molecule=p,
                            resolution=IMP.atom.ALL_RESOLUTIONS)
        psel=s.get_selected_particles()
        color=IMP.display.Color(col[0],col[1],col[2])
        for part in psel:
            IMP.display.Colored(part).set_color(color)

    o=IMP.pmi.output.Output()
    o.init_rmf(rmf_file,[prots])
    o.write_rmf(rmf_file)
    o.close_rmf(rmf_file)
    return rmf_file

    

def write_pdb(rmf_file):
    pdb_file = rmf_file.split('.rmf3')[0]+'.pdb'
    mm = IMP.Model()
    prots = IMP.pmi.analysis.get_hiers_from_rmf(mm, 0, rmf_file)
    print(rmf_file, prots)
    prot=prots[0]
    for p,col in color_dict.items():
        s=IMP.atom.Selection(prots,
                            molecule=p,
                            resolution=IMP.atom.ALL_RESOLUTIONS)
        psel=s.get_selected_particles()
        color=IMP.display.Color(col[0],col[1],col[2])
        for part in psel:
            IMP.display.Colored(part).set_color(color)

    if not prots:
        raise ValueError("Cannot read hiearchy from rmf")
    o = IMP.pmi.output.Output()
    o.init_pdb(pdb_file, prot)
    o.write_pdb(pdb_file)
    del o, mm

#################################
########### MAIN ################
#################################
top_dir =  '../free_noGDP/PMI_analysis_singlenode_52cores_output'
clusters = ['GSMs_cl0']
samples = ['sample_A']
for cluster in clusters:
    for sample in samples:
        rmfdir = top_dir+'/'+cluster+'/'+sample+'/'
        rmf_files = glob.glob(rmfdir + '*.rmf3')
        for rmf_file in rmf_files:
            color_rmf(rmf_file)
            write_pdb(color_rmf(rmf_file))





