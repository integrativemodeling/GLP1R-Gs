"""
#############################################
##  IMP script for GLP1R G protein complex modeling
##
#############################################
#
# Short modeling script using Crosslinking data to model GLP1R G protein complex
#
# Authors: Liping Sun
#
# References: Papers where this data is shown...
#
"""
import IMP
import RMF
import IMP.atom
import IMP.core
import IMP.algebra
import IMP.container
import IMP.rmf
import os
import sys
import time

'''
import DLFCN as dl
sys.setdlopenflags(dl.RTLD_NOW | dl.RTLD_GLOBAL)
'''
import IMP.pmi
import IMP.pmi.topology
import IMP.pmi.dof
import IMP.pmi.macros
import IMP.pmi.restraints
import IMP.pmi.restraints.em
import IMP.pmi.restraints.basic
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.restraints.crosslinking
import IMP.pmi.tools
import IMP.pmi.samplers
import IMP.pmi.output

# set time
start=time.time()

#---------------------------
# Define Input Files
#---------------------------
datadirectory = "../data/"
topology_file = datadirectory + "topology_free.txt"
outdirectory = "output"

# cross links, using the Maximum Calpha-Calpha distances as the cutoff - Angstrom
BS3_all = datadirectory + 'BS3.csv'
DSG_all = datadirectory + 'DSG.csv'
EDC_all = datadirectory + 'EDC.csv'
KArGO_all = datadirectory + 'KArGO.csv'

BS3_EUD = 35
DSG_EUD = 30
EDC_EUD = 25
KArGO_EUD = 40

#--------------------------
# Set MC Sampling Parameters
#--------------------------
#num_frames = 20000
num_frames = 150000
if '--test' in sys.argv:
    num_frames = 200
num_mc_steps = 10

#--------------------------
# Create movers and setup weights
#--------------------------
# rigid body movement params
rb_max_trans = 1.0 # 1-5, 1.0
rb_max_rot = 0.01 # 0-0.3, 0.01

# flexible bead movement
bead_max_trans = 2.00 # 1-5, 2.0

# super rigid body movement params
srb_max_trans = 1.0 # 1-5, 1.0
srb_max_rot = 0.01 # 0-0.3, 0.01

# Restraint weights
# one may set different values for inter and intra links depending on the specific purposes, maybe 1~10
# I used maximum Capha-Calpha distance as the cutoff so I set it as 7.5
xl_weight = 7.5 
mem_weight = 10

#--------------------------------
# Build the Model Representation
#--------------------------------
# Initialize model
m = IMP.Model()

# Create list of components from topology file
topology = IMP.pmi.topology.TopologyReader(topology_file,
                                           pdb_dir=datadirectory,
                                           fasta_dir=datadirectory
                                           )
# Get a list of domains
#domains = topology.get_components()

bs = IMP.pmi.macros.BuildSystem(m)
bs.add_state(topology)
representation, dof = bs.execute_macro(max_rb_trans=rb_max_trans,
                                       max_rb_rot=rb_max_rot,
                                       max_bead_trans=bead_max_trans,
                                       max_srb_trans=srb_max_trans,
                                       max_srb_rot=srb_max_rot)

#--------------------------
# Define Degrees of Freedom
#--------------------------
# Add default mover parameters to simulation
outputobjects = []  # reporter objects (for stat files)
sampleobjects = []  # sampling objects

#-----------------------------------
# Define Scoring Function Components
#-----------------------------------
# Here we are defining a number of restraints on our system.
#  For all of them we call add_to_model() so they are incorporated into scoring
# We also add them to the outputobjects list, so they are reported in stat
# files

# Excluded Volume Restraint
#  To speed up this expensive restraint, one might operate it at resolution 10-20. Here we have a small protein complex, so we do resolution = 1

sf = IMP.core.RestraintsScoringFunction(IMP.pmi.tools.get_restraint_set(m))
print("liping0", sf.evaluate(False))

humanglp1r = []
crs = []

moldict = bs.get_molecules()[0]
for molname in moldict:
    print(molname)
    for mol in moldict[molname]:

        cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(
            mol, scale=2.0, label=molname) # usually we use a scale of 2.0
        cr.add_to_model()
        outputobjects.append(cr)
        sampleobjects.append(cr)
        crs.append(cr)

        humanglp1r.append(mol)

sf = IMP.core.RestraintsScoringFunction(IMP.pmi.tools.get_restraint_set(m))
print("liping0", sf.evaluate(False))
ev1 = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(included_objects=humanglp1r,
                                                              resolution=1) # the resolution of the coarsest particle is usually chosen, slack is 10 by default in ExcludedVolumeRestraint
ev1.add_to_model()
ev1.set_label('ExcludedVolume')
outputobjects.append(ev1)

sf = IMP.core.RestraintsScoringFunction(IMP.pmi.tools.get_restraint_set(m))
print("liping1", sf.evaluate(False))

#--------------------------
# Membrane Restraint
#--------------------------

# we select residue 167 as the center of the membrane, TMD2
#p = IMP.atom.Selection(representation,molecule="GLP1R").get_selected_particles()[154]
#z_center = IMP.isd.Nuisance.setup_particle(p)
#z_center.set_nuisance(0.0)
mr1 = IMP.pmi.restraints.basic.MembraneRestraint(representation,
                                                 #z_nuisance = z_center.get_particle_index(),
                                                 center = 121.22, # average z coordinate of CA in center residues of all 7 helices
                                                 thickness=29.4, # thickness of hydrophobic lipid tails
                                                 softness=3.0, # refered to basic.py and test_membrane_restraint.py in imp_develop/
                                                 plateau=0.0000000001, # refered to basic.py and test_membrane_restraint.py in imp_develop/
                                                 objects_above=[(0,119,"GLP1R"),(178,200,"GLP1R"),(265,282,"GLP1R"),
                                                               (350,358,"GLP1R")],
                                                 objects_inside=[(120,142,"GLP1R"),(156,177,"GLP1R"),(201,226,"GLP1R"),
                                                                (245,264,"GLP1R"),(283,306,"GLP1R"),(326,349,"GLP1R"),(359,381,"GLP1R")],
                                                 objects_below=[(143,155,"GLP1R"),(227,244,"GLP1R"),(307,325,"GLP1R"),
                                                               (382,440,"GLP1R"),(1,394,"GA"), (1,350,"GB"),
                                                               (1,71,"GG"),(1,129,"NB")],
                                                 weight=mem_weight)
                                                 #)

#print(mr1.get_particles_inside())
#print(mr1.get_particles_above())
#print(mr1.get_particles_below())

mr1.add_to_model()
#mr1.set_label('Membrane')
outputobjects.append(mr1)

# For visualization purposes
mr1.create_membrane_density()

sf = IMP.core.RestraintsScoringFunction(IMP.pmi.tools.get_restraint_set(m))
print("liping2", sf.evaluate(False))

#--------------------------
# Cross linking restraints
#--------------------------
# Crosslinks - dataset
#  To use this restraint we have to first define the data format
#  Here assuming that it's a CSV file with column names that may need to change
# Other options include the linker length and the slope (for nudging
# components together)
kw = IMP.pmi.io.crosslink.CrossLinkDataBaseKeywordsConverter()
kw.set_unique_id_key("id")
kw.set_protein1_key("prot1")
kw.set_protein2_key("prot2")
kw.set_residue1_key("res1")
kw.set_residue2_key("res2")
kw.set_id_score_key(None)

#############
#BS3
#############

# Medium + High Confidence Intermolecular crosslinks
BS3 = IMP.pmi.io.crosslink.CrossLinkDataBase(kw)
BS3.create_set_from_file(BS3_all)

x_BS3 = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
    root_hier=representation,
    CrossLinkDataBase=BS3,
    length=BS3_EUD,
    label="BS3",
    resolution=1.0,
    slope=0.02)
x_BS3.rs.set_weight(xl_weight)
x_BS3.add_to_model()
sampleobjects.append(x_BS3)
outputobjects.append(x_BS3)
dof.get_nuisances_from_restraint(x_BS3)

sf = IMP.core.RestraintsScoringFunction(IMP.pmi.tools.get_restraint_set(m))
print("liping3", sf.evaluate(False))

#############
#DSG
#############
# Medium + High Confidence Intermolecular crosslinks
DSG = IMP.pmi.io.crosslink.CrossLinkDataBase(kw)
DSG.create_set_from_file(DSG_all)

x_DSG = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
    root_hier=representation,
    CrossLinkDataBase=DSG,
    length=DSG_EUD,
    label="DSG",
    resolution=1.0,
    slope=0.02)
x_DSG.rs.set_weight(xl_weight)
x_DSG.add_to_model()
sampleobjects.append(x_DSG)
outputobjects.append(x_DSG)
dof.get_nuisances_from_restraint(x_DSG)

sf = IMP.core.RestraintsScoringFunction(IMP.pmi.tools.get_restraint_set(m))
print("liping3", sf.evaluate(False))

#############
#EDC
#############
# Medium + High Confidence Intermolecular crosslinks
EDC = IMP.pmi.io.crosslink.CrossLinkDataBase(kw)
EDC.create_set_from_file(EDC_all)

x_EDC = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
    root_hier=representation,
    CrossLinkDataBase=EDC,
    length=EDC_EUD,
    label="EDC",
    resolution=1.0,
    slope=0.02)
x_EDC.rs.set_weight(xl_weight)
x_EDC.add_to_model()
sampleobjects.append(x_EDC)
outputobjects.append(x_EDC)
dof.get_nuisances_from_restraint(x_EDC)

sf = IMP.core.RestraintsScoringFunction(IMP.pmi.tools.get_restraint_set(m))
print("liping3", sf.evaluate(False))

#############
#KArGO
#############
# Medium + High Confidence Intermolecular crosslinks
KArGO = IMP.pmi.io.crosslink.CrossLinkDataBase(kw)
KArGO.create_set_from_file(KArGO_all)

x_KArGO = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
    root_hier=representation,
    CrossLinkDataBase=KArGO,
    length=KArGO_EUD,
    label="KArGO",
    resolution=1.0,
    slope=0.02)
x_KArGO.rs.set_weight(xl_weight)
x_KArGO.add_to_model()
sampleobjects.append(x_KArGO)
outputobjects.append(x_KArGO)
dof.get_nuisances_from_restraint(x_KArGO)

sf = IMP.core.RestraintsScoringFunction(IMP.pmi.tools.get_restraint_set(m))
print("liping3", sf.evaluate(False))

###Shuffling for random initial conformations 
IMP.pmi.tools.shuffle_configuration(humanglp1r)
dof.optimize_flexible_beads(200)

#--------------------------
# Monte-Carlo Sampling
#--------------------------
# This object defines all components to be sampled as well as the sampling
# protocol
mc1 = IMP.pmi.macros.ReplicaExchange0(m,
                                root_hier=representation,
                                monte_carlo_sample_objects=dof.get_movers(),
                                output_objects=outputobjects,
                                crosslink_restraints=sampleobjects,
                                monte_carlo_temperature=1.0,

                                simulated_annealing=True,
                                simulated_annealing_minimum_temperature=1.0,
                                simulated_annealing_maximum_temperature=1.5,
                                simulated_annealing_minimum_temperature_nframes=200,
                                simulated_annealing_maximum_temperature_nframes=20,

                                replica_exchange_minimum_temperature=1.5,
                                replica_exchange_maximum_temperature=2.5,

                                number_of_best_scoring_models=1,
                                monte_carlo_steps=num_mc_steps,
                                number_of_frames=num_frames,
                                global_output_directory=outdirectory)

# Start Sampling
mc1.execute
# sampling time
end=time.time()
print("Sampling finished succesfully")
print('running time={} s'.format(end-start))
