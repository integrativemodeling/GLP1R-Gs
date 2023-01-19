"""
#############################################
##  Modeller script for GLP-1R-Gs complex modeling
##
#############################################
#
# Short modeling script to complete density missing regions for the em structure of the GLP1R G protein complex
#
# Authors: Liping Sun
#
# References: Papers where this data is shown...
#
"""
from modeller import *
from modeller.automodel import *
import sys

env = environ(rand_seed=int(sys.argv[1]))
'''
6au6A @1.7 - R-Value Work: 0.168 
6eg8I @2.8 - R-Value Work: 0.217 
3sn6A @3.2 - R-Value Work: 0.226 
6b3jA @3.3
5vaiR @4.1  
'''

env.io.atom_files_directory = ['.', '../../atom_files']

# override the select_atoms method to select only the missing residues
# Please note that the residue numbers and chain IDs refer to the built model, not to the template(s)
# Modeller always names the model residues consistently, counting up from 1.
'''
# The same thing from chain A (required for multi-chain models):
# return selection(self.residue_range(’1:A’, ’2:A’))
# Residues 4, 6, 10:
# return selection(self.residues[’4’], self.residues[’6’],
#                  self.residues[’10’])
# All residues except 1-5:
# return selection(self) - selection(self.residue_range(’1’, ’5’))
'''
class MyModel(automodel):
                        
     def select_atoms(self):
        return selection(self.residue_range('1:A', '10:A'),
                         self.residue_range('65:A', '86:A'),
                         self.residue_range('255:A', '261:A'),
                             self.residue_range('395:B', '396:B'),
                             self.residue_range('735:C', '740:C'),
                             self.residue_range('797:C', '805:C'),
                             self.residues['806:D'],
                             self.residue_range('933:D', '934:D'),
                             self.residues['965:E'],
                             self.residue_range('966:F', '970:F'),
                             self.residue_range('1072:F', '1077:F'),
                             self.residue_range('1281:F', '1285:F'),
                             self.residue_range('1366:F', '1405:F'))
'''       
    def select_atoms(self):
        return selection(self.residue_range('1:A', '10:A'),
                         self.residue_range('65:A', '86:A'),
                         self.residue_range('255:A', '261:A'),
                             self.residue_range('1:B', '2:B'),
                             self.residue_range('1:G', '6:G'),
                             self.residue_range('63:G', '71:G'),
                             self.residues['1:N'],
                             self.residue_range('128:N', '129:N'),
                             self.residues['31:P'],
                             self.residue_range('1:R', '5:R'),
                             self.residue_range('107:R', '112:R'),
                             self.residue_range('316:R', '320:R'),
                             self.residue_range('401:R', '440:R'))
'''  
# create the automodel object
a = MyModel(env, alnfile='alignment.ali', # PIR alignment filename
              knowns=('6x18A','6x18B','6x18G','6x18N','6x18P','6x18R'), # codes of the templates
              sequence='humanglp1r',
              assess_methods=(assess.DOPE, assess.GA341, assess.normalized_dope)  # GA341 and DOPE model assessment
              ) # code of the target

# To get an approximate model very quickly (to get a rough idea of what it looks like, or to confirm that the alignment is reasonable)
# a.very_fast()

a.starting_model = 1 # index of the first model
a.ending_model = 20   # index of the last model
                     # determines how many models to calculate

#superpose the structures, this determines the number of the last model to build
#a.initial_malign3d = True
         
#If set to True, then all of the generated models (and loop models, if using loopmodel) are fit to the templates, and written out with the fit.pdb extension.
#a.final_malign3d = True

# a.auto_align() # get an automatic alignment, should not be used unless the unless the sequences are so similar (> 50%) that the calculated alignment is likely to be correct

# Very thorough VTFM optimization:
a.library_schedule = autosched.slow
a.max_var_iterations = 300

# Thorough MD optimization:
a.md_level = refine.slow

# Repeat the whole cycle 2 times and do not stop unless obj.func. > 1E6
a.repeat_optimization = 2
a.max_molpdf = 1e6

# a.set_output_model_format("MMCIF")  # request mmCIF rather than PDB outputs
if '--test' in sys.argv: a.ending_model = 0
a.make()

# cluster(cluster cut=1.5)
#a.cluster(cluster_cut=1.5)

#-----------------------------------------------
#Accessing output data after modeling is complete based on the DOPE score
#-----------------------------------------------

# Get a list of all successfully built models from a.outputs
ok_models = [x for x in a.outputs if x['failure'] is None]

# Rank the models by DOPE score
key = 'DOPE score'
if sys.version_info[:2] == (2,3):
    # Python 2.3’s sort doesn’t have a ’key’ argument
    ok_models.sort(lambda a,b: cmp(a[key], b[key]))
else:
    ok_models.sort(key=lambda a: a[key])
    
# Get top model
m = ok_models[0]

print("Top model: %s (DOPE score %.3f)" % (m['name'], m[key]))

