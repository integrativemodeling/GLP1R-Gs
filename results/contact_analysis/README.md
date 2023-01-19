# Scripts for analyzing residue contacts and buried surface area of the GLP-1 receptor-Gs complex

### Compute residue contacts in the GLP-1R-Gs interface

We use pymol to compute residue contacts in the GLP-1R-Gs interface.

A residue pair (one from GLP-1R, the other from a Gs subunit) is defined as a residue-residue contact if and only if the residue pair can potentially form a hydrogen bond, a salt bridge, or a hydrophobic contact according to the following criteria: 

1. A hydrogen bond is defined when the distance between a probable donor and acceptor with an electronegative atom (N, O, S) is smaller than or equal to 3.5 Å. 
2. A salt bridge is defined when the distance between a positively charged residue (K, R, H) and a negatively charged residue (D, E) is smaller than or equal to 4.0 Å. 
3. A hydrophobic contact is defined when the distance between two hydrophobic residues is within 4.0 Å, or the distance between the atoms in one residue and the center of the benzene ring in another aromatic residue is smaller than or equal to 4.0Å, or the distance between the centers of two benzene ring is smaller than or equal to 5.0 Å.
    
The distance cut-offs used here are generally accepted empirical values according to previous GPCR structural studies. The interface contacts are also obtained from the complete cryo-EM structure of the GLP-1R-Gs complex using the same criteria. Detailed descriptions of the contacts can be found in the manuscript.

#### Dependency:
* Pymol 2.1 (https://pymolwiki.org/index.php/Main_Page)


#### Usage:

```Shell

conda create --no-default-packages -n myenv python=3.6 pymol
source activate myenv 
python ./residue_contact/residue_contact.py
python ./residue_contact/aroma_count.py

```
### Compute buried surface area in the GLP-1R-Gs interface

To compute the buried surface area in the GLP-1R-Gs interface, we use the program FreeSASA 2.0, through the Sharke-Rupley algorithm with a probe radius of 1.4 Å, and refer to definition of interface area using solvent accessible surface area (SASA) as in https://www.ebi.ac.uk/pdbe/pisa/. 

The calculation of the SASA of the interface is calculated as follows:
If T is the total SASA of the combined molecular components A:B where ":" is the interface between A and B, then you can see this total excludes the interface area. If TA and TB are the surface area of each molecular component A and B, then the missing area for A:B due to ":" is the same on A and B in isolation, so the measure of contact ":" missing in the A:B complex is counted twice - once in TA and once in TB. Thus, the interface Area = (TA + TB - T) / 2.0

We compute the surface area of GLP-1R-Gs interface, GLP-1R-Gs-ICL1 interface, GLP-1R-Gs-ICL2 interface and GLP-1R-Gs-ICL3 interface for comparison.


#### Dependency:
* FreeSASA 2.0 (https://freesasa.github.io/)

#### Usage:

```Shell

pip install freesasa
python ./interface_area/compute_surface_area.py

```
