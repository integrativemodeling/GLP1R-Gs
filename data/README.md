# Representation and data used for the integrative structure determination of the GLP-1 receptor-Gs complex

Master data directory; used by the modeling script to generate the ensemble of solutions of the Glp-1R-Gs complex.

## List of files and directories:
`data` contains all relevant data

## Sequence and structures files
- `humanglp1r.fasta` contains all the sequences of the Glp-1R-Gs complex.
- `6x18.pdb` is the EM structure of the Glp-1R-Gs complex (PDB code 6X18; [article](10.1016/j.molcel.2020.09.020))
- `GLP1_z.pdb, GLP1R_z.pdb, GA_noGDP_z.pdb, GB_z.pdb, GG_z.pdb, and NB_z.pdb` are the complete models of all subunits of the Glp-1R-Gs complex derived from `humanglp1r_run_48_renr.pdb`, which is the complete model of the Glp-1R-Gs complex, obtained using MODELLER (see [humanglp1r](./humanglp1r))

## Cross-linking datasets for the GLP-1 receptor-Gs complex
- `BS3.csv`: BS3 - 40 intermolecular and 188 intramolecular cross-links for the GLP-1 receptor-Gs complex
- `DSG.csv`: DSG - 35 intermolecular and 168 intramolecular cross-links for the GLP-1 receptor-Gs complex
- `EDC`: EDC - 15 intermolecular and 111 intramolecular cross-links for the GLP-1 receptor-Gs complex
- `KArGO`: KArGO - 2 intermolecular and 38 intramolecular cross-links for the GLP-1 receptor-Gs complex

### Summary of the Glp-1R-Gs cross-linking dataset
![](./Glp-1R-Gs_XL.pdf)

## Topology files used for modeling: defining rigid bodies and flexible beads regions.
- `topology_free.txt` is the topology file for the Glp-1R-Gs complex.

## Information

_Author(s)_: Liping Sun

_License_: [LGPL](http://www.gnu.org/licenses/old-licenses/lgpl-2.1.html).
This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2 of the License, or (at your option) any later version.

_Publications_:
