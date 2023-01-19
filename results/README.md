# Results for  the GLP-1 receptor-Gs complex

Directory for the main results of the modeling pipelines for the GLP-1R-Gs complex that follows the four stages of integrative modeling. (see README.md in parent directory)

## List of files and directories:

`results` contains the results reported for the structure of the GLP-1R-Gs complexs

- `integrative_structures` contains the results of integrative modeling for the GLP-1R-Gs complex
    * [Integrative Structures GLP-1R-Gs](./integrative_structures/localization_probability_densities/) contains the localization probability densities of the integrative structures of GLP-1R-G,s at a precision of 16.9Å
    * [Sampling Precision](./sampling_precision) contains the results of the sampling exhaustiveness test
- The centroid structure of the ensemble

- `xlink_satisfaction` contains the analysis of the cross-links on the integrative structures of the GLP-1R-Gs complex

- `backmap` contains the scripts to backmap coarse-grained integrative structures to atomic ones
- rmf2pdb.py is the script to generate .pdb files of the coarse-grained integrative structures from .rmf files
- ca2all.py is the script to generate .pdb files of atomic integrative structures from coarse-grained ones

- `contact_analysis` contains the scripts to compute residue contacts and buried surface area of the GLP-1R-Gs complex

- `representative_structures` contains .pdb files of seven representative structures

_Author(s)_:  Liping Sun

_License_: [LGPL](http://www.gnu.org/licenses/old-licenses/lgpl-2.1.html).
This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2 of the License, or (at your option) any later version.

_Publications_:
