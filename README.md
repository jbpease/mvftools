![logo](https://github.com/jbpease/mvftools/blob/master/doc/logo.png)

# mvftools

## Authors

* James B. Pease ([website](http://peaselab.github.io)([@jamesbpease](https://twitter.com/jamesbpease/))
* Ben K. Rosenzweig

## Description
Multisample Variant Format (MVF), which is designed for compact storage and efficient analysis of multi-genome and multi-transcriptome datasets.  The programs provided in MVFtools support this format, both with conversion utilities, filtering and transformation programs, and data analysis and visualization modules.  MVF format is designed specifically for biological data analysis, since sequence data is encoded based on the information content at a particular aligned sequence site.  This contextual encoding allows for rapid computation of phylogenetic and population genetic analyses, and small file sizes that enable data sharing and distribution.

## Citation Information

Pease JB and BK Rosenzweig. 2018. "Encoding Data Using Biological Principles: the Multisample Variant Format for Phylogenomics and Population Genomics" *IEEE/ACM Transactions on Computational Biology and Bioinformatics*. 15(4):1231-1238.  http://www.dx.doi.org/10.1109/tcbb.2015.2509997

Please also include the URL [https://www.github.com/jbpease/mvftools] in your methods section where the program is referenced.

(Note this paper was originally published online in 2015, but did not receive final citation page numbering until 2018.  You may see older citations as 2015, which is the same paper.)

## Manual

Please see the full manual at (https://github.com/jbpease/mvftools/blob/master/MVFtools.pdf)


## Requirements

MVFtools requires:
  * **Python 3.x** (will not work on Python 2.x), but can be run on any operating system.
  * BioPython ([https://biopython.org/])
  * Numpy ([https://www.numpy.org/])
  * Scipy ([https://www.scipy.org/])

*Optionally required for certain modules:*
  * RAxML (8.x recommended; [https://sco.h-its.org/exelixis/web/software/
raxml/index.html])
  * PAML ([http://abacus.gene.ucl.ac.uk/software/paml.html])


## FAQ and Questions/Comments
See the manual above and visit the Google Groups site for FAQs and to ask question:
https://groups.google.com/forum/#!forum/mvftools

## Examples of papers that have used MVFtools

* [http://dx.doi.org/10.1111/mec.13679]
* [http://dx.doi.org/10.1111/mec.13610]
* [http://dx.doi.org/10.1371/journal.pbio.1002379]
* [https://doi.org/10.1038/s41467-018-04963-6]
* [https://doi.org/10.1093/gbe/evy227]
* [https://doi.org/10.1093/molbev/msy198]
* [https://doi.org/10.1534/genetics.118.301831]


## License
This file is part of MVFtools.

MVFtools is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

MVFtools is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with Foobar.  If not, see (http://www.gnu.org/licenses/).


