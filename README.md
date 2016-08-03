# mvftools
###Multisample Variant Format ToolKit
**Update 2016-09-02: MVFtools now requires Python 3.x**
**Update 2016-07-15: the mvftools-dev repository has merged with mvftools under the branch *develop*, see this branch for the latest features**

## Description:
This repository contains the specification for the Multisample Variant Format (MVF), which is designed for compact storage and efficient analysis of multi-genome and multi-transcriptome datasets.  The programs provided in MVFtools support this format, both with conversion utilities, filtering and transformation programs, and data analysis and visualization modules.  MVF format is designed specifically for biological data analysis, since sequence data is encoded based on the information content at a particular aligned sequence site.  This contextual encoding allows for rapid computation of phylogenetic and population genetic analyses, and small file sizes that enable data sharing and distribution.

##Citation Information:
Pease JB and BK Rosenzweig. 2016. "Encoding Data Using Biological Principles: the Multisample Variant Format for Phylogenomics and Population Genomics" *IEEE/ACM Transactions on Computational Biology and Bioinformatics*. In press. http://www.dx.doi.org/10.1109/tcbb.2015.2509997

## Requirements:
MVFtools requires **Python 3.x** (will not work on Python 2.x), but can be run on any operating system.

Programs within MVFtools have specific dependencies as follows (\* denotes non-Python program):
* mvf_chromoplot.py: [PIL](http://www.pythonware.com/products/pil/), [Scipy](http://www.scipy.org/)
* mvf_window_tree.py: [Biopython](http://www.biopython.org/) 1.6+, [Numpy](http://www.numpy.org/), [Scipy](http://www.scipy.org/), [RAxML](http://sco.h-its.org/exelixis/web/software/raxml/index.html/)\* 

## Installation:
Simply clone MVFtools to use the scripts, unless other dependencies are required (see above). All external depencies should be installed as recommended by their individual documentation.

##Contributors:
* James B. Pease ([website](http://peaselab.github.io)([@jamesbpease](https://twitter.com/jamesbpease/))
* Ben K. Rosenzweig

##FAQ and Questions/Comments:
Visit the Google Groups site for FAQs and to ask question:
https://groups.google.com/forum/#!forum/mvftools

##Other Studies using MVFtools:
* http://dx.doi.org/10.1111/mec.13679
* http://dx.doi.org/10.1111/mec.13610
* http://dx.doi.org/10.1371/journal.pbio.1002379

## License
This file is part of MVFtools.

MVFtools is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

MVFtools is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with Foobar.  If not, see (http://www.gnu.org/licenses/).


