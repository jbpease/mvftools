# mvftools
###Multisample Variant Format ToolKit

## Description:
This repository contains the specification for the Multisample Variant Format (MVF), which is designed for compact storage and efficient analysis of multi-genome and multi-transcriptome datasets.  The programs provided in MVFtools support this format, both with conversion utilities, filtering and transformation programs, and data analysis and visualization modules.  MVF format is designed specifically for biological data analysis, since sequence data is encoded based on the information content at a particular aligned sequence site.  This contextual encoding allows for rapid computation of phylogenetic and population genetic analyses, and small file sizes that enable data sharing and distribution.

## Requirements:
MVFtools requires Python 2.7+ or Python 3.x, but can be run on any operating system.

Programs within MVFtools have specific dependencies as follows (\* denotes non-Python program):
* mvf_chromoplot.py: [PIL](http://www.pythonware.com/products/pil/), [Scipy](http://www.scipy.org/)
* mvf_window_tree.py: [Biopython](http://www.biopython.org/) 1.6+, [Numpy](http://www.numpy.org/), [Scipy](http://www.scipy.org/), [RAxML](http://sco.h-its.org/exelixis/web/software/raxml/index.html/)\* 

## Installation:
Simply clone MVFtools to use the scripts, unless other dependencies are required (see above). All external depencies should be installed as recommended by their individual documentation.

##Contributors:
* James B. Pease ([website](http://pages.iu.edu/~jbpease/))([@jamesbpease](https://twitter.com/jamesbpease/))
* Ben K. Rosenzweig

## Version History:

2015-02-01: First Public Release
2015-02-26: Patch upgrade for multiple modules.  Featured upgrades are efficiency increases to mvfbase and major overhaul of mvf_window_tree.


## License
This file is part of MVFtools.

MVFtools is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

MVFtools is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with Foobar.  If not, see (http://www.gnu.org/licenses/).

## FAQ

1. What is an example of a simple MVF workflow using VCF?
1.1. After aligning your reads and sorting them using Using [SAMtools and BCFtools](http://www.htslib.org/) and 1+ *sorted* .bam files, this command is used to generate a multi-sample VCF. 

    samtools mpileup -uD -f $GENOME.fasta $SPECIESA.sorted.bam $SPECIESB.sorted.bam ... $SPECIESN.sorted.bam | bcftools view -cg - > $INTERMED.vcf
    
    vcf2mvf.py $INTERMED.vcf --out $OUTPUT.mvf --labelreplace SPECIESA:SPECA SPECIESB:SPECB ... SPECIESN:SPECN --maskqual 0 --lowqual 20 --maskdepth 1 --lowdepth 3 --quiet --chromindex chr01 chr02 chr03 ...

The '--labelpreplace' amd '--chromindex' flags are optional, as VCF files should be converted automatically.

