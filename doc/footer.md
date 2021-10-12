# Retired Modules

* JoinMVF > Use ConcatenateMVF instead.
* CheckMVF > Use VerifyMVF instead.

# Version History

**v.0.6.2**

2021-10-12: Added 'randomallele' action in FilterMVF. Fixed issue with 'collapsemerge' in FilterMVF. Added ability to label columns in 'collapsemerge'.  Fixed program errors in TranslateMVF and MVF2Fasta.  Moved the insertion of commands used to generate MVF to notes lines. **Functionailty change**: The column specifications in the collapse functions are now done with respect to the original columns. You no longer need to adjust the numeric columns in later actions to account for changes in prior actions. **Functionality change**: RAxML-ng is now the default for InferTrees.  You can still run RAxML 8.2 by setting the options manually. 

**v.0.6.1**

2021-04-20: Updates to fix problems with the ConvertMAF2MVF Module.  Note that commas now separate the --sample-tags and --ref-tag is now required.

**v.0.6.0**

*2020-12-20: Major Update* - Major update to the back-end of MVF for speed and stability.  Adds support for faster file access by creating an MVF index, fixes to the filtering modules, changes to the specification of sample and contig labels to improve usability. Major improvements to the MVFTranslate module.  Added support for tetraploid and hexaploid VCF files. Many other fixes


**v.0.5.4**

2019-07-16: Adds the CalcAllCharacterCountPerSample function, continued upgrades to the screen output to provide more realtime information.  Other small fixes to the CalcPairwiseDistances module in ambiguous character mode.


**v.0.5.3**

2019-07-14: Critial update, strongly recommend updating to this version.  Major efficiency fix in the base iteration modules.  Several key bug fixes implemented in FilterMVF, MergeMVF.  Enhanced support for ambiguous sequences and polyploids in several modules including CalcPairwiseDistance.  Restructuring of FilterMVF for cleaner syntax.


**v.0.5.2**

2019-06-29: Added MergeMVF to join several files together (still experimental, use with caution).  JoinMVF is now called ConcatenateMVF to avoid confusion.  CheckMVF changed to VerifyMVF to make it more clear.  ConvertVCF2MVF now has experimental support for tetraploid and hexaploid VCF files through the --ploidy flag.  Other small fixes to the software and manual to update issues with the VCF interpreter.

**v.0.5.1**

2018-02-01: Changes to the --sample and --outgroup arguments for some calculations into separate --sample-indices and --sample-labels arguments.  This fixes an issue where if the sample labels are numerical they are misinterpreted when specified at the command line. All sample/outgroup indices or labels should be specified as a single comma-separated list.

**v.0.5.0**

*2017-11-27 - Major Upgrade*: Change to single-command structure

**v.2017-06-25**

*Major Upgrade*: Full manual documentation added, standardization and cleanup of paramaters and upgrades and bugfixes throughout.

**v.2017-05-18**

Fixes to VCF conversion for compatibility

**v.2017-04-10**

Added MVF-to-Phylip output conversion ``mvf2phy``

**v.2017-03-25**

Multiple bug fixes, merged and removed the development instance

**v.2016-02-15**

Fix to vcf2mvf for VCF with truncated entries

**v.2016-10-25**

Efficiency upgrades for mvfbase entry iteration.

**v.2016-09-10**

Minor fixes to gz reading and MVF chromoplot shading

**v.2016-08-02**

Python3 conversion, integrate analysis_base

**v.2016-01-11**

fix for dna ambiguity characters

**v.2016-01-01**

Python3 compatiblity fix

**v.2015-12-31**

Header changesand cleanup

**v.2015-12-15**

Python3 compatibilty fix

**v.2015-09-04**

Small style fixes

**v.2015-06-09**

MVF1.2.1 upgrade

**v.2015-02-26**

Efficiency upgrades for iterators

**v.2015-02-01**

First Public Release

# License
MVFtools is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
MVFtools is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with MVFtools.  If not, see <http://www.gnu.org/licenses/>.
