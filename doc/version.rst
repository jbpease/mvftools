===============
Version History
===============
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

=======
License
=======
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
