.. AnnotateMVF:

AnnotateMVF
===========

Description
-----------
None

Parameters
----------

``-h/--help``
^^^^^^^^^^^^^

**Description:** show this help message and exit

**Type:** boolean flag



``--mvf`` (required)
^^^^^^^^^^^^^^^^^^^^

**Description:** Input MVF file.

**Type:** file path; **Default:** None



``--out`` (required)
^^^^^^^^^^^^^^^^^^^^

**Description:** Output file

**Type:** file path; **Default:** None



``--filter-annotation/--filterannotation``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** Skip entries in the GFF file that contain this string in their 'Notes'

**Type:** None; **Default:** None



``--gff``
^^^^^^^^^

**Description:** Input gff annotation file.

**Type:** file path; **Default:** None



``--line-buffer/--linebuffer``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** Number of entries to store in memory at a time.

**Type:** integer; **Default:** 100000



``--nongenic-margin/--nongenicmargin``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** for --unnanotated-mode, only retain positions that are this number of bp away from an annotated region boundary

**Type:** integer; **Default:** 0



``--nongenic-mode/--nongenicmode``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** Instead of returning annotated genes, return the non-genic regions without without changing contigs or coordinates

**Type:** boolean flag



``--quiet``
^^^^^^^^^^^

**Description:** Suppress screen output.

**Type:** boolean flag


.. ConvertFasta2MVF:

ConvertFasta2MVF
================

Description
-----------
None

Parameters
----------

``-h/--help``
^^^^^^^^^^^^^

**Description:** show this help message and exit

**Type:** boolean flag



``--fasta`` (required)
^^^^^^^^^^^^^^^^^^^^^^

**Description:** input FASTA file(s)

**Type:** None; **Default:** None



``--out`` (required)
^^^^^^^^^^^^^^^^^^^^

**Description:** output MVF file

**Type:** None; **Default:** None



``--contig-by-file/--contigbyfile``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** Contigs are designated by separate files.

**Type:** boolean flag



``--contig-field/--contigfield``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** When headers are split by --field-sep, the 0-based index of the contig id.

**Type:** integer; **Default:** None



``--contig-ids/--contigids``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** manually specify one or more contig ids as ID:NAME

**Type:** None; **Default:** None



``--field-sep/--fieldsep``
^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** FASTA field separator; assumes '>database accession locus' format

**Type:** None; **Default:** None

**Choices:** ['TAB', 'SPACE', 'DBLSPACE', 'COMMA', 'MIXED', 'PIPE', 'AT', 'UNDER', 'DBLUNDER']


``--flavor``
^^^^^^^^^^^^

**Description:** type of file [dna] or protein

**Type:** None; **Default:** dna

**Choices:** ['dna', 'protein']


``--manual-coord/--manualcoord``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** manually specify reference coordinates for each file in the format CONTIGID:START..STOP, ...

**Type:** None; **Default:** None



``--overwrite``
^^^^^^^^^^^^^^^

**Description:** USE WITH CAUTION: force overwrite of outputs

**Type:** boolean flag



``--quiet``
^^^^^^^^^^^

**Description:** Suppress screen output.

**Type:** boolean flag



``--read-buffer/--readbuffer``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** number of lines to hold in READ buffer

**Type:** integer; **Default:** 100000



``--ref-label/--reflabel``
^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** label for reference sample

**Type:** None; **Default:** REF



``--sample-field/--samplefield``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** when headers are split by --field-sep, the 0-based index of the sample id

**Type:** integer; **Default:** None



``--sample-replace/--samplereplace``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** one or more TAG:NEWLABEL or TAG, items, if TAG found in sample label, replace with NEW (or TAG if NEW not specified) NEW and TAG must each be unique

**Type:** None; **Default:** None



``--write-buffer/--writebuffer``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** number of lines to hold in WRITE buffer

**Type:** integer; **Default:** 100000


.. ConvertMAF2MVF:

ConvertMAF2MVF
==============

Description
-----------
None

Parameters
----------

``-h/--help``
^^^^^^^^^^^^^

**Description:** show this help message and exit

**Type:** boolean flag



``--maf`` (required)
^^^^^^^^^^^^^^^^^^^^

**Description:** input MAF file

**Type:** file path; **Default:** None



``--out`` (required)
^^^^^^^^^^^^^^^^^^^^

**Description:** output MVF file

**Type:** file path; **Default:** None



``--sample-tags/--sampletags`` (required)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** one or more TAG:NEWLABEL or TAG, items, if TAG found in sample label, replace with NEW (or TAG if NEW not specified) NEW and TAG must each be unique.

**Type:** None; **Default:** None



``--line-buffer/--linebuffer``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** Number of entries to store in memory at a time.

**Type:** integer; **Default:** 100000



``--mvf-ref-label/--mvfreflabel``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** new label for reference sample (default='REF')

**Type:** None; **Default:** REF



``--overwrite``
^^^^^^^^^^^^^^^

**Description:** USE WITH CAUTION: force overwrite of outputs

**Type:** boolean flag



``--quiet``
^^^^^^^^^^^

**Description:** Suppress screen output.

**Type:** boolean flag



``--ref-tag/--reftag``
^^^^^^^^^^^^^^^^^^^^^^

**Description:** old reference tag

**Type:** None; **Default:** None


.. ConvertMVF2Fasta:

ConvertMVF2Fasta
================

Description
-----------
None

Parameters
----------

``-h/--help``
^^^^^^^^^^^^^

**Description:** show this help message and exit

**Type:** boolean flag



``--mvf`` (required)
^^^^^^^^^^^^^^^^^^^^

**Description:** Input MVF file.

**Type:** file path; **Default:** None



``--out`` (required)
^^^^^^^^^^^^^^^^^^^^

**Description:** Output path of FASTA file.

**Type:** file path; **Default:** None



``--buffer``
^^^^^^^^^^^^

**Description:** size (Mbp) of write buffer for each sample

**Type:** integer; **Default:** 10



``--label-type/--labeltype``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** Long labels with all metadata or short ids

**Type:** None; **Default:** long

**Choices:** ('long', 'short')


``--output-data/--outputdata``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** Output dna, rna or prot data.

**Type:** None; **Default:** None

**Choices:** ('dna', 'rna', 'prot')


``--quiet``
^^^^^^^^^^^

**Description:** Suppress screen output.

**Type:** boolean flag



``--regions``
^^^^^^^^^^^^^

**Description:** Path of a plain text file containing one more lines with entries 'contigid,stop,start' (one per line, inclusive coordinates) all data will be returned if left blank.

**Type:** file path; **Default:** None



``--samples``
^^^^^^^^^^^^^

**Description:** Specify comma-separated list of samples, Leave blank for all samples.

**Type:** None; **Default:** None



``--temp_dir/--tempdir``
^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** directory to write temporary fasta files

**Type:** None; **Default:** .


.. ConvertMVF2Phylip:

ConvertMVF2Phylip
=================

Description
-----------
None

Parameters
----------

``-h/--help``
^^^^^^^^^^^^^

**Description:** show this help message and exit

**Type:** boolean flag



``--mvf`` (required)
^^^^^^^^^^^^^^^^^^^^

**Description:** Input MVF file.

**Type:** file path; **Default:** None



``--out`` (required)
^^^^^^^^^^^^^^^^^^^^

**Description:** Output Phylip file.

**Type:** file path; **Default:** None



``--buffer``
^^^^^^^^^^^^

**Description:** size (bp) of write buffer for each sample

**Type:** integer; **Default:** 100000



``--contigs``
^^^^^^^^^^^^^

**Description:** Specify comma-separated list of contigs.

**Type:** None; **Default:** None



``--label-type/--labeltype``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** Long labels with all metadata or short ids

**Type:** None; **Default:** short

**Choices:** ('long', 'short')


``--output-data/--outputdata``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** Output dna, rna or prot data.

**Type:** None; **Default:** None

**Choices:** ('dna', 'rna', 'prot')


``--partition``
^^^^^^^^^^^^^^^

**Description:** Output a CSV partitions file with RAxMLformatting for use in partitioned phylogenetic methods.

**Type:** boolean flag



``--quiet``
^^^^^^^^^^^

**Description:** Suppress screen output.

**Type:** boolean flag



``--regions``
^^^^^^^^^^^^^

**Description:** Path of a plain text file containing one more lines with entries 'contigid,stop,start' (one per line, inclusive coordinates) all data will be returned if left blank.

**Type:** file path; **Default:** None



``--samples``
^^^^^^^^^^^^^

**Description:** Specify comma-separated list of samples, Leave blank for all samples.

**Type:** None; **Default:** None



``--temp_dir/--tempdir``
^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** directory to write temporary fasta files

**Type:** None; **Default:** .


.. ConvertVCF2MVF:

ConvertVCF2MVF
==============

Description
-----------
None

Parameters
----------

``-h/--help``
^^^^^^^^^^^^^

**Description:** show this help message and exit

**Type:** boolean flag



``--out`` (required)
^^^^^^^^^^^^^^^^^^^^

**Description:** output MVF file

**Type:** None; **Default:** None



``--alleles-from/--allelesfrom``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** get additional alignment columns
                from INFO fields (:-separated)

**Type:** None; **Default:** None



``--contig-ids/--contigids``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** manually specify one or more contig ids as ID;VCFLABE;MVFLABEL, note that VCFLABEL must match EXACTLY the contig string labels in the VCF file

**Type:** None; **Default:** None



``--field-sep/--fieldsep``
^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** VCF field separator (default='TAB')

**Type:** None; **Default:** TAB

**Choices:** ['TAB', 'SPACE', 'DBLSPACE', 'COMMA', 'MIXED']


``--line-buffer/--linebuffer``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** Number of entries to store in memory at a time.

**Type:** integer; **Default:** 100000



``--low-depth/--lowdepth``
^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** below this read depth coverage, convert to lower case set to 0 to disable

**Type:** integer; **Default:** 3



``--low-qual/--lowqual``
^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** below this quality convert to lower case set to 0 to disable

**Type:** integer; **Default:** 20



``--mask-depth/--maskdepth``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** below this read depth mask with N/n

**Type:** integer; **Default:** 1



``--mask-qual/--maskqual``
^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** low quality cutoff, bases replaced by N/- set to 0 to disable

**Type:** integer; **Default:** 3



``--no-autoindex/--noautoindex``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** do not automatically index contigs from the VCF

**Type:** boolean flag



``--out-flavor/--outflavor``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** choose output MVF flavor to include quality scores and/or indels

**Type:** None; **Default:** dna

**Choices:** ['dna', 'dnaqual', 'dnaqual-indel', 'dna-indel']


``--overwrite``
^^^^^^^^^^^^^^^

**Description:** USE WITH CAUTION: force overwrite of outputs

**Type:** boolean flag



``--qual``
^^^^^^^^^^

**Description:** Include Phred genotype quality (GQ) scores

**Type:** boolean flag



``--quiet``
^^^^^^^^^^^

**Description:** Suppress screen output.

**Type:** boolean flag



``--ref-label/--reflabel``
^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** label for reference sample (default='REF')

**Type:** None; **Default:** REF



``--sample-replace/--samplereplace``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** one or more TAG:NEWLABEL or TAG, items, if TAG found in sample label, replace with NEW (or TAG if NEW not specified) NEW and TAG must each be unique

**Type:** None; **Default:** None



``--vcf``
^^^^^^^^^

**Description:** VCF input file

**Type:** file path; **Default:** None


.. CalcCharacterCount:

CalcCharacterCount
==================

Description
-----------
None

Parameters
----------

``-h/--help``
^^^^^^^^^^^^^

**Description:** show this help message and exit

**Type:** boolean flag



``--mvf`` (required)
^^^^^^^^^^^^^^^^^^^^

**Description:** Input MVF file.

**Type:** file path; **Default:** None



``--out`` (required)
^^^^^^^^^^^^^^^^^^^^

**Description:** Output file

**Type:** file path; **Default:** None



``--base-match/--basematch``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** String of bases to match (i.e. numerator).

**Type:** None; **Default:** None



``--base-total/--basetotal``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** String of bases for total (i.e. denominator).

**Type:** None; **Default:** None



``--contigs``
^^^^^^^^^^^^^

**Description:** Specify comma-separated list of contigs.

**Type:** None; **Default:** None



``--mincoverage``
^^^^^^^^^^^^^^^^^

**Description:** Mininum sample coverage for sites.

**Type:** integer; **Default:** None



``--quiet``
^^^^^^^^^^^

**Description:** Suppress screen output.

**Type:** boolean flag



``--samples``
^^^^^^^^^^^^^

**Description:** Specify comma-separated list of samples, Leave blank for all samples.

**Type:** None; **Default:** None


.. CalcDstatCombinations:

CalcDstatCombinations
=====================

Description
-----------
None

Parameters
----------

``-h/--help``
^^^^^^^^^^^^^

**Description:** show this help message and exit

**Type:** boolean flag



``--mvf`` (required)
^^^^^^^^^^^^^^^^^^^^

**Description:** Input MVF file.

**Type:** file path; **Default:** None



``--out`` (required)
^^^^^^^^^^^^^^^^^^^^

**Description:** Output file

**Type:** file path; **Default:** None



``--contigs``
^^^^^^^^^^^^^

**Description:** Specify comma-separated list of contigs.

**Type:** None; **Default:** None



``--quiet``
^^^^^^^^^^^

**Description:** Suppress screen output.

**Type:** boolean flag



``--samples``
^^^^^^^^^^^^^

**Description:** Specify comma-separated list of samples, Leave blank for all samples.

**Type:** None; **Default:** None


.. CalcPairwiseDistances:

CalcPairwiseDistances
=====================

Description
-----------
None

Parameters
----------

``-h/--help``
^^^^^^^^^^^^^

**Description:** show this help message and exit

**Type:** boolean flag



``--mvf`` (required)
^^^^^^^^^^^^^^^^^^^^

**Description:** Input MVF file.

**Type:** file path; **Default:** None



``--out`` (required)
^^^^^^^^^^^^^^^^^^^^

**Description:** Output file

**Type:** file path; **Default:** None



``--mincoverage``
^^^^^^^^^^^^^^^^^

**Description:** Mininum sample coverage for sites.

**Type:** integer; **Default:** None



``--quiet``
^^^^^^^^^^^

**Description:** Suppress screen output.

**Type:** boolean flag



``--samples``
^^^^^^^^^^^^^

**Description:** Specify comma-separated list of samples, Leave blank for all samples.

**Type:** None; **Default:** None


.. CalcPatternCount:

CalcPatternCount
================

Description
-----------
None

Parameters
----------

``-h/--help``
^^^^^^^^^^^^^

**Description:** show this help message and exit

**Type:** boolean flag



``--mvf`` (required)
^^^^^^^^^^^^^^^^^^^^

**Description:** Input MVF file.

**Type:** file path; **Default:** None



``--out`` (required)
^^^^^^^^^^^^^^^^^^^^

**Description:** Output file

**Type:** file path; **Default:** None



``--mincoverage``
^^^^^^^^^^^^^^^^^

**Description:** Mininum sample coverage for sites.

**Type:** integer; **Default:** None



``--quiet``
^^^^^^^^^^^

**Description:** Suppress screen output.

**Type:** boolean flag



``--samples``
^^^^^^^^^^^^^

**Description:** Specify comma-separated list of samples, Leave blank for all samples.

**Type:** None; **Default:** None


.. CalcSampleCoverage:

CalcSampleCoverage
==================

Description
-----------
None

Parameters
----------

``-h/--help``
^^^^^^^^^^^^^

**Description:** show this help message and exit

**Type:** boolean flag



``--mvf`` (required)
^^^^^^^^^^^^^^^^^^^^

**Description:** Input MVF file.

**Type:** file path; **Default:** None



``--out`` (required)
^^^^^^^^^^^^^^^^^^^^

**Description:** Output file

**Type:** file path; **Default:** None



``--contigs``
^^^^^^^^^^^^^

**Description:** Specify comma-separated list of contigs.

**Type:** None; **Default:** None



``--quiet``
^^^^^^^^^^^

**Description:** Suppress screen output.

**Type:** boolean flag



``--samples``
^^^^^^^^^^^^^

**Description:** Specify comma-separated list of samples, Leave blank for all samples.

**Type:** None; **Default:** None


.. CheckMVF:

CheckMVF
========

Description
-----------
None

Parameters
----------

``-h/--help``
^^^^^^^^^^^^^

**Description:** show this help message and exit

**Type:** boolean flag



``--mvf`` (required)
^^^^^^^^^^^^^^^^^^^^

**Description:** Input MVF file.

**Type:** file path; **Default:** None



``--quiet``
^^^^^^^^^^^

**Description:** Suppress screen output.

**Type:** boolean flag


.. FilterMVF:

FilterMVF
=========

Description
-----------
None

Parameters
----------

``-h/--help``
^^^^^^^^^^^^^

**Description:** show this help message and exit

**Type:** boolean flag



``--mvf`` (required)
^^^^^^^^^^^^^^^^^^^^

**Description:** Input MVF file.

**Type:** file path; **Default:** None



``--out`` (required)
^^^^^^^^^^^^^^^^^^^^

**Description:** Output file

**Type:** file path; **Default:** None



``--actions``
^^^^^^^^^^^^^

**Description:** set of actions:args to perform, note these are done in order as listed

**Type:** None; **Default:** None



``--labels``
^^^^^^^^^^^^

**Description:** use sample labels instead of indices

**Type:** boolean flag



``--line-buffer/--linebuffer``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** Number of entries to store in memory at a time.

**Type:** integer; **Default:** 100000



``--more-help/--morehelp``
^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** prints full module list and descriptions

**Type:** boolean flag



``--overwrite``
^^^^^^^^^^^^^^^

**Description:** USE WITH CAUTION: force overwrite of outputs

**Type:** boolean flag



``--quiet``
^^^^^^^^^^^

**Description:** Suppress screen output.

**Type:** boolean flag



``--test``
^^^^^^^^^^

**Description:** manually input a line for testing

**Type:** None; **Default:** None



``--test-nchar/--textnchar``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** total number of samples for test string

**Type:** integer; **Default:** None



``--verbose``
^^^^^^^^^^^^^

**Description:** report every line (for debugging)

**Type:** boolean flag


.. InferGroupSpecificAllele:

InferGroupSpecificAllele
========================

Description
-----------
None

Parameters
----------

``-h/--help``
^^^^^^^^^^^^^

**Description:** show this help message and exit

**Type:** boolean flag



``--mvf`` (required)
^^^^^^^^^^^^^^^^^^^^

**Description:** Input MVF file.

**Type:** file path; **Default:** None



``--out`` (required)
^^^^^^^^^^^^^^^^^^^^

**Description:** Output file

**Type:** file path; **Default:** None



``--all-sample-trees/--allsampletrees``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** Makes trees from all samples instead of only the most complete sequence from each species

**Type:** boolean flag



``--allele-groups/--allelegroups``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** GROUP1:LABEL,LABEL GROUP2:LABEL,LABEL 

**Type:** None; **Default:** None



``--branch-lrt/--branchlrt``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** Specify the output file for and turn on the RAxML-PAML format LRT test scan for selection on the target branch in addition to the basic patterns scan

**Type:** file path; **Default:** None



``--chi-test/--chitest``
^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** Input two number values for expected Nonsynonymous and Synonymous expected values.

**Type:** None; **Default:** None



``--codeml-path/--codemlpath``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** Full path for PAML codeml executable.

**Type:** file path; **Default:** codeml



``--contigs``
^^^^^^^^^^^^^

**Description:** Specify comma-separated list of contigs.

**Type:** None; **Default:** None



``--end-contig/--endcontig``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** Numerical id for the ending contig.

**Type:** integer; **Default:** 100000000



``--gff``
^^^^^^^^^

**Description:** Input gff annotation file.

**Type:** file path; **Default:** None



``--mincoverage``
^^^^^^^^^^^^^^^^^

**Description:** Mininum sample coverage for sites.

**Type:** integer; **Default:** None



``--num-target-species/--targetspec``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** Specify the minimum number of taxa in the target set that are required to conduct analysis

**Type:** integer; **Default:** 1



``--outgroup``
^^^^^^^^^^^^^^

**Description:** Specify sample name with which to root trees.

**Type:** None; **Default:** None



``--output-align/--outputalign``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** Output alignment to this file path in phylip format.

**Type:** None; **Default:** None



``--paml-tmp/--pamltmp``
^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** path for temporary folder for PAML output files

**Type:** file path; **Default:** pamltmp



``--quiet``
^^^^^^^^^^^

**Description:** Suppress screen output.

**Type:** boolean flag



``--raxml-path/--raxmlpath``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** Full path to RAxML program executable.

**Type:** file path; **Default:** raxml



``--samples``
^^^^^^^^^^^^^

**Description:** Specify comma-separated list of samples, Leave blank for all samples.

**Type:** None; **Default:** None



``--species-groups/--speciesgroups``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** None

**Type:** None; **Default:** None



``--start-contig/--startcontig``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** Numerical ID for the starting contig.

**Type:** integer; **Default:** 0



``--target``
^^^^^^^^^^^^

**Description:** Specify the taxa labels that define the target lineage-specific branch to be tested.

**Type:** None; **Default:** None



``--use-labels/--uselabels``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** Use contig labels instead of IDs in output.

**Type:** boolean flag


.. InferTree:

InferTree
=========

Description
-----------
None

Parameters
----------

``-h/--help``
^^^^^^^^^^^^^

**Description:** show this help message and exit

**Type:** boolean flag



``--mvf`` (required)
^^^^^^^^^^^^^^^^^^^^

**Description:** Input MVF file.

**Type:** file path; **Default:** None



``--out`` (required)
^^^^^^^^^^^^^^^^^^^^

**Description:** Output file

**Type:** file path; **Default:** None



``--bootstrap``
^^^^^^^^^^^^^^^

**Description:** turn on rapid bootstrapping for RAxML and perform specified number of replicates

**Type:** integer; **Default:** None



``--choose-allele/--chooseallele/--hapmode``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** Chooses how heterozygous alleles are handled. (none=no splitting (default); randomone=pick one allele randomly (recommended); randomboth=pick two alleles randomly, but keep both; major=pick the more common allele; minor=pick the less common allele; majorminor= pick the major in 'a' and minor in 'b'

**Type:** None; **Default:** none

**Choices:** ['none', 'randomone', 'randomboth', 'major', 'minor', 'majorminor']


``--contigs``
^^^^^^^^^^^^^

**Description:** Specify comma-separated list of contigs.

**Type:** None; **Default:** None



``--duplicate-seq/--duplicateseq``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** dontuse=remove duplicate sequences prior to RAxML tree inference, then add them to the tree manually as zero-branch-length sister taxa; keep=keep in for RAxML tree inference (may cause errors for RAxML); remove=remove entirely from alignment

**Type:** None; **Default:** dontuse

**Choices:** ['dontuse', 'keep', 'remove']


``--min-depth/--mindepth``
^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** minimum number of alleles per site

**Type:** integer; **Default:** 4



``--min-seq-coverage/--minseqcoverage``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** proportion of total alignment a sequencemust cover to be retianed [0.1]

**Type:** float; **Default:** 0.1



``--min-sites/--minsites``
^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** minimum number of sites 

**Type:** integer; **Default:** 100



``--output-contig-labels/--outputcontiglabels``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** Output will use contig labels instead of id numbers.

**Type:** boolean flag



``--output-empty/--outputempty``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** Include entries of windows with no data in output.

**Type:** boolean flag



``--quiet``
^^^^^^^^^^^

**Description:** Suppress screen output.

**Type:** boolean flag



``--raxml-model/--raxmlmodel``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** choose RAxML model

**Type:** None; **Default:** GTRGAMMA



``--raxml-opts/--raxmlopts``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** specify additional RAxML arguments as a double-quotes encased string

**Type:** None; **Default:** 



``--raxml-outgroups/--raxmloutgroups``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** Comma-separated list of outgroup taxon labels to use in RAxML.

**Type:** None; **Default:** None



``--raxml-path/--raxmlpath``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** RAxML path for manual specification.

**Type:** None; **Default:** raxml



``--root-with/--rootwith``
^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** Comma-separated list of taxon labels to root trees with after RAxML

**Type:** None; **Default:** None



``--samples``
^^^^^^^^^^^^^

**Description:** Specify comma-separated list of samples, Leave blank for all samples.

**Type:** None; **Default:** None



``--temp-dir/--tempdir``
^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** Temporary directory path

**Type:** file path; **Default:** ./raxmltemp



``--temp-prefix/--tempprefix``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** Temporary file prefix

**Type:** None; **Default:** mvftree


.. JoinMVF:

JoinMVF
=======

Description
-----------
None

Parameters
----------

``-h/--help``
^^^^^^^^^^^^^

**Description:** show this help message and exit

**Type:** boolean flag



``--mvf`` (required)
^^^^^^^^^^^^^^^^^^^^

**Description:** One or more mvf files.

**Type:** file path; **Default:** None



``--out`` (required)
^^^^^^^^^^^^^^^^^^^^

**Description:** Output file

**Type:** file path; **Default:** None



``--line-buffer/--linebuffer``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** Number of entries to store in memory at a time.

**Type:** integer; **Default:** 100000



``--main_header_file/--mainheaderfile``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** Output file will use same headers as this input file (default=first in list).

**Type:** None; **Default:** None



``--new-contigs/--newcontigs``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** By default, contigs are matched between files using their text labels in the header. Use this option to turn matching off and treat each file's contigs as distinct.

**Type:** boolean flag



``--newsamples``
^^^^^^^^^^^^^^^^

**Description:** By default, samples are matched between files using their text labels in the header. Use this option to turn matching off and treat each file's sample columns as distinct.

**Type:** boolean flag



``--overwrite``
^^^^^^^^^^^^^^^

**Description:** USE WITH CAUTION: force overwrite of outputs

**Type:** boolean flag



``--quiet``
^^^^^^^^^^^

**Description:** Suppress screen output.

**Type:** boolean flag


.. PlotChromoplot:

PlotChromoplot
==============

Description
-----------
None

Parameters
----------

``-h/--help``
^^^^^^^^^^^^^

**Description:** show this help message and exit

**Type:** boolean flag



``--mvf`` (required)
^^^^^^^^^^^^^^^^^^^^

**Description:** Input MVF file.

**Type:** file path; **Default:** None



``--outgroup`` (required)
^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** 1 or more outgroups to use for quartets

**Type:** None; **Default:** None



``--samples`` (required)
^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** 3 or more taxa to use for quartets

**Type:** None; **Default:** None



``--colors``
^^^^^^^^^^^^

**Description:** three colors to use for chromoplot

**Type:** None; **Default:** None

**Choices:** {'lgrey': (250, 250, 250), 'dgrey': (192, 192, 192), 'black': (0, 0, 0), 'white': (255, 255, 255), 'red': (192, 0, 0), 'orange': (217, 95, 2), 'yellow': (192, 192, 0), 'green': (0, 192, 0), 'blue': (0, 0, 192), 'teal': (27, 158, 119), 'puce': (117, 112, 179), 'purple': (192, 0, 192), 'none': ()}


``--contigs``
^^^^^^^^^^^^^

**Description:** Enter the ids of one or more contigs in the order they will appear in the chromoplot. (defaults to all ids in order present in MVF)

**Type:** None; **Default:** None



``--empty-mask/--emptymask``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** Mask empty regions with this color.

**Type:** None; **Default:** none

**Choices:** {'lgrey': (250, 250, 250), 'dgrey': (192, 192, 192), 'black': (0, 0, 0), 'white': (255, 255, 255), 'red': (192, 0, 0), 'orange': (217, 95, 2), 'yellow': (192, 192, 0), 'green': (0, 192, 0), 'blue': (0, 0, 192), 'teal': (27, 158, 119), 'puce': (117, 112, 179), 'purple': (192, 0, 192), 'none': ()}


``--info-track/--infotrack``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** Include an additional coverage information track that will show empty, uninformative, and informative loci. (Useful for ranscriptomes/RAD or other reduced sampling.

**Type:** boolean flag



``--majority``
^^^^^^^^^^^^^^

**Description:** Plot only 100% shading in the majority track  rather than shaded proportions in all tracks.

**Type:** boolean flag



``--out-prefix/--outprefix``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** Output prefix (not required).

**Type:** None; **Default:** None



``--plot-type/--plottype``
^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** PNG image (default) or graph via matplotlib (experimental)

**Type:** None; **Default:** image

**Choices:** ['graph', 'image']


``--quiet``
^^^^^^^^^^^

**Description:** Suppress screen output.

**Type:** boolean flag



``--xscale``
^^^^^^^^^^^^

**Description:** Width (in number of pixels) for each window

**Type:** integer; **Default:** 1



``--yscale``
^^^^^^^^^^^^

**Description:** Height (in number of pixels) for each track

**Type:** integer; **Default:** 20


.. TranslateMVF:

TranslateMVF
============

Description
-----------
None

Parameters
----------

``-h/--help``
^^^^^^^^^^^^^

**Description:** show this help message and exit

**Type:** boolean flag



``--mvf`` (required)
^^^^^^^^^^^^^^^^^^^^

**Description:** Input MVF file.

**Type:** file path; **Default:** None



``--out`` (required)
^^^^^^^^^^^^^^^^^^^^

**Description:** Output file

**Type:** file path; **Default:** None



``--filter-annotation/--filterannotation``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** skip GFF entries with text matching this in their 'Notes' field

**Type:** None; **Default:** None



``--gff``
^^^^^^^^^

**Description:** Input GFF3 file. If GFF3 not provided, alignments are assumed to be in-frame coding sequences.

**Type:** file path; **Default:** None



``--line-buffer/--linebuffer``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** Number of entries to store in memory at a time.

**Type:** integer; **Default:** 100000



``--output-data/--outputdata``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** protein=single data column of protein alleles; codon=four columns with: protein frame1 frame2 frame3

**Type:** None; **Default:** codon

**Choices:** ['protein', 'codon']


``--overwrite``
^^^^^^^^^^^^^^^

**Description:** USE WITH CAUTION: force overwrite of outputs

**Type:** boolean flag



``--quiet``
^^^^^^^^^^^

**Description:** Suppress screen output.

**Type:** boolean flag


