Program Parameter Descriptions
##############################

.. fasta2mvf:

fasta2mvf
=========

Description
-----------

This program is used to convert a FASTA file into MVF format.


Parameters
----------

``-h/--help``
^^^^^^^^^^^^^

**Description:** show this help message and exit

**Type:** boolean flag



``-i/--fasta`` (required)
^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** input FASTA file(s)

**Type:** None; **Default:** None



``-o/--out`` (required)
^^^^^^^^^^^^^^^^^^^^^^^

**Description:** output MVF file

**Type:** None; **Default:** None



``-c/--contig-ids/--contigids``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** manually specify one or more contig ids as ID:NAME

**Type:** None; **Default:** None



``--contig-by-file/--contigbyfile``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** Contigs are designated by separate files.

**Type:** boolean flag



``--contig-field/--contigfield``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** When headers are split by --field-sep, the 0-based index of the contig id.

**Type:** integer; **Default:** None



``-f/--flavor``
^^^^^^^^^^^^^^^

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



``-s/--sample-replace/--samplereplace``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** one or more TAG:NEWLABEL or TAG, items, if TAG found in sample label, replace with NEW (or TAG if NEW not specified) NEW and TAG must each be unique

**Type:** None; **Default:** None



``--sample-field/--samplefield``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** when headers are split by --field-sep, the 0-based index of the sample id

**Type:** integer; **Default:** None



``-B/--read-buffer/--readbuffer``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** number of lines to hold in READ buffer

**Type:** integer; **Default:** 100000



``-F/--field-sep/--fieldsep``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** FASTA field separator; assumes '>database accession locus' format

**Type:** None; **Default:** None

**Choices:** ['TAB', 'SPACE', 'DBLSPACE', 'COMMA', 'MIXED', 'PIPE', 'AT', 'UNDER', 'DBLUNDER']


``-R/--ref-label/--reflabel``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** label for reference sample

**Type:** None; **Default:** REF



``-W/--write-buffer/--writebuffer``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** number of lines to hold in WRITE buffer

**Type:** integer; **Default:** 100000


.. maf2mvf:

maf2mvf
=======

Description
-----------

This program analyzes a DNA MVF alignment using the modules specified below,
use the --morehelp option for additional module information.


Parameters
----------

``-h/--help``
^^^^^^^^^^^^^

**Description:** show this help message and exit

**Type:** boolean flag



``-i/--maf`` (required)
^^^^^^^^^^^^^^^^^^^^^^^

**Description:** input MAF file

**Type:** file path; **Default:** None



``-o/--out`` (required)
^^^^^^^^^^^^^^^^^^^^^^^

**Description:** output MVF file

**Type:** file path; **Default:** None



``-s/--sample-tags/--sampletags`` (required)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** one or more TAG:NEWLABEL or TAG, items, if TAG found in sample label, replace with NEW (or TAG if NEW not specified) NEW and TAG must each be unique.

**Type:** None; **Default:** None



``--overwrite``
^^^^^^^^^^^^^^^

**Description:** None

**Type:** boolean flag



``-B/--line-buffer/--linebuffer``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** number of lines to hold in read/write buffer

**Type:** integer; **Default:** 100000



``-M/--mvf-ref-label/--mvfreflabel``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** new label for reference sample (default='REF')

**Type:** None; **Default:** REF



``-R/--ref-tag/--reftag``
^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** old reference tag

**Type:** None; **Default:** None


.. mvf2dump:

mvf2dump
========

Description
-----------

This program exports the entirety of an MVF to FASTA format,
with many fewer options than mvf2fasta.py.  This is designed
to export large MVF files faster, but with less specific
formatting and region-finding options.


Parameters
----------

``-h/--help``
^^^^^^^^^^^^^

**Description:** show this help message and exit

**Type:** boolean flag



``-i/--mvf`` (required)
^^^^^^^^^^^^^^^^^^^^^^^

**Description:** Input MVF file.

**Type:** file path; **Default:** None



``-o/--outprefix`` (required)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** Target FASTA file

**Type:** file path; **Default:** None



``-d/--outdata``
^^^^^^^^^^^^^^^^

**Description:** output dna, rna or prot data

**Type:** None; **Default:** None

**Choices:** ('dna', 'rna', 'prot')


``--quiet``
^^^^^^^^^^^

**Description:** suppress screen output

**Type:** boolean flag



``-s/--samples``
^^^^^^^^^^^^^^^^

**Description:** One or more taxon labels, leave blank for all

**Type:** None; **Default:** None



``-t/--tmpdir``
^^^^^^^^^^^^^^^

**Description:** directory to write temporary fasta files

**Type:** file path; **Default:** .



``-B/--buffer``
^^^^^^^^^^^^^^^

**Description:** size (Mbp) of write buffer for each sample

**Type:** integer; **Default:** 10


.. mvf2fasta:

mvf2fasta
=========

Description
-----------

This program takes an MVF file and converts the data to a FASTA file


Parameters
----------

``-h/--help``
^^^^^^^^^^^^^

**Description:** show this help message and exit

**Type:** boolean flag



``-i/--mvf`` (required)
^^^^^^^^^^^^^^^^^^^^^^^

**Description:** Input MVF file.

**Type:** file path; **Default:** None



``-o/--out`` (required)
^^^^^^^^^^^^^^^^^^^^^^^

**Description:** target FASTA file

**Type:** file path; **Default:** None



``-r/--regions`` (required)
^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** A file path to a plain-text file withone region per line formatted asformatted as: contigid,start,stop(coordinates are inclusive)

**Type:** None; **Default:** None



``-d/--outdata``
^^^^^^^^^^^^^^^^

**Description:** Output dna, rna or prot data.

**Type:** None; **Default:** None

**Choices:** ('dna', 'rna', 'prot')


``-l/--labeltype``
^^^^^^^^^^^^^^^^^^

**Description:** Long labels with all metadata or short ids

**Type:** None; **Default:** long

**Choices:** ('long', 'short')


``--quiet``
^^^^^^^^^^^

**Description:** suppress screen output

**Type:** boolean flag



``-s/--samples``
^^^^^^^^^^^^^^^^

**Description:** One or more taxon labels, leave blank for all

**Type:** None; **Default:** None



``-t/--tmpdir``
^^^^^^^^^^^^^^^

**Description:** directory to write temporary fasta files

**Type:** None; **Default:** .



``-B/--buffer``
^^^^^^^^^^^^^^^

**Description:** size (Mbp) of write buffer for each sample

**Type:** integer; **Default:** 10


.. mvf2phy:

mvf2phy
=======

Description
-----------

This program is used to export MVF data to Phylip format.


Parameters
----------

``-h/--help``
^^^^^^^^^^^^^

**Description:** show this help message and exit

**Type:** boolean flag



``-i/--mvf`` (required)
^^^^^^^^^^^^^^^^^^^^^^^

**Description:** Input MVF file.

**Type:** file path; **Default:** None



``-o/--out`` (required)
^^^^^^^^^^^^^^^^^^^^^^^

**Description:** Output Phylip file.

**Type:** file path; **Default:** None



``-d/--outdata``
^^^^^^^^^^^^^^^^

**Description:** Output dna, rna or prot data.

**Type:** None; **Default:** None

**Choices:** ('dna', 'rna', 'prot')


``-p/--partition``
^^^^^^^^^^^^^^^^^^

**Description:** Output a CSV partitions file with RAxMLformatting for use in partitioned phylogenetic methods.

**Type:** boolean flag



``--quiet``
^^^^^^^^^^^

**Description:** suppress screen output

**Type:** boolean flag



``-r/--region``
^^^^^^^^^^^^^^^

**Description:** Path of a plain text file containing one more lines with entries 'contigid,stop,start' (one per line, inclusive coordinates) all data will be returned if left blank.

**Type:** file path; **Default:** None



``-s/--samples``
^^^^^^^^^^^^^^^^

**Description:** One or more taxon labels, leave blank for all.

**Type:** None; **Default:** None



``-t/--tmpdir``
^^^^^^^^^^^^^^^

**Description:** directory to write temporary fasta files

**Type:** None; **Default:** .



``-B/--buffer``
^^^^^^^^^^^^^^^

**Description:** size (bp) of write buffer for each sample

**Type:** integer; **Default:** 100000



``-L/--labeltype``
^^^^^^^^^^^^^^^^^^

**Description:** Long labels with all metadata or short ids

**Type:** None; **Default:** short

**Choices:** ('long', 'short')

.. mvf_analyze_codon:

mvf_analyze_codon
=================

Description
-----------

This program analyzes a codon MVF using several analysis modules.
Run 'python3 mvf_analyze_codon.py --morehelp' for details on
module functions.


Parameters
----------

module
^^^^^^

**Description:** None

**Type:** None; **Default:** None

**Choices:** ('Coverage', 'GroupUniqueAlleleWindow', 'PiDiversityWindow', 'PairwiseNS')


``-h/--help``
^^^^^^^^^^^^^

**Description:** show this help message and exit

**Type:** boolean flag



``--all-sample-trees/--allsampletrees``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** (GroupUniqueAlleleWindow) Makes trees from all samples instead of only the most complete sequence from each species

**Type:** boolean flag



``--allele-groups/--allelegroups``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** GROUP1:LABEL,LABEL GROUP2:LABEL,LABEL
                                (GroupUniqueAlleleWindow)

**Type:** None; **Default:** None



``--branchlrt``
^^^^^^^^^^^^^^^

**Description:** (GroupUniqueAlleleWindow) Specify the output file for and turn on the RAxML-PAML format LRT test scan for selection on the target branch in addition to the basic patterns scan

**Type:** file path; **Default:** None



``-c/--contigs``
^^^^^^^^^^^^^^^^

**Description:** List of space-separated contig ids.

**Type:** None; **Default:** None



``-g/--gff``
^^^^^^^^^^^^

**Description:** GFF3 file for use in annotation

**Type:** None; **Default:** None



``-i/--mvf``
^^^^^^^^^^^^

**Description:** Input MVF file.

**Type:** file path; **Default:** None



``-m/--mincoverage``
^^^^^^^^^^^^^^^^^^^^

**Description:** Minimum number of samples with alleles needed to use site for analysis.

**Type:** integer; **Default:** None



``--morehelp``
^^^^^^^^^^^^^^

**Description:** Get additional information on modules.

**Type:** boolean flag



``--num-target-species/--targetspec``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** (GroupUniqueAlleleWindow) Specify the minimum number of taxa in the target set that are required to conduct analysis

**Type:** integer; **Default:** 1



``-o/--out``
^^^^^^^^^^^^

**Description:** output file

**Type:** file path; **Default:** None



``--output-align/--outputalign``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** (GroupUniqueAlleleWindow) Output alignment to this file path in phylip format.

**Type:** None; **Default:** None



``--pamltmp``
^^^^^^^^^^^^^

**Description:** path for temporary folder for PAML output files

**Type:** file path; **Default:** pamltmp



``-s/--samples``
^^^^^^^^^^^^^^^^

**Description:** List of space-separated sample names.

**Type:** None; **Default:** None



``--species-groups/--speciesgroups``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** None

**Type:** None; **Default:** None



``--target``
^^^^^^^^^^^^

**Description:** (GroupUniqueAlleleWindow) Specify the taxa labels that define the target lineage-specific branch to be tested.

**Type:** None; **Default:** None



``-w/--windowsize``
^^^^^^^^^^^^^^^^^^^

**Description:** Window size in bp, use -1 for whole contig.

**Type:** integer; **Default:** -1



``-x/--chi-test/--chitest``
^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** (GroupUniqueAlleleWindow,PairwiseDNDS)Input two number values for expected Nonsynonymous and Synonymous expected values. 

**Type:** None; **Default:** None



``-E/--endcontig``
^^^^^^^^^^^^^^^^^^

**Description:** Numerical id for the ending contig.

**Type:** integer; **Default:** 100000000



``-L/--uselabels``
^^^^^^^^^^^^^^^^^^

**Description:** Use contig labels instead of IDs in output.

**Type:** boolean flag



``-O/--outgroup``
^^^^^^^^^^^^^^^^^

**Description:** (GroupUniqueAlleleWindow) Specify sample name with which to root trees.

**Type:** None; **Default:** None



``-P/--codemlpath``
^^^^^^^^^^^^^^^^^^^

**Description:** Full path for PAML codeml executable.

**Type:** file path; **Default:** codeml



``-S/--startcontig``
^^^^^^^^^^^^^^^^^^^^

**Description:** Numerical ID for the starting contig.

**Type:** integer; **Default:** 0



``-X/--raxmlpath``
^^^^^^^^^^^^^^^^^^

**Description:** Full path to RAxML program executable.

**Type:** file path; **Default:** raxml


.. mvf_analyze_dna:

mvf_analyze_dna
===============

Description
-----------

This program analyzes a DNA MVF alignment using the modules specified below,
use the --morehelp option for additional module information.


Parameters
----------

module
^^^^^^

**Description:** analysis module to run

**Type:** None; **Default:** None

**Choices:** ('BaseCountWindow', 'Coverage', 'DstatComb', 'PairwiseDistance', 'PairwiseDistanceWindow', 'PatternCount', 'PatternList')


``-h/--help``
^^^^^^^^^^^^^

**Description:** show this help message and exit

**Type:** boolean flag



``-i/--mvf`` (required)
^^^^^^^^^^^^^^^^^^^^^^^

**Description:** Input MVF file.

**Type:** file path; **Default:** None



``-o/--out`` (required)
^^^^^^^^^^^^^^^^^^^^^^^

**Description:** output file

**Type:** file path; **Default:** None



``--base-match``
^^^^^^^^^^^^^^^^

**Description:** [BaseCountWindow] string of bases to match (i.e. numerator).

**Type:** None; **Default:** None



``--base-total``
^^^^^^^^^^^^^^^^

**Description:** [BaseCountWindow] string of bases for total (i.e. denominator).

**Type:** None; **Default:** None



``-c/--contigs``
^^^^^^^^^^^^^^^^

**Description:** limit analyses to these contigs

**Type:** None; **Default:** None



``-m/--mincoverage``
^^^^^^^^^^^^^^^^^^^^

**Description:** mininum sample coverage for site

**Type:** integer; **Default:** None



``--morehelp``
^^^^^^^^^^^^^^

**Description:** get additional information on modules

**Type:** boolean flag



``-s/--samples``
^^^^^^^^^^^^^^^^

**Description:** limit analyses to these samples

**Type:** None; **Default:** None



``-w/--windowsize``
^^^^^^^^^^^^^^^^^^^

**Description:** window size, use -1 to use whole contigs

**Type:** integer; **Default:** 100000


.. mvf_annotate:

mvf_annotate
============

Description
-----------

This program takes a DNA MVF alignment and annotates the output into
gene boudaries.


Parameters
----------

``-h/--help``
^^^^^^^^^^^^^

**Description:** show this help message and exit

**Type:** boolean flag



``-g/--gff``
^^^^^^^^^^^^

**Description:** Input gff annotation file.

**Type:** file path; **Default:** None



``-i/--mvf``
^^^^^^^^^^^^

**Description:** Input MVF file.

**Type:** file path; **Default:** None



``-o/--out``
^^^^^^^^^^^^

**Description:** Output annotated MVF file

**Type:** file path; **Default:** None



``--overwrite``
^^^^^^^^^^^^^^^

**Description:** USE WITH CAUTION: force overwrite of outputs

**Type:** boolean flag



``--quiet``
^^^^^^^^^^^

**Description:** suppress progress meter

**Type:** boolean flag



``-B/--linebuffer``
^^^^^^^^^^^^^^^^^^^

**Description:** Number of entries to store in memory at a time.

**Type:** integer; **Default:** 100000



``-F/--filter_annotation``
^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** Skip entries in the GFF file that contain this string in their 'Notes'

**Type:** None; **Default:** None



``-M/--nongenic-margin``
^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** for --unnanotated-mode, only retain positions that are this number of bp away from an annotated region boundary

**Type:** integer; **Default:** 0



``-N/--nongenic-mode``
^^^^^^^^^^^^^^^^^^^^^^

**Description:** Instead of returning annotated genes, return the non-genic regions without without changing contigs or coordinates

**Type:** boolean flag


.. mvf_check:

mvf_check
=========

Description
-----------

This program checks an MVF file for inconsistencies or errors


Parameters
----------

mvf
^^^

**Description:** Input MVF file.

**Type:** file path; **Default:** None



``-h/--help``
^^^^^^^^^^^^^

**Description:** show this help message and exit

**Type:** boolean flag


.. mvf_chromoplot:

mvf_chromoplot
==============

Description
-----------

This program creates a chromoplot from an MVF alignment.
A chromoplot shows a genome-wide diagram of different
evolutionary histories for a given quartet of taxa.


Parameters
----------

``-h/--help``
^^^^^^^^^^^^^

**Description:** show this help message and exit

**Type:** boolean flag



``-i/--mvf`` (required)
^^^^^^^^^^^^^^^^^^^^^^^

**Description:** Input MVF file.

**Type:** file path; **Default:** None



``-s/--samples`` (required)
^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** 3 or more taxa to use for quartets

**Type:** None; **Default:** None



``-G/--outgroup`` (required)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** 1 or more outgroups to use for quartets

**Type:** None; **Default:** None



``-c/--contigs``
^^^^^^^^^^^^^^^^

**Description:** Enter the ids of one or more contigs in the order they will appear in the chromoplot. (defaults to all ids in order present in MVF)

**Type:** None; **Default:** None



``-o/--outprefix``
^^^^^^^^^^^^^^^^^^

**Description:** Output prefix (not required).

**Type:** None; **Default:** None



``-q/--quiet``
^^^^^^^^^^^^^^

**Description:** suppress all output messages

**Type:** boolean flag



``-w/--windowsize``
^^^^^^^^^^^^^^^^^^^

**Description:** None

**Type:** integer; **Default:** 100000



``-x/--xscale``
^^^^^^^^^^^^^^^

**Description:** Width (in number of pixels) for each window

**Type:** integer; **Default:** 1



``-y/--yscale``
^^^^^^^^^^^^^^^

**Description:** Height (in number of pixels) for each track

**Type:** integer; **Default:** 20



``-C/--colors``
^^^^^^^^^^^^^^^

**Description:** three colors to use for chromoplot

**Type:** None; **Default:** None

**Choices:** {'lgrey': (250, 250, 250), 'dgrey': (192, 192, 192), 'black': (0, 0, 0), 'white': (255, 255, 255), 'red': (192, 0, 0), 'orange': (217, 95, 2), 'yellow': (192, 192, 0), 'green': (0, 192, 0), 'blue': (0, 0, 192), 'teal': (27, 158, 119), 'puce': (117, 112, 179), 'purple': (192, 0, 192), 'none': ()}


``-E/--emptymask``
^^^^^^^^^^^^^^^^^^

**Description:** Mask empty regions with this color.

**Type:** None; **Default:** none

**Choices:** {'lgrey': (250, 250, 250), 'dgrey': (192, 192, 192), 'black': (0, 0, 0), 'white': (255, 255, 255), 'red': (192, 0, 0), 'orange': (217, 95, 2), 'yellow': (192, 192, 0), 'green': (0, 192, 0), 'blue': (0, 0, 192), 'teal': (27, 158, 119), 'puce': (117, 112, 179), 'purple': (192, 0, 192), 'none': ()}


``-I/--infotrack``
^^^^^^^^^^^^^^^^^^

**Description:** Include an additional coverage information track that will show empty, uninformative, and informative loci. (Useful for ranscriptomes/RAD or other reduced sampling.

**Type:** boolean flag



``-M/--majority``
^^^^^^^^^^^^^^^^^

**Description:** Plot only 100% shading in the majority track  rather than shaded proportions in all tracks.

**Type:** boolean flag



``-P/--plottype``
^^^^^^^^^^^^^^^^^

**Description:** PNG image (default) or graph via matplotlib (experimental)

**Type:** None; **Default:** image

**Choices:** ['graph', 'image']

.. mvf_filter:

mvf_filter
==========

Description
-----------

This program filters an MVF alignment using the modules specified below,
use the --morehelp option for additional module information.


Parameters
----------

``-h/--help``
^^^^^^^^^^^^^

**Description:** show this help message and exit

**Type:** boolean flag



``-a/--actions``
^^^^^^^^^^^^^^^^

**Description:** set of actions:args to perform, note these are done in order as listed

**Type:** None; **Default:** None



``-i/--mvf``
^^^^^^^^^^^^

**Description:** Input MVF file.

**Type:** file path; **Default:** None



``-l/--labels``
^^^^^^^^^^^^^^^

**Description:** use sample labels instead of indices

**Type:** boolean flag



``--morehelp``
^^^^^^^^^^^^^^

**Description:** prints full module list and descriptions

**Type:** boolean flag



``-o/--out``
^^^^^^^^^^^^

**Description:** Output MVF file

**Type:** file path; **Default:** None



``--overwrite``
^^^^^^^^^^^^^^^

**Description:** USE WITH CAUTION: force overwrite of outputs

**Type:** boolean flag



``-q/--quiet``
^^^^^^^^^^^^^^

**Description:** suppress progress meter

**Type:** boolean flag



``--test``
^^^^^^^^^^

**Description:** manually input a line for testing

**Type:** None; **Default:** None



``--test-nchar``
^^^^^^^^^^^^^^^^

**Description:** total number of samples for test string

**Type:** integer; **Default:** None



``-B/--linebuffer``
^^^^^^^^^^^^^^^^^^^

**Description:** number of lines to write at once to MVF

**Type:** integer; **Default:** 100000



``-V/--verbose``
^^^^^^^^^^^^^^^^

**Description:** report every line (for debugging)

**Type:** boolean flag


.. mvf_join:

mvf_join
========

Description
-----------

This program checks an MVF file for inconsistencies or errors


Parameters
----------

``-h/--help``
^^^^^^^^^^^^^

**Description:** show this help message and exit

**Type:** boolean flag



``-i/--mvf`` (required)
^^^^^^^^^^^^^^^^^^^^^^^

**Description:** One or more mvf files.

**Type:** file path; **Default:** None



``-o--out`` (required)
^^^^^^^^^^^^^^^^^^^^^^

**Description:** Output mvf file.

**Type:** file path; **Default:** None



``-c/--newcontigs``
^^^^^^^^^^^^^^^^^^^

**Description:** By default, contigs are matched between files using their text labels in the header. Use this option to turn matching off and treat each file's contigs as distinct.

**Type:** boolean flag



``--overwrite``
^^^^^^^^^^^^^^^

**Description:** USE WITH CAUTION: force overwrite of outputs

**Type:** boolean flag



``--quiet``
^^^^^^^^^^^

**Description:** suppress progress meter

**Type:** boolean flag



``-s/--newsamples``
^^^^^^^^^^^^^^^^^^^

**Description:** By default, samples are matched between files using their text labels in the header. Use this option to turn matching off and treat each file's sample columns as distinct.

**Type:** boolean flag



``-B/--linebuffer``
^^^^^^^^^^^^^^^^^^^

**Description:** number of entries to write in a block

**Type:** integer; **Default:** 100000



``-M/--main_header_file``
^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** Output file will use same headers as this input file (default=first in list).

**Type:** None; **Default:** None


.. mvf_translate:

mvf_translate
=============

Description
-----------

This program translates a DNA MVF file into a codon or protein MVF file
using a GFF3 annotation file.


Parameters
----------

``-h/--help``
^^^^^^^^^^^^^

**Description:** show this help message and exit

**Type:** boolean flag



``-i/--mvf`` (required)
^^^^^^^^^^^^^^^^^^^^^^^

**Description:** Input MAF file

**Type:** file path; **Default:** None



``-o/--out`` (required)
^^^^^^^^^^^^^^^^^^^^^^^

**Description:** Output MVF file

**Type:** None; **Default:** None



``-g/--gff``
^^^^^^^^^^^^

**Description:** Input GFF3 file. If GFF3 not provided, alignments are assumed to be in-frame coding sequences.

**Type:** file path; **Default:** None



``--overwrite``
^^^^^^^^^^^^^^^

**Description:** USE WITH CAUTION: force overwrite of outputs

**Type:** boolean flag



``--quiet``
^^^^^^^^^^^

**Description:** suppress progress meter

**Type:** boolean flag



``-t--outtype``
^^^^^^^^^^^^^^^

**Description:** protein=single data column of protein alleles; codon=four columns with: protein frame1 frame2 frame3

**Type:** None; **Default:** codon

**Choices:** ['protein', 'codon']


``-B/--line-buffer/--linebuffer``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** number of entries to write in a block

**Type:** integer; **Default:** 100000



``-F/--filter-annotation``
^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** skip GFF entries with text matching this in their 'Notes' field

**Type:** None; **Default:** None


.. mvf_window_tree:

mvf_window_tree
===============

Description
-----------

This program makes phylogenies from individual genomic windows of
a DNA MVF alignment (Requires: BioPython).


Parameters
----------

``-h/--help``
^^^^^^^^^^^^^

**Description:** show this help message and exit

**Type:** boolean flag



``-i/--mvf`` (required)
^^^^^^^^^^^^^^^^^^^^^^^

**Description:** Input MVF file.

**Type:** file path; **Default:** None



``-o/--out`` (required)
^^^^^^^^^^^^^^^^^^^^^^^

**Description:** Tree list output text file.

**Type:** file path; **Default:** None



``-b/--bootstrap``
^^^^^^^^^^^^^^^^^^

**Description:** turn on rapid bootstrapping for RAxML and perform specified number of replicates

**Type:** integer; **Default:** None



``-c/--contigs``
^^^^^^^^^^^^^^^^

**Description:** Contig ids to use in analysis (default=all)

**Type:** None; **Default:** None



``-d/--duplicateseq``
^^^^^^^^^^^^^^^^^^^^^

**Description:** dontuse=remove duplicate sequences prior to RAxML tree inference, then add them to the tree manually as zero-branch-length sister taxa; keep=keep in for RAxML tree inference (may cause errors for RAxML); remove=remove entirely from alignment

**Type:** None; **Default:** dontuse

**Choices:** ['dontuse', 'keep', 'remove']


``-e/--outputempty``
^^^^^^^^^^^^^^^^^^^^

**Description:** Include entries of windows with no data in output.

**Type:** boolean flag



``-g/--raxml-outgroups/--raxml_outgroups``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** Outgroups taxon labels to use in RAxML.

**Type:** None; **Default:** None



``-m/--raxml_model``
^^^^^^^^^^^^^^^^^^^^

**Description:** choose RAxML model

**Type:** None; **Default:** GTRGAMMA



``--outputcontiglabels``
^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** Output will use contig labels instead of id numbers.

**Type:** boolean flag



``--quiet``
^^^^^^^^^^^

**Description:** suppress screen output

**Type:** boolean flag



``-r/--rootwith``
^^^^^^^^^^^^^^^^^

**Description:** Root output trees with these taxa after RAxML.

**Type:** None; **Default:** None



``-s/--samples``
^^^^^^^^^^^^^^^^

**Description:** One or more taxon labels (default=all)

**Type:** None; **Default:** None



``--tempdir``
^^^^^^^^^^^^^

**Description:** Temporary directory path

**Type:** file path; **Default:** ./raxmltemp



``--tempprefix``
^^^^^^^^^^^^^^^^

**Description:** Temporary file prefix

**Type:** None; **Default:** mvftree



``-w/--windowsize``
^^^^^^^^^^^^^^^^^^^

**Description:** specify genomic region size, or use -1 for whole contig

**Type:** integer; **Default:** 10000



``-A/--choose_allele/--hapmode``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** Chooses how heterozygous alleles are handled. (none=no splitting (default); randomone=pick one allele randomly (recommended); randomboth=pick two alleles randomly, but keep both; major=pick the more common allele; minor=pick the less common allele; majorminor= pick the major in 'a' and minor in 'b'

**Type:** None; **Default:** none

**Choices:** ['none', 'randomone', 'randomboth', 'major', 'minor', 'majorminor']


``-C/--minseqcoverage``
^^^^^^^^^^^^^^^^^^^^^^^

**Description:** proportion of total alignment a sequence
                                must cover to be retianed [0.1]

**Type:** float; **Default:** 0.1



``-D/--mindepth``
^^^^^^^^^^^^^^^^^

**Description:** minimum number of alleles per site

**Type:** integer; **Default:** 4



``-M/--minsites``
^^^^^^^^^^^^^^^^^

**Description:** minimum number of sites 

**Type:** integer; **Default:** 100



``-R/--raxmlopts``
^^^^^^^^^^^^^^^^^^

**Description:** specify additional RAxML arguments as a double-quotes encased string

**Type:** None; **Default:** 



``-X/--raxmlpath``
^^^^^^^^^^^^^^^^^^

**Description:** RAxML path for manual specification.

**Type:** None; **Default:** raxml


.. vcf2mvf:

vcf2mvf
=======

Description
-----------

MVFtools: Multisample Variant Format Toolkit
James B. Pease and Ben K. Rosenzweig
http://www.github.org/jbpease/mvftools


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



``--vcf`` (required)
^^^^^^^^^^^^^^^^^^^^

**Description:** input VCF file

**Type:** file path; **Default:** None



``--allelesfrom``
^^^^^^^^^^^^^^^^^

**Description:** get additional alignment columns
                                from INFO fields (:-separated)

**Type:** None; **Default:** None



``--contigids``
^^^^^^^^^^^^^^^

**Description:** manually specify one or more contig ids
                                 as ID;VCFLABE;MVFLABEL, note that
                                 VCFLABEL must match EXACTLY the contig string
                                 labels in the VCF file

**Type:** None; **Default:** None



``--fieldsep``
^^^^^^^^^^^^^^

**Description:** VCF field separator (default='TAB')

**Type:** None; **Default:** TAB

**Choices:** ['TAB', 'SPACE', 'DBLSPACE', 'COMMA', 'MIXED']


``--linebuffer``
^^^^^^^^^^^^^^^^

**Description:** number of lines to hold in read/write buffer

**Type:** integer; **Default:** 100000



``--lowdepth``
^^^^^^^^^^^^^^

**Description:** below this read depth coverage, convert to lower case set to 0 to disable

**Type:** integer; **Default:** 3



``--lowqual``
^^^^^^^^^^^^^

**Description:** below this quality convert to lower case
                                set to 0 to disable

**Type:** integer; **Default:** 20



``--maskdepth``
^^^^^^^^^^^^^^^

**Description:** below this read depth mask with N/n

**Type:** integer; **Default:** 1



``--maskqual``
^^^^^^^^^^^^^^

**Description:** low quality cutoff, bases replaced by N/-
                             set to 0 to disable

**Type:** integer; **Default:** 3



``--no_autoindex``
^^^^^^^^^^^^^^^^^^

**Description:** do not automatically index contigs from the VCF

**Type:** boolean flag



``--outflavor``
^^^^^^^^^^^^^^^

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



``--reflabel``
^^^^^^^^^^^^^^

**Description:** label for reference sample (default='REF')

**Type:** None; **Default:** REF



``--samplereplace``
^^^^^^^^^^^^^^^^^^^

**Description:** one or more TAG:NEWLABEL or TAG, items,
                                if TAG found in sample label, replace with
                                NEW (or TAG if NEW not specified)
                                NEW and TAG must each be unique

**Type:** None; **Default:** None


