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


