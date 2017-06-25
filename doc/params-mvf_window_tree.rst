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


