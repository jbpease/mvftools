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


