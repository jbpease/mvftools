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


