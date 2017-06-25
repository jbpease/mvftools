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


