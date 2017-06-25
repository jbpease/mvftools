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


