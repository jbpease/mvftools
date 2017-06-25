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


