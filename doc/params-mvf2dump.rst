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


