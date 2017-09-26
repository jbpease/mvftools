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


