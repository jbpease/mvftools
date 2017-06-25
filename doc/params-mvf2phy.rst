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

