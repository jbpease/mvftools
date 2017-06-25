.. mvf_join:

mvf_join
========

Description
-----------

This program checks an MVF file for inconsistencies or errors


Parameters
----------

``-h/--help``
^^^^^^^^^^^^^

**Description:** show this help message and exit

**Type:** boolean flag



``-i/--mvf`` (required)
^^^^^^^^^^^^^^^^^^^^^^^

**Description:** One or more mvf files.

**Type:** file path; **Default:** None



``-o--out`` (required)
^^^^^^^^^^^^^^^^^^^^^^

**Description:** Output mvf file.

**Type:** file path; **Default:** None



``-c/--newcontigs``
^^^^^^^^^^^^^^^^^^^

**Description:** By default, contigs are matched between files using their text labels in the header. Use this option to turn matching off and treat each file's contigs as distinct.

**Type:** boolean flag



``--overwrite``
^^^^^^^^^^^^^^^

**Description:** USE WITH CAUTION: force overwrite of outputs

**Type:** boolean flag



``--quiet``
^^^^^^^^^^^

**Description:** suppress progress meter

**Type:** boolean flag



``-s/--newsamples``
^^^^^^^^^^^^^^^^^^^

**Description:** By default, samples are matched between files using their text labels in the header. Use this option to turn matching off and treat each file's sample columns as distinct.

**Type:** boolean flag



``-B/--linebuffer``
^^^^^^^^^^^^^^^^^^^

**Description:** number of entries to write in a block

**Type:** integer; **Default:** 100000



``-M/--main_header_file``
^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** Output file will use same headers as this input file (default=first in list).

**Type:** None; **Default:** None


