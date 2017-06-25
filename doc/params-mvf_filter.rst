.. mvf_filter:

mvf_filter
==========

Description
-----------

This program filters an MVF alignment using the modules specified below,
use the --morehelp option for additional module information.


Parameters
----------

``-h/--help``
^^^^^^^^^^^^^

**Description:** show this help message and exit

**Type:** boolean flag



``-a/--actions``
^^^^^^^^^^^^^^^^

**Description:** set of actions:args to perform, note these are done in order as listed

**Type:** None; **Default:** None



``-i/--mvf``
^^^^^^^^^^^^

**Description:** Input MVF file.

**Type:** file path; **Default:** None



``-l/--labels``
^^^^^^^^^^^^^^^

**Description:** use sample labels instead of indices

**Type:** boolean flag



``--morehelp``
^^^^^^^^^^^^^^

**Description:** prints full module list and descriptions

**Type:** boolean flag



``-o/--out``
^^^^^^^^^^^^

**Description:** Output MVF file

**Type:** file path; **Default:** None



``--overwrite``
^^^^^^^^^^^^^^^

**Description:** USE WITH CAUTION: force overwrite of outputs

**Type:** boolean flag



``-q/--quiet``
^^^^^^^^^^^^^^

**Description:** suppress progress meter

**Type:** boolean flag



``--test``
^^^^^^^^^^

**Description:** manually input a line for testing

**Type:** None; **Default:** None



``--test-nchar``
^^^^^^^^^^^^^^^^

**Description:** total number of samples for test string

**Type:** integer; **Default:** None



``-B/--linebuffer``
^^^^^^^^^^^^^^^^^^^

**Description:** number of lines to write at once to MVF

**Type:** integer; **Default:** 100000



``-V/--verbose``
^^^^^^^^^^^^^^^^

**Description:** report every line (for debugging)

**Type:** boolean flag


