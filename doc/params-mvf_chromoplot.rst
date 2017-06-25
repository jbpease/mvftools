.. mvf_chromoplot:

mvf_chromoplot
==============

Description
-----------

This program creates a chromoplot from an MVF alignment.
A chromoplot shows a genome-wide diagram of different
evolutionary histories for a given quartet of taxa.


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



``-s/--samples`` (required)
^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** 3 or more taxa to use for quartets

**Type:** None; **Default:** None



``-G/--outgroup`` (required)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:** 1 or more outgroups to use for quartets

**Type:** None; **Default:** None



``-c/--contigs``
^^^^^^^^^^^^^^^^

**Description:** Enter the ids of one or more contigs in the order they will appear in the chromoplot. (defaults to all ids in order present in MVF)

**Type:** None; **Default:** None



``-o/--outprefix``
^^^^^^^^^^^^^^^^^^

**Description:** Output prefix (not required).

**Type:** None; **Default:** None



``-q/--quiet``
^^^^^^^^^^^^^^

**Description:** suppress all output messages

**Type:** boolean flag



``-w/--windowsize``
^^^^^^^^^^^^^^^^^^^

**Description:** None

**Type:** integer; **Default:** 100000



``-x/--xscale``
^^^^^^^^^^^^^^^

**Description:** Width (in number of pixels) for each window

**Type:** integer; **Default:** 1



``-y/--yscale``
^^^^^^^^^^^^^^^

**Description:** Height (in number of pixels) for each track

**Type:** integer; **Default:** 20



``-C/--colors``
^^^^^^^^^^^^^^^

**Description:** three colors to use for chromoplot

**Type:** None; **Default:** None

**Choices:** {'lgrey': (250, 250, 250), 'dgrey': (192, 192, 192), 'black': (0, 0, 0), 'white': (255, 255, 255), 'red': (192, 0, 0), 'orange': (217, 95, 2), 'yellow': (192, 192, 0), 'green': (0, 192, 0), 'blue': (0, 0, 192), 'teal': (27, 158, 119), 'puce': (117, 112, 179), 'purple': (192, 0, 192), 'none': ()}


``-E/--emptymask``
^^^^^^^^^^^^^^^^^^

**Description:** Mask empty regions with this color.

**Type:** None; **Default:** none

**Choices:** {'lgrey': (250, 250, 250), 'dgrey': (192, 192, 192), 'black': (0, 0, 0), 'white': (255, 255, 255), 'red': (192, 0, 0), 'orange': (217, 95, 2), 'yellow': (192, 192, 0), 'green': (0, 192, 0), 'blue': (0, 0, 192), 'teal': (27, 158, 119), 'puce': (117, 112, 179), 'purple': (192, 0, 192), 'none': ()}


``-I/--infotrack``
^^^^^^^^^^^^^^^^^^

**Description:** Include an additional coverage information track that will show empty, uninformative, and informative loci. (Useful for ranscriptomes/RAD or other reduced sampling.

**Type:** boolean flag



``-M/--majority``
^^^^^^^^^^^^^^^^^

**Description:** Plot only 100% shading in the majority track  rather than shaded proportions in all tracks.

**Type:** boolean flag



``-P/--plottype``
^^^^^^^^^^^^^^^^^

**Description:** PNG image (default) or graph via matplotlib (experimental)

**Type:** None; **Default:** image

**Choices:** ['graph', 'image']

