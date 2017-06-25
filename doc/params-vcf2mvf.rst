.. vcf2mvf:

vcf2mvf
=======

Description
-----------

MVFtools: Multisample Variant Format Toolkit
James B. Pease and Ben K. Rosenzweig
http://www.github.org/jbpease/mvftools


Parameters
----------

``-h/--help``
^^^^^^^^^^^^^

**Description:** show this help message and exit

**Type:** boolean flag



``--out`` (required)
^^^^^^^^^^^^^^^^^^^^

**Description:** output MVF file

**Type:** None; **Default:** None



``--vcf`` (required)
^^^^^^^^^^^^^^^^^^^^

**Description:** input VCF file

**Type:** file path; **Default:** None



``--allelesfrom``
^^^^^^^^^^^^^^^^^

**Description:** get additional alignment columns
                                from INFO fields (:-separated)

**Type:** None; **Default:** None



``--contigids``
^^^^^^^^^^^^^^^

**Description:** manually specify one or more contig ids
                                 as ID;VCFLABE;MVFLABEL, note that
                                 VCFLABEL must match EXACTLY the contig string
                                 labels in the VCF file

**Type:** None; **Default:** None



``--fieldsep``
^^^^^^^^^^^^^^

**Description:** VCF field separator (default='TAB')

**Type:** None; **Default:** TAB

**Choices:** ['TAB', 'SPACE', 'DBLSPACE', 'COMMA', 'MIXED']


``--linebuffer``
^^^^^^^^^^^^^^^^

**Description:** number of lines to hold in read/write buffer

**Type:** integer; **Default:** 100000



``--lowdepth``
^^^^^^^^^^^^^^

**Description:** below this read depth coverage, convert to lower case set to 0 to disable

**Type:** integer; **Default:** 3



``--lowqual``
^^^^^^^^^^^^^

**Description:** below this quality convert to lower case
                                set to 0 to disable

**Type:** integer; **Default:** 20



``--maskdepth``
^^^^^^^^^^^^^^^

**Description:** below this read depth mask with N/n

**Type:** integer; **Default:** 1



``--maskqual``
^^^^^^^^^^^^^^

**Description:** low quality cutoff, bases replaced by N/-
                             set to 0 to disable

**Type:** integer; **Default:** 3



``--no_autoindex``
^^^^^^^^^^^^^^^^^^

**Description:** do not automatically index contigs from the VCF

**Type:** boolean flag



``--outflavor``
^^^^^^^^^^^^^^^

**Description:** choose output MVF flavor to include quality scores and/or indels

**Type:** None; **Default:** dna

**Choices:** ['dna', 'dnaqual', 'dnaqual-indel', 'dna-indel']


``--overwrite``
^^^^^^^^^^^^^^^

**Description:** USE WITH CAUTION: force overwrite of outputs

**Type:** boolean flag



``--qual``
^^^^^^^^^^

**Description:** Include Phred genotype quality (GQ) scores

**Type:** boolean flag



``--reflabel``
^^^^^^^^^^^^^^

**Description:** label for reference sample (default='REF')

**Type:** None; **Default:** REF



``--samplereplace``
^^^^^^^^^^^^^^^^^^^

**Description:** one or more TAG:NEWLABEL or TAG, items,
                                if TAG found in sample label, replace with
                                NEW (or TAG if NEW not specified)
                                NEW and TAG must each be unique

**Type:** None; **Default:** None


