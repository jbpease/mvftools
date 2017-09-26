======================================
MVF Format Specification (version 1.2)
======================================

Version History
===============

v1.1.1
------
Codons and Proteins accommodated

v1.2
----
Dot masking, multi-line header, adoption of "X" in place of "N" for nucleotides, support for non-reference aligned sequences.

MVF General Notes and Usage
===========================

General Features
----------------
MVF is primarily intended for site-wise analyses in phylogenomics and population genomics. MVF is formatted to contain one aligned site per line, but contains only allelic information, therefore MVF most closely mimics VCF files in formatting, but resembles MAF format in informational content,  Additionally, MVF uses special formatting to lower file sizes and speed up filtering and analysis.  MVF can readily be adapted from other common sequence formats including VCF, FSATA, and MAF.  MVF is also designed to be able to accommodate readily store other information for phylogenomic projects, including tree topologies and sample metadata.

Native Gzip read/write
----------------------

MVF is designed to work natively with GZIP compression and uses a formatting that attempts to strike a balance between fast filtering, easy visual inspection, while using character patterns that create a good Gzip compression ratio. As long as any input or output file path ends with exactly ".gz", all MVF scripts will natively read/write to gzip-compressed files.

General Notes on Filtering
--------------------------

MVF was specifically designed as a "vertical" format for rapid filtering of *sites* in large-scale phylogenomic analyses. (rather than being "horizontal" to visually show alignment) Therefore, the following should be noted to take advantage of MVF formatting for rapid filtering (i.e. with grep/zgrep).

* ``#`` is present iff. the line is in the header
* ``@`` is present iff. the position is non-reference
* ``X`` is present in the allele string iff. the positon has ambiguity data
* ``#:`` can quickly filter by chromosome
* ``:#`` can quickly filter by coordinate numbers
* Allele strings with one or two characters have full sample coverage (no gaps)
* Allele strings with ``@[any]+`` have coverage=1, ``[not@][any]+`` have coverage=2 
* One or two-character allele strings, or notation with ``[any]+`` CANNOT contain homoplasy or synapomorphy (by definition).

Header Specification
====================

All header lines begin with one or more ``#`` and contain single-space separated fields.

MVF declaration line
--------------------
First header line always starts with ``##mvf``, followed by required metadata fields:

   * version=1.2
   * mvftype=[dna, protein, codon]
     
and optionally:

   * an arbitrary number of metadata fields in key=value format ('mvftype' and 'version' not allowed as key)

Sample information
------------------
Sample information (columns) header lines are specified by:

* line starts with ``#s`` ("s" for sample) with no leading spaces
* LABEL (must be unique, no spaces)
* an arbitrary number of metadata fields in key=value format ('label' not allowed as key)

The first entry should be the reference sequence (if aligned to reference) or can be any sequence in the case of non-reference-aligned de novo alignment).

Contig information
------------------

Contig information header lines are specified by:

* line starts with ``#c`` ("c" for contig)
* CONTIG_ID (must be unique, alpha-numeric strong recommended, must not contain ``*:;,@!+`` or spaces)
* ``label=[NAME]`` (recommended by not required to be unique, no spaces allowed)
* ``len=[LENGTH]`` (integer > 0, or zero for unknown)
* ``ref=[0/1]``, indicates if contig is reference-based (=1) or not (=0)
* an arbitrary number of metadata fields in key=value format ("label", "len", and "ref" not allowed as key)

Tree information
----------------
Tree information may (optionally) be specified in header lines by:

* line starts with ``#t`` ("t" for tree/topology)
* ``TREE_ID=[###]` (must be unique, alpha-numeric)
* ``TOPOLOGY=[tree_String]`` in Newick/Phylip/parenthetical format (must end with ';')
* an arbitrary number of metadata fields in key=value format
To take full advantage of MVF tree storage, use the same sample labels as in the ``#s`` header lines
	
Notes
-----
General project notes may (optionally) be specified in the header lines by:

* line starts with ``#n`` ("n" for notes)
* Text is unstructured and is not necessarily formatted as metadata
	
Example Header
--------------
::
  ##mvf version=1.2 mvftype=[MVFTYPE]
  #s SAMPLE0 meta0=somevalue meta1=0 ...
  #s SAMPLE1 meta0=somethingele meta1=1 ...
  #s SAMPLE2 meta0=somesome meta1=0 ...
  ...
  #c 0 label=CONTIG0 length=100 ref=1 meta0=somevalue ...
  #c 1 label=CONTIG1 length=200 ref=0 meta0=someother ...
  ...
  #t 0 ((SAMPLE0,SAMPLE1),SAMPLE2); model=GTRGAMMA software=RAxML
  #t 1 ((SAMPLE2,SAMPLE0),SAMPLE1); model=GTRGAMMA software=RAxML partition=chrom1
  ...
  #n Notes on this project.


Entry Specification
===================

.. note:: all examples show an MVF entry with REF and four samples

Entries are structured as two space-separated columns:

``ID:POSITION	ALLELES [ALLELES ALLELES ...]``

  * ``ID:POSITION`` = chromosomal id matching the first element of a contig in the ``#c`` header element
  * ``POSITION`` = 1-based position on the contig with matching ``CONTIG_ID``
  * ``ALLELES`` = one or more records of alleles at reference-based location specified by ``ID:POSITION`` and matching the formatting below

For mvftype=codon
-----------------
* Allele columns are ``PROTEIN DNA1 DNA2 DNA3`` where the three DNA columns represent three codon positions in collated form
* Position is the position of the lowest numbered codon position (regardless of transcript strand) and ``DNA1/2/3`` codon columns are given in order to match the protein (again regardless of transcript orientation)

Allele formatting
-----------------

.. note:: all examples show an MVF entry with five samples.

For reference-anchored contigs, the first allele is assumed to be the "reference" allele by default. Each entry must either (1) contain the same number of characters as sample labels specified in the header or (2) use one of the special cases in the section below.

``ATCTG`` =  (REF is 'A' samples 1&3 are 'T', sample 2 is 'C', sample 4 is 'G')

Special cases
-------------

.. note:: all examples show an MVF entry with five samples

Invariant sites
---------------

When all alleles are both present (non-gap) and all the same, this is represented by a single base.

  ``A = AAAAA``

Monoallelic non-reference samples 
---------------------------------

When all alleles in the samples (non-REF) are the same but differ from REF, this is represented by two bases.

  ``AT = ATTTT``
  ``Aa = Aaaaa``

Single-variant sites
--------------------

When only one of the samples varies from the others, this is specified as:

::

  [reference_base, majority_base, "+", unique_base, unique_position]

This is useful shorthand for both sites with one a single base that differs and samples with only one sample represented.  When the site only has coverage via one sample (i.e. all other bases are empty, the '-' is omitted from the second position.

  ``AC+T2 = ACTCC``
  ``AA+C2 = AACAA``
  ``-+A2  = --A--``
  ``A+A2  = A-A--``
  ``A+a2  = A-a--``
  ``A+C2  = A-C--``

Non-reference aligned sites 
---------------------------
Added in MVF v.1.2, this facilitates using MVF for non-reference aligned sequences (e.g. aligned sets of orthologs from de novo assembled transcripts). These non-reference-anchored alignments can comprise the entire MVF file or be included in addition to reference-aligned contigs. Non-reference-contigs in their header entry should include the keyword "nonref" (see Section 1.3). Contigs labels and coordinates are labelled the same as reference-based entries. To denote that the sequence is non-reference and not simply a deletion in the reference, the character "@" should be the first character of the alignment.  In the case an entirely non-reference MVF, all contigs can be labelled as "nonref," but one sequence should be chosen as the reference for the purposes of the allele
string.  When this sequence is not present, ``@`` is still used.

  ``@AATT   = -AATT``
  ``@A+T3   = -A-T-``
  ``@-+A3   = ---A-``

Character encoding
==================

Nucleotide Notation
-------------------

* Standard IUPAC nucleotide codes are used: ``ACGT``, and ``U`` for uracil in RNA
* Standard IUPAC bialleic ambiguity codes ``KMRSWY`` are used also.
* Current MVF formatting does NOT allow triallelic ambiguity codes (``BDHV``), which are converted to ambiguous (``X``) instead.
* Current MVF formatting does NOT recognize rare symbols (``ISOX``, or ``Phi``)
* Ambiguous nucleotide is denoted by ``X`` instead of standard ``N`` 
  
Amino Acid Notation
-------------------

* Standard IUPAC amino acid codes are used: ``ACDEFGHIKLMNPQRSTVWY``
* Standard stop codon symbol ``*`` is used
* Currently the ambiguous/rare symbols are not recognized (``BZ``)

Use of ``X`` for ambiguous nucleotides and amino acids
------------------------------------------------------

In standard notation, "``N``" is used for an ambiguous nucleotide, which could be any of A/C/G/T.  
However, in amino acid notation ``N`` stands for "Asparagine" and is a valid character, while ``X`` is used for an ambiguous amino acid.
MVF v1.2 adopts ``X`` as unified ambiguity character for both nucleotides and proteins for MVF files for two purposes:
1. To creates a unified ambiguity character for MVF codon files for faster processing
2. To allow fast filtering of ambiguous lines
Also note that while 'X' in expanded IUPAC notation refers to 'xanthosine,' MVF currently does not support rare nucleotides.
.. note:: In all conversion utilities that export from MVF format to another file format conversion to the standard "N"/"X" for ambiguous nucleotides/amino acids should ALWAYS be implemented.
