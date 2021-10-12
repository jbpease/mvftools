% MVFtools
  <br>![](logo.png)
%
% 2021-10-12


<link rel="stylesheet" type="text/css" media="all" href="mdmanstyle1.css">

---

***Version 0.6.2***

---

# Introduction

## Authors

>James B. Pease ([http://www.peaselab.org](http://www.peaselab.org))

>Benjamin K. Rosenzweig 

**Contributors:**

>Roddra Johnson

>Ellen Weinheimer


## What is MVFtools?
Multisample Variant Format (MVF) is designed for compact storage and efficient analysis of aligned multi-genome and multi-transcriptome datasets.  The programs provided in MVFtools support this format, including data conversion, filtering and transformation, and data analysis and visualization modules.  MVF format is designed specifically for biological data analysis, and sequence data is encoded based on the information content at a particular aligned sequence site.  This contextual encoding allows for rapid computation of phylogenetic and population genetic analyses, and small file sizes that enable unified and compact data sharing and distribution.


## How do I cite this?
Pease JB and BK Rosenzweig. 2018. "Encoding Data Using Biological Principles: the Multisample Variant Format for Phylogenomics and Population Genomics" *IEEE/ACM Transactions on Computational Biology and Bioinformatics*. 15(4):1231-1238.  [http://www.dx.doi.org/10.1109/tcbb.2015.2509997]()

Please also include the URL [https://www.github.com/peaselab/mvftools]() in your methods section where the program is referenced.

(Note this paper was originally published online in 2015, but did not receive final citation page numbering until 2018.  You may see older citations as 2015, which is the same paper.)

---

# Getting Started

## Requirements
* Python 3.x: [https://www.python.org/downloads/]()
* Biopython 1.6+: [http://www.biopython.org/]()
* Scipy: [http://www.scipy.org/]()
* Numpy [http://www.numpy.org/]()
* Matplotlib [http://www.matplotlib.org/]()

### Additional Requirements for Some Modules
* RAxML-ng; [https://github.com/amkozlov/raxml-ng]())
* RAxML: 8.2+.X; [https://sco.h-its.org/exelixis/web/software/raxml/index.html]())
* PAML: [http://abacus.gene.ucl.ac.uk/software/paml.html]()

## Installation
No installation is required, mvftools scripts should work as long as Python3 is installed.  The repository can be cloned or downloaded as a .zip file from GitHub.

```
git clone https://www.github.com/jbpease/mvftools
```
---

## Preparing your data

**Sequence Alignment**

MVF files can be created from VCF, FASTA, and MAF files using the `ConvertVCF2MVF`, `ConvertFasta2MVF`, or `ConvertMAF2MVF` commands, respectively.  Once converted to MVF format, analyses and manipulations can be carried out using the rest of the commands in MVFtools.


## Basic usage examples

**Case #1: Generate phylogenies from 100kb windows using a VCF data**
```
python3 mvftools.py ConvertVCF2MVF --vcf DATA.vcf --mvf DATA.mvf
python3 mvftools.py InferWindowTree --mvf DATA.mvf --out WINDOWTREES.txt \
--windowsize 100000
```
**Case #2: Convert a large FASTA file, then generate window-based counts for DFOIL/D-statistic introgression testing from the first five samples**::
```
python3 mvftools.py ConvertFasta2MVF --fasta DATA.fasta --mvf DATA.mvf
python3 mvftools.py CalcPatternCount --mvf DATA.mvf --out PATTERNS.txt \ 
--windowsize 100000 --samples 0,1,2,3,4
```
The file is now ready to use as an input file for with dfoil:
([http://www.github.com/jbpease/dfoil](http://www.github.com/jbpease/dfoil)).

**Case #3: Convert a VCF file, then generate window-based counts for DFOIL/D-statistic introgression testing from the first five samples**::
```
python3 mvftools.py ConvertVCF2MVF --vcf INPUT.vcf --mvf DATA.mvf
python3 mvftools.py CalcPatternCount --mvf DATA.mvf --out PATTERNS.txt \ 
--windowsize 100000 --samples 0,1,2,3,4
```
The file is now ready to use as an input file for with dfoil:
([http://www.github.com/jbpease/dfoil](http://www.github.com/jbpease/dfoil)).


---

# MVF Format Specification (version 1.2)

## MVF General Notes and Usage

### General Features
MVF is primarily intended for site-wise analyses in phylogenomics and population genomics. MVF is formatted to contain one aligned site per line, but contains only allelic information, therefore MVF most closely mimics VCF files in formatting, but resembles MAF format in informational content,  Additionally, MVF uses special formatting to lower file sizes and speed up filtering and analysis.  MVF can readily be adapted from other common sequence formats including VCF, FSATA, and MAF.  MVF is also designed to be able to accommodate readily store other information for phylogenomic projects, including tree topologies and sample metadata.

### Native Gzip read/write

MVF is designed to work natively with GZIP compression and uses a formatting that attempts to strike a balance between fast filtering, easy visual inspection, while using character patterns that create a good Gzip compression ratio. As long as any input or output file path ends with exactly ".gz", all MVF scripts will natively read/write to gzip-compressed files.

### General Notes on Filtering

MVF was specifically designed as a "vertical" format for rapid filtering of *sites* in large-scale phylogenomic analyses. (rather than being "horizontal" to visually show alignment) Therefore, the following should be noted to take advantage of MVF formatting for rapid filtering (i.e. with grep/zgrep).

* `#` is present iff. the line is in the header
* `@` is present iff. the position is non-reference
* `X` is present in the allele string iff. the positon has ambiguity data
* `#:` can quickly filter by chromosome
* `:#` can quickly filter by coordinate numbers
* Allele strings with one or two characters have full sample coverage (no gaps)
* Allele strings with `@[any]+` have coverage=1, `[not@][any]+` have coverage=2 
* One or two-character allele strings, or notation with `[any]+` CANNOT contain homoplasy or synapomorphy (by definition).

## Header Specification

All header lines begin with one or more `#` and contain single-space separated fields.

### MVF declaration line
First header line always starts with `##mvf`, followed by required metadata fields:

* version=1.2
* mvftype=[dna, protein, codon]
     
and optionally:

* an arbitrary number of metadata fields in `key=value` format (`mvftype` and `version` not allowed as key)

### Sample information
Sample information (columns) header lines are specified by:

* line starts with `#s` ("s" for sample) with no leading spaces
* LABEL (must be unique, no spaces)
* an arbitrary number of metadata fields in key=value format ('label' not allowed as key)

The first entry should be the reference sequence (if aligned to reference) or can be any sequence in the case of non-reference-aligned de novo alignment).

### Contig information

Contig information header lines are specified by:

* line starts with `#c` ("c" for contig)
* CONTIG_ID (must be unique, alpha-numeric strong recommended, must not contain `*:;,@!+` or spaces)
* `label=[NAME]` (recommended by not required to be unique, no spaces allowed)
* `len=[LENGTH]` (integer > 0, or zero for unknown)
* `ref=[0/1]`, indicates if contig is reference-based (=1) or not (=0)
* an arbitrary number of metadata fields in key=value format ("label", "len", and "ref" not allowed as key)

### Tree information
Tree information may (optionally) be specified in header lines by:

* line starts with `#t` ("t" for tree/topology)
* `TREE_ID=[###]` (must be unique, alpha-numeric)
* `TOPOLOGY=[tree_String]` in Newick/Phylip/parenthetical format (must end with ';')
* an arbitrary number of metadata fields in key=value format

To take full advantage of MVF tree storage, use the same sample labels as in the `#s` header lines
	
### Notes
General project notes may (optionally) be specified in the header lines by:

* line starts with `#n` ("n" for notes)
* Text is unstructured and is not necessarily formatted as metadata
	
### Example Header

```
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
```

## Entry Specification

***Note: all examples show an MVF entry with REF and four samples.***

Entries are structured as two space-separated columns:

`ID:POSITION	ALLELES [ALLELES ALLELES ...]`

  * `ID:POSITION` = chromosomal id matching the first element of a contig in the `#c` header element
  * `POSITION` = 1-based position on the contig with matching `CONTIG_ID`
  * `ALLELES` = one or more records of alleles at reference-based location specified by `ID:POSITION` and matching the formatting below

### For mvftype=codon
* Allele columns are `PROTEIN DNA1 DNA2 DNA3` where the three DNA columns represent three codon positions in collated form
* Position is the position of the lowest numbered codon position (regardless of transcript strand) and `DNA1/2/3` codon columns are given in order to match the protein (again regardless of transcript orientation)

## Allele formatting

***Note: all examples show an MVF entry with five samples.***

For reference-anchored contigs, the first allele is assumed to be the "reference" allele by default. Each entry must either (1) contain the same number of characters as sample labels specified in the header or (2) use one of the special cases in the section below.

`ATCTG` =  (REF is 'A' samples 1&3 are 'T', sample 2 is 'C', sample 4 is 'G')

## Special cases

### Invariant sites

When all alleles are both present (non-gap) and all the same, this is represented by a single base.

  `A = AAAAA`

### Monoallelic non-reference samples 

When all alleles in the samples (non-REF) are the same but differ from REF, this is represented by two bases.

  `AT = ATTTT`

  `Aa = Aaaaa`

### Single-variant sites

When only one of the samples varies from the others, this is specified as:

  `[reference_base, majority_base, "+", unique_base, unique_position]`

This is useful shorthand for both sites with one a single base that differs and samples with only one sample represented.  When the site only has coverage via one sample (i.e. all other bases are empty, the '-' is omitted from the second position.

  `AC+T2 = ACTCC`

  `AA+C2 = AACAA`

  `-+A2  = --A--`

  `A+A2  = A-A--`

  `A+a2  = A-a--`

  `A+C2  = A-C--`

### Non-reference aligned sites 
Added in MVF v.1.2, this facilitates using MVF for non-reference aligned sequences (e.g. aligned sets of orthologs from de novo assembled transcripts). These non-reference-anchored alignments can comprise the entire MVF file or be included in addition to reference-aligned contigs. Non-reference-contigs in their header entry should include the keyword "nonref" (see Section 1.3). Contigs labels and coordinates are labelled the same as reference-based entries. To denote that the sequence is non-reference and not simply a deletion in the reference, the character "@" should be the first character of the alignment.  In the case an entirely non-reference MVF, all contigs can be labelled as "nonref," but one sequence should be chosen as the reference for the purposes of the allele
string.  When this sequence is not present, `@` is still used.

  `@AATT   = -AATT`

  `@A+T3   = -A-T-`

  `@-+A3   = ---A-`

## Character encoding

### Nucleotide Notation

* Standard IUPAC nucleotide codes are used: `ACGT`, and `U` for uracil in RNA
* Standard IUPAC bialleic ambiguity codes `KMRSWY` are used also.
* <span style="color:blue;">Standard IUPAC triallelic ambiguity codes (`BDHV`) are allowed in convsersion from FASTA and MAF, or when a polyploid VCF is converted.  For diploid VCFs, triallelic sites are converted to ambiguous (`X`) instead. (*Changed in 1.2.1*)</span>
* Current MVF formatting does NOT recognize rare symbols (`ISOX`, or `Phi`)
* Ambiguous nucleotide is denoted by `X` instead of standard `N` 
  
### Amino Acid Notation

* Standard IUPAC amino acid codes are used: `ACDEFGHIKLMNPQRSTVWY`
* Standard stop codon symbol `*` is used
* Currently the ambiguous/rare symbols are not recognized (`BZ`)

### Use of `X` for ambiguous nucleotides and amino acids

In standard notation, "`N`" is used for an ambiguous nucleotide, which could be any of A/C/G/T.  
However, in amino acid notation `N` stands for "Asparagine" and is a valid character, while `X` is used for an ambiguous amino acid.

The MVF Standard since v1.2 adopts `X` as unified ambiguity character for both nucleotides and proteins for MVF files for two purposes:

1. To creates a unified ambiguity character for MVF codon files for faster processing

2. To allow fast filtering of ambiguous lines

Also note that while "`X`" in expanded IUPAC notation refers to "xanthosine," MVF currently does not support rare nucleotides.

***NOTE: In all conversion utilities that export from MVF format to another file format conversion to the standard "N"/"X" for ambiguous nucleotides/amino acids should be implemented.***

---

## Examples of the same data in MVF Format and other formats

### MVF Format

```
  ##mvf sourceformat=fasta version=1.2 mvftype=dna ncol=5 
  #s Hsapiens
  #s Ptroglodytes 
  #s Ppaniscus
  #s Ggorilla 
  #s Mmusculus
  #c 1 label=Chromosome1 length=248956422
  #n Note: This is an example file showing data formatting 
  1:100 A
  1:101 A
  1:102 A
  1:103 T
  1:104 TT+C4
  1:105 GC
  1:106 A+A4
  1:107 AATTA
  1:108 AC+G4
```
### FASTA Format
``` 
  >Hsapiens gi:1234 geneid:GeneOfInterest chrom:1 start:100 end:108
  AAATTGAAA
  >Ptroglodytes geneid:GeneOfInterest 
  AAATTC-AC
  >Ppaniscus geneid:GeneOfInterest 
  AAATTC-TC
  >Ggorilla geneid:GeneOfInterest
  AAATTC-TC
  >Mmusculus geneid:GeneOfInterest 
  AAATCCAAG
```
### VCF Format
```  
  ##fileformat=VCFv4.1
  ##samtoolsVersion=0.1.19-44428cd
  ##reference=hg19.fa
  ##contig=<ID=Chromosome1,length=248956422>
  ##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw read depth">
  ##INFO=<ID=DP4,Number=4,Type=Integer,Description="# high-quality ref-forward bases, ref-reverse, alt-forward and alt-reverse bases">
  ##INFO=<ID=MQ,Number=1,Type=Integer,Description="Root-mean-square mapping quality of covering reads">
  ##INFO=<ID=FQ,Number=1,Type=Float,Description="Phred probability of all samples being the same">
  ##INFO=<ID=AF1,Number=1,Type=Float,Description="Max-likelihood estimate of the first ALT allele frequency (assuming HWE)">
  ##INFO=<ID=AC1,Number=1,Type=Float,Description="Max-likelihood estimate of the first ALT allele count (no HWE assumption)">
  ##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
  ##INFO=<ID=IS,Number=2,Type=Float,Description="Maximum number of reads supporting an indel and fraction of indel reads">
  ##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes for each ALT allele, in the same order as listed">
  ##INFO=<ID=G3,Number=3,Type=Float,Description="ML estimate of genotype frequencies">
  ##INFO=<ID=HWE,Number=1,Type=Float,Description="Chi^2 based HWE test P-value based on G3">
  ##INFO=<ID=CLR,Number=1,Type=Integer,Description="Log ratio of genotype likelihoods with and without the constraint">
  ##INFO=<ID=UGT,Number=1,Type=String,Description="The most probable unconstrained genotype configuration in the trio">
  ##INFO=<ID=CGT,Number=1,Type=String,Description="The most probable constrained genotype configuration in the trio">
  ##INFO=<ID=PV4,Number=4,Type=Float,Description="P-values for strand bias, baseQ bias, mapQ bias and tail distance bias">
  ##INFO=<ID=INDEL,Number=0,Type=Flag,Description="Indicates that the variant is an INDEL.">
  ##INFO=<ID=PC2,Number=2,Type=Integer,Description="Phred probability of the nonRef allele frequency in group1 samples being larger (,smaller) than in group2.">
  ##INFO=<ID=PCHI2,Number=1,Type=Float,Description="Posterior weighted chi^2 P-value for testing the association between group1 and group2 samples.">
  ##INFO=<ID=QCHI2,Number=1,Type=Integer,Description="Phred scaled PCHI2.">
  ##INFO=<ID=PR,Number=1,Type=Integer,Description="# permutations yielding a smaller PCHI2.">
  ##INFO=<ID=QBD,Number=1,Type=Float,Description="Quality by Depth: QUAL/#reads">
  ##INFO=<ID=RPB,Number=1,Type=Float,Description="Read Position Bias">
  ##INFO=<ID=MDV,Number=1,Type=Integer,Description="Maximum number of high-quality nonRef reads in samples">
  ##INFO=<ID=VDB,Number=1,Type=Float,Description="Variant Distance Bias (v2) for filtering splice-site artefacts in RNA-seq data. Note: this version may be broken.">
  ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
  ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
  ##FORMAT=<ID=GL,Number=3,Type=Float,Description="Likelihoods for RR,RA,AA genotypes (R=ref,A=alt)">
  ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="# high-quality bases">
  ##FORMAT=<ID=DV,Number=1,Type=Integer,Description="# high-quality non-reference bases">
  ##FORMAT=<ID=SP,Number=1,Type=Integer,Description="Phred-scaled strand bias P-value">
  ##FORMAT=<ID=PL,Number=G,Type=Integer,Description="List of Phred-scaled genotype likelihoods">
  #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Ptroglodytes Ppaniscus	Ggorilla	Mmusculus
  ch01	100	.	A	.	30	.	DP=5;AF1=0;AC1=0;DP4=5,0,0,0;MQ=20;FQ=-23.4	PL:DP	0/0:0,6,40:2:4	0/0:0,6,40:2:4	0/0:0,6,40:2:4	0/0:0,6,40:2:4
  ch01	101	.	A	.	30	.	DP=5;AF1=0;AC1=0;DP4=5,0,0,0;MQ=20;FQ=-23.4	PL:DP	0/0:0,6,40:2:4	0/0:0,6,40:2:4	0/0:0,6,40:2:4	0/0:0,6,40:2:4
  ch01	102	.	A	.	30	.	DP=5;AF1=0;AC1=0;DP4=5,0,0,0;MQ=20;FQ=-23.4	PL:DP	0/0:0,6,40:2:4	0/0:0,6,40:2:4	0/0:0,6,40:2:4	0/0:0,6,40:2:4
  ch01	103	.	T	.	32	.	DP=5;AF1=0;AC1=0;DP4=5,0,0,0;MQ=20;FQ=-23.4	PL:DP	0/0:0,6,40:2:4	0/0:0,6,40:2:4	0/0:0,6,40:2:4	0/0:0,6,40:2:4
  ch01	104	.	T	C	7.61	.	DP=2;VDB=6.720000e-02;AF1=1;AC1=58;DP4=0,0,1,1;MQ=20;FQ=-23.8	GT:PL:DP:GQ	0/0:0,6,40:2:4	0/0:0,6,40:2:4	0/0:0,6,40:2:4	1/1:38,6,0:2:4
  ch01	105	.	G	C	32.1	.	DP=5;AF1=0;AC1=0;DP4=5,0,0,0;MQ=20;FQ=-23.4	PL:DP	0/0:0,6,40:2:4	0/0:0,6,40:2:4	0/0:0,6,40:2:4	1/1:38,6,0:2:4
  ch01	106	.	A	.	30	.	DP=5;AF1=0;AC1=0;DP4=5,0,0,0;MQ=20;FQ=-23.4	PL:DP	0:0	0:0	0:0	0/0:0,6,40:2:4
  ch01	107	.	A	T	24.4	.	DP=5;AF1=1;AC1=58;DP4=0,0,1,0;MQ=20;FQ=-23.4	PL:DP	0/0:0,6,40:2:4	1/1:38,6,0:2:4	1/1:38,6,0:2:4	0/0:0,6,40:2:4
  ch01	108	.	A	C,G	999	.	DP=52;VDB=6.361343e-02;RPB=-1.264051e+00;AF1=0.9325;AC1=54;DP4=0,2,20,26;MQ=20;FQ=-16.1;PV4=0.5,1,1,1	GT:PL:DP:GQ	1/1:20,3,0,20,3,20:1:11	1/1:36,6,0,36,6,36:2:13	1/1:36,6,0,36,6,36:2:13	1/1:95,95,95,18,18,0:6:8
```

---

# Program Parameters

