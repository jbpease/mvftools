#########################################################
Examples of the same data in MVF Format and other formats
#########################################################

MVF Format
==========

::

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

FASTA Format
============
::
  
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

VCF Format
===========

::
  
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

