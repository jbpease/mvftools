---

## CalcAllCharacterCountPerSample
***Calculates the count of different character typesin an MVF file***

**Parameters**


``--mvf`` (required) = Input MVF file. (type=file path, default=None)

``--out`` (required) = Output file (type=file path, default=None)

``--contig-ids/--contigids`` = Specify comma-separated list of contig short ids. Must match exactly. Do not use with --contig-labels. (type=None, default=None)

``--contig-labels/--contiglabels`` = Specify comma-separated list of contig full labels. Must match exactly. Do not use with --contig-ids (type=None, default=None)

``--mincoverage`` = Mininum sample coverage for sites. (type=integer, default=None)


``--quiet`` = Suppress screen output. (flag, default=False)


``--sample-indices/--sampleindices`` = Specify comma-separated list of sample numerical indices (first sample is 0). Leave blank for all samples. Do not use with --sample_labels. (type=None, default=None)

``--sample-labels`` = Specify comma-separated list of sample labels. Labels must be exact (case-sensitive). Leave blank for all samples.Do not use with --sample_indicies. (type=None, default=None)


``--windowsize`` = Set integer window size. Use 0 for whole file. Use -1 for whole contigs.  (flag, default=100000)


---

## CalcCharacterCount
***Calculates the count of different character typesin an MVF file***

**Parameters**


``--mvf`` (required) = Input MVF file. (type=file path, default=None)

``--out`` (required) = Output file (type=file path, default=None)

``--base-match/--basematch`` = String of bases to match (i.e. numerator). (type=None, default=None)

``--base-total/--basetotal`` = String of bases for total (i.e. denominator). (type=None, default=None)

``--contig-ids/--contigids`` = Specify comma-separated list of contig short ids. Must match exactly. Do not use with --contig-labels. (type=None, default=None)

``--contig-labels/--contiglabels`` = Specify comma-separated list of contig full labels. Must match exactly. Do not use with --contig-ids (type=None, default=None)

``--mincoverage`` = Mininum sample coverage for sites. (type=integer, default=None)


``--quiet`` = Suppress screen output. (flag, default=False)


``--sample-indices/--sampleindices`` = Specify comma-separated list of sample numerical indices (first sample is 0). Leave blank for all samples. Do not use with --sample_labels. (type=None, default=None)

``--sample-labels`` = Specify comma-separated list of sample labels. Labels must be exact (case-sensitive). Leave blank for all samples.Do not use with --sample_indicies. (type=None, default=None)


``--windowsize`` = Set integer window size. Use 0 for whole file. Use -1 for whole contigs.  (flag, default=100000)


---

## CalcDstatCombinations
***Calculates all D-statistics for all combinations ofspecified taxa in an MVF file.***

**Parameters**


``--mvf`` (required) = Input MVF file. (type=file path, default=None)

``--out`` (required) = Output file (type=file path, default=None)

``--contig-ids/--contigids`` = Specify comma-separated list of contig short ids. Must match exactly. Do not use with --contig-labels. (type=None, default=None)

``--contig-labels/--contiglabels`` = Specify comma-separated list of contig full labels. Must match exactly. Do not use with --contig-ids (type=None, default=None)

``--outgroup-indices/--outgroupindices`` = Specify comma-separated list of outgroup sample numerical indices (first column is 0). Leave blank for all samples. Do not use with --outgroup_labels. (type=None, default=None)

``--outgroup-labels/--outgrouplabels`` = Specify comma-separated list of outgroup sample labels. Labels must be exact (case-sensitive). Leave blank for all samples.Do not use with --outgroup_indicies. (type=None, default=None)


``--quiet`` = Suppress screen output. (flag, default=False)


``--sample-indices/--sampleindices`` = Specify comma-separated list of 3 or more sample numerical indices (first sample is 0). Leave blank for all samples. Do not use with --sample_labels. (type=None, default=None)

``--sample-labels`` = Specify comma-separated list of 3 or more sample labels. Labels must be exact (case-sensitive). Leave blank for all samples.Do not use with --sample_indicies. (type=None, default=None)

---

## CalcPairwiseDistances
***Calculates pairwise sequence distances for combinations ofspecified taxa in an MVF file.***

**Parameters**


``--mvf`` (required) = Input MVF file. (type=file path, default=None)

``--out`` (required) = Output file (type=file path, default=None)

``--ambig`` = By default, ambiguous nucleotides are excluded.  This option will include sets of ambiguous characters by randomly choosing one of the options for: RYMKWS ('random2') or RYMKWS+BDHV ('random3') (type=None, default=None)
Choices: ('random2', 'random3')

``--data-type/--datatype`` = Data type to compare.(This option is only needed for codon  MVF files, others will default.) (type=None, default=None)
Choices: ('dna', 'prot')


``--emit-counts`` = output additional file that presents the raw counts of pairwise patterns for each sample pair tested for each window (flag, default=False)


``--mincoverage`` = Mininum sample coverage for sites. (type=integer, default=None)


``--quiet`` = Suppress screen output. (flag, default=False)


``--sample-indices/--sampleindices`` = Specify comma-separated list of 2 or more sample numerical indices (first sample is 0). Leave blank for all samples. Do not use with --sample_labels. (type=None, default=None)

``--sample-labels`` = Specify comma-separated list of 2 or more sample labels. Labels must be exact (case-sensitive). Leave blank for all samples.Do not use with --sample_indicies. (type=None, default=None)


``--windowsize`` = Set integer window size. Use 0 for whole file. Use -1 for whole contigs.  (flag, default=100000)


---

## CalcPatternCount
***Counts biallelic site pattersn (AB-patterns) forspecified combinations of taxa in an MVF file.***

**Parameters**


``--mvf`` (required) = Input MVF file. (type=file path, default=None)

``--out`` (required) = Output file (type=file path, default=None)

``--mincoverage`` = Mininum sample coverage for sites. (type=integer, default=None)


``--output-lists`` = None (flag, default=False)



``--quiet`` = Suppress screen output. (flag, default=False)


``--sample-indices/--sampleindices`` = Specify comma-separated list of sample numerical indices (first sample is 0). Leave blank for all samples. Do not use with --sample_labels. (type=None, default=None)

``--sample-labels`` = Specify comma-separated list of sample labels. Labels must be exact (case-sensitive). Leave blank for all samples.Do not use with --sample_indicies. (type=None, default=None)


``--windowsize`` = Set integer window size. Use 0 for whole file. Use -1 for whole contigs.  (flag, default=100000)


---

## CalcSampleCoverage
***Counts per-contig coverage forspecified sample columns in an MVF file.***

**Parameters**


``--mvf`` (required) = Input MVF file. (type=file path, default=None)

``--out`` (required) = Output file (type=file path, default=None)

``--contig-ids/--contigids`` = Specify comma-separated list of contig short ids. Must match exactly. Do not use with --contig-labels. (type=None, default=None)

``--contig-labels/--contiglabels`` = Specify comma-separated list of contig full labels. Must match exactly. Do not use with --contig-ids (type=None, default=None)


``--quiet`` = Suppress screen output. (flag, default=False)


``--sample-indices/--sampleindices`` = Specify comma-separated list of sample numerical indices (first sample is 0). Leave blank for all samples. Do not use with --sample_labels. (type=None, default=None)

``--sample-labels`` = Specify comma-separated list of sample labels. Labels must be exact (case-sensitive). Leave blank for all samples.Do not use with --sample_indicies. (type=None, default=None)

---

## ConcatenateMVF
***Combine non-overlapping contigs from one or more MVF files into asingle MVF file.  This does NOT merge columns.Use MergeMVF to merge sample columns from multiple files.***

**Parameters**


``--mvf`` (required) = One or more mvf files. (type=file path, default=None)

``--out`` (required) = Output file (type=file path, default=None)

``--line-buffer/--linebuffer`` = Number of entries to store in memory at a time. (type=integer, default=100000)

``--main_header_file/--mainheaderfile`` = Output file will use same headers as this input file (default=first in list). (type=None, default=None)


``--new-contigs/--newcontigs`` = By default, contigs are matched between files using their text labels in the header. Use this option to turn matching off and treat each file's contigs as distinct. (flag, default=False)



``--newsamples`` = By default, samples are matched between files using their text labels in the header. Use this option to turn matching off and treat each file's sample columns as distinct. (flag, default=False)



``--overwrite`` = USE WITH CAUTION: force overwrite of outputs (flag, default=False)



``--quiet`` = Suppress screen output. (flag, default=False)


---

## ConvertFasta2MVF
***Converts a FASTA file to MVF format***

**Parameters**


``--fasta`` (required) = input FASTA file(s) (type=None, default=None)

``--out`` (required) = output MVF file (type=None, default=None)


``--contig-by-file/--contigbyfile`` = Contigs are designated by separate files. (flag, default=False)


``--contig-field/--contigfield`` = When headers are split by --field-sep, the 0-based index of the contig id. (type=integer, default=None)

``--contig-ids/--contigids`` = manually specify one or more contig ids as ID:LABEL (type=None, default=None)

``--field-sep/--fieldsep`` = FASTA field separator; assumes '>database accession locus' format (type=None, default=None)
Choices: ['TAB', 'SPACE', 'DBLSPACE', 'COMMA', 'MIXED', 'PIPE', 'AT', 'UNDER', 'DBLUNDER']

``--flavor`` = type of file [dna] or protein (type=None, default=dna)
Choices: ['dna', 'protein']

``--manual-coord/--manualcoord`` = manually specify reference coordinates for each file in the format CONTIGID:START..STOP, ... (type=None, default=None)


``--overwrite`` = USE WITH CAUTION: force overwrite of outputs (flag, default=False)



``--quiet`` = Suppress screen output. (flag, default=False)


``--read-buffer/--readbuffer`` = number of lines to hold in READ buffer (type=integer, default=100000)

``--ref-label/--reflabel`` = label for reference sample (type=None, default=REF)

``--sample-field/--samplefield`` = when headers are split by --field-sep, the 0-based index of the sample id (type=integer, default=None)

``--sample-replace/--samplereplace`` = one or more TAG:NEWLABEL or TAG, items, if TAG found in sample label, replace with NEW (or TAG if NEW not specified) NEW and TAG must each be unique (type=None, default=None)

``--write-buffer/--writebuffer`` = number of lines to hold in WRITE buffer (type=integer, default=100000)

---

## ConvertMAF2MVF
***Converts a MAF file to a MVF file***

**Parameters**


``--maf`` (required) = input MAF file (type=file path, default=None)

``--out`` (required) = output MVF file (type=file path, default=None)

``--ref-tag/--reftag`` (required) = Specify which TAG in --sample-tags is the reference genome. (type=None, default=None)

``--sample-tags/--sampletags`` (required) = One or more TAG:NEW or TAG, items separated by commas.Each TAG is partial text-matched to the sample labels in the MAF. For example, hsap18.chr1 and hsap18.chr2 would be matched tag 'hsap18'. If :NEW is added, then the MVF sample will be labeled NEW.  Otherwise, the sample will be labeled simply TAG. (type=None, default=None)

``--line-buffer/--linebuffer`` = Number of entries to store in memory at a time. (type=integer, default=100000)

``--mvf-ref-label/--mvfreflabel`` = new label for reference sample (default='REF') (type=None, default=REF)


``--overwrite`` = USE WITH CAUTION: force overwrite of outputs (flag, default=False)



``--quiet`` = Suppress screen output. (flag, default=False)


---

## ConvertMVF2Fasta
***Converts an MVF file to a FASTA file***

**Parameters**


``--mvf`` (required) = Input MVF file. (type=file path, default=None)

``--out`` (required) = Output path of FASTA file. (type=file path, default=None)

``--buffer`` = size (Mbp) of write buffer for each sample (type=integer, default=10)


``--gene-mode`` = None (flag, default=False)


``--label-type/--labeltype`` = Long labels with all metadata or short ids (type=None, default=long)
Choices: ('long', 'short')

``--output-data/--outputdata`` = Output dna, rna or prot data. (type=None, default=None)
Choices: ('dna', 'rna', 'prot')


``--quiet`` = Suppress screen output. (flag, default=False)


``--regions`` = Path of a plain text file containing one more lines with entries 'contigid,stop,start' (one per line, inclusive coordinates) all data will be returned if left blank. (type=file path, default=None)

``--sample-indices/--sampleindices`` = Specify comma-separated list of sample numerical indices (first sample is 0). Leave blank for all samples. Do not use with --sample_labels. (type=None, default=None)

``--sample-labels`` = Specify comma-separated list of sample labels. Labels must be exact (case-sensitive). Leave blank for all samples.Do not use with --sample_indicies. (type=None, default=None)

``--temp_dir/--tempdir`` = directory to write temporary fasta files (type=None, default=.)

---

## ConvertMVF2FastaGene
***Converts an MVF file to a set ofFASTA files per gene***

**Parameters**


``--mvf`` (required) = Input MVF file. (type=file path, default=None)

``--output-dir`` (required) = Output directory of FASTA files. (type=file path, default=None)

``--buffer`` = size (Mbp) of write buffer for each sample (type=integer, default=10)

``--choose-allele/--chooseallele`` = Chooses how heterozygous alleles are handled. (none=no splitting (default); random1=pick one allele randomly (type=None, default=none)
Choices: ['none', 'random1']


``--ignore-strand`` = Do not read strand info from contigs (flag, default=False)


``--output-data/--outputdata`` = Output dna, rna or prot data. (type=None, default=None)
Choices: ('dna', 'rna', 'prot')


``--quiet`` = Suppress screen output. (flag, default=False)


``--sample-indices/--sampleindices`` = Specify comma-separated list of sample numerical indices (first sample is 0). Leave blank for all samples. Do not use with --sample_labels. (type=None, default=None)

``--sample-labels`` = Specify comma-separated list of sample labels. Labels must be exact (case-sensitive). Leave blank for all samples.Do not use with --sample_indicies. (type=None, default=None)

``--temp_dir/--tempdir`` = directory to write temporary fasta files (type=None, default=.)

---

## ConvertMVF2Phylip
***Converts an MVF file to a Phylip file***

**Parameters**


``--mvf`` (required) = Input MVF file. (type=file path, default=None)

``--out`` (required) = Output Phylip file. (type=file path, default=None)

``--buffer`` = size (bp) of write buffer for each sample (type=integer, default=100000)

``--label-type/--labeltype`` = Long labels with all metadata or short ids (type=None, default=short)
Choices: ('long', 'short')

``--output-data/--outputdata`` = Output dna, rna or prot data. (type=None, default=None)
Choices: ('dna', 'rna', 'prot')


``--partition`` = Output a CSV partitions file with RAxMLformatting for use in partitioned phylogenetic methods. (flag, default=False)



``--quiet`` = Suppress screen output. (flag, default=False)


``--regions`` = Path of a plain text file containing one more lines with entries 'contigid,stop,start' (one per line, inclusive coordinates) all data will be returned if left blank. (type=file path, default=None)

``--sample-indices/--sampleindices`` = Specify comma-separated list of sample numerical indices (first sample is 0). Leave blank for all samples. Do not use with --sample_labels. (type=None, default=None)

``--sample-labels`` = Specify comma-separated list of sample labels. Labels must be exact (case-sensitive). Leave blank for all samples.Do not use with --sample_indicies. (type=None, default=None)

``--temp_dir/--tempdir`` = directory to write temporary fasta files (type=None, default=.)

---

## ConvertVCF2MVF
***Converts a VCF file to an MVF file***

**Parameters**


``--out`` (required) = output MVF file (type=None, default=None)

``--alleles-from/--allelesfrom`` = get additional alignment columns
                from INFO fields (:-separated) (type=None, default=None)

``--contig-ids/--contigids`` = manually specify one or more contig ids as ID;VCFLABE;MVFLABEL, note that VCFLABEL must match EXACTLY the contig string labels in the VCF file (type=None, default=None)

``--field-sep/--fieldsep`` = VCF field separator (default='TAB') (type=None, default=TAB)
Choices: ['TAB', 'SPACE', 'DBLSPACE', 'COMMA', 'MIXED']


``--filter-nonref-empty`` = Do not output entries that are masked
                        or empty for all samples (not the reference). (flag, default=False)


``--line-buffer/--linebuffer`` = Number of entries to store in memory at a time. (type=integer, default=100000)

``--low-depth/--lowdepth`` = below this read depth coverage, convert to lower case set to 0 to disable (type=integer, default=3)

``--low-qual/--lowqual`` = below this quality convert to lower case set to 0 to disable (type=integer, default=20)

``--mask-depth/--maskdepth`` = below this read depth mask with N/n (type=integer, default=1)

``--mask-qual/--maskqual`` = low quality cutoff, bases replaced by N/- set to 0 to disable (type=integer, default=3)


``--no-autoindex/--noautoindex`` = do not automatically index contigs from the VCF (flag, default=False)


``--out-flavor/--outflavor`` = choose output MVF flavor to include quality scores and/or indels (type=None, default=dna)
Choices: ['dna', 'dnaqual', 'dnaqual-indel', 'dna-indel']


``--overwrite`` = USE WITH CAUTION: force overwrite of outputs (flag, default=False)


``--ploidy`` = Use for hexaploid and tetraploid
                        (Experimental, use with caution (type=integer, default=2)
Choices: (2, 4, 6)


``--qual`` = Include Phred genotype quality (GQ) scores (flag, default=False)



``--quiet`` = Suppress screen output. (flag, default=False)


``--ref-label/--reflabel`` = label for reference sample (default='REF') (type=None, default=REF)

``--sample-replace/--samplereplace`` = one or more TAG:NEWLABEL or TAG, items, if TAG found in sample label, replace with NEW (or TAG if NEW not specified) NEW and TAG must each be unique (type=None, default=None)


``--skip-contig-label-check`` = When there are many contigs
                        skip checking for repeat labels.
                        (use with caution). (flag, default=False)


``--vcf`` = VCF input file (type=file path, default=None)


``--verbose`` = Output excessive data to screen for debugging (flag, default=False)


---

## FilterMVF
***Filter an MVF file using various parameters.***

**Parameters**


``--actions`` = set of actions:args to perform, note these are done in order as listed (type=None, default=None)


``--labels`` = use sample labels instead of indices (flag, default=False)


``--line-buffer/--linebuffer`` = Number of entries to store in memory at a time. (type=integer, default=100000)


``--more-help/--morehelp`` = prints full module list and descriptions (flag, default=False)


``--mvf`` = Input MVF file. (type=file path, default=None)

``--out`` = Output file (type=file path, default=None)


``--overwrite`` = USE WITH CAUTION: force overwrite of outputs (flag, default=False)



``--quiet`` = Suppress screen output. (flag, default=False)



``--retain-empty/--retainempty`` = keep empty entries during filtering (flag, default=False)


``--test`` = manually input a line for testing (type=None, default=None)

``--test-nchar/--textnchar`` = total number of samples for test string (type=integer, default=None)


``--verbose`` = report every line (for debugging) (flag, default=False)


---

## InferGroupSpecificAllele
***Infer Group-specific alleles using PAML.***

**Parameters**


``--mvf`` (required) = Input MVF file. (type=file path, default=None)

``--out`` (required) = Output file (type=file path, default=None)


``--all-sample-trees/--allsampletrees`` = Makes trees from all samples instead of only the most complete sequence from each species (flag, default=False)


``--allele-groups/--allelegroups`` = GROUP1:LABEL,LABEL GROUP2:LABEL,LABEL  (type=None, default=None)

``--branch-lrt/--branchlrt`` = Specify the output file for and turn on the RAxML-PAML format LRT test scan for selection on the target branch in addition to the basic patterns scan (type=file path, default=None)

``--chi-test/--chitest`` = Input two number values for expected Nonsynonymous and Synonymous expected values. (type=None, default=None)

``--codeml-path/--codemlpath`` = Full path for PAML codeml executable. (type=file path, default=codeml)

``--end-contig/--endcontig`` = Numerical id for the ending contig. (type=integer, default=100000000)

``--gff`` = Input gff annotation file. (type=file path, default=None)

``--mincoverage`` = Mininum sample coverage for sites. (type=integer, default=None)

``--num-target-species/--targetspec`` = Specify the minimum number of taxa in the target set that are required to conduct analysis (type=integer, default=1)

``--outgroup`` = Specify sample name with which to root trees. (type=None, default=None)

``--output-align/--outputalign`` = Output alignment to this file path in phylip format. (type=None, default=None)

``--paml-tmp/--pamltmp`` = path for temporary folder for PAML output files (type=file path, default=pamltmp)


``--quiet`` = Suppress screen output. (flag, default=False)


``--raxml-path/--raxmlpath`` = Full path to RAxML program executable. (type=file path, default=raxml)

``--species-groups/--speciesgroups`` = None (type=None, default=None)

``--start-contig/--startcontig`` = Numerical ID for the starting contig. (type=integer, default=0)

``--target`` = Specify the taxa labels that define the target lineage-specific branch to be tested. (type=None, default=None)


``--use-labels/--uselabels`` = Use contig labels instead of IDs in output. (flag, default=False)



``--verbose`` = additional screen output (flag, default=False)



``--windowsize`` = Set integer window size. Use 0 for whole file. Use -1 for whole contigs.  (flag, default=100000)


---

## InferTree
***Infer phylogenies for various windows or contigs in anMVF file.***

**Parameters**


``--mvf`` (required) = Input MVF file. (type=file path, default=None)

``--out`` (required) = Output file (type=file path, default=None)

``--bootstrap`` = turn on rapid bootstrapping for RAxML and perform specified number of replicates (type=integer, default=None)

``--choose-allele/--chooseallele/--hapmode`` = Chooses how heterozygous alleles are handled. (none=no splitting (default); randomone=pick one allele randomly (recommended); randomboth=pick two alleles randomly, but keep both; major=pick the more common allele (type=None, default=none)
Choices: ['none', 'randomone', 'randomboth']


``--collapse-polytomies/--collapsepolytomies`` = Collapses internal branches with length 0to polytomies.  Off by default, so arbitrarytopological resolutions in trees may bemaintained if sequences are highly similar. (flag, default=False)


``--contig-ids/--contigids`` = Specify comma-separated list of contig short ids. Must match exactly. Do not use with --contig-labels. (type=None, default=None)

``--contig-labels/--contiglabels`` = Specify comma-separated list of contig full labels. Must match exactly. Do not use with --contig-ids (type=None, default=None)

``--duplicate-seq/--duplicateseq`` = dontuse=remove duplicate sequences prior to RAxML tree inference, then add them to the tree manually as zero-branch-length sister taxa; keep=keep in for RAxML tree inference (may cause errors for RAxML); remove=remove entirely from alignment (type=None, default=dontuse)
Choices: ['dontuse', 'keep', 'remove']

``--engine`` = Choose a phylogenetic inference 'engine' application. The defaultis 'raxml-ng'. (type=None, default=raxml-ng)
Choices: ('raxml', 'raxml-ng')

``--engine-opts/--engineopts/--raxml-opts/--raxmlopts`` = specify additional RAxML arguments as a double-quotes encased string (type=None, default=)

``--engine-path/--enginepath/--raxml-path/--raxmlpath`` = manually specify the path of the phylogenetic engine. (type=None, default=raxml-ng)

``--min-depth/--mindepth`` = minimum number of alleles per site (type=integer, default=4)

``--min-seq-coverage/--minseqcoverage`` = proportion of total alignment a sequencemust cover to be retianed [0.1] (type=float, default=0.1)

``--min-sites/--minsites`` = minimum number of sites  (type=integer, default=100)

``--model/--model/--raxml-model`` = choose model of sequence evolution. defaults are GTRGAMMA for RAxML, or GTR+G for RAxML-ng. (type=None, default=None)


``--output-contig-labels/--outputcontiglabels`` = Output will use contig labels instead of id numbers. (flag, default=False)



``--output-empty/--outputempty`` = Include entries of windows with no data in output. (flag, default=False)



``--quiet`` = Suppress screen output. (flag, default=False)


``--raxml-outgroups/--raxmloutgroups`` = Comma-separated list of outgroup taxon labels to use in RAxML. (type=None, default=None)

``--root-with/--rootwith`` = Comma-separated list of taxon labels to root trees with after RAxML (type=None, default=None)

``--sample-indices/--sampleindices`` = Specify comma-separated list of sample numerical indices (first sample is 0). Leave blank for all samples. Do not use with --sample_labels. (type=None, default=None)

``--sample-labels`` = Specify comma-separated list of sample labels. Labels must be exact (case-sensitive). Leave blank for all samples.Do not use with --sample_indicies. (type=None, default=None)

``--temp-dir/--tempdir`` = Temporary directory path (type=file path, default=./raxmltemp)

``--temp-prefix/--tempprefix`` = Temporary file prefix (type=None, default=mvftree)


``--windowsize`` = Set integer window size. Use 0 for whole file. Use -1 for whole contigs.  (flag, default=100000)


---

## MergeMVF
***Combines columns from multiple MVF files into a single output MVF(this is a newer module, use with caution!)***

**Parameters**


``--mvf`` (required) = One or more mvf files. (type=file path, default=None)

``--out`` (required) = Output file (type=file path, default=None)

``--line-buffer/--linebuffer`` = Number of entries to store in memory at a time. (type=integer, default=100000)

``--main_header_file/--mainheaderfile`` = Output file will use same headers as this input file (default=first in list). (type=None, default=None)


``--new-contigs/--newcontigs`` = By default, contigs are matched between files using their text labels in the header. Use this option to turn matching off and treat each file's contigs as distinct. (flag, default=False)



``--newsamples`` = By default, samples are matched between files using their text labels in the header. Use this option to turn matching off and treat each file's sample columns as distinct. (flag, default=False)



``--overwrite`` = USE WITH CAUTION: force overwrite of outputs (flag, default=False)



``--quiet`` = Suppress screen output. (flag, default=False)



``--skip-index/--skipindex`` = Skip index because index exists (flag, default=False)


---

## PlotChromoplot
***Plot a Chromoplot from an MVF file for all combinationsof the specified samples.***

**Parameters**


``--mvf`` (required) = Input MVF file. (type=file path, default=None)

``--colors`` = three colors to use for chromoplot (type=None, default=None)
Choices: {'lgrey': (250, 250, 250), 'dgrey': (192, 192, 192), 'black': (0, 0, 0), 'white': (255, 255, 255), 'red': (192, 0, 0), 'orange': (217, 95, 2), 'yellow': (192, 192, 0), 'green': (0, 192, 0), 'blue': (0, 0, 192), 'teal': (27, 158, 119), 'puce': (117, 112, 179), 'purple': (192, 0, 192), 'none': ()}

``--contig-ids/--contigids/--contigs`` = Enter the labels of one or more contigs in the order they will appear in the chromoplot (as comma-separated list)(defaults to all ids in order present in MVF) (type=None, default=None)

``--contig-labels/--contiglabels`` = Enter the ids of one or more contigs in the order they will appear in the chromoplot (as comma-separated list)(defaults to all ids in order present in MVF) (type=None, default=None)

``--empty-mask/--emptymask`` = Mask empty regions with this color. (type=None, default=none)
Choices: {'lgrey': (250, 250, 250), 'dgrey': (192, 192, 192), 'black': (0, 0, 0), 'white': (255, 255, 255), 'red': (192, 0, 0), 'orange': (217, 95, 2), 'yellow': (192, 192, 0), 'green': (0, 192, 0), 'blue': (0, 0, 192), 'teal': (27, 158, 119), 'puce': (117, 112, 179), 'purple': (192, 0, 192), 'none': ()}


``--info-track/--infotrack`` = Include an additional coverage information track that will show empty, uninformative, and informative loci. (Useful for ranscriptomes/RAD or other reduced sampling. (flag, default=False)



``--majority`` = Plot only 100% shading in the majority track  rather than shaded proportions in all tracks. (flag, default=False)


``--out-prefix/--outprefix`` = Output prefix (not required). (type=None, default=None)

``--outgroup-indices/--outgroupindices`` = Specify comma-separated list of 1 or more outgroup sample numerical indices (first column is 0). Leave blank for all samples. Do not use with --outgroup_labels. (type=None, default=None)

``--outgroup-labels/--outgrouplabels`` = Specify comma-separated list of 1 or more outgroup sample labels. Labels must be exact (case-sensitive). Leave blank for all samples.Do not use with --outgroup_indicies. (type=None, default=None)

``--plot-type/--plottype`` = PNG image (default) or graph via matplotlib (experimental) (type=None, default=image)
Choices: ['graph', 'image']


``--quiet`` = Suppress screen output. (flag, default=False)


``--sample-indices/--sampleindices`` = Specify comma-separated list of 3 or more sample numerical indices (first sample is 0). Leave blank for all samples. Do not use with --sample_labels. (type=None, default=None)

``--sample-labels`` = Specify comma-separated list of 3 or more sample labels. Labels must be exact (case-sensitive). Leave blank for all samples.Do not use with --sample_indicies. (type=None, default=None)


``--windowsize`` = Set integer window size. Use 0 for whole file. Use -1 for whole contigs.  (flag, default=100000)


``--xscale`` = Width (in number of pixels) for each window (type=integer, default=1)

``--yscale`` = Height (in number of pixels) for each track (type=integer, default=20)

---

## TranslateMVF
***Annotates a chromosomal MVF file with new contigboundaries based on genes/features from a GFF file.***

**Parameters**


``--mvf`` (required) = Input MVF file. (type=file path, default=None)

``--out`` (required) = Output file (type=file path, default=None)

``--filter-annotation/--filterannotation`` = skip GFF entries with text matching this in their 'Notes' field (type=None, default=None)

``--gene-pattern/--genepattern`` = Gene name pattern finder when interpreting GFF/GTF.  Use %% in place of gene name. (type=None, default=gene_id "%")

``--gff`` = Input gff annotation file. (type=file path, default=None)

``--line-buffer/--linebuffer`` = Number of entries to store in memory at a time. (type=integer, default=100000)

``--non-genic-margin/--nongenicmargin`` = For –-non-genic-mode, pad the boundaries of unannotated regions by this amount. (type=integer, default=0)


``--non-genic-mode/--nongenicmode`` = Instead of returning annotated genes, return the non-genic regions without changing contigs or coordinates. (flag, default=False)


``--output-data/--outputdata`` = dna=single data column of dna alleles; protein=single data column of protein alleles; codon=four columns with: protein frame1 frame2 frame3 (type=None, default=codon)
Choices: ['dna', 'protein', 'codon']


``--overwrite`` = USE WITH CAUTION: force overwrite of outputs (flag, default=False)



``--quiet`` = Suppress screen output. (flag, default=False)


``--require-annotation/--requireannotation`` = require GFF entries with text matching this in their 'Notes' field (type=None, default=None)


``--retain-contigs`` = maintain original contig numbering (flag, default=False)



``--retain-coords`` = maintain original coordinates (flag, default=False)


---

## VerifyMVF
***Checks an MVF file for errors.***

**Parameters**


``--mvf`` (required) = Input MVF file. (type=file path, default=None)


``--quiet`` = Suppress screen output. (flag, default=False)


---

## LegacyAnnotateMVF
***This is deprecated now, but maintained for legacy functions.Use TranslateMVF with --output-data dna to annotate regions.Annotates a chromosomal MVF file with new contigboundaries based on genes/features from a GFF file.***

**Parameters**


``--mvf`` (required) = Input MVF file. (type=file path, default=None)

``--out`` (required) = Output file (type=file path, default=None)

``--filter-annotation/--filterannotation`` = Skip entries in the GFF file that contain this string in their 'Notes' (type=None, default=None)

``--gene-pattern/--genepattern`` = Gene name pattern finder when interpreting GFF/GTF.  Use % in place of gene name. (type=None, default=gene_id "%")

``--gff`` = Input gff annotation file. (type=file path, default=None)

``--line-buffer/--linebuffer`` = Number of entries to store in memory at a time. (type=integer, default=100000)

``--nongenic-margin/--nongenicmargin`` = for --nongenic-mode, only retain positions that are this number of bp away from an annotated region boundary (type=integer, default=0)


``--nongenic-mode/--nongenicmode`` = Instead of returning annotated genes, return the non-genic regions without without changing contigs or coordinates (flag, default=False)



``--overwrite`` = USE WITH CAUTION: force overwrite of outputs (flag, default=False)



``--quiet`` = Suppress screen output. (flag, default=False)


---

## LegacyTranslateMVF
***Note this is deprecated now, but maintained for legacy function.Use TranslateMVF with '--output-data protein'or '--output-data codon.'Translate a DNA MVF to a protein or codon MVF***

**Parameters**


``--mvf`` (required) = Input MVF file. (type=file path, default=None)

``--out`` (required) = Output file (type=file path, default=None)

``--filter-annotation/--filterannotation`` = skip GFF entries with text matching this in their 'Notes' field (type=None, default=None)

``--gff`` = Input GFF3 file. If GFF3 not provided, alignments are assumed to be in-frame coding sequences. (type=file path, default=None)

``--line-buffer/--linebuffer`` = Number of entries to store in memory at a time. (type=integer, default=100000)

``--output-data/--outputdata`` = protein=single data column of protein alleles; codon=four columns with: protein frame1 frame2 frame3 (type=None, default=codon)
Choices: ['protein', 'codon']


``--overwrite`` = USE WITH CAUTION: force overwrite of outputs (flag, default=False)


``--parent-gene-pattern/--parentgenepattern`` = Parent genes prefix when interpretingGFF files.  For GFF3 files, 'gene:' is standard, but for older or custom GFF files this may vary.  Use 'none' to make empty. (type=None, default=gene_id "%")


``--quiet`` = Suppress screen output. (flag, default=False)


``--require-annotation/--requireannotation`` = require GFF entries with text matching this in their 'Notes' field (type=None, default=None)


``--verbose`` = Output excessive data to screen for debugging (flag, default=False)


