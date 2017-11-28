******************
mvf_filter modules
******************

GENERAL NOTES
=============
mvf_filter is a script that processes an MVF file using a variety
of modules that can be used in any combination of orders.  There
are three types of actions: 

* Transformations: alter the character strings and may remove empty entries
* Filters: remove entries that meet specific criteria
* Location: remove entries based on their genomic location

Modules can be used in any order and as many as you like.  However,
this means that when multiple transformations are used any changes
to the column numbering must be accounted for.  For example,
if you want to remove columns 3 and then 5, you have to specify this as 
"columns:0,1,2,4,5 columns:0,1,2,3" since after the first 
transformation column 5 would become the new column 4.


allelegroup 
=============
This filter requires that all members of each group contain
valid alleles.  The groups are specified by a series of colon-separated
groups of comma-separate columns.
EXAMPLE ACTION: allelegroup:1,2,3:4,5,6
EXAMPLE #1 AA-AATA --> *retained* (first and second group both have alleles)
EXAMPLE #2 A-X-ATA --> *filtered out* (first group does not have valid alleles)
EXAMPLE #3 AACC--- --> *filtered out* (second group does not have valid alleles)
a

collapsepriority
==================
This transformation will combine the alleles from several 
columns using a priority ranked order. This is useful for collapsing 
low-coverage samples into a single combined sample column.
The columns  are specified after the colon using comma-separated integers 
(or text labels with the --labels option).
EXAMPLE ACTION: collapsepriority:2,3,4 
EXAMPLE #1 ABCDE --> ABC   (column 3 present, so column 3 used)
EXAMPLE #2 AB-DE --> ABD   (column 3 is a gap, so column 4 used)
EXAMPLE #3 ABX-E --> ABE   (column 3 is ambig, 4 is gap, so column 5 used.


collapsemerge
==================
This transformation combines alleles from several columns 
into a single representative allele. This is useful for 
combining haplotypes or population samples. The columns 
are specified after the colon using comma-separated integers 
(or text labels with the --labels option).
EXAMPLE ACTION: collapsemerge:2,3,4
EXAMPLE #1 AACAA --> AAM (CAA becomes ambiguity code 'M')
EXAMPLE #2 AACAG --> AAX (CAG would be 'V'. However, X is used since triallelic is not allowed in MVF.
EXAMPLE #3 AAT-T --> AAT (both non-gap columns are 'T' so T is just used.

columns
=========
This transformation returns only the specified columns.
The columns are specified after the colon using comma-separated integers 
(or text labels with the --labels option).
EXAMPLE ACTION: columns:1,3
EXAMPLE #1 ABCDE --> BD (columns 1 and 3 are returned)
EXAMPLE #2 A-C-E --> [filtered out] (Since there is no data in columns 1 and 3.

maskchar
=========
This transformation will replace the specified character(s) with "X".
Characters to be masked are specified after the column 
as a comma-separated list of single characters.
EXAMPLE ACTION: maskchar:K,M
EXAMPLE #1: AAKA --> AAXA
EXAMPLE #2: AAMX --> AAXX                                                                                               


masklower
===========
This transformation will replace all lower case characters with "X".
This takes no paramters.
EXAMPLE ACTION: masklower
EXAMPLE #1: AaTa --> AXTX
EXAMPLE #2: aaaa --> XXXX
 
mincoverage
=============
This filter will remove entries with fewer non-gap/ambiguous alleles 
than the specified cutoff. This is useful before conducting scans
(such as phylogenetic scans or chromoplots ) that require a minimum 
number of taxa.  The action is specified by a single integer after 
the colon. 
EXAMPLE ACTION: mincoverage:3
EXAMPLE #1: A--A --> *filtered out* (coverage = 2)
EXAMPLE #2: AA-A --> *retained* (coverage = 3)

"notchar
=========
This filter will remove entries with any of the specifed characters.
This can be useful for removing entries with ambiguous characters 
or missing data.  Note that these are *case sensitive* so lower-case 
characters should be entered alongside upper-case when both are 
filtered.  The action is specified by one or more comma-separated 
characters after the colon.
EXAMPLE ACTION: notchar:X,K,M
EXAMPLE #1: AK-X --> *filtered out* (contains K and X)
EXAMPLE #2: AA-A --> *retained* (contains none of specific characters)

promotelower
==============
This transformation will change all lower-case characters to upper-case.
This takes no paramters.
EXAMPLE ACTION: promotelower
EXAMPLE #1: AaTa --> AATA
EXAMPLE #2: aaaa --> AAAA

removelower
=============
This transformation will change all lower-case characters to gaps.
This action takes no paramters.

::

  EXAMPLE ACTION: promotelower
  EXAMPLE #1: AaTa --> A-T-
  EXAMPLE #2: aaaa --> ----

removechar
============
This transformation will change all instances of the specified
characters to gaps. Characters are *case sensitive*. The action is 
specified by one or more comma-separated characters after the colon.

::
  EXAMPLE ACTION: removechar:a
  EXAMPLE #1: AaTa --> A-T-
  EXAMPLE #2: aaaa --> ----


reqallchar
============
This filter will remove entries that do no contain all of the specified 
characters. Characters are *case sensitive*. The action is 
specified by one or more comma-separated characters after the colon.

::
  EXAMPLE ACTION: reqallchar:A,K
  EXAMPLE #1: AaTa --> *filtered out* (contains "A" but not "K")
  EXAMPLE #2: aKaa --> *filtered out* (contains "K" and "a" but not "A")
  EXAMPLE #3: AKAT --> *retained*


reqcontig
=========
This location filter removes entries not on the specified contig.
The action is specified by a numerical contig id after the colon.

::
 EXAMPLE ACTION: reqcontig:1
 EXAMPLE #1: 1:100 AAA --> *retained* 
 EXAMPLE #2: 2:110 AAA --> *filtered out*
 EXAMPLE #3: X:101 AAA --> *filtered out*
 

reqinformative
==============
This filter removes sites without at least two instances of
at least two alleles (phylogenetically informative sites). 
This action takes no paramters.

::
 EXAMPLE ACTION: reqinformative
 EXAMPLE #1: AATA --> *filtered out* (only one "T")
 EXAMPLE #2: ATTA --> *retained* (contains "A" and "T" twice)
 EXAMPLE #3: ATCA --> *filtered out* (only one each of "T" and "C")


reqinvariant
============
This filter removes variant sites (not including gaps or ambiguities)
This action takes no paramters.

 ::
  EXAMPLE ACTION: reqinvariant
  EXAMPLE #1: AATA --> *filteredout* 
  EXAMPLE #2: AAAA --> *retained*
  EXAMPLE #3: AA-A --> *retained
  EXAMPLE #3: AAXA --> *retained

reqregion
=========
This location filter removes entries not on the specified contig
within in the specified bounds.
The action is specified by a numerical contig id, then start and 
stop coordinates (inclusive) after the colon.
 
 ::
  EXAMPLE ACTION: reqregion:1,101,110
  EXAMPLE #1: 1:100 AAA --> *filtered out*
  EXAMPLE #2: 1:110 AAA --> *retained* 
  EXAMPLE #3: 2:101 AAA --> *filtered out*


reqonechar
==========
This filter will remove entries that do no contain at least
one of the of the specified  characters. Characters are 
*case sensitive*. The action is specified by one or more 
comma-separated characters after the colon.
 ::
  EXAMPLE ACTION: reqonechar:A,K
  EXAMPLE #1: AaTa --> *retained* 
  EXAMPLE #2: CTCC --> *filtered out* 
  EXAMPLE #3: aaTC --> *filtered out*



reqsample
=========
This filter requires that the given sample(s) be a non-gap/ambiguous
allele. The action is specified by one or more
comma-separated integer column indices after the colon.
EXAMPLE ACTION: reqample:1,2
EXAMPLE #1: AAAA --> *retained* 
EXAMPLE #2: A-AA --> *filtered out* 
EXAMPLE #3: AA-A --> *filtered out*

reqvariant
==========
This filter removes invariant sites.
This action takes no paramters.
 ::
  EXAMPLE ACTION: reqinvariant
  EXAMPLE #1: AATA --> *retained* 
  EXAMPLE #2: AAAA --> *filtered out*
  EXAMPLE #3: AA-A --> *filtered out*
  EXAMPLE #4: AAXA --> *filtered out*


reqnonrefsample
===============
This filter removes sites with no non-reference information.
This action takes no paramters.
 ::
  EXAMPLE ACTION: reqnonrefsample
  EXAMPLE #1: AATA --> *retained* 
  EXAMPLE #2: A--A --> *retained*
  EXAMPLE #3: A--- --> *filtered out*

