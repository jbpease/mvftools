import sys
from mvfbase import MultiVariantFile

mvf = MultiVariantFile(sys.argv[1], 'read')

#for entry in mvf.iterentries(contig_ids=['57881']):
#    print(entry)

print(mvf.get_header())
