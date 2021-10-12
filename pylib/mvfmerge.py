# -*- coding: utf-8 -*-
"""
This program checks an MVF file for inconsistencies or errors

MVFtools: Multisample Variant Format Toolkit
James B. Pease and Ben K. Rosenzweig
http://www.github.org/jbpease/mvftools

If you use this software please cite:
Pease, James B. and Benjamin K. Rosenzweig. 2018.
"Encoding Data Using Biological Principles: the Multisample Variant Format
for Phylogenomics and Population Genomics"
IEEE/ACM Transactions on Computational Biology and Bioinformatics. 15(4) 1231–1238.
http://www.dx.doi.org/10.1109/tcbb.2015.2509997

This file is part of MVFtools.

MVFtools is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
MVFtools is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MVFtools.  If not, see <http://www.gnu.org/licenses/>.
"""

from pylib.mvfbase import MultiVariantFile


class MvfTransformer():
    """MVF Transformer Object
       Creates a Data Structure to translate sample and contig names
       across mulitple MVF files
       Arguments:
           labels: dict of sample labels
           contigs: dict of contig info

    """

    def __init__(self, labels=None, contigs=None):
        self.contigs = contigs if contigs else {}
        self.labels = labels if labels else {}
        self.labels_rev = (dict((v, k) for k, v in labels.items()) if
                           labels else {})

    def set_label(self, localindex, consensusindex):
        """Set label transform
        """
        self.labels[consensusindex] = localindex
        self.labels_rev[localindex] = consensusindex
        return ''

    def set_contig(self, localid, consensusid):
        """Set contigid transform
        """
        self.contigs[consensusid] = localid
        return ''


def verify_mvf(args):
    """Main method"""
    args.qprint("Running VerifyMVF")
    mvf = MultiVariantFile(args.mvf, 'read')
    contigs = mvf.metadata['contigs']
    ncol = mvf.metadata['ncol']
    previous_location = (None, None)
    if mvf.metadata['mvftype'] in ('dna', 'protein'):
        if mvf.metadata['mvftype'] == 'protein':
            valid_bases = 'ACDEFGHIKLMNPQRSTVWY'
            valid_characters = 'ACDEFGHIKLMNPQRSTVWYX-'
        else:
            valid_bases = 'ATGCKMRSWY'
            valid_characters = 'ATGCKMRSWYX-'
        for contigid, pos, allelesets in mvf:
            alleles = allelesets[0]
            nonref = False
            if alleles[0] == '@':
                alleles = alleles[1:]
                nonref = True
            errmsg = []
            #  CHECK ALLELES
            if len(alleles) == 1:
                if alleles in 'X-':
                    errmsg.append("no data")
                elif alleles not in valid_bases:
                    errmsg.append("invalid alleles")
            elif len(alleles) == 2:
                if alleles[0] == alleles[1]:
                    errmsg.append("invalid format")
                elif alleles[0] in '-' and not nonref:
                    errmsg.append("empty reference")
                elif (alleles[0] not in valid_characters or
                      alleles[1] not in valid_characters):
                    errmsg.append("invalid alleles")
            elif alleles[1] == '+':
                if alleles[2] == '-':
                    errmsg.append("invalid format")
                elif alleles[0] == '-' and not nonref:
                    errmsg.append("empty reference")
                elif (alleles[0] not in valid_characters or
                      alleles[2] not in valid_characters):
                    errmsg.append("invalid alleles")
                elif int(alleles[3:]) > ncol:
                    errmsg.append("invalid sample number")
            elif alleles[2] == '+':
                if alleles[0] == alleles[1] and alleles[0] == alleles[3]:
                    errmsg.append("invalid format")
                elif any(alleles[x] not in valid_characters
                         for x in (0, 1, 3)):
                    errmsg.append("invalid alleles")
                elif int(alleles[4:]) > ncol:
                    errmsg.append("invalid sample number")
            else:
                if alleles[0] in '-' and not nonref:
                    errmsg.append("empty reference")
                if alleles[0] in '-':
                    errmsg.append("empty reference")
                if any(x not in valid_characters for x in alleles):
                    errmsg.append("invalid alleles")
            #  CHECK POSITION
            if contigid not in contigs:
                errmsg.append("invalid contigid")
            elif pos > contigs[contigid]['length']:
                errmsg.append("invalid position on contig")
            elif contigid != previous_location[0]:
                previous_location = (contigid, pos)
            elif pos <= previous_location[1]:
                errmsg.append("position out of order")
            #  PRINT MESSAGES
            if errmsg:
                print(contigid, pos, allelesets, errmsg)
    elif mvf.metadata['mvftype'] == 'codon':
        args.qprint("codon checking coming soon")
    return ''


def concatenate_mvf(args):
    """Main method"""
    args.qprint("Running ConcatenateMVF")
    concatmvf = MultiVariantFile(args.out, 'write', overwrite=args.overwrite)
    # Copy the first file's metadata
    if args.main_header_file:
        if args.main_header_file not in args.mvf:
            raise RuntimeError("{} not found in files".format(
                args.main_header_file))
        args.main_header_file = args.mvf.index(args.main_header_file)
    else:
        args.main_header_file = 0
    first_mvf = MultiVariantFile(args.mvf[args.main_header_file], 'read')
    concatmvf.metadata = first_mvf.metadata.copy()
    # Open each MVF file, read headers to make unified header
    transformers = []
    for mvfname in args.mvf:
        # This will create a dictionary of samples{old:new}, contigs{old:new}
        transformer = MvfTransformer()
        mvf = MultiVariantFile(mvfname, 'read')
        mvf.reset_max_contig()
        for i, label in enumerate(mvf.get_sample_ids()):
            if label not in concatmvf.get_sample_ids():
                concatmvf.metadata['labels'].append(label)
                concatmvf.metadata['samples'][
                    concatmvf.metadata['labels'].index(label)] = {
                        'label': label}
            if concatmvf.metadata['labels'].index(label) != i:
                transformer.set_label(
                    i, concatmvf.metadata['labels'].index(label))
        for contigid, contigdata in iter(mvf.metadata['contigs'].items()):
            if contigdata['label'] not in [
                    concatmvf.metadata['contigs'][x]['label']
                    for x in concatmvf.metadata['contigs']]:
                newid = (contigid if
                         contigid not in concatmvf.metadata['contigs'] else
                         concatmvf.get_next_contig_id())
                concatmvf.metadata['contigs'][newid] = contigdata
            else:
                for concatid, concatdata in (
                        concatmvf.metadata['contigs'].items()):
                    if contigdata['label'] == concatdata['label']:
                        newid = concatid
                        break
            if newid != contigid:
                transformer.set_contig(contigid, newid)
        transformers.append(transformer)
    # Write output header
    concatmvf.notes.append(args.command_string) 
    concatmvf.write_data(concatmvf.get_header())
    # Now loop through each file
    entries = []
    nentries = 0
    for ifile, mvfname in enumerate(args.mvf):
        args.qprint("Processing {} ...\n".format(mvfname))
        transformer = transformers[ifile]
        mvf = MultiVariantFile(mvfname, 'read')
        for contigid, pos, allelesets in mvf.iterentries(decode=False):
            if transformer.labels:
                allelesets = [mvf.decode(x) for x in allelesets]
                for j, alleles in enumerate(allelesets):
                    allelesets[j] = concatmvf.encode(''.join([
                        x in transformer.labels and
                        alleles[transformer.labels[x]] or alleles[x]
                        for x in range(len(alleles))]))
            if transformer.contigs:
                contigid = (transformer['contigs'][contigid] if
                            contigid in transformer['contigs'] else
                            contigid)
            entries.append((contigid, pos, allelesets))
            nentries += 1
            if nentries == args.line_buffer:
                concatmvf.write_entries(entries)
                entries = []
                nentries = 0
        if entries:
            concatmvf.write_entries(entries)
            entries = []
            nentries = 0
    return ''


def merge_mvf(args):
    """Main method"""
    args.qprint("Running MergeMVF")
    if any(fpath.endswith('.gz') for fpath in args.mvf):
        print("WARNING! Running MergeMVF with gzipped input files is "
              "extremely slow and strongly discouraged.")
    concatmvf = MultiVariantFile(args.out, 'write', overwrite=args.overwrite)
    # Copy the first file's metadata
    args.qprint("Reading First File and Establishing Output")
    if args.main_header_file:
        if args.main_header_file not in args.mvf:
            raise RuntimeError("{} not found in files".format(
                args.main_header_file))
        args.main_header_file = args.mvf.index(args.main_header_file)
    else:
        args.main_header_file = 0
    first_mvf = MultiVariantFile(args.mvf[args.main_header_file], 'read')
    concatmvf.copy_header(first_mvf)
    # Open each MVF file, read headers to make unified header
    transformers = []
    mvfmetadata = []
    inputfiles = []
    for mvfname in args.mvf:
        args.qprint("Reading headers from {}".format(mvfname))
        # This will create a dictionary of samples{old:new}, contigs{old:new}
        args.qprint("Processing Headers and Indexing: {}".format(mvfname))
        transformer = MvfTransformer()
        mvf = MultiVariantFile(mvfname, 'read',
                               contigindex=(not args.skip_index))
        if args.skip_index:
            mvf.read_index_file()
        mvf.reset_max_contig()
        mvfmetadata.append(mvf.metadata)
        for i, sid in enumerate(mvf.get_sample_ids()):
            if sid not in concatmvf.get_sample_ids():
                new_sindex = concatmvf.max_sample_index + 0
                concatmvf.max_sample_index += 1
                concatmvf.sample_indices.append(new_sindex)
                concatmvf.sample_ids.append(sid)
                concatmvf.sample_data[new_sindex] = {}
                concatmvf.sample_data[new_sindex]['id'] = sid
                concatmvf.sample_id_to_index[sid] = new_sindex
            transformer.set_label(
                i, concatmvf.sample_id_to_index[sid])
        for cindex in mvf.contig_indices:
            if (mvf.contig_data[cindex]['label'] not in
                    concatmvf.contig_label_to_index):
                new_cindex = (mvf.contig_data[cindex]['id'] if
                              mvf.contig_data[cindex]['id'] not in
                              concatmvf.contig_ids else
                              concatmvf.get_next_contig_index())
                concatmvf.contig_data[new_cindex] = (
                    mvf.contig_data[cindex].copy())
            else:
                new_cindex = concatmvf.contig_label_to_index[
                    mvf.contig_data[cindex]['label']]
            transformer.set_contig(cindex, new_cindex)
        transformers.append(transformer)
        inputfiles.append(mvf)
    # Write output header
    args.qprint("Writing headers to merge output")
    concatmvf.reset_max_sample()
    concatmvf.notes.append(args.command_string) 
    concatmvf.write_data(concatmvf.get_header())
    # Now loop through each file
    blank_entry = '-' * len(concatmvf.sample_indices)
    for cons_contig in concatmvf.contig_indices:
        contig_merged_entries = {}
        args.qprint("Merging Contig Index: {}".format(cons_contig))
        for ifile, mvffile in enumerate(inputfiles):
            if cons_contig not in transformers[ifile].contigs:
                continue
            localcontig = transformers[ifile].contigs[cons_contig]
            if 'idx' not in mvffile.contig_data[localcontig]:
                print("not found")
                continue
            for _, pos, allelesets in mvffile.itercontigentries(
                    localcontig, decode=True):
                if pos not in contig_merged_entries:
                    contig_merged_entries[pos] = blank_entry[:]
                for j, base in enumerate(allelesets[0]):
                    xcoord = transformers[ifile].labels_rev[j]
                    if contig_merged_entries[pos][xcoord] != '-':
                        if contig_merged_entries[pos][xcoord] == base:
                            continue
                        if base in '-X':
                            continue
                        raise RuntimeError(
                            ("Merging columns have two different bases: "
                             "{} {} {}").format(
                                 pos, contig_merged_entries[pos][xcoord], base))
                    contig_merged_entries[pos] = (
                        contig_merged_entries[pos][:xcoord] + base +
                        contig_merged_entries[pos][xcoord+1:])
        if contig_merged_entries:
            concatmvf.write_entries(((cons_contig, coord, (entry,))
                                     for coord, entry in
                                     sorted(contig_merged_entries.items())
                                     ), encoded=False)
        args.qprint("Entries written for contig {}: {}".format(
            cons_contig, len(contig_merged_entries)))
    return ''
