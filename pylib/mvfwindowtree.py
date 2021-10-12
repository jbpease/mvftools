#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This program makes phylogenies from individual genomic windows of
a DNA MVF alignment (Requires: BioPython).

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

import os
import subprocess
from random import randint
from datetime import datetime
from io import StringIO
from itertools import combinations
from Bio import Phylo
from pylib.mvfbase import MultiVariantFile, same_window
from pylib.mvfbiolib import MvfBioLib
MLIB = MvfBioLib()


class WindowData():
    """Container for window data
        Arguments:
            labels: sequence labels
            hapsplt: split mode (none, randomone, randomboth, majorminor)
            seqs: list of sequence strings
    """
    def __init__(self, window_params=None, seqs=None):
        self.contigname = window_params.get('contigname', 'CONTIG')
        self.windowstart = window_params.get('windowstart', -1)
        self.windowsize = window_params.get('windowsize', -1)
        self.labels = window_params.get('labels', '')
        self.seqs = [[] for _ in range(len(self.labels))]

    def append_alleles(self, alleles, mindepth=1):
        """Add alleles to the window
            Arguments:
                alleles: list of allele strings
        """
        site_depth = mindepth + 0
        i = 0
        while site_depth and i < len(alleles):
            if alleles[i] in 'AaTtGgCcUu':
                site_depth -= 1
            i += 1
        if not site_depth:
            for j, allele in enumerate(alleles):
                try:
                    self.seqs[j].append(allele)
                except:
                    raise RuntimeError(alleles, len(self.seqs), len(alleles))
        return ''

    def prepare_alignment(self, params):
        """Checks and prepares alignment for RAxML
            Arguments:
                params: dict job parameters

            Returns list of any duplicates and empty entry when all checks pass
            In the case of errors, returns entry with information
        """
        # CHECK OVERALL ALIGNMENT DEPTH
        if len(self.seqs) < params.get('mindepth', 4):
            return ('', {'status': 'few',
                         'comment': "Not enough taxa {} (< {}).".format(
                             len(self.seqs), params.get('mindepth', 4)),
                         'contig': self.contigname,
                         'windowstart': self.windowstart,
                         'windowsize': self.windowsize,
                         'alignlength': (
                             len(self.seqs[0]) if len(self.seqs) > 0 else 0),
                         'aligndepth': len(self.seqs) or 0})
        # CHECK FOR OVERALL ALIGNMENT LENGTH
        if len(self.seqs[0]) < params.get('minsites', 0):
            return ('', {'status': 'short',
                         'comment': "Not enough sites {} (< {})".format(
                             len(self.seqs[0]), params.get('minsites', 0)),
                         'contig': self.contigname,
                         'windowstart': self.windowstart,
                         'windowsize': self.windowsize,
                         'alignlength': len(self.seqs[0]),
                         'aligndepth': len(self.seqs)})
        # CHECK FOR EMPTY SEQUENCES
        self.remove_empty_sequences()
        if not self.seqs:
            return ('', {'status': 'few.proc',
                         'comment': 'no sequences left after filtering',
                         'contig': self.contigname,
                         'windowstart': self.windowstart,
                         'windowsize': self.windowsize,
                         'alignlength': (self.seqs and
                                         len(self.seqs[0]) or 0),
                         'aligndepth': len(self.seqs) or 0})
        # CHECK FOR EMPTY SITES
        self.remove_empty_sites()
        if not self.seqs[0]:
            return ('', {'status': 'short.proc',
                         'comment': 'no sites left after filtering',
                         'contig': self.contigname,
                         'windowstart': self.windowstart,
                         'windowsize': self.windowsize,
                         'alignlength': (self.seqs and
                                         len(self.seqs[0]) or 0),
                         'aligndepth': len(self.seqs) or 0})
        # CHECK FOR MINIMUM SEQUENCE COVERAGE
        if params.get('minseqcoverage', 0):
            self.remove_low_coverage_sequences(
                mincov=params.get('minseqcoverage', 0))
            if not self.seqs:
                return ('', {'status': 'spotty.proc',
                             'comment': 'no sequences left after filtering',
                             'contig': self.contigname,
                             'windowstart': self.windowstart,
                             'windowsize': self.windowsize,
                             'alignlength': (self.seqs and
                                             len(self.seqs[0]) or 0),
                             'aligndepth': len(self.seqs) or 0})
        # CHECK AGAIN FOR EMPTY SITES
        self.remove_empty_sites()
        if len(self.seqs[0]) < params.get('minsites', 0):
            return ('', {'status': 'short.proc',
                         'comment': 'no sites left after filtering',
                         'contig': self.contigname,
                         'windowstart': self.windowstart,
                         'windowsize': self.windowsize,
                         'alignlength': len(self.seqs[0]),
                         'aligndepth': len(self.seqs)})
        # CHECK FOR IDENTICAL DUPLICATE SEQUENCES BY BASE COMPOSITION
        duplicates = ''
        if params.get('duplicateseq', 'dontuse') in ('dontuse', 'remove'):
            duplicates = self.remove_duplicates(params.get('duplicateseq',
                                                           'dontuse'))
        # CHECK AGAIN FOR OVERALL ALIGNMENT DEPTH
        if len(self.seqs) < params.get('few.dup', 4):
            return ('', {'status': 'mindepth',
                         'comment': ("Not enough sequences after filtering"
                                     "{} (< {}).").format(
                                         len(self.seqs),
                                         params.get('few.dup', 4)),
                         'contig': self.contigname,
                         'windowstart': self.windowstart,
                         'windowsize': self.windowsize,
                         'alignlength': (
                             len(self.seqs[0]) if len(self.seqs) > 0 else 0),
                         'aligndepth': len(self.seqs) or 0})

        # IF ALL OTHER CHECKS PASS, RETURN DUPLICATE LIST AND EMPTY ENTRY
        return (duplicates, {})

    def maketree_raxml(self, params):
        """Runs RAXML, and processes output
            Arguments:
                params: dict job parameters
        """
        # Prepare alignment and check for problems
        (duplicates, entry) = self.prepare_alignment(params)
        if entry:
            return entry
        # Alignment Has Passed Preliminary Checks, ready for RAxML run
        # Establish Temporary Directory
        jobname = "{}.{}-{}".format(
            params.get("tempprefix", "mvftree"),
            datetime.now().strftime('%Y-%m-%d-%H-%M-%S'),
            randint(100000, 999999))
        temp_filepath = os.path.abspath(
            "{}/{}_temp.phy".format(params['tempdir'], jobname))
        # Temporarily Shorten Labels for Use in Phylip
        temp_labels = {'encode': {}, 'decode': {}}
        for labellen in (len(x) > 10 for x in self.labels):
            if labellen:
                for j, oldlabel in enumerate(self.labels):
                    newlabel = 's{}'.format(j)
                    temp_labels['encode'][oldlabel] = newlabel
                    temp_labels['decode'][newlabel] = oldlabel
                    self.labels[j] = newlabel
                break
        params['outgroups'] = (
            [temp_labels['encode'].get(x, x) for x in params['rootwith']] if
            params['rootwith'] is not None else [])
        self.write_phylip(temp_filepath)
        try:
            if params['engine'] == 'raxml':
                run_raxml(temp_filepath, jobname, params)
                if params.get('bootstrap', 0):
                    tree = Phylo.read('RAxML_bipartitions.' + jobname, 'newick',
                                      comments_are_confidence=True)
                else:
                    tree = Phylo.read('RAxML_bestTree.' + jobname, 'newick')
            elif params['engine'] == 'raxml-ng':
                run_raxml_ng(temp_filepath, jobname, params)
                if params.get('bootstrap', 0):
                    tree = Phylo.read(jobname + '.raxml.support', 'newick',
                                      comments_are_confidence=True)
                else:
                    tree = Phylo.read(jobname + '.raxml.bestTree', 'newick')
        except IOError:
            return {'status': 'IOerror',
                    'contig': self.contigname,
                    'windowstart': self.windowstart,
                    'windowsize': self.windowsize,
                    'alignlength': len(self.seqs[0]),
                    'aligndepth': len(self.seqs)}
        if temp_labels['decode']:
            self.labels = [temp_labels.get(x, x) for x in self.labels]
        for node in tree.get_terminals():
            node.name = temp_labels['decode'].get(node.name, node.name)
            if node.name in duplicates:
                name_list = [node.name] + duplicates[node.name]
                node.split(n=len(duplicates[node.name]) + 1,
                           branch_length=0.00000)
                node.name = None
                for j, subnode in enumerate(node.get_terminals()):
                    subnode.name = name_list[j]
        tree.ladderize()
        full_tree = tree.__format__('newick').replace(
            ':1.00000;', ';').replace("'", '')
        for node in tree.find_clades():
            node.branch_length = None
            if not node.is_terminal:
                node.name = ''
        if params['rootwith'] is not None:
            params['rootwith'] = [x for x in params['rootwith']
                                  if x in self.labels]
        topology = ladderize_alpha_tree(
            #tree.__format__('newick'),
            full_tree,
            rootwith=params.get('rootwith', self.labels[0]),
            prune=params.get('prune', []),
            collapse_polytomies=params.get('collapse_polytomies', False)
            )
        return {'status': 'ok',
                'duplicates': bool(duplicates),
                'contig': self.contigname,
                'windowstart': self.windowstart,
                'windowsize': self.windowsize,
                'tree': full_tree.strip(),
                'topology': topology.strip(),
                'alignlength': len(self.seqs[0]),
                'aligndepth': len(self.seqs)}

    def write_phylip(self, outpath):
        """Write Phylip"""
        with open(outpath, 'w') as alignfile:
            alignfile.write("{} {}\n".format(len(self.seqs),
                                             len(self.seqs[0])))
            for j in range(len(self.seqs)):
                alignfile.write("{} {}\n".format(self.labels[j],
                                                 ''.join(self.seqs[j])))
        return ''

    def remove_empty_sequences(self):
        """Remove sequences with no data"""
        removal_list = []
        for j in range(len(self.seqs)):
            if all([x in 'Xx-' for x in self.seqs[j]]):
                removal_list.append(j)
        for j in sorted(removal_list, reverse=True):
            del self.seqs[j]
            del self.labels[j]
        return ''

    def remove_low_coverage_sequences(self, mincov=0.):
        """Removes lowcoverage sequences"""
        removal_list = []
        if mincov:
            for i, seq in enumerate(self.seqs):
                cov = (float(sum([int(x in 'AaTtGgCcUu') for x in seq])) /
                       (float(len(seq))))
                if cov < mincov:
                    removal_list.append(i)
        for j in sorted(removal_list, reverse=True):
            del self.seqs[j]
            del self.labels[j]
        return ''

    def remove_empty_sites(self):
        """Remove sites without data"""
        removal_list = []
        for j in range(len(self.seqs[0])):
            isempty = True
            for base in (seq[j] for seq in self.seqs):
                if base in 'AaTtGgCcUu':
                    isempty = False
                    break
            if isempty:
                removal_list.append(j)
        if removal_list:
            self.seqs = [[base for (j, base) in enumerate(seq)
                          if j not in removal_list]
                         for seq in self.seqs]
        return ''

    def remove_duplicates(self, mode='dontuse'):
        """Remove duplicate sequences"""
        duplicates = {}
        remove_indices = set([])
        for i, j in combinations(range(len(self.seqs)), 2):
            duplicate = True
            for k, base in enumerate(self.seqs[i]):
                if base != self.seqs[j][k]:
                    duplicate = False
                    break
            if duplicate:
                remove_indices.update([j])
                if mode == 'dontuse':
                    if self.labels[i] not in duplicates:
                        duplicates[self.labels[i]] = []

                    duplicates[self.labels[i]].append(self.labels[j])
        for j in sorted(remove_indices, reverse=True):
            del self.seqs[j]
            del self.labels[j]
        return duplicates


class OutputFile():
    """Set up Output File
        Params:
            path: file path
            headers: list of header elements
    """

    def __init__(self, path, headers):
        self.headers = headers
        self.path = os.path.abspath(path)
        self.write_headers()

    def write_headers(self):
        """Write headers to file"""
        with open(self.path, 'w') as outfile:
            outfile.write('#' + '\t'.join(self.headers) + "\n")
        return ''

    def write_entry(self, entry):
        """Writes entry to file
            Arguments:
                entry: dict of values with keys matching header
        """
        with open(self.path, 'a') as outfile:
            outfile.write("\t".join([str(entry.get(k, '.'))
                                     for k in self.headers]) + "\n")
        return ''


def verify_raxml(params):
    """verify raxml path"""
    try:
        out = str(subprocess.check_output([params['engine_path'], "-v"]))
        if out.find("RAxML") == -1:
            raise RuntimeError("Program not found at path:{}\n{}".format(
                params['engine_path'], out))
    except FileNotFoundError:
        raise RuntimeError("Program not found at path: {}".format(
            params['engine_path']),)
    return ''


def run_raxml(filename, jobname, params):
    """Runs RAxML
        Parameters:
            filename: temporary phy filepath
            jobname: unique jobid for RAxML
            engine_path: path to RAxML program (default=raxmlHPC)
            outgroup: list of outgroup labels
            threads: multithreading tasks (experimental)
            bootstrap: int bootstrap replicates (default=none)
            optbranch: run -k option for RAxML
            opts: any additional options (as a string)
    """
    outgroups = params.get('outgroups', '') or []
    logpath = params.get('logpath', jobname + '.log')
    log_file = open(logpath, 'a')
    cmd = [params.get('engine_path', 'raxmlHPC'),
           '-m', params.get('model', 'GTRGAMMA'),
           '-n', jobname,
           '-s', filename,
           '-p', str(randint(1, 100000000))]
    if outgroups:
        cmd.extend(['-o', ','.join(outgroups)])
    if params.get('bootstrap', 0):
        if params.get('optbranch', 0):
            cmd.extend(['-k', '-#', str(params['bootstrap']),
                        '-f', 'a', '-x', str(randint(1, 100000000))])
        else:
            cmd.extend(['-#', str(params['bootstrap']), '-f', 'a', '-x',
                        str(randint(1, 100000000))])
    cmd.append(params.get('engine_opts', ''))
    process = subprocess.Popen(' '.join(cmd), shell=True,
                               stdout=log_file, stderr=log_file)
    process.communicate()
    log_file.close()
    return ''


def run_raxml_ng(filename, jobname, params):
    """Runs RAxML-ng
        Parameters:
            filename: temporary phy filepath
            jobname: unique jobid for RAxML
            engine_path: path to RAxML program (default=raxmlHPC)
            bootstrap: int bootstrap replicates (default=none)
            opts: any additional options (as a string)
    """
    outgroups = params.get('outgroups', '') or []
    logpath = params.get('logpath', jobname + '.log')
    log_file = open(logpath, 'a')
    cmd = [params.get('engine_path', 'raxml-ng'),
           '--model', params.get('model', 'GTR+G'),
           '--prefix', jobname,
           '--msa', filename
           ]
    if outgroups:
        cmd.extend(['--outgroup', ','.join(outgroups)])
    if params.get('bootstrap', 0):
        cmd.extend(['-all', '--bs-trees', str(params['bootstrap'])])
    cmd.append(params.get('engine_opts', ''))
    process = subprocess.Popen(' '.join(cmd), shell=True,
                               stdout=log_file, stderr=log_file)
    process.communicate()
    log_file.close()
    return ''


def ladderize_alpha_tree(treestring, prune=None, rootwith=None,
                         collapse_polytomies=False):
    """Ladderizes and Alphabetizes the Tree to create a unique
       string for each topology
    """
    tree0 = Phylo.read(StringIO(treestring), 'newick')
    for node in tree0.get_terminals():
        if prune:
            if any([x in node.name for x in prune]):
                tree0.prune(node)
    if collapse_polytomies is True:
        tree0.collapse_all(lambda c: c.branch_length == 0.0 
                           and c != tree0.root)
    for node in tree0.find_clades():
        node.branch_length = None
    for node in tree0.get_nonterminals(order="postorder"):
        if not sum([int(not x.is_terminal) for x in node.clades]):
            node.clades.sort(key=lambda c: sorted([
                x.name for x in c.get_terminals()])[0])
        else:
            node.clades.sort(key=lambda c: c.name)
    tree0.ladderize()
    if rootwith:
        if len(rootwith) == 1:
            for node in tree0.get_terminals():
                if node.name == rootwith[0]:
                    tree0.root_with_outgroup(node)
        else:
            for node in tree0.get_nonterminals():
                if all([subnode.name in rootwith for subnode
                        in node.get_terminals()]):
                    tree0.root_with_outgroup(node)
                    break
    tree_string = tree0.__format__('newick').replace(
        ':1.00000', '').replace(":0.00000", "").rstrip()
    return tree_string


def hapsplit(alleles, mode):
    """Process Alleles into Haplotypes"""
    if all([x not in 'RYMKWS' for x in alleles]):
        return (alleles if mode in ['major', 'minor', 'randomone'] else
                ''.join([base*2 for base in alleles]))
    if mode in ['major', 'minor', 'majorminor']:
        hapleles = ''.join([MLIB.hapsplit(x, mode=mode) for x in alleles])
        counts = sorted([(hapleles.count(x), x) for x in set(hapleles)],
                        reverse=True)
        order = [x[1] for x in counts]
        newalleles = []
        for base in alleles:
            if base in 'RYMKWS':
                newalleles.extend([
                    x for x in order if x in MLIB.hapsplit(base, mode=mode)])
            else:
                newalleles.extend([base, base])
        if mode == 'major':
            alleles = ''.join([x[0] for x in newalleles])
        elif mode == 'minor':
            alleles = ''.join([x[1] for x in newalleles])
        elif mode == 'majorminor':
            alleles = ''.join(newalleles)
    elif mode == 'randomone':
        alleles = ''.join([MLIB.hapsplit(x, mode=mode) for x in alleles])
    elif mode == 'randomboth':
        alleles = ''.join([MLIB.hapsplit(x, mode=mode) for x in alleles])
    return alleles


def infer_window_tree(args):
    """Main method"""
    args.qprint("Running InferTree")
    # ESTABLISH FILE OBJECTS
    mvf = MultiVariantFile(args.mvf, 'read')
    args.qprint("Read MVF File: {}".format(args.mvf))
    # Set up contig ids
    if args.contig_ids is not None:
        contig_ids = args.contig_ids[0].split(",")
    elif args.contig_labels is not None:
        contig_ids = mvf.get_contig_ids(
            labels=args.contig_labels[0].split(","))
    else:
        contig_ids = mvf.get_contig_ids()
    treefile = OutputFile(
        args.out,
        headers=['contig', 'windowstart', 'windowsize', 'tree',
                 'topology', 'topoid',
                 # 'templabels', ### USED FOR DEBUGGING ###
                 'alignlength', 'aligndepth', 'status'])
    topofile = OutputFile(args.out + '.counts',
                          headers=['rank', 'topology', 'count'])
    if args.sample_indices is not None:
        sample_indices = [int(x) for x in
                          args.sample_indices[0].split(",")]
    elif args.sample_labels is not None:
        sample_indices = mvf.get_sample_indices(
            ids=args.sample_labels[0].split(","))
    else:
        sample_indices = mvf.get_sample_indices()
    if not os.path.exists(args.temp_dir):
        os.mkdir(args.temp_dir)
    os.chdir(args.temp_dir)
    # SETUP PARAMS
    main_labels = mvf.get_sample_ids(sample_indices)
    if args.choose_allele in ['randomboth', 'majorminor']:
        main_labels = [label + x for x in ['a', 'b'] for label in main_labels]
    params = {
        'bootstrap': args.bootstrap,
        'chooseallele': args.choose_allele,
        'collapse_polytomies': args.collapse_polytomies,
        'duplicateseq': args.duplicate_seq,
        'engine': args.engine,
        'engine_path': args.engine_path,
        'engine_opts': args.engine_opts,
        'mindepth': args.min_depth,
        'minseqcoverage': args.min_seq_coverage,
        'minsites': args.min_sites,
        'model': args.model,
        'outgroups': (args.raxml_outgroups 
                      if args.raxml_outgroups is not None
                      else None),
        'rootwith': (args.root_with.split(',')
                     if args.root_with is not None
                    else []),
        'tempdir': args.temp_dir,
        'tempprefix': args.temp_prefix,
        'windowsize': args.windowsize,
        }
    # DEFAULT MODEL
    if params['model'] is None:
        if params['engine'] == 'raxml':
            params['model'] = 'GTRGAMMA'
        elif params['engine'] == 'raxml-ng':
            params['model'] = "GTR+G"
    # WINDOW START INTERATION
    verify_raxml(params)
    args.qprint("RAxML Found.")
    current_contig = None
    current_position = 0
    window_data = None
    # skip_contig = False
    topo_ids = {}
    topo_counts = {}
    args.qprint("Prcocessing Records")
    windowsizename = "window size={}".format(args.windowsize)
    if windowsizename == "window size=-1":
        windowsizename = "whole contig"
    elif windowsizename == "window size=0":
        windowsizename = "whole genome"
        window_data = WindowData(window_params={
            'contigname': 'all',
            "windowstart": 0,
            "windowsize": 0,
            "labels": main_labels[:]})
    for contig, pos, allelesets in mvf.iterentries(
            contig_ids=contig_ids, subset=sample_indices,
            no_invariant=False, no_ambig=False, no_gap=False, decode=True):
        # if current_contig == contig:
        #     if skip_contig is True:
        #         args.qprint("Skipping contig: {}".format(current_contig))
        #         continue
        if not same_window((current_contig, current_position),
                           (contig, pos), args.windowsize):
            # skip_contig = False
            if window_data is not None:
                args.qprint(("Making tree for {} "
                             "at contig {} position {}").format(
                                 windowsizename,
                                 current_contig,
                                 current_position))
                entry = window_data.maketree_raxml(params)
                if entry['status'] != 'ok':
                    if args.output_empty:
                        treefile.write_entry(entry)
                    # if args.windowsize != -1:
                    #     skip_contig = True
                    args.qprint(
                        "TREE REJECTED with error code: {} ({})".format(
                            entry['status'], entry.get('comment', "None")))
                else:
                    args.qprint("Tree completed.")
                    topo = entry["topology"]
                    topo_counts[topo] = topo_counts.get(topo, 0) + 1
                    if topo not in topo_ids:
                        topo_ids[topo] = (max(topo_ids.values()) + 1
                                          if topo_ids else 0)
                    entry["topoid"] = topo_ids[topo]
                    treefile.write_entry(entry)
                current_position = current_position + args.windowsize if (
                    contig == current_contig and args.windowsize > 0) else 0
            current_contig = contig[:]
            window_data = None
            window_data = WindowData(window_params={
                'contigname': (mvf.get_contig_labels(ids=current_contig) if
                               args.output_contig_labels is not None else
                               current_contig[:]),
                "windowstart": ('-1' if args.windowsize == -1
                                else current_position + 0),
                "windowsize": args.windowsize,
                "labels": main_labels[:]})
        # ADD ALLELES
        if mvf.flavor == 'dna':
            if args.choose_allele != 'none':
                allelesets[0] = hapsplit(allelesets[0], args.choose_allele)
            window_data.append_alleles(allelesets[0], mindepth=args.min_depth)
        elif mvf.flavor == 'codon':
            for i in (1, 2, 3):
                if args.choose_allele != 'none':
                    allelesets[i] = hapsplit(allelesets[i], args.choose_allele)
                window_data.append_alleles(allelesets[i], mindepth=args.min_depth)
    # LAST LOOP
    if window_data:
        entry = window_data.maketree_raxml(params)
        if entry['status'] != 'ok':
            if args.output_empty:
                treefile.write_entry(entry)
        else:
            topo = entry["topology"]
            topo_counts[topo] = topo_counts.get(topo, 0) + 1
            if topo not in topo_ids:
                topo_ids[topo] = (
                    max(topo_ids.values()) + 1 if topo_ids else 0)
            entry["topoid"] = topo_ids[topo]
            treefile.write_entry(entry)
        window_data = None
    # END WINDOW ITERATION
    topo_list = sorted([(v, k) for k, v in topo_counts.items()],
                       reverse=True)
    for rank, [value, topo] in enumerate(topo_list):
        topofile.write_entry({'rank': rank, 'count': value, 'topology': topo})
    return ''
