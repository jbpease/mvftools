# -*- coding: utf-8 -*-
"""

MVFtools: Multisample Variant Format Toolkit
http://www.github.org/jbpease/mvftools

mvf_window_tree: RAxML trees from sequence regions encoded by MVF format
@author: James B. Pease
@author: Ben K. Rosenzweig

Version: 2015-02-01 - First Public Release

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

from __future__ import print_function
import sys, argparse, os, subprocess
from mvfbase import  MultiVariantFile
from random import randint
from itertools import combinations
from mvfbiolib import HAPSPLIT
from Bio import Phylo
from datetime import datetime
from StringIO import StringIO

class WindowData(object):
    """Container for window data
        Arguments:
            labels: sequence labels
            hapsplt: split mode (none, randomone, randomboth, majorminor)
            seqs: list of sequence strings
    """
    def __init__(self, labels=None, hapmode="none", seqs=None):
        if hapmode in ['randomboth', 'majorminor']:
            self.labels = [label + x for x in ['a', 'b'] for label in labels]
        else:
            self.labels = labels[:]
        if not seqs:
            self.seqs = [[] for _ in xrange(len(self.labels))]
        self.hapmode = hapmode

    def append_alleles(self, alleles):
        """Add alleles to the window
            Arguments:
                alleles: list of allele strings
        """
        for j, allele in enumerate(alleles):
            try:
                self.seqs[j].append(allele)
            except:
                raise RuntimeError(alleles, len(self.seqs), len(alleles))

    def maketree_raxml(self, **kwargs):
        """Makes alignment, runs RAXML, and processes output
            Arguments:
                kwargs: passthrough arguments
        """
        jobname = 'mvftree-' + datetime.now().strftime(
            '%Y-%m-%d-%H-%M-%S-') + str(randint(100000, 999999))
        if not os.path.exists("raxmltemp"):
            os.mkdir("raxmltemp")
        os.chdir('raxmltemp')
        temp_labels = {}
        rev_labels = {}
        if len(self.seqs[0]) < kwargs.get('minsites', 0):
            os.chdir('..')
            return 'nodata'
        duplicates = ''
        if not kwargs.get('keep_duplicate_sequences', False):
            duplicates = self.remove_duplicates()
        for j in xrange(len(self.labels)):
            if len(self.labels[j]) > 10:
                temp_labels['s' + str(j)] = self.labels[j]
                rev_labels[self.labels[j]] = 's' + str(j)
                self.labels[j] = 's' + str(j)
        kwargs['outgroups'] = kwargs.get('outgroups', False) and [
            rev_labels.get(x, x) for x in kwargs['outgroups']] or []

        temp_file = os.path.abspath(os.curdir + "/" + jobname + "_temp.phy")
        self.write_reduced_phylip(temp_file)
        try:
            run_raxml(temp_file, jobname, **kwargs)
            #REGULAR TREE
            if kwargs.get('bootstrap', 0):
                tree = Phylo.read('RAxML_bipartitions.' + jobname, 'newick',
                                  comments_are_confidence=True)
            else:
                tree = Phylo.read('RAxML_bestTree.' + jobname, 'newick')
        except IOError:
            os.chdir('..')
            return 'nodata'
        for node in tree.get_terminals():
            kwargs['outgroup'] = temp_labels.get(node.name, node.name)
            node.name = temp_labels.get(node.name, node.name)
            if node.name in duplicates:
                name_list = [node.name] + duplicates[node.name]
                node.split(
                    n=len(duplicates[node.name]) + 1, branch_length=0.00000)
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
        if kwargs['outgroups']:
            rootwith = kwargs['outgroups'][0]
        else:
            rootwith = self.labels[0]
        topology = ladderize_alpha_tree(tree.__format__('newick'),
                                        rootwith=rootwith,
                                        prune=kwargs.get('prune', []))
        os.chdir('..')
        return {'tree': full_tree.strip(),
                'topology': topology.strip(),
                'templabels': ';'.join(["{}={}".format(k, v)
                                        for (k, v) in temp_labels.items()]),
                'alignlength': len(self.seqs[0]),
                'aligndepth': len(self.seqs)}

    def write_reduced_phylip(self, outpath):
        """Write Phylip"""
        with open(outpath, 'w') as alignfile:
            alignfile.write("{} {}\n".format(len(self.seqs),
                                             len(self.seqs[0])))
            for j in xrange(len(self.seqs)):
                alignfile.write("{} {}\n".format(self.labels[j],
                                                 ''.join(self.seqs[j])))
        return ''

    def remove_duplicates(self):
        """Remove duplicate sequences"""
        duplicates = {}
        remove_indices = set([])
        for i, j in combinations(range(len(self.seqs)), 2):
            if self.seqs[i] == self.seqs[j]:
                remove_indices.update([j])
                if self.labels[i] in duplicates:
                    duplicates[self.labels[i]].append(self.labels[j])
                else:
                    duplicates[self.labels[i]] = [self.labels[j]]
        for j in sorted(remove_indices, reverse=True):
            del self.seqs[j]
            del self.labels[j]
        return duplicates

def run_raxml(filename, jobname, **kwargs):
    """Runs RAxML
        Parameters:
            filename: temporary phy filepath
            jobname: unique jobid for RAxML
            raxmlpath: path to RAxML program (default=raxmlHPC)
            outgroup: list of outgroup labels
            threads: multithreading tasks (experimental)
            bootstrap: int bootstrap replicates (default=none)
            optbranch: run -k option for RAxML
    """
    raxml_exec = 'raxmlpath' in kwargs and kwargs['raxmlpath'] or 'raxmlHPC'
    model = 'model' in kwargs and kwargs['model'] or 'GTRGAMMA'
    outgroups = kwargs.get('outgroups', '') or []
    logpath = 'logpath' in kwargs and kwargs['logpath'] or (jobname + '.log')
    log_file = open(logpath, 'a')
    arr_cmd = [raxml_exec, '-m', model, '-n', jobname,
               '-s', filename, '-p', str(randint(1000000, 999999999))]
    if outgroups:
        arr_cmd.extend(['-o', ','.join(outgroups)])
    if kwargs.get('threads', 1) > 1:
        arr_cmd.extend(['-T', str(kwargs['threads'])])
    if kwargs.get('bootstrap', 0):
        if kwargs.get('optbranch', 0):
            arr_cmd.extend(['-k', '-#',
                            str(kwargs['bootstrap']), '-f', 'a', '-x',
                            str(randint(1, 100000000))])
        else:
            arr_cmd.extend(['-#', str(kwargs['bootstrap']), '-f', 'a', '-x',
                            str(randint(1, 100000000))])
    process = subprocess.Popen(' '.join(arr_cmd),
                               shell=True, stdout=log_file, stderr=log_file)
    process.communicate()
    log_file.close()
    return ''


def ladderize_alpha_tree(treestring, prune=None, rootwith=None):
    """Ladderizes and Alphabetizes the Tree to create a unique
    string for each topology"""
    tree0 = Phylo.read(StringIO(treestring), 'newick')
    for node in tree0.get_terminals():
        if prune:
            if any(x in node.name for x in prune):
                tree0.prune(node)
            elif rootwith in node.name:
                tree0.root_with_outgroup(node)
    for node in tree0.find_clades():
        node.branch_length = 1.0
    for node in tree0.get_nonterminals(order="postorder"):
        if not sum([int(not x.is_terminal) for x in node.clades]):
            node.clades.sort(key=lambda c: sorted([
                x.name for x in c.get_terminals()])[0])
        else:
            node.clades.sort(key=lambda c: c.name)
    tree0.ladderize()
    tree_string = tree0.__format__('newick').replace(
        ':1.00000', '').rstrip()
    return tree_string

def hapsplit(alleles, mode):
    """Process Alleles into Haplotypes"""
    if all(x not in 'RYMKWS' for x in alleles):
        if mode in ['major', 'minor', 'randomone']:
            return alleles
        elif mode in ['majorminor', 'randomboth']:
            return ''.join([base*2 for base in alleles])
    else:
        if mode in ['major', 'minor', 'majorminor']:
            hapleles = ''.join([HAPSPLIT[x] for x in alleles])
            counts = sorted([(hapleles.count(x), x) for x in set(hapleles)],
                            reverse=True)
            order = [x[1] for x in counts]
            newalleles = []
            for base in alleles:
                if base in 'RYMKWS':
                    newalleles.extend([x for x in order if x in HAPSPLIT[base]])
                else:
                    newalleles.extend([base, base])
            if mode == 'major':
                alleles = ''.join([x[0] for base in newalleles])
            elif mode == 'minor':
                alleles = ''.join([x[1] for x in newalleles])
            elif mode == 'majorminor':
                alleles = ''.join([x for x in newalleles])
        elif mode == 'randomone':
            alleles = ''.join([HAPSPLIT[x][randint(0, 1)] for x in alleles])
        elif mode == 'randomboth':
            randx = randint(0, 1)
            alleles = ''.join([HAPSPLIT[x][randx] + HAPSPLIT[x][1 - randx]
                               for x in alleles])
        return alleles

def main(arguments=sys.argv[1:]):
    """Main MVF Treemaker"""
    parser = argparse.ArgumentParser(description="""
    Process MVF into alignment""")
    parser.add_argument("--mvf", help="inputmvf")
    parser.add_argument("--out", help="alignment file")
    parser.add_argument("--regions", nargs='*',
                        help="one or more regions id,start,stop (inclusive)")
    parser.add_argument("--samples", nargs='*',
                        help="one or more taxon labels, leave blank for all")
    parser.add_argument("--outgroups", nargs="*")
    parser.add_argument("--contigs", nargs='*')
    parser.add_argument("--hapmode", default="none",
                        choices=["none", "randomone", "randomboth",
                                 "major", "minor", "majorminor"])
    parser.add_argument("--windowsize", type=int, default=10000)
    parser.add_argument("--minsites", type=int, default=100)
    parser.add_argument("--raxmlpath",
                        help="manually specify RAxML path")
    parser.add_argument("--keep_duplicate_sequences", action="store_true")
    parser.add_argument("--quiet", action="store_true",
                        help="suppress screen output")
    parser.add_argument("-v", "--version", action="store_true",
                        help="display version information")
    args = parser.parse_args(args=arguments)
    if args.version:
        print("Version 2015-02-01: Initial Public Release")
        sys.exit()
    args.contigs = args.contigs or []
    args.outgroups = args.outgroups or []
    regions = []
    if args.regions:
        try:
            regions = [tuple([x.split(',')[0], int(x.split(',')[1]),
                              int(x.split(',')[2])])
                       for x in args.regions]
            for j, (contig, rmin, rmax) in enumerate(regions):
                if rmin > rmax:
                    args.complement = True
                    regions[j] = (contig, rmax, rmin)
            regions.sort()

        except:
            raise RuntimeError(args.regions, "invalid regions")
    mvf = MultiVariantFile(args.mvf, 'read')
    sample_cols = args.samples and mvf.get_sample_indices(args.samples) or []
    labels = mvf.get_sample_labels(sample_cols)
    ## WINDOW START
    headers = ['contig', 'winstart', 'winlen', 'tree', 'topology', 'topoid',
               'templabels', 'alignlength', 'aligndepth']
    with open(args.out, 'w') as outfile:
        outfile.write("\t".join(headers) + "\n")
        current_contig = ''
        window_start = 0
        window = None
        topo_ids = {}
        topo_counts = {}
        for contig, pos, allelesets in mvf.iterentries(
                contigs=args.contigs, subset=sample_cols, quiet=args.quiet,
                no_invariant=True, no_ambig=True, no_gap=True, decode=True):
            if regions:
                if not any(rcontig == contig and rstart <= pos <= rstop
                           for (rcontig, rstart, rstop) in regions):
                    continue
            if not current_contig:
                current_contig = contig[:]
                window = WindowData(labels=labels, hapmode=args.hapmode)
                window_start = 0
            elif (pos > window_start + args.windowsize
                  or contig != current_contig):
                entry = window.maketree_raxml(
                    outgroups=args.outgroups, minsites=args.minsites,
                    raxmlpath=args.raxmlpath,
                    keep_duplicate_sequences=args.keep_duplicate_sequences)
                if entry != 'nodata':
                    entry.update([("contig", current_contig),
                                  ("winstart", window_start),
                                  ("winlen", args.windowsize)])
                    topo = entry["topology"]
                    topo_counts[topo] = topo_counts.get(topo, 0) + 1
                    if topo not in topo_ids:
                        topo_ids[topo] = (
                            topo_ids and max(topo_ids.values()) + 1 or 0)
                    entry["topoid"] = topo_ids[topo]
                    outfile.write("\t".join([str(entry[k])
                                             for k in headers]) + "\n")
                if contig != current_contig:
                    window_start = 0
                    current_contig = contig[:]
                else:
                    window_start += args.windowsize
                window = WindowData(labels=labels, hapmode=args.hapmode)
            if args.hapmode != 'none':
                allelesets[0] = hapsplit(allelesets[0], args.hapmode)
            elif any(x in allelesets[0] for x in 'RYWSKM'):
                continue
            window.append_alleles(allelesets[0])
        entry = window.maketree_raxml(
            outgroups=args.outgroups, minsites=args.minsites,
            raxmlpath=args.raxmlpath,
            keep_duplicate_sequences=args.keep_duplicate_sequences)
        if entry != 'nodata':
            entry.update([("contig", current_contig),
                          ("winstart", window_start),
                          ("winlen", args.windowsize)])
            topo = entry["topology"]
            topo_counts[topo] = topo_counts.get(topo, 0) + 1
            if topo not in topo_ids:
                topo_ids[topo] = topo_ids and max(topo_ids.values()) + 1 or 0
            entry["topoid"] = topo_ids[topo]
            outfile.write("\t".join([str(entry[k])
                                     for k in headers]) + "\n")
    window = None
    ## END WINDOW ITERATION
    topo_list = sorted([(v, k) for k, v in topo_counts.iteritems()],
                       reverse=True)
    with open(args.out + ".counts", 'w') as countfile:
        for rank, [value, topo] in enumerate(topo_list):
            countfile.write("#{}\t{}\t{}\n".format(rank, value, topo))
            topo_ids[topo] = str(rank)
    return ''


if __name__ == "__main__":
    main()
