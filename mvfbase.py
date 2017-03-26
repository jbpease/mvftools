# -*- coding: utf-8 -*-
"""
MVFtools: Multisample Variant Format Toolkit
http://www.github.org/jbpease/mvftools

If you use this software please cite:
Pease JB and BK Rosenzweig. 2016.
"Encoding Data Using Biological Principles: the Multisample Variant Format
for Phylogenomics and Population Genomics"
IEEE/ACM Transactions on Computational Biology and Bioinformatics. In press.
http://www.dx.doi.org/10.1109/tcbb.2015.2509997

MVFbase: Base class MVF handler and functions
@author: James B. Pease
@author: Ben K. Rosenzweig

version: 2015-02-01 - First Public Release
version: 2015-02-26 - Efficiency upgrades for iterators
version: 2015-06-09 - MVF1.2.1 upgrade
version: 2015-09-04 - Small style fixes
version: 2015-12-15 - Python3 compatibilty fix
version: 2015-12-31 - Header and cleanup
version:  2016-01-01 - Python3 compatiblity fix
version:  2016-01-11 - fix for dna ambiguity characters
version: 2016-08-02 - Python3 conversion, integrate analysis_base
@version: 2016-09-10 - Minor fixes to gz reading

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
import sys
import gzip
from itertools import groupby

# Math Functions


def is_int(num):
    """Checks if num is an integer, accepts strings as well"""
    try:
        return not bool(float(num) % 1)
    except ValueError:
        return False
    except TypeError:
        return False


def interpret_param(string):
    """Check if value is a bool, then int, then float, or returns string"""
    if string.lower() in ['true', 't', 'yes']:
        return True
    elif string.lower() in ['false', 'f', 'no']:
        return False
    else:
        try:
            if not bool(float(string) % 1):
                return int(string)
            raise ValueError
        except ValueError:
            try:
                return float(string)
            except ValueError:
                return string


def abpattern(num, digits=0):
    return bin(num)[2:].zfill(digits).replace('0', 'A').replace('1', 'B')


# FASTA FILE HANDLER
def fasta_iter(fasta_name):
    """
        given a fasta file. yield tuples of header, sequence
        Adapted from https://github.com/brentp
    """
    filehandler = open(fasta_name, 'rt')
    faiter = (x[1] for x in groupby(filehandler, lambda line: line[0] == ">"))
    for header in faiter:
        header = next(header)[1:].strip()
        seq = "".join(s.strip() for s in next(faiter))
        yield header, seq


# MVF Class Object

class MultiVariantFile(object):
    """Multisample Variant Format Handler
    Object Structure:
        path = file path (converted to absolute path)
        filemode = read/write mode
        entrystart: auto-generated file coordinate where entries start
        metadata = dict of associated data including (but not limited to):
            contigs: dict[id] = dict(metadata)
            isgzip: boolean if file is gzip-compressed
            flavor: [dna, codon, protein, dnaqual, dnaindel]
            samples: dict[index] = dict(sample_info)
            ncol: auto-generated int(number of samples)
            labels: auto-generated tuple of sample labels
    """

    def __init__(self, path, filemode, **kwargs):
        self.path = os.path.abspath(path)
        self.metadata = {'flavor': 'dna'}
        self.metadata['labels'] = []
        self.metadata['samples'] = {}
        self.metadata['contigs'] = {}
        if filemode not in ('read', 'r', 'rb', 'write', 'w', 'wb'):
            raise RuntimeError("Invalid filemode {}".format(filemode))
        self.filemode = filemode
        self.entrystart = 0
        # Check for Gzip and establish file object
        self.metadata['isgzip'] = (self.path.endswith(".gz") or
                                   kwargs.get('isgzip', False))
        # READ MODE
        if filemode in ('read', 'r', 'rb'):
            if os.path.exists(self.path):
                filehandler = (self.metadata.get('isgzip', False) and
                               gzip.open(self.path, 'rt') or
                               open(self.path, 'rt'))
                # Process header lines
                header_lines = []
                self.entrystart = filehandler.tell()
                line = filehandler.readline()
                while line:
                    if not line.startswith("#"):
                        break
                    header_lines.append(line.rstrip())
                    self.entrystart = filehandler.tell()
                    line = filehandler.readline()
                self._process_header(header_lines)
                # Establish number of columns
                self.metadata['ncol'] = len(self.metadata['labels'])
            else:
                raise IOError("MVF path {} not found!".format(self.path))
        # WRITE MODE
        elif filemode in ('write', 'w', 'wb'):
            if os.path.exists(self.path) and not kwargs.get('overwrite',
                                                            False):
                raise IOError(
                    """MVF path {} already exists, use --overwrite
                    to replace""".format(self.path))
            filehandler = (self.metadata.get('isgzip', False) and
                           gzip.open(self.path, 'wt') or
                           open(self.path, 'wt'))
            self.metadata['ncol'] = kwargs.get('ncol', 2)
        filehandler.close()
        self.metadata['flavor'] = (kwargs.get('flavor', False) or
                                   self.metadata.get('flavor', 'dna'))

    def _process_header(self, headerlines):
        """Processes header lines into metadata
            Arguments:
                headerlines: list of unprocessed header strings
        """
        sample_index = 0
        for line in headerlines:
            try:
                entry = line.split()
                if entry[0].startswith('##mvf'):
                    for elem in entry[1:]:
                        if '=' in elem:
                            elem = elem.split('=')
                            self.metadata[elem[0]] = interpret_param(elem[1])
                        else:
                            self.metadata[elem] = True
                elif entry[0].startswith('#s'):
                    self.metadata['samples'][sample_index] = {
                        'label': entry[1]}
                    self.metadata['labels'].append(entry[1])
                    for elem in entry[2:]:
                        if '=' in elem:
                            elem = elem.split('=')
                            self.metadata['samples'][sample_index][elem[0]] = (
                                interpret_param(elem[1]))
                        else:
                            self.metadata['samples'][sample_index][elem] = True
                    sample_index += 1
                elif entry[0].startswith('#c'):
                    contigid = entry[1]
                    self.metadata['contigs'][contigid] = {}
                    for elem in entry[2:]:
                        if '=' in elem:
                            elem = elem.split('=')
                            self.metadata['contigs'][contigid][elem[0]] = (
                                interpret_param(elem[1]))
                        else:
                            self.metadata['contigs'][contigid][elem] = True
                else:
                    raise ValueError
            except ValueError:
                sys.stderr.write("Skipping invalid header line '{}'\n".format(
                    entry))
        return ''

    def get_sample_indices(self, labels=None):
        """Get indices for a specified named group or set of labels
           Arguments:
                labels: list/set of str or str; default= all
        """
        if labels is None:
            return range(len(self.metadata['labels']))
        try:
            if hasattr(labels, '__iter__'):
                return [self.metadata['labels'].index(x) for x in labels]
            else:
                return self.metadata['labels'].index(labels)
        except IndexError:
            raise IndexError(labels, "contains invalid label")

    def get_sample_labels(self, indices=None):
        """Get labels for the specified named indices
            Arguments:
                indices: list/set of int or int; default= all
        """
        if indices is None:
            return self.metadata['labels']
        try:
            if hasattr(indices, '__iter__'):
                return [self.metadata['labels'][x] for x in indices]
            else:
                return self.metadata['labels'][indices]
        except IndexError:
            raise IndexError(indices, "contains invalid label index")

    def get_contig_id(self, label):
        """Returns contig id given contig label
        """
        for contigid, contig in self.metadata['contigs'].items():
            if contig['label'] == label:
                return contigid
        raise IndexError("contig '{}' not found".format(label))

    def get_contig_ids(self):
        """Returns all contig ids
        """
        return self.metadata['contigs'].keys()

    def get_contig_label(self, contigid):
        """Returns contig label given contig id
        """
        if contigid in self.metadata['contigs']:
            return self.metadata['contigs'][contigid]['label']
        raise IndexError("contig '{}' not found".format(contigid))

    def get_contig_labels(self):
        """Returns contig labels as a list
        """
        return [self.metadata['contigs'][x]['label']
                for x in self.metadata['contigs']]

    def get_next_contig_id(self):
        """Returns the (highest integer id) + 1 or 0"""
        maxid = 0
        for contigid in self.metadata['contigs']:
            if is_int(contigid):
                maxid = max(maxid, int(contigid))
        return str(maxid + 1)

    def __iter__(self, quiet=False):
        """Simple entry iterator
           Returns (str(chrom), int(pos), list(allele entries))
        """
        linecount = 0
        if self.metadata['isgzip']:
            filehandler = gzip.open(self.path, 'rt')
        else:
            filehandler = open(self.path, 'rt')
        filehandler.seek(self.entrystart)
        for line in filehandler:
            if line.strip() == '':
                    continue
            try:
                arr = line.rstrip().split()
                loc = arr[0].split(':')
                yield (loc[0], int(loc[1]), arr[1:])
            except:
                raise RuntimeError(
                    "Error processing MVF at line# {} = {} ".format(
                        linecount, line))
        filehandler.close()

    def iterentries(self, decode=True, contigs=None, no_invariant=False,
                    no_gap=False, no_ambig=False, no_nonref=False,
                    onlyalleles=False, quiet=False, subset=None):
        """
        Fully-optioned iterator for MVF entries with filtering
        Returns (str(chrom), int(pos), list(allele entries))
                 or list(allele entries) with 'onlyalleles'

        Arguments:
            contigs: list of contig ids to include (default=all)
            decode:fully decode the allele sets (T/F)
            no_invariant: set to false to skip invariant sites
            no_ambig: set to false to skip positions with 'N'
            no_gap: set to false to skip positions with '-'
            no_nonref: set to false to skip nonref contigs
            onlyalleles: return only list of alleles
            quiet: suppress progress meter (default=False)
            subset: list of column indices

        Note: for codons, filters must apply to all allele strings
        Note: using subset without decode returns encoded subset
        """
        if not contigs:
            if no_nonref:
                contigs = sorted([x for x in self.metadata['contigs']
                                  if self.metadata['contigs'][x].get(
                                      'ref', False)])
            else:
                contigs = sorted(self.metadata['contigs'].keys())
        else:
            contigs = [x if x in self.metadata['contigs'] else
                       self.get_contig_id(x)
                       for x in contigs]
        subset = subset or ''
        current_contigid = ''
        linecount = 0
        if self.metadata['isgzip']:
            filehandler = gzip.open(self.path, 'rb')
        else:
            filehandler = open(self.path, 'r')
        filehandler.seek(self.entrystart)
        for line in filehandler:
            try:
                if self.metadata['isgzip']:
                    arr = line.decode().rstrip().split()
                else:
                    arr = line.rstrip().split()
                loc = str(arr[0]).split(':')
                contigid = loc[0]
                pos = int(loc[1])
                allelesets = arr[1:]
                if contigid != current_contigid:
                    if current_contigid in contigs:
                        contigs.remove(current_contigid)
                        if not contigs:
                            break
                    current_contigid = contigid[:]
                if contigid not in contigs:
                    continue
                if subset:
                    try:
                        allelesets = [''.join([alleles[j] for j in subset])
                                      for alleles in [self.decode(x)
                                                      for x in allelesets]]
                    except:
                        raise RuntimeError(allelesets)
                if no_gap:
                    if any(x in allelesets[0] for x in '-@'):
                        continue
                if no_ambig:
                    if any('X' in x for x in allelesets):
                        continue
                if no_invariant:
                    if all(len(x) == 1 for x in allelesets):
                        continue
                    elif len(allelesets[0]) > 1 and (
                            allelesets[0][1] == '+' and
                            allelesets[0][2] == allelesets[0][0]):
                        continue
                    elif subset:
                        if all(x == allelesets[0][0]
                               for x in allelesets[0][1:]):
                            continue

                if subset and not decode:
                    allelesets = [self.encode(x) for x in allelesets]
                if decode and not subset:
                    allelesets = [self.decode(x) for x in allelesets]
                if onlyalleles:
                    yield allelesets
                else:
                    yield (contigid, pos, allelesets)
            except:
                raise RuntimeError(
                    "Error processing MVF at line# {} = {} ".format(
                        linecount, line))
        filehandler.close()

    def get_header(self):
        """Returns formatted header string (with final newline)
        """
        header = ["##mvf version=1.2 flavor={} {}".format(
            self.metadata['flavor'], ' '.join([
                "{}={}".format(k, v) for (k, v) in sorted(
                    self.metadata.items())
                if k not in ('contigs', 'samples', 'flavor', 'version',
                             'isgzip', 'labels')]))]
        header.extend(["#s {} {}".format(
            self.metadata['samples'][x]['label'],
            ' '.join(["{}={}".format(k, v) for (k, v) in (
                sorted(self.metadata['samples'][x].items())) if k != 'label']))
            for x in range(len(self.metadata['samples']))])
        contigs = [is_int(k) and (int(k), v) or (k, v)
                   for (k, v) in self.metadata['contigs'].items()]
        header.extend(["#c {} label={} length={} {}".format(
            cid, cdata['label'], cdata['length'],
            ' '.join(["{}={}".format(k, v) for k, v in (
                sorted(cdata.items())) if k not in ['length', 'label']]))
            for cid, cdata in (sorted(contigs))])
        return '\n'.join(header) + '\n'

    def decode(self, alleles):
        """Decode entry into full-length alleles
            Arguments:
                alleles = encoded allele string
        """
        ref = True
        ncol = self.metadata['ncol']
        if alleles.startswith('@'):
            alleles = alleles[1:]
            ref = False
            ncol -= 1
        if len(alleles) == 1:
            alleles = alleles[0] * ncol
        elif len(alleles) == 2:
            alleles = alleles[0] + (alleles[1] * (ncol - 1))
        elif alleles[1] == '+':
            newalleles = [alleles[0]] + ['-'] * (ncol - 1)
            newalleles[int(alleles[3:])] = alleles[2]
            alleles = ''.join(newalleles)
        elif alleles[2] == '+':
            tmp = int(alleles[4:])
            alleles = "{}{}{}{}".format(alleles[0],
                                        alleles[1]*(tmp - 1),
                                        alleles[3],
                                        alleles[1]*(ncol - tmp - 1))
        if not ref:
            return '@' + alleles
        return alleles

    def encode(self, alleles):
        """Internal copy of encode_mvfstring
            Arguments:
                alleles: unencoded allele string
        """
        if self.metadata['flavor'] == 'dna':
            return encode_mvfstring(alleles).replace(
                'N', 'X').replace('n', 'X')
        else:
            return encode_mvfstring(alleles)

    def write_data(self, data):
        """Writes datastring to the MVF file
            Argunments:
                data: string datastream
        """
        if self.metadata['isgzip'] or self.path.endswith('.gz'):
            with gzip.open(self.path, 'at') as outfile:
                outfile.write(data)
        else:
            with open(self.path, 'at') as outfile:
                outfile.write(data)

    def write_entries(self, entries, encoded=True):
        """Write MVF entries
            Arguments:
                entries: list of entry tuples (contigid, pos, alleles)
                encoded: entries have been pre-encoded (default=True)
        """
        self.write_data('\n'.join(["{}:{} {}".format(
            entry[0], entry[1], (
                ' '.join([encoded and x or encode_mvfstring(x)
                          for x in entry[2]])))
                for entry in entries]) + '\n')
        return ''


def encode_mvfstring(alleles):
    """Encode full-length alleles to MVF short form
    """
    denovo = False
    if alleles.startswith('@'):
        alleles = alleles[1:]
        denovo = True
    if all(x == alleles[1] for x in alleles[2:]):
        if alleles[0] == alleles[1]:
            alleles = alleles[0]
        else:
            alleles = alleles[0:2]
    elif alleles[1:].count(alleles[1]) == len(alleles) - 2:
        pos = [(j, x) for j, x in enumerate(alleles[1:]) if x != alleles[1]]
        alleles = "{}{}+{}{}".format(
            alleles[0], alleles[1] != '-' and alleles[1] or '',
            pos[0][1], pos[0][0] + 1)
    elif all(x == alleles[2] for x in alleles[3:]):
        alleles = "{}{}+{}1".format(
            alleles[0], alleles[2] != '-' and alleles[2] or '', alleles[1])
    if denovo:
        return '@' + alleles
    return alleles

# ANALYSIS BACKEND


class Counter(dict):
    """dict subclass for counter with integer values"""

    def add(self, key, val=1):
        """Adds integer value to key, default = 1)"""
        self[key] = self.get(key, 0) + val

    def subtract(self, key, val=1):
        """Subtracts integer value to key, default = 1)"""
        self[key] = self.get(key, 0) - val
        return ''

    def iter_sorted(self, percentage=True, n_values=None, sort='desc'):
        """Iterates values sorting in a list"""
        total = float(sum(self.values()))
        if not total:
            yield "Nothing counted!"
        entries = sorted([(v, k) for (k, v) in self.items()],
                         reverse=(sort != 'asc'))
        for (val, k) in entries[:n_values]:
            if percentage:
                yield "{} {} {}%\n".format(
                    k, val, round(float(val) * 100/total, 2))
            else:
                yield "{} {}\n".format(k, val)


class OutputFile(object):
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
        with open(self.path, 'wt') as outfile:
            outfile.write('#' + '\t'.join(self.headers) + "\n")
        return ''

    def write_entry(self, entry):
        """Writes entry to file
            Arguments:
                entry: dict of values with keys matching header
        """
        with open(self.path, 'at') as outfile:
            outfile.write("\t".join([str(k not in entry and '.' or entry[k])
                                     for k in self.headers]) + "\n")
        return ''


class AnalysisModule(object):
    """General Functions for Analysis Modules
        Params:
            data = data dict
            params = parameters dict
            **kwargs = all other keywords are added to params as key=value

    """
    def __init__(self, data=None, params=None, **kwargs):
        self.data = data or {}
        self.params = params or {}
        self.params.update([(k, v) for k, v in kwargs.items()
                            if k not in ['data', 'params']])


if __name__ == ("__main__"):
    print("""MVF base handler library v. 2016-08-02, please run one of the
          other MVFtools scripts to access these functions""")
