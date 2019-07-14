#! -*- coding: utf-8 -*-

"""
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
import sys
import gzip
from itertools import groupby

# ==== Math Functions ====

def is_int(num):
    """Checks if num is an integer, accepts strings as well"""
    try:
        return not bool(float(num) % 1)
    except ValueError:
        return False
    except TypeError:
        return False


def zerodiv(numer, denom):
    """Divide but return zero if denominator is zero
    """
    return 0 if denom == 0 else numer/denom


def interpret_param(string):
    """Check if value is a bool, then int, then float, or returns string"""
    if string.lower() in ['true', 't', 'yes']:
        return True
    if string.lower() in ['false', 'f', 'no']:
        return False
    try:
        if not bool(float(string) % 1):
            return int(string)
        raise ValueError
    except ValueError:
        try:
            return float(string)
        except ValueError:
            return string


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


def same_window(coords1, coords2, windowsize):
    """ coords1/coords1 = a tuple or list with (contig, position)
        windowsize = the windowsize, 0=whole file (always True), -1 contigs
    """
    if windowsize == 0:
        return True
    if windowsize > 0:
        if coords1[0] != coords2[0]:
            return False
        if coords2[1] > coords1[1] + windowsize:
            return False
        return True
    return coords1[0] == coords2[0]


def mixed_sorter(elem):
    """sort a mix of strings and integers
    """
    if isinstance(elem[0], int) is True:
        return (str(elem[0]).zfill(20), elem[1])
    return elem


# MVF Class Object

class MultiVariantFile():
    """Multisample Variant Format Handler
    Object Structure:
        path = file path (converted to absolute path)
        filemode = read/write mode
        entrystart: auto-generated file coordinate where entries start
        flavor: [dna, codon, protein, dnaqual, dnaindel]
        metadata = dict of associated data including (but not limited to):
            contigs: dict[id] = dict(metadata)
            isgzip: boolean if file is gzip-compressed
            samples: dict[index] = dict(sample_info)
            ncol: auto-generated int(number of samples)
            labels: auto-generated tuple of sample labels
    """

    def __init__(self, path, filemode, **kwargs):
        self.path = os.path.abspath(path)
        self.flavor = 'dna'
        self.metadata = {'labels': []}
        self.metadata['samples'] = {}
        self.metadata['contigs'] = {}
        self.metadata['trees'] = []
        self.metadata['notes'] = []
        self.metadata['maxcontigid'] = 0
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
                filehandler = (gzip.open(self.path, 'rt') if
                               self.metadata.get('isgzip', False)
                               else open(self.path, 'rt'))
                # Process header lines
                header_lines = []
                self.entrystart = filehandler.tell()
                line = filehandler.readline()
                while line:
                    if not line.startswith("#"):
                        previous_contig = line.strip().split(":")[0]
                        break
                    header_lines.append(line.rstrip())
                    self.entrystart = filehandler.tell()
                    line = filehandler.readline()
                index_mvf = kwargs.get('contigindex', False)
                if os.path.exists(self.path + ".idx"):
                    # Checks if index file is newer than mvf file
                    if os.path.getmtime(self.path + '.idx') > (
                            os.path.getmtime(self.path)):
                        print("Index file is newer than source, skipping indexing...")
                        index_mvf = False
                if index_mvf is True:
                    line = filehandler.readline()
                    with open(self.path + ".idx", "w") as idxfile:
                        idxfile.write(previous_contig + "\t" +
                                      str(self.entrystart) + "\n")
                        while line:
                            coord = filehandler.tell()
                            line = filehandler.readline()
                            if line.split(":")[0] != previous_contig:
                                if ":" in line:
                                    idxfile.write(line.split(":")[0] + "\t" +
                                                  str(coord) + "\n")
                                    previous_contig = line.split(":")[0]
                self._process_header(header_lines)
                if kwargs.get('contigindex', False) is True:
                    self.read_index_file()
                # Establish number of columns
                self.metadata['ncol'] = len(self.metadata['labels'])
            else:
                raise IOError("MVF path {} not found!".format(self.path))
        # WRITE MODE
        elif filemode in ('write', 'w', 'wb'):
            if os.path.exists(self.path) and (
                    not kwargs.get('overwrite', False)):
                raise IOError(
                    """MVF path {} already exists, use --overwrite
                    to replace""".format(self.path))
            filehandler = (gzip.open(self.path, 'wt') if
                           self.metadata.get('isgzip', False) else
                           open(self.path, 'wt'))
            self.metadata['ncol'] = kwargs.get('ncol', 2)
        filehandler.close()
        self.flavor = kwargs['flavor'] if 'flavor' in kwargs else self.flavor

    def _process_header(self, headerlines):
        """Processes header lines into metadata
            Arguments:
                headerlines: list of unprocessed header strings
        """
        sample_index = 0
        for line in headerlines:
            try:
                entry = line.split()
                # Primary MVF header line
                if entry[0].startswith('##mvf'):
                    for elem in entry[1:]:
                        if '=' not in elem:
                            self.metadata[elem] = True
                            continue
                        elem = elem.split('=')
                        if elem[0] == "flavor":
                            self.flavor = elem[1]
                            if self.flavor not in (
                                    'dna', 'rna', 'prot', 'codon'):
                                raise(RuntimeError(
                                    "ERROR: flavor '{}' is not valid!".format(
                                        self.flavor)))
                        self.metadata[elem[0]] = (
                            interpret_param(elem[1]))
                # Sample column information header lines
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
                # Contig information header lines
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
                # Tree/phylogeny header lines
                elif entry[0].startswith("#t"):
                    self.metadata['trees'].append(line[2:].strip())
                # Notes header lines
                elif entry[0].startswith("#n"):
                    self.metadata['notes'].append(line[2:].strip())
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
            return self.metadata['labels'][indices]
        except IndexError:
            raise IndexError(indices, "contains invalid label index")

    def get_contig_ids(self, labels=None):
        """Returns contig id given contig label
        """
        if labels is None:
            return [x for x in self.metadata['contigs']]
        try:
            if isinstance(labels, (list, tuple, set)):
                return [str(x) for x in self.metadata['contigs']
                        if self.metadata['contigs'][x]['label'] in labels]
            if isinstance(labels, (str, int)):
                return [x for x in self.metadata['contigs']
                        if self.metadata['contigs'][x]['label'] == str(labels)]
            raise TypeError("contig labels not correct datatype")
        except IndexError:
            raise IndexError("contig labels '{}' not found".format(labels))

    def get_contig_labels(self, ids=None):
        """Returns contig label given contig id
        """
        if ids is None:
            return [self.metadata['contigs'][x]['label']
                    for x in self.metadata['contigs']]
        try:
            if isinstance(ids, (list, tuple, set)):
                return [self.metadata['contigs'][x]['label'] for x in ids]
            if isinstance(ids, (str, int)):
                return self.metadata['contigs'][str(ids)]['label']
            raise TypeError("contig labels not correct datatype")
        except IndexError:
            raise IndexError("contig ids '{}' not found".format(ids))

    def get_contig_reverse_dict(self):
        """Returns a dict[label] = id for contigs
        """
        return dict((self.metadata['contigs'][x]['label'], x)
                    for x in self.metadata['contigs'])

    def reset_max_contig_id(self):
        """Requeries the max contig id after modification
        """
        maxid = 0
        if self.metadata['contigs']:
            maxid = max([int(contigid) if is_int(contigid) else 0 for
                         contigid in self.metadata['contigs']])
            self.metadata['maxcontigid'] = maxid
        return ''

    def reset_ncol(self):
        """Resets the ncol metadata after modifying columns
        """
        self.metadata['ncol'] = len(self.metadata['samples'])
        return ''

    def get_next_contig_id(self):
        """Returns the (highest integer id) + 1 or 0"""
        self.metadata['maxcontigid'] += 1
        return str(self.metadata['maxcontigid'])

    def read_index_file(self):
        """Reads the mvf.idx file with contig coordinates
        """
        with open(self.path + ".idx") as idxfile:
            for line in idxfile:
                entry = line.rstrip().split("\t")
                self.metadata['contigs'][entry[0]]['idx'] = int(entry[1])
        return ''

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
            except ValueError:
                raise RuntimeError(
                    "Error processing MVF at line# {} = {} ".format(
                        linecount, line))
        filehandler.close()


    def itercontigentries(self, target_contig, decode=True, no_invariant=False,
                          no_gap=False, no_ambig=False,
                          onlyalleles=False, subset=None):
        """
        Fully-optioned iterator for MVF entries that only returns a
        specific contig.  Intended to be used with indexed mode.
        Returns (str(chrom), int(pos), list(allele entries))
                 or list(allele entries) with 'onlyalleles'

        Arguments:
            decode:fully decode the allele sets (T/F)
            no_invariant: set to false to skip invariant sites
            no_ambig: set to false to skip positions with 'N'
            no_gap: set to false to skip positions with '-'
            onlyalleles: return only list of alleles
            quiet: suppress progress meter (default=False)
            subset: list of column indices

        Note: for codons, filters must apply to all allele strings
        Note: using subset without decode returns encoded subset
        """
        subset = subset or ''
        linecount = 0
        if self.metadata['isgzip']:
            filehandler = gzip.open(self.path, 'rb')
        else:
            filehandler = open(self.path, 'r')
        filehandler.seek(self.metadata['contigs'][target_contig]['idx'])
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
                if contigid != target_contig:
                    break
                if subset:
                    try:
                        allelesets = [''.join([alleles[j] for j in subset])
                                      for alleles in [self.decode(x)
                                                      for x in allelesets]]
                    except IndexError:
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


    def iterentries(self, decode=True, contigs=None, no_invariant=False,
                    no_gap=False, no_ambig=False, no_nonref=False,
                    onlyalleles=False, subset=None):
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
        if contigs is None:
            if no_nonref:
                contigs = sorted([x for x in self.metadata['contigs']
                                  if self.metadata['contigs'][x].get(
                                      'ref', False)])
            # This was turned off in order to speed up
            #else:
               # contigs = sorted(self.metadata['contigs'].keys())
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
                # The below line was added to speed up checking when you
                # are looking at all contigs with a large number
                if contigs is not None:
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
                    except IndexError:
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
            self.flavor, ' '.join([
                "{}={}".format(k, v) for (k, v) in sorted(
                    self.metadata.items())
                if k not in ('contigs', 'samples', 'flavor', 'version',
                             'isgzip', 'labels', 'trees', 'notes')]))]
        header.extend(["#s {} {}".format(
            self.metadata['samples'][x]['label'],
            ' '.join(["{}={}".format(k, v) for (k, v) in (
                sorted(self.metadata['samples'][x].items())) if k != 'label']))
                       for x in range(len(self.metadata['samples']))])
        contigs = [(int(k) if is_int(k) else k, v)
                   for (k, v) in self.metadata['contigs'].items()]
        header.extend(["#c {} label={} length={} {}".format(
            cid, cdata['label'], cdata['length'],
            ' '.join(["{}={}".format(k, v) for k, v in (
                sorted(cdata.items(), key=mixed_sorter))
                      if k not in ['length', 'label', 'idx']]))
                       for cid, cdata in sorted(contigs, key=mixed_sorter)])
        if self.metadata["trees"]:
            header.extend(["#t {}".format(x) for x in self.metadata["trees"]])
        if self.metadata["notes"]:
            header.extend(["#n {}".format(x) for x in self.metadata["notes"]])
        return '\n'.join(header) + '\n'

    def decode(self, alleles):
        """Decode entry into full-length alleles
            Arguments:
                alleles = encoded allele string
        """
        ncol = self.metadata['ncol']
        return decode_mvfstring(alleles, ncol)

    def encode(self, alleles):
        """Internal copy of encode_mvfstring
            Arguments:
                alleles: unencoded allele string
        """
        if self.flavor == 'dna':
            return encode_mvfstring(alleles).replace(
                'N', 'X').replace('n', 'X')
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
            entry[0], entry[1], ' '.join([
                x if encoded else encode_mvfstring(x) for x in entry[2]]))
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
            alleles[0], alleles[1] if alleles[1] != '-' else '',
            pos[0][1], pos[0][0] + 1)
    elif all(x == alleles[2] for x in alleles[3:]):
        alleles = "{}{}+{}1".format(
            alleles[0], alleles[2] if alleles[2] != '-' else '', alleles[1])
    if denovo:
        return '@' + alleles
    return alleles


def decode_mvfstring(alleles, ncol):
    """Decode entry into full-length alleles
        Arguments:
            alleles = encoded allele string
    """
    ref = True
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
        with open(self.path, 'wt') as outfile:
            outfile.write('#' + '\t'.join(self.headers) + "\n")
        return ''

    def write_entry(self, entry):
        """Writes entry to file
            Arguments:
                entry: dict of values with keys matching header
        """
        with open(self.path, 'at') as outfile:
            outfile.write("\t".join([str('.' if k not in entry else entry[k])
                                     for k in self.headers]) + "\n")
        return ''

    def write(self, data):
        """Writes data to file
        """
        with open(self.path, 'at') as outfile:
            outfile.write(data)
        return ''


class AnalysisModule():
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
    print("This is the MVF base function library. "
          "Please run one of the other MVFtools scripts to access these "
          "functions")
