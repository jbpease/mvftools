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
from copy import deepcopy

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
        sample_indices:
    """

    def __init__(self, path, filemode, **kwargs):
        self.path = os.path.abspath(path)
        if filemode not in ('read', 'r', 'rb', 'write', 'w', 'wb'):
            raise RuntimeError("Invalid filemode {}".format(filemode))
        self.filemode = filemode
        self.flavor = 'dna'
        self.metadata = {}
        self.sample_indices = []
        self.sample_ids = []
        self.max_sample_index = 0
        self.sample_id_to_index = {}
        self.sample_data = {}
        self.contig_indices = []
        self.contig_ids = []
        self.contig_id_to_index = {}
        self.contig_labels = []
        self.contig_label_to_index = {}
        self.contig_data = {}
        self.max_contig_index = 0
        self.max_contig_id = 0
        self.trees = []
        self.notes = []
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
                self._process_header(header_lines)
                self.max_sample_index = len(self.sample_indices)
                index_mvf = kwargs.get('contigindex', False)
                if index_mvf is True:
                    if os.path.exists(self.path + ".idx"):
                        # Checks if index file is newer than mvf file
                        if os.path.getmtime(self.path + '.idx') > (
                                os.path.getmtime(self.path)):
                            print("Index file is newer than source, "
                                  "skipping indexing...")
                            index_mvf = False
                if index_mvf is True:
                    line = filehandler.readline()
                    with open(self.path + ".idx", "w") as idxfile:
                        idxfile.write("{}\t{}\n".format(previous_contig,
                                                        self.entrystart))
                        while line:
                            coord = filehandler.tell()
                            line = filehandler.readline()
                            if ":" in line:
                                if line.split(":")[0] != previous_contig:
                                    idxfile.write('{}\t{}\n'.format(
                                        line.split(":")[0],
                                        coord))
                                    previous_contig = line.split(":")[0]
                if kwargs.get('contigindex', False) is True:
                    self.read_index_file()
                # Establish number of columns
                self.max_sample_index = len(self.sample_indices)
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
            self.max_sample_index = kwargs.get('ncol', 2)
        filehandler.close()
        # Final check on flavor after processing headers
        self.flavor = kwargs['flavor'] if 'flavor' in kwargs else self.flavor

    def _process_header(self, headerlines):
        """Processes header lines into metadata
            Arguments:
                headerlines: list of unprocessed header strings
        """
        sample_index = 0
        contig_index = 0
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
                    self.sample_indices.append(sample_index)
                    self.sample_id_to_index[entry[1]] = sample_index
                    self.sample_data[sample_index] = {}
                    self.sample_data[sample_index]['id'] = entry[1]
                    self.sample_ids.append(entry[1])
                    for elem in entry[2:]:
                        if '=' in elem:
                            elem = elem.split('=')
                            self.sample_data[sample_index][elem[0]] = (
                                interpret_param(elem[1]))
                        else:
                            self.sample_data[sample_index][elem] = True
                    sample_index += 1
                # Contig information header lines
                elif entry[0].startswith('#c'):
                    self.contig_indices.append(contig_index)
                    self.contig_id_to_index[entry[1]] = contig_index
                    self.contig_data[contig_index] = {}
                    self.contig_data[contig_index]['id'] = entry[1]
                    self.contig_ids.append(entry[1])
                    for elem in entry[2:]:
                        if '=' in elem:
                            elem = elem.split('=')
                            self.contig_data[contig_index][elem[0]] = (
                                interpret_param(elem[1]))
                        else:
                            self.contig_data[contig_index][elem] = True
                    if 'label' in self.contig_data[contig_index]:
                        self.contig_label_to_index[
                            self.contig_data[contig_index]['label']] = (
                                contig_index + 0)
                    self.contig_labels.append(
                        self.contig_data[contig_index].get('label', None))
                    contig_index += 1
                # Tree/phylogeny header lines
                elif entry[0].startswith("#t"):
                    self.trees.append(line[2:].strip())
                # Notes header lines
                elif entry[0].startswith("#n"):
                    self.notes.append(line[2:].strip())
                else:
                    raise ValueError
            except ValueError:
                sys.stderr.write("Skipping invalid header line '{}'\n".format(
                    entry))
        # self.sample_indices = tuple(self.sample_indices)
        # self.contig_indices = tuple(self.contig_indices)
        return ''

    def __repr__(self):
        ret = ("mvf.metadata = {}\n"
               "mvf.flavor = {}\n"
               "------------\n"
               "mvf.sample_indices = {}\n"
               "mvf.sample_ids = {}\n"
               "mvf.sample_data = {}\n"
               "mvf.max_sample_index = {}\n"
               "mvf.sample_id_to_index = {}\n"
               "------------\n"
               "mvf.contig_indices = {}\n"
               "mvf.contig_ids = {}\n"
               "mvf.contig_id_to_index = {}\n"
               "mvf.contig_labels = {}\n"
               "mvf.contig_label_to_index = {}\n"
               "mvf.contig_data = {}\n"
               "mvf.max_contig_index = {}\n"
               "mvf.max_contig_id = {}\n"
               "------------\n"
               "mvf.trees = {}\n"
               "------------\n"
               "mvf.notes = {}\n"
               "------------\n"
               "mvf.entrystart = {}\n"
               ).format(
                   self.metadata,
                   self.flavor,
                   self.sample_indices,
                   self.sample_ids,
                   self.sample_data,
                   self.max_sample_index,
                   self.sample_id_to_index,
                   self.contig_indices,
                   self.contig_ids,
                   self.contig_id_to_index,
                   self.contig_labels,
                   self.contig_label_to_index,
                   self.contig_data,
                   self.max_contig_index,
                   self.max_contig_id,
                   self.trees,
                   self.notes,
                   self.entrystart
                   )
        return ret

    def copy_headers_from(self, othermvf):
        """Copy headers from one MVF to another"""
        self.filemode = othermvf.filemode
        self.flavor = othermvf.flavor
        self.metadata = othermvf.metadata.copy()
        self.sample_indices = othermvf.sample_indices[:]
        self.sample_ids = othermvf.sample_ids[:]
        self.max_sample_index = othermvf.max_sample_index + 0
        self.sample_id_to_index = othermvf.sample_id_to_index.copy()
        self.sample_data = deepcopy(othermvf.sample_data)
        self.contig_indices = othermvf.contig_indices[:]
        self.contig_ids = othermvf.contig_ids[:]
        self.contig_id_to_index = othermvf.contig_id_to_index.copy()
        self.contig_labels = othermvf.contig_labels[:]
        self.contig_label_to_index = othermvf.contig_label_to_index.copy()
        self.contig_data = deepcopy(othermvf.contig_data)
        self.max_contig_index = othermvf.max_contig_index + 0
        self.max_contig_id = othermvf.max_contig_id + 0
        self.trees = othermvf.trees[:]
        self.notes = othermvf.notes[:]

    def get_sample_indices(self, ids=None):
        """Get indices for a specified named group or set of labels
           Arguments:
                labels: list/set of str or str; default= all
        """
        if ids is None:
            return tuple(self.sample_indices)
        try:
            if isinstance(ids, str):
                return self.sample_id_to_index[ids]
            if hasattr(ids, '__iter__'):
                return [self.sample_id_to_index[x] for x in ids]
            return self.sample_id_to_index[ids]
        except IndexError:
            raise IndexError(ids, "contains invalid ids")

    def get_sample_ids(self, indices=None):
        """Get labels for the specified named indices
            Arguments:
                indices: list/set of int or int; default= all
        """
        if indices is None:
            return tuple(self.sample_ids)
        try:
            if hasattr(indices, '__iter__'):
                return [self.sample_data[x]['id']for x in indices]
            return self.sample_data[indices]['id']
        except IndexError:
            raise IndexError(indices, "contains invalid index")

    def get_contig_indices(self, ids=None, labels=None):
        """Returns contig id given contig label
        """
        if labels is None and ids is None:
            return self.contig_indices
        if labels is not None and ids is not None:
            raise RuntimeError("cannot specify both ids and labels for "
                               "get_contig_indices()")
        if ids is not None:
            try:
                if isinstance(ids, str):
                    return self.contig_id_to_index[ids]
                if hasattr(ids, '__iter__'):
                    return [self.contig_id_to_index[x] for x in ids]
                return self.contig_id_to_index[ids]
            except KeyError:
                raise KeyError("contig id(s) '{}' not found".format(ids))
        if labels is not None:
            try:
                if isinstance(labels, str):
                    return self.contig_label_to_index[labels]
                if hasattr(labels, '__iter__'):
                    return [self.contig_label_to_index[x] for x in labels]
                return self.contig_label_to_index[labels]
            except KeyError:
                raise KeyError("contig label(s) '{}' not found".format(labels))
        return ''

    def get_contig_ids(self, indices=None, labels=None):
        """Returns contig id given contig label
        """
        if labels is None and indices is None:
            return tuple(self.contig_ids)
        if labels is not None and indices is not None:
            raise RuntimeError("cannot specify both indices and labels for "
                               "get_contig_ids()")
        if indices is not None:
            try:
                if hasattr(indices, '__iter__'):
                    return tuple(self.contig_data[x]['id'] for x in indices)
                return (self.contig_data[indices]['id'], )
            except KeyError:
                raise KeyError("contig index '{}' not found".format(indices))
        elif labels is not None:
            try:
                if isinstance(labels, str):
                    return self.contig_data[
                        self.contig_label_to_index[labels]]['id']
                if hasattr(labels, '__iter__'):
                    return [
                        self.contig_data[self.contig_label_to_index[x]]['id']
                        for x in labels]
                return self.contig_data[
                    self.contig_label_to_index[labels]]['id']
            except KeyError:
                raise KeyError("contig label(s) '{}' not found".format(labels))
        return ''

    def get_contig_labels(self, indices=None, ids=None):
        """Returns contig id given contig label
        """
        if ids is None and indices is None:
            return tuple(self.contig_labels)
        if ids is not None and indices is not None:
            raise RuntimeError("cannot specify both indices and labels for "
                               "get_contig_ids()")
        if indices is not None:
            try:
                if hasattr(indices, '__iter__'):
                    return [self.contig_data[x]['label'] for x in indices]
                return [self.contig_data[indices]['label'], ]
            except KeyError:
                raise KeyError("contig index '{}' not found".format(indices))
        if ids is not None:
            try:
                if isinstance(ids, str):
                    return self.contig_data[self.contig_id_to_index[ids]]['label']
                if hasattr(ids, '__iter__'):
                    return [self.contig_data[
                        self.contig_id_to_index[x]]['label']
                            for x in ids]
                return self.contig_data[self.contig_id_to_index[ids]]['label']
            except KeyError:
                raise KeyError("contig id(s) '{}' not found".format(ids))
        return ''

    def reset_max_contig(self):
        """Requeries the max contig id after modification
        """
        maxid = 0
        if self.contig_ids:
            maxid = max([int(contigid) if is_int(contigid) else 0 for
                         contigid in self.contig_ids])
        self.max_contig_id = maxid
        self.max_contig_index = len(self.contig_indices)
        return ''

    def reset_max_sample(self):
        """Resets the ncol metadata after modifying columns
        """
        self.max_sample_index = len(self.sample_indices)
        self.metadata['ncol'] = self.max_sample_index + 0
        return ''

    def get_next_contig_id(self):
        """Returns the (highest integer id) + 1 or 0"""
        self.max_contig_id += 1
        return str(self.max_contig_id)

    def get_next_contig_index(self):
        """Returns the (highest integer id) + 1 or 0"""
        self.max_contig_index += 1
        return self.max_contig_index

    def read_index_file(self):
        """Reads the mvf.idx file with contig coordinates
           Index File stores contig cordinates based on id, not index
        """
        with open(self.path + ".idx") as idxfile:
            for line in idxfile:
                entry = line.rstrip().split("\t")
                self.contig_data[
                    self.contig_id_to_index[entry[0]]]['idx'] = int(entry[1])
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

    def itercontigentries(self, target_contig_index,
                          decode=True, no_invariant=False,
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
        if 'idx' not in self.contig_data[target_contig_index]:
            return
        filehandler.seek(self.contig_data[target_contig_index]['idx'])
        target_id = self.contig_data[target_contig_index]['id']
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
                if contigid != target_id:
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
                    if len(allelesets[0]) > 1 and (
                            allelesets[0][1] == '+' and
                            allelesets[0][2] == allelesets[0][0]):
                        continue
                    if subset:
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
            except Exception as exc:
                raise RuntimeError(
                    "Error processing MVF at line# {} = {}\n{}".format(
                        linecount, line, exc))
        filehandler.close()

    def iterentries(self, decode=True,
                    contig_ids=None,
                    contig_labels=None,
                    contig_indices=None,
                    no_invariant=False,
                    no_gap=False, no_ambig=False,
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
            onlyalleles: return only list of alleles
            quiet: suppress progress meter (default=False)
            subset: list of column indices

        Note: for codons, filters must apply to all allele strings
        Note: using subset without decode returns encoded subset
        """
        contigs = None
        if contig_labels is not None:
            contigs = [self.contig_data[self.contig_label_to_index[x]]['id']
                       for x in contig_labels]
        elif contig_ids is not None:
            contigs = list(contig_ids)
        elif contig_indices is not None:
            contigs = [self.contig_data[x]['id'] for x in contig_indices]
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
                    if len(allelesets[0]) > 1 and (
                            allelesets[0][1] == '+' and
                            allelesets[0][2] == allelesets[0][0]):
                        continue
                    if subset:
                        if all(x == allelesets[0][0]
                               for x in allelesets[0][1:]):
                            continue
                if subset and decode is False:
                    allelesets = [self.encode(x) for x in allelesets]
                if decode is True and not subset:
                    allelesets = [self.decode(x) for x in allelesets]
                if onlyalleles:
                    yield allelesets
                else:
                    yield (contigid, pos, allelesets)
            except Exception as exc:
                raise RuntimeError(
                    "Error processing MVF at line# {} = {}\n{}".format(
                        linecount, line, exc))
        filehandler.close()

    def get_header(self):
        """Returns formatted header string (with final newline)
        """
        header = ["##mvf version=1.2 flavor={} {}".format(
            self.flavor, ' '.join([
                "{}={}".format(k, v) for (k, v) in self.metadata.items()
                if k not in ('version', 'isgzip', 'labels', 'flavor',
                             'trees', 'notes')]))]
        if self.notes:
            header.extend(["#n {}".format(x) for x in self.notes])
        if self.trees:
            header.extend(["#t {}".format(x) for x in self.trees])
        header.extend(["#s {} {}".format(
            self.sample_data[x]['id'],
            ' '.join(["{}={}".format(k, v) for (k, v) in (
                self.sample_data[x].items()) if k != 'id']))
                       for x in self.sample_indices])
        header.extend(["#c {} label={} length={} {}".format(
            self.contig_data[x]['id'], self.contig_data[x]['label'],
            self.contig_data[x]['length'],
            ' '.join(["{}={}".format(k, v) for k, v in (
                self.contig_data[x].items())
                      if k not in ['length', 'label', 'idx', 'id']]))
                       for x in self.contig_indices])
        return '\n'.join(header) + '\n'

    def copy_header(self, mvfold):
        """Returns formatted header string (with final newline)
        """
        self.flavor = mvfold.flavor
        self.metadata = mvfold.metadata.copy()
        self.sample_indices = mvfold.sample_indices[:]
        self.sample_ids = mvfold.sample_ids[:]
        self.max_sample_index = len(self.sample_indices)
        self.sample_id_to_index = mvfold.sample_id_to_index.copy()
        self.sample_data = mvfold.sample_data.copy()
        self.contig_indices = mvfold.contig_indices[:]
        self.max_contig_index = max(self.contig_indices)
        self.contig_ids = mvfold.contig_ids[:]
        self.contig_id_to_index = mvfold.contig_id_to_index.copy()
        self.contig_labels = mvfold.contig_labels.copy()
        self.contig_label_to_index = mvfold.contig_label_to_index.copy()
        self.contig_data = mvfold.contig_data.copy()
        self.trees = mvfold.trees
        self.notes = mvfold.notes

    def decode(self, alleles):
        """Decode entry into full-length alleles
            Arguments:
                alleles = encoded allele string
        """
        return decode_mvfstring(alleles, self.max_sample_index)

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
        if entries:
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

    def write_entry(self, entry, defaultvalue="."):
        """Writes entry to file
            Arguments:
                entry: dict of values with keys matching header
        """
        with open(self.path, 'at') as outfile:
            outfile.write("\t".join([str(defaultvalue if
                                         k not in entry else entry[k])
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
