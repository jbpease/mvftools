#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This program creates a chromoplot from an MVF alignment.
A chromoplot shows a genome-wide diagram of different
evolutionary histories for a given quartet of taxa.
"""

import os
import sys
import argparse
from itertools import combinations
from mvfbase import MultiVariantFile
from scipy.stats import chi2

_LICENSE = """
If you use this software please cite:
Pease JB and BK Rosenzweig. 2016.
"Encoding Data Using Biological Principles: the Multisample Variant Format
for Phylogenomics and Population Genomics"
IEEE/ACM Transactions on Computational Biology and Bioinformatics. In press.
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


class Counter(dict):
    """Extension of dictionary functions for a sinple counter"""

    def _sorted_items_list(self):
        """Generate a list of tuples with sorted values"""
        return sorted([(v, k) for k, v in self.items()], reverse=True)

    def add(self, key, val=1):
        """Add item to counter, creates key if DNE"""
        self[key] = self.get(key, 0) + val
        return ''

    def get_maximum(self, target='value'):
        """Get maximum value or key"""
        sortitems = self._sorted_items_list()
        if sortitems[0][0] == sortitems[1][0]:
            return -1
        if target == 'key':
            return sortitems[0][1]
        else:
            return sortitems[0][0]

    def get_ranked_keys(self, mode='descending'):
        """Get keys in order of values"""
        sortitems = self._sorted_items_list()
        if mode == 'centered':
            return [sortitems[2][1], sortitems[0][1], sortitems[1][1]]
        return [x[1] for x in sortitems]


class Pallette(object):
    """Class for color management for chromoplots"""

    def __init__(self, basecolors=None):
        self.colornames = {
            'lgrey': (250, 250, 250), 'dgrey': (192, 192, 192),
            'black': (0, 0, 0), 'white': (255, 255, 255),
            'red': (192, 0, 0), 'orange': (217, 95, 2),
            'yellow': (192, 192, 0), 'green': (0, 192, 0),
            'blue': (0, 0, 192), 'teal': (27, 158, 119),
            'puce': (117, 112, 179), 'purple': (192, 0, 192),
            'none': ()}
        self.basecolors = (basecolors
                           if basecolors is not None
                           else ['puce', 'orange', 'teal'])

    def color_str(self, color):
        """Return ASCII color codes for numerical trio"""
        return (self.colornames[color][0],
                self.colornames[color][1],
                self.colornames[color][2],
                255)

    def solid_colortracks(self, values, infotrack=False):
        """Majority single-color tracks"""
        ret = [(self.color_str(self.basecolors[j]) if x > 0 else (
            self.color_str('white'))) for j, x in enumerate(values)] + ([
                self.color_str('dgrey')] * infotrack)
        return ret

    def shaded_colortracks(self, values, infotrack=False):
        """Make shaded color using relative intensities"""
        total = float(sum(values))
        if not values or not total:
            colortracks = [self.color_str('white')]*(3 + int(infotrack))
        else:
            colortracks = [
                tuple([int(255 - ((255 - x) * (values[j]/total)))
                       for x in self.colornames[self.basecolors[j]]] + [255])
                for j in range(len(values))]
            if infotrack:
                colortracks.append(self.color_str('dgrey'))
        return colortracks


class Chromoplot(object):
    """Chromomplot data container and processor"""

    def __init__(self, **kwargs):
        self.data = kwargs.get('data', 0) or {}
        self.counts = kwargs.get('counts', 0) or {}
        self.pallette = kwargs.get('pallette', 0) or Pallette()
        self.params = kwargs.get('params', 0) or {}
        self.params.update([(k, v) for k, v in kwargs.items()
                            if k not in ['data', 'params', 'counts'
                                         'pallette']])
        self.datalog = kwargs.get('datalog', 0) or []
        if 'data' not in kwargs:
            for param in ['contigs', 'windowsize']:
                if param not in self.params:
                    raise RuntimeError("'{}' not in params".format(param))
            for contigid, _, length in self.params['contigs']:
                self.counts[contigid] = {}
                self.data[contigid] = dict([
                    (start, dict())
                    for start in range(
                        0, int(length // self.params['windowsize']) + 1)])
        self.params['winlogpath'] = self.params.get(
            'winlogpath', self.params['outpath'] + '.windows.log')
        self.params['totallogpath'] = self.params.get(
            'totallogpath', self.params['outpath'] + '.summary.log')

    def add_data(self, contig, window, site_code, val=1):
        """Adds site code data to chromoplot object"""
        if contig not in self.data:
            self.data[contig] = {}
        if window not in self.data[contig]:
            self.data[contig][window] = {}
        self.data[contig][window][site_code] = self.data[contig][window].get(
            site_code, 0) + val
        if contig not in self.counts:
            self.counts[contig] = {}
        self.counts[contig][site_code] = self.counts[contig].get(
            site_code, 0) + val
        return ''

    def interpret_colors(self, valuelist):
        """Create color from list of values"""
        colortracks = []
        ntracks = self.params['ntracks']
        infotrack = self.params['infotrack']
        reorder = []
        if self.params['track_order'] != [1, 2, 3]:
            reorder = [[1, 2, 3].index(x) for x in self.params['track_order']]
        for entry in valuelist:
            if -1 in entry:
                colortracks.append([self.pallette.color_str('lgrey')] * (
                    ntracks + int(infotrack)))
            else:
                if reorder:
                    entry = [entry[x] for x in reorder] + [ntracks]*infotrack
                if self.params['majority']:
                    colortracks.append(self.pallette.solid_colortracks(
                        entry, infotrack=infotrack))
                else:
                    colortracks.append(self.pallette.shaded_colortracks(
                        entry, infotrack=infotrack))
        return b''.join([b''.join(bytes(y) for y in x) *
                         self.params['yscale'] for x in zip(*colortracks)])

    def parse_count_trio(self, codes):
        """Parse Site Counts from Trios"""
        majority_count = 0
        if codes == 'nodata':
            values = [-1]*(3+int(self.params['infotrack']))
        elif all(x not in codes for x in (12, 10, 6)):
            if any(x in codes for x in (8, 4, 2)):
                values = [0]*3 + [0.5] * self.params['infotrack']
            else:
                values = [0]*3 + [1] * self.params['infotrack']
        else:
            tree0 = codes.get(12, 0)  # 12 = 1100 = BBAA
            tree1 = codes.get(10, 0)  # 10 = 1010 = BABA
            tree2 = codes.get(6, 0)   # 6  = 0110 = ABBA
            if tree0 > tree1 and tree0 > tree2:
                values = [1, 0, 0] + [0.5] * self.params['infotrack']
                majority_count = 1
            elif tree1 > tree0 and tree1 > tree2:
                values = [0, 1, 0] + [0.5] * self.params['infotrack']
                majority_count = 2
            elif tree2 > tree0 and tree2 > tree1:
                values = [0, 0, 1] + [0.5] * self.params['infotrack']
                majority_count = 3
            else:
                values = [0, 0, 0] + [0.5] * self.params['infotrack']
            if not self.params['majority']:
                values = [tree0, tree1, tree2]
                total = sum(values)
                if total:
                    values = [float(x) / total for x in values]
        return majority_count, values

    def plot_chromoplot(self):
        """Make Chromoplot for count-based trio"""
        maxlen = max([x[2] for x in self.params['contigs']])
        width = int(maxlen // self.params['windowsize']) + 1
        self.params['ntracks'] = [0, 0, 0, 0, 3, 15][4]
        # total_windows = (maxlen // self.params['windowsize'] + 1) * len(
        #    self.data)
        majority_counts = Counter([(1, 0), (2, 0), (3, 0)])
        contig_ab_values = {}
        self.write_window_log(headermode=True)

        for contig, _, _ in self.params['contigs']:

            i = 0
            contig_ab_values[contig] = []
            while i < width:
                window_codes = self.data[contig].get(i, 'nodata')
                if self.params['ntracks'] == 3:
                    majority_call, trio = self.parse_count_trio(window_codes)
                    contig_ab_values[contig].append(trio)
                    if window_codes != 'nodata':
                        if majority_call:
                            majority_counts.add(majority_call)
                        self.datalog.append((window_codes, contig,
                                             i * self.params['windowsize']))
                i += 1
        self.params['track_order'] = majority_counts.get_ranked_keys(
            mode='centered')

        for entry in self.datalog:
            self.write_window_log(window_codes=entry[0], contig=entry[1],
                                  pos=entry[2],
                                  track_order=self.params['track_order'])
        if self.params['plottype'] == "graph":
            from matplotlib import pyplot as plt
            import numpy as np
            maxlen = max([len(contig_ab_values[x]) for x in contig_ab_values])
            fig, axl = plt.subplots()
            cvals = np.empty((3*len(self.params['contigs']), maxlen))
            for i, (contig, _, _) in enumerate(self.params['contigs']):
                c_ab_vals = contig_ab_values[contig]
                for j in self.params['track_order']:
                    for k, valk in enumerate(c_ab_vals[j]):
                        cvals[i + k, j] = valk
            axl.pcolormesh(cvals)
            plt.savefig(self.params['outpath'])
        else:
            image_data = bytes()
            for icontig, (contig, _, _) in enumerate(self.params['contigs']):
                image_data += self.interpret_colors(contig_ab_values[contig])
                if icontig != len(self.params['contigs']) - 1:
                    image_data = image_data + b'\x00\x00\x00\xff'*width
            height = len(image_data) / (width * (
                1 + self.params['ntracks'] + int(self.params['infotrack'])))
            write_png(self.params['outpath'],
                      image_data, width, height)
        return ''

    def write_window_log(self, window_codes=None, contig=None, pos=None,
                         headermode=False, track_order=None):
        """Write entry for window data
            Arguments:
                window_codes: list of data from window
                contig: string contig name
                position: int position on contig
                headermode: invoke write headers instead of data
                track_order: list of track orders to keep majority in center
        """
        if not track_order:
            track_order = [1, 2, 3]
        headers = ['taxon1', 'taxon2', 'taxon3', 'outgroup',
                   'contig', 'winstart', 'windowcall', 'totalsites',
                   'ambiguous', 'nonpolar', 'triallelic', 'gap',
                   'BBAA', 'ABBA', 'BABA', 'Dleft', 'Dright',
                   'Dsites', 'Dstat', 'Pvalue',
                   'pBBAA', 'pABBA', 'pBABA', 'Dorder']
        abcodes = ['BBAA', 'BABA', 'ABBA']
        if headermode:
            with open(self.params['winlogpath'], 'w') as logfile:
                logfile.write('\t'.join(headers) + "\n")
        else:
            with open(self.params['winlogpath'], 'a') as logfile:
                # 12 = BBAA = tree0 = 1
                # 10 = BABA = tree1 = 2
                # 06 - ABBA = tree2 = 3
                leftcount = window_codes.get(
                    [12, 10, 6][track_order[0] - 1], 0)
                rightcount = window_codes.get(
                    [12, 10, 6][track_order[2] - 1], 0)
                dval, pval = dcalc(leftcount, rightcount)
                total_ab = sum([window_codes.get(x, 0) for x in [6, 10, 12]])
                d_order = leftcount + rightcount and "{}-{}".format(
                    abcodes[track_order[0] - 1],
                    abcodes[track_order[2] - 1]) or 'NA'
                logfile.write("{}\n".format('\t'.join([
                    str(x) for x in list(self.params['labels']) +
                    [contig, pos, abcodes[track_order[1] - 1],
                     sum(window_codes.values())] +
                    [window_codes.get(x, 0) for x in [
                        'ambiguous', 'nonpolar', 'triallelic',
                        'gap', 12, 6, 10]] +
                    [leftcount, rightcount, total_ab, dval, pval] +
                    [(round(float(window_codes.get(x, 0))/total_ab, 3)
                      if total_ab > 0 else 0) for x in [6, 10, 12]] +
                    [d_order]
                    ])))
        return ''

    def write_total_log(self):
        """Write summary log data
        """
        nwindows = {'TOTAL': 0}
        contigorder = [x[0] for x in self.params['contigs']]
        self.counts['TOTAL'] = {}
        contigorder.append('TOTAL')
        for contig in self.counts:
            if contig == 'TOTAL':
                continue
            nwindows[contig] = len(self.data[contig])
            for code in self.counts[contig]:
                self.counts['TOTAL'][code] = (
                    self.counts['TOTAL'].get(code, 0)) + (
                        self.counts[contig].get(code, 0))
            nwindows['TOTAL'] += len(self.data[contig])
        with open(self.params['totallogpath'], 'w') as logfile:
            logfile.write("{}\n".format("\t".join(
                ['taxon1', 'taxon2', 'taxon3', 'outgroup',
                 'contig', 'nwindow', 'totalsites', 'absites',
                 'ambiguous', 'nonpolar', 'triallelic', 'gap',
                 'BBAA', 'ABBA', 'BABA', 'Dleft', 'Dright',
                 'Dsites', 'Dstat', 'Pvalue',
                 'pBBAA', 'pABBA', 'pBABA', 'Dorder'
                 ])))
            for contig in contigorder:
                counts = self.counts[contig]
                total = float(sum([counts.get(x, 0) for x in [
                    0, 2, 4, 8, 6, 10, 12, 14, 'invariant']]))
                abba = counts.get(6, 0)
                baba = counts.get(10, 0)
                bbaa = counts.get(12, 0)
                total_ab = float(abba + baba + bbaa)
                leftcount = counts.get(
                    [12, 10, 6][self.params['track_order'][0] - 1], 0)
                rightcount = counts.get(
                    [12, 10, 6][self.params['track_order'][2] - 1], 0)
                d_order = "{}-{}".format(
                    ['BBAA', 'BABA', 'ABBA'][
                        self.params['track_order'][0] - 1],
                    ['BBAA', 'BABA', 'ABBA'][
                        self.params['track_order'][2] - 1])
                dval, pval = dcalc(leftcount, rightcount)
                logfile.write('{}\n'.format("\t".join([
                    str(x) for x in list(self.params['labels']) +
                    [str(contig).zfill(2), nwindows[contig],
                     sum(counts.values()), total] +
                    [counts.get(x, 0) for x in [
                        'ambiguous', 'nonpolar', 'triallelic', 'gap',
                        12, 6, 10]] +
                    [leftcount, rightcount, total_ab, dval, pval] +
                    [(round(float(counts.get(x, 0))/total_ab, 3)
                      if total_ab > 0 else 0)
                     for x in [6, 10, 12]] + [d_order]
                    ])))
        return ''


def dcalc(abba, baba):
    """Calculate the D-statistic and Chi2 P-value
    """
    dval = (float(abba - baba)/(abba + baba)) if abba + baba > 0 else 0
    pval = chi2_test(abba, baba)[1] if abba + baba > 0 else 1
    return dval, pval


def chi2_test(val0, val1):
    """Calculate Pearson Chi-Squared for two values that should be equal
    """
    try:
        chisq = float((val0 - val1)**2) / float(val0 + val1)
        if chisq == 0:
            return (0, 1)
        pval = 1.0 - chi2.cdf(chisq, 1)
        return (chisq, pval)
    except ZeroDivisionError:
        return (0, 1)


def counter_print(dict_counter, reverse=True, percentage=True, n_values=None):
    """Prints a Dictionary with Integer Values, with sorting and percentages
    """
    total = float(sum(dict_counter.values()))
    if not total:
        return "Nothing counted!"
    arr_counter = sorted([(v, k) for (k, v) in dict_counter.items()],
                         reverse=reverse)
    for (val, k) in arr_counter[:n_values]:
        if percentage:
            return "{} {} {}%".format(
                k, val, round(float(val) * 100/total, 2))
        else:
            return "{} {}".format(k, val)
    return ''


def write_png(filepath, buf, width, height):
    """Writes PNG image"""
    import zlib
    import struct
    with open(filepath, 'wb') as imgfile:
        height = int(height)
        width_byte_4 = width * 4
        raw_data = b"".join(bytearray([0]) + buf[span:span + width_byte_4]
                            for span in range(
                                (height - 1) * width * 4, -1, - width_byte_4))
        print(len(raw_data))
        print(width, height, width * height)

        def png_pack(png_tag, data):
            """PNG packaging"""
            chunk_head = png_tag + data
            return struct.pack("!I", len(data)) + chunk_head + struct.pack(
                "!I", 0xFFFFFFFF & zlib.crc32(chunk_head))

        imgfile.write(b"".join([b'\x89PNG\r\n\x1a\n',
                                png_pack(b'IHDR', struct.pack(
                                    "!2I5B", width, height, 8, 6, 0, 0, 0)),
                                png_pack(b'IDAT', zlib.compress(raw_data, 9)),
                                png_pack(b'IEND', b'')]))


def generate_argparser():
    pallette = Pallette()
    parser = argparse.ArgumentParser(
        prog="mvf_chromoplot.py",
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog=_LICENSE)
    parser.add_argument("-i", "--mvf", help="Input MVF file.", required=True,
                        type=os.path.abspath)
    parser.add_argument("-o", "--outprefix",
                        help="Output prefix (not required).")
    parser.add_argument("-s", "--samples", nargs='*', required=True,
                        help="3 or more taxa to use for quartets")
    parser.add_argument("-G", "--outgroup", nargs='*', required=True,
                        help="1 or more outgroups to use for quartets")
    parser.add_argument("-w", "--windowsize", type=int, default=100000)
    parser.add_argument("-c", "--contigs", nargs='*',
                        help=("Enter the ids of one or more contigs in the "
                              "order they will appear in the chromoplot. "
                              "(defaults to all ids in order present in MVF)"))
    parser.add_argument("-M", "--majority", action="store_true",
                        help=("Plot only 100% shading in the majority track "
                              " rather than shaded proportions in all "
                              "tracks."))
    parser.add_argument("-I", "--infotrack", action="store_true",
                        help=("Include an additional coverage information "
                              "track that will show empty, uninformative, "
                              "and informative loci. (Useful for "
                              "ranscriptomes/RAD or other reduced sampling."))
    parser.add_argument("-E", "--emptymask", choices=pallette.colornames,
                        default="none",
                        help="Mask empty regions with this color.")
    parser.add_argument("-y", "--yscale", default=20, type=int,
                        help="Height (in number of pixels) for each track")
    parser.add_argument("-x", "--xscale", default=1, type=int,
                        help="Width (in number of pixels) for each window")
    parser.add_argument("-C", "--colors", nargs=3, choices=pallette.colornames,
                        help="three colors to use for chromoplot")
    parser.add_argument("-q", "--quiet", action="store_true",
                        help="suppress all output messages")
    parser.add_argument("-P", "--plottype", choices=["graph", "image"],
                        default="image",
                        help=("PNG image (default) "
                              "or graph via matplotlib (experimental)"))
    parser.add_argument("-v", "--version", action="version",
                        version="2017-06-24",
                        help="display version information")
    return parser


def main(arguments=None):
    """Main method"""
    arguments = arguments if arguments is not None else sys.argv[1:]
    pallette = Pallette()
    parser = generate_argparser()
    args = parser.parse_args(args=arguments)
    if args.colors is None:
        pallette.basecolors = args.colors
    # Establish MVF and parse chromosome information
    if args.quiet is False:
        print("Reading MVF...")
    mvf = MultiVariantFile(args.mvf, 'read')
    if args.quiet is False:
        print("Parsing headers...")
    if args.contigs is None:
        contignames = args.contigs
    else:
        contignames = [mvf.metadata['contigs'][contigid]['label']
                       for contigid in mvf.metadata['contigs']]
        for i, contigname in enumerate(contignames):
            try:
                contignames[i] = int(contigname)
            except:
                pass
        contignames = [str(x) for x in sorted(contignames)]
    if args.quiet is False:
        print("Plotting chromoplot for contigs: {}".format(
            ",".join(contignames)))
    master_contigs = []
    for contigname in contignames:
        contig_found = False
        for contigid in mvf.metadata['contigs']:
            if (contigname == contigid or
                    contigname == mvf.metadata['contigs'][contigid]['label']):
                master_contigs.append((
                    contigid,
                    mvf.metadata['contigs'][contigid]['label'],
                    mvf.metadata['contigs'][contigid]['length']))
                contig_found = True
        if contig_found:
            continue
        raise RuntimeError(contigname, "not found in MVF contig ids or labels")
    quartets = [(x, y, z, outgroup) for x, y, z in
                combinations(args.samples, 3) for outgroup in args.outgroup]
    # Begin iterations
    for quartet in quartets:
        if args.quiet is False:
            print("Beginning quartet {}".format(",".join(quartet)))
        params = {'contigs': [[str(x), y, z] for [x, y, z] in master_contigs],
                  'outpath': (args.outprefix or '_'.join(quartet)) + ".png",
                  'labels': quartet,
                  'windowsize': args.windowsize,
                  'majority': args.majority,
                  'infotrack': args.infotrack,
                  'yscale': args.yscale,
                  'xscale': args.xscale,
                  'quiet': args.quiet,
                  'plottype': args.plottype}
        chromoplot = Chromoplot(params=params, pallette=pallette)
        quartet_indices = mvf.get_sample_indices(labels=quartet)
        current_contig = ''
        for contig, pos, allelesets in mvf.iterentries(
                subset=quartet_indices, decode=True,
                contigs=[str(x[1]) for x in master_contigs]):
            if contig != current_contig:
                if args.quiet is False:
                    print("Starting contig {}".format(contig))
                    current_contig = contig[:]
            alleles = allelesets[0]
            if '-' in alleles:
                site_code = 'gap'
            elif any(x not in 'ATGCatgc' for x in alleles):
                site_code = 'ambiguous'
            elif alleles[3] not in alleles[:3]:
                site_code = 'nonpolar'
            elif len(set(alleles)) > 2:
                site_code = 'triallelic'
            else:
                site_code = sum([2**(3-j) * (alleles[j] != alleles[3])
                                 for j in range(3)])
            chromoplot.add_data(str(contig),
                                int(pos // args.windowsize), site_code)
        if not args.quiet:
            print("Writing image...")
        chromoplot.plot_chromoplot()

        if not args.quiet:
            print("Writing log...")
        chromoplot.write_total_log()
    return ''


if __name__ == "__main__":
    main()
