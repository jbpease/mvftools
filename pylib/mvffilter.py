#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This program filters an MVF alignment using the modules specified below,
use the --more-help option for additional module information.

MVFtools: Multisample Variant Format Toolkit
James B. Pease and Ben K. Rosenzweig
http://www.github.org/jbpease/mvftools

If you use this software please cite:
Pease, James B. and Benjamin K. Rosenzweig. 2018.
"Encoding Data Using Biological Principles: the Multisample Variant Format
for Phylogenomics and Population Genomics"
IEEE/ACM Transactions on Computational Biology and Bioinformatics.
15(4) 1231–1238.
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

import sys
from pylib.mvfbase import MultiVariantFile, encode_mvfstring
from pylib.mvfbiolib import MvfBioLib

MLIB = MvfBioLib()

# Note action modules are designed in these types:
# filter = returns boolean False to filter OUT, True to retain line
# transform = applies a transformation function to the allele string
# location = applies a transformation to ht
#
# Actions are applied IN ORDER, meaning that filters can be applied
# before and/or after transforms, the order of actions is important
#
# For new modules, the preferred order of 'mvfenc' checks is:
# full, invar, onecov, onevar, refvar
# (descending order of most common frequency)


def get_linetype(alleles):
    """Determines the line type from allele string
    """
    if not alleles:
        return 'empty'
    if len(alleles) == 1:
        linetype = ('empty' if alleles[0] == '-' else
                    'invar')
    elif len(alleles) == 2:
        linetype = ('empty' if alleles[0:1] == '--' else
                    'refvar')
    elif alleles[1] == '+':
        linetype = ('empty' if all(alleles[x] == '-' for x in (0, 2)) else
                    'onecov')
    elif alleles[2] == '+':
        linetype = ('empty' if all(alleles[x] == '-' for x in (0, 1, 3)) else
                    'onevar')
    else:
        linetype = ('empty' if all(x == '-' for x in alleles) else
                    'full')
    return linetype


def make_module(modulename, ncol, optargs=None):
    """Generate Modules"""

    # ALLELEGROUP
    def allelegroup(entry, mvfenc):
        """Only retain sites where all members of the specific
           group contain alleles
        """
        if mvfenc == 'full':
            allele_groups = [set(entry[x] for x in y) -
                             set('X-') for y in optargs]
            if not all(allele_groups):
                return False
        return False

    # COLLAPSEMERGE
    def collapsemerge(entry, mvfenc):
        """Samples merged completely, uses ambiguity codes
           for heterozygous alleles"""
        if mvfenc == 'full':
            newbase = MLIB.merge_bases([entry[x] for x in optargs[0]])
            return ''.join([(j == optargs[0][0] and newbase) or
                            (j not in optargs[0] and entry[j]) or
                            '' for j in range(len(entry))])
        if mvfenc == 'invar':
            return entry
        if mvfenc == 'onecov':
            num = int(entry[3:])
            if optargs[0][0] == 0 and num in optargs[0]:
                return ''
            return "{}{}".format(entry[0:3], num - len([
                x for x in optargs[0] if x < num]))
        if mvfenc == 'onevar':
            num = int(entry[4:])
            if optargs[0][0] == 0:
                if num in optargs[0]:
                    if len(optargs[0]) > 2:
                        entry = "{}{}".format(MLIB.merge_bases([
                            entry[0], entry[1], entry[3]]), entry[1])
                    else:
                        entry = "{}{}".format(MLIB.merge_bases([
                            entry[0], entry[1]]), entry[2:])
            elif optargs[0][0] == num:
                entry = "{}{}{}".format(
                    entry[0:3], MLIB.merge_bases([entry[1], entry[3]]),
                    num)
            elif num in optargs:
                entry = "{}{}{}".format(
                    entry[0:3], MLIB.merge_bases([entry[1], entry[3]]),
                    optargs[0][0])
            else:
                entry = "{}{}".format(
                    entry[:4], 
                    num - len([x for x in optargs[0] if x < num]))
            return entry
        if mvfenc == 'refvar':
            if optargs[0][0] == 0:
                return "{}{}".format(MLIB.merge_bases(entry[0:2]),
                                     entry[1])
        return entry

    # COLLAPSEPRIORITY
    def collapsepriority(entry, mvfenc):
        """sample alleles combined,
           using a priority list if gap encountered"""
        if mvfenc == 'full':
            colbases = [entry[x] for x in optargs[0]
                        if entry[x] not in 'X-'] + ['-']
            return ''.join([(j not in optargs[0] and entry[j]) or
                            (j == optargs[0][0] and colbases[0]) or
                            '' for j in range(len(entry))])
        if mvfenc == 'invar':
            return entry
        if mvfenc == 'onecov':
            num = int(entry[3:])
            if optargs[0][0] == 0 and num in optargs[0]:
                return ''
            return "{}{}".format(entry[0:3], num - len([
                x for x in optargs[0] if x < num]))
        if mvfenc == 'onevar':
            num = int(entry[4:])
            if optargs[0][0] == 0:
                if num in optargs[0]:
                    entry = entry[0:1]
            elif optargs[0][0] == num:
                entry = "{}{}".format(entry[0:4], num - len([
                    x for x in optargs[0] if x < num]))
            elif num in optargs[0]:
                entry = entry[0:2]
            return entry
        if mvfenc == 'refvar':
            if 0 in optargs[0] and optargs[0][0] != 0:
                return entry[1]
            if optargs[0][0] == 0 and len(optargs[0]) == ncol:
                return entry[0]
        return entry

    # COLUMNS
    def columns(entry, mvfenc):
        """return only these sample columns (arg=1,2,3...)"""
        if mvfenc == 'full':
            return ''.join([entry[j] for j in optargs[0]])
        if mvfenc == 'invar':
            return entry
        if mvfenc == 'onecov':
            num = int(entry[3:])
            if num not in optargs[0]:
                if 0 in optargs[0]:
                    return '{}-'.format(entry[0])
                return "-"
            if list(sorted(optargs[0])) == [0, num]:
                return '{}{}'.format(entry[0], entry[2])
            return '{}{}'.format(entry[0:3],
                                 len([x < num for x in optargs[0]]))
        if mvfenc == 'onevar':
            if optargs[0] == [0]:
                return entry[0]
            num = int(entry[4:])
            print(optargs, num, entry)
            if num not in optargs[0]:
                if 0 not in optargs[0]:
                    return entry[1]
                return entry[0:2]
            if list(sorted(optargs[0])) == [0, num]:
                return '{}{}'.format(entry[0], entry[3])
            return '{}{}'.format(entry[0:4],
                                 len([x < num for x in optargs[0]]))
        if mvfenc == 'refvar':
            if optargs[0] == [0]:
                return entry[0]
            return entry if 0 in optargs[0] else entry[1]
        return entry

    # MASKCHAR
    def maskchar(entry, mvfenc):
        """replace specified characters with 'X'"""
        if mvfenc in ('refvar', 'full'):
            return ''.join("X" if x in optargs[0] else x for x in entry)
        if mvfenc == 'invar':
            return 'X' if entry in optargs[0] else entry
        if mvfenc == 'onecov':
            return "{}+{}{}".format(
                entry[0] in optargs[0] and 'X' or entry[0],
                entry[2] in optargs[0] and 'X' or entry[2],
                entry[3:])
        if mvfenc == 'onevar':
            if entry[1] in optargs[0] and entry[3] in optargs[0]:
                if entry[0] in optargs[0]:
                    return 'X'
                return '{}X'.format(entry[0])
            return "{}{}+{}{}".format(
                entry[0] in optargs[0] and 'X' or entry[0],
                entry[1] in optargs[0] and 'X' or entry[1],
                entry[3] in optargs[0] and 'X' or entry[3],
                entry[4:])
        return entry

    # MASKLOWER
    def masklower(entry, mvfenc):
        """turn lower case to 'X'"""
        if mvfenc == 'invar':
            return 'X' if entry.islower() else entry
        if mvfenc in ('refvar', 'full'):
            return ''.join(['X' if x.islower() is True
                            else x for x in entry])
        if mvfenc == 'onecov':
            return "{}+{}{}".format(
                'X' if entry[0].islower() else entry[0],
                'X' if entry[2].islower() else entry[2],
                entry[3:])
        if mvfenc == 'onevar':
            return "{}{}{}{}".format(
                'X' if entry[0].islower() else entry[0],
                'X' if entry[1].islower() else entry[1],
                'X' if entry[3].islower() else entry[2],
                entry[3:])
        return entry

    # MINCOVERAGE
    def mincoverage(entry, mvfenc):
        """minimum sample coverage"""
        # print(mvfenc)
        if mvfenc == 'full':
            zcnt = 0
            for allele in entry:
                if allele not in 'X-':
                    zcnt += 1
                    if zcnt == optargs[0][0]:
                        return True
            return False
        if mvfenc == 'invar':
            return entry not in 'X-'
        if mvfenc == 'onecov':
            return optargs[0][0] <= 2
        if mvfenc == 'onevar':
            if entry[1] == '-' and optargs[0][0] >= 2:
                return False
            if entry[3] == '-' and optargs[0][0] > ncol - 1:
                return False
            return True
        if mvfenc == 'refvar':
            if entry[0] in 'X-':
                if entry[1] in 'X-':
                    return False
                if optargs[0][0] >= ncol - 1:
                    return False
            if entry[1] in 'X-' and optargs[0][0] > 1:
                return False
            return True
        return False

    # NOTCHAR
    def notchar(entry, mvfenc):
        """filter out if any specified character present in entry
        """
        if mvfenc == 'full':
            return all([x not in optargs[0] for x in entry])
        if mvfenc in ('invar', 'refvar'):
            return all([x not in optargs[0] for x in entry])
        if mvfenc == 'onecov':
            return all([entry[x] not in optargs[0] for x in (0, 2)])
        if mvfenc == 'onevar':
            return all([entry[x] not in optargs[0] for x in (0, 1, 3)])
        return False

    # PROMOTELOWER
    def promotelower(entry, mvfenc):
        """turn lower case to upper case
        """
        if mvfenc in ('full', 'refvar', 'onecov', 'invar'):
            return entry.upper()
        if mvfenc == 'onevar':
            if entry[1].upper() == entry[3].upper():
                return entry[0:2].upper()
            return entry.upper()
        return entry.upper()

    # RANDOMALLELE
    def randomallele(entry, mvfenc):
        """replaces ambiguous alleles (RYMKWSBDHV) with a randomly selected
           allele"""
        # if all(x not in MLIB.validchars['dnaambig23'] for x in entry):
        #     return entry
        if mvfenc == 'full':
            return ''.join(
                (
                    MLIB.randomnuc(x)
                    if x in MLIB.validchars['dnaambig23']
                    else x
                )
                    for x in entry
            )
        if mvfenc == 'refvar':
            return '{}{}'.format(
                     (MLIB.randomnuc(entry[0])
                      if entry[0] in MLIB.validchars['dnaambig23']
                      else entry[0]),
                     (''.join(MLIB.randomnuc(entry[1])
                              for _ in range(ncol - 1))
                      if entry[1] in MLIB.validchars['dnaambig23']
                      else entry[1])
            )
        if mvfenc == 'invar':
            return (''.join(MLIB.randomnuc(entry) for _ in range(ncol))
                    if entry in MLIB.validchars['dnaambig23']
                    else entry)
        if mvfenc == 'onecov':
            return "{}+{}{}".format(
                (MLIB.randomnuc(entry[0])
                 if entry[0] in MLIB.validchars['dnaambig23']
                 else entry[0]),
                (MLIB.randomnuc(entry[2])
                 if entry[2] in MLIB.validchars['dnaambig23']
                 else entry[2]),
                entry[3:])
        if mvfenc == 'onevar':
            varloc = int(entry[4:])
            return "{}{}".format(
                (MLIB.randomnuc(entry[0])
                 if entry[0] in MLIB.validchars['dnaambig23']
                 else entry[0]),
                ''.join((MLIB.randomnuc(entry[3])
                         if entry[3] in MLIB.validchars['dnaambig23']
                         else entry[3])
                        if j == varloc else (
                            MLIB.randomnuc(entry[1])
                            if entry[1] in MLIB.validchars['dnaambig23']
                            else entry[1])
                        for j in range(ncol - 1)))
        return entry

    # REMOVECHAR
    def removechar(entry, mvfenc):
        """replace specified characters with '-'"""
        if mvfenc in ('refvar', 'full'):
            if all([(x in optargs[0] or x == '-') for x in entry]):
                return ''
            return ''.join(['-' if x in optargs[0] else x for x in entry])
        if mvfenc == 'invar':
            return entry if entry not in optargs[0] else ''
        if mvfenc == 'onecov':
            if entry[2] in optargs[0] or entry[2] == '-':
                if entry[0] in optargs[0] or entry[0] == '-':
                    return ''
                return "{}-".format(entry[0])
            if entry[0] in optargs[0] or entry[0] == '-':
                return '-{}'.format(entry[1:])
        elif mvfenc == 'onevar':
            if all([(entry[x] in optargs[0] or
                     entry[x] == '-') for x in (1, 3)]):
                if entry[0] in optargs[0] or entry[0] == '-':
                    return ''
                return "{}-".format(entry[0])
            return '{}{}+{}{}'.format(
                ('-' if (entry[0] in optargs[0] or entry[0] == '-')
                 else entry[0]),
                ((entry[1] not in optargs[0] and entry[1] != '-') and
                 entry[1] or ''),
                ('-' if (entry[3] in optargs[0] or entry[3] == '-') else
                 entry[3]),
                entry[4:])
        return entry

    # REMOVELOWER
    def removelower(entry, mvfenc):
        """turn lower case to '-'
        """
        if mvfenc in ('refvar', 'full'):
            entry = ''.join([x if (x.isupper() or x == '-')
                             else '-' for x in entry])
            if all(x == '-' for x in entry):
                return ''
        elif mvfenc == 'invar':
            return entry if entry.isupper() else ''
        elif mvfenc == 'onecov':
            if entry[2].islower() or entry[2] == '-':
                if entry[0].islower() or entry[0] == '-':
                    return ''
                return "{}-".format(entry[0])
            if entry[0].islower() or entry[0] == '-':
                return "-{}".format(entry[1:])
        elif mvfenc == 'onevar':
            if all([(entry[x].islower() or entry[x] == '-')
                    for x in (1, 3)]):
                if entry[0].islower() or entry[0] == '-':
                    return ''
                return "{}-".format(entry[0])
            return '{}{}+{}{}'.format(
                ('-' if (entry[0].islower() or entry[0] == '-')
                 else entry[0]),
                (entry[1] if (entry[1].isupper() or entry[1] == '-')
                 else ''),
                ('-' if (entry[3].islower() or entry[3] == '-')
                 else entry[3]),
                entry[4:])
        return entry

    # REQALLCHAR
    def reqallchar(entry, mvfenc):
        """require all of the specified characters appear in the entry
        """
        if mvfenc in ('full', 'invar', 'refvar'):
            return all([x in entry for x in optargs[0]])
        if mvfenc == 'onecov':
            return all([x in (entry[0], entry[2]) for x in optargs[0]])
        if mvfenc == 'onevar':
            return all([x in (entry[0], entry[1], entry[3])
                        for x in optargs[0]])
        return False

    # REQCONTIG
    def reqcontig(entry):
        """return sites in ID,START,STOP (inclusive)
        """
        return entry[0] in optargs[0]

    # REQINFORMATIVE
    def reqinformative(entry, mvfenc):
        """only retain parsimony informative sites (2+ alleles in 2+ samples)
        """
        if mvfenc == 'full':
            return len([x for x in set(entry.upper())
                        if entry.upper().count(x) > 1 and
                        x not in 'X-']) > 1
        return False

    # REQINVARIANT
    def reqinvariant(entry, mvfenc):
        """only retain invariant sites"""
        if mvfenc in ('full', 'onevar', 'refvar'):
            return False
        if mvfenc == 'invar':
            return True
        if mvfenc == 'onecov':
            return entry[0].upper() == entry[2].upper()
        return False

    # REQNONREFSAMPLE
    def reqnonrefsample(entry, mvfenc):
        """Returns entries only where one non-reference sample
           has an allele (not X or -)
        """
        if mvfenc == 'full':
            for i in range(1, len(entry)):
                if entry[i] not in 'X-':
                    return True
        elif mvfenc == 'invar':
            if entry not in 'X-':
                return True
        elif mvfenc == 'onecov':
            if entry[2] not in 'X-':
                return True
        elif mvfenc == 'onevar':
            if entry[1] not in 'X-' or entry[3] not in 'X-':
                return True
        elif mvfenc == 'refvar':
            if entry[1] not in 'X-':
                return True
        return False

    # REQONECHAR
    def reqonechar(entry, mvfenc):
        """require one of the specified characters appear in entry
        """
        if mvfenc in ('full', 'invar', 'refvar'):
            return any([x in entry for x in optargs[0]])
        if mvfenc == 'onecov':
            return any([x in (entry[0], entry[2]) for x in optargs[0]])
        if mvfenc == 'onevar':
            return any([x in (entry[0], entry[1], entry[3])
                        for x in optargs[0]])
        return False

    # REQREGION
    def reqregion(entry):
        """return sites in ID,START,STOP (inclusive)"""
        return (entry[0] == optargs[0][0] and
                optargs[0][1] <= entry[1] <= optargs[0][2])

    # REQSAMPLE
    def reqsample(entry, mvfenc):
        """require specific samples to have alleles (not X/-)
        """
        if mvfenc == 'full':
            return all([entry[x] not in 'X-' for x in optargs[0]])
        if mvfenc == 'invar':
            return entry not in 'X-'
        if mvfenc == 'onecov':
            return (all([x in [0, int(entry[3:])] for x in optargs[0]])
                    and (0 in optargs[0]
                         and entry[0] not in 'X-')
                    and (int(entry[3:]) in optargs[0]
                         and entry[2] not in 'X'))
        if mvfenc == 'onevar':
            return not((entry[0] in 'X-'
                        and 0 in optargs[0])
                       or (entry[3] in 'X-'
                           and (int(entry[4:]) in optargs[0]))
                       or (entry[1] in 'X-'
                           and (set(optargs[0]) - set([0, int(entry[4:])])))
                       )
        if mvfenc == 'refvar':
            return not((entry[0] in 'X-' and 0 in optargs[0])
                       or (entry[1] in 'X-'
                           and set(optargs[0]) - {0, })
                       )

        return False

    # REQVARIANT
    def reqvariant(entry, mvfenc):
        """only retain variable sites
        """
        if mvfenc == 'full':
            return len(set(entry.upper()) - set('X-')) > 1
        if mvfenc == 'invar':
            return False
        if mvfenc == 'onecov':
            return (entry[0].upper() != entry[2].upper() if
                    entry[0] not in 'X-' and entry[2] not in 'X-' else
                    False)
        if mvfenc == 'onevar':
            return (len(set([entry[0], entry[1], entry[2]])) > 1 if
                    entry[2] not in 'X-' and entry[1] not in 'X-' and
                    entry[0] not in 'X-' else False)
        if mvfenc == 'refvar':
            return entry[1] not in 'X-'
        return False

    module_toc = {
        'allelegroup': ('filter', allelegroup),
        'collapsemerge': ('transform', collapsemerge),
        'collapsepriority': ('transform', collapsepriority),
        'columns': ('transform', columns),
        'maskchar': ('transform', maskchar),
        'masklower': ('transform', masklower),
        'mincoverage': ('filter', mincoverage),
        'notchar': ('filter', notchar),
        'promotelower': ('transform', promotelower),
        'randomallele': ('transform', randomallele),
        'removechar': ('transform', removechar),
        'removelower': ('transform', removelower),
        'reqallchar': ('filter', reqallchar),
        'reqcontig': ('location', reqcontig),
        'reqinformative': ('filter', reqinformative),
        'reqinvariant': ('filter', reqinvariant),
        'reqnonrefsample': ('filter', reqnonrefsample),
        'reqonechar': ('filter', reqonechar),
        'reqregion': ('location', reqregion),
        'reqsample': ('filter', reqsample),
        'reqvariant': ('filter', reqvariant),
        }

    if modulename == 'getnames':
        return tuple(sorted(x for x in module_toc))

    return (modulename, module_toc[modulename][0],
            module_toc[modulename][1], optargs)

# END OF MODULE DEFINITION

# MODULENAMES = ['allelegroup', 'collapsemerge', 'collapsepriority',
#               'columns', 'maskchar', 'masklower', 'mincoverage',
#               'notchar', 'promotelower', 'randomallele',
#               'removechar', 'removelower',
#               'reqallchar', 'reqcontig', 'reqinformative',
#               'reqinvariant', 'reqonechar', 'reqregion',
#               'reqsample', 'reqvariant', 'reqnonrefsample']


MODULENAMES = make_module('getnames', 0)


def modulehelp():
    """Extended Help for Modules
    """
    for modulename in sorted(MODULENAMES):
        modulename, _, module, _ = make_module(
            modulename, 0, optargs='')
        description = module.__doc__.replace("\n", " ")
        while "  " in description:
            description = description.replace("  ", " ")
        print(modulename, description.strip())
    return ''

# HELP Generator


def build_actionset(moduleargs, ncol):
    """Create action set modules
        Arguments:
            moduleargs: arguments for using the module
            ncol: int number of columns in the base MVF
    """
    actionset = []
    if moduleargs is None:
        raise RuntimeError("ERROR: No modules selected")
    for module in moduleargs:
        if ':' in module:
            modargs = module.split(':')
            if modargs[0] not in MODULENAMES:
                raise RuntimeError("Module {} not found".format(modargs[0]))
            for i in range(1, len(modargs)):
                modargs[i] = (modargs[i].split(',') if ',' in modargs[i] else
                              [modargs[i]])
            if modargs[0] == 'reqregion':
                for i in range(1, len(modargs)):
                    if len(modargs[i]) != 3:
                        raise RuntimeError((
                            "ERROR: Region specification '{}' "
                            "invalid!").format(modargs[i]))
                    modargs[i] = (modargs[i][0], int(modargs[i][1]),
                                  int(modargs[i][2]))
            elif modargs[0] == 'mincoverage':
                modargs[1][0] = int(modargs[1][0])
                if modargs[1][0] > ncol:
                    raise RuntimeError((
                        "ERROR: Minimum columns specified ({}) is "
                        "greater than number of MVF total columns"
                        "({}).").format(modargs[1][0], ncol))
            elif modargs[0] in ['columns', 'collapsemerge', 'reqsample']:
                for i in range(1, len(modargs)):
                    try:
                        modargs[i] = [int(x) for x in modargs[i]]
                    except ValueError:
                        continue
            actionset.append(make_module(modargs[0],
                                         ncol, optargs=modargs[1:]))
        else:
            actionset.append(make_module(module, ncol))
            if module not in MODULENAMES:
                raise RuntimeError("Module {} not found".format(module))
    return actionset


def filter_mvf(args):
    """Main method"""
    args.qprint("Running FilterMVF")
    if args.more_help is True:
        modulehelp()
        sys.exit()
    if args.mvf is None and args.test is None:
        raise RuntimeError("No input file specified with --mvf")
    if args.out is None and args.test is None:
        raise RuntimeError("No output file specified with --out")
    # Establish Input MVF
    if args.test is not None:
        ncol = args.test_nchar or len(args.test.split()[1])
    else:
        mvf = MultiVariantFile(args.mvf, 'read')
        ncol = mvf.metadata['ncol']
    args.qprint("Input MVF read with {} columns.".format(ncol))
    # Create Actionset
    if args.labels:
        for i in range(len(args.actions)):
            action = args.actions[i]
            arr = action.split(':')
            if arr[0] in ('collapsepriority', 'collapsemerge'):
                arr[1] = ','.join([
                    str(mvf.sample_id_to_index[x])
                    for x in arr[1].split(',')])
            if arr[0] in ('columns', 'allelegroup', 
                          'notmultigroup', 'reqsample'):
                for j in range(1, len(arr)):
                    arr[j] = ','.join([
                        str(mvf.sample_id_to_index[x])
                        for x in arr[j].split(',')])
            args.actions[i] = ':'.join(arr)
    removed_columns = set([])
    for i in range(len(args.actions)):
        action = args.actions[i]
        arr = action.split(':')
        if arr[0] in ('collapsepriority', 'collapsemerge'):
            tmp_arr = arr[1][:]
            arr[1] = ','.join([
                str(int(x) - len([y for y in removed_columns if y < int(x)]))
                for x in arr[1].split(',')])
            removed_columns.update([int(x) for x in tmp_arr.split(',')[1:]])
            print(arr)
            print(removed_columns)
        if arr[0] in ('columns', 'allelegroup', 
                      'notmultigroup', 'reqsample'):
            for j in range(1, len(arr)):
                arr[j] = ','.join([
                    str(int(x) - len([y for y in removed_columns if y < int(x)]))
                    for x in arr[j].split(',')])
        args.actions[i] = ':'.join(arr)
            
            
    actionset = build_actionset(args.actions, ncol)
    args.qprint("Actions established.")
    args.qprint(actionset)
    # TESTING MODE
    if args.test:
        loc, alleles = args.test.split()
        linefail = False
        transformed = False
        # invar = invariant (single character)
        # refvar (all different than reference, two chars)
        # onecov (single coverage, + is second character)
        # onevar (one variable base, + is third character)
        # full = full alleles (all chars)
        if args.verbose:
            print(alleles)
        linetype = get_linetype(alleles)
        sys.stdout.write("MVF Encoding type '{}' detected\n".format(linetype))
        for actionname, actiontype, actionfunc, actionarg in actionset:
            sys.stdout.write("Applying action {} ({}): ".format(
                actionname, actiontype))
            if actiontype == 'filter':
                if not actionfunc(alleles, linetype):
                    linefail = True
                    sys.stdout.write("Filter Fail\n")
                    break
                sys.stdout.write("Filter Pass\n")
            elif actiontype == 'transform':
                transformed = True
                alleles = actionfunc(alleles, linetype)
                linetype = get_linetype(alleles)
                if linetype == 'empty':
                    linefail = True
                    sys.stdout.write("Transform removed all alleles\n")
                    break
                sys.stdout.write("Transform result {}\n".format(alleles))
            elif actiontype == 'location':
                loc = loc.split(':')
                loc[1] = int(loc[1])
                if actionfunc(loc) is False:
                    linefail = True
                    sys.stdout.write("Location Fail\n")
                    break
                sys.stdout.write("Location Pass\n")
        if linefail is False:
            if transformed:
                if linetype == 'full':
                    alleles = encode_mvfstring(alleles)
                if alleles:
                    test_output = "{}\t{}\n".format(loc, alleles)
                    sys.stdout.write("Final output = {}\n".format(
                        test_output))
                else:
                    sys.stdout.write("Transform removed all alleles\n")
            else:
                sys.stdout.write("No changes applied\n")
                sys.stdout.write("Final output = {}\n".format(args.test))
        sys.exit()
    # MAIN MODE
    # Set up file handler
    outmvf = MultiVariantFile(args.out, 'write', overwrite=args.overwrite)
    outmvf.copy_headers_from(mvf)

    removed_indices = set([])
    # reprocess header if actions are used that filter columns
    if any(x == y[0] for x in ('columns', 'collapsepriority', 'collapsemerge')
           for y in actionset):
        for actionname, actiontype, actionfunc, actionarg in actionset:
            if actionname == 'columns':
                if args.labels:
                    oldindices = [outmvf.sample_id_to_index[int(x)]
                                  for x in actionarg[0]]
                else:
                    oldindices = [int(x) for x in actionarg[0]]
            elif actionname in ('collapsepriority', 'collapsemerge'):
                actionarg[0] = [x - len([y for y in removed_indices if y < x])
                                 for x in actionarg[0]]
                oldindices = [x for x in outmvf.sample_indices
                              if x not in actionarg[0][1:]]
            outmvf.sample_ids = outmvf.get_sample_ids(oldindices)
            outmvf.sample_data = dict(
                (i, outmvf.sample_data[oldindices[i]])
                for i, _ in enumerate(oldindices))

            if actionname in ('collapsepriority', 'collapsemerge'):
                if len(actionarg) == 2:
                    outmvf.sample_data[actionarg[0][0]]['id'] = actionarg[1][0]
                    outmvf.sample_ids[actionarg[0][0]] = actionarg[1][0]
            outmvf.sample_indices = list(range(len(oldindices)))
    Youtmvf.metadata['ncol'] = len(outmvf.sample_indices)
    outmvf.notes.append(args.command_string)
    outmvf.write_data(outmvf.get_header())
    args.qprint("Output MVF established.")
    # End header editing
    linebuffer = []
    nbuffer = 0
    args.qprint("Processing Entries.")
    write_total = 0
    for chrom, pos, allelesets in mvf.iterentries(decode=False):
        linefail = False
        transformed = False
        # invar = invariant (single character)
        # refvar (all different than reference, two chars)
        # onecov (single coverage, + is second character)
        # onevar (one variable base, + is third character)
        # full = full alleles (all chars)
        alleles = allelesets[0]
        linetype = get_linetype(alleles)
        if linetype == 'empty':
            continue
        if args.verbose is True:
            sys.stdout.write(" {} {} ".format(alleles, linetype))
        for actionname, actiontype, actionfunc, _ in actionset:
            if actiontype == 'filter':
                linefail = not actionfunc(alleles, linetype)
            elif actiontype == 'transform':
                transformed = True
                alleles = actionfunc(alleles, linetype)
                linetype = get_linetype(alleles)
                linefail = linetype == 'empty'
            elif actiontype == 'location':
                linefail = not actionfunc([chrom, pos])
            if linefail:
                break
        if linefail is False:
            if transformed:
                if linetype == 'full':
                    alleles = mvf.encode(alleles)
                if not alleles:
                    linefail = True
            nbuffer += 1
            linebuffer.append((chrom, pos, (alleles,)))
            if args.verbose:
                sys.stdout.write("{}\n".format(alleles))
            if nbuffer == args.line_buffer:
                write_total += args.line_buffer
                args.qprint("{} entries written. Total written: {}.".format(
                    args.line_buffer, write_total))
                outmvf.write_entries(linebuffer)
                linebuffer = []
                nbuffer = 0
        elif args.verbose:
            sys.stdout.write("FAIL\n")
    if linebuffer:
        outmvf.write_entries(linebuffer)
        write_total += len(linebuffer)
        args.qprint("{} entries written. Total written: {}.".format(
            args.line_buffer, write_total))
        linebuffer = []
    return ''
