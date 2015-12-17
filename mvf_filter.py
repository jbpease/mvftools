#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
MVFtools: Multisample Variant Format Toolkit
http://www.github.org/jbpease/mvftools

mvf_filter_dev: Testing for new Modules for Filtering/Transformation

@author: James B. Pease
@author: Ben K. Rosenzweig

version: 2015-02-16 Release
version: 2015-06-09 Fixes and updates for 1.2.1
@version: 2015-09-04 Style fixes and upgrades

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


from __future__ import print_function
import sys
import argparse
from copy import deepcopy
from itertools import combinations
from time import time
from mvfbase import MultiVariantFile, encode_mvfstring


def get_linetype(alleles):
    """Determines the line type from allele string
    """
    if not alleles:
        return 'empty'
    if len(alleles) == 1:
        linetype = 'invar'
        if alleles[0] == '-':
            linetype = 'empty'
    elif len(alleles) == 2:
        linetype = 'refvar'
        if alleles[0:1] == '--':
            linetype = 'empty'
    elif alleles[1] == '+':
        linetype = 'onecov'
        if alleles[0] == '-' and alleles[2] == '-':
            linetype = 'empty'
    elif alleles[2] == '+':
        linetype = 'onevar'
        if all([alleles[x] == '-' for x in (0, 1, 3)]):
            linetype = 'empty'
    else:
        linetype = 'full'
        if all([x == '-' for x in alleles]):
            return 'empty'
    return linetype


def make_module(modulename, ncol, optargs=None):
    """Generate Modules"""

    # ALLELEGROUP
    if modulename == 'allelegroup':
        moduletype = 'filter'

        def allelegroup(entry, mvfenc):
            """Only retain sites where all members of the specific
               group contain alleles
            """
            if mvfenc == 'full':
                allele_groups = [set([entry[x] for x in y]) -
                                 set('X-') for y in optargs]
                if not all(allele_groups):
                    return False
                if any(grp0 & grp1 for grp0, grp1 in
                       combinations(allele_groups, 2)):
                    return False
                return True
            else:
                return False
    # COLLAPSEPRIORITY
    elif modulename == "collapsepriority":
        moduletype = "transform"

        def collapsepriority(entry, mvfenc):
            """sample alleles combined,
               using a priority list if gap encountered"""
            if mvfenc == 'full':
                colbases = [entry[x] for x in optargs[0]
                            if entry[x] not in 'NX-'] + ['-']
                return ''.join([(j not in optargs[0] and entry[j]) or
                                (j == optargs[0][0] and colbases[0]) or
                                '' for j in range(len(entry))])
            elif mvfenc == 'invar':
                return entry
            elif mvfenc == 'onecov':
                num = int(entry[3:])
                if optargs[0][0] == 0 and num in optargs[0]:
                    return ''
                return "{}{}".format(entry[0:3], num - len([
                    x for x in optargs if x < num]))
            elif mvfenc == 'onevar':
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
            elif mvfenc == 'refvar':
                if 0 in optargs[0] and optargs[0][0] is not 0:
                    return entry[1]
                elif optargs[0][0] == 0 and len(optargs[0]) == ncol:
                    return entry[0]
                return entry
    # COLUMNS
    elif modulename == 'columns':
        moduletype = 'transform'

        def columns(entry, mvfenc):
            """return only these sample columns (arg=1,2,3...)"""
            if mvfenc == 'full':
                return ''.join([entry[j] for j in optargs[0]])
            elif mvfenc == 'invar':
                return entry
            elif mvfenc == 'onecov':
                num = int(entry[3:])
                if num not in optargs[0]:
                    return '{}-'.format(entry[0])
                elif list(sorted(optargs[0])) == [0, num]:
                    return '{}{}'.format(entry[0], entry[2])
                return '{}{}'.format(entry[0:3],
                                     len([x < num for x in optargs[0]]))
            elif mvfenc == 'onevar':
                if optargs[0] == [0]:
                    return entry[0]
                num = int(entry[4:])
                if num not in optargs[0]:
                    return entry[0:2]
                elif list(sorted(optargs[0])) == [0, num]:
                    return '{}{}'.format(entry[0], entry[3])
                return '{}{}'.format(entry[0:4],
                                     len([x < num for x in optargs[0]]))
            elif mvfenc == 'refvar':
                if optargs[0] == [0]:
                    return entry[0]
                else:
                    return 0 in optargs[0] and entry or entry[1]

    # MASKCHAR
    elif modulename == "maskchar":
        moduletype = "transform"

        def maskchar(entry, mvfenc):
            """replace specified characters with 'X'"""
            if mvfenc in ('refvar', 'full'):
                return ''.join([x in optargs[0] and 'X' or x for x in entry])
            elif mvfenc == 'invar':
                return entry in optargs[0] and 'X' or entry
            elif mvfenc == 'onecov':
                return "{}+{}{}".format(
                    entry[0] in optargs[0] and 'X' or entry[0],
                    entry[2] in optargs[0] and 'X' or entry[2],
                    entry[3:])
            elif mvfenc == 'onevar':
                if entry[1] in optargs[0] and entry[3] in optargs[0]:
                    if entry[0] in optargs[0]:
                        return 'X'
                    else:
                        return '{}X'.format(entry[0])
                return "{}{}+{}{}".format(
                    entry[0] in optargs[0] and 'X' or entry[0],
                    entry[1] in optargs[0] and 'X' or entry[1],
                    entry[3] in optargs[0] and 'X' or entry[3],
                    entry[4:])

    # MASKLOWER
    elif modulename == 'masklower':
        moduletype = 'transform'

        def masklower(entry, mvfenc):
            """turn lower case to 'X'"""
            if mvfenc == 'invar':
                return entry.islower() and 'X' or entry
            elif mvfenc in ('refvar', 'full'):
                return ''.join([x.isupper() and 'X' or x for x in entry])
            elif mvfenc == 'onecov':
                return "{}+{}{}".format(
                    entry[0].islower() and 'X' or entry[0],
                    entry[2].islower() and 'X' or entry[2],
                    entry[3:])
            elif mvfenc == 'onevar':
                return "{}{}+{}{}".format(
                    entry[0].islower() and 'X' or entry[0],
                    entry[1].islower() and 'X' or entry[1],
                    entry[3].islower() and 'X' or entry[2],
                    entry[4:])

    # MINCOVERAGE
    elif modulename == "mincoverage":
        moduletype = 'filter'

        def mincoverage(entry, mvfenc):
            """minimum sample coverage"""
            if mvfenc == 'full':
                return sum(
                    [int(x not in 'NX-') for x in entry]) >= optargs[0][0]
            if mvfenc in ('invar', 'refvar'):
                return True
            elif mvfenc == 'onecov':
                return optargs[0][0] <= 2
            elif mvfenc == 'onevar':
                if entry[1] == '-' and optargs[0][0] >= 2:
                    return False
                elif entry[3] == '-' and optargs[0][0] > ncol - 1:
                    return False
                return True

    # NOTCHAR
    elif modulename == "notchar":
        moduletype = 'filter'

        def notchar(entry, mvfenc):
            """filter out if any specified character present in entry"""
            if mvfenc == 'full':
                return all([x not in optargs[0] for x in entry])
            elif mvfenc in ('invar', 'refvar'):
                return all([x not in optargs[0] for x in entry])
            elif mvfenc == 'onecov':
                return all([entry[x] not in optargs[0] for x in (0, 2)])
            elif mvfenc == 'onevar':
                return all([entry[x] not in optargs[0] for x in (0, 1, 3)])

    # PROMOTELOWER
    elif modulename == 'promotelower':
        moduletype = 'transform'

        def promotelower(entry, mvfenc):
            """turn lower case to upper case"""
            if mvfenc in ['full', 'refvar', 'onecov', 'invar']:
                return entry.upper()
            elif mvfenc == 'onevar':
                if entry[1].upper() == entry[3].upper():
                    return entry[0:2].upper()
                else:
                    return entry.upper()

    # REMOVELOWER
    elif modulename == 'removelower':
        moduletype = 'transform'

        def removelower(entry, mvfenc):
            """turn lower case to '-'"""
            if mvfenc in ('refvar', 'full'):
                entry = ''.join([(x.isupper() or x == '-') and
                                 x or '-' for x in entry])
                if all(x == '-' for x in entry):
                    return ''
                return entry
            elif mvfenc == 'invar':
                return entry.isupper() and entry or ''
            elif mvfenc == 'onecov':
                if entry[2].islower() or entry[2] == '-':
                    if entry[0].islower() or entry[0] == '-':
                        return ''
                    else:
                        return "{}-".format(entry[0])
                if entry[0].islower() or entry[0] == '-':
                    return "-{}".format(entry[1:])
                return entry
            elif mvfenc == 'onevar':
                if all([(entry[x].islower() or entry[x] == '-')
                        for x in (1, 3)]):
                    if entry[0].islower() or entry[0] == '-':
                        return ''
                    else:
                        return "{}-".format(entry[0])
                return '{}{}+{}{}'.format(
                    ((entry[0].islower() or entry[0] == '-') and
                     '-' or entry[0]),
                    ((entry[1].isupper() or entry[1] != '-') and
                     entry[1] or ''),
                    ((entry[3].islower() or entry[3] == '-') and
                     '-' or entry[3]),
                    entry[4:])

    # REMOVECHAR
    elif modulename == "removechar":
        moduletype = "transform"

        def removechar(entry, mvfenc):
            """"replace specified characters with '-'"""
            if mvfenc in ('refvar', 'full'):
                if all([x in optargs[0] or x == '-' for x in entry]):
                    return ''
                return ''.join([x in optargs[0] and '-' or x for x in entry])
            elif mvfenc == 'invar':
                return entry not in optargs[0] and entry or ''
            elif mvfenc == 'onecov':
                if entry[2] in optargs[0] or entry[2] == '-':
                    if entry[0] in optargs[0] or entry[0] == '-':
                        return ''
                    else:
                        return "{}-".format(entry[0])
                if entry[0] in optargs[0] or entry[0] == '-':
                    return '-{}'.format(entry[1:])
                return entry
            elif mvfenc == 'onevar':
                if all([(entry[x] in optargs[0] or
                         entry[x] == '-') for x in (1, 3)]):
                    if entry[0] in optargs[0] or entry[0] == '-':
                        return ''
                    else:
                        return "{}-".format(entry[0])
                return '{}{}+{}{}'.format(
                    ((entry[0] in optargs[0] or entry[0] == '-') and '-' or
                     entry[0]),
                    ((entry[1] not in optargs[0] and entry[1] != '-') and
                     entry[1] or ''),
                    ((entry[3] in optargs[0] or entry[3] == '-') and
                     '-' or entry[3]),
                    entry[4:])

    # REQALLCHAR
    elif modulename == "reqallchar":
        moduletype = "filter"

        def reqallchar(entry, mvfenc):
            """require all of the specified characters appear in the entry"""
            if mvfenc in ('full', 'invar', 'refvar'):
                return all([x in entry for x in optargs[0]])
            elif mvfenc == 'onecov':
                return all([x in (entry[0], entry[2]) for x in optargs[0]])
            elif mvfenc == 'onevar':
                return all([x in (entry[0], entry[1], entry[3])
                            for x in optargs[0]])

    # REQCONTIG
    elif modulename == 'reqcontig':
        moduletype = 'location'

        def reqcontig(entry):
            """return sites in ID,START,STOP (inclusive)"""
            return entry[0] in optargs[0]

    # REQINFORMATIVE
    elif modulename == 'reqinformative':
        moduletype = 'filter'

        def reqinformative(entry, mvfenc):
            """only retain informative sites (2+ alleles in 2+ samples)"""
            if mvfenc == 'full':
                return len([x for x in set(entry.upper())
                            if entry.upper().count(x) > 1 and
                            x not in 'NX-']) > 1
            return False

    # REQINVARIANT
    elif modulename == 'reqinvariant':
        moduletype = 'filter'

        def reqinvariant(entry, mvfenc):
            """only retain invariant sites"""
            if mvfenc in ('full', 'onevar', 'refvar'):
                return False
            elif mvfenc == 'invar':
                return True
            elif mvfenc == 'onecov':
                return entry[0].upper() == entry[2].upper()

    # REQREGION
    elif modulename == "reqregion":
        moduletype = 'location'

        def reqregion(entry):
            """return sites in ID,START,STOP (inclusive)"""
            return (entry[0] == optargs[0][0] and
                    optargs[0][1] <= entry[1] <= optargs[0][2])

    # REQONECHAR
    elif modulename == "reqonechar":
        moduletype = 'filter'

        def reqonechar(entry, mvfenc):
            """require one of the specified characters appear in entry
            """
            if mvfenc in ('full', 'invar', 'refvar'):
                return any([x in entry for x in optargs[0]])
            elif mvfenc == 'onecov':
                return any([x in (entry[0], entry[2]) for x in optargs[0]])
            elif mvfenc == 'onevar':
                return any([x in (entry[0], entry[1], entry[3])
                            for x in optargs[0]])

    # REQSAMPLE
    elif modulename == "reqsample":
        moduletype = 'filter'

        def reqsample(entry, mvfenc):
            """require specific samples to be present
            """
            if mvfenc == 'full':
                return all([entry[x] not in 'NX-' for x in optargs[0]])
            elif mvfenc == 'invar':
                return True
            elif mvfenc == 'onecov':
                return all([x in [0, int(entry[3:])] for x in optargs[0]])
            elif mvfenc == 'onevar':
                return not (entry[3] == '-' and int(entry[4:]) in optargs[0])

    # REQVARIANT
    elif modulename == 'reqvariant':
        moduletype = 'filter'

        def reqvariant(entry, mvfenc):
            """only retain variable sites
            """
            if mvfenc == 'full':
                return len(set(entry.upper()) - set('NX-')) > 1
            elif mvfenc == 'invar':
                return False
            elif mvfenc == 'onecov':
                return (entry[0].upper() != entry[2].upper() and
                        entry[0] not in 'NX-' and entry[2] not in 'NX-')
            elif mvfenc == 'onevar':
                return (entry[0].upper() == entry[1].upper() and
                        entry[2] in 'NX-')
            elif mvfenc == 'refvar':
                return True

    elif modulename == 'reqnonrefsample':
        moduletype = "filter"

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

    return (modulename, moduletype, eval(modulename), optargs)

# END OF MODULE DEFINITION

MODULENAMES = ['allelegroup', 'collapsepriority',
               'columns', 'maskchar', 'masklower', 'mincoverage',
               'notchar', 'promotelower', 'removechar', 'removelower',
               'reqallchar', 'reqcontig', 'reqinformative',
               'reqinvariant', 'reqonechar', 'reqregion',
               'reqsample', 'reqvariant', 'reqnonrefsample']


def modulehelp():
    """Extended Help for Modules
    """
    for modulename in sorted(MODULENAMES):
        modulename, modtype, module, _ = make_module(
            modulename, 0, optargs='')
        print("{}: {} ({})".format(modulename, module.__doc__, modtype))
    sys.exit()
    return ''

# HELP Generator


def build_actionset(moduleargs, ncol):
    """Create action set modules
        Arguments:
            moduleargs: arguments for using the module
            ncol: int number of columns in the base MVF
    """
    actionset = []
    for module in moduleargs:
        if ':' in module:
            modargs = module.split(':')
            if modargs[0] not in MODULENAMES:
                raise RuntimeError("Module {} not found".format(modargs[0]))
            for i in range(1, len(modargs)):
                modargs[i] = (',' in modargs[i] and modargs[i].split(',') or
                              [modargs[i]])
            print(modargs)
            if modargs[0] == 'region':
                modargs[0][1] = int(modargs[0][1])
                modargs[0][2] = int(modargs[0][2])
            elif modargs[0] == 'mincoverage':
                if int(modargs[0][0]) > ncol:
                    raise RuntimeError()
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


def main(arguments=sys.argv[1:]):
    """Main method for mvf_filter"""
    parser = argparse.ArgumentParser(description="""
        Filters and Transforms MVF files""")
    parser.add_argument("--mvf", help="input MVF file")
    parser.add_argument("--out", help="output MVF file")
    parser.add_argument("--actions", nargs='*',
                        help=("set of actions:args to perform,"
                              " note these are done in order as listed"))
    parser.add_argument("--labels", action="store_true",
                        help="use sample labels instead of indices")
    parser.add_argument("--test", help="manually input a line for testing")
    parser.add_argument("--testnchar", type=int,
                        help="total number of samples for test string")
    parser.add_argument("--modulehelp", action="store_true",
                        help="prints full module list and descriptions")
    parser.add_argument("--linebuffer", type=int, default=100000,
                        help="number of lines to write at once to MVF")
    parser.add_argument("--verbose", action="store_true",
                        help="report every line (for debugging)")
    parser.add_argument("--overwrite", action="store_true",
                        help="USE WITH CAUTION: force overwrite of outputs")
    parser.add_argument("--quiet", action="store_true",
                        help="suppress progress meter")
    parser.add_argument("-v", "--version", action="store_true",
                        help="display version information")
    args = parser.parse_args(args=arguments)
    if args.version:
        print("Version 2015-09-04")
        sys.exit()
    args = parser.parse_args(args=arguments)
    time0 = time()
    if args.modulehelp:
        modulehelp()
    if not args.mvf and not args.test:
        raise RuntimeError("No input file specified with --mvf")
    if not args.out and not args.test:
        raise RuntimeError("No output file specified with --outs")
    if not args.actions:
        raise RuntimeError("No --actions specified!")
    # Establish Input MVF
    if args.test:
        ncol = args.testnchar or len(args.test)
    else:
        mvf = MultiVariantFile(args.mvf, 'read')
        ncol = mvf.metadata['ncol']
    # Create Actionset
    if args.labels:
        labels = mvf.get_sample_labels()[:]
        for i in range(len(args.actions)):
            action = args.actions[i]
            arr = action.split(':')
            if arr[0] in ('columns', 'collapsepriority',
                          'allelegroup', 'notmultigroup'):
                for j in range(1, len(arr)):
                    arr[j] = ','.join([
                        str(labels.index(x)) for x in arr[j].split(',')])
            args.actions[i] = ':'.join(arr)
    actionset = build_actionset(args.actions, ncol)
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
                else:
                    sys.stdout.write("Filter Pass\n")
            elif actiontype == 'transform':
                transformed = True
                alleles = actionfunc(alleles, linetype)
                linetype = get_linetype(alleles)
                if linetype == 'empty':
                    linefail = True
                    sys.stdout.write("Transform removed all alleles\n")
                    break
                else:
                    sys.stdout.write("Transform result {}\n".format(alleles))
            elif actiontype == 'location':
                if not actionfunc([int(x) for x in loc.split(':')]):
                    linefail = True
                    sys.stdout.write("Location Fail\n")
                    break
                else:
                    sys.stdout.write("Location Pass\n")
        if not linefail:
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
    outmvf.metadata = deepcopy(mvf.metadata)
    # reprocess header if actions are used that filter columns
    if any(x == y[0] for x in ('columns', 'collapsepriority')
           for y in actionset):
        labels = outmvf.metadata['labels'][:]
        for actionname, actiontype, actionfunc, actionarg in actionset:
            if actionname == 'columns':
                labels = [labels[x] for x in actionarg]
            elif actionname == 'collapsepriority':
                labels = [labels[x] for x in range(len(labels))
                          if x not in actionarg[1:]]
        oldindicies = mvf.get_sample_indices(labels)
        newsamples = {}
        for i, _ in enumerate(labels):
            newsamples[i] = mvf.metadata['samples'][oldindicies[i]]
        outmvf.metadata['samples'] = newsamples.copy()
        outmvf.metadata['labels'] = labels[:]
    outmvf.write_data(outmvf.get_header())
    # End header editing
    linebuffer = []
    nbuffer = 0
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
        if args.verbose:
            sys.stdout.write(" {} {}".format(alleles, linetype))
        for actionname, actiontype, actionfunc, actionargs in actionset:
            if actiontype == 'filter':
                if not actionfunc(alleles, linetype):
                    linefail = True
            elif actiontype == 'transform':
                transformed = True
                alleles = actionfunc(alleles, linetype)
                linetype = get_linetype(alleles)
                if linetype == 'empty':
                    linefail = True
            elif actiontype == 'location':
                if not actionfunc([chrom, pos]):
                    linefail = True
            if linefail:
                break
        if not linefail:
            if transformed:
                if linetype == 'full':
                    alleles = mvf.encode(alleles)
                if not alleles:
                    linefail = True
        if not linefail:
            nbuffer += 1
            linebuffer.append((chrom, pos, (alleles,)))
            if args.verbose:
                sys.stdout.write("{}\n".format(alleles))
            if nbuffer == args.linebuffer:
                outmvf.write_entries(linebuffer)
                linebuffer = []
                nbuffer = 0
        elif args.verbose:
            sys.stdout.write("FAIL\n")
    if linebuffer:
        outmvf.write_entries(linebuffer)
        linebuffer = []
    if not args.quiet:
        print("Completed in {} seconds".format(time() - time0))
    return ''

if __name__ == "__main__":
    main()