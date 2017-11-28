# -*- coding: utf-8 -*-
"""
mvfpaml - PAML library for MVFtools
MVFtools: Multisample Variant Format Toolkit
James B. Pease and Ben K. Rosenzweig
http://www.github.org/jbpease/mvftools

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


import subprocess
import os
import re
from random import randint
from Bio import Phylo


CTLFILE = """seqfile = INPUTPATH
outfile = OUTPUTPATH
noisy = 0
verbose = 0
runmode = -2
seqtype = 1
CodonFreq = 2
aaDist = 0
aaRatefile = dat/wag.dat
model = 0
NSsites = 0
icode = 0
Mgene = 0
fix_kappa = 0
kappa = 2
fix_omega = 0
omega = 0.4
fix_alpha = 1
alpha = 0.
Malpha = 0
ncatG = 8
getSE = 0
RateAncestor = 0
Small_Diff = .5e-6
cleandata = 0
method = 0
fix_blength = 2
"""

BSNULLCTLFILE = """seqfile = INPUTPATH
outfile = OUTPUTPATH
treefile = TREEPATH
noisy = 0
verbose = 0
runmode = 0
seqtype = 1
CodonFreq = 2
aaDist = 0
aaRatefile = dat/wag.dat
model = 2
NSsites = 2
icode = 0
Mgene = 0
fix_kappa = 0
kappa = 2
fix_omega = 1
omega = 0
fix_alpha = 1
alpha = 0.
Malpha = 0
ncatG = 8
getSE = 0
RateAncestor = 0
Small_Diff = .5e-6
cleandata = 0
method = 0
fix_blength = 2
"""

BSTESTCTLFILE = """seqfile = INPUTPATH
outfile = OUTPUTPATH
treefile = TREEPATH
noisy = 0
verbose = 0
runmode = 0
seqtype = 1
CodonFreq = 2
aaDist = 0
aaRatefile = dat/wag.dat
model = 2
NSsites = 2
icode = 0
Mgene = 0
fix_kappa = 0
kappa = 2
fix_omega = 0
omega = 1
fix_alpha = 1
alpha = 0.
Malpha = 0
ncatG = 8
getSE = 0
RateAncestor = 0
Small_Diff = .5e-6
cleandata = 0
method = 0
fix_blength = 2
"""

REGEX_NP1 = re.compile(r"np:\s(.*?)\)")
REGEX_LNL = re.compile(r"\):\s*(.*?)\s*\+")


def paml_pwcalc_dnds(seqs, pamltmp='pamltmp', codemlpath="codeml"):
    logfile = open("{}.log".format(pamltmp), 'w')
    ctlpath = "{}.ctl".format(pamltmp)
    seqpath = "{}.phy".format(pamltmp)
    outpath = "{}.out".format(pamltmp)
    if len(seqs) < 2 or len(seqs[0]) < 5:
        logfile.close()
        return (0, 0, 0)
    with open(seqpath, 'w') as infile:
        infile.write("   {}    {}\n".format(len(seqs), len(seqs[0]) * 3))
        for i, seq in enumerate(seqs):
            infile.write("{}  {}\n".format(i, ''.join(seq)))
    with open(ctlpath, 'w') as ctlfile:
        ctlfile.write(CTLFILE.replace(
            "INPUTPATH", seqpath).replace("OUTPUTPATH", outpath))
    subcmd = [codemlpath, ctlpath]
    proc = subprocess.Popen(' '.join(subcmd), shell=True, stdout=logfile,
                            stderr=logfile)
    print('Excecuting {}'.format(' '.join(subcmd)))
    proc.communicate()
    logfile.close()
    with open(outpath, 'r') as infile:
        for line in infile:
            if line.startswith("t="):
                arr = line.rstrip().split()
                if arr[10] == 'nan' or arr[10] == 'nan':
                    return 0., 0.
                return float(arr[10]), float(arr[13])
    return 0., 0.


def parse_branchsite(filepath=None, errorstate='ok'):
    entry = {'likelihood': -1,
             'nparam': -1,
             'fgdnds0': -1,
             'fgdnds1': -1,
             'fgdnds2a': -1,
             'fgdnds2b': -1,
             'bgdnds0': -1,
             'bgdnds1': -1,
             'bgdnds2a': -1,
             'bgdnds2b': -1,
             'dndstree': '.',
             'errorstate': 'ok',
            }
    if errorstate != 'ok':
        entry['errorstate'] = errorstate
        return entry
    with open(filepath) as infile:
        line = infile.readline()
        i = 0
        while line:
            i += 1
            if line.startswith("tree length"):
                entry['tree_length'] = line.rstrip().split()[3],
                j = 3
                while j:
                    infile.readline()
                    j -= 1
                line = infile.readline()
                entry['dndstree'] = line.strip()
            elif line.startswith('lnL'):
                entry['nparam'] = REGEX_NP1.findall(line)[0]
                entry['likelihood'] = float(REGEX_LNL.findall(line)[0])
            elif line.startswith("foreground w"):
                fgdnds = line.rstrip().split()
                entry['fgdnds0'] = fgdnds[2]
                entry['fgdnds1'] = fgdnds[3]
                entry['fgdnds2a'] = fgdnds[4]
                entry['fgdnds2b'] = fgdnds[5]
            elif line.startswith("background w"):
                bgdnds = line.rstrip().split()
                entry['bgdnds0'] = bgdnds[2]
                entry['bgdnds1'] = bgdnds[3]
                entry['bgdnds2a'] = bgdnds[4]
                entry['bgdnds2b'] = bgdnds[5]
            line = infile.readline()
    return entry


def paml_branchsite(seqs, labels, pamltmp="pamltmp",
                    species=None, speciesrev=None, mincov=8,
                    codemlpath="codeml", raxmlpath="raxmlHPC",
                    target=None, targetspec=None,
                    allsampletrees=False,
                    outgroup=None):
    target_taxa = set(target)
    logfile = open("{}.log".format(pamltmp), 'w')
    ctlpath = "{}.ctl".format(pamltmp)
    seqpath = "{}.phy".format(pamltmp)
    outpath = "{}.out".format(pamltmp)
    remx = []
    if len(seqs[0]) < 5:
        cleanup(pamltmp)
        return (parse_branchsite(errorstate='lowcov'),
                parse_branchsite(errorstate='lowcov'), '.')
    if not allsampletrees:
        best_species_rep = dict.fromkeys(species, -1)
        for i in range(len(seqs)):
            if i not in speciesrev:
                remx.append(i)
                continue
            ngap = sum([int(all(x in 'XN-' for x in y) is True)
                        for y in seqs[i]])
            if ngap == len(seqs[i]):
                remx.append(i)
            elif best_species_rep[speciesrev[i]] == -1:
                best_species_rep[speciesrev[i]] = (i, ngap + 0)
            elif ngap > best_species_rep[speciesrev[i]][1]:
                remx.append(i)
            else:
                remx.append(best_species_rep[speciesrev[i]][0])
            best_species_rep[speciesrev[i]] = (i, ngap + 0)
        for i in sorted(remx, reverse=True):
            del seqs[i]
            del labels[i]
    if len(seqs) < mincov or len(seqs[0]) < 5:
        logfile.close()
        cleanup(pamltmp)
        return (parse_branchsite(errorstate='lowcov'),
                parse_branchsite(errorstate='lowcov'), '.')

    try:
        with open(seqpath, 'w') as infile:
            infile.write("   {}    {}\n".format(len(seqs), len(seqs[0]) * 3))
            for i, seq in enumerate(seqs):
                infile.write("{}  {}\n".format(labels[i], ''.join(seq)))
        subcmd = [raxmlpath, "-s", seqpath, "-n", pamltmp, "-m" "GTRGAMMA",
                  '--no-bfgs', '-p', str(randint(100000000, 999999999))]
        print(' '.join(subcmd))
        proc = subprocess.Popen(' '.join(subcmd), shell=True, stdout=logfile,
                                stderr=logfile)
        proc.communicate()
        tree = Phylo.read('RAxML_bestTree.' + pamltmp, 'newick')
        treestr = tree.__format__("newick")
    except Exception as exception:
        cleanup(pamltmp)
        return (parse_branchsite(errorstate="raxmlfail"),
                parse_branchsite(errorstate="raxmlfail"), '.')
    try:
        bestnode = (0, None)
        tree.root_with_outgroup(outgroup)
        for node in tree.get_nonterminals():
            subnodes = set([x.name for x in node.get_terminals()])
            if not subnodes - target_taxa:
                if len(subnodes) > bestnode[0]:
                    bestnode = (len(subnodes), node)
        if bestnode[0] >= targetspec:
            bestnode[1].name = '#1'
        else:
            cleanup(pamltmp)
            return (parse_branchsite(errorstate='nobranch'),
                    parse_branchsite(errorstate='nobranch'), treestr)

        Phylo.write(tree, seqpath + '.tree.nwk', 'newick')
        treestr = tree.__format__("newick")
        print(tree.__format__('newick'))
    except Exception as exception:
        cleanup(pamltmp)
        return (parse_branchsite(errorstate="treeproc"),
                parse_branchsite(errorstate="treeproc"), treestr)
#    try:
    with open(ctlpath, 'w') as ctlfile:
        ctlfile.write(BSNULLCTLFILE.replace(
            "INPUTPATH", seqpath).replace("OUTPUTPATH", outpath).replace(
                "TREEPATH", seqpath + '.tree.nwk'))
    subcmd = [codemlpath, ctlpath]
    proc = subprocess.Popen(' '.join(subcmd), shell=True, stdout=logfile,
                            stderr=logfile)
    proc.communicate()
    pamlnull = parse_branchsite(outpath)
    print("pamlnull done ({})".format(pamlnull))
    with open(ctlpath + '.b', 'w') as ctlfile:
        ctlfile.write(BSTESTCTLFILE.replace(
            "INPUTPATH", seqpath).replace(
                "OUTPUTPATH", outpath + '.b').replace(
                    "TREEPATH", seqpath + '.tree.nwk'))
    subcmd = [codemlpath, ctlpath + '.b']
    proc = subprocess.Popen(' '.join(subcmd), shell=True, stdout=logfile,
                            stderr=logfile)
    proc.communicate()
    pamltest = parse_branchsite(outpath + '.b')
    print("pamltest done ({})".format(pamltest))
#    except:
#        cleanup(pamltmp)
#        return (parse_branchsite(errorstate='pamlfail'),
#                parse_branchsite(errorstate='pamlfail'), treestr)
    logfile.close()
    cleanup(pamltmp)
    return pamlnull, pamltest, treestr


def cleanup(pamltmp):
    """Remove temporary RAXML and PAML files
    """
    for filename in [x.format(pamltmp) for x in (
            'RAxML_info.{}', 'RAxML_parsimonyTree.{}',
            'RAxML_bestTree.{}', 'RAxML_log.{}',
            'RAxML_result.{}')]:
        try:
            os.remove(filename)
        except Exception as exception:
            pass
    return ''


if __name__ == ("__main__"):
    print("""MVF PAML wrappers, please run one of the
             other MVFtools scripts to access these functions""")
