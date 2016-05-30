"""
Parse bootstrap NJ outputs from MEGA
"""
import sys
import os
from glob import glob
from Bio import Phylo



# Nj trees generated with 1000 bootstrap samples
files = glob('../data/*.tn93-nj(2).nwk')
bootcuts = [0.8, 0.9, 0.95]
distcuts = [0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.0075, 0.008, 0.009, 0.01, 0.0125, 0.015, 0.02]

outfile = open('../results/nj-clades.csv', 'w')
outfile.write('treatment,rep,boot.cut,dist.cut,true.neg,false.pos,false.neg,true.pos\n')

for f in files:
    path, filename = os.path.split(f)
    treatment, rep, _, _ = filename.split('.')
    if treatment == 'Control':
        continue

    rep = int(rep.split('_')[0])

    tree = Phylo.read(f, 'newick')
    tree.root_at_midpoint()
    tipnames = set([tip.name for tip in tree.get_terminals()])

    nodes = tree.get_nonterminals(order='preorder')

    for bootcut in bootcuts:
        for distcut in distcuts:
            clustered = set()

            for node in nodes:
                if node.confidence < bootcut:
                    continue

                bl = [n.branch_length for n in node.find_clades(branch_length=True)]
                mean_bl = sum(bl) / len(bl)
                if mean_bl > distcut:
                    continue

                my_tips = node.get_terminals()
                for tip in my_tips:
                    clustered.add(tip.name)

            complement = tipnames.difference(clustered)  # not clustered

            true_neg = sum(map(lambda x: '_0_' in x, complement))
            false_pos = sum(map(lambda x: '_0_' in x, clustered))
            false_neg = len(complement) - true_neg
            true_pos = len(clustered) - false_pos

            outfile.write('%s,%d,%f,%f,%d,%d,%d,%d\n' % (
                treatment, rep, bootcut, distcut, true_neg, false_pos, false_neg, true_pos
            ))
            outfile.flush()

outfile.close()
