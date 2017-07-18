import sys
import os
import patrist
from glob import glob
from Bio import Phylo

# ML trees from FastTree2
files = glob('../data/*.*_TRUE.fas.tree')
cutoffs = [
    0.004, 0.006, 0.008, 0.01, 0.012, 0.014, 0.016, 0.018, 0.02, 0.022, 0.024, 0.026, 0.028,
    0.03, 0.035, 0.04, 0.045, 0.05, 0.075, 0.1
]

outfile = open('../results/patrist.csv', 'w')
outfile.write('treatment,rep,cutoff,true.neg,false.pos,false.neg,true.pos\n')

for f in files:
    # parse filename
    path, filename = os.path.split(f)
    treatment, rep, _, _ = filename.split('.')
    rep = int(rep.split('_')[0])

    tree = Phylo.read(f, 'newick')
    tipnames = set([tip.name for tip in tree.get_terminals()])
    for cutoff in cutoffs:
        print treatment, rep, cutoff  # monitor progress

        # find all tips that are within cutoff distance of another tip
        clustered = set()
        result = patrist.find_short_edges(tree=tree, cutoff=cutoff)
        for key, dist in result.iteritems():
            for tipname in key:
                clustered.add(tipname)

        complement = tipnames.difference(clustered)  # not clustered

        true_neg = sum(map(lambda x: '_0_' in x, complement))
        false_pos = sum(map(lambda x: '_0_' in x, clustered))
        false_neg = len(complement) - true_neg
        true_pos = len(clustered) - false_pos

        outfile.write('%s,%d,%f,%d,%d,%d,%d\n' % (
            treatment, rep, cutoff, true_neg, false_pos, false_neg, true_pos
        ))
        outfile.flush()

outfile.close()
