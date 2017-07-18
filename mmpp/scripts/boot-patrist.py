import sys
import os
import patrist
from glob import glob
from Bio import Phylo
from mpi4py import MPI

my_rank = MPI.COMM_WORLD.Get_rank()
nprocs = MPI.COMM_WORLD.Get_size()

# ML trees from FastTree2
files = glob('../data2/*.*_TRUE.bootstrap.nwk')
cutoffs = [
    0.004, 0.006, 0.008, 0.01, 0.012, 0.014, 0.016, 0.018, 0.02, 0.022, 0.024, 0.026, 0.028,
    0.03, 0.035, 0.04, 0.045, 0.05
]  # WARNING: performance is determined by the largest cutoff value
THRESHOLD = 0.8


if my_rank == 0:
    outfile = open('../results2/patrist-boot.csv', 'w')
    outfile.write('treatment,rep,cutoff,true.neg,false.pos,false.neg,true.pos\n')

for f in files:
    # parse filename
    path, filename = os.path.split(f)
    treatment, rep, _, _ = filename.split('.')
    rep = int(rep.split('_')[0])

    trees = list(Phylo.parse(f, 'newick'))
    tipnames = set([tip.name for tip in trees[0].get_terminals()])

    # distribute trees across nodes
    sys.stdout.write('Process %d of %d analyzing trees\n' % (my_rank, nprocs))
    result = []
    for tn, tree in enumerate(trees):
        if tn % nprocs != my_rank:
            continue
        result.extend(patrist.find_short_edges(tree=tree, cutoff=max(cutoffs), returnlist=True))

    sys.stdout.write('Process %d of %d computing minimum distances\n' % (my_rank, nprocs))
    counts = dict([(cutoff, {}) for cutoff in cutoffs])
    for ii, (tip1, tip2, dist) in enumerate(result):
        for tipname in [tip1, tip2]:
            for cutoff in cutoffs:
                if tipname not in counts[cutoff]:
                    counts[cutoff].update({tipname: 0})
                counts[cutoff][tipname] += int(dist < cutoff)

    all_counts = MPI.COMM_WORLD.gather(counts, root=0)

    if my_rank == 0:
        sys.stdout.write('Collating counts\n')

        # collate counts
        collated_counts = dict([(cutoff, {}) for cutoff in cutoffs])
        for counts in all_counts:
            for cutoff in cutoffs:
                for tipname, count in counts[cutoff].iteritems():
                    if tipname not in collated_counts[cutoff]:
                        collated_counts[cutoff].update({tipname: 0})
                    collated_counts[cutoff][tipname] += count

        # output to file
        for cutoff in cutoffs:
            clustered = set()
            for tipname, count in collated_counts[cutoff].iteritems():
                if count >= THRESHOLD * len(trees):
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

    MPI.COMM_WORLD.Barrier()

if my_rank == 0:
    outfile.close()

MPI.Finalize()
