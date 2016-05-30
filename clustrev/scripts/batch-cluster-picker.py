"""
Batch process with Cluster Picker
"""
import os
import subprocess
from glob import glob

devnull = open(os.devnull, 'w')

def clusterpick (fasta_file, tree_file, init_support, boot_support, max_dist, max_size):
    """ Call Cluster Picker and parse outputs """
    p = subprocess.check_call(map(str, [
        'java', '-jar', '../bin/ClusterPicker_1.2.4.jar', 
        fasta_file, tree_file, init_support, boot_support, max_dist, max_size
    ]), stdout=devnull)
    
    outfile = fasta_file.replace('.fas', '.fas_clusterPicks_list.txt')
    res = dict([(x, 0) for x in ['tp', 'fp', 'tn', 'fn']])
    with open(outfile, 'rU') as handle:
        _ = handle.next()  # skip first line
        for line in handle:
            label, cluster = line.strip('\n').split('\t')
            if '_1_' in label:
                if cluster == '-1':
                    res['fn'] += 1
                else:
                    res['tp'] += 1
            else:
                if cluster == '-1':
                    res['tn'] += 1
                else:
                    res['fp'] += 1
    return res


boots = [0.9, 0.95, 0.99]
dists = [0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.04, 0.05, 0.075, 0.1]
max_size = 500

files = glob('../data/*.*_TRUE.fas')

outfile = open('../results/cluster-picker.csv', 'w')
outfile.write('treatment,rep,boot.cut,dist.cut,true.neg,false.pos,false.neg,true.pos\n')

for f in files:
    path, filename = os.path.split(f)
    treatment, rep, _ = filename.split('.')
    rep = int(rep.split('_')[0])
    
    treefile = f.replace('.fas', '.fas.tree')
    for boot in boots:
        for dist in dists:
            print treatment, rep, boot, dist
            res = clusterpick(f, treefile, boot, boot, dist, max_size)
            outfile.write('%s,%d,%.3f,%.3f,%d,%d,%d,%d\n' % (
                treatment, rep, boot, dist, res['tn'], res['fp'], res['fn'], res['tp']
            ))
            outfile.flush()

outfile.close()
