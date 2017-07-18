"""
Batch process with phyloPart
"""
import os
import subprocess
from glob import glob
import tempfile

tmp = os.path.join(tempfile.gettempdir(), 'phylopart.out')
devnull = open(os.devnull, 'w')

def phylopart (tree_file, max_dist, outfile):
    """ Call phylopart and parse outputs """
    p = subprocess.check_call(map(str, [
        'java', '-jar', '-Djava.awt.headless=true', '../bin/phyloPart.jar', tree_file, max_dist, '-o'+outfile
    ]), stdout=devnull)
    
    res = dict([(x, 0) for x in ['tp', 'fp', 'tn', 'fn']])
    with open(outfile, 'rU') as handle:
        _ = handle.next()  # skip first line
        for line in handle:
            cluster, bootstrap, label, path, median, size = line.strip('\n').split(';')
            if '_1_' in label:
                if cluster == '0':
                    res['fn'] += 1
                else:
                    res['tp'] += 1
            else:
                if cluster == '0':
                    res['tn'] += 1
                else:
                    res['fp'] += 1
    return res


dists = [0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.04, 0.05, 0.075, 0.1]
dists += [0.001, 0.002, 0.003, 0.004, 0.0075, 0.0125]
dists.sort()

# which trees to analyze
files = glob('../data2/*.*_TRUE.fas.tree')  # ML from fasttree2

outfile = open('../results2/phylopart.csv', 'a')
#outfile.write('treatment,rep,dist.cut,true.neg,false.pos,false.neg,true.pos\n')

for f in files:
    path, filename = os.path.split(f)
    treatment, rep = filename.split('.')[:2]
    rep = int(rep.split('_')[0])
    
    for dist in dists:
        print treatment, rep, dist
        res = phylopart(f, dist, tmp)
        outfile.write('%s,%d,%.3f,%d,%d,%d,%d\n' % (
            treatment, rep, dist, res['tn'], res['fp'], res['fn'], res['tp']
        ))
        outfile.flush()

outfile.close()
