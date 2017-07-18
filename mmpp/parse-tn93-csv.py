"""
Parse TN93 CSV outputs.
Find the shortest genetic distance for each individual.
"""
import os
from glob import glob

def parse(handle):
    for line in handle:
        id1, id2, dist = map(lambda x: x.strip(), line.strip('\n').split(','))
        dist = float(dist)
        for idn in [id1, id2]:
            if idn not in mindist:
                mindist.update({idn: 100.0})
            if dist < mindist[idn]:
                mindist[idn] = dist
    return mindist


files = glob('../data2/*.tn93.csv')

outfile = open('../results2/tn93.csv', 'w')
outfile.write('treatment,rep,label,distance\n')

for f in files:
    # parse filename
    path, filename = os.path.split(f)
    treatment, rep, _, _ = filename.split('.')
    rep = int(rep.split('_')[0])
    print treatment, rep

    mindist = {}
    handle = open(f, 'rU')
    _ = handle.next()  # skip header

    # get the shortest distance for each individual
    mindist = parse(handle)

    # output to file
    for label, dist in mindist.iteritems():
        outfile.write('%s,%d,%s,%s\n' % (treatment, rep, label, dist))

outfile.close()
