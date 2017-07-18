"""
Trailing whitespace in sequence headers causing problems for ClusterPicker
"""

import sys
from glob import glob
#from Bio import Phylo


files = glob('../data2/*.fas.rapidnj.nwk')

for f in files:
    print f
    with open(f, 'rU') as handle:
        line = handle.readline()
    outfile = f.replace('.rapidnj.', '.rapidnj2.')
    with open(outfile, 'w') as handle:
        handle.write(line.replace(' ', '').replace("'", ''))

sys.exit()

# remove from FASTA
files = glob('../data/*.*_TRUE.fas')

for f in files:
    print f
    with open(f, 'rU') as handle:
        lines = handle.readlines()
    
    # OVERWRITE!
    with open(f, 'w') as handle:
        for line in lines:
            handle.write(line.strip('\n').strip(' ')+'\n')

