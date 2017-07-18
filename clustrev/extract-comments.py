"""
Extract node annotations from comment fields in BEAST/MASTER NEXUS files 
and output in CSV format.

Terminal node label has the format:
1[&type="I",location="0",reaction="Sampling",time=0.2458663625884157]:0.12195191891543124
Internal node label has the format:
...)[&type="I",location="0",reaction="Infection",time=0.15089530647673044]:0.011186556452893681
"""

import re
import sys
import os
from glob import glob

files = glob('/Users/art/git/papers/clustrev/data/*.nexus')

outfile = open('extract-comments.csv', 'w')
outfile.write('treatment,rep,is.internal,type,location,reaction,time,branch.length\n')


for f in files:
    path, filename = os.path.split(f)
    treatment = filename.split('.')[0]
    
    handle = open(f, 'rU')
    tree_count = 0
    for line in handle:
        if not line.startswith('tree'):
            continue
        
        # parse comment fields for each node
        matches = re.finditer('\[[^]]*\]:[-.E0-9]*', line)  # re.MatchObject
        for m in matches:
            left, right = m.span()
            is_internal = (line[left-1] == ')')
    
            comment, branch_length = m.group().split(':')
            tokens = comment.strip('[&]').split(',')
            values = dict([token.replace('"','').split('=') for token in tokens])
    
            outfile.write(','.join(map(str, [treatment, 
                tree_count, is_internal, values['type'], values['location'],
                values['reaction'], values['time'], branch_length
            ])))
            outfile.write('\n')
            
        tree_count += 1  # update count

outfile.close()
