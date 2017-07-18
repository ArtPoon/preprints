"""
Parse and modify BEAST/MASTER NEXUS tree outputs by transferring the contents of the 
comment strings to the node labels.  The purpose of this is to make these outputs 
more readable by BioPython.

Terminal node label has the format:
1[&type="I",location="0",reaction="Sampling",time=0.2458663625884157]:0.12195191891543124
Internal node label has the format:
...)[&type="I",location="0",reaction="Infection",time=0.15089530647673044]:0.011186556452893681
"""
import re
import sys

from glob import glob


def parse_nexus(handle):
    """
    Generator that parses each line of NEXUS input and updates 
    node names with information extracted from comments.
    Output only the Newick tree strings.
    """
    for line in handle:
        node_count = 0  # for indexing internal nodes
        
        items = line.split()
        if len(items) < 3:
            continue
        line = items[-1]  # exclude the tree prefix
        
        has_tree = False
        while True:
            m = re.search('\)*[0-9]*\[[^]]*\]', line)  # node name
            if m is None:
                break
                
            has_tree = True
            node_name = m.group()
            left, right = m.span()
        
            is_internal = node_name.startswith(')')
            tip_index = None
            if is_internal:
                node_name = node_name.lstrip(')')  # remove bracket
            else:
                m2 = re.match('[0-9]+', node_name)
                if m2 is not None:
                    tip_index = m2.group()
                    node_name = node_name[len(tip_index):]
            
            tokens = node_name.strip('[]').split(',')
            values = dict([token.lstrip('&').replace('"','').split('=') for token in tokens])
            
            
            insert = '%s%d_%s_%s_%s' % (
                ')N' if is_internal else 'T', 
                node_count if is_internal else int(tip_index),
                values.get('type', ''),
                values.get('location', ''),
                values.get('reaction', '')
            )
            if is_internal:
                insert = ')'
            
            line = line[:left] + insert + line[right:]
            if is_internal:
                node_count += 1
        
        if has_tree:
            yield line


nexus = glob('../data2/*.nexus')

for f in nexus:
    print (f)
    with open(f, 'rU') as handle:
        outfile = open(f.replace('.nexus', '.labeled.nwk'), 'w')
        for i, line in enumerate(parse_nexus(handle)):
            print (i)
            outfile.write(line+'\n')
        outfile.close()
    
