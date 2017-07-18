"""
Batch INDELible
Write control.txt file and execute indelible to simulate sequences.
"""
import sys
import os
from glob import glob
import re
import subprocess
import tempfile

from mpi4py import MPI

my_rank = MPI.COMM_WORLD.Get_rank()
nprocs = MPI.COMM_WORLD.Get_size()
devnull = open('/dev/null', 'w')
tempdir = tempfile.gettempdir()

## simulation parameters

# expected number of substitutions at a codon site in the entire tree
scaling_factor = 15.0
"""
On examining subtype B sequences from BC database, I've found the mean TN93
distance to be about 0.055 (IQR 0.047-0.062).
I experimented with different scaling factors to try to get to this mean TN93 
among simulated sequences.
12.3 gave 0.0475200
15.0 gave 0.059120 (IQR 0.050610-0.069950)  <- going with this one
20.0 was definitely too high
"""


# Gamma(shape=1.5, rate=3), histogram with 50 breaks
omegas = [0.0510204, 0.1530612, 0.2551020, 0.3571428, 0.4591836, 0.5612244, 0.6632652, 0.7653060, 0.8673468,
0.9693876, 1.0714284, 1.1734692, 1.2755100, 1.3775508, 1.4795916, 1.5816324, 1.6836732, 1.7857140,
1.8877548, 1.9897956, 2.0918364, 2.1938772, 2.2959180, 2.3979588, 2.4999996, 2.6020404, 2.7040812,
2.8061220, 2.9081628, 3.0102036, 3.1122444, 3.2142852, 3.3163260, 3.4183668, 3.5204076, 3.6224484,
3.7244892, 3.8265300, 3.9285708, 4.0306116, 4.1326524, 4.2346932, 4.3367340, 4.4387748, 4.5408156,
4.6428564, 4.7448972, 4.8469380, 4.9489788, 5.0510196]

prop = [
    1.063764e-01, 1.464867e-01, 1.401630e-01, 1.223917e-01, 1.023007e-01, 8.333079e-02, 6.673162e-02, 
    5.279576e-02, 4.139382e-02, 3.222714e-02, 2.495008e-02, 1.922789e-02, 1.476163e-02, 1.129629e-02,
    8.620596e-03, 6.562952e-03, 4.985996e-03, 3.780965e-03, 2.862470e-03, 2.163923e-03, 1.633688e-03,
    1.231906e-03, 9.279287e-04, 6.982665e-04, 5.249686e-04, 3.943506e-04, 2.960034e-04, 2.220245e-04,
    1.664245e-04, 1.246710e-04, 9.333909e-05, 6.984377e-05, 5.223631e-05, 3.904912e-05, 2.917804e-05,
    2.179306e-05, 1.627075e-05, 1.214322e-05, 9.059520e-06, 6.756626e-06, 5.037501e-06, 3.754636e-06,
    2.797655e-06, 2.084012e-06, 1.551998e-06, 1.155507e-06, 8.600995e-07, 6.400653e-07, 4.762155e-07
]

kappa = 8.0

# absolute path to write control file
#control_file = os.path.join(os.getcwd(), 'control.txt')

def sci_to_float(match):
    """
    Convert scientific notation string to float.
    """
    return format(float(match.group()), '.12f')


def indelible (tree_string, control_file, outfile):
    """
    Write control file and do a system call to run INDELIBLE.
    :param: handle is an open file to write the INDELIBLE control text
    """
    
    # process the tree string
    ts = tree_string.rstrip('-E0123456789.;:\n')
    
    # INDELIBLE can't parse scientific notation, e.g., "2.4998751172983225E-5"
    ts = re.sub('[0-9]+\.[0-9]+E-[0-9]+', sci_to_float, ts)
    
    with open(control_file, 'w') as handle:
        # write minimal contents of INDELible control file
        handle.write('[TYPE] CODON 1\n')
        handle.write('[SETTINGS]\n[output] FASTA\n[printrates] FALSE\n[randomseed] 41\n')
        handle.write('[MODEL] M3\n[submodel] %f\n' % kappa)

        prop_string = ' '.join(map(lambda p: '%1.9f' % p, prop))
        omega_string = ' '.join(map(str, omegas))

        handle.write(' %s\n %s\n' % (prop_string, omega_string))

        #handle.write('[indelmodel] NB 0.2 4\n')
        handle.write('[indelmodel] LAV 1.5 4\n')
        handle.write('[indelrate] 0.001\n')

        handle.write('[TREE] bigtree %s;\n' % ts)
        handle.write('[treelength] %1.5f\n' % scaling_factor)
        handle.write('[PARTITIONS] partitionname\n')
        #handle.write('  [bigtree M3 %d]\n' % n_sites)  # treename name rootlength
        handle.write("  [bigtree M3 hxb2_pol.txt]\n")
        handle.write('[EVOLVE] partitionname 1 %s\n' % (outfile, ))
    
    p = subprocess.check_call(['indelible', control_file], stdout=devnull)
    #os.system('/Users/art/src/INDELibleV1.03/src/indelible')  # can also use lldb


# iterate over simulated trees
newicks = glob('../data2/*.labeled.nwk')

for f in newicks:
    path, filename = os.path.split(f)
    treatment = filename.split('.')[0]
    
    with open(f, 'rU') as f:
        for rep, tree in enumerate(f):
            if rep % nprocs != my_rank:
                continue
            if len(tree.strip()) < 100:
                sys.stdout.write('process %d of %d skipping truncated tree\n' % (my_rank, nprocs))
                continue
            
            sys.stdout.write('process %d of %d starting task %s %d\n' % (my_rank, nprocs, treatment, rep))
            
            control_file = os.path.join(tempdir, 'control_%s_%d.txt' % (treatment, rep))
            outfile = os.path.join(path, '%s.%d' % (treatment, rep))
            indelible(tree, control_file, outfile)

MPI.COMM_WORLD.Barrier()
MPI.Finalize()
