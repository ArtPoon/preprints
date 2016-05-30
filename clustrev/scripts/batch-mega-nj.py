"""
Batch process MEGA neighbor-joining analysis with bootstrap 
sampling using MPI.
"""
import sys
import subprocess
from glob import glob
from mpi4py import MPI

my_rank = MPI.COMM_WORLD.Get_rank()
nprocs = MPI.COMM_WORLD.Get_size()
devnull = open('/dev/null', 'w')

maofile = 'infer_NJ_nucleotide.mao'

files = glob('../data/*.*_TRUE.fas')

for i, f in enumerate(files):
    if i % nprocs != my_rank:
        continue
    sys.stdout.write('process %d of %d starting task %s\n' % (my_rank, nprocs, f))
    outfile = f.replace('.fas', '.tn93-nj.nwk')
    #os.system('megacc -a %s -d %s -o %s' % (maofile, f, outfile))
    p = subprocess.check_call(['megacc', '-a', maofile, '-d', f, '-o', outfile], stdout=devnull)
    

MPI.COMM_WORLD.Barrier()
MPI.Finalize()
