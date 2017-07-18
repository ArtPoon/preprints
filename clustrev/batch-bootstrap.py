import os
import sys
from glob import glob
import bootstrap
import tempfile
from mpi4py import MPI

my_rank = MPI.COMM_WORLD.Get_rank()
nprocs = MPI.COMM_WORLD.Get_size()

files = glob('../data/fast*.*_TRUE.fas')
tmpfile = os.path.join(tempfile.gettempdir(), '%d.fa' % (my_rank,))

for i, f in enumerate(files):
    if i % nprocs != my_rank:
        continue

    sys.stdout.write('process %d of %d starting task %s\n' % (my_rank, nprocs, f))

    outfile = open(f.replace('.fas', '.bootstrap.nwk'), 'w')
    bootstrap.bootstrap(f, 100, outfile, tmpfile)
    outfile.close()

MPI.COMM_WORLD.Barrier()
MPI.Finalize()
