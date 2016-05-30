"""
Use Sergei's fast TN93 code.
"""
import os
from glob import glob
import subprocess

files = glob('../data/*.*_TRUE.fas')

for f in files:
    path, filename = os.path.split(f)
    print filename
    outfile = f.replace('.fas', '.tn93.csv')
    p = subprocess.check_call(['tn93', '-o', outfile, '-t', '0.5', f])
