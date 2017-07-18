import os
import argparse
import subprocess
import tempfile
import random


def convert_fasta (handle):
    result = []
    sequence = ''
    for line in handle:
        if line.startswith('$'): # skip header line
            continue
        elif line.startswith('>') or line.startswith('#'):
            if len(sequence) > 0:
                result.append([h,sequence])
                sequence = ''   # reset
            h = line.strip('>#\n')
        else:
            sequence += line.strip(' \n')

    result.append([h,sequence]) # handle last entry

    return result


def transpose_fasta (fasta):
    # some checks to make sure the right kind of object is being sent
    if type(fasta) is not list:
        return None
    if type(fasta[0]) is not list or len(fasta[0]) != 2:
        return None

    n_columns = len(fasta[0][1])
    res = []
    for c in range(n_columns):
        res.append ( [ s[c] for h, s in fasta ] )

    return res


def untranspose_fasta(tfasta):
    nseq = len(tfasta[0])
    res = [ '' for s in range(nseq) ]
    for col in tfasta:
        for i in range(nseq):
            res[i] += col[i]
    return res



def bootstrap(fasta, reps, outfile, tmpfile):
    with open(fasta, 'rU') as handle:
        fasta = convert_fasta(handle)
        seqlen = len(fasta[0][1])
        tfasta = transpose_fasta(fasta)

    for rep in range(reps):
        print rep, reps
        sample = [random.sample(tfasta, 1)[0] for _ in range(seqlen)]
        seqs = untranspose_fasta(sample)

        # export bootstrap alignment
        handle = open(tmpfile, 'w')
        for i, s in enumerate(seqs):
            handle.write('>%s\n%s\n' % (fasta[i][0], s))
        handle.close()

        # build tree
        output = subprocess.check_output(
           ['fasttree2', '-quiet', '-gtr', '-nosupport', '-nt'],
           stdin=open(tmpfile, 'rU')
        )

        outfile.write(output)
        outfile.flush()




def main():
    parser = argparse.ArgumentParser(description='Generate a set of trees with bootstrap sampling')
    parser.add_argument('fasta', help='<input> FASTA containing multiple sequence alignment')
    parser.add_argument('outfile', type=argparse.FileType('w'), help='<output> File to write Newick strings.')
    parser.add_argument('-n', type=int, default=100, help='Number of bootstrap samples to generate.')
    parser.add_argument('-tmp', default=os.path.join(tempfile.gettempdir(), 'bootstrap.fa'),
                        help='<optional> File to write bootstrap alignment.')

    args = parser.parse_args()
    bootstrap(args.fasta, args.n, args.outfile, args.tmp)


if __name__ == '__main__':
    main()
