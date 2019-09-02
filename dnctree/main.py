import argparse
import sys

from alv.io import guess_format, read_alignment
from dnctree import divide_n_conquer_tree
from dnctree.msa import MSA

def cmd_line_args():
    ap = argparse.ArgumentParser(description='Divide-and-conquer phylogenetic distance-based tree inference. Currently, the only implemented sequence evolution model is Poisson.')
    ap.add_argument('infile', help='Multiple sequence alignment in a standard format (FASTA, Phylip, Nexus, Clustal, or Stockholm format).')
    ap.add_argument('-t', '--seqtype', default='guess', choices=['aa', 'dna', 'rna', 'guess'],
                    help='Set the sequence type to expect. The default is to guess the input type.')
    ap.add_argument('-f', '--format', default='guess',
                    choices=['guess', 'fasta', 'clustal', 'nexus', 'phylip', 'stockholm'],
                    help="Specify what sequence type to assume. Be specific if the file is not recognized automatically. Default: %(default)s")
    ap.add_argument('--max_clade_size', type=float, default=0.5,
                    help='Stop sampling triples when the largest subclade is this fraction of the number of taxa.')
    ap.add_argument('--max_n_attempts', type=int, default=100,
                    help='Make at most this many attempts. Only applies when using --random.')
    ap.add_argument('--base_case_size', default=100, type=int,
                    help='Just run NJ when a subproblem has less than this many taxa. Default: %(default)s')
    ap.add_argument('--verbose', action='store_true',
                    help='Output progress information')

    ap.add_argument('--first_triple', nargs=3, metavar='taxa',
                    help='Give three taxa to induce first subproblems.')

    args = ap.parse_args()
    return args


def main():
    args = cmd_line_args()

    try:
        if args.format == 'guess':
            inputformat = guess_format(args.infile)
        else:
            inputformat = args.format
        alv_msa, x = read_alignment(args.infile, args.seqtype, inputformat, None, None)
    except Exception as e:
        print('Error in dnctree:', e, file=sys.stderr)
        sys.exit(1)

    msa = MSA(alv_msa)

    try:
        t = divide_n_conquer_tree(msa, args.max_n_attempts, args.max_clade_size, args.base_case_size, args.first_triple, args.verbose)
        print(t)
    except Exception as e:
        print('Error in dnctree', file=sys.stderr)
        print(e, file=sys.stderr)
        sys.exit(2)
