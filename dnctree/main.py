import argparse
import sys
import traceback

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

    ap.add_argument('--verbose', action='store_true',
                    help='Output progress information')

    ap.add_argument('-w', '--supress-warnings', action='store_true',
                    help='Do not warn about sequence pairs not sharing columns or sequences being completely different.')

    group = ap.add_argument_group('Export options', 'You need to understand the dnctree algorithm to tweak these options in a meaningful way.')
    group.add_argument('--max-clade-size', type=float, default=0.5, metavar='float',
                       help='Stop sampling triples when the largest subclade is this fraction of the number of taxa. Default: %(default)s')
    group.add_argument('--max-n-attempts', type=int, default=10, metavar='int',
                       help='Make at most this many attempts. Default: %(default)s')
    group.add_argument('--base-case-size', default=100, type=int, metavar='int',
                       help='When a subproblem has at most this many taxa, full NJ is run. Default: %(default)s')

    group.add_argument('--first-triple', nargs=3, metavar='taxa',
                       help='Give three taxa to induce first subproblems.')

    args = ap.parse_args()
    return args


def check_args(args):
    '''
    Impose constraints on arguments. Some badly chosen parameters result in program exit, 
    others simply a correction to a valid parameter value.
    '''

    if args.max_clade_size <= 0.01:
        sys.exit('Error: --max-clade-size cannot be that small')
    elif args.max_clade_size > 1.0:
        print('Warning: --max-clade-size cannot be that large, changing to 1.0', file=sys.stderr)
        args.max_clade_size = 1.0

    if args.base_case_size < 3:
        sys.exit('Error: --base-case-size must be at least 3.')

    if args.max_n_attempts < 1:
        print('Warning: --max-n-attempts cannot be smaller than 1. Setting that parameter to 1 and continuing.', file=sys.stderr)
        args.max_n_attempts = 1


def main():
    try:
        args = cmd_line_args()
        check_args(args)

        if args.format == 'guess':
            inputformat = guess_format(args.infile)
        else:
            inputformat = args.format
        alv_msa, x = read_alignment(args.infile, args.seqtype, inputformat, None, None)
    except KeyboardInterrupt:
        sys.exit()
    except Exception as e:
        print('Error when reading data in dnctree:', e, file=sys.stderr)
        sys.exit(1)

    verbosity = []
    msa = MSA(alv_msa)
    if args.verbose:
        n = len(msa.taxa())
        print(f'Input MSA has {n} taxa.', file=sys.stderr)

        verbosity.append('verbose')
        if args.supress_warnings:
            verbosity.append('supress_warnings')

    try:
        t = divide_n_conquer_tree(msa, args.max_n_attempts, args.max_clade_size, args.base_case_size, args.first_triple, verbosity)
        print(t)
    except AssertionError:
        sys.exit(f'A bug has occured. Please report an issue on http://github.com/arvestad/dnctree, and include sample data.')
    except KeyboardInterrupt:
        sys.exit()
    except Exception as e:
        print('Error in dnctree:', e, file=sys.stderr)
        sys.exit(2)
