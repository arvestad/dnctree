import argparse
import json
import os
import sys
import time
import traceback

from alv.exceptions import AlvPossibleFormatError, AlvEmptyAlignment
from alv.io import guess_format, read_alignment
from dnctree.pahmm_adaptor import pahmm_available
from dnctree import divide_n_conquer_tree, choose_distance_function
from dnctree.algtesting import run_alg_testing
from dnctree.version import __version__ as dnctree_version

if pahmm_available():
    from dnctree.msa import MSA, MSApaHMM
    from pahmm import PAHMMError
else:
    from dnctree.msa import MSA 
    class PAHMMError(Exception):
        pass


aa_models = ['WAG', 'LG', 'VT', 'JTT', 'Dayhoff', 'cpREV']
aa_default_model = aa_models[0]
pahmm_aa_models = ['WAG', 'LG', 'JTT']

dna_models =['JC', 'K2P']
dna_default_model = dna_models[0]
pahmm_dna_models = ['GTR', 'HKY85']

def cmd_line_args():
    ap = argparse.ArgumentParser(description='Divide-and-conquer phylogenetic distance-based tree inference.')
    ap.add_argument('infile', help='Multiple sequence alignment in a standard format '
                                   '(FASTA, Phylip, Nexus, Clustal, or Stockholm format).')
    ap.add_argument('-t', '--seqtype', default='guess', choices=['aa', 'dna', 'rna', 'guess'],
                    help='Set the sequence type to expect. The default is to guess the input type.')
    ap.add_argument('-m', '--model', default='guess', choices=['guess', 'kimura'] + aa_models + dna_models,
                    help=f'Choose one of the named models or let dnctree guess based on the sequence type. '
                         f'Default is {aa_default_model}.')
    ap.add_argument('-f', '--format', default='guess',
                    choices=['guess', 'fasta', 'clustal', 'nexus', 'phylip', 'stockholm'],
                    help="Specify what sequence type to assume. "
                         "Be specific if the file is not recognized automatically. Default: %(default)s")
    ap.add_argument('--pahmm', metavar='pahmm_model', choices=pahmm_aa_models+pahmm_dna_models,
                    help='Use paHMM library to calculate distances using one of the given models.'
                    'Note that paHMM implements a different set of models. Model parameters for'
                    'the DNA models, HKY85 and GTR, are inferred.')
    ap.add_argument('-j', '--json-output', action='store_true',
                    help='Output in JSON format, as an object with fields. Key object is "tree".')

    info = ap.add_argument_group('Diagnostic output')
    info.add_argument('-i', '--info', action='store_true',
                      help='Show some basic info about input and output.')
    info.add_argument('--verbose', action='store_true',
                      help='Output progress information')

    info.add_argument('-w', '--supress-warnings', action='store_true',
                      help='Do not warn about sequence pairs not sharing columns '
                           'or sequences being completely different.')

    info.add_argument('-v', '--version', action='version', version = f'dnctree {dnctree_version}')

    group = ap.add_argument_group('Expert options', 'You need to understand the dnctree '
                                                    'algorithm to tweak these options in a meaningful way.')
    group.add_argument('--base-case-size', default=100, type=int, metavar='int',
                       help='When a subproblem has at most this many taxa, full NJ is run. Default: %(default)s')
    group.add_argument('--max-n-attempts', type=int, default=1, metavar='int',
                       help='Make at most this many attempts. Default: %(default)s')
    group.add_argument('--max-clade-size', type=float, default=0.5, metavar='float',
                       help='Stop sampling triples when the largest subclade is this fraction of the '
                            'number of taxa. Default: %(default)s')
    group.add_argument('--first-triple', nargs=3, metavar='taxa',
                       help='Give three taxa to induce first subproblems.')
    group.add_argument('--secret-developer-options', action='store_true',
                       help='This option shows how to access secret expert developer options')

    if 'DNCTREE_TESTING' in os.environ:
        algtesting = ap.add_argument_group('Development options', 'Secret option setup.')
        algtesting.add_argument('--alg-testing', type=float,
                                help='Enables algorithm evaluation. The infile is read as model tree, '
                                     'defining distance, and the parameter to this option is the randomised error.')
        algtesting.add_argument('--alg-testing-base-case-sizes', default='5,10',
                                help='Write a comma-separated list of base-case sizes')
        algtesting.add_argument('--alg-testing-nj', action='store_true',
                                help='Compare with NJ. This option is dependent on --alg-testing.')

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
        print('Warning: --max-n-attempts cannot be smaller than 1. '
              'Setting that parameter to 1 and continuing.', file=sys.stderr)
        args.max_n_attempts = 1


def main():
    args = None

    try:
        args = cmd_line_args()
        check_args(args)

        verbosity = []
        if args.info:
            verbosity.append('info')
        if args.verbose:
            verbosity.append('verbose')
        if args.supress_warnings:
            verbosity.append('supress_warnings')

        if 'DNCTREE_TESTING' in os.environ:
            if args.alg_testing:
                run_alg_testing(args, verbosity)
                sys.exit()

        if args.secret_developer_options:
            print('Set the environment variable DNCTREE_TESTING to enable some additional developer options.')
            sys.exit(0)

    except KeyboardInterrupt:
        sys.exit()

    try:
        start_time = time.perf_counter()
        if args.pahmm:
            if not pahmm_available():
                print("Could not use '--pahmm' because the pahmm library is not available.", file=sys.stderr)
                print("To install pahmm, run 'python3 -m pip install pahmm'.", file=sys.stderr)
                sys.exit(1)
            model=args.pahmm
            model_name = args.pahmm
            seq_data = MSApaHMM.from_file(args.infile, model, verbosity)
        else:
            if args.format == 'guess':
                inputformat = guess_format(args.infile)
            else:
                inputformat = args.format
            alv_msa, _ = read_alignment(args.infile, args.seqtype, inputformat, None, None)
            seq_data = MSA(alv_msa)
            model_name = None
            if args.model != 'guess':
                model_name = args.model
    except KeyboardInterrupt:
        sys.exit()
    except FileNotFoundError:
        sys.exit(f'Could not read file "{args.infile}"')
    except AlvPossibleFormatError:
        sys.exit(f'File "{args.infile}" does not look like an alignment')
    except PAHMMError as e:
        sys.exit(f'Error from paHMM: {e}')
    except Exception as e:
        print('Error when reading data in dnctree:', e, file=sys.stderr)
        traceback.print_exc()
        sys.exit(1)

    try:
        if args.verbose:
            n = len(seq_data.taxa())
            print(f'Input has {n} taxa.', file=sys.stderr)

        t, aux_info = divide_n_conquer_tree(seq_data, model_name=model_name,
                                  max_n_attempts=args.max_n_attempts,
                                  max_clade_size=args.max_clade_size,
                                  base_case_size=args.base_case_size,
                                  first_triple=args.first_triple,
                                  verbose=verbosity)

        stop_time = time.perf_counter()
        if args.json_output:
            aux_info['computing-time'] = stop_time - start_time
            print(make_json_string(t, aux_info, args))
        else:
            print(t)
    except AssertionError:
        sys.exit(f'A bug has occured. Please report an issue on '
                 f'http://github.com/arvestad/dnctree, and include sample data.')
    except KeyboardInterrupt:
        sys.exit()
    except Exception as e:
        print('Error in dnctree:', e, file=sys.stderr)
        traceback.print_exc()
        sys.exit(2)


def make_json_string(t, aux_info, args):
    info = dict()
    info['version'] = f'dnctree {dnctree_version}'
    info['tree'] = str(t)
    info['infile'] = args.infile
    info['aligned'] = not(args.pahmm)
    info['base-case-size'] = args.base_case_size

    for key, val in aux_info.items():
        info[key] = val
    return json.dumps(info, indent=4)



if __name__ == "__main__":
    main()
