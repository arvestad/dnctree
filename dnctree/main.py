import argparse
import json
import logging
import os
import sys
import time
import traceback

import fastphylo
from alv.exceptions import AlvPossibleFormatError, AlvEmptyAlignment
from alv.io import guess_format, read_alignment
from dnctree import divide_n_conquer_tree
from dnctree.algtesting import run_alg_testing
from dnctree.version import __version__ as dnctree_version
from dnctree.msa import MSA
from dnctree.api import dna_models

# Formats fastphylo can read directly (fast C++ parser). Anything else
# (Nexus, Clustal) falls back to the alv/BioPython pipeline.
FASTPHYLO_FORMATS = ('fasta', 'phylip', 'stockholm')

FASTPHYLO_SEQTYPES = {
    'aa': fastphylo.SeqType.Protein,
    'dna': fastphylo.SeqType.DNA,
    'rna': fastphylo.SeqType.RNA,
}


aa_models = ['WAG', 'LG', 'VT', 'JTT', 'Dayhoff', 'cpREV']
aa_default_model = aa_models[0]

def cmd_line_args():
    ap = argparse.ArgumentParser(description='Divide-and-conquer phylogenetic distance-based tree inference.')
    ap.add_argument('infile', help='Multiple sequence alignment in a standard format '
                                   '(FASTA, Phylip, Nexus, Clustal, or Stockholm format).')
    ap.add_argument('-t', '--seqtype', default='guess', choices=['aa', 'dna', 'rna', 'guess'],
                    help='Set the input sequence type. The default is to guess the input type.')
    ap.add_argument('-m', '--model', default='guess', choices=['guess', 'kimura'] + aa_models + dna_models,
                    help=f'Choose one of the named models or let dnctree guess based on the sequence type. '
                         f'Default is {aa_default_model}.')
    ap.add_argument('-f', '--format', default='guess',
                    choices=['guess', 'fasta', 'clustal', 'nexus', 'phylip', 'stockholm'],
                    help="Specify what sequence type to assume. "
                         "Be specific if the file is not recognized automatically. Default: %(default)s")
    simple_group = ap.add_mutually_exclusive_group()
    simple_group.add_argument('-s', '--simple', action='store_true',
                    help='Use the simple algorithm, which uses three taxa and sorts the rest into three '
                    'subproblems, joined on a centroid vertex with no distances of its own. Worse quality '
                    'than the standard algorithm.')
    simple_group.add_argument('-so', '--simple-old', action='store_true',
                    help='Like --simple, but joined on a vertex whose distances to the rest of the taxa '
                    'are estimated - the original version of the simple algorithm.')
    simple_group.add_argument('-M', '--multi', action='store_true',
                    help='Use dnc_tree_multi, a variant of the standard (core-tree) algorithm with a '
                    'choice of non-core-taxa partition method - see --partition, --core-tree-size.')
    ap.add_argument('-j', '--json-output', action='store_true',
                    help='Output in JSON format, as an object with fields. Key object is "tree".')

    info = ap.add_argument_group('Logging')
    info.add_argument('-l', '--log-level', choices=['quiet', 'progress', 'verbose'],
                      default='quiet',
                      help="Level 'quiet', default, shows warnings and errors."
                      "At level 'progress' you get progress information and more."
                      "Level 'verbose' shows you more algorithmic details and is probably more than you want to see.")
    info.add_argument('--log-file',
                      help='Name a file to send logging info to.')

    info.add_argument('-w', '--supress-warnings', action='store_true',
                      help='Do not warn about sequence pairs not sharing columns '
                           'or sequences being completely different.')

    info.add_argument('-v', '--version', action='version', version = f'dnctree {dnctree_version}')

    group = ap.add_argument_group('Expert options', 'You need to understand the dnctree '
                                                    'algorithm to tweak these options in a meaningful way.')
    group.add_argument('-b', '--base-case-size', default=100, type=int, metavar='int',
                       help='When a subproblem has at most this many taxa, full NJ is run. Default: %(default)s')
    group.add_argument('-c', '--core-tree-size', type=int, default=None, metavar='int',
                       help='Size of the randomly sampled core taxa set, used by the core-tree and '
                            'multi algorithms. Defaults to --base-case-size if not given.')
    group.add_argument('-p', '--partition', choices=['quartet', 'tree'], default='quartet',
                       help="Non-core-taxa partition method for the multi algorithm (see --multi): "
                            "'quartet' (default) uses a 4-point quartet test against one representative "
                            "taxon per clade; 'tree' builds a fresh NJ tree per non-core taxon (much more "
                            "expensive). No effect without --multi.")
    group.add_argument('--base-case-method', choices=['nj', 'bionj'], default='bionj',
                       help="Which method the multi algorithm's terminal base case uses (see --multi). "
                            "Default: %(default)s. No effect without --multi.")
    group.add_argument('--max-n-attempts', type=int, default=1, metavar='int',
                       help='Make at most this many attempts. Default: %(default)s')
    group.add_argument('--max-clade-size', type=float, default=0.5, metavar='float',
                       help='Stop sampling triples when the largest subclade is this fraction of the '
                            'number of taxa. Default: %(default)s')
    group.add_argument('--first-triple', nargs=3, metavar='taxa',
                       help='Give three taxa to induce first subproblems when using --simple. No effect'
                            'at all without --simple.')
    group.add_argument('--secret-developer-options', action='store_true',
                       help='This option shows how to access secret expert developer options')

    if 'DNCTREE_TESTING' in os.environ:
        algtesting = ap.add_argument_group('Development options', 'Secret option setup.')
        algtesting.add_argument('--alg-testing', type=float,
                                help='Enables algorithm evaluation. The infile is read as model tree, '
                                     'defining distance, and the parameter to this option is the randomised error.')
        algtesting.add_argument('--alg-testing-algorithm', choices=['simple', 'simple_old', 'core-tree', 'multi'],
                                default='core-tree',
                                help='Decide which algorithm to test')
        algtesting.add_argument('--alg-testing-base-case-sizes', default='5,10',
                                help='Write a comma-separated list of base-case sizes')
        algtesting.add_argument('--alg-testing-nj', action='store_true',
                                help='Compare with NJ. This option is dependent on --alg-testing.')
        algtesting.add_argument('--alg-testing-err-distribution', default='uniform', choices=['uniform', 'normal', 'evil'],
                                help='Decide wether to use U[-eps, eps], N(0, eps), or {-0.99*eps, 0.99*eps} for the error term.'
                                     'I.e., uniform or normal distribution around zero, and eps is the parameter to --alg-testing,'
                                     'or "evil", meaning the we are always close to worst error.'
                                     'Note that distances are not allowed to be zero, so the error term is conditioned.')

    args = ap.parse_args()
    return args


def check_args(args):
    '''
    Impose constraints on arguments. Some badly chosen parameters result in program exit,
    others simply a correction to a valid parameter value.
    '''

    try:
        if args.secret_developer_options:
            print('Set the environment variable DNCTREE_TESTING to enable some additional developer options.')
            sys.exit(0)

        if args.max_clade_size <= 0.01:
            sys.exit('Error: --max-clade-size cannot be that small')
        elif args.max_clade_size > 1.0:
            logging.warning('--max-clade-size cannot be that large, changing to 1.0')
            args.max_clade_size = 1.0

        if args.base_case_size < 3:
            sys.exit('Error: --base-case-size must be at least 3.')

        if args.core_tree_size is not None and args.core_tree_size < 3:
            sys.exit('Error: --core-tree-size must be at least 3.')

        if args.max_n_attempts < 1:
            logging.warning('--max-n-attempts cannot be smaller than 1. '
                            'Setting that parameter to 1 and continuing.')
            args.max_n_attempts = 1
    except KeyboardInterrupt:
        sys.exit()

        
def main_load_data_and_model(args):
    '''
    Returns: the sequence data (an MSA) and the name of the model to use.

    There are several ways this function can go wrong due to user input and strange data, so
    we catch several different exceptions and give appropriate (?) feedback.
    '''
    try:
        logging.info(f'Reading data from {args.infile}')
        if args.format == 'guess':
            inputformat = guess_format(args.infile)
        else:
            inputformat = args.format

        if inputformat in FASTPHYLO_FORMATS:
            fp_seq_type = FASTPHYLO_SEQTYPES.get(args.seqtype)  # None -> auto-detect
            fp_alignment = fastphylo.read(args.infile, format=inputformat, seq_type=fp_seq_type)
            seq_data = MSA.from_fastphylo(fp_alignment)
        else:
            alv_msa, _ = read_alignment(args.infile, args.seqtype, inputformat, None, None)
            seq_data = MSA.from_alv(alv_msa)
        model_name = None
        if args.model != 'guess':
            model_name = args.model
        return seq_data, model_name
    except KeyboardInterrupt:
        sys.exit(17)
    except FileNotFoundError:
        logging.critical(f'Could not read file "{args.infile}"')
        sys.exit(f'Could not read file "{args.infile}"')
    except AlvPossibleFormatError:
        logging.critical(f'File "{args.infile}" does not look like an alignment')
        sys.exit(f'File "{args.infile}" does not look like an alignment')
    except Exception as e:
        logging.critical(f'Error when reading data in dnctree: {e}')
        sys.exit(1)


def logging_details(args):
    if args.log_level == 'quiet':
        level = logging.WARNING
    elif args.log_level == 'progress':
        level = logging.INFO
    elif args.log_level == 'verbose':
        level = logging.DEBUG
    else:
        sys.exit(f'Not a logging level: {args.log_level}')

    if args.log_file:
        logging.basicConfig(filename=args.log_file, level=level, format='%(levelname)s %(asctime)s: %(message)s')
    else:
        logging.basicConfig(stream=sys.stderr, level=level, format='%(levelname)s %(asctime)s: %(message)s')

       
def main():
    args = cmd_line_args()
    logging_details(args)
    check_args(args)

    if 'DNCTREE_TESTING' in os.environ and args.alg_testing:
        try:
            run_alg_testing(args)
        except KeyboardInterrupt:
            sys.exit()
        sys.exit()

    start_time = time.perf_counter()
    seq_data, model_name = main_load_data_and_model(args)

    try:
        n = len(seq_data.taxa())
        logging.info(f'Input has {n} taxa.')

        if args.simple:
            algorithm = 'simple'
        elif args.simple_old:
            algorithm = 'simple_old'
        elif args.multi:
            algorithm = 'multi'
        else:
            algorithm = 'core-tree'

        t, aux_info = divide_n_conquer_tree(seq_data, model_name=model_name,
                                  max_n_attempts=args.max_n_attempts,
                                  max_clade_size=args.max_clade_size,
                                  base_case_size=args.base_case_size,
                                  core_size=args.core_tree_size,
                                  base_case_method=args.base_case_method,
                                  partition_method=args.partition,
                                  algorithm=algorithm)

        logging.info('Done.')
        stop_time = time.perf_counter()
        if args.json_output:
            aux_info['computing-time (s)'] = stop_time - start_time
            print(make_json_string(t, aux_info, args))
        else:
            print(t)
    except AssertionError:
        sys.exit(f'A bug has occured. Please report an issue on '
                 f'http://github.com/arvestad/dnctree, and include sample data.')
    except KeyboardInterrupt:
        sys.exit(17)
    except Exception as e:
        logging.critical(f'Error in dnctree: {e}')
        traceback.print_exc()
        sys.exit(2)


def make_json_string(t, aux_info, args):
    info = dict()
    info['version'] = f'dnctree {dnctree_version}'
    info['tree'] = str(t)
    info['infile'] = args.infile
    info['aligned'] = True
    info['base-case-size'] = args.base_case_size
    info['core-tree-size'] = args.core_tree_size if args.core_tree_size is not None else args.base_case_size

    for key, val in aux_info.items():
        info[key] = val
    return json.dumps(info, indent=4)



if __name__ == "__main__":
    main()
