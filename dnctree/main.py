import argparse
import json
import logging
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
                    help='Set the input sequence type. The default is to guess the input type.')
    ap.add_argument('-m', '--model', default='guess', choices=['guess', 'kimura'] + aa_models + dna_models,
                    help=f'Choose one of the named models or let dnctree guess based on the sequence type. '
                         f'Default is {aa_default_model}.')
    ap.add_argument('-f', '--format', default='guess',
                    choices=['guess', 'fasta', 'clustal', 'nexus', 'phylip', 'stockholm'],
                    help="Specify what sequence type to assume. "
                         "Be specific if the file is not recognized automatically. Default: %(default)s")
    ap.add_argument('-s', '--simple', action='store_true',
                    help='Use the simple algorithm, which uses three taxa and sorts the rest into three'
                    'subproblems. Worse quality than the standard algorithm.')
    ap.add_argument('--pahmm', metavar='pahmm_model', choices=pahmm_aa_models+pahmm_dna_models,
                    help='Experimental: Use paHMM library to calculate distances using one of the given models.'
                    'Note that paHMM implements a different set of models. Model parameters for'
                    'the DNA models, HKY85 and GTR, are inferred.')
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
    group.add_argument('--base-case-size', default=100, type=int, metavar='int',
                       help='When a subproblem has at most this many taxa, full NJ is run. Default: %(default)s')
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
        algtesting.add_argument('--alg-testing-algorithm', choices=['simple', 'core-tree'],
                                default='core-tree',
                                help='Decide which algorithm to test')
        algtesting.add_argument('--alg-testing-base-case-sizes', default='5,10',
                                help='Write a comma-separated list of base-case sizes')
        algtesting.add_argument('--alg-testing-nj', action='store_true',
                                help='Compare with NJ. This option is dependent on --alg-testing.')
        algtesting.add_argument('--alg-testing-err-distribution', default='uniform', choices=['uniform', 'normal'],
                                help='Decide wether to use U[-eps, eps] or N(0, eps) for the error term.'
                                     'I.e., uniform or normal distribution around zero, and eps is the parameter to --alg-testing.'
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

        if args.max_n_attempts < 1:
            logging.warning('--max-n-attempts cannot be smaller than 1. '
                            'Setting that parameter to 1 and continuing.')
            args.max_n_attempts = 1

        if 'DNCTREE_TESTING' in os.environ:
            if args.alg_testing:
                run_alg_testing(args)
                sys.exit()
    except KeyboardInterrupt:
        sys.exit()

        
def main_load_data_and_model(args):
    '''
    The model is in part decided by the data. In case we want to try PaHMM, we have to deal with
    that in a special way.

    Returns: the sequence data (typically an MSA) and the name of the model to use.

    There are several ways this function can go wrong due to user input and strange data, so
    we catch several different exceptions and give appropriate (?) feedback.
    '''
    try:
        logging.info(f'Reading data from {args.infile}')
        if args.pahmm:
            if not pahmm_available():
                logging.critical("Could not use '--pahmm' because the pahmm library is not available."
                                "To install pahmm, run 'python3 -m pip install pahmm'.")
                sys.exit(1)
            model=args.pahmm
            model_name = args.pahmm
            seq_data = MSApaHMM.from_file(args.infile, model)
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
        return seq_data, model_name
    except KeyboardInterrupt:
        sys.exit(17)
    except FileNotFoundError:
        logging.critical(f'Could not read file "{args.infile}"')
        sys.exit(f'Could not read file "{args.infile}"')
    except AlvPossibleFormatError:
        logging.critical(f'File "{args.infile}" does not look like an alignment')
        sys.exit(f'File "{args.infile}" does not look like an alignment')
    except PAHMMError as e:
        logging.critical(f'Error from paHMM: {e}')
        sys.exit(f'Error from paHMM: {e}')
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

    start_time = time.perf_counter()
    seq_data, model_name = main_load_data_and_model(args)

    try:
        n = len(seq_data.taxa())
        logging.info(f'Input has {n} taxa.')

        t, aux_info = divide_n_conquer_tree(seq_data, model_name=model_name,
                                  max_n_attempts=args.max_n_attempts,
                                  max_clade_size=args.max_clade_size,
                                  base_case_size=args.base_case_size,
                                  simple_alg=args.simple)

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
    info['aligned'] = not(args.pahmm)
    info['base-case-size'] = args.base_case_size

    for key, val in aux_info.items():
        info[key] = val
    return json.dumps(info, indent=4)



if __name__ == "__main__":
    main()
