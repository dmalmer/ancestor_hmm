
import argparse
import math
import numpy

from collections import namedtuple
from datetime import datetime
from time import time

from prob import calc_recomb_rate, calc_new_emit_p
from util import atof, create_grid_range, get_emit_key, get_states, read_recomb_rates, read_SNPs, read_SVs
from hmm import expectation_maximization, viterbi


# Arguments
def read_arguments():
    parser = argparse.ArgumentParser(add_help=True, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input-file', help='Input SNP data file (BED file format)', type=str, required=True)
    parser.add_argument('-d', '--desc-strain', help='Name of descendant in input SNP file', type=str,
                        required=True)
    parser.add_argument('-o', '--output-dir', help='Directory for output files (if not specified, output directory will '
                                                   'default to the same directory as the input file)', type=str)

    parser.add_argument('-t', '--trans-in-p', help='Starting trans-in probability (can be a value or range of values in '
                                                   'the form "[x-y]", which is divided into parts with the -gs flag)',
                                                   type=atof, default=.64)
    parser.add_argument('-e', '--emit-same-p', help='Starting emit-same probability (can be a value or range of values '
                                                    'in the form "[x-y]", which is divided into parts with the -gs flag)',
                                                    type=atof, default=.99)
    parser.add_argument('-m', '--max-iter', help='Maximum number of EM iterations', type=int, default=50)
    parser.add_argument('-c', '--prob-dist-cutoff', help='Probability distance cutoff to end EM loop', type=float,
                        default=.001)

    parser.add_argument('-p', '--parallel', help='Run viterbi algorithm over each chromosome in parallel',
                        action='store_true')

    parser.add_argument('-r', '--recomb-rates-file', help='Input file with recombination rates to be used as priors for '
                                                          'transition probabilities', type=str)
    parser.add_argument('-a', '--adjust-recomb', help='Multiplier to adjust expected number of recombinations (can be a '
                                                      'value or range of values in the form "[x-y]", which is divided '
                                                      'into parts with the -gs flag)', type=atof, default=1.)
    parser.add_argument('-u', '--unk-cutoff', help='Cutoff for fraction of Unk SNPs required for an ancestor block to '
                                                   'be relabeled as Unk (can be a value or range of values in the form '
                                                   '"[x-y]", which is divided into parts with the -gs flag)', type=atof,
                                                   default=1.)

    parser.add_argument('-ep', '--effective-pop', help='Effective population (N_e) used in recombination rate '
                                                       'calculations', type=int, default=1)
    parser.add_argument('-ng', '--num-generations', help='Estimated number of generations between ancestors and '
                                                         'descendant used in recombation rate calculations', type=int,
                                                         default=1)

    parser.add_argument('-si', '--sv-insertions-file', help='Input file for insertion structural variants used to score '
                                                            'HMM results', type=str)
    parser.add_argument('-sd', '--sv-deletions-file', help='Input file for deletion structural variants used to score '
                                                            'HMM results', type=str)

    parser.add_argument('-gs', '--grid-size', help='Number of items to divide a range of input values into', type=int,
                        default=2)

    parser.add_argument('-ii', '--ignore-inconsistent', help='Ignore SNPs that are present in the ancestors, but not in '
                                                             'the descendants. This can be useful when using SNP data '
                                                             'that doesn\'t fully cover the genome (eg. SNP chip data).',
                                                             action='store_true')

    parser.add_argument('-ad', '--append-date', help='Append date to output filename', action='store_true')
    parser.add_argument('-ap', '--append-params', help='Append string to output filename based on the input parameters',
                        action='store_true')
    parser.add_argument('-w', '--write-iter', help='Calculate scores and write to output file at each iteration within '
                                                   'the EM loop', action='store_true')
    parser.add_argument('-v', '--verbose', help='Verbose', action='store_true')

    return parser.parse_args()


# Main method
if __name__ == '__main__':
    # Start timer
    time_start = time()

    # Constants
    STATE_RGBS = {'Unk': '128,128,128', 'A': '0,153,0', 'AKR': '51,102,255', 'BALBc': '255,255,51',
                  'C3HHe': '255,153,51', 'C57BL6N': '102,0,204', 'DBA2': '255,51,51', 'IBA': '153,255,255'}
    
    # Read in arguments
    args = read_arguments()
    if '/' in args.input_file:
        filename_in = args.input_file.rsplit('/', 1)[1]
        output_dir = args.output_dir if args.output_dir else args.input_file.rsplit('/', 1)[0] + '/'
    else:
        filename_in = args.input_file
        output_dir = args.output_dir if args.output_dir else './'

    # Read in data
    SNPs_by_chr = read_SNPs(args.input_file)
    states = get_states(SNPs_by_chr, args.desc_strain, True)

    recomb_rate_dict = {}
    use_recomb_rates = False
    if args.recomb_rates_file:
        recomb_rate_dict = read_recomb_rates(args.recomb_rates_file)
        use_recomb_rates = True

    desc_ins_by_chr = {}
    desc_del_by_chr = {}
    anc_ins_by_chr = {}
    anc_del_by_chr = {}
    if args.sv_insertions_file or args.sv_deletions_file:
        if not (args.sv_insertions_file and args.sv_deletions_file):
            raise Exception('Passed in either an SV insertions (-si) or SV deletions file (-sd), but missing the other')
        desc_ins_by_chr, desc_del_by_chr, anc_ins_by_chr, anc_del_by_chr = \
            read_SVs(args.sv_insertions_file, args.sv_deletions_file, states, args.desc_strain)

    # Use a dict for EM inputs to clean up function calls
    # I would like to use a namedtuple here, but parallelpython doesn't seem to like them
    start_params = {
        'desc_strain': args.desc_strain,
        'states': states,
        'use_recomb_rates': use_recomb_rates,
        'recomb_rate_dict': recomb_rate_dict,
        'max_iter': args.max_iter,
        'prob_dist_cutoff': args.prob_dist_cutoff,
        'parallel': args.parallel,
        'effective_pop': args.effective_pop,
        'num_generations': args.num_generations,
        'desc_ins_by_chr': desc_ins_by_chr,
        'desc_del_by_chr': desc_del_by_chr,
        'anc_ins_by_chr': anc_ins_by_chr,
        'anc_del_by_chr': anc_del_by_chr,
        'ignore_inconsistent': args.ignore_inconsistent,
        'write_iter': args.write_iter,
        'filename_in': filename_in,
        'output_dir': output_dir,
        'state_rgbs': STATE_RGBS,
        'time_start': time_start,
        'verbose': args.verbose
    }

    # Create lists of input parameters for grid searching (lists will be of size 1 if there is no range of values)
    trans_in_p_range = [args.trans_in_p]
    emit_same_p_range = [args.emit_same_p]
    adjust_recomb_range = [args.adjust_recomb]
    unk_cutoff_range = [args.unk_cutoff]

    # Set parameter ranges for grid searches
    if not isinstance(args.trans_in_p, float):
        trans_in_p_range = create_grid_range(args.trans_in_p, args.grid_size)
        use_auto_str = True
    if not isinstance(args.emit_same_p, float):
        emit_same_p_range = create_grid_range(args.emit_same_p, args.grid_size)
        use_auto_str = True
    if not isinstance(args.adjust_recomb, float):
        adjust_recomb_range = create_grid_range(args.adjust_recomb, args.grid_size)
        use_auto_str = True
    if not isinstance(args.unk_cutoff, float):
        unk_cutoff_range = create_grid_range(args.unk_cutoff, args.grid_size)
        use_auto_str = True

    # Kickoff EM loop across all ranges of parameters
    job_server = None
    vit_func = None
    if args.parallel:
        import pp
        job_server = pp.Server()
        vit_func = pp.Template(job_server, viterbi, depfuncs=(calc_recomb_rate, get_emit_key),
                               modules=('numpy', 'math'))
    
    for trans_in_p in trans_in_p_range:
        for emit_same_p in emit_same_p_range:
            for adjust_recomb in adjust_recomb_range:
                for unk_cutoff in unk_cutoff_range:
                    append_str = ''
                    if args.append_params:
                        append_str += '_%.2ft-%.2fe-%.2fu-%.2fa' % (trans_in_p, emit_same_p, unk_cutoff, adjust_recomb)
                    if args.append_date:
                        append_str += datetime.now().strftime('_%y-%m-%d_%H-%M')
                    expectation_maximization(SNPs_by_chr, trans_in_p, emit_same_p, adjust_recomb, unk_cutoff, append_str,
                                             start_params, job_server, vit_func)

    if args.parallel:
        job_server.wait()
