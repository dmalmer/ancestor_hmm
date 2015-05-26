
import argparse
import math
import numpy
import pp
import sys

from datetime import datetime
from math import log
from os import environ
from time import time

from hmm_output import write_ancestors, write_scores, write_statistics
from hmm_prob import calc_new_emit_p, calc_new_trans_p, calc_recomb_rate, prob_dist, reclassify_ibd_and_unk, score_results
from hmm_util import atof, create_grid_range, get_emit_key, get_states, prob_tuples, read_recomb_rates, read_SNPs, read_SVs


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

    parser.add_argument('-ad', '--append-date', help='Append date to output filename', action='store_true')
    parser.add_argument('-ap', '--append-params', help='Append string to output filename based on the input parameters',
                        action='store_true')
    parser.add_argument('-w', '--write-iter', help='Calculate scores and write to output file at each iteration within '
                                                   'the EM loop', action='store_true')
    parser.add_argument('-v', '--verbose', help='Verbose', action='store_true')

    return parser.parse_args()


# Viterbi algorithm
def viterbi(SNPs, states, trans_p, emit_p, desc_strain, recomb_rate_dict, effective_pop, num_generations,
            recomb_adjustment, use_recomb_rates, verbose):
    # Initialize
    prob_nodes = numpy.zeros(len(SNPs), dtype={'names': states, 'formats': ['f8']*len(states)})
    ancestors_by_state = {}
    recomb_index = None  # Initially set to None so calc_recomb_rate uses a special case when called for the first time

    # Start probabilities
    for s in states:
        # For each state, the emission probability is either emit_p[state] or emit_p[~state]
        emit_key = get_emit_key(s, SNPs[0][3], desc_strain)
        # Probability of a given state at SNPs[0] is emit prob of stat
        prob_nodes[0][s] = emit_p[s][emit_key]
        ancestors_by_state[s] = [s]

    # Rest of probabilities
    for i in range(1, len(SNPs)):
        #if verbose and i % 1000 == 0:
        #    print 'i = ' + str(i)
        new_ancestors_by_state = {}

        # Calculate recombination rate
        expected_recombs = 1.
        if use_recomb_rates and len(recomb_rate_dict[SNPs[i][0]]) > 0:
            expected_recombs, recomb_index = calc_recomb_rate(int(SNPs[i-1][1]), int(SNPs[i][1]), recomb_index,
                                                              recomb_rate_dict[SNPs[i][0]], effective_pop,
                                                              num_generations)
            expected_recombs *= recomb_adjustment

        # At every SNP, find probabilities for each state
        for curr_state in states:
            # For each state, the emission probability is either emit_p[state] or emit_p[~state]
            emit_key = get_emit_key(curr_state, SNPs[i][3], desc_strain)

            # The probability of a given state for a given SNP is the maximum out of ((prob prev_state) *
            #  (prob trans prev_state -> curr_state) * (prob emit curr_state)) for all previous states
            state_probabilities = []
            for prev_state in states:
                state_probabilities.append((prob_nodes[i-1][prev_state] + trans_p[prev_state][curr_state] +
                                            emit_p[curr_state][emit_key], prev_state))

            (recomb_prob, orig_prob, prev_state) = max([(prob + math.log(expected_recombs), prob, prev_s)
                                                        if curr_state != prev_s else (prob, prob, prev_s)
                                                        for prob, prev_s in state_probabilities])

            # Keep track of probabilities in prob_nodes and ancestors in new_ancestors_by_state for each curr_state
            prob_nodes[i][curr_state] = orig_prob
            new_ancestors_by_state[curr_state] = ancestors_by_state[prev_state] + [curr_state]

        # Update ancestors with additional iteration
        ancestors_by_state = new_ancestors_by_state

    # Find maximum-likelihood ancestors for the current chromosome
    (prob, best_state) = max((prob_nodes[len(SNPs)-1][s], s) for s in states)
    return ancestors_by_state[best_state]


# Expectation-Maximization loop
def expectation_maximization(states, trans_in_p, emit_same_p, adjust_recomb, unk_cutoff, use_recomb_rates,
                             desc_ins_by_chr, desc_del_by_chr, anc_ins_by_chr, anc_del_by_chr, args,
                             output_dir, append_str, job_server, vit_func):
    # Start, transition, and emission probabilities
    trans_p = {s_outer: {s_inner: log(trans_in_p) if s_inner == s_outer else log((1 - trans_in_p) / 6)
                         for s_inner in states} for s_outer in states}
    emit_p = {s: {s: log(emit_same_p), '~'+s: log(1 - emit_same_p)} for s in states}

    # EM loop
    tot_prob_dist = 10.
    run_count = 0
    while tot_prob_dist > PROB_DIST_CUTOFF and run_count < int(args.max_iter):
        tot_prob_dist = 0.

        print '---- Run %i ----' % run_count
        sys.stdout.flush()

        # Run viterbi to find maximum likelihood path
        ancestors_by_chr = {}
        if args.parallel:
            jobs = []
            for curr_chr, SNPs in SNPs_by_chr.items():
                job = vit_func.submit(SNPs, states, trans_p, emit_p, args.desc_strain, recomb_rate_dict, args.effective_pop,
                                      args.num_generations, adjust_recomb, use_recomb_rates, args.verbose)
                jobs.append((curr_chr, job))

            for curr_chr, job in jobs:
                ancestors_by_chr[curr_chr] = job()
            job_server.wait()

            if args.verbose:
                job_server.print_stats()
        else:
            for curr_chr, SNPs in SNPs_by_chr.items():
                ancestors_by_chr[curr_chr] = viterbi(SNPs, states, trans_p, emit_p, args.desc_strain, recomb_rate_dict,
                                                     args.effective_pop, args.num_generations, adjust_recomb,
                                                     use_recomb_rates, args.verbose)

        # Recalculate transition and emission probabilities
        new_trans_p = calc_new_trans_p(ancestors_by_chr, states)
        tot_prob_dist += prob_dist(trans_p, new_trans_p)

        new_emit_p = calc_new_emit_p(ancestors_by_chr, SNPs_by_chr, states, args.desc_strain, emit_same_p)
        tot_prob_dist += prob_dist(emit_p, new_emit_p)

        if args.verbose:
            print 'Transition probabilities'
            print '  Before:'
            tuples_dict = prob_tuples(trans_p)
            for s in states:
                print '    %s ->' % s,
                for tup in tuples_dict[s]:
                    print '%s: %.2f,' % tup,
                print ''
            print '  After:'
            tuples_dict = prob_tuples(new_trans_p)
            for s in states:
                print '    %s ->' % s,
                for tup in tuples_dict[s]:
                    print '%s: %.2f,' % tup,
                print ''

            print '\nEmission probabilities'
            print '  Before:'
            tuples_dict = prob_tuples(emit_p)
            for s in states:
                print '    %s ->' % s,
                for tup in tuples_dict[s]:
                    print '%s: %.2f,' % tup,
                print ''
            print '  After:'
            tuples_dict = prob_tuples(new_emit_p)
            for s in states:
                print '    %s ->' % s,
                for tup in tuples_dict[s]:
                    print '%s: %.2f,' % tup,
                print ''

            print 'Total probability distance: %.3f' % tot_prob_dist

        # If set, write out results at each iteration
        if args.write_iter:
            # Score results
            final_score = -1
            if len(desc_ins_by_chr) > 0:
                hits_by_chr, misses_by_chr, all_scores_by_chr = score_results(ancestors_by_chr, SNPs_by_chr,
                                                                              desc_ins_by_chr, desc_del_by_chr,
                                                                              anc_ins_by_chr, anc_del_by_chr)
                try:
                    final_score = sum(hits_by_chr.values())/float(sum(misses_by_chr.values()))
                except ZeroDivisionError:
                    final_score = -1

                # Output scores
                write_scores(output_dir, filename_in, append_str, all_scores_by_chr)

            # Output results (results during EM iterations will not have reclassified IBD and Unk)
            write_ancestors(output_dir, filename_in, append_str, ancestors_by_chr, SNPs_by_chr, STATE_RGBS)
            write_statistics(output_dir, filename_in, append_str, ancestors_by_chr, SNPs_by_chr,
                             (use_recomb_rates, trans_in_p, emit_same_p, trans_p, emit_p), final_score,
                             run_count, time() - time_start, tot_prob_dist)

        trans_p = new_trans_p
        emit_p = new_emit_p

        run_count += 1

    # Reclassify segments where segment likely came from an unsequenced ancestor or where SNP counts for multiple
    #  ancestors are identical
    ancestors_by_chr = reclassify_ibd_and_unk(ancestors_by_chr, SNPs_by_chr, args.desc_strain, unk_cutoff)

    print '\nTotal time (min): ' + str((time() - time_start)/60)
    print 'Total runs: ' + str(run_count)

    # Score results
    final_score = -1
    if len(desc_ins_by_chr) > 0:
        hits_by_chr, misses_by_chr, all_scores_by_chr = score_results(ancestors_by_chr, SNPs_by_chr, desc_ins_by_chr,
                                                                      desc_del_by_chr, anc_ins_by_chr, anc_del_by_chr)
        try:
            final_score = sum(hits_by_chr.values())/float(sum(misses_by_chr.values()))
        except ZeroDivisionError:
            final_score = -1

        print '\nChromosome scores:'
        for curr_chr in hits_by_chr.keys():
            try:
                chr_score = hits_by_chr[curr_chr]/float(misses_by_chr[curr_chr])
            except ZeroDivisionError:
                chr_score = -1
            print ' %s - Hits: %i, Misses: %i, Ratio: %.3f' % (curr_chr, hits_by_chr[curr_chr], misses_by_chr[curr_chr],
                                                               chr_score)
        print 'Final score: %.3f' % (final_score)

        # Output scores
        write_scores(output_dir, filename_in, append_str, all_scores_by_chr)

    # Output final results after reclassifying IBD and Unk
    write_ancestors(output_dir, filename_in, append_str, ancestors_by_chr, SNPs_by_chr, STATE_RGBS)
    write_statistics(output_dir, filename_in, append_str, ancestors_by_chr, SNPs_by_chr, (use_recomb_rates,
                     trans_in_p, emit_same_p, trans_p, emit_p), final_score, run_count, time() - time_start,
                     tot_prob_dist)


# Main method
if __name__ == '__main__':
    # Start timer
    time_start = time()

    # Constants
    MAX_EMIT_SAME_RATE = .99  # maximum allowed emit-same rate
    PROB_DIST_CUTOFF = .01  # prob dist threshold for ending EM loop

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
        job_server = pp.Server()
        vit_func = pp.Template(job_server, viterbi, depfuncs=(calc_recomb_rate, get_emit_key),
                               modules=('numpy', 'math'))
    
    for trans_in_p in trans_in_p_range:
        for emit_same_p in emit_same_p_range:
            for adjust_recomb in adjust_recomb_range:
                for unk_cutoff in unk_cutoff_range:
                    append_str = ''
                    if args.append_params:
                        append_str += '_%.2ft-%.2fe-%.2fu-%.2fa' % (trans_in_p, emit_same_p, adjust_recomb, unk_cutoff)
                    if args.append_date:
                        append_str += datetime.now().strftime('_%y-%m-%d_%H-%M')
                    expectation_maximization(states, trans_in_p, emit_same_p, adjust_recomb, unk_cutoff, use_recomb_rates,
                                             desc_ins_by_chr, desc_del_by_chr, anc_ins_by_chr, anc_del_by_chr,
                                             args, output_dir, append_str, job_server, vit_func)

    if args.parallel:
        job_server.wait()
