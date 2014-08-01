
import argparse
import pp
import sys

from math import log
from numpy import zeros
from os import environ
from time import time

from hmm_output import write_ancestors_to_file, write_confidence_interval, \
                       write_indentical_by_ancestor, write_statistics
from hmm_prob import calc_confidence_intervals, calc_new_emit_p, calc_new_trans_p, label_identical_ancestors, prob_dist
from hmm_util import calc_recomb_rate, count_hotspots, log_add_pair, read_hotspots_data, read_recomb_rates_data, \
                     read_SNP_data

#-----------
# Arguments
#-----------
def read_arguments():
    args = argparse.ArgumentParser(add_help=False)
    args.add_argument('-i', '--input-file', help='Input SNP data file', required=True)

    args.add_argument('-t', '--trans-in-p', help='Starting trans-in probability', default=.94)
    args.add_argument('-e', '--emit-same-p', help='Starting emit-same probability', default=.99)
    args.add_argument('-m', '--max-iter', help='Maximum number of EM iterations', default=50)

    args.add_argument('-r', '--use-recomb-rates', help='Use recombination rates', default=True)
    args.add_argument('-s', '--use-snp-dist', help='Use SNP distance', default=False)
    args.add_argument('-h', '--use-hotspots', help='Use hotspots', default=False)
    args.add_argument('-f', '--hotspot-fi', help='Fold increase per hotspot', default=40)

    args.add_argument('-d', '--append-date', help='Append date to output filename', default=True)
    args.add_argument('-o', '--append-str', help='Append string to output filename', default='')
    args.add_argument('-v', '--verbose', help='Verbose', default=False)

    return args.parse_args()


#-------------------
# Viterbi algorithm
#-------------------
def viterbi(SNPs, states, trans_p, emit_p, input_strain, fi_per_hotspot, hotspot_dict, recomb_rate_dict,
            effective_pop, num_generations, use_hotspots, use_SNP_dist, use_recomb_rates, verbose):
    # Initialize
    prob_nodes = zeros(len(SNPs), dtype={'names': states, 'formats': ['f8']*len(states)})
    ancestors_by_state = {}
    recomb_index = None # Initially set to None so calc_recomb_rate uses a special case when called for the first time

    # Start probabilities
    for s in states:
        # For each state, the emission probability is either emit_p[state] or emit_p[~state]
        emit_key = s
        if s == 'Unk':
            if SNPs[0][3] != input_strain:
                emit_key = '~' + s
        elif s not in SNPs[0][3].split('_'):
            emit_key = '~' + s
        # Probability of a given state at SNPs[0] is emit prob of stat
        prob_nodes[0][s] = emit_p[s][emit_key]
        ancestors_by_state[s] = [s]

    # Rest of probabilities
    for i in range(1, len(SNPs)):
        if verbose and i % 1000 == 0:
            print 'i = ' + str(i)
        new_ancestors_by_state = {}

        # Use of additional data flags
        hotspots_count = 0
        if use_hotspots:
            hotspots_count = count_hotspots(SNPs[i][0], int(SNPs[i-1][1]), int(SNPs[i][1]), hotspot_dict)

        SNP_dist = 1
        if use_SNP_dist:
            SNP_dist = max((int(SNPs[i][1]) - int(SNPs[i-1][1])) / 100, 1)

        expected_recombs = 1.
        if use_recomb_rates and len(recomb_rate_dict[SNPs[i][0]]) > 0:
            expected_recombs, recomb_index = calc_recomb_rate(int(SNPs[i-1][1]), int(SNPs[i][1]), recomb_index,
                                                              recomb_rate_dict[SNPs[i][0]], effective_pop,
                                                              num_generations)

        # At every SNP, find probabilities for each state
        for curr_state in states:
            # For each state, the emission probability is either emit_p[state] or emit_p[~state]
            emit_key = curr_state
            if curr_state == 'Unk':
                if SNPs[i][3] != input_strain:
                    emit_key = '~' + curr_state
            elif curr_state not in SNPs[i][3].split('_'):
                emit_key = '~' + curr_state

            # The probability of a given state for a given SNP is the maximum out of ((prob prev_state) *
            #  (prob trans prev_state -> curr_state) * (prob emit curr_state)) for all previous states
            state_probabilities = []
            for prev_state in states:
                curr_hotspot_fi = 1
                if prev_state == curr_state:
                    curr_trans_p = trans_p[prev_state][curr_state] * SNP_dist
                else:
                    curr_trans_p = 0.
                    for j in range(1,SNP_dist):
                        # hopefully I can eventually remove this check to speed things up (should always be true)
                        if curr_trans_p < trans_p[prev_state][prev_state]*j:
                            raise Exception('log_add: curr_trans_p < trans_p[prev_state][prev_state]*j, need to add check')
                        curr_trans_p = log_add_pair(curr_trans_p, trans_p[prev_state][prev_state]*j)
                    curr_trans_p += trans_p[prev_state][curr_state]

                    # Only apply hotspot fold increase to transition probabilities from one state to a different state
                    curr_hotspot_fi = max(fi_per_hotspot * hotspots_count, 1)

                # Calculate probabilities
                state_probabilities.append((prob_nodes[i-1][prev_state] + curr_trans_p + log(curr_hotspot_fi) +
                                            emit_p[curr_state][emit_key], prev_state))

            (recomb_prob, orig_prob, prev_state) = max([(prob + log(expected_recombs), prob, prev_s)
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

#-------------
# Main method
#-------------
if __name__ == "__main__":
    # Start timer
    time_start = time()

    # Constants
    WORKING_DIR = '../'
    try:
        if environ['HOSTNAME'][-6:] == '.local':
            WORKING_DIR = '/Users/dama9282/AncestorInference/'
    except KeyError:
        pass

    STATES = ('Unk', 'A', 'ARK', 'BALBc', 'C3HHe', 'C57BL6N', 'DBA2')
    STATE_RGBS = {'Unk': '128,128,128', 'A': '0,153,0', 'ARK': '51,51,255', 'BALBc': '255,255,51',
                  'C3HHe': '255,153,51', 'C57BL6N': '102,0,204', 'DBA2': '255,51,51', 'IBA': '153,255,255'}

    EFFECTIVE_POP = 1  # effective population (N_e) for recombination rate calculations
    NUM_GENERATIONS = 25  # number of generations between ancestors and ILS/ISS strains

    # Read in arguments
    args = read_arguments()
    input_strain = args.input_file.strip().rsplit('/',1)[1].split('_')[0] \
                   if args.input_file.strip().rsplit('/',1)[1].split('_')[0] != 'TEST' else 'ISS' #ILS or ISS

    # Read in data
    SNPs_by_chr = read_SNP_data(args.input_file)
    hotspot_dict = {}
    if args.use_hotspots:
        hotspot_dict = read_hotspots_data(WORKING_DIR + 'data/mouse_hotspots.csv')
    recomb_rate_dict = {}
    if args.use_recomb_rates:
        recomb_rate_dict = read_recomb_rates_data(WORKING_DIR + 'data/mouse_recomb_rates.csv')

    # Start, transition, and emission probabilities
    trans_p = {s_outer: {s_inner: log(args.trans_in_p) if s_inner == s_outer else log((1 - args.trans_in_p) / 6) \
                         for s_inner in STATES} for s_outer in STATES}
    emit_p = {s: {s: log(args.emit_same_p), '~'+s: log(1 - args.emit_same_p)} for s in STATES}

    # Set up parallel python server for viterbi function
    job_server = pp.Server()
    vit_func = pp.Template(job_server, viterbi, depfuncs=(count_hotspots, calc_recomb_rate, log_add_pair),
                           modules=('from numpy import zeros', 'from math import e, log'))

    # Run algorithms
    tot_prob_dist = 10.
    run_count = 0
    while tot_prob_dist > 0.01 and run_count < int(args.max_iter):
        tot_prob_dist = 0.

        print '---- Run %i ----' % run_count
        sys.stdout.flush()

        jobs = []
        for curr_chr, SNPs in SNPs_by_chr.items():
            # Run viterbi to find maximum likelihood path
            job = vit_func.submit(SNPs, STATES, trans_p, emit_p, input_strain, args.hotspot_fi,
                                                        hotspot_dict, recomb_rate_dict, EFFECTIVE_POP, NUM_GENERATIONS,
                                                        args.use_hotspots, args.use_snp_dist, args.use_recomb_rates,
                                                        args.verbose)
            jobs.append((curr_chr, job))

        ancestors_by_chr = {}
        for curr_chr, job in jobs:
            ancestors_by_chr[curr_chr] = job()

        print 'before wait:'
        print job_server.print_stats()
        job_server.wait()
        print 'after wait:'
        print job_server.print_stats()

        # Label segments where SNP counts for multiple ancestors are identical
        ancestors_by_chr = label_identical_ancestors(ancestors_by_chr, SNPs_by_chr, input_strain)

        # Recalculate transition and emission probabilities
        new_trans_p = calc_new_trans_p(ancestors_by_chr, STATES)
        tot_prob_dist += prob_dist(trans_p, new_trans_p)

        new_emit_p = calc_new_emit_p(ancestors_by_chr, SNPs_by_chr, STATES, input_strain, args.emit_same_p)
        tot_prob_dist += prob_dist(emit_p, new_emit_p)

        if args.verbose:
            print 'Transition probabilities'
            print '  Before: ' + str(trans_p)
            print '  After:  ' + str(new_trans_p)

            print 'Emission probabilities'
            print '  Before: ' + str(emit_p)
            print '  After:  ' + str(new_emit_p)

            print 'Total probability distance: %.3f' % tot_prob_dist

        trans_p = new_trans_p
        emit_p = new_emit_p

        run_count += 1

    #confidence_intervals = calc_confidence_intervals(ancestors, SNPs)

    print 'Total time (min): ' + str((time() - time_start)/60)
    print 'Total runs: ' + str(run_count)

    filename_in = args.input_file.rsplit('/', 1)[1]

    #write_indentical_by_ancestor(WORKING_DIR, filename_in, unique_output_name, identical_by_anc)
    #write_confidence_interval(WORKING_DIR, filename_in, unique_output_name, confidence_intervals)

    write_ancestors_to_file(WORKING_DIR, filename_in, args.append_str, ancestors_by_chr, SNPs_by_chr, STATE_RGBS)
    write_statistics(WORKING_DIR, filename_in, args.append_str, ancestors_by_chr, SNPs_by_chr, (args.use_recomb_rates,
                     args.trans_in_p, args.emit_same_p, trans_p, emit_p), run_count, time() - time_start)
