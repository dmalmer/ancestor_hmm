
import sys

from math import log
from numpy import loadtxt, zeros
from os import environ
from time import time

from hmm_output import print_ancestors, write_ancestors_to_file, write_confidence_interval, \
                       write_indentical_by_ancestor, write_statistics
from hmm_prob import calc_confidence_intervals, calc_identical_by_ancestor, calc_new_emit_p, calc_new_trans_p, prob_dist
from hmm_util import calc_recomb_rate, count_hotspots, log_add_pair, read_hotspots_data, read_recomb_rates_data

WORKING_DIR = '../'
try:
    if environ['HOSTNAME'][-6:] == '.local':
        WORKING_DIR = '/Users/dama9282/AncestorInference/'
except KeyError:
    pass

#-------------------
# Viterbi algorithm
#-------------------
def viterbi(SNPs, states, start_p, trans_p, emit_p, fi_per_hotspot, hotspot_dict, recomb_rate_dict, effective_pop,
            num_generations, input_group, use_hotspots, use_SNP_dist, use_recomb_rates, verbose):
    # Initialize
    prob_nodes = zeros(len(SNPs), dtype={'names': states, 'formats': ['f8']*len(states)})
    ancestors_by_state = {}
    recomb_index = None # Initially set to None so calc_recomb_rate uses a special case when called for the first time

    # SNPs[0] probabilities
    for s in states:
        # For each state, the emission probability is either emit_p[state] or emit_p[~state]
        emit_key = s
        if s == 'Unk':
            if SNPs[0][3] != input_group:
                emit_key = '~' + s
        elif s not in SNPs[0][3].split('_'):
            emit_key = '~' + s
        # Probability of a given state at SNPs[0] is the (start prob of state) * (emit prob of state)
        prob_nodes[0][s] = start_p[s] + emit_p[s][emit_key]
        ancestors_by_state[s] = [s]

    # Rest of SNP probabilities
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
                if SNPs[i][3] != input_group:
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

            # print 'expected recombs:'
            # print expected_recombs
            # print 'state_probs:'
            # print state_probabilities
            # print 'max(state_probs):'
            # print max(state_probabilities)
            # print [(prob + log(expected_recombs), prev_s) if curr_state != prev_s else (prob, prev_s) for prob, prev_s in state_probabilities]
            # print max([(prob + log(expected_recombs), prob, prev_s) if curr_state != prev_s else (prob, prob, prev_s) for prob, prev_s in state_probabilities])
            (recomb_prob, orig_prob, prev_state) = max([(prob + log(expected_recombs), prob, prev_s)
                                                        if curr_state != prev_s else (prob, prob, prev_s)
                                                        for prob, prev_s in state_probabilities])

            # Keep track of probabilities in prob_nodes and ancestors in new_ancestors_by_state for each curr_state
            prob_nodes[i][curr_state] = orig_prob
            new_ancestors_by_state[curr_state] = ancestors_by_state[prev_state] + [curr_state]

        # Update ancestors with additional iteration
        ancestors_by_state = new_ancestors_by_state

    # Find maximum-likelihood ancestors
    (prob, best_state) = max((prob_nodes[len(SNPs)-1][s], s) for s in states)

    return ancestors_by_state[best_state], prob_nodes


#-------------
# Main method
#-------------
if __name__ == "__main__":
    # Start timer
    time_start = time()

    # Input/output names
    input_group = sys.argv[1].strip().rsplit('/',1)[1].split('_')[0] if \
        sys.argv[1].strip().split('/')[1].split('_')[0] != 'TEST' else 'ISS' #ILS or ISS
    unique_output_name = ''

    # Starting settings
    #  The probabilities are translated into log space later on
    def_start_p = 1/.7
    def_trans_in_p = .64
    def_trans_out_p = (1 - def_trans_in_p) / 6
    def_emit_same_p = .95
    def_emit_other_p = 1 - def_emit_same_p

    fi_per_hotspot = 40  # fold increase per hotspot
    effective_pop = 1  # effective population (N_e) for recombination rate calculations
    num_generations = 25  # number of generations between ancestors and ILS/ISS strains
    iba_cutoff = .90 # how close probabilities need to be to be considered "identical by ancestor"
    use_hotspots = False
    use_SNP_dist = False
    use_recomb_rates = True

    verbose = True
    max_run_count = 1

    # Read in SNP data
    SNPs = loadtxt(sys.argv[1], dtype='string')
    print 'Input file length: ' + str(len(SNPs))

    # Read in hotspot data
    hotspot_dict = read_hotspots_data(WORKING_DIR + 'data/mouse_hotspots.csv')

    # Read in recombination rate data
    recomb_rate_dict = {}
    if use_recomb_rates:
        recomb_rate_dict = read_recomb_rates_data(WORKING_DIR + 'data/mouse_recomb_rates.csv')

    # States
    states = ('Unk', 'A', 'ARK', 'BALBc', 'C3HHe', 'C57BL6N', 'DBA2')
    state_RGBs = {'Unk': '128,128,128', 'A': '0,153,0', 'ARK': '51,51,255', 'BALBc': '255,255,51',
                  'C3HHe': '255,153,51', 'C57BL6N': '102,0,204', 'DBA2': '255,51,51'}

    # Start, transition, and emission probabilities
    start_p = {s: log(def_start_p) for s in states}
    trans_p = {s_outer: {s_inner: log(def_trans_in_p) if s_inner == s_outer else log(def_trans_out_p) for s_inner in states} for s_outer in states}
    emit_p = {s: {s: log(def_emit_same_p), '~'+s: log(def_emit_other_p)} for s in states}

    # Run algorithms
    tot_prob_dist = 10.
    run_count = 0
    while tot_prob_dist > 0.01 and run_count < max_run_count:
        tot_prob_dist = 0.

        print '---- Run %i ----' % run_count
        sys.stdout.flush()

        ancestors, prob_nodes = viterbi(SNPs, states, start_p, trans_p, emit_p, fi_per_hotspot, hotspot_dict,
                                        recomb_rate_dict, effective_pop, num_generations, input_group, use_hotspots,
                                        use_SNP_dist, use_recomb_rates, verbose)

        new_trans_p = calc_new_trans_p(ancestors, states)

        tot_prob_dist += prob_dist(trans_p, new_trans_p)
        trans_p = new_trans_p

        new_emit_p = calc_new_emit_p(ancestors, SNPs, states, input_group, def_emit_same_p, def_emit_other_p)

        tot_prob_dist += prob_dist(emit_p, new_emit_p)
        emit_p = new_emit_p

        if verbose:
            print 'Transition probabilities'
            print '  Before: ' + str(trans_p)
            print '  After:  ' + str(new_trans_p)

            print 'Emission probabilities'
            print '  Before: ' + str(emit_p)
            print '  After:  ' + str(new_emit_p)

            print 'Total probability distance: %.3f' % tot_prob_dist

        run_count += 1

    confidence_intervals = calc_confidence_intervals(ancestors, SNPs, prob_nodes)
    identical_by_anc = calc_identical_by_ancestor(ancestors, SNPs, prob_nodes, states, iba_cutoff)

    print 'Total time (min): ' + str((time() - time_start)/60)
    print 'Total runs: ' + str(run_count)

    if verbose:
        print_ancestors(ancestors, SNPs, 'Viterbi')

    filename_in = sys.argv[1].rsplit('/', 1)[1]

    write_indentical_by_ancestor(WORKING_DIR, filename_in, unique_output_name, identical_by_anc)
    write_confidence_interval(WORKING_DIR, filename_in, unique_output_name, confidence_intervals)

    write_ancestors_to_file(WORKING_DIR, filename_in, unique_output_name, ancestors, SNPs, state_RGBs)
    write_statistics(WORKING_DIR, filename_in, unique_output_name, ancestors, SNPs, (def_trans_in_p, def_trans_out_p,
                     def_emit_same_p, def_emit_other_p, fi_per_hotspot, use_hotspots, use_SNP_dist, use_recomb_rates),
                     run_count, time() - time_start)
