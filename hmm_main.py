
import sys

from collections import defaultdict
from math import log
from numpy import array, loadtxt, zeros

from hmm_output import print_ancestors, write_ancestors_to_file, write_statistics
from hmm_prob import calc_new_emit_p, calc_new_trans_p_and_hs_fi, prob_dist
from hmm_util import count_hotspots, log_add


#-------------------
# Viterbi algorithm
#-------------------
def viterbi(SNPs, states, start_p, trans_p, emit_p, fold_increase_per_hotspot, hotspot_dict, input_group,
            use_hotspots, use_SNP_dist):
    # Initialize
    prob_nodes = zeros(len(SNPs), dtype={'names': states, 'formats': ['f8']*len(states)})
    ancestors_by_state = {}

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
        if i % 1000 == 0:
            print 'i = ' + str(i)
        new_ancestors_by_state = {}

        hotspots_count = count_hotspots(SNPs[i][0], int(SNPs[i-1][1]), int(SNPs[i][1]), hotspot_dict)
        SNP_dist = max((int(SNPs[i][1]) - int(SNPs[i-1][1])) / 100, 1)

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
                # If set, incorporate recombination hotspot data and SNP distance
                hotspot_fi = 1
                if prev_state == curr_state:
                    curr_trans_p = trans_p[prev_state][curr_state] * SNP_dist if \
                        use_SNP_dist else trans_p[prev_state][curr_state]
                else:
                    if use_SNP_dist:
                        curr_trans_p = 0.
                        for j in range(1,SNP_dist):
                            # hopefully I can eventually remove this check to speed things up (should always be true)
                            if curr_trans_p < trans_p[prev_state][prev_state]*j:
                                raise Exception('log_add: curr_trans_p < trans_p[prev_state][prev_state]*j, need to add check')
                            curr_trans_p = log_add(curr_trans_p, trans_p[prev_state][prev_state]*j)
                        curr_trans_p += trans_p[prev_state][curr_state]
                    else:
                        curr_trans_p = trans_p[prev_state][curr_state]

                    # Only apply hotspot fold increase to transition probabilities from one state to a different state
                    if use_hotspots:
                        hotspot_fi = max(fold_increase_per_hotspot * hotspots_count, 1)

                # Calculate probabilities
                state_probabilities.append((prob_nodes[i-1][prev_state] + curr_trans_p + log(hotspot_fi) +
                                            emit_p[curr_state][emit_key], prev_state))

            (prob, prev_state) = max(state_probabilities)

            # Keep track of probabilities in prob_nodes and ancestors in new_ancestors_by_state for each curr_state
            prob_nodes[i][curr_state] = prob
            new_ancestors_by_state[curr_state] = ancestors_by_state[prev_state] + [curr_state]

        # Update ancestors with additional iteration
        ancestors_by_state = new_ancestors_by_state

    # Find maximum-likelihood ancestors
    (prob, best_state) = max((prob_nodes[len(SNPs)-1][s], s) for s in states)

    return ancestors_by_state[best_state]


#-------------
# Main method
#-------------
if __name__ == "__main__":
    # Input/output names
    input_group = sys.argv[1].strip().rsplit('/',1)[1].split('_')[0] if \
        sys.argv[1].strip().split('/')[1].split('_')[0] != 'TEST' else 'ISS' #ILS or ISS
    unique_output_name = ''

    # Starting settings
    #  --The probabilities need to be translated into log space
    def_start_p = 1/.7
    def_trans_in_p = .9999
    def_trans_out_p = (1 - def_trans_in_p) / 6
    def_emit_same_p = .95
    def_emit_other_p = 1 - def_emit_same_p
    def_fold_increase_per_hotspot = 1
    use_hotspots = True
    use_SNP_dist = True

    # Read in SNP data
    SNPs = loadtxt(sys.argv[1], dtype='string')
    print 'Input file length: ' + str(len(SNPs))

    # Read in hotspot data
    #  defaultdict list that begins with a 0 element rather than starting empty
    #  All even indexes are the pos start of cold (not hot) spots, all odd indexes are the pos start of hot spots
    hotspot_dict = defaultdict(lambda: [0])
    with open('data/mouse_hotspots.csv', 'r') as f:
        f.readline()
        line = f.readline()
        while line != '':
            splits = line.split(',')
            hotspot_dict['chr'+str(splits[0])].extend([int(splits[1]), int(splits[2])+1])
            line = f.readline()
    #  Convert to numpy arrays
    for k, v in hotspot_dict.items():
        hotspot_dict[k] = array(v)

    # States
    states = ('Unk', 'A', 'ARK', 'BALBc', 'C3HHe', 'C57BL6N', 'DBA2')

    # Start, transition, and emission probabilities
    #  Equal start probability for each state
    #  Transition to same state = 0.46, to other state = 0.09
    #  Emit same state = 0.95, other state = 0.05 (labeled '~State')
    start_p = {s: log(def_start_p) for s in states}
    trans_p = {s_outer: {s_inner: log(def_trans_in_p) if s_inner == s_outer else log(def_trans_out_p) for s_inner in states} for s_outer in states}
    emit_p = {s: {s: log(def_emit_same_p), '~'+s: log(def_emit_other_p)} for s in states}

    # Hotspot fold increase
    fold_increase_per_hotspot = def_fold_increase_per_hotspot

    # Run algorithms
    tot_prob_dist = 10.
    run_count = 0
    while tot_prob_dist > 0.01 and run_count < 1:
        tot_prob_dist = 0.
        print '\n---- RUN %i ----' % run_count
        ancestors = viterbi(SNPs, states, start_p, trans_p, emit_p, fold_increase_per_hotspot, hotspot_dict,
                            input_group, use_hotspots, use_SNP_dist)

        new_trans_p, new_hs_fi = calc_new_trans_p_and_hs_fi(ancestors, SNPs, states, hotspot_dict)

        print '\nHotspot fold increase'
        print '  Before: ' + str(fold_increase_per_hotspot)
        print '  After:  ' + str(new_hs_fi)
        fold_increase_per_hotspot = new_hs_fi

        print 'Transition probabilities'
        print '  Before: ' + str(trans_p)
        print '  After:  ' + str(new_trans_p)
        tot_prob_dist += prob_dist(trans_p, new_trans_p)
        trans_p = new_trans_p

        new_emit_p = calc_new_emit_p(ancestors, SNPs, states, input_group, def_emit_same_p, def_emit_other_p)

        print 'Emission probabilities'
        print '  Before: ' + str(emit_p)
        print '  After:  ' + str(new_emit_p)
        tot_prob_dist += prob_dist(emit_p, new_emit_p)
        emit_p = new_emit_p

        print 'Total probability distance: %.3f' % tot_prob_dist

        run_count += 1

    print_ancestors(ancestors, SNPs, 'Viterbi')
    write_ancestors_to_file(sys.argv[1], unique_output_name, ancestors, SNPs)

    write_statistics(sys.argv[1], ancestors, SNPs, (def_trans_in_p, def_trans_out_p, def_emit_same_p,
                     def_emit_other_p, def_fold_increase_per_hotspot, use_hotspots, use_SNP_dist), run_count,
                     fold_increase_per_hotspot)

    print trans_p