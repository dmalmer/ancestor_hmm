
import sys

from itertools import tee, izip
from numpy import array, loadtxt, zeros
from math import *
from collections import defaultdict


#----------------
# Helper methods
#----------------
# Return a list of pairwise elements (taken from https://docs.python.org/2/library/itertools.html#recipes)
def pairwise(iterable):
    # s -> (s0,s1), (s1,s2), (s2, s3), ...
    a, b = tee(iterable)
    next(b, None)
    return izip(a, b)


# Count the number of hotspots between two chromosome positions
def count_hotspots(chromosome, pos_start, pos_end, hotspot_dict):
    # All even indexes are the pos start of cold (not hot) spots, all odd indexes are the pos start of hot spots
    i_start = len(hotspot_dict[chromosome][pos_start-hotspot_dict[chromosome] > 0]) - 1
    i_end = len(hotspot_dict[chromosome][pos_end-hotspot_dict[chromosome] > 0]) - 1

    hs_count = ceil((i_end-i_start)/2.)

    # if both the start and end indexes are on hotspots, we need to add 1 more to the count
    if i_start % 2 == 1 and i_end % 2 == 1:
        hs_count += 1

    return hs_count


# Generator to loop over only unique ancestors
def unique_ancestors(ancestors, SNPs):
    curr_i = 0
    while curr_i < len(ancestors):
        start_i = curr_i
        curr_i += 1
        while curr_i < len(ancestors) and ancestors[curr_i-1] == ancestors[curr_i]:
            curr_i += 1

        yield (SNPs[start_i,0], SNPs[start_i,1], SNPs[curr_i-1,2], ancestors[start_i])



#----------------
# Output methods
#----------------
def print_ancestors(ancestors, SNPs, title):
    print ''
    print title
    cur = ''
    for i in range(len(ancestors)):
        if cur != ancestors[i]:
            print str(SNPs[i][1]) + ': ' + ancestors[i]
            cur = ancestors[i]


def write_ancestors_to_file(ancestors, SNPs):
    extension = sys.argv[1].rsplit('.', 1)[1]
    out_file = open(sys.argv[1].split('.' + extension)[0] + '_hmm-out' + unique_output_name + '.' + extension, 'w')

    out_len = 0
    for chromosome, pos_start, pos_end, ancestor in unique_ancestors(ancestors, SNPs):
        out_file.write(chromosome + '\t' + pos_start + '\t' + pos_end + '\t' + ancestor + '\n')
        out_len += 1
    out_file.close()

    print '\nOutput file length: ' + str(out_len)
    print 'Percentage change: ' + str(float(len(SNPs)-out_len)/float(len(SNPs)))


#---------------------
# Probability methods
#---------------------
def calc_new_trans_p_and_hs_fi(ancestors, SNPs, states, hotspot_dict):
    tot_hotspots = 0
    tot_hs_trans = 0
    tot_hs_not_trans = 0
    #set minimum transition counts to 1 so we don't have any 0 probabilities
    trans_counts = {s_outer: {s_inner: 1 for s_inner in states} for s_outer in states}
    for ((anc_prev, anc_curr), (pos_prev, pos_curr), (chr_prev, chr_curr)) in \
            zip(pairwise(ancestors), pairwise(SNPs[:,1]), pairwise(SNPs[:,0])):
        if chr_prev == chr_curr:
            hotspots_count = count_hotspots(chr_curr, int(pos_prev), int(pos_curr), hotspot_dict)
            if hotspots_count > 0:
                tot_hotspots += hotspots_count
                if anc_prev != anc_curr:
                    tot_hs_trans += 1
                else:
                    tot_hs_not_trans += 1

            trans_counts[anc_prev][anc_curr] += 1

    print 'tot_hotspots:'
    print tot_hotspots
    print 'tot_hs_trans:'
    print tot_hs_trans
    print 'tot_hs_not_trans:'
    print tot_hs_not_trans

    new_trans_p = {}
    for s_outer in states:
        new_trans_p[s_outer] = {}
        for s_inner in states:
            tot_trans = float(sum(trans_counts[s_outer].values()))
            new_trans_p[s_outer][s_inner] = log(trans_counts[s_outer][s_inner] / tot_trans)

    new_hs_fi = max(float(tot_hs_trans)/max(tot_hs_not_trans, 1), 1.)

    return new_trans_p, new_hs_fi


def calc_new_emit_p(ancestors, SNPs, states, input_group, def_emit_same_p, def_emit_other_p):
    # SNP_counts[<state>] = total number of <state> appearances
    # SNP_counts[<~state>] = total number of <~state> appearances
    SNP_counts = defaultdict(int)
    # ancestor_counts[<state>] = total number of calculated <state> when <state> is observed
    # ancestor_counts[<~state>] = total number of calculated <state> when <state> is NOT observed
    ancestor_counts = defaultdict(int)

    for i in range(len(SNPs)):
        # SNP_counts
        for s in states:
            SNP_key = s
            if s == 'Unk':
                if SNPs[i][3] != input_group:
                    SNP_key = '~' + s
            elif s not in SNPs[i][3].split('_'):
                SNP_key = '~' + s
            SNP_counts[SNP_key] += 1

        # ancestor_counts
        if ancestors[i] in SNPs[i][3].split('_'):
            ancestor_counts[ancestors[i]] += 1
        else:
            ancestor_counts['~'+ancestors[i]] += 1

    new_emit_p = {}
    for s in states:
        new_emit_p[s] = {}
        # if we don't observe or calculate a state, use the default probabilities
        if SNP_counts[s] == 0 or ancestor_counts[s] == 0:
            new_emit_p[s][s] = log(def_emit_same_p)
            new_emit_p[s]['~'+s] = log(def_emit_other_p)
        else:
            # normalize <state> and <~state> probabilities to one
            state_p = float(ancestor_counts[s]) / SNP_counts[s]
            notstate_p = float(ancestor_counts['~'+s]) / SNP_counts['~'+s]
            normalizer = 1 / (state_p + notstate_p)
            # don't allow probabilities of 1.0 and 0.0
            new_emit_p[s][s] = log(min(normalizer * state_p, .99))
            new_emit_p[s]['~'+s] = log(max(normalizer * notstate_p, .01))

    return new_emit_p


def prob_dist(old_probs, new_probs):
    tot_dist = 0.
    tot_probs = 0
    for key in old_probs.keys():
        for old_p, new_p in izip(old_probs[key].values(), new_probs[key].values()):
            tot_dist += abs(old_p - new_p)
            tot_probs += 1

    # return average prob dist
    return tot_dist / tot_probs


#-------------------
# Viterbi algorithm
#-------------------
def viterbi(SNPs, states, start_p, trans_p, emit_p, fold_increase_per_hotspot, hotspot_dict, input_group):
    # initialize
    prob_nodes = zeros(len(SNPs), dtype={'names': states, 'formats': ['f8']*len(states)})
    ancestors_by_state = {}

    # SNPs[0] probabilities
    for s in states:
        # for each state, the emission probability is either emit_p[state] or emit_p[~state]
        emit_key = s
        if s == 'Unk':
            if SNPs[0][3] != input_group:
                emit_key = '~' + s
        elif s not in SNPs[0][3].split('_'):
            emit_key = '~' + s
        # probability of a given state at SNPs[0] is the (start prob of state) * (emit prob of state)
        prob_nodes[0][s] = start_p[s] + emit_p[s][emit_key]
        ancestors_by_state[s] = [s]

    # rest of SNP probabilities
    for i in range(1, len(SNPs)):
        if i % 10000 == 0:
            print 'i = ' + str(i)
        new_ancestors_by_state = {}

        hotspots_count = count_hotspots(SNPs[i][0], int(SNPs[i-1][1]), int(SNPs[i][1]), hotspot_dict)

        # at every SNP, find probabilities for each state
        for curr_state in states:
            # for each state, the emission probability is either emit_p[state] or emit_p[~state]
            emit_key = curr_state
            if curr_state == 'Unk':
                if SNPs[i][3] != input_group:
                    emit_key = '~' + curr_state
            elif curr_state not in SNPs[i][3].split('_'):
                emit_key = '~' + curr_state

            # the probability of a given state for a given SNP is the maximum
            #  out of (prob prev_state) * (prob trans prev_state -> curr_state) * (prob emit curr_state)
            #  for all previous states
            state_probabilities = []
            for prev_state in states:
                # incorporate recombination hotspot data
                hotspot_fi = 1
                if prev_state != curr_state:
                    # only apply hotspot fold increase to transition probabilities from one state to a different state
                    hotspot_fi = max(fold_increase_per_hotspot * hotspots_count, 1)

                state_probabilities.append((prob_nodes[i-1][prev_state] + trans_p[prev_state][curr_state] +
                                            log(hotspot_fi) + emit_p[curr_state][emit_key], prev_state))

            (prob, prev_state) = max(state_probabilities)

            # keep track of probabilities in prob_nodes and ancestors in new_ancestors_by_state for each curr_state
            prob_nodes[i][curr_state] = prob
            new_ancestors_by_state[curr_state] = ancestors_by_state[prev_state] + [curr_state]

        # update ancestors with additional iteration
        ancestors_by_state = new_ancestors_by_state

    # find maximum-likelihood ancestors
    (prob, best_state) = max((prob_nodes[len(SNPs)-1][s], s) for s in states)

    return ancestors_by_state[best_state]


#-------------
# Main method
#-------------
if __name__ == "__main__":
    # Input/output names
    input_group = sys.argv[1].strip().rsplit('/',1)[1].split('_')[0] if sys.argv[1].strip().split('/')[1].split('_')[0] != 'TEST' else 'ISS' #ILS or ISS
    unique_output_name = '-NEWPROB_testtest'

    # Starting probabilities
    #  --these need to be translated into log space
    def_start_p = 1/.7
    def_trans_in_p = .46
    def_trans_out_p = .09
    def_emit_same_p = .95
    def_emit_other_p = .05
    def_fold_increase_per_hotspot = 40

    # Read in SNP data
    SNPs = loadtxt(sys.argv[1], dtype='string')
    print 'Input file length: ' + str(len(SNPs))

    # Read in hotspot data
    #  defaultdict list that begins with a 0 element rather than starting empty
    #  all even indexes are the pos start of cold (not hot) spots, all odd indexes are the pos start of hot spots
    hotspot_dict = defaultdict(lambda: [0])
    with open('data/mouse_hotspots.csv', 'r') as f:
        f.readline()
        line = f.readline()
        while line != '':
            splits = line.split(',')
            hotspot_dict['chr'+str(splits[0])].extend([int(splits[1]), int(splits[2])+1])
            line = f.readline()
    #  convert to numpy arrays
    for k, v in hotspot_dict.items():
        hotspot_dict[k] = array(v)

    #print hotspot_dict['chr1'][:9]
    #print ['cold ', 'hot ', 'cold', 'hot ', 'cold', 'hot ', 'cold', 'hot ', 'cold']
    #print hotspots_count('chr1', 3020000, 3095000, hotspot_dict)

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
    i = 0
    while tot_prob_dist > 0.01 and i < 2:
        tot_prob_dist = 0.
        print '\n---- RUN %i ----' % i
        ancestors = viterbi(SNPs, states, start_p, trans_p, emit_p, fold_increase_per_hotspot, hotspot_dict, input_group)

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

        i += 1

    print_ancestors(ancestors, SNPs, 'Viterbi')
    write_ancestors_to_file(ancestors, SNPs)
