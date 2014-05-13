
import sys

from itertools import tee, izip
from numpy import loadtxt, zeros
from math import *
from collections import defaultdict

#-------------------------------------
# Global variables and helper methods
#-------------------------------------
# Default probabilities
#  --these need to be translated into log space
def_start_p = 1/.7
def_trans_in_p = .46
def_trans_out_p = .09
def_emit_same_p = .95
def_emit_other_p = .05

# Return a list of pairwise elements (taken from https://docs.python.org/2/library/itertools.html#recipes)
def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return izip(a, b)

#----------------
# Output methods
#----------------
def print_path(path, title):
    print ''
    print title
    cur = ''
    for i in range(len(path)):
        if cur != path[i]:
            print str(obs[i][1]) + ': ' + path[i]
            cur = path[i]
    print ''

def output_path(path):
    extension = sys.argv[1].rsplit('.', 1)[1]
    out_file = open(sys.argv[1].split('.' + extension)[0] + '_hmm-out' + unique_output_name + '.' + extension,'w')
    cur = ''
    out_len = 0
    for i in range(len(path)):
        if cur != path[i]:
            if i > 0:
                #write ending position and ancestor of prev line
                out_file.write(obs[i-1][2] + '\t' + path[i-1] + '\n')
            #write chromosome and starting position of current line
            out_file.write(obs[i][0] + '\t' + obs[i][1] + '\t')
            cur = path[i]
            out_len += 1
    #write ending position and ancestor of last line
    out_file.write(obs[len(path)-1][2] + '\t' + path[len(path)-1] + '\n')
    out_file.close()

    print 'Output file length: ' + str(out_len)
    print 'Percentage change: ' + str(float(len(obs)-out_len)/float(len(obs)))

#---------------------
# Probability methods
#---------------------
def calc_new_trans_p(path, states):
    #set minimum transition counts to 1 so we don't have any 0 probabilities
    trans_counts = {s_outer: {s_inner: 1 for s_inner in states} for s_outer in states}
    for prev, curr in pairwise(path):
        trans_counts[prev][curr] += 1

    new_trans_p = {}
    for s_outer in states:
        new_trans_p[s_outer] = {}
        for s_inner in states:
            tot_trans = float(sum(trans_counts[s_outer].values()))
            new_trans_p[s_outer][s_inner] = log(trans_counts[s_outer][s_inner] / tot_trans)

    return new_trans_p

def calc_new_emit_p(path, obs, states, input_group):
    # obs_counts[<state>] = total number of <state> appearances
    # obs_counts[<~state>] = total number of <~state> appearances
    obs_counts = defaultdict(int)
    # path_counts[<state>] = total number of calculated <state> when <state> is observed
    # path_counts[<~state>] = total number of calculated <state> when <state> is NOT observed
    path_counts = defaultdict(int)

    for i in range(len(obs)):
        # obs_counts
        for s in states:
            obs_key = s
            if s == 'Unk':
                if obs[i][3] != input_group:
                    obs_key = '~' + s
            elif s not in obs[i][3].split('_'):
                obs_key = '~' + s
            obs_counts[obs_key] += 1

        # path_counts
        if path[i] in obs[i][3].split('_'):
            path_counts[path[i]] += 1
        else:
            path_counts['~'+path[i]] += 1

    new_emit_p = {}
    for s in states:
        new_emit_p[s] = {}
        # if we don't observe or calculate a state, use the default probabilities
        if obs_counts[s] == 0 or path_counts[s] == 0:
            new_emit_p[s][s] = log(def_emit_same_p)
            new_emit_p[s]['~'+s] = log(def_emit_other_p)
        else:
            # normalize <state> and <~state> probabilities to one
            state_p = float(path_counts[s]) / obs_counts[s]
            notstate_p = float(path_counts['~'+s]) / obs_counts['~'+s]
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
def viterbi(obs, states, start_p, trans_p, emit_p, input_group):
    # intialize
    V = zeros(len(obs), dtype={'names':states, 'formats':['f8']*len(states)})
    path = {}

    # obs[0] probabilities
    for s in states:
        emit_key = s
        if s == 'Unk':
            if obs[0][3] != input_group:
                emit_key = '~' + s
        elif s not in obs[0][3].split('_'):
            emit_key = '~' + s
        V[0][s] = start_p[s] + emit_p[s][emit_key]
        path[s] = [s]

    # run viterbi
    for i in range(1, len(obs)):
        if i % 10000 == 0:
            print 'i = ' + str(i)
        new_path = {}

        # at every observation, find probabilities for each state
        for curr_state in states:
            # for each state, the emission probability is either emit_p[state] or emit_p[~state]
            emit_key = curr_state
            if curr_state == 'Unk':
                if obs[i][3] != input_group:
                    emit_key = '~' + curr_state
            elif curr_state not in obs[i][3].split('_'):
                emit_key = '~' + curr_state

            # the probability of a given state for a given observation is the maximum
            #  out of (prob prev state) * (prob trans prev state -> cur state) * (prob emit cur state)
            #  for all previous states
            (prob, prev_state) = max((V[i-1][prev_state] + trans_p[prev_state][curr_state] + emit_p[curr_state][emit_key], prev_state) for prev_state in states)

            # keep track of probabilities in V and paths in new_path for each currState
            V[i][curr_state] = prob
            new_path[curr_state] = path[prev_state] + [curr_state]

        # update path with additional iteration
        path = new_path

    # find maximum-likelihood path
    (prob, state) = max((V[i][s], s) for s in states)

    return path[state], V

#-------------
# Main method
#-------------
if __name__ == "__main__":
    # Global variables
    input_group = sys.argv[1].strip().rsplit('/',1)[1].split('_')[0] if sys.argv[1].strip().split('/')[1].split('_')[0] != 'TEST' else 'ISS' #ILS or ISS
    unique_output_name = '-NEWPROB2'

    # Read in SNP data
    obs = loadtxt(sys.argv[1], dtype='string')
    print 'Input file length: ' + str(len(obs))

    # States
    states = ('Unk', 'A', 'ARK', 'BALBc', 'C3HHe', 'C57BL6N', 'DBA2')

    # Start, transition, and emission probabilities
    #  Equal start probability for each state
    #  Transition to same state = 0.46, to other state = 0.09
    #  Emit same state = 0.95, other state = 0.05 (labeled '~State')
    start_p = {s: log(def_start_p) for s in states}
    trans_p = {s_outer: {s_inner: log(def_trans_in_p) if s_inner == s_outer else log(def_trans_out_p) for s_inner in states} for s_outer in states}
    emit_p = {s: {s: log(def_emit_same_p), '~'+s: log(def_emit_other_p)} for s in states}

    # Run algorithms
    print ''
    tot_prob_dist = 10.
    i = 0
    while tot_prob_dist > 0.01 and i < 10:
        tot_prob_dist = 0.
        print '-----RUN %i-----' % i
        path, V = viterbi(obs, states, start_p, trans_p, emit_p, input_group)

        #print_path(path, 'Viterbi')
        #output_path(path)

        new_trans_p = calc_new_trans_p(path, states)
        print '\ntransition probabilities:'
        print trans_p
        print new_trans_p
        tot_prob_dist += prob_dist(trans_p, new_trans_p)
        print ''
        trans_p = new_trans_p

        new_emit_p = calc_new_emit_p(path, obs, states, input_group)
        print 'emission probabilities:'
        print emit_p
        print new_emit_p
        tot_prob_dist += prob_dist(emit_p, new_emit_p)
        print ''
        emit_p = new_emit_p

        print 'total prob distance:'
        print tot_prob_dist

        print ''

        i += 1
