
import sys
from numpy import loadtxt, zeros
from math import *
from collections import defaultdict

# Output
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

# Baum-Welch algorithm for expectation maximization
# **** THIS DOES NOT WORK ****
#  I tried to write this mostly based on the wikipedia B-W page, but after a few
#   generations it quickly ends up with transition probabilities to other states 
#   higher than the transition probabilities to the same state, so something is wrong.
#  I really recommend scrapping this entirely a re-writing with a cleaner 
#   foward-backward approach
def em():
    #initialize sums and maximums to zero
    trans_sums = {s_outer: {s_inner: 0. for s_inner in states} for s_outer in states}
    trans_max = 0.

    emit_sums = {s: {s: 0., '~'+s: 0.} for s in states}
    emit_max = 0.

    #find new transition probabilities:
    #  at each observation, find the prob of transitioning from the prev obs to the current obs for each state transition
    #  sum these state transition prob, and
    #  keep track of the highest prob state path
    max_seq = defaultdict(float)
    prev_obs = obs[0][3]
    for i in range(1, 100000):#len(obs)):
        if i % 10000 == 0:
            print i
        cur_obs = obs[i][3]
        max_prob = -9999.
        for s in states:
            cur_e = s
            if s == 'Unk':
                if cur_obs != input_group:
                    cur_e = '~' + s
            elif s not in cur_obs.split('_'):
                cur_e = '~' + s
            for p in states:
                prev_e = p
                if p == 'Unk':
                    if prev_obs != input_group:
                        prev_e = '~' + p
                elif p not in prev_obs.split('_'):
                    prev_e = '~' + p
                trans_prob = start_p[p] + trans_p[p][s] + emit_p[s][cur_e] + emit_p[p][prev_e]
                #print 's = ' + s + ', p = ' + p + ': ' + str(trans_prob)
                trans_sums[p][s] += trans_prob
                if trans_prob > max_prob:
                    max_prob = trans_prob
        trans_max += max_prob
        prev_obs = cur_obs

    #normalize transitions of each state to 1 in base 10, then move back to log space
    print 'before:'
    for t in trans_p.keys():
        print t + ': ' + str(trans_p[t])
    for k in trans_sums.keys():
        sum = 0.
        for s in states:
            sum += exp(-1*trans_sums[k][s]/trans_max)
        print k + ' sum: ' + str(sum)
        new = 0.
        for s in states:
            trans_p[k][s] = (-1*trans_sums[k][s]/trans_max)/sum
            print 'k=' + k + ',s=' + s +': ' + str(trans_p[k][s])
            new += exp(-1*trans_sums[k][s]/trans_max)/sum
        print new
    print 'after:'
    for t in trans_p.keys():
        print t + ': ' + str(trans_p[t])

    print ''
    for m in max_seq:
        print m + ': ' + str(max_seq[m])

# Viterbi algorithm
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
        if i%10000 == 0:
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

    return path[state]

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
    start_p = {s: log(1/7.) for s in states}
    trans_p = {s_outer: {s_inner: log(0.46) if s_inner == s_outer else log(0.09) for s_inner in states} for s_outer in states}
    emit_p = {s: {s: log(0.95), '~'+s: log(0.05)} for s in states}

    # Run algorithms
    path = viterbi(obs, states, start_p, trans_p, emit_p, input_group)

    print_path(path, 'Viterbi')
    output_path(path)
