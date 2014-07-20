

from collections import defaultdict
from itertools import izip
from math import log

from hmm_util import ancestor_blocks, pairwise


def calc_new_trans_p(ancestors, states):
    # set minimum transition counts to 1 so we don't have any 0 probabilities
    trans_counts = {s_outer: {s_inner: 1 for s_inner in states} for s_outer in states}
    for anc_prev, anc_curr in pairwise(ancestors):
            trans_counts[anc_prev][anc_curr] += 1

    new_trans_p = {}
    for s_outer in states:
        new_trans_p[s_outer] = {}
        for s_inner in states:
            tot_trans = float(sum(trans_counts[s_outer].values()))
            new_trans_p[s_outer][s_inner] = log(trans_counts[s_outer][s_inner] / tot_trans)

    return new_trans_p


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


def posterior_probabilities(ancestors, SNPs, prob_nodes, states, iba_cutoff):
    identical_by_anc = []
    confidence_interval = []

    for chromosome, pos_start, pos_end, ancestor, prob_section in ancestor_blocks(ancestors, SNPs, prob_nodes=prob_nodes):
        print ancestor
        print prob_section

        sum_probs = 0.
        cur_iba = []
        for s in states:
            if s != ancestor:
                if sum(prob_section[ancestor])/sum(prob_section[s]) > iba_cutoff:
                    cur_iba.append(s)
            sum_probs += sum(prob_section[s])

        if any(cur_iba):
            identical_by_anc.append(chromosome, pos_start, pos_end, '_'.join(cur_iba))
        confidence_interval.append((chromosome, pos_start, pos_end, 1./sum(prob_section[ancestor])/sum_probs))

    return (identical_by_anc, confidence_interval)


def prob_dist(old_probs, new_probs):
    tot_dist = 0.
    tot_probs = 0
    for key in old_probs.keys():
        for old_p, new_p in izip(old_probs[key].values(), new_probs[key].values()):
            tot_dist += abs(old_p - new_p)
            tot_probs += 1

    # return average prob dist
    return tot_dist / tot_probs

