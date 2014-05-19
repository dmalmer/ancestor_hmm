

from collections import defaultdict
from itertools import izip
from math import log

from hmm_util import pairwise, count_hotspots


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

    new_trans_p = {}
    for s_outer in states:
        new_trans_p[s_outer] = {}
        for s_inner in states:
            tot_trans = float(sum(trans_counts[s_outer].values()))
            new_trans_p[s_outer][s_inner] = log(trans_counts[s_outer][s_inner] / tot_trans)

    #new_hs_fi = max(float(tot_hs_trans)/max(tot_hs_not_trans, 1), 1.)
    new_hs_fi = 1

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

