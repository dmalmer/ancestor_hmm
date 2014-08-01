

from collections import defaultdict
from itertools import izip
from math import log, e
from numpy import mean

from hmm_util import ancestor_blocks, log_add_list, pairwise


def calc_new_trans_p(ancestors_by_chr, states):
    # set minimum transition counts to 1 so we don't have any 0 probabilities
    trans_counts = {s_outer: {s_inner: 1 for s_inner in states} for s_outer in states}
    for ancestors in ancestors_by_chr.values():
        for anc_prev_grp, anc_curr_grp in pairwise(ancestors):
            for anc_prev in anc_prev_grp.split('_'):
                for anc_curr in anc_curr_grp.split('_'):
                    trans_counts[anc_prev][anc_curr] += 1

    new_trans_p = {}
    for s_outer in states:
        new_trans_p[s_outer] = {}
        for s_inner in states:
            tot_trans = float(sum(trans_counts[s_outer].values()))
            new_trans_p[s_outer][s_inner] = log(trans_counts[s_outer][s_inner] / tot_trans)

    return new_trans_p


def calc_new_emit_p(ancestors_by_chr, SNPs_by_chr, states, input_group, def_emit_same_p, def_emit_other_p):
    # SNP_counts[<state>] = total number of <state> appearances
    # SNP_counts[<~state>] = total number of <~state> appearances
    SNP_counts = defaultdict(int)
    # ancestor_counts[<state>] = total number of calculated <state> when <state> is observed
    # ancestor_counts[<~state>] = total number of calculated <state> when <state> is NOT observed
    ancestor_counts = defaultdict(int)

    for curr_chr, ancestors in ancestors_by_chr.items():
        for i in range(len(SNPs_by_chr[curr_chr])):
            # SNP_counts
            for s in states:
                SNP_key = s
                if s == 'Unk':
                    if SNPs_by_chr[curr_chr][i][3] != input_group:
                        SNP_key = '~' + s
                elif s not in SNPs_by_chr[curr_chr][i][3].split('_'):
                    SNP_key = '~' + s
                SNP_counts[SNP_key] += 1

            # ancestor_counts
            for anc in ancestors[i].split('_'):
                if anc in SNPs_by_chr[curr_chr][i][3].split('_'):
                    ancestor_counts[anc] += 1
                else:
                    ancestor_counts['~'+anc] += 1

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

    # Take mean of each state to fix emission probabilities across states
    mean_same = mean([e ** new_emit_p[s][s] for s in states])
    mean_other = mean([e ** new_emit_p[s]['~'+s] for s in states])
    new_emit_p = {s: {s: log(mean_same), '~'+s: log(mean_other)} for s in states}

    return new_emit_p


def calc_confidence_intervals(ancestors, SNPs, prob_nodes):
    confidence_intervals = []
    for anc, SNP, p_node in zip(ancestors, SNPs, prob_nodes):
        p_sum = log_add_list(list(p_node))
        # (chromosome, start_pos, confidence)
        confidence_intervals.append((SNP[0], SNP[1], e ** (p_node[anc] - p_sum)))
    return confidence_intervals


def label_identical_ancestors(ancestors_by_chr, SNPs_by_chr, input_group):

    new_ancestors_by_chr = {}
    for curr_chr in ancestors_by_chr.keys():
        new_ancestors = []
        for chromosome, pos_start, pos_end, ancestor, SNPs_section in ancestor_blocks(ancestors_by_chr[curr_chr], \
                                                                                      SNPs_by_chr[curr_chr], \
                                                                                      return_SNPs=True):
            SNP_counts = defaultdict(int)
            for SNP in SNPs_section:
                if SNP[3] == input_group:
                    SNP_counts['Unk'] += 1
                else:
                    for anc in SNP[3].split('_')[1:]:
                        SNP_counts[anc] += 1

            indent_ancestors = [k for k,v in SNP_counts.items() if v >= SNP_counts[ancestor]]

            new_ancestors.extend(['_'.join(indent_ancestors)] * len(SNPs_section))

        new_ancestors_by_chr[curr_chr] = new_ancestors

    return new_ancestors_by_chr

def prob_dist(old_probs, new_probs):
    tot_dist = 0.
    tot_probs = 0
    for key in old_probs.keys():
        for old_p, new_p in izip(old_probs[key].values(), new_probs[key].values()):
            tot_dist += abs(old_p - new_p)
            tot_probs += 1

    # return average prob dist
    return tot_dist / tot_probs

