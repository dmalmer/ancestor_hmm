

from collections import defaultdict
from itertools import izip
from math import log, e

from hmm_util import ancestor_blocks, log_add_list, pairwise


def calc_new_trans_p(ancestors_by_chr, states):
    # set minimum transition counts to 1 so we don't have any 0 probabilities
    trans_counts = {s_outer: {s_inner: 1 for s_inner in states} for s_outer in states}
    for ancestors in ancestors_by_chr.values():
        for anc_prev, anc_curr in pairwise(ancestors):
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
            if ancestors[i] in SNPs_by_chr[curr_chr][i][3].split('_'):
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


def calc_confidence_intervals(ancestors, SNPs, prob_nodes):
    confidence_intervals = []
    for anc, SNP, p_node in zip(ancestors, SNPs, prob_nodes):
        p_sum = log_add_list(list(p_node))
        # (chromosome, start_pos, confidence)
        confidence_intervals.append((SNP[0], SNP[1], e ** (p_node[anc] - p_sum)))
    return confidence_intervals


def calc_identical_by_ancestor(ancestors, SNPs, trans_p, emit_p, states, IBA_cutoff):
    prev_anc = ''
    identical_by_anc = []
    for chromosome, pos_start, pos_end, ancestor, SNPs_section in ancestor_blocks(ancestors, SNPs, return_SNPs=True):
        for s in states:
            pass
        print SNPs_section

        prev_anc = ancestor

def calc_identical_by_ancestor_old(ancestors, SNPs, prob_nodes, states, IBA_cutoff):

    identical_by_anc = []
    for chromosome, pos_start, pos_end, ancestor, prob_section in ancestor_blocks(ancestors, SNPs, prob_nodes=prob_nodes):

        diff_ratios = defaultdict(list)
        #print 'p_nodes:'
        for p_node in prob_section:
            anc_prob = p_node[ancestor]
            for s in states:
                if s != ancestor:
                    diff_ratios[s].append(e ** (p_node[s] - anc_prob))
        IBA = [s for s in states if s != ancestor and sum(diff_ratios[s])/len(diff_ratios[s]) > IBA_cutoff]

        # if pos_start == '90310154':
        #     print 'at pos = 90310154:'
        #     print 'anc_prob = ' + str(prob_section[ancestor])
        #     for s in states:
        #         print 'state = ' + str(s)
        #         print prob_section[s]
        #         print diff_ratios[s]
        #         print sum(diff_ratios[s])
        #         print len(diff_ratios[s])
        #
        # if pos_start == '63749154':
        #     print 'C3HHe, pos ' + str(pos_start)
        #     print prob_section[ancestor]
        #     print prob_section['A']
        #     print diff_ratios['A']
        #     print [sum(diff_ratios[s])/len(diff_ratios[s]) for s in states if s != ancestor]
        #
        # if any(IBA):
        #     print 'IBAs!!'
        #     print IBA
        #     identical_by_anc.append((chromosome, pos_start, pos_end, '%s_%s' % (ancestor, '_'.join(IBA))))

        #
        # sum_probs = 0.
        # cur_iba = []
        # anc_prob_sum = log_add_list(prob_section[ancestor])
        # print 'section:'
        # print prob_section[ancestor]
        # print 'anc_prob_sum:'
        # print anc_prob_sum
        # print 'state sums:'
        # for s in states:
        #     if s != ancestor:
        #         print s
        #         print 'log space:'
        #         print log_add_list(prob_section[s])
        #         print anc_prob_sum - log_add_list(prob_section[s])
        #
        #         print 'base 10:'
        #         print e ** log_add_list(prob_section[s])
        #         print (e ** anc_prob_sum) / (e ** log_add_list(prob_section[s]))
        #
        #         if sum(prob_section[ancestor])/sum(prob_section[s]) > iba_cutoff:
        #             cur_iba.append(s)
        #     sum_probs += sum(prob_section[s])
        #
        # if any(cur_iba):
        #     identical_by_anc.append((chromosome, pos_start, pos_end, '_'.join(cur_iba)))
        #
        # # conf_int = prob of flipping a coin [1/7] - prob of choice [prob_anc / sum(prob_every_state)]
        # confidence_interval.append((chromosome, pos_start, pos_end, (1./7) - sum(prob_section[ancestor])/sum_probs))

        #if sum(prob_section[ancestor])/sum_probs < .1:
        #    print prob_section
        #    print sum(prob_section[ancestor])/sum_probs

    return identical_by_anc


def prob_dist(old_probs, new_probs):
    tot_dist = 0.
    tot_probs = 0
    for key in old_probs.keys():
        for old_p, new_p in izip(old_probs[key].values(), new_probs[key].values()):
            tot_dist += abs(old_p - new_p)
            tot_probs += 1

    # return average prob dist
    return tot_dist / tot_probs

