
from collections import defaultdict
from itertools import izip
from math import log, e
from numpy import mean

from hmm_util import ancestor_blocks, log_add_list, pairwise


# Calculate a confidence interval for each classified haplotype block
def calc_confidence_intervals(ancestors, SNPs, prob_nodes):
    confidence_intervals = []
    for anc, SNP, p_node in zip(ancestors, SNPs, prob_nodes):
        p_sum = log_add_list(list(p_node))
        # (chromosome, start_pos, confidence)
        confidence_intervals.append((SNP[0], SNP[1], e ** (p_node[anc] - p_sum)))
    return confidence_intervals


# Calculate new transition probabilities
def calc_new_trans_p(ancestors_by_chr, states):
    # Count number of transitions from one state to another
    #  (set minimum transition counts to 1 so we don't have any 0 probabilities)
    trans_counts = {s_outer: {s_inner: 1 for s_inner in states} for s_outer in states}
    for ancestors in ancestors_by_chr.values():
        for anc_prev_grp, anc_curr_grp in pairwise(ancestors):
            for anc_prev in anc_prev_grp.split('_'):
                for anc_curr in anc_curr_grp.split('_'):
                    trans_counts[anc_prev][anc_curr] += 1

    # New transition rate of s1->s2 is count(s1->s2)/count(s1->*)
    new_trans_p = {}
    for s_outer in states:
        new_trans_p[s_outer] = {}
        for s_inner in states:
            tot_trans = float(sum(trans_counts[s_outer].values()))
            new_trans_p[s_outer][s_inner] = log(trans_counts[s_outer][s_inner] / tot_trans)

    return new_trans_p


# Calculate new emission probabilities
def calc_new_emit_p(ancestors_by_chr, SNPs_by_chr, states, input_strain, max_emit_same_p):
    # SNP counts
    #  SNP_counts[<state>] = total number of <state> appearances
    #  SNP_counts[<~state>] = total number of <~state> appearances
    SNP_counts = defaultdict(int)

    # Ancestor counts
    #  ancestor_counts[<state>] = total number of calculated <state> when <state> is observed
    #  ancestor_counts[<~state>] = total number of calculated <state> when <state> is NOT observed
    ancestor_counts = defaultdict(int)

    for curr_chr, ancestors in ancestors_by_chr.items():
        for SNP, anc in zip(SNPs_by_chr[curr_chr], ancestors):
            # Count SNP observations
            for s in states:
                SNP_key = s
                if s == 'Unk':
                    if SNP[3] != input_strain:
                        SNP_key = '~' + s
                elif s not in SNP[3].split('_'):
                    SNP_key = '~' + s
                SNP_counts[SNP_key] += 1

            # Count ancestor classifications
            for a in anc.split('_'):
                if a in SNP[3].split('_'):
                    ancestor_counts[a] += 1
                else:
                    ancestor_counts['~'+a] += 1

    # New emission rates of <state> and <~state> are ancestor_counts[<state>]/SNP_counts[<state>] and
    #  ancestor_counts[<~state>]/SNP_counts[<~state>] respectively, normalized to 1.0
    new_emit_p = {}
    for s in states:
        new_emit_p[s] = {}
        if SNP_counts[s] == 0 or ancestor_counts[s] == 0:
            # If we don't observe or calculate a state, use the max emission probabilities
            new_emit_p[s][s] = log(max_emit_same_p)
            new_emit_p[s]['~'+s] = log(1 - max_emit_same_p)
        else:
            # Normalize <state> and <~state> probabilities to 1.0
            state_p = float(ancestor_counts[s]) / SNP_counts[s]
            notstate_p = float(ancestor_counts['~'+s]) / SNP_counts['~'+s]
            normalizer = 1 / (state_p + notstate_p)
            # Don't allow probabilities of 1.0 and 0.0
            new_emit_p[s][s] = log(min(normalizer * state_p, max_emit_same_p))
            new_emit_p[s]['~'+s] = log(max(normalizer * notstate_p, 1 - max_emit_same_p))

    # Take mean of each state to fix emission probabilities across states
    mean_same = mean([e ** new_emit_p[s][s] for s in states])
    mean_other = mean([e ** new_emit_p[s]['~'+s] for s in states])
    new_emit_p = {s: {s: log(mean_same), '~'+s: log(mean_other)} for s in states}

    return new_emit_p


# Calculate expected number of recombinations between two SNPs
def calc_recomb_rate(SNP_start, SNP_end, recomb_main_i, recomb_map, effective_pop, num_generations):
    # Find recomb_map starting position
    if recomb_main_i is None:
        recomb_main_i = 0
        while recomb_main_i < len(recomb_map) and int(recomb_map[recomb_main_i][0]) < SNP_start:
            recomb_main_i += 1
        recomb_start_i = max(recomb_main_i-1, 0)
    else:
        recomb_start_i = recomb_main_i - 1

    # Quick check to make sure recomb_start_i >= 0 (should be, as recomb_index should always be >=1 in else statement)
    #  --remove this later
    if recomb_start_i < 0:
        raise Exception('recomb_start_i should never be less than 0')

    # Find recomb_map ending position
    while (recomb_main_i < len(recomb_map) and int(recomb_map[recomb_main_i][0]) < SNP_end) or recomb_main_i == 0:
        recomb_main_i += 1
    recomb_end_i = recomb_main_i

    # Calc recomb rates
    #  First, special case for SNPs between adjacent genetic markers
    if recomb_end_i - recomb_start_i == 1:
        expected_recombs = ((((SNP_end - SNP_start) / 1000.) * recomb_map[recomb_start_i][1]) / (4 * effective_pop)) * \
                            num_generations

    #  Otherwise, calc with all recomb rates between SNPs
    else:
        # Proportional rate of first SNP
        expected_recombs = ((((recomb_map[recomb_start_i+1][0] - SNP_start) / 1000.) * recomb_map[recomb_start_i][1]) /
                            (4 * effective_pop)) * num_generations

        # Rates in the middle
        for i in range(recomb_start_i+1, recomb_end_i-1):
            expected_recombs += ((((recomb_map[i+1][0] - recomb_map[i][0]) / 1000.) * recomb_map[i][1]) /
                                 (4 * effective_pop)) * num_generations

        # Proportional rate of second SNP
        expected_recombs += ((((SNP_end - recomb_map[recomb_end_i-1][0]) / 1000.) * recomb_map[recomb_end_i-1][1]) /
                             (4 * effective_pop)) * num_generations

    # Expected_recombs must be > 0 in order to convert to log space
    return max(expected_recombs, .00000001), recomb_main_i


# Calculate the distance between old and new transition or emission rates
def prob_dist(old_probs, new_probs):
    # Probability distance = abs(log-space prob of old t/e rate - log-space prob of new t/e rate)
    tot_dist = 0.
    tot_probs = 0
    for key in old_probs.keys():
        for old_p, new_p in izip(old_probs[key].values(), new_probs[key].values()):
            tot_dist += abs(old_p - new_p)
            tot_probs += 1

    # Return average prob dist
    return tot_dist / tot_probs


# Reclassify all haplotype blocks that are identical by descent or likely from one of the two unknown strains
def reclassify_ibd_and_unk(ancestors_by_chr, SNPs_by_chr, input_strain, unk_cutoff):
    # An ancestor is IBD to the classified ancestor if it has as many or more SNPs in a haplotype block
    #  as the classified ancestor
    # Alternatively, a haplotype block is likely from one of the unsequenced ancestors if a high percentage of SNPs are
    #  labeled ILS or ISS
    new_ancestors_by_chr = {}
    for curr_chr in ancestors_by_chr.keys():
        new_ancestors = []
        for chromosome, pos_start, pos_end, ancestor, SNPs_section in ancestor_blocks(ancestors_by_chr[curr_chr],
                                                                                      SNPs_by_chr[curr_chr],
                                                                                      return_SNPs=True):
            # Count the number of SNPs from each ancestor
            SNP_counts = defaultdict(int)
            for SNP in SNPs_section:
                if SNP[3] == input_strain:
                    SNP_counts['Unk'] += 1
                else:
                    for anc in SNP[3].split('_')[1:]:
                        SNP_counts[anc] += 1

            # First, check if haplotype block comes from an unsequenced ancestor
            if int(SNP_counts['Unk']) > 2 and SNP_counts['Unk'] >= unk_cutoff * SNP_counts[ancestor]:
                new_ancestors.extend(['Unk'] * len(SNPs_section))

            # If not, check if classified ancestor is IBD
            else:
                # Find all other ancestors with counts >= the classified ancestor count
                indent_ancestors = [k for k,v in SNP_counts.items() if v >= SNP_counts[ancestor]]

                # Add all IBD ancestors to back of classified ancestor
                new_ancestors.extend(['_'.join(indent_ancestors)] * len(SNPs_section))

        new_ancestors_by_chr[curr_chr] = new_ancestors

    return new_ancestors_by_chr
