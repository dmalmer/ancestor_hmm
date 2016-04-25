
from collections import defaultdict
from itertools import izip
from math import log, e
from numpy import mean

from util import ancestor_blocks, get_emit_key, pairwise


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
def calc_new_emit_p(ancestors_by_chr, SNPs_by_chr, states, desc_strain, max_emit_same_p):
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
                SNP_key = get_emit_key(s, SNP[3], desc_strain)
                SNP_counts[SNP_key] += 1

            # Count ancestor classifications
            for a in anc.split('_'):
                anc_key = get_emit_key(a, SNP[3], desc_strain)
                ancestor_counts[anc_key] += 1

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
def reclassify_ibd_and_unk(ancestors_by_chr, SNPs_by_chr, desc_strain, use_unknown, unk_cutoff):
    # An ancestor is IBD to the classified ancestor if it has as many or more SNPs in a haplotype block
    #  as the classified ancestor
    # Alternatively, a haplotype block is likely from one of the unsequenced ancestors if a high percentage of SNPs are
    #  labeled the descendant strain
    new_ancestors_by_chr = {}
    for curr_chr in ancestors_by_chr.keys():
        new_ancestors = []
        for pos_start, pos_end, ancestor, SNPs_section in ancestor_blocks(ancestors_by_chr[curr_chr],
                                                                          SNPs_by_chr[curr_chr], return_SNPs=True):
            # Count the number of SNPs from each ancestor
            SNP_counts = defaultdict(int)
            for SNP in SNPs_section:
                if use_unknown and SNP[3] == desc_strain:
                    SNP_counts['Unk'] += 1
                elif desc_strain in SNP[3]:
                    for anc in SNP[3].split('_'):
                        if anc != desc_strain:
                            SNP_counts[anc] += 1

            # First, check if haplotype block comes from an unsequenced ancestor 
            if use_unknown and (int(SNP_counts['Unk']) > 2 and SNP_counts['Unk'] >= unk_cutoff * SNP_counts[ancestor]):
                new_ancestors.extend(['Unk'] * len(SNPs_section))

            # If not, check if classified ancestor is IBD
            elif ancestor != 'Unk':
                # Find all other ancestors with counts >= the classified ancestor count
                indent_ancestors = [k for k,v in SNP_counts.items() if v >= SNP_counts[ancestor] and k != 'Unk']
                # Weird case where entire section is not consistent with descendant strain so counts are 0, very rare
                if len(indent_ancestors) == 0:
                    indent_ancestors = [ancestor]

                # Add all IBD ancestors to back of classified ancestor
                new_ancestors.extend(['_'.join(indent_ancestors)] * len(SNPs_section))

        new_ancestors_by_chr[curr_chr] = new_ancestors

    return new_ancestors_by_chr


# Score results of HMM using structural variants
def score_results(ancestors_by_chr, SNPs_by_chr, desc_ins_by_chr, desc_del_by_chr, anc_ins_by_chr, anc_del_by_chr):
    # For every ancestor block, check for descendant strain structural variants that overlaps the region
    # For every descendant SV, check for individual ancestor SVs that overlap
    # For every set of the ancestor SVs covering the same region
    #  If one of the ancestors is the classified ancestor, score a hit
    #  If none of the ancestors is the classified ancestor, score a miss
    all_scores_by_chr = defaultdict(list)
    for curr_chr, ancestors in ancestors_by_chr.items():
        # first score insertions, then score deletions
        for desc_SVs, anc_SVs, SV_type in ((desc_ins_by_chr[curr_chr], anc_ins_by_chr[curr_chr], 'Ins'),
                                           (desc_del_by_chr[curr_chr], anc_del_by_chr[curr_chr], 'Del')):
            d_sv_ind = 0
            a_sv_ind = 0
            # loop over HMM classification blocks
            for pos_start, pos_end, ancestor in ancestor_blocks(ancestors, SNPs_by_chr[curr_chr]):
                # can't score Unk classifications
                if ancestor == 'Unk':
                    continue

                # if the descendant SV list is empty, we're done
                if len(desc_SVs) == 0:
                    break
                # walk to earliest descendant SV overlapping the current HMM block
                while d_sv_ind < len(desc_SVs) and desc_SVs[d_sv_ind][1] < int(pos_start):
                    d_sv_ind += 1
                #  if we've looped through all the descendant SV's, we're done
                if d_sv_ind == len(desc_SVs):
                    break
                #  if the start of the current descendant SV is past the end of the HMM block, skip (no overlap)
                if desc_SVs[d_sv_ind][0] > int(pos_end):
                    continue

                # otherwise, the descendant SV and HMM block must overlap
                # loop over every descendant SV that overlaps this HMM block
                #  increment a tmp variable for this walk as the current descendant SV may overlap multiple HMM blocks
                d_tmp_ind = d_sv_ind
                while d_tmp_ind < len(desc_SVs) and desc_SVs[d_tmp_ind][0] < int(pos_end):
                    # if the ancestor SV list is empty, score a miss
                    if len(anc_SVs) == 0:
                        # only append miss to all_scores if the current desc_SV hasn't been scored already (hits overwrite misses)
                        if len(all_scores_by_chr[curr_chr]) == 0:
                            all_scores_by_chr[curr_chr].append((pos_start, pos_end, SV_type, 'Miss'))
                        elif pos_start != all_scores_by_chr[curr_chr][-1][0]:
                            all_scores_by_chr[curr_chr].append((pos_start, pos_end, SV_type, 'Miss'))
                        d_tmp_ind += 1
                        continue
                    # walk to earliest ancestor SV overlapping the current descendant SV
                    while a_sv_ind < len(anc_SVs) and anc_SVs[a_sv_ind][1] < desc_SVs[d_tmp_ind][0]:
                        a_sv_ind += 1
                    #  if we've looped through all the ancestor SV's, score a miss
                    if a_sv_ind == len(anc_SVs):
                        # only append miss to all_scores if the current desc_SV hasn't been scored already (hits overwrite misses)
                        if len(all_scores_by_chr[curr_chr]) == 0:
                            all_scores_by_chr[curr_chr].append((pos_start, pos_end, SV_type, 'Miss'))
                        elif pos_start != all_scores_by_chr[curr_chr][-1][0]:
                            all_scores_by_chr[curr_chr].append((pos_start, pos_end, SV_type, 'Miss'))
                        d_tmp_ind += 1
                        continue
                    #  if start of current ancestor SV is past the end of the overlap region, score a miss (no overlap)
                    if anc_SVs[a_sv_ind][0] > desc_SVs[d_tmp_ind][1]:
                        # only append miss to all_scores if the current desc_SV hasn't been scored already (hits overwrite misses)
                        if len(all_scores_by_chr[curr_chr]) == 0:
                            all_scores_by_chr[curr_chr].append((pos_start, pos_end, SV_type, 'Miss'))
                        elif pos_start != all_scores_by_chr[curr_chr][-1][0]:
                            all_scores_by_chr[curr_chr].append((pos_start, pos_end, SV_type, 'Miss'))
                        d_tmp_ind += 1
                        continue
                    # walk over every ancestor SV that overlaps the current descendant SV
                    #  increment a tmp variable for this walk as ancestor SVs can overlap multiple desc SVs
                    a_tmp_ind = a_sv_ind
                    while a_tmp_ind < len(anc_SVs) and anc_SVs[a_tmp_ind][0] < desc_SVs[d_tmp_ind][1]:
                        if anc_SVs[a_tmp_ind][1] > desc_SVs[d_tmp_ind][0]:
                            # check if any ancestor SVs match the HMM output
                            for anc in ancestor.split('_'):
                                if anc in anc_SVs[a_tmp_ind][2]:
                                    # if the current desc_SV has already been scored, overwrite with a hit
                                    if len(all_scores_by_chr[curr_chr]) == 0:
                                        all_scores_by_chr[curr_chr].append((pos_start, pos_end, SV_type, 'Hit'))
                                    elif pos_start == all_scores_by_chr[curr_chr][-1][0]:
                                        all_scores_by_chr[curr_chr][-1] = (pos_start, pos_end, SV_type, 'Hit')
                                    else:
                                        all_scores_by_chr[curr_chr].append((pos_start, pos_end, SV_type, 'Hit'))
                                    break
                            # clever way to break out of both for loops (taken from stackoverflow)
                            else:
                                a_tmp_ind += 1
                                continue
                            break
                        a_tmp_ind += 1
                    # if we reach this else statement, the while loop never called break, so must be a miss
                    else:
                        # only append miss to all_scores if the current desc_SV hasn't been scored already (hits overwrite misses)
                        if len(all_scores_by_chr[curr_chr]) == 0:
                            all_scores_by_chr[curr_chr].append((pos_start, pos_end, SV_type, 'Miss'))
                        elif pos_start != all_scores_by_chr[curr_chr][-1][0]:
                            all_scores_by_chr[curr_chr].append((pos_start, pos_end, SV_type, 'Miss'))
                    d_tmp_ind += 1

    # Count total hits and misses
    hits_by_chr = defaultdict(int)
    misses_by_chr = defaultdict(int)
    for curr_chr, all_scores in all_scores_by_chr.items():
        for score in all_scores:
            if score[3] == 'Hit':
                hits_by_chr[curr_chr] += 1
            elif score[3] == 'Miss':
                misses_by_chr[curr_chr] += 1

    return hits_by_chr, misses_by_chr, all_scores_by_chr



