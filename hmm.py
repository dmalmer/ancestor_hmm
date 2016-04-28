
import math
import numpy

from math import log
from time import time

from prob import calc_new_emit_p, calc_new_trans_p, calc_recomb_rate, prob_dist, reclassify_ibd_and_unk, score_results
from util import get_emit_key, prob_tuples, write_ancestors, write_scores, write_statistics


# Viterbi algorithm
def viterbi(SNPs, trans_p, emit_p, adjust_recomb, sp):
    # Initialize
    prob_nodes = numpy.zeros(len(SNPs), dtype={'names': sp['states'], 'formats': ['f8']*len(sp['states'])})
    ancestors_by_state = {}
    recomb_index = None  # Initially set to None so calc_recomb_rate uses a special case when called for the first time

    # Start probabilities
    for s in sp['states']:
        # For each state, the emission probability is either emit_p[state] or emit_p[~state]
        emit_key = get_emit_key(s, SNPs[0][3], sp['desc_strain'])
        # Probability of a given state at SNPs[0] is emit prob of stat
        prob_nodes[0][s] = emit_p[s][emit_key]
        ancestors_by_state[s] = [s]

    # Rest of probabilities
    for i in range(1, len(SNPs)):
        #if sp['verbose'] and i % 1000 == 0:
        #    print 'i = ' + str(i)
        new_ancestors_by_state = {}

        # Calculate recombination rate
        expected_recombs = 1.
        if sp['use_recomb_rates'] and len(sp['recomb_rate_dict'][SNPs[i][0]]) > 0:
            expected_recombs, recomb_index = calc_recomb_rate(int(SNPs[i-1][1]), int(SNPs[i][1]), recomb_index,
                                                              sp['recomb_rate_dict'][SNPs[i][0]], sp['effective_pop'],
                                                              sp['num_generations'])
            expected_recombs *= adjust_recomb

        # At every SNP, find probabilities for each state
        for curr_state in sp['states']:
            # For each state, the emission probability is either emit_p[state] or emit_p[~state]
            emit_key = get_emit_key(curr_state, SNPs[i][3], sp['desc_strain'])

            # The probability of a given state for a given SNP is the maximum out of ((prob prev_state) *
            #  (prob trans prev_state -> curr_state) * (prob emit curr_state)) for all previous states
            state_probabilities = []
            for prev_state in sp['states']:
                state_probabilities.append((prob_nodes[i-1][prev_state] + trans_p[prev_state][curr_state] +
                                            emit_p[curr_state][emit_key], prev_state))

            (recomb_prob, orig_prob, prev_state) = max([(prob + math.log(expected_recombs), prob, prev_s)
                                                        if curr_state != prev_s else (prob, prob, prev_s)
                                                        for prob, prev_s in state_probabilities])

            # Keep track of probabilities in prob_nodes and ancestors in new_ancestors_by_state for each curr_state
            prob_nodes[i][curr_state] = orig_prob
            new_ancestors_by_state[curr_state] = ancestors_by_state[prev_state] + [curr_state]

        # Update ancestors with additional iteration
        ancestors_by_state = new_ancestors_by_state

    # Find maximum-likelihood ancestors for the current chromosome
    (prob, best_state) = max((prob_nodes[len(SNPs)-1][s], s) for s in sp['states'])
    return ancestors_by_state[best_state]


# Expectation-Maximization loop
def expectation_maximization(SNPs_by_chr, trans_in_p, emit_same_p, adjust_recomb, use_unknown, unk_cutoff, append_str, sp,
                             job_server, vit_func):
    # Start, transition, and emission probabilities
    trans_p = {s_outer: {s_inner: log(trans_in_p) if s_inner == s_outer else log((1 - trans_in_p) / 6)
                         for s_inner in sp['states']} for s_outer in sp['states']}
    emit_p = {s: {s: log(emit_same_p), '~'+s: log(1 - emit_same_p)} for s in sp['states']}

    # EM loop
    tot_prob_dist = 10.
    run_count = 0
    while tot_prob_dist > sp['prob_dist_cutoff'] and run_count < int(sp['max_iter']):
        print '---- EM loop %i ----' % run_count
        tot_prob_dist = 0.

        # Run viterbi to find maximum likelihood path
        ancestors_by_chr = {}
        if sp['parallel']:
            jobs = []
            for curr_chr, SNPs in SNPs_by_chr.items():
                job = vit_func.submit(SNPs, trans_p, emit_p, adjust_recomb, sp)
                jobs.append((curr_chr, job))

            for curr_chr, job in jobs:
                ancestors_by_chr[curr_chr] = job()
            job_server.wait()

            if sp['verbose']:
                job_server.print_stats()
        else:
            for curr_chr, SNPs in SNPs_by_chr.items():
                ancestors_by_chr[curr_chr] = viterbi(SNPs, trans_p, emit_p, adjust_recomb, sp)

        # Recalculate transition and emission probabilities
        new_trans_p = calc_new_trans_p(ancestors_by_chr, sp['states'])
        tot_prob_dist += prob_dist(trans_p, new_trans_p)

        new_emit_p = calc_new_emit_p(ancestors_by_chr, SNPs_by_chr, sp['states'], sp['desc_strain'], emit_same_p)
        tot_prob_dist += prob_dist(emit_p, new_emit_p)

        if sp['verbose']:
            print 'Transition probabilities'
            print '  Before:'
            tuples_dict = prob_tuples(trans_p)
            for s in sp['states']:
                print '    %s ->' % s,
                for tup in tuples_dict[s]:
                    print '%s: %.2f,' % tup,
                print ''
            print '  After:'
            tuples_dict = prob_tuples(new_trans_p)
            for s in sp['states']:
                print '    %s ->' % s,
                for tup in tuples_dict[s]:
                    print '%s: %.2f,' % tup,
                print ''

            print '\nEmission probabilities'
            print '  Before:'
            tuples_dict = prob_tuples(emit_p)
            for s in sp['states']:
                print '    %s ->' % s,
                for tup in tuples_dict[s]:
                    print '%s: %.2f,' % tup,
                print ''
            print '  After:'
            tuples_dict = prob_tuples(new_emit_p)
            for s in sp['states']:
                print '    %s ->' % s,
                for tup in tuples_dict[s]:
                    print '%s: %.2f,' % tup,
                print ''

            print 'Total probability distance: %.3f' % tot_prob_dist

        # If set, write out results at each iteration
        if sp['write_iter']:
            # Score results
            final_score = -1
            if len(sp['desc_ins_by_chr']) > 0:
                hits_by_chr, misses_by_chr, all_scores_by_chr = score_results(ancestors_by_chr, SNPs_by_chr,
                                                                              sp['desc_ins_by_chr'], sp['desc_del_by_chr'],
                                                                              sp['anc_ins_by_chr'], sp['anc_del_by_chr'])
                try:
                    final_score = sum(hits_by_chr.values())/float(sum(misses_by_chr.values()))
                except ZeroDivisionError:
                    final_score = -1

                # Output scores
                write_scores(sp['output_dir'], sp['filename_in'], append_str, all_scores_by_chr)

            # Output results (results during EM iterations will not have reclassified IBD and Unk)
            write_ancestors(sp['output_dir'], sp['filename_in'], append_str, ancestors_by_chr, SNPs_by_chr, sp['state_rgbs'])
            write_statistics(sp['output_dir'], sp['filename_in'], append_str, sp['states'], ancestors_by_chr, SNPs_by_chr,
                             (sp['use_recomb_rates'], trans_in_p, emit_same_p, trans_p, emit_p), final_score,
                             run_count, time() - sp['time_start'], tot_prob_dist)

        trans_p = new_trans_p
        emit_p = new_emit_p

        run_count += 1

    # Reclassify segments where segment likely came from an unsequenced ancestor or where SNP counts for multiple
    #  ancestors are identical
    ancestors_by_chr = reclassify_ibd_and_unk(ancestors_by_chr, SNPs_by_chr, sp['desc_strain'], use_unknown, unk_cutoff)

    print '\nTotal time (min): ' + str((time() - sp['time_start'])/60)
    print 'Total runs: ' + str(run_count)

    # Score results
    final_score = -1
    if len(sp['desc_ins_by_chr']) > 0:
        hits_by_chr, misses_by_chr, all_scores_by_chr = score_results(ancestors_by_chr, SNPs_by_chr, sp['desc_ins_by_chr'],
                                                                      sp['desc_del_by_chr'], sp['anc_ins_by_chr'],
                                                                      sp['anc_del_by_chr'])
        try:
            final_score = sum(hits_by_chr.values())/float(sum(misses_by_chr.values()))
        except ZeroDivisionError:
            final_score = -1

        print '\nChromosome scores:'
        for curr_chr in hits_by_chr.keys():
            try:
                chr_score = hits_by_chr[curr_chr]/float(misses_by_chr[curr_chr])
            except ZeroDivisionError:
                chr_score = -1
            print ' %s - Hits: %i, Misses: %i, Ratio: %.3f' % (curr_chr, hits_by_chr[curr_chr], misses_by_chr[curr_chr],
                                                               chr_score)
        print 'Final score: %.3f' % (final_score)

        # Output scores
        write_scores(sp['output_dir'], sp['filename_in'], append_str, all_scores_by_chr)

    # Output final results after reclassifying IBD and Unk
    write_ancestors(sp['output_dir'], sp['filename_in'], append_str, ancestors_by_chr, SNPs_by_chr, sp['state_rgbs'])
    write_statistics(sp['output_dir'], sp['filename_in'], append_str, sp['states'], ancestors_by_chr, SNPs_by_chr, 
                     (sp['use_recomb_rates'], trans_in_p, emit_same_p, trans_p, emit_p), final_score, run_count, 
                     time() - sp['time_start'], tot_prob_dist)
