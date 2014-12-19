
from collections import defaultdict
from datetime import datetime
from os.path import isfile

from hmm_util import ancestor_blocks, natural_keys


# Write ancestor classifications to a .bed file
def write_ancestors(working_dir, filename_in, unique_output_name, ancestors_by_chr, SNPs_by_chr, state_RGBs):
    with open(working_dir + '/results/' + filename_in.rsplit('.', 1)[0] + '_hmm-out' + unique_output_name +
                    '.' + filename_in.rsplit('.', 1)[1], 'w') as f_out:
        for curr_chr in sorted(ancestors_by_chr.keys(), key=natural_keys):
            for pos_start, pos_end, ancestor in ancestor_blocks(ancestors_by_chr[curr_chr], SNPs_by_chr[curr_chr]):
                color_key = ancestor
                if '_' in color_key:
                    color_key = 'IBA'
                f_out.write('{0}\t{1}\t{2}\t{3}\t0\t+\t{1}\t{2}\t{4}\n'.format(
                    curr_chr,  # {0} - chromosome
                    pos_start,  # {1} - start pos
                    pos_end,  # {2} - end pos
                    ancestor,  # {3} - ancestor label
                    state_RGBs[color_key],  # {4} - label color
                ))


# Write hits and misses to a .bed file
def write_scores(working_dir, filename_in, unique_output_name, all_scores_by_chr):
    with open(working_dir + '/results/' + filename_in.rsplit('.', 1)[0] + '_scores' + unique_output_name +
                    '.' + filename_in.rsplit('.', 1)[1], 'w') as f_out:
        for curr_chr in sorted(all_scores_by_chr.keys(), key=natural_keys):
            for score in all_scores_by_chr[curr_chr]:
                color = '51,255,51' if score[4] == 'Hit' else '255,51,51'
                f_out.write('{0}\t{1}\t{2}\t{3}\t0\t+\t{1}\t{2}\t{4}\n'.format(
                    curr_chr,  # {0} - chromosome
                    score[0],  # {1} - start pos
                    score[1],  # {2} - end pos
                    score[4],  # {3} - 'Hit' or 'Miss'
                    color  # {4} - green for hit, red for miss
                ))


# Write statistics of each run out to a file
def write_statistics(working_dir, filename_in, unique_output_name, ancestors_by_chr, SNPs_by_chr, starting_params,
                     final_score, run_count, tot_run_time, final_prob_dist):
    len_before = 0
    for SNPs in SNPs_by_chr.values():
        len_before += len(SNPs)
    len_after = 0
    anc_counts = defaultdict(int)

    for curr_chr in ancestors_by_chr.keys():
        for pos_start, pos_end, ancestor in ancestor_blocks(ancestors_by_chr[curr_chr], SNPs_by_chr[curr_chr]):
            anc_counts[ancestor] += 1
            len_after += 1

    line = '{0}\t{1:.3f}\t{2}\t{3}\t{4:.2%}\t{5:.3f}\t{6:.3f}\t{7:.3f}\t{8:.3f}\t{9:.3f}\t{10:.3f}\t{11:.3f}\t'.format(
                datetime.now().strftime('%m/%d/%y-%H:%M'),  # {0} - date and time of run
                final_score, # {1} - final score
                len_before,  # {2} - input file length
                len_after,  # {3} - output file length
                float(len_before-len_after)/len_before,  # {4} - percentage difference
                anc_counts['A']/float(len_after),  # {5}
                anc_counts['AKR']/float(len_after),  # {6}
                anc_counts['BALBc']/float(len_after),  # {7}
                anc_counts['C3HHe']/float(len_after),  # {8}
                anc_counts['C57BL6N']/float(len_after),  # {9}
                anc_counts['DBA2']/float(len_after),  # {10}
                anc_counts['Unk']/float(len_after)  # {11}
            )
    line += '{0}\t{1:.2f}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7:.4f}\n'.format(
                run_count,  # {0}
                tot_run_time,  # {1}
                starting_params[0],  # {2} - use_recomb_rates
                starting_params[1],  # {3} - def_trans_in_p
                starting_params[2],  # {4} - def_emit_same_p
                starting_params[3],  # {5} - final trans_p
                starting_params[4],  # {6} - final emit_p
                final_prob_dist  # {7} - final probability distance
            )

    filename_out = working_dir + '/results/' + filename_in.rsplit('.', 1)[0] + '_stats' + unique_output_name + \
                   '.txt'
    if not isfile(filename_out):
        line = 'Datetime\tFinal_score\tStart_len\tFinal_len\t%_diff\t%_A\t%_AKR\t%_BALBc\t%_C3HHe\t%_C57BL6N\t' + \
               '%_DBA2\t%_Unknown\tRun_count\tTotal_run_time(s)\tUse_recomb_rates\tStart_trans_in\t' + \
               'Start_emit_same\tFinal_trans_p\tFinal_emit_p\n' + line

    with open(filename_out, 'a') as f:
        f.write(line)
