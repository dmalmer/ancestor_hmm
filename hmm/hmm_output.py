
from collections import defaultdict
from datetime import datetime
from os.path import isfile

from hmm_util import ancestor_blocks


# Write ancestor classifications to a .bed file
def write_ancestors_to_file(working_dir, filename_in, unique_output_name, ancestors_by_chr, SNPs_by_chr, state_RGBs):
    out_file = open(working_dir + '/results/' + filename_in.rsplit('.', 1)[0] + '_hmm-out' + unique_output_name +
                    '.' + filename_in.rsplit('.', 1)[1], 'w')
    out_len = 0
    for curr_chr in sorted(ancestors_by_chr.keys()):
        for pos_start, pos_end, ancestor in ancestor_blocks(ancestors_by_chr[curr_chr], SNPs_by_chr[curr_chr]):
            color_key = ancestor
            if '_' in color_key:
                color_key = 'IBA'
            out_file.write(curr_chr + '\t' + pos_start + '\t' + pos_end + '\t' + ancestor + '\t0\t+\t' + pos_start +
                           '\t' + pos_end + '\t' + state_RGBs[color_key] + '\n')
            out_len += 1
    out_file.close()


# Write statistics of each run out to a file
def write_statistics(working_dir, filename_in, unique_output_name, ancestors_by_chr, SNPs_by_chr, starting_params,
                     run_count, tot_run_time):
    len_before = 0
    for SNPs in SNPs_by_chr.values():
        len_before += len(SNPs)
    len_after = 0
    anc_counts = defaultdict(int)

    for curr_chr in ancestors_by_chr.keys():
        for pos_start, pos_end, ancestor in ancestor_blocks(ancestors_by_chr[curr_chr], SNPs_by_chr[curr_chr]):
            anc_counts[ancestor] += 1
            len_after += 1

    line = '{0}\t{1}\t{2}\t{3:.2%}\t{4:.3f}\t{5:.3f}\t{6:.3f}\t{7:.3f}\t{8:.3f}\t{9:.3f}\t{10}\t'.format(
                datetime.now().strftime('%m/%d/%y-%H:%M'),  # {0} - date and time of run
                len_before,  # {1} - input file length
                len_after,  # {2} - output file length
                float(len_before-len_after)/len_before,  # {3} - percentage difference
                anc_counts['A']/float(len_after),  # {4}
                anc_counts['AKR']/float(len_after),  # {5}
                anc_counts['BALBc']/float(len_after),  # {6}
                anc_counts['C3HHe']/float(len_after),  # {7}
                anc_counts['C57BL6N']/float(len_after),  # {8}
                anc_counts['DBA2']/float(len_after),  # {9}
                anc_counts['Unk']/float(len_after)  # {10}
            )
    line += '{0}\t{1:.2f}\t{2}\t{3}\t{4}\t{5}\t{6}\n'.format(
                run_count,  # {0}
                tot_run_time,  # {1}
                starting_params[0],  # {2} - use_recomb_rates
                starting_params[1],  # {3} - def_trans_in_p
                starting_params[2],  # {4} - def_emit_same_p
                starting_params[3],  # {5} - final trans_p
                starting_params[4]  # {6} - final emit_p
            )

    filename_out = working_dir + '/results/' + filename_in.rsplit('.', 1)[0] + '_stats' + unique_output_name + \
                   '.txt'
    if not isfile(filename_out):
        line = 'Datetime\tStart_len\tFinal_len\t%_diff\t%_A\t%_AKR\t%_BALBc\t%_C3HHe\t%_C57BL6N\t%_DBA2\t' + \
               '%_Unknown\tRun_count\tTotal_run_time(s)\tUse_recomb_rates\tStart_trans_in\tStart_emit_same\t' + \
               'Final_trans_p\tFinal_emit_p\n' + line

    with open(filename_out, 'a') as f:
        f.write(line)
