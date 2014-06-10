
from collections import defaultdict
from os.path import isfile

from hmm_util import ancestor_blocks


def print_ancestors(ancestors, SNPs, title):
    print ''
    print title
    cur = ''
    for i in range(len(ancestors)):
        if cur != ancestors[i]:
            print str(SNPs[i][1]) + ': ' + ancestors[i]
            cur = ancestors[i]


def write_ancestors_to_file(filename_in, unique_output_name, ancestors, SNPs):
    extension = filename_in.rsplit('.', 1)[1]
    out_file = open(filename_in.split('.' + extension)[0] + '_hmm-out' + unique_output_name + '.' + extension, 'w')

    out_len = 0
    for chromosome, pos_start, pos_end, ancestor in ancestor_blocks(ancestors, SNPs):
        out_file.write(chromosome + '\t' + pos_start + '\t' + pos_end + '\t' + ancestor + '\n')
        out_len += 1
    out_file.close()

    print '\nOutput file length: ' + str(out_len)
    print 'Percentage change: ' + str(float(len(SNPs)-out_len)/float(len(SNPs)))


def write_statistics(filename_in, ancestors, SNPs, starting_params, run_count, final_hs_fi):
    len_before = len(SNPs)
    len_after = 0
    anc_counts = defaultdict(int)

    for chromosome, pos_start, pos_end, ancestor in ancestor_blocks(ancestors, SNPs):
        anc_counts[ancestor] += 1
        len_after += 1

    line = '%i\t%i\t%.3f\t%.5f\t%.5f\t%.2f\t%.2f\t%.1f\t%.1f\t%s\t%s\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%i\n' % \
           (len_before, len_after, float(len_before-len_after)/len_before, starting_params[0], starting_params[1],
            starting_params[2], starting_params[3], starting_params[4], starting_params[5], starting_params[6],
            final_hs_fi, anc_counts['A']/float(len_after), anc_counts['ARK']/float(len_after),
            anc_counts['BALBc']/float(len_after), anc_counts['C3HHe']/float(len_after),
            anc_counts['C57BL6N']/float(len_after), anc_counts['DBA2']/float(len_after),
            anc_counts['Unk']/float(len_after), run_count)

    filename_out = '/'.join(filename_in.split('/')[:-1]) + '/STATS_' + filename_in.split('/')[-1]
    if not isfile(filename_out):
        line = 'Start_len\tFinal_len\t%_diff\tTrans_in\tTrans_out\tEmit_same\tEmit_other\tHS_FI\tUse_HS\tUse_SNP_dist\t + \
        %_A\t%_ARK\t%_BALBc\t%_C3HHe\t%_C57BL6N\t%_DBA2\t%_Unknown\tRun_count\n' + line

    with open(filename_out, 'a') as f:
        f.write(line)