
from collections import defaultdict
from os.path import isfile

from hmm_util import ancestor_blocks

import os

def print_ancestors(ancestors, SNPs, title):
    print ''
    print title
    cur = ''
    for i in range(len(ancestors)):
        if cur != ancestors[i]:
            print str(SNPs[i][1]) + ': ' + ancestors[i]
            cur = ancestors[i]


def write_ancestors_to_file(WORKING_DIR, filename_in, unique_output_name, ancestors, SNPs, state_RGBs):
    out_file = open(WORKING_DIR + '/results/' + filename_in.rsplit('.', 1)[0] + '_hmm-out' + unique_output_name +
                    '.' + filename_in.rsplit('.', 1)[1], 'w')

    out_len = 0
    for chromosome, pos_start, pos_end, ancestor in ancestor_blocks(ancestors, SNPs):
        out_file.write(chromosome + '\t' + pos_start + '\t' + pos_end + '\t' + ancestor + '\t0\t+\t' + pos_start + '\t' +
                        pos_end + '\t' + state_RGBs[ancestor] + '\n')
        out_len += 1
    out_file.close()

    print '\nOutput file length: ' + str(out_len)
    print 'Percentage change: ' + str(float(len(SNPs)-out_len)/float(len(SNPs)))


def write_confidence_interval(WORKING_DIR, filename_in, unique_output_name, confidence_intervals):
    out_file = open(WORKING_DIR + '/results/' + filename_in.rsplit('.', 1)[0] + '_conf-int' + unique_output_name +
                    '.wig', 'w')
    out_file.write('track type=wiggle_0 graphType=line viewLimits=0:1\n')

    cur_chrom = ''
    for next_chrom, pos_start, conf in confidence_intervals:
        if next_chrom != cur_chrom:
            cur_chrom = next_chrom
            out_file.write('variableStep chrom=%s\n' % cur_chrom)
        out_file.write('%s\t%f\n' % (pos_start, conf))


def write_indentical_by_ancestor(WORKING_DIR, filename_in, unique_output_name, identical_by_anc):
    out_file = open(WORKING_DIR + '/results/' + filename_in.rsplit('.', 1)[0] + '_IBA' + unique_output_name +
                    '.bed', 'w')

    for chromosome, pos_start, pos_end, iba in identical_by_anc:
        out_file.write('%s\t%s\t%s\t%s\n' % (chromosome, pos_start, pos_end, iba))


def write_statistics(WORKING_DIR, filename_in, unique_output_name, ancestors, SNPs, starting_params, run_count, tot_run_time):
    len_before = len(SNPs)
    len_after = 0
    anc_counts = defaultdict(int)

    for chromosome, pos_start, pos_end, ancestor in ancestor_blocks(ancestors, SNPs):
        anc_counts[ancestor] += 1
        len_after += 1

    line = '%i\t%i\t%.3f\t%.5f\t%.5f\t%.2f\t%.2f\t%.1f\t%s\t%s\t%s\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%i\t%i\n' % \
           (len_before, len_after, float(len_before-len_after)/len_before, starting_params[0], starting_params[1],
            starting_params[2], starting_params[3], starting_params[4], starting_params[5], starting_params[6],
            starting_params[7], anc_counts['A']/float(len_after), anc_counts['ARK']/float(len_after),
            anc_counts['BALBc']/float(len_after), anc_counts['C3HHe']/float(len_after),
            anc_counts['C57BL6N']/float(len_after), anc_counts['DBA2']/float(len_after),
            anc_counts['Unk']/float(len_after), run_count, tot_run_time)

    filename_out = WORKING_DIR + '/results/' + filename_in.rsplit('.', 1)[0] + '_stats' + unique_output_name + \
                    '.txt'
    if not isfile(filename_out):
        line = 'Start_len\tFinal_len\t%_diff\tTrans_in\tTrans_out\tEmit_same\tEmit_other\tHS_FI\tUse_HS\tUse_SNP_dist\t + \
        Use_recomb_rate\t%_A\t%_ARK\t%_BALBc\t%_C3HHe\t%_C57BL6N\t%_DBA2\t%_Unknown\tRun_count\tTotal_run_time(s)\n' + line

    with open(filename_out, 'a') as f:
        f.write(line)
