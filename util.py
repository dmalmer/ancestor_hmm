
import re
from collections import defaultdict
from datetime import datetime
from itertools import tee, izip
from math import e, log
from numpy import arange, array
from os.path import isfile


# Generator to loop over blocks of sequential ancestors
def ancestor_blocks(ancestors, SNPs, return_SNPs=False):
    curr_i = 0
    while curr_i < len(ancestors):
        start_i = curr_i
        curr_i += 1
        while curr_i < len(ancestors) and ancestors[curr_i-1] == ancestors[curr_i]:
            curr_i += 1

        if return_SNPs:
            yield (SNPs[start_i, 1], SNPs[curr_i-1, 2], ancestors[start_i], SNPs[start_i:curr_i])
        else:
            yield (SNPs[start_i, 1], SNPs[curr_i-1, 2], ancestors[start_i])


# Convert string to int (if possible)
def atoi(text):
    try:
        return int(text)
    except ValueError:
        return text


# Convert string to float (if possible)
def atof(text):
    try:
        return float(text)
    except ValueError:
        return text


# Create range of values for grid search
def create_grid_range(input_params, grid_size):
    start, stop = [float(p) for p in input_params.strip('[](){}<>').split('-')]
    grid_range = list(arange(start, stop, (stop-start)/(grid_size-1)))
    grid_range.append(stop)

    return grid_range


# Return the emit key
def get_emit_key(state, SNPs_str, desc_strain):
    if state == 'Unk' and SNPs_str == desc_strain:
        return state
    SNPs = SNPs_str.split('_')
    if (desc_strain in SNPs and state in SNPs) or (desc_strain not in SNPs and state not in SNPs):
        return state
    else:
        return '~' + state


# Get all unique states from SNP data
def get_states(SNPs_by_chr, desc_strain, use_unk):
    states = []
    for v in SNPs_by_chr.values():
        for SNPs in v:
            for anc in SNPs[3].split('_'):
                if anc != desc_strain and anc not in states:
                    states.append(anc)

    if use_unk:
        states.append('Unk')
    
    return states


# Find log(A+B) when A and B are in log-space
#  (taken from https://facwiki.cs.byu.edu/nlp/index.php/Log_Domain_Computations)
def log_add_pair(log_A, log_B):
    # ***Note***: log_A needs to be greater than log_B, but I'm removing setting the max here to speed things up
    return log_A + log(1. + e ** (log_B - log_A))


# Find sum of a list of numbers in log-space
def log_add_list(log_list):
    log_list.sort(reverse=True)
    log_sum = log_list[0]
    for l in log_list[1:]:
        log_sum = log_add_pair(log_sum, l)
    return log_sum


# Key to sort by numerics within strings (ie. chr2 comes before chr11, etc.)
#  (taken from http://nedbatchelder.com/blog/200712/human_sorting.html)
def natural_keys(text):
    # Use: alist.sort(key=natural_keys) sorts in human order
    return [ atoi(c) for c in re.split('(\d+)', text) ]


# Return a list of pairwise elements
#  (taken from https://docs.python.org/2/library/itertools.html#recipes)
def pairwise(iterable):
    # s -> (s0,s1), (s1,s2), (s2, s3), ...
    a, b = tee(iterable)
    next(b, None)
    return izip(a, b)


# Return dictionary of probability tuples to be used in print statements
def prob_tuples(probs):
    prob_dict = defaultdict(list)
    for s in probs.keys():
        for k, v in probs[s].items():
            prob_dict[s].append((k, e ** v))

    return prob_dict


# Read in recombination rates data
def read_recomb_rates(filename):
    # Divide recombination rate data into dictionary by chromosome
    recomb_rate_dict = defaultdict(list)
    with open(filename, 'r') as f:
        f.readline()
        line = f.readline()
        while line != '':
            splits = line.strip().split(',')
            # Convert DNA position from Kb to b
            recomb_rate_dict[splits[0]].append([float(splits[1])*1000, float(splits[2])])
            line = f.readline()

    # Convert to numpy arrays
    for k, v in recomb_rate_dict.items():
        recomb_rate_dict[k] = array(v, dtype=float)

    return recomb_rate_dict


# Read in SNP data from .bed file
def read_SNPs(filename):
    # Divide SNP data into dictionary by chromosome
    SNPs_dict = defaultdict(list)
    with open(filename, 'r') as f:
        line = f.readline().strip()
        while line != '':
            cols = line.split('\t')
            SNPs_dict[cols[0]].append(cols)

            line = f.readline().strip()

    # Convert to numpy arrays
    for k, v in SNPs_dict.items():
        SNPs_dict[k] = array(v)

    return SNPs_dict


# Read in structural variant files
def read_SVs(anc_ins_filename, anc_del_filename, states, desc_strain):
    # Divide SV data into dictionaries by chromosome
    #  Read individual insertions, separate descendant insertions from ancestor insertions
    desc_ins_by_chr = defaultdict(list)
    anc_ins_by_chr = defaultdict(list)
    with open(anc_ins_filename, 'r') as f:
        line = f.readline().strip()
        while line != '':
            cols = line.split('\t')

            ancestors = []
            has_desc = False
            for s in cols[3].split('_'):
                if s in states:
                    ancestors.append(s)
                elif s == desc_strain:
                    has_desc = True
            if len(ancestors) > 0:
                anc_ins_by_chr[cols[0]].append((int(cols[1]), int(cols[2]), ancestors))
            if has_desc:
                desc_ins_by_chr[cols[0]].append((int(cols[1]), int(cols[2])))

            line = f.readline().strip()

    #  Sort by start pos
    for insertions in desc_ins_by_chr.values():
        insertions.sort(key=lambda x: x[0])
    for insertions in anc_ins_by_chr.values():
        insertions.sort(key=lambda x: x[0])

    #  Read individual deletions, separate descendant deletions from ancestor deletions
    desc_del_by_chr = defaultdict(list)
    anc_del_by_chr = defaultdict(list)
    with open(anc_del_filename, 'r') as f:
        line = f.readline().strip()
        while line != '':
            cols = line.split('\t')

            ancestors = []
            has_desc = False
            for s in cols[3].split('_'):
                if s in states:
                    ancestors.append(s)
                elif s == desc_strain:
                    has_desc = True
            if len(ancestors) > 0:
                anc_del_by_chr[cols[0]].append((int(cols[1]), int(cols[2]), ancestors))
            if has_desc:
                desc_del_by_chr[cols[0]].append((int(cols[1]), int(cols[2])))

            line = f.readline().strip()

    #  Sort by start pos
    for deletions in desc_del_by_chr.values():
        deletions.sort(key=lambda x: x[0])
    for deletions in anc_del_by_chr.values():
        deletions.sort(key=lambda x: x[0])

    return (desc_ins_by_chr, desc_del_by_chr, anc_ins_by_chr, anc_del_by_chr)


# Write ancestor classifications to a .bed file
def write_ancestors(working_dir, filename_in, append_str, ancestors_by_chr, SNPs_by_chr, state_RGBs):
    extension = '' if '.' not in filename_in else '.' + filename_in.rsplit('.', 1)[1]
    with open(working_dir + filename_in.rsplit('.', 1)[0] + '_hmm-out' + append_str + extension, 'w') as f_out:
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
def write_scores(working_dir, filename_in, append_str, all_scores_by_chr):
    extension = '' if '.' not in filename_in else '.' + filename_in.rsplit('.', 1)[1]
    with open(working_dir + filename_in.rsplit('.', 1)[0] + '_scores' + append_str + extension, 'w') as f_out:
        for curr_chr in sorted(all_scores_by_chr.keys(), key=natural_keys):
            for score in all_scores_by_chr[curr_chr]:
                color = '51,255,51' if score[3] == 'Hit' else '255,51,51'
                f_out.write('{0}\t{1}\t{2}\t{3}\t0\t+\t{1}\t{2}\t{4}\n'.format(
                    curr_chr,  # {0} - chromosome
                    score[0],  # {1} - start pos
                    score[1],  # {2} - end pos
                    score[3],  # {3} - 'Hit' or 'Miss'
                    color  # {4} - green for hit, red for miss
                ))


# Write statistics of each run out to a file
def write_statistics(working_dir, filename_in, append_str, ancestors_by_chr, SNPs_by_chr, starting_params,
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

    line = '{0}\t{1:.3f}\t{2:.4f}\t{3}\t{4}\t{5:.2%}\t'.format(
                datetime.now().strftime('%m/%d/%y-%H:%M'),  # {0} - date and time of run
                final_score, # {1} - final score
                final_prob_dist, # {2} - final probability distance
                len_before,  # {3} - input file length
                len_after,  # {4} - output file length
                float(len_before-len_after)/len_before,  # {5} - percentage difference
            )
    line += '{0:.3f}\t{1:.3f}\t{2:.3f}\t{3:.3f}\t{4:.3f}\t{5:.3f}\t{6:.3f}\t'.format(
                anc_counts['A']/float(len_after),  # {0}
                anc_counts['AKR']/float(len_after),  # {1}
                anc_counts['BALBc']/float(len_after),  # {2}
                anc_counts['C3HHe']/float(len_after),  # {3}
                anc_counts['C57BL6N']/float(len_after),  # {4}
                anc_counts['DBA2']/float(len_after),  # {5}
                anc_counts['Unk']/float(len_after)  # {6}
            )
    line += '{0}\t{1:.2f}\t{2}\t{3}\t{4}\t{5}\t{6}\n'.format(
                run_count,  # {0}
                tot_run_time,  # {1}
                starting_params[0],  # {2} - use_recomb_rates
                starting_params[1],  # {3} - def_trans_in_p
                starting_params[2],  # {4} - def_emit_same_p
                starting_params[3],  # {5} - final trans_p
                starting_params[4],  # {6} - final emit_p
            )

    filename_out = working_dir + filename_in.rsplit('.', 1)[0] + '_stats' + append_str + '.txt'
    if not isfile(filename_out):
        header = '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\t{15}\t{16}\t{17}\t' \
                 '{18}\t{19}\n'.format(
                    'Datetime',
                    'Final_score',
                    'Final_prob_dist',
                    'Start_len',
                    'Final_len',
                    '%_diff',
                    '%_A',
                    '%_AKR',
                    '%_BALBc',
                    '%_C3HHe',
                    '%_C57BL6N',
                    '%_DBA2',
                    '%_Unknown',
                    'Run_count',
                    'Total_run_time(s)',
                    'Use_recomb_rates',
                    'Start_trans_in',
                    'Start_emit_same',
                    'Final_trans_p',
                    'Final_emit_p'
                )
        line = header + line

    with open(filename_out, 'a') as f:
        f.write(line)
