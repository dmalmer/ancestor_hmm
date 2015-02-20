
import re
from collections import defaultdict
from itertools import tee, izip
from math import e, log
from numpy import arange, array


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
def get_emit_key(state, SNPs_str, input_strain):
    if state == 'Unk' and SNPs_str == input_strain:
        return state
    SNPs = SNPs_str.split('_')
    if (input_strain in SNPs and state in SNPs) or (input_strain not in SNPs and state not in SNPs):
        return state
    else:
        return '~' + state


# Get all unique states from SNP data
def get_states(SNPs_by_chr, input_strain, use_unk):
    states = []
    for v in SNPs_by_chr.values():
        for SNPs in v:
            for anc in SNPs[3].split('_'):
                if anc != input_strain and anc not in states:
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
def read_SVs(strain_SV_filename, anc_ins_filename, anc_del_filename):
    # Divide SV data into dictionaries by chromosome
    #  Read ISS or ILS SVs
    strain_SVs_by_chr = defaultdict(list)
    with open(strain_SV_filename, 'r') as f:
        f.readline()
        line = f.readline()
        while line != '':
            cols = line.split('\t')
            strain_SVs_by_chr[cols[0]].append([int(cols[1]), int(cols[2]), cols[3]])

            line = f.readline()

    for SVs in strain_SVs_by_chr.values():
        SVs.sort()

    #  Read individual ancestor insertions
    anc_ins_by_chr = defaultdict(list)
    with open(anc_ins_filename, 'r') as f:
        line = f.readline().strip()
        while line != '':
            cols = line.split('\t')
            anc_ins_by_chr[cols[0]].append([int(cols[1]), int(cols[2]), cols[3]])

            line = f.readline().strip()

    for insertions in anc_ins_by_chr.values():
        insertions.sort()

    #  Read individual ancestor deletions
    anc_del_by_chr = defaultdict(list)
    with open(anc_del_filename, 'r') as f:
        line = f.readline().strip()
        while line != '':
            cols = line.split('\t')
            anc_del_by_chr[cols[0]].append([int(cols[1]), int(cols[2]), cols[3]])

            line = f.readline().strip()

    for deletions in anc_del_by_chr.values():
        deletions.sort()

    return (strain_SVs_by_chr, anc_ins_by_chr, anc_del_by_chr)
