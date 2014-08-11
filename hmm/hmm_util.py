
from collections import defaultdict
from itertools import tee, izip
from math import e, log
from numpy import array


# Generator to loop over blocks of sequential ancestors
def ancestor_blocks(ancestors, SNPs, return_SNPs=False):
    curr_i = 0
    while curr_i < len(ancestors):
        start_i = curr_i
        curr_i += 1
        while curr_i < len(ancestors) and ancestors[curr_i-1] == ancestors[curr_i]:
            curr_i += 1

        if return_SNPs:
            yield (SNPs[start_i, 0], SNPs[start_i, 1], SNPs[curr_i-1, 2], ancestors[start_i], SNPs[start_i:curr_i])
        else:
            yield (SNPs[start_i, 0], SNPs[start_i, 1], SNPs[curr_i-1, 2], ancestors[start_i])


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


# Return a list of pairwise elements
#  (taken from https://docs.python.org/2/library/itertools.html#recipes)
def pairwise(iterable):
    # s -> (s0,s1), (s1,s2), (s2, s3), ...
    a, b = tee(iterable)
    next(b, None)
    return izip(a, b)


# Return dictionary of probability tuples to be used in print statements
def prob_tuples(probs):
    ret_dict = {}
    for s in probs.keys():
        prob_list = [s]
        for k, v in probs[s].items():
            prob_list.extend([k, e ** v])
        ret_dict[s] = tuple(prob_list)

    return ret_dict


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
