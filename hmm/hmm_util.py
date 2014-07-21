
from collections import defaultdict
from itertools import tee, izip
from math import ceil, e, log
from numpy import array


# Generator to loop over only unique ancestors
def ancestor_blocks(ancestors, SNPs, prob_nodes=[]):
    curr_i = 0
    while curr_i < len(ancestors):
        start_i = curr_i
        curr_i += 1
        while curr_i < len(ancestors) and ancestors[curr_i-1] == ancestors[curr_i]:
            curr_i += 1

        if any(prob_nodes):
            yield (SNPs[start_i,0], SNPs[start_i,1], SNPs[curr_i-1,2], ancestors[start_i], prob_nodes[start_i:curr_i])
        else:
            yield (SNPs[start_i,0], SNPs[start_i,1], SNPs[curr_i-1,2], ancestors[start_i])


def calc_recomb_rate(SNP_start, SNP_end, recomb_main_i, recomb_map, effective_pop, num_generations):
    # Find recomb_map starting position
    if recomb_main_i is None:
        recomb_main_i = 0
        while recomb_main_i < len(recomb_map) and int(recomb_map[recomb_main_i][0]) < SNP_start:
            recomb_main_i += 1
        recomb_start_i = max(recomb_main_i-1, 0)
    else:
        recomb_start_i = recomb_main_i - 1

    # Quick check to make sure recomb_start_i >= 0 (should be as recomb_index should always be >=1 in else statement)
    #  remove this later
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
        expected_recombs = ((((recomb_map[recomb_start_i+1][0] - SNP_start) / 1000.) * recomb_map[recomb_start_i][1]) / \
                            (4 * effective_pop)) * num_generations

        # Rates in the middle
        for i in range(recomb_start_i+1, recomb_end_i-1):
            expected_recombs += ((((recomb_map[i+1][0] - recomb_map[i][0]) / 1000.) * recomb_map[i][1]) / \
                                 (4 * effective_pop)) * num_generations

        # Proportional rate of second SNP
        expected_recombs += ((((SNP_end - recomb_map[recomb_end_i-1][0]) / 1000.) * recomb_map[recomb_end_i-1][1]) / \
                             (4 * effective_pop)) * num_generations

    return expected_recombs, recomb_main_i


# Count the number of hotspots between two chromosome positions
def count_hotspots(chromosome, pos_start, pos_end, hotspot_dict):
    # All even indexes are the pos start of cold (not hot) spots, all odd indexes are the pos start of hot spots
    i_start = len(hotspot_dict[chromosome][pos_start-hotspot_dict[chromosome] > 0]) - 1
    i_end = len(hotspot_dict[chromosome][pos_end-hotspot_dict[chromosome] > 0]) - 1

    hs_count = ceil((i_end-i_start)/2.)

    # if both the start and end indexes are on hotspots, we need to add 1 more to the count
    if i_start % 2 == 1 and i_end % 2 == 1:
        hs_count += 1

    return hs_count


# Find log(A+B) when A and B are in log-space
#  (taken from https://facwiki.cs.byu.edu/nlp/index.php/Log_Domain_Computations)
def log_add_pair(log_A, log_B):
    # ***Note***: log_A needs to be greater than log_B, but I'm removing setting the max here to speed things up
    return log_A + log(1. + e ** (log_B - log_A))


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


# Read in hotspot data
def read_hotspots_data(filename):
    # defaultdict list that begins with a 0 element rather than starting empty
    # All even indexes are the pos start of cold (not hot) spots, all odd indexes are the pos start of hot spots
    hotspot_dict = defaultdict(lambda: [0])
    with open(filename, 'r') as f:
        f.readline()
        line = f.readline()
        while line != '':
            splits = line.split(',')
            hotspot_dict['chr'+str(splits[0])].extend([int(splits[1]), int(splits[2])+1])
            line = f.readline()
    # Convert to numpy arrays
    for k, v in hotspot_dict.items():
        hotspot_dict[k] = array(v)

    return hotspot_dict


# Read in recombination rates data
def read_recomb_rates_data(filename):
    recomb_rate_dict = defaultdict(list)
    with open(filename, 'r') as f:
        f.readline()
        line = f.readline()
        while line != '':
            splits = line.strip().split(',')
            recomb_rate_dict[splits[0]].append([float(splits[1])*1000, float(splits[2])])
            line = f.readline()

    # Convert to numpy arrays
    for k, v in recomb_rate_dict.items():
        recomb_rate_dict[k] = array(v, dtype=float)

    return recomb_rate_dict


