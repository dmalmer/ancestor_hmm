
from itertools import tee, izip
from math import ceil, e, log


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
def log_add(log_A, log_B):
    log_max = max(log_A, log_B)
    log_min = min(log_A, log_B)
    return log_max + log(1. + e ** (log_min - log_max))

# Return a list of pairwise elements
#  (taken from https://docs.python.org/2/library/itertools.html#recipes)
def pairwise(iterable):
    # s -> (s0,s1), (s1,s2), (s2, s3), ...
    a, b = tee(iterable)
    next(b, None)
    return izip(a, b)


# Generator to loop over only unique ancestors
def unique_ancestors(ancestors, SNPs):
    curr_i = 0
    while curr_i < len(ancestors):
        start_i = curr_i
        curr_i += 1
        while curr_i < len(ancestors) and ancestors[curr_i-1] == ancestors[curr_i]:
            curr_i += 1

        yield (SNPs[start_i,0], SNPs[start_i,1], SNPs[curr_i-1,2], ancestors[start_i])
