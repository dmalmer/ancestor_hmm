
from collections import defaultdict
from operator import itemgetter


def clean_label(label):
    return '_'.join(sorted(list(set(label[:-1].split('_')))))


insertions_by_chr = defaultdict(lambda : defaultdict(str))
deletions_by_chr = defaultdict(lambda : defaultdict(str))

# read in anc indels
for s in ('A', 'AKR', 'BALBc', 'C3HHe', 'DBA2'):
    with open('data/' + s + '_insertions.bed', 'r') as f:
        line = f.readline().strip()
        while line != '':
            cols = line.split('\t')
            insertions_by_chr[cols[0]][(int(cols[1]), int(cols[2]))] += s + '_'

            line = f.readline().strip()

    with open('data/' + s + '_deletions.bed', 'r') as f:
        line = f.readline().strip()
        while line != '':
            cols = line.split('\t')
            deletions_by_chr[cols[0]][(int(cols[1]), int(cols[2]))] += s + '_'

            line = f.readline().strip()

# write to new file, merging overlaps
with open('data/ancestor_insertions.bed', 'w') as f:
    for curr_chr in ('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12',
                     'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chrX', 'chrY'):
        sorted_dict = sorted(insertions_by_chr[curr_chr].iteritems(), key=itemgetter(0))
        i = 0
        while i < len(sorted_dict):
            pos_start = sorted_dict[i][0][0]
            pos_end = sorted_dict[i][0][1]
            label = sorted_dict[i][1]
            i += 1
            while i < len(sorted_dict) and sorted_dict[i][0][0] < pos_end:
                pos_start = min(sorted_dict[i][0][0], pos_start)
                pos_end = max(sorted_dict[i][0][1], pos_end)
                label += sorted_dict[i][1]
                i += 1
            f.write(curr_chr + '\t' + str(pos_start) + '\t' + str(pos_end) + '\t' + clean_label(label) + '\n')

with open('data/ancestor_deletions.bed', 'w') as f:
    for curr_chr in ('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12',
                     'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chrX', 'chrY'):
        sorted_dict = sorted(deletions_by_chr[curr_chr].iteritems(), key=itemgetter(0))
        i = 0
        while i < len(sorted_dict):
            pos_start = sorted_dict[i][0][0]
            pos_end = sorted_dict[i][0][1]
            label = sorted_dict[i][1]
            i += 1
            while i < len(sorted_dict) and sorted_dict[i][0][0] < pos_end:
                pos_start = min(sorted_dict[i][0][0], pos_start)
                pos_end = max(sorted_dict[i][0][1], pos_end)
                label += sorted_dict[i][1]
                i += 1
            f.write(curr_chr + '\t' + str(pos_start) + '\t' + str(pos_end) + '\t' + clean_label(label) + '\n')
