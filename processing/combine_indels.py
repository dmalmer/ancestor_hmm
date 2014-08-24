
from collections import defaultdict
from operator import itemgetter

insertions_by_chr = defaultdict(lambda : defaultdict(str))
deletions_by_chr = defaultdict(lambda : defaultdict(str))

for s in ('A', 'ARK', 'BALBc', 'C3HHe', 'DBA2'):
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

with open('data/ancestor_insertions.bed', 'w') as f:
    for curr_chr in ('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12',
                     'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chrX', 'chrY'):
        sorted_dict = sorted(insertions_by_chr[curr_chr].iteritems(), key=itemgetter(0))
        for item in sorted_dict:
            f.write(curr_chr + '\t' + str(item[0][0]) + '\t' + str(item[0][1]) + '\t' + str(item[1][:-1]) + '\n')

with open('data/ancestor_deletions.bed', 'w') as f:
    for curr_chr in ('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12',
                     'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chrX', 'chrY'):
        sorted_dict = sorted(deletions_by_chr[curr_chr].iteritems(), key=itemgetter(0))
        for item in sorted_dict:
            f.write(curr_chr + '\t' + str(item[0][0]) + '\t' + str(item[0][1]) + '\t' + str(item[1][:-1]) + '\n')
