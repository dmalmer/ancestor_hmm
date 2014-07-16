
#run with:
# python bed_stats.py ISS unique
# python bed_stats.py ISS full
# python bed_stats.py ILS unique
# python bed_stats.py ILS full

import sys

if __name__ == "__main__":

    UNIQUE_NAME = '_94in'

    strain = sys.argv[1]

    if sys.argv[2].lower() == 'full':
        snp_type = 'Full'
    elif sys.argv[2].lower() == 'unique':
        snp_type = 'temporary'

    total = 0

    print 'Classifications per chromosome:'

    #chr 1 through 19
    for i in range(1,20):
        cur_chr = 'chr'+str(i)

        f_in = open('../results/' + strain + '_' + snp_type + '_sorted_' + cur_chr + '_hmm-out' + UNIQUE_NAME + '.bed', 'r')

        lines = f_in.readlines()
        print cur_chr + ': ' + str(len(lines))

        total += len(lines)

    #chrX
    cur_chr = 'chrX'

    f_in = open('../results/' + strain + '_' + snp_type + '_sorted_' + cur_chr + '_hmm-out' + UNIQUE_NAME + '.bed', 'r')

    lines = f_in.readlines()
    print cur_chr + ': ' + str(len(lines))

    total += len(lines)

    #chrY
    cur_chr = 'chrY'

    f_in = open('../results/' + strain + '_' + snp_type + '_sorted_' + cur_chr + '_hmm-out' + UNIQUE_NAME + '.bed', 'r')

    lines = f_in.readlines()
    print cur_chr + ': ' + str(len(lines))

    total += len(lines)

    print '\nTotal classifications: ' + str(total)
