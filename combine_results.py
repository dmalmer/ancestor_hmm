
#run with:
# python combine_results.py ISS unique
# python combine_results.py ISS full
# python combine_results.py ILS unique
# python combine_results.py ILS full

import sys

if __name__ == "__main__":
    strain = sys.argv[1]

    if sys.argv[2].lower() == 'full':
        snp_type = 'Full'
    elif sys.argv[2].lower() == 'unique':
        snp_type = 'temporary'

    f_out = open('results/' + strain + '_' + snp_type + '_sorted_hmm-out.bed', 'w')

    #chr 1 through 19
    for i in range(1,20):
        f_in = open('results/' + strain + '_' + snp_type + '_sorted_chr' + str(i) + '_hmm-out.bed', 'r')
        f_out.write(f_in.read())

    #chrX
    f_in = open('results/' + strain + '_' + snp_type + '_sorted_chrX_hmm-out.bed', 'r')
    f_out.write(f_in.read())

    #chrY
    f_in = open('results/' + strain + '_' + snp_type + '_sorted_chrY_hmm-out.bed', 'r')
    f_out.write(f_in.read())