
#run with:
# python create_beds_by_ancestor.py ISS unique
# python create_beds_by_ancestor.py ILS unique

# python create_beds_by_ancestor.py ISS full
# python create_beds_by_ancestor.py ILS full

# python create_beds_by_ancestor.py ISS sh-nsh
# python create_beds_by_ancestor.py ILS sh-nsh

# python create_beds_by_ancestor.py ISS filt-anc
# python create_beds_by_ancestor.py ILS filt-anc


import sys

if __name__ == "__main__":
    strain = sys.argv[1]

    print sys.argv[1]
    print sys.argv[2]

    if sys.argv[2].lower() == 'full':
        snp_type = 'Full'
    elif sys.argv[2].lower() == 'unique':
        snp_type = 'temporary'
    elif sys.argv[2].lower() == 'sh-nsh':
        snp_type = 'sh-nsh_full'
    elif sys.argv[2].lower() == 'filt-anc':
        snp_type = 'sh-nsh_filt-anc_full'
    elif sys.argv[2].lower() == 'rmc57':
        snp_type = 'sh-nsh_filt-anc_rmC57_full'

    anc_RGBs = {strain: '128,128,128', 'A': '0,153,0', 'AKR': '51,102,255', 'BALBc': '255,255,51', 'C3HHe': '255,153,51',
                'C57BL6N': '102,0,204', 'DBA2': '255,51,51', 'not_shared': '0,0,0'}

    f_in = open('../data/' + strain + '_' + snp_type + '.bed', 'r')
    f_outs = {
               strain:    open('../igv/' + strain + '_' + snp_type + '_' + strain + '_sorted_black-not-shared.bed', 'w'),
               'A':       open('../igv/' + strain + '_' + snp_type + '_' + 'A' + '_sorted_black-not-shared.bed', 'w'),
               'AKR':     open('../igv/' + strain + '_' + snp_type + '_' + 'AKR' + '_sorted_black-not-shared.bed', 'w'),
               'BALBc':   open('../igv/' + strain + '_' + snp_type + '_' + 'BALBc' + '_sorted_black-not-shared.bed', 'w'),
               'C3HHe':   open('../igv/' + strain + '_' + snp_type + '_' + 'C3HHe' + '_sorted_black-not-shared.bed', 'w'),
               'C57BL6N': open('../igv/' + strain + '_' + snp_type + '_' + 'C57BL6N' + '_sorted_black-not-shared.bed', 'w'),
               'DBA2':    open('../igv/' + strain + '_' + snp_type + '_' + 'DBA2' + '_sorted_black-not-shared.bed', 'w')
            }

    line = f_in.readline().strip()
    while line != '':
        cols = line.split('\t')
        ancestors = cols[3].split('_')

        if strain in ancestors:
            for anc in ancestors:
                f_outs[anc].write('\t'.join(cols[0:3]) + '\t' + anc + '\t0\t+\t' + cols[1] + '\t' + cols[2] + '\t'
                                 + anc_RGBs[anc] + '\n')
        else:
            for anc in ancestors:
                f_outs[anc].write('\t'.join(cols[0:3]) + '\t' + anc + '\t0\t+\t' + cols[1] + '\t' + cols[2] + '\t'
                                 + anc_RGBs['not_shared'] + '\n')

        line = f_in.readline().strip()
