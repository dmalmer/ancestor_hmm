
#run with:
# python bed_plots.py ISS unique
# python bed_plots.py ISS full
# python bed_plots.py unique
# python bed_plots.py full

import sys
import pylab
from collections import defaultdict

if __name__ == "__main__":

    UNIQUE_NAME = '_94in'

    strain = sys.argv[1]

    if sys.argv[2].lower() == 'full':
        snp_type = 'Full'
    elif sys.argv[2].lower() == 'unique':
        snp_type = 'temporary'

    #chr3	3140748	4057786	ARK

    mice_distrb_indv = {} #dictionary of defaultdict's
    mice_distrb_aggr = defaultdict(float)

    chromosomes = [
            'chr1',
            'chr2',
            'chr3',
            'chr4',
            'chr5',
            'chr6',
            'chr7',
            'chr8',
            'chr9',
            'chr10',
            'chr11',
            'chr12',
            'chr13',
            'chr14',
            'chr15',
            'chr16',
            'chr17',
            'chr18',
            'chr19',
            'chrX',
            'chrY'
        ]


    total_recomb = 0
    for cur_chr in chromosomes:
        f_in = open('../results/' + strain + '_' + snp_type + '_sorted_' + cur_chr + '_hmm-out' + UNIQUE_NAME + '.bed','r')
        lines = f_in.readlines()

        mice_distrb_indv[cur_chr] = defaultdict(float)
        for l in lines:
            spl = l.strip().split('\t')
            mice_distrb_indv[cur_chr][spl[3]] += (int(spl[2])-int(spl[1]))

            mice_distrb_aggr[spl[3]] += (int(spl[2])-int(spl[1]))

        total_recomb += len(lines)

    print 'total recombinations: ' + str(total_recomb)

    # Hard code colors to match WashU genome viewer
    mice_colors = {
        'Unk'	    : '#707070', #gray
        'DBA2'      : 'r', #red
        'A'         : 'g', #green
        'ARK'       : 'b', #blue
        'BALBc'     : '#FFFF00', #yellow
        'C3HHe'     : '#FF9900', #orange
        'C57BL6N'   : '#800080' #purple
    }

    # Genome-wide
    pylab.figure(figsize=(6,6))
    pylab.pie(mice_distrb_aggr.values(), explode=None, labels=mice_distrb_aggr.keys(),
              autopct='%1.1f%%', shadow=True, colors=map(lambda x: mice_colors[x], mice_distrb_aggr.keys()))
    pylab.title(sys.argv[2][0].upper() + sys.argv[2][1:].lower() + ' ' + strain + ' - HMM - Genome-Wide Ancestor Distributions')

    # Individual chromosomes
    pylab.figure()
    i = 0
    for chr in chromosomes:
        pylab.subplot(3,7,i)
        i += 1
        pylab.pie(mice_distrb_indv[chr].values(), explode=None, labels=mice_distrb_indv[chr].keys(),
                  autopct='%1.1f%%', shadow=True, colors=map(lambda x: mice_colors[x], mice_distrb_indv[chr].keys()))
        pylab.title(chr)
    pylab.suptitle(sys.argv[2][0].upper() + sys.argv[2][1:].lower() + ' ' + strain + ' - HMM - Chromosome-Wide Ancestor Distributions')
    #pylab.legend(miceColors.keys(),loc=6)

    pylab.show()
