
#run with:
# python bed_plots.py ISS unique <unique-name>
# python bed_plots.py ILS unique <unique-name>

# python bed_plots.py ISS full <unique-name>
# python bed_plots.py ILS full <unique-name>

# python bed_plots.py ISS sh-nsh <unique-name>
# python bed_plots.py ILS sh-nsh <unique-name>

# python bed_plots.py ISS filt-anc <unique-name>
# python bed_plots.py ILS filt-anc <unique-name>

import sys
import pylab
from collections import defaultdict

if __name__ == "__main__":

    strain = sys.argv[1].upper()
    
    if sys.argv[2].lower() == 'full':
        snp_type = 'Full'
    elif sys.argv[2].lower() == 'unique':
        snp_type = 'temporary'
    elif sys.argv[2].lower() == 'sh-nsh':
        snp_type = 'sh-nsh_full'
    elif sys.argv[2].lower() == 'filt-anc':
        snp_type = 'sh-nsh_filt-anc_full'
    else:
        raise Exception('Bad inputs. See comments at top of script for correct usage.')

    unique_name = sys.argv[3] if len(sys.argv) > 3 else ''

    mice_distrb_indv = {} #dictionary of defaultdict's
    mice_distrb_aggr = defaultdict(float)

    #separate counters for ibd regions
    mice_distrb_indv_ibd = {} #dictionary of defaultdict's
    mice_distrb_aggr_ibd = defaultdict(float)

    chromosomes = []

    total_recomb = 0
    with open('../results/' + strain + '_' + snp_type + '_hmm-out' + unique_name + '.bed','r') as f:
        line = f.readline()

        curr_chr = line.split('\t')[0]
        chromosomes.append(curr_chr)
        mice_distrb_indv[curr_chr] = defaultdict(float)
        mice_distrb_indv_ibd[curr_chr] = defaultdict(float)
        while line != '':
            spl = line.strip().split('\t')
            if spl[0] != curr_chr:
                curr_chr = spl[0]
                chromosomes.append(curr_chr)
                mice_distrb_indv[curr_chr] = defaultdict(float)
                mice_distrb_indv_ibd[curr_chr] = defaultdict(float)

            if '_' in spl[3]:
                mice_distrb_indv[curr_chr]['IBD'] += (int(spl[2])-int(spl[1])) 
                mice_distrb_aggr['IBD'] += (int(spl[2])-int(spl[1]))
                for anc in spl[3].split('_'):
                    mice_distrb_indv_ibd[curr_chr][anc] += (int(spl[2])-int(spl[1])) 
                    mice_distrb_aggr_ibd[anc] += (int(spl[2])-int(spl[1]))
            else:
                mice_distrb_indv[curr_chr][spl[3]] += (int(spl[2])-int(spl[1]))
                mice_distrb_aggr[spl[3]] += (int(spl[2])-int(spl[1]))

            line = f.readline()
            total_recomb += 1

    print 'total recombinations: ' + str(total_recomb)

    # Hard code colors to match WashU genome viewer
    mice_colors = {
        'Unk'       : '#707070', #gray
        'DBA2'      : 'r', #red
        'A'         : 'g', #green
        'AKR'       : '#3366FF', #light blue
        'BALBc'     : '#FFFF00', #yellow
        'C3HHe'     : '#FF9900', #orange
        'C57BL6N'   : '#800080', #purple
        'IBD'       : '#99FFFF' #teal
    }

    # Genome-wide
    #   Non-IBD
    pylab.figure(figsize=(6,6), facecolor='white')
    patches, texts = pylab.pie(mice_distrb_aggr.values(), explode=None, labels=mice_distrb_aggr.keys(), #autopct='%1.1f%%',
              shadow=True, colors=map(lambda x: mice_colors[x], mice_distrb_aggr.keys()), labeldistance=0.65)
    pylab.title(sys.argv[2][0].upper() + sys.argv[2][1:].lower() + ' ' + strain + ' - HMM - Genome-Wide Ancestor Distributions')
    for t in texts:
        t.set_horizontalalignment('center')
    
    #   IBD
    pylab.figure(figsize=(6,6), facecolor='white')
    patches, texts = pylab.pie(mice_distrb_aggr_ibd.values(), explode=None, labels=mice_distrb_aggr_ibd.keys(), #autopct='%1.1f%%',
              shadow=True, colors=map(lambda x: mice_colors[x], mice_distrb_aggr_ibd.keys()), labeldistance=0.65)
    pylab.title(sys.argv[2][0].upper() + sys.argv[2][1:].lower() + ' ' + strain + ' - IBD - HMM - Genome-Wide Ancestor Distributions')
    for t in texts:
        t.set_horizontalalignment('center')

    # Individual chromosomes
    #   Non-IBD
    pylab.figure(facecolor='white')
    i = 1
    for curr_chr in chromosomes:
        pylab.subplot(3,7,i)
        i += 1
        pylab.pie(mice_distrb_indv[curr_chr].values(), explode=None, #labels=mice_distrb_indv[curr_chr].keys(), autopct='%1.1f%%',  
                  shadow=True, colors=map(lambda x: mice_colors[x], mice_distrb_indv[curr_chr].keys()))
        pylab.title(curr_chr)
    pylab.suptitle(sys.argv[2][0].upper() + sys.argv[2][1:].lower() + ' ' + strain + ' - HMM - Chromosome-Wide Ancestor Distributions')
    #pylab.legend(miceColors.keys(),loc=6)
    
    #   IBD
    pylab.figure(facecolor='white')
    i = 1
    for curr_chr in chromosomes:
        pylab.subplot(3,7,i)
        i += 1
        pylab.pie(mice_distrb_indv_ibd[curr_chr].values(), explode=None, #labels=mice_distrb_indv_ibd[curr_chr].keys(), autopct='%1.1f%%',  
                  shadow=True, colors=map(lambda x: mice_colors[x], mice_distrb_indv_ibd[curr_chr].keys()))
        pylab.title(curr_chr)
    pylab.suptitle(sys.argv[2][0].upper() + sys.argv[2][1:].lower() + ' ' + strain + ' - IBD - HMM - Chromosome-Wide Ancestor Distributions')
    #pylab.legend(miceColors.keys(),loc=6)

    pylab.show()
