
import pylab
import matplotlib
matplotlib.use('Agg')
from collections import defaultdict

QTL_regions = defaultdict(list)
with open('../data/LORE_QTL.bed', 'r') as f:
    for line in f:
        cols = line.strip().split('\t')
        QTL_regions[cols[0]].append((int(cols[1]), int(cols[2]), cols[0] + '_' + cols[1][0] + '-' + cols[2][0]))

QTL_dist = defaultdict(int)
QTL_per_reg = defaultdict(lambda: defaultdict(int))
with open('../results/ILS_sh-nsh_filt-anc_rmC57_full_hmm-out_rmC57t64u9a15pw.bed', 'r') as f:
    for line in f:
        cols = line.strip().split('\t')
        res_s = int(cols[1])
        res_e = int(cols[2])
        for qtl_s, qtl_e, reg in QTL_regions[cols[0]]:
            if res_s > qtl_s and res_e < qtl_e:
                for anc in cols[3].split('_'):
                    QTL_dist[anc] += res_e - res_s
                    QTL_per_reg[reg][anc] += res_e - res_s
            elif res_s < qtl_s and (res_e > qtl_s and res_e < qtl_e):
                for anc in cols[3].split('_'):
                    QTL_dist[anc] += res_e - qtl_s
                    QTL_per_reg[reg][anc] += res_e - qtl_s
            elif (res_s > qtl_s and res_s < qtl_e) and res_e > qtl_e:
                for anc in cols[3].split('_'):
                    QTL_dist[anc] += qtl_e - res_s
                    QTL_per_reg[reg][anc] += qtl_e - res_s

print 'Total distributions of ancestors in QTL regions:'
tot = 0.
for anc, ct in QTL_dist.items():
    tot += ct
print '  total area: ' + str(tot)
for anc, ct in QTL_dist.items():
    print '  %s: %f' % (anc, ct/tot)

mice_colors = {
    'Unk'       : '#707070', #gray
    'DBA2'      : 'r', #red
    'A'         : 'g', #green
    'AKR'       : '#3366FF', #light blue
    'BALBc'     : '#FFFF00', #yellow
    'C3HHe'     : '#FF9900', #orange
    'C57BL6N'   : '#800080', #purple
}

pylab.figure(figsize=(6,6), facecolor='white')
pylab.pie(QTL_dist.values(), explode=None, labels=QTL_dist.keys(), autopct='%1.1f%%',  
            shadow=True, colors=map(lambda x: mice_colors[x], QTL_dist.keys()))
pylab.title('Ancestor distribution per ancestor')
pylab.savefig('QTL_distribution_per_ancestor.png')


pylab.figure(facecolor='white')
i = 1
for reg in QTL_per_reg:
    pylab.subplot(3,5,i)
    i += 1
    pylab.pie(QTL_per_reg[reg].values(), explode=None, #labels=QTL_per_reg[reg].keys(), autopct='%1.1f%%',  
                shadow=True, colors=map(lambda x: mice_colors[x], QTL_per_reg[reg].keys()))
    pylab.title(reg)
#pylab.title('Ancestor distribution per ancestor')
pylab.savefig('QTL_distribution_per_region.png')

