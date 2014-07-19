
from numpy import mean, std
from numpy import genfromtxt

# find std_dev and mean
rates = genfromtxt('../data/mouse_recomb_rates.csv', delimiter=',', skiprows=1, usecols=2)
rates_std_dev = std(rates)
rates_mean = mean(rates)
#print rates_std_dev
#print rates_mean

# find normalizing number after removing > and < 2 std_dev's from mean
rates = [r for r in rates if r < rates_mean + 2*rates_std_dev and r > rates_mean - 2*rates_std_dev]
rates_min = min(rates)
rates_max = max(rates)
#print rates_min
#print rates_max

normal_num = 1000/(rates_max - rates_min)

# create wig file with heatmap of recomb rates
#  --anything > 2 std_dev's from mean, change to max (1000)
outside_ct = 0
with open('../data/mouse_recomb_rates.wig', 'w') as f_out:
	f_out.write('track type=wiggle_0 graphType=heatmap\n')

	with open('../data/mouse_recomb_rates.csv', 'r') as f_in:
		next(f_in)
		chrom = ''

		for line in f_in:
			splits = line.split(',')

			new_chrom = splits[0]
			if new_chrom != chrom:
				chrom = new_chrom
				f_out.write('variableStep chrom=%s\n' % chrom)

			if float(splits[2]) > rates_mean + 2*rates_std_dev:
				f_out.write('%i\t%i\n' % (int(float(splits[1])*1000), 1000))
				outside_ct += 1
			else:
				f_out.write('%i\t%i\n' % (int(float(splits[1])*1000), int(float(splits[2]) * normal_num)))

print outside_ct