
from scipy.interpolate import spline
from matplotlib import pyplot
from numpy import linspace, mean, std

c = []
x = []
y = []

with open('../data/mouse_recomb_rates.csv','r') as f:
    f.readline()
    
    line = f.readline().strip().split(',')
    curr_chr = line[0]
    adjustment = -1 * int(float(line[1]) * 1000)
    while line != ['']:
        if line[0] != curr_chr:
            curr_chr = line[0]
            adjustment = x[-1] - int(float(line[1]) * 1000) + 10
        c.append(curr_chr)
        x.append(int(float(line[1]) * 1000) + adjustment)
        y.append(float(line[2]))
        line = f.readline().strip().split(',')

ws = 100
y_smooth = [mean(y[max(0, i-ws):min(len(y)-1, i+ws)]) for i in range(len(y))]

pyplot.plot(x, y_smooth)
pyplot.show()
