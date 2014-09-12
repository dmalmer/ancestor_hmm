
from scipy.interpolate import spline
from matplotlib import pyplot
from numpy import linspace, mean, std

c = []
x = []
y = []

with open('../data/mouse_recomb_rates.csv','r') as f:
    f.readline()
    
    line = f.readline().strip().split(',')
    curr_chr = int(line[0][3:])
    adjustment = -1 * int(float(line[1]) * 1000)
    while line != ['']:
        if int(line[0][3:]) != curr_chr:
            c[-1] = -.015
            curr_chr = int(line[0][3:])
            adjustment = x[-1] - int(float(line[1]) * 1000) + 10
        #c.append((curr_chr + 1) * .1)
        c.append(-.001)
        x.append(int(float(line[1]) * 1000) + adjustment)
        y.append(float(line[2]))
        line = f.readline().strip().split(',')
    c[-1] = -.015
    c[0] = -.015

ws = 100
y_smooth = [mean(y[max(0, i-ws):min(len(y)-1, i+ws)]) for i in range(len(y))]

pyplot.figure(facecolor='white')

pyplot.plot(x, y_smooth, color='#004444')
pyplot.plot(x, c, color='k')

pyplot.xlim(0, max(x))
pyplot.ylim(-.05,.15)

pyplot.show()
