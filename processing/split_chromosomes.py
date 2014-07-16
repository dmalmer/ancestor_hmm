
import sys

filename = sys.argv[1].rsplit('.',1)[0]
ext = sys.argv[1].rsplit('.',1)[1]

f_in = open(filename + '.' + ext,'r')

line = f_in.readline()
cur_chr = line.split('\t')[0]
f_out = open(filename + '_' + cur_chr + '.' + ext,'w')
while line != '':
    if line.split('\t')[0] == cur_chr:
        f_out.write(line)
    else:
        f_out.close()
        cur_chr = line.split('\t')[0]
        f_out = open(filename + '_' + cur_chr + '.' + ext,'w')
        f_out.write(line)
    line = f_in.readline()

