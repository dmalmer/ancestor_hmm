
#PBS -S /bin/sh
#PBS -V

#PBS -q long
#PBS -l nodes=1:ppn=1
#PBS -l mem=256mb

#PBS -m ae
#PBS -M daniel.malmer@colorado.edu

### #PBS -o /Users/dama9282/AncestorInference/output/test.o
#### PBS -j oe

python2.7 /Users/dama9282/AncestorInference/hmm_main.py /Users/dama9282/AncestorInference/data/${STRAIN}_temporary_sorted_chr${CHR}.txt
