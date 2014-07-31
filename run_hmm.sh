
#PBS -S /bin/sh
#PBS -V

#PBS -q long
#PBS -l nodes=1:ppn=1
#PBS -l mem=256mb

#PBS -m ae
#PBS -M daniel.malmer@colorado.edu

python2.7 /Users/dama9282/AncestorInference/hmm/hmm_main.py /Users/dama9282/AncestorInference/data/${STRAIN}_temporary_sorted${CHR}.bed
