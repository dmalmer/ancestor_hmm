
#PBS -S /bin/sh
#PBS -V

#PBS -q long
#PBS -l nodes=1:ppn=1
#PBS -l mem=1gb

#PBS -m ae
#PBS -M daniel.malmer@colorado.edu

python2.7 /Users/dama9282/AncestorInference/hmm/hmm_main.py -i /Users/dama9282/AncestorInference/data/${STRAIN}_temporary_sorted${CHR}.bed -r -p
