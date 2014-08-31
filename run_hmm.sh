
#PBS -S /bin/sh
#PBS -V

#PBS -q long
#PBS -l nodes=1:ppn=16
#PBS -l mem=8gb
#PBS -l walltime=200:00:00

#PBS -m bae
#PBS -M daniel.malmer@colorado.edu

python2.7 /Users/dama9282/AncestorInference/hmm/hmm_main.py -i /Users/dama9282/AncestorInference/data/${STRAIN}_full_sorted${CHR}.bed -r -p -v ${A} ${U} ${O}
