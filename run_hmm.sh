
#PBS -S /bin/sh
#PBS -V

###PBS -q long8gb

#PBS -l nodes=1:ppn=20
#PBS -l pmem=8gb
#PBS -l walltime=120:00:00

#PBS -m ae
#PBS -M daniel.malmer@colorado.edu

python2.7 /Users/dama9282/AncestorInference/hmm/hmm_main.py \
    -i /Users/dama9282/AncestorInference/data/${STRAIN}_sh-nsh_filt-anc_full.bed \
    -r -p -w -v ${P}
