
cd output

iss_c()
{
	qsub -v STRAIN=ISS,CHR= -N hmm_ISS_comp ../run_hmm.sh
}

ils_c()
{
	qsub -v STRAIN=ILS,CHR= -N hmm_ILS_comp ../run_hmm.sh
}

iss_i()
{
	qsub -v STRAIN=ISS,CHR=_chr1 -N hmm_ISS_1 ../run_hmm.sh
	qsub -v STRAIN=ISS,CHR=_chr2 -N hmm_ISS_2 ../run_hmm.sh
	qsub -v STRAIN=ISS,CHR=_chr3 -N hmm_ISS_3 ../run_hmm.sh
	qsub -v STRAIN=ISS,CHR=_chr4 -N hmm_ISS_4 ../run_hmm.sh
	qsub -v STRAIN=ISS,CHR=_chr5 -N hmm_ISS_5 ../run_hmm.sh
	qsub -v STRAIN=ISS,CHR=_chr6 -N hmm_ISS_6 ../run_hmm.sh
	qsub -v STRAIN=ISS,CHR=_chr7 -N hmm_ISS_7 ../run_hmm.sh
	qsub -v STRAIN=ISS,CHR=_chr8 -N hmm_ISS_8 ../run_hmm.sh
	qsub -v STRAIN=ISS,CHR=_chr9 -N hmm_ISS_9 ../run_hmm.sh
	qsub -v STRAIN=ISS,CHR=_chr10 -N hmm_ISS_10 ../run_hmm.sh
	qsub -v STRAIN=ISS,CHR=_chr11 -N hmm_ISS_11 ../run_hmm.sh
	qsub -v STRAIN=ISS,CHR=_chr12 -N hmm_ISS_12 ../run_hmm.sh
	qsub -v STRAIN=ISS,CHR=_chr13 -N hmm_ISS_13 ../run_hmm.sh
	qsub -v STRAIN=ISS,CHR=_chr14 -N hmm_ISS_14 ../run_hmm.sh
	qsub -v STRAIN=ISS,CHR=_chr15 -N hmm_ISS_15 ../run_hmm.sh
	qsub -v STRAIN=ISS,CHR=_chr16 -N hmm_ISS_16 ../run_hmm.sh
	qsub -v STRAIN=ISS,CHR=_chr17 -N hmm_ISS_17 ../run_hmm.sh
	qsub -v STRAIN=ISS,CHR=_chr18 -N hmm_ISS_18 ../run_hmm.sh
	qsub -v STRAIN=ISS,CHR=_chr19 -N hmm_ISS_19 ../run_hmm.sh
	qsub -v STRAIN=ISS,CHR=_chrX -N hmm_ISS_X ../run_hmm.sh
	qsub -v STRAIN=ISS,CHR=_chrY -N hmm_ISS_Y ../run_hmm.sh
}

ils_i()
{
	qsub -v STRAIN=ILS,CHR=_chr1 -N hmm_ILS_1 ../run_hmm.sh
	qsub -v STRAIN=ILS,CHR=_chr2 -N hmm_ILS_2 ../run_hmm.sh
	qsub -v STRAIN=ILS,CHR=_chr3 -N hmm_ILS_3 ../run_hmm.sh
	qsub -v STRAIN=ILS,CHR=_chr4 -N hmm_ILS_4 ../run_hmm.sh
	qsub -v STRAIN=ILS,CHR=_chr5 -N hmm_ILS_5 ../run_hmm.sh
	qsub -v STRAIN=ILS,CHR=_chr6 -N hmm_ILS_6 ../run_hmm.sh
	qsub -v STRAIN=ILS,CHR=_chr7 -N hmm_ILS_7 ../run_hmm.sh
	qsub -v STRAIN=ILS,CHR=_chr8 -N hmm_ILS_8 ../run_hmm.sh
	qsub -v STRAIN=ILS,CHR=_chr9 -N hmm_ILS_9 ../run_hmm.sh
	qsub -v STRAIN=ILS,CHR=_chr10 -N hmm_ILS_10 ../run_hmm.sh
	qsub -v STRAIN=ILS,CHR=_chr11 -N hmm_ILS_11 ../run_hmm.sh
	qsub -v STRAIN=ILS,CHR=_chr12 -N hmm_ILS_12 ../run_hmm.sh
	qsub -v STRAIN=ILS,CHR=_chr13 -N hmm_ILS_13 ../run_hmm.sh
	qsub -v STRAIN=ILS,CHR=_chr14 -N hmm_ILS_14 ../run_hmm.sh
	qsub -v STRAIN=ILS,CHR=_chr15 -N hmm_ILS_15 ../run_hmm.sh
	qsub -v STRAIN=ILS,CHR=_chr16 -N hmm_ILS_16 ../run_hmm.sh
	qsub -v STRAIN=ILS,CHR=_chr17 -N hmm_ILS_17 ../run_hmm.sh
	qsub -v STRAIN=ILS,CHR=_chr18 -N hmm_ILS_18 ../run_hmm.sh
	qsub -v STRAIN=ILS,CHR=_chr19 -N hmm_ILS_19 ../run_hmm.sh
	qsub -v STRAIN=ILS,CHR=_chrX -N hmm_ILS_X ../run_hmm.sh
	qsub -v STRAIN=ILS,CHR=_chrY -N hmm_ILS_Y ../run_hmm.sh
}

test_run()
{
	qsub -v STRAIN=TEST,CHR= -N hmm_TEST ../run_hmm.sh
}

if [[ $1 == "ALL_C" ]]
then
	iss_c
	ils_c
elif [[ $1 == "ISS_C" ]]
then
	iss_c
elif [[ $1 == "ILS_C" ]]
then
	ils_c
elif [[ $1 == "ALL_I" ]]
then
	iss_i
	ils_i
elif [[ $1 == "ISS_I" ]]
then
	iss_i
elif [[ $1 == "ILS_I" ]]
then
	ils_i
elif [[ $1 == "TEST" ]]
then
    test_run
fi
