

iss()
{
	qsub -v STRAIN=ISS,CHR=1 -N hmm_ISS_1 run_hmm.sh
	qsub -v STRAIN=ISS,CHR=2 -N hmm_ISS_2 run_hmm.sh
	qsub -v STRAIN=ISS,CHR=3 -N hmm_ISS_3 run_hmm.sh
	qsub -v STRAIN=ISS,CHR=4 -N hmm_ISS_4 run_hmm.sh
	qsub -v STRAIN=ISS,CHR=5 -N hmm_ISS_5 run_hmm.sh
	qsub -v STRAIN=ISS,CHR=6 -N hmm_ISS_6 run_hmm.sh
	qsub -v STRAIN=ISS,CHR=7 -N hmm_ISS_7 run_hmm.sh
	qsub -v STRAIN=ISS,CHR=8 -N hmm_ISS_8 run_hmm.sh
	qsub -v STRAIN=ISS,CHR=9 -N hmm_ISS_9 run_hmm.sh
	qsub -v STRAIN=ISS,CHR=10 -N hmm_ISS_10 run_hmm.sh
	qsub -v STRAIN=ISS,CHR=11 -N hmm_ISS_11 run_hmm.sh
	qsub -v STRAIN=ISS,CHR=12 -N hmm_ISS_12 run_hmm.sh
	qsub -v STRAIN=ISS,CHR=13 -N hmm_ISS_13 run_hmm.sh
	qsub -v STRAIN=ISS,CHR=14 -N hmm_ISS_14 run_hmm.sh
	qsub -v STRAIN=ISS,CHR=15 -N hmm_ISS_15 run_hmm.sh
	qsub -v STRAIN=ISS,CHR=16 -N hmm_ISS_16 run_hmm.sh
	qsub -v STRAIN=ISS,CHR=17 -N hmm_ISS_17 run_hmm.sh
	qsub -v STRAIN=ISS,CHR=18 -N hmm_ISS_18 run_hmm.sh
	qsub -v STRAIN=ISS,CHR=19 -N hmm_ISS_19 run_hmm.sh
	qsub -v STRAIN=ISS,CHR=X -N hmm_ISS_X run_hmm.sh
	qsub -v STRAIN=ISS,CHR=Y -N hmm_ISS_Y run_hmm.sh
}

ils()
{
	qsub -v STRAIN=ILS,CHR=1 -N hmm_ILS_1 run_hmm.sh
	qsub -v STRAIN=ILS,CHR=2 -N hmm_ILS_2 run_hmm.sh
	qsub -v STRAIN=ILS,CHR=3 -N hmm_ILS_3 run_hmm.sh
	qsub -v STRAIN=ILS,CHR=4 -N hmm_ILS_4 run_hmm.sh
	qsub -v STRAIN=ILS,CHR=5 -N hmm_ILS_5 run_hmm.sh
	qsub -v STRAIN=ILS,CHR=6 -N hmm_ILS_6 run_hmm.sh
	qsub -v STRAIN=ILS,CHR=7 -N hmm_ILS_7 run_hmm.sh
	qsub -v STRAIN=ILS,CHR=8 -N hmm_ILS_8 run_hmm.sh
	qsub -v STRAIN=ILS,CHR=9 -N hmm_ILS_9 run_hmm.sh
	qsub -v STRAIN=ILS,CHR=10 -N hmm_ILS_10 run_hmm.sh
	qsub -v STRAIN=ILS,CHR=11 -N hmm_ILS_11 run_hmm.sh
	qsub -v STRAIN=ILS,CHR=12 -N hmm_ILS_12 run_hmm.sh
	qsub -v STRAIN=ILS,CHR=13 -N hmm_ILS_13 run_hmm.sh
	qsub -v STRAIN=ILS,CHR=14 -N hmm_ILS_14 run_hmm.sh
	qsub -v STRAIN=ILS,CHR=15 -N hmm_ILS_15 run_hmm.sh
	qsub -v STRAIN=ILS,CHR=16 -N hmm_ILS_16 run_hmm.sh
	qsub -v STRAIN=ILS,CHR=17 -N hmm_ILS_17 run_hmm.sh
	qsub -v STRAIN=ILS,CHR=18 -N hmm_ILS_18 run_hmm.sh
	qsub -v STRAIN=ILS,CHR=19 -N hmm_ILS_19 run_hmm.sh
	qsub -v STRAIN=ILS,CHR=X -N hmm_ILS_X run_hmm.sh
	qsub -v STRAIN=ILS,CHR=Y -N hmm_ILS_Y run_hmm.sh
}

if [[ $1 == "ALL" ]]
then
	iss
	ils
elif [[ $1 == "ISS" ]]
then
	iss
elif [[ $1 == "ILS" ]]
then
	ils
fi
