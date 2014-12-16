
cd output

iss_c()
{
	qsub -v STRAIN=ISS -N hmm_ISS_comp ../run_hmm.sh
}

ils_c()
{
	qsub -v STRAIN=ILS -N hmm_ILS_comp ../run_hmm.sh
}

iss_a()
{
	qsub -v STRAIN=ISS,P="-a 1.5 -u .6 -o _a15u6" -N hmm_ISS_a15u6 ../run_hmm.sh
	qsub -v STRAIN=ISS,P="-a 1.5 -u .8 -o _a15u8" -N hmm_ISS_a15u8 ../run_hmm.sh
	qsub -v STRAIN=ISS,P="-a 1. -u .6 -o _a10u6" -N hmm_ISS_a10u6 ../run_hmm.sh
	qsub -v STRAIN=ISS,P="-a 1. -u .8 -o _a10u8" -N hmm_ISS_a10u8 ../run_hmm.sh
	qsub -v STRAIN=ISS,P="-a .5 -u .6 -o _a5u6" -N hmm_ISS_a5u6 ../run_hmm.sh
	qsub -v STRAIN=ISS,P="-a .5 -u .8 -o _a5u8" -N hmm_ISS_a5u8 ../run_hmm.sh
}

ils_a()
{
	qsub -v STRAIN=ILS,P="-a 1.5 -u .6 -o _a15u6" -N hmm_ILS_a15u6 ../run_hmm.sh
	qsub -v STRAIN=ILS,P="-a 1.5 -u .8 -o _a15u8" -N hmm_ILS_a15u8 ../run_hmm.sh
	qsub -v STRAIN=ILS,P="-a 1. -u .6 -o _a10u6" -N hmm_ILS_a10u6 ../run_hmm.sh
	qsub -v STRAIN=ILS,P="-a 1. -u .8 -o _a10u8" -N hmm_ILS_a10u8 ../run_hmm.sh
	qsub -v STRAIN=ILS,P="-a .5 -u .6 -o _a5u6" -N hmm_ILS_a5u6 ../run_hmm.sh
	qsub -v STRAIN=ILS,P="-a .5 -u .8 -o _a5u8" -N hmm_ILS_a5u8 ../run_hmm.sh
}

iss_p()
{
	qsub -v STRAIN=ISS,P="-pw -a 1.5 -u .6 -o _a15u6pw" -N hmm_ISS_a15u6pw ../run_hmm.sh
	qsub -v STRAIN=ISS,P="-pw -a 1.5 -u .8 -o _a15u8pw" -N hmm_ISS_a15u8pw ../run_hmm.sh
	qsub -v STRAIN=ISS,P="-pw -a 1. -u .6 -o _a10u6pw" -N hmm_ISS_a10u6pw ../run_hmm.sh
	qsub -v STRAIN=ISS,P="-pw -a 1. -u .8 -o _a10u8pw" -N hmm_ISS_a10u8pw ../run_hmm.sh
	qsub -v STRAIN=ISS,P="-pw -a .5 -u .6 -o _a5u6pw" -N hmm_ISS_a5u6pw ../run_hmm.sh
	qsub -v STRAIN=ISS,P="-pw -a .5 -u .8 -o _a5u8pw" -N hmm_ISS_a5u8pw ../run_hmm.sh
}

ils_p()
{
	qsub -v STRAIN=ILS,P="-pw -a 1.5 -u .6 -o _a15u6pw" -N hmm_ILS_a15u6pw ../run_hmm.sh
	qsub -v STRAIN=ILS,P="-pw -a 1.5 -u .8 -o _a15u8pw" -N hmm_ILS_a15u8pw ../run_hmm.sh
	qsub -v STRAIN=ILS,P="-pw -a 1. -u .6 -o _a10u6pw" -N hmm_ILS_a10u6pw ../run_hmm.sh
	qsub -v STRAIN=ILS,P="-pw -a 1. -u .8 -o _a10u8pw" -N hmm_ILS_a10u8pw ../run_hmm.sh
	qsub -v STRAIN=ILS,P="-pw -a .5 -u .6 -o _a5u6pw" -N hmm_ILS_a5u6pw ../run_hmm.sh
	qsub -v STRAIN=ILS,P="-pw -a .5 -u .8 -o _a5u8pw" -N hmm_ILS_a5u8pw ../run_hmm.sh
}

iss_u()
{
	qsub -v STRAIN=ISS,P="-a 1. -u .6 -o _a10u6" -N hmm_ISS_a10u6 ../run_hmm.sh
	qsub -v STRAIN=ISS,P="-a 1. -u .7 -o _a10u7" -N hmm_ISS_a10u7 ../run_hmm.sh
	qsub -v STRAIN=ISS,P="-a 1. -u .8 -o _a10u8" -N hmm_ISS_a10u8 ../run_hmm.sh
	qsub -v STRAIN=ISS,P="-a 1. -u .9 -o _a10u9" -N hmm_ISS_a10u9 ../run_hmm.sh
}

ils_u()
{
	qsub -v STRAIN=ILS,P="-a 1. -u .6 -o _a10u6" -N hmm_ILS_a10u6 ../run_hmm.sh
	qsub -v STRAIN=ILS,P="-a 1. -u .7 -o _a10u7" -N hmm_ILS_a10u7 ../run_hmm.sh
	qsub -v STRAIN=ILS,P="-a 1. -u .8 -o _a10u8" -N hmm_ILS_a10u8 ../run_hmm.sh
	qsub -v STRAIN=ILS,P="-a 1. -u .9 -o _a10u9" -N hmm_ILS_a10u9 ../run_hmm.sh
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
elif [[ $1 == "ALL_A" ]]
then
	iss_a
	ils_a
elif [[ $1 == "ISS_A" ]]
then
	iss_a
elif [[ $1 == "ILS_A" ]]
then
	ils_a
elif [[ $1 == "ALL_U" ]]
then
	iss_u
	ils_u
elif [[ $1 == "ISS_U" ]]
then
	iss_u
elif [[ $1 == "ILS_U" ]]
then
	ils_u
elif [[ $1 == "ALL_P" ]]
then
	iss_p
	ils_p
elif [[ $1 == "ISS_P" ]]
then
	iss_p
elif [[ $1 == "ILS_P" ]]
then
	ils_p
elif [[ $1 == "TEST" ]]
then
    test_run
fi
