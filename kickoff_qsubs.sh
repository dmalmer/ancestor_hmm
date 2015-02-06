
cd output

iss_c()
{
	qsub -v STRAIN=ISS -N hmm_ISS_comp ../run_hmm.sh
}

ils_c()
{
	qsub -v STRAIN=ILS -N hmm_ILS_comp ../run_hmm.sh
}

iss_g()
{
	qsub -v STRAIN=ISS,P="-t [.4-.99] -e [.98-.9999] -a [.2-2] -u [.2-.95] -g 10 -m 1" -N hmm_ISS_grid ../run_hmm.sh
}

ils_g()
{
	qsub -v STRAIN=ILS,P="-t [.4-.99] -e [.98-.9999] -a [.2-2] -u [.2-.95] -g 10 -m 1" -N hmm_ILS_grid ../run_hmm.sh
}

iss_a()
{
	qsub -v STRAIN=ISS,P="-a 1.5 -u .6 -o _a15u6 -m 20" -N hmm_ISS_a15u6 ../run_hmm.sh
	qsub -v STRAIN=ISS,P="-a 1.5 -u .8 -o _a15u8 -m 20" -N hmm_ISS_a15u8 ../run_hmm.sh
	qsub -v STRAIN=ISS,P="-a 1. -u .6 -o _a10u6 -m 20" -N hmm_ISS_a10u6 ../run_hmm.sh
	qsub -v STRAIN=ISS,P="-a 1. -u .8 -o _a10u8 -m 20" -N hmm_ISS_a10u8 ../run_hmm.sh
	qsub -v STRAIN=ISS,P="-a .5 -u .6 -o _a5u6 -m 20" -N hmm_ISS_a5u6 ../run_hmm.sh
	qsub -v STRAIN=ISS,P="-a .5 -u .8 -o _a5u8 -m 20" -N hmm_ISS_a5u8 ../run_hmm.sh
}

ils_a()
{
	qsub -v STRAIN=ILS,P="-a 1.5 -u .6 -o _a15u6 -m 20" -N hmm_ILS_a15u6 ../run_hmm.sh
	qsub -v STRAIN=ILS,P="-a 1.5 -u .8 -o _a15u8 -m 20" -N hmm_ILS_a15u8 ../run_hmm.sh
	qsub -v STRAIN=ILS,P="-a 1. -u .6 -o _a10u6 -m 20" -N hmm_ILS_a10u6 ../run_hmm.sh
	qsub -v STRAIN=ILS,P="-a 1. -u .8 -o _a10u8 -m 20" -N hmm_ILS_a10u8 ../run_hmm.sh
	qsub -v STRAIN=ILS,P="-a .5 -u .6 -o _a5u6 -m 20" -N hmm_ILS_a5u6 ../run_hmm.sh
	qsub -v STRAIN=ILS,P="-a .5 -u .8 -o _a5u8 -m 20" -N hmm_ILS_a5u8 ../run_hmm.sh
}

iss_e()
{
	qsub -v STRAIN=ISS,P="-pw -e .96 -a 1.5 -o _e96a15pw -m 20" -N hmm_ISS_e96a15pw ../run_hmm.sh
	qsub -v STRAIN=ISS,P="-pw -e .96 -a 1.5 -o _e96a15pw -m 20" -N hmm_ISS_e96a15pw ../run_hmm.sh
	qsub -v STRAIN=ISS,P="-pw -e .96 -a 1. -o _e96a10pw -m 20" -N hmm_ISS_e96a10pw ../run_hmm.sh
	qsub -v STRAIN=ISS,P="-pw -e .96 -a 1. -o _e96a10pw -m 20" -N hmm_ISS_e96a10pw ../run_hmm.sh
	qsub -v STRAIN=ISS,P="-pw -e .96 -a .5 -o _e96a5pw -m 20" -N hmm_ISS_e96a5pw ../run_hmm.sh
	qsub -v STRAIN=ISS,P="-pw -e .96 -a .5 -o _e96a5pw -m 20" -N hmm_ISS_e96a5pw ../run_hmm.sh
	
    qsub -v STRAIN=ISS,P="-pw -e .94 -a 1.5 -o _e94a15pw -m 20" -N hmm_ISS_e94a15pw ../run_hmm.sh
	qsub -v STRAIN=ISS,P="-pw -e .94 -a 1.5 -o _e94a15pw -m 20" -N hmm_ISS_e94a15pw ../run_hmm.sh
	qsub -v STRAIN=ISS,P="-pw -e .94 -a 1. -o _e94a10pw -m 20" -N hmm_ISS_e94a10pw ../run_hmm.sh
	qsub -v STRAIN=ISS,P="-pw -e .94 -a 1. -o _e94a10pw -m 20" -N hmm_ISS_e94a10pw ../run_hmm.sh
	qsub -v STRAIN=ISS,P="-pw -e .94 -a .5 -o _e94a5pw -m 20" -N hmm_ISS_e94a5pw ../run_hmm.sh
	qsub -v STRAIN=ISS,P="-pw -e .94 -a .5 -o _e94a5pw -m 20" -N hmm_ISS_e94a5pw ../run_hmm.sh
}

ils_e()
{
	qsub -v STRAIN=ILS,P="-pw -e .96 -a 1.5 -o _e96a15pw -m 20" -N hmm_ILS_e96a15pw ../run_hmm.sh
	qsub -v STRAIN=ILS,P="-pw -e .96 -a 1.5 -o _e96a15pw -m 20" -N hmm_ILS_e96a15pw ../run_hmm.sh
	qsub -v STRAIN=ILS,P="-pw -e .96 -a 1. -o _e96a10pw -m 20" -N hmm_ILS_e96a10pw ../run_hmm.sh
	qsub -v STRAIN=ILS,P="-pw -e .96 -a 1. -o _e96a10pw -m 20" -N hmm_ILS_e96a10pw ../run_hmm.sh
	qsub -v STRAIN=ILS,P="-pw -e .96 -a .5 -o _e96a5pw -m 20" -N hmm_ILS_e96a5pw ../run_hmm.sh
	qsub -v STRAIN=ILS,P="-pw -e .96 -a .5 -o _e96a5pw -m 20" -N hmm_ILS_e96a5pw ../run_hmm.sh

	qsub -v STRAIN=ISS,P="-pw -e .94 -a 1.5 -o _e94a15pw -m 20" -N hmm_ISS_e94a15pw ../run_hmm.sh
	qsub -v STRAIN=ISS,P="-pw -e .94 -a 1.5 -o _e94a15pw -m 20" -N hmm_ISS_e94a15pw ../run_hmm.sh
	qsub -v STRAIN=ISS,P="-pw -e .94 -a 1. -o _e94a10pw -m 20" -N hmm_ISS_e94a10pw ../run_hmm.sh
	qsub -v STRAIN=ISS,P="-pw -e .94 -a 1. -o _e94a10pw -m 20" -N hmm_ISS_e94a10pw ../run_hmm.sh
	qsub -v STRAIN=ISS,P="-pw -e .94 -a .5 -o _e94a5pw -m 20" -N hmm_ISS_e94a5pw ../run_hmm.sh
	qsub -v STRAIN=ISS,P="-pw -e .94 -a .5 -o _e94a5pw -m 20" -N hmm_ISS_e94a5pw ../run_hmm.sh
}

iss_t()
{
	qsub -v STRAIN=ISS,P="-pw -t .82 -u .9 -a 1.5 -o _t82u9a15pw -m 20" -N hmm_ISS_t82u9a15pw ../run_hmm.sh
	qsub -v STRAIN=ISS,P="-pw -t .82 -u .9 -a 1.5 -o _t82u9a15pw -m 20" -N hmm_ISS_t82u9a15pw ../run_hmm.sh
	qsub -v STRAIN=ISS,P="-pw -t .82 -u .9 -a 1. -o _t82u9a10pw -m 20" -N hmm_ISS_t82u9a10pw ../run_hmm.sh
	qsub -v STRAIN=ISS,P="-pw -t .82 -u .9 -a 1. -o _t82u9a10pw -m 20" -N hmm_ISS_t82u9a10pw ../run_hmm.sh
	qsub -v STRAIN=ISS,P="-pw -t .82 -u .9 -a .5 -o _t82u9a5pw -m 20" -N hmm_ISS_t82u9a5pw ../run_hmm.sh
	qsub -v STRAIN=ISS,P="-pw -t .82 -u .9 -a .5 -o _t82u9a5pw -m 20" -N hmm_ISS_t82u9a5pw ../run_hmm.sh
	
	qsub -v STRAIN=ISS,P="-pw -t .64 -u .9 -a 1.5 -o _t64u9a15pw -m 20" -N hmm_ISS_t64u9a15pw ../run_hmm.sh
	qsub -v STRAIN=ISS,P="-pw -t .64 -u .9 -a 1.5 -o _t64u9a15pw -m 20" -N hmm_ISS_t64u9a15pw ../run_hmm.sh
	qsub -v STRAIN=ISS,P="-pw -t .64 -u .9 -a 1. -o _t64u9a10pw -m 20" -N hmm_ISS_t64u9a10pw ../run_hmm.sh
	qsub -v STRAIN=ISS,P="-pw -t .64 -u .9 -a 1. -o _t64u9a10pw -m 20" -N hmm_ISS_t64u9a10pw ../run_hmm.sh
	qsub -v STRAIN=ISS,P="-pw -t .64 -u .9 -a .5 -o _t64u9a5pw -m 20" -N hmm_ISS_t64u9a5pw ../run_hmm.sh
	qsub -v STRAIN=ISS,P="-pw -t .64 -u .9 -a .5 -o _t64u9a5pw -m 20" -N hmm_ISS_t64u9a5pw ../run_hmm.sh
	
    qsub -v STRAIN=ISS,P="-pw -t .52 -u .9 -a 1.5 -o _t52u9a15pw -m 20" -N hmm_ISS_t52u9a15pw ../run_hmm.sh
	qsub -v STRAIN=ISS,P="-pw -t .52 -u .9 -a 1.5 -o _t52u9a15pw -m 20" -N hmm_ISS_t52u9a15pw ../run_hmm.sh
	qsub -v STRAIN=ISS,P="-pw -t .52 -u .9 -a 1. -o _t52u9a10pw -m 20" -N hmm_ISS_t52u9a10pw ../run_hmm.sh
	qsub -v STRAIN=ISS,P="-pw -t .52 -u .9 -a 1. -o _t52u9a10pw -m 20" -N hmm_ISS_t52u9a10pw ../run_hmm.sh
	qsub -v STRAIN=ISS,P="-pw -t .52 -u .9 -a .5 -o _t52u9a5pw -m 20" -N hmm_ISS_t52u9a5pw ../run_hmm.sh
	qsub -v STRAIN=ISS,P="-pw -t .52 -u .9 -a .5 -o _t52u9a5pw -m 20" -N hmm_ISS_t52u9a5pw ../run_hmm.sh
}

ils_t()
{
	qsub -v STRAIN=ILS,P="-pw -t .82 -u .9 -a 1.5 -o _t82u9a15pw -m 20" -N hmm_ILS_t82u9a15pw ../run_hmm.sh
	qsub -v STRAIN=ILS,P="-pw -t .82 -u .9 -a 1.5 -o _t82u9a15pw -m 20" -N hmm_ILS_t82u9a15pw ../run_hmm.sh
	qsub -v STRAIN=ILS,P="-pw -t .82 -u .9 -a 1. -o _t82u9a10pw -m 20" -N hmm_ILS_t82u9a10pw ../run_hmm.sh
	qsub -v STRAIN=ILS,P="-pw -t .82 -u .9 -a 1. -o _t82u9a10pw -m 20" -N hmm_ILS_t82u9a10pw ../run_hmm.sh
	qsub -v STRAIN=ILS,P="-pw -t .82 -u .9 -a .5 -o _t82u9a5pw -m 20" -N hmm_ILS_t82u9a5pw ../run_hmm.sh
	qsub -v STRAIN=ILS,P="-pw -t .82 -u .9 -a .5 -o _t82u9a5pw -m 20" -N hmm_ILS_t82u9a5pw ../run_hmm.sh

	qsub -v STRAIN=ILS,P="-pw -t .64 -u .9 -a 1.5 -o _t64u9a15pw -m 20" -N hmm_ILS_t64u9a15pw ../run_hmm.sh
	qsub -v STRAIN=ILS,P="-pw -t .64 -u .9 -a 1.5 -o _t64u9a15pw -m 20" -N hmm_ILS_t64u9a15pw ../run_hmm.sh
	qsub -v STRAIN=ILS,P="-pw -t .64 -u .9 -a 1. -o _t64u9a10pw -m 20" -N hmm_ILS_t64u9a10pw ../run_hmm.sh
	qsub -v STRAIN=ILS,P="-pw -t .64 -u .9 -a 1. -o _t64u9a10pw -m 20" -N hmm_ILS_t64u9a10pw ../run_hmm.sh
	qsub -v STRAIN=ILS,P="-pw -t .64 -u .9 -a .5 -o _t64u9a5pw -m 20" -N hmm_ILS_t64u9a5pw ../run_hmm.sh
	qsub -v STRAIN=ILS,P="-pw -t .64 -u .9 -a .5 -o _t64u9a5pw -m 20" -N hmm_ILS_t64u9a5pw ../run_hmm.sh

	qsub -v STRAIN=ISS,P="-pw -t .52 -u .9 -a 1.5 -o _t52u9a15pw -m 20" -N hmm_ISS_t52u9a15pw ../run_hmm.sh
	qsub -v STRAIN=ISS,P="-pw -t .52 -u .9 -a 1.5 -o _t52u9a15pw -m 20" -N hmm_ISS_t52u9a15pw ../run_hmm.sh
	qsub -v STRAIN=ISS,P="-pw -t .52 -u .9 -a 1. -o _t52u9a10pw -m 20" -N hmm_ISS_t52u9a10pw ../run_hmm.sh
	qsub -v STRAIN=ISS,P="-pw -t .52 -u .9 -a 1. -o _t52u9a10pw -m 20" -N hmm_ISS_t52u9a10pw ../run_hmm.sh
	qsub -v STRAIN=ISS,P="-pw -t .52 -u .9 -a .5 -o _t52u9a5pw -m 20" -N hmm_ISS_t52u9a5pw ../run_hmm.sh
	qsub -v STRAIN=ISS,P="-pw -t .52 -u .9 -a .5 -o _t52u9a5pw -m 20" -N hmm_ISS_t52u9a5pw ../run_hmm.sh
}

iss_p()
{
	qsub -v STRAIN=ISS,P="-pw -a 1.5 -u .6 -o _a15u6pw -m 20" -N hmm_ISS_a15u6pw ../run_hmm.sh
	qsub -v STRAIN=ISS,P="-pw -a 1.5 -u .8 -o _a15u8pw -m 20" -N hmm_ISS_a15u8pw ../run_hmm.sh
	qsub -v STRAIN=ISS,P="-pw -a 1. -u .6 -o _a10u6pw -m 20" -N hmm_ISS_a10u6pw ../run_hmm.sh
	qsub -v STRAIN=ISS,P="-pw -a 1. -u .8 -o _a10u8pw -m 20" -N hmm_ISS_a10u8pw ../run_hmm.sh
	qsub -v STRAIN=ISS,P="-pw -a .5 -u .6 -o _a5u6pw -m 20" -N hmm_ISS_a5u6pw ../run_hmm.sh
	qsub -v STRAIN=ISS,P="-pw -a .5 -u .8 -o _a5u8pw -m 20" -N hmm_ISS_a5u8pw ../run_hmm.sh
}

ils_p()
{
	qsub -v STRAIN=ILS,P="-pw -a 1.5 -u .6 -o _a15u6pw -m 20" -N hmm_ILS_a15u6pw ../run_hmm.sh
	qsub -v STRAIN=ILS,P="-pw -a 1.5 -u .8 -o _a15u8pw -m 20" -N hmm_ILS_a15u8pw ../run_hmm.sh
	qsub -v STRAIN=ILS,P="-pw -a 1. -u .6 -o _a10u6pw -m 20" -N hmm_ILS_a10u6pw ../run_hmm.sh
	qsub -v STRAIN=ILS,P="-pw -a 1. -u .8 -o _a10u8pw -m 20" -N hmm_ILS_a10u8pw ../run_hmm.sh
	qsub -v STRAIN=ILS,P="-pw -a .5 -u .6 -o _a5u6pw -m 20" -N hmm_ILS_a5u6pw ../run_hmm.sh
	qsub -v STRAIN=ILS,P="-pw -a .5 -u .8 -o _a5u8pw -m 20" -N hmm_ILS_a5u8pw ../run_hmm.sh
}

iss_u()
{
	qsub -v STRAIN=ISS,P="-a 1. -u .6 -o _a10u6 -m 20" -N hmm_ISS_a10u6 ../run_hmm.sh
	qsub -v STRAIN=ISS,P="-a 1. -u .7 -o _a10u7 -m 20" -N hmm_ISS_a10u7 ../run_hmm.sh
	qsub -v STRAIN=ISS,P="-a 1. -u .8 -o _a10u8 -m 20" -N hmm_ISS_a10u8 ../run_hmm.sh
	qsub -v STRAIN=ISS,P="-a 1. -u .9 -o _a10u9 -m 20" -N hmm_ISS_a10u9 ../run_hmm.sh
}

ils_u()
{
	qsub -v STRAIN=ILS,P="-a 1. -u .6 -o _a10u6 -m 20" -N hmm_ILS_a10u6 ../run_hmm.sh
	qsub -v STRAIN=ILS,P="-a 1. -u .7 -o _a10u7 -m 20" -N hmm_ILS_a10u7 ../run_hmm.sh
	qsub -v STRAIN=ILS,P="-a 1. -u .8 -o _a10u8 -m 20" -N hmm_ILS_a10u8 ../run_hmm.sh
	qsub -v STRAIN=ILS,P="-a 1. -u .9 -o _a10u9 -m 20" -N hmm_ILS_a10u9 ../run_hmm.sh
}

test_run()
{
	qsub -v STRAIN=TEST,P="-p" -N hmm_TEST ../run_hmm.sh
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
elif [[ $1 == "ALL_G" ]]
then
	iss_g
	ils_g
elif [[ $1 == "ISS_G" ]]
then
	iss_g
elif [[ $1 == "ILS_G" ]]
then
	ils_g
elif [[ $1 == "ALL_E" ]]
then
	iss_e
	ils_e
elif [[ $1 == "ISS_E" ]]
then
	iss_e
elif [[ $1 == "ILS_E" ]]
then
	ils_e
elif [[ $1 == "ALL_T" ]]
then
	iss_t
	ils_t
elif [[ $1 == "ISS_T" ]]
then
	iss_t
elif [[ $1 == "ILS_T" ]]
then
	ils_t
elif [[ $1 == "TEST" ]]
then
    test_run
fi
