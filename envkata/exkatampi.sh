#!/bin/sh

MaxSubBand=7
MaxCPU=7

gTubeType=$1
gSubBand=$2


SubmitJob () {
    local TubeType=$1
    local SubBand=$2
    local LOG="exkatampi-${TubeType}-e${SubBand}${SubBand}.log"
    local JOBNAME="E${SubBand}${SubBand}${TubeType}"
#    bsub -m "hpc12 hpc11" -o $LOG -J $JOBNAME -n $MaxCPU "mpijob mpirun ./exkatampi.out $TubeType $SubBand"
    bsub -o $LOG -J $JOBNAME -n $MaxCPU "mpijob mpirun ./exkatampi2.out $TubeType $SubBand"
#    bsub -o $LOG -J $JOBNAME -n $MaxCPU "mpijob mpirun ./kappa124-2.out $TubeType $SubBand"
}


SelectSubBand () {
    local TubeType=$1
    local SubBand=$2
    case $SubBand in
	all)
	    n=1
	    while [ $n -le $MaxSubBand ]
            do
	      SubmitJob $TubeType $n
	      n=`expr $n + 1`
	    done
	    ;;
	*)
	    if [ $SubBand -le $MaxSubBand ]; then
		SubmitJob $TubeType $SubBand
	    else
		echo "error: subband number is wrong."
		exit
	    fi
	    ;;
    esac
}


SelectTubeType () {
    local TubeType=$1
    local SubBand=$2
    case $TubeType in
	s0|s1|s2)
	    SelectSubBand $TubeType $SubBand
	    ;;
	all)
	    SelectSubBand s0 $SubBand
	    SelectSubBand s1 $SubBand
	    SelectSubBand s2 $SubBand
	    ;;
	*)
	    echo "error: tube type is wrong."
	    exit
	    ;;
    esac
}


SelectTubeType $gTubeType $gSubBand
