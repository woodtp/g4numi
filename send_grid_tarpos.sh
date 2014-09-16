#!/bin/bash

# $Header: /cvs/projects/numi-beam-sim/numi-beam-sim/g4numi/Attic/send_grid_tarpos.sh,v 1.1.2.1 2014/09/16 16:37:36 laliaga Exp $

# for interactive running use GRID=""
#export GRID=""
# for production account, non group-writable, use GRID="-g"
#export GRID="-g"
# for submitting as myself, use GRID="-g --use_gftp"
export GRID="-g --opportunistic --use_gftp"

#############################
#General variables:
export DOWATER="false"
export WATERCM="3"
export DOIMPWT="true"
export POT="20000" 
export NJOBS="5"
export RUN="5"

BEAMS=( "le010z185i" "le010z185i" "le010z185i" "le010z185i" "le010z-185i" "le010z-185i" "le010z-185i" "le010z000i" "le100z200i" "le100z200i" "le100z-200i" "le100z-200i" "le250z200i" "le250z200i")
MODES=( "1" "7" "9" "13" "0" "5" "10" "6" "2" "11" "3" "12" "4" "8")

len=${#BEAMS[*]}  
echo "len"

for (( i=0; i<len; i++ )) 
do
  export BEAMCONFIG=${BEAMS[${i}]}
  export PLAYLIST=${MODES[${i}]}
  export OUTDIR="/minerva/data/users/$USER/flux/TEST/${BEAMCONFIG}/playlist${PLAYLIST}"
  export LOGDIR="${OUTDIR}/log"
  export LOGFILE="${LOGDIR}/g4numi_${BEAMCONFIG}_${PLAYLIST}_${RUN}_\${PROCESS}.log"
  mkdir -p $OUTDIR
  chmod 01775 $OUTDIR
  mkdir -p $LOGDIR
  chmod 01775 $LOGDIR

  echo "Output to $OUTDIR"
  export G4NUMIAREA=`pwd`
  echo "Will use code and scripts from $G4NUMIAREA"  

  jobsub $GRID -N $NJOBS -dG4NUMI $OUTDIR -e G4NUMIAREA  -e BEAMCONFIG -e DOWATER -e WATERCM -e POT -e RUN -e DOIMPWT -e PLAYLIST -L $LOGFILE $BEAMSIM/g4numi_job.sh  
done
