#!/bin/bash

export GRID="--environment=LD_LIBRARY_PATH --OS=SL6 -g --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC --use_gftp --role=Analysis -G nova "
#############################
#General:
export POT="1000" 
export NJOBS="5"
export RUN="6"
#Beam:
export BEAM_SPOTSIZE_X="1.3"
export BEAM_SPOTSIZE_Y="1.3"
export BEAM_POSITION_X="0"
export BEAM_POSITION_Y="0"
#Target:
export DO_TARGET_WATER="false"
export TARGET_WATER_CM="3"
export TARGET_POSITION_X="0"
export TARGET_POSITION_Y="0"
export TARGET_POSITION_Z="-143.3"
#Horns:
export HORN_WATER_MM="1"
export DO_HORN1_OLD_GEOMETRY="false"
export DO_HORN1_FINE_SEGMENTATION="false"
export HORN1_POSITION_X="0"
export HORN1_POSITION_Y="0"
export HORN1_POSITION_Z="3"
export HORN2_POSITION_X="0"
export HORN2_POSITION_Y="0"

BEAMS=( "me000z200i" "me000z-200i" )
MODES=( "minervame" "minervame" )


len=${#BEAMS[*]}  
echo "Number of beam modes to run: ${len}"

for (( i=0; i<len; i++ )) 
do
  export BEAMCONFIG=${BEAMS[${i}]}
  export PLAYLIST=${MODES[${i}]}
  export OUTDIR="/pnfs/nova/scratch/users/$USER/test_g4numi/${BEAMCONFIG}/${PLAYLIST}"
  export LOGDIR="${OUTDIR}/log"
  export LOGFILE="${LOGDIR}/g4numi_${BEAMCONFIG}_${PLAYLIST}_${RUN}_\${PROCESS}.log"
  mkdir -p $OUTDIR
  chmod 01775 $OUTDIR
  mkdir -p $LOGDIR
  chmod 01775 $LOGDIR

  echo "Output to $OUTDIR"
  export G4NUMIAREA=`pwd`
  echo "Will use code and scripts from $G4NUMIAREA"  

  jobsub_submit $GRID -N $NJOBS -dG4NUMI $OUTDIR -e G4NUMIAREA  -e HOME -e BEAMCONFIG -e BEAM_SPOTSIZE_X -e BEAM_SPOTSIZE_Y -e BEAM_POSITION_X -e BEAM_POSITION_Y -e DO_TARGET_WATER -e TARGET_WATER_CM -e TARGET_POSITION_X -e TARGET_POSITION_Y -e TARGET_POSITION_Z -e DO_HORN1_OLD_GEOMETRY -e DO_HORN1_FINE_SEGMENTATION -e HORN_WATER_MM -e HORN1_POSITION_X -e HORN1_POSITION_Y -e HORN1_POSITION_Z -e HORN2_POSITION_X -e HORN2_POSITION_Y -e POT -e RUN -e PLAYLIST -L $LOGFILE file://$BEAMSIM/g4numi_job.sh

done
