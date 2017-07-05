#!/bin/bash

# $Header: /cvs/projects/numi-beam-sim/numi-beam-sim/g4numi/Attic/g4numi_grid_submit.sh,v 1.1.2.8 2017/07/05 17:13:33 bmesserl Exp $

export GRID="--environment=LD_LIBRARY_PATH --OS=SL6 -g --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC --role=Analysis -G minerva "
export MEMORY="--memory 500MB "

########################### g4numi configuration ###########################
export DO_HORN1_OLD_GEOMETRY="false"       #"false" the new horn is considered default (01/2016)
export DO_HORN1_FINE_SEGMENTATION="false"  #"false" works for both old and new horn geometry
export DO_TARGET_WATER_LAYER="false"       #"false" TARGET water layer
export IMPORTANCE_WEIGHTING="true"         #"true"
export TARGET_WATER_CM="3"                 #"3" TARGET water layer
export HORN1_POSITION_X="0"                #"0.0" #cm
export HORN1_POSITION_Y="0"                #"0.0" #cm
export HORN2_POSITION_X="0"                #"0.0" #cm
export HORN_WATER_MM="1.0"                 #"1.0" #mm
export TARGET_POSITION_X="0.0"             #"0.0" #cm
export TARGET_POSITION_Y="0.0"             #"0.0" #cm
export TARGET_POSITION_Z="-143.3"          #"-143.3" #cm
export BEAM_POSITION_X="0.0"               #"0" #m
export BEAM_POSITION_Y="0.0"               #"0" #m
export BEAM_SPOTSIZE_X="1.4"               #"1.4" #mm. obtained from typical ME-era ACNET values (~1.3),
export BEAM_SPOTSIZE_Y="1.4"               #"1.4" #mm. properly scaled by 1.08 (see Hylen numi ops 02/16/16)
export POT="400000"                        #"400000"
export NJOBS="1020"                        #"1020"
export RUN="1"
export BEAMCONFIG="me000z200i"             #the whole beamconfig thing is not used for ME. Should still work for LE.
export PLAYLIST="minervame"
export OUTDIR="/pnfs/minerva/persistent/users/${USER}/flux/test/"
export LOGFILE="${OUTDIR}/g4numi_${BEAMCONFIG}_${RUN}_\${PROCESS}.log"
mkdir -p $OUTDIR
# make it writable for minervaana
chmod 01775 $OUTDIR
echo "Output to $OUTDIR"
export G4NUMIAREA=`pwd`
echo "Will use code and scripts from $G4NUMIAREA"

############################### submit the job #############################

jobsub_submit $GRID $MEMORY -N $NJOBS -dG4NUMI $OUTDIR -e G4NUMIAREA  -e BEAMCONFIG -e DO_TARGET_WATER_LAYER -e TARGET_WATER_CM -e HORN1_POSITION_X -e HORN1_POSITION_Y -e HORN2_POSITION_X -e TARGET_POSITION_X -e TARGET_POSITION_Y -e TARGET_POSITION_Z -e HORN_WATER_MM -e BEAM_POSITION_X -e BEAM_POSITION_Y -e BEAM_SPOTSIZE_X -e BEAM_SPOTSIZE_Y -e DO_HORN1_OLD_GEOMETRY -e DO_HORN1_FINE_SEGMENTATION -e POT -e RUN -e IMPORTANCE_WEIGHTING -e PLAYLIST -L $LOGFILE file://$BEAMSIM/g4numi_job.sh
