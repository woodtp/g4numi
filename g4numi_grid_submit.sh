#!/bin/bash

# $Header: /cvs/projects/numi-beam-sim/numi-beam-sim/g4numi/Attic/g4numi_grid_submit.sh,v 1.1.2.3 2014/02/14 20:16:05 kordosky Exp $

# for interactive running use GRID=""
#export GRID=""
# for production account, non group-writable, use GRID="-g"
#export GRID="-g"
# for submitting as myself, use GRID="-g --use_gftp"
export GRID="-g --use_gftp"
#export GRID="-g --opportunistic"

########################### g4numi configuration ###########################
export BEAMCONFIG="le010z185i"
export DOWATER="false"
export WATERCM="3"
export DOIMPWT="false"
export POT="500000" 
export NJOBS="10"
export RUN="1"
export OUTDIR="/minerva/data/users/$USER/flux/targ_pos/"
export LOGFILE="${OUTDIR}/g4numi_${BEAMCONFIG}_${RUN}_\${PROCESS}.log"
mkdir -p $OUTDIR
# make it writable for minervaana
chmod 01775 $OUTDIR
echo "Output to $OUTDIR"
export G4NUMIAREA=`pwd`
echo "Will use code and scripts from $G4NUMIAREA"

############################### submit the job #############################

jobsub $GRID -N $NJOBS -dG4NUMI $OUTDIR -e G4NUMIAREA  -e BEAMCONFIG -e DOWATER -e WATERCM -e POT -e RUN -e DOIMPWT -L $LOGFILE $BEAMSIM/g4numi_job.sh

