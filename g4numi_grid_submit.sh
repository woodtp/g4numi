#!/bin/bash

# $Header: /cvs/projects/numi-beam-sim/numi-beam-sim/g4numi/Attic/g4numi_grid_submit.sh,v 1.1.2.6 2015/02/25 20:32:33 drut1186 Exp $

# for interactive running use GRID=""
#export GRID=""
# for production account, non group-writable, use GRID="-g"
#export GRID="-g"
# for submitting as myself, use GRID="-g --use_gftp"
#export for the old job submission method kept here just in case
#export GRID="--environment=LD_LIBRARY_PATH --OS=SL6 -g --opportunistic --use_gftp "
export GRID="--environment=LD_LIBRARY_PATH --OS=SL6 -g --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC --use_gftp --role=Analysis -G minerva "


########################### g4numi configuration ###########################
export BEAMCONFIG="le010z185i"
export PLAYLIST="minerva13"
export DOWATER="false"
export WATERCM="3"
export DOIMPWT="true"
export POT="1000" 
export NJOBS="2"
export RUN="1"
export OUTDIR="/minerva/data/users/$USER/flux/sl6test/"
export LOGFILE="${OUTDIR}/g4numi_${BEAMCONFIG}_${RUN}_\${PROCESS}.log"
mkdir -p $OUTDIR
# make it writable for minervaana
chmod 01775 $OUTDIR
echo "Output to $OUTDIR"
export G4NUMIAREA=`pwd`
echo "Will use code and scripts from $G4NUMIAREA"

############################### submit the job #############################

#old submission method kept here in case 
#jobsub $GRID -N $NJOBS -dG4NUMI $OUTDIR -e G4NUMIAREA  -e BEAMCONFIG -e DOWATER -e WATERCM -e POT -e RUN -e DOIMPWT -e PLAYLIST -L $LOGFILE $BEAMSIM/g4numi_job.sh

jobsub_submit $GRID -N $NJOBS -dG4NUMI $OUTDIR -e G4NUMIAREA  -e BEAMCONFIG -e DOWATER -e WATERCM -e POT -e RUN -e DOIMPWT -e PLAYLIST -L $LOGFILE file://$BEAMSIM/g4numi_job.sh
