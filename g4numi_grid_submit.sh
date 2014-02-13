#!/bin/bash

# $Header: /cvs/projects/numi-beam-sim/numi-beam-sim/g4numi/Attic/g4numi_grid_submit.sh,v 1.1.2.2 2014/02/13 21:12:03 kordosky Exp $

# for interactive running use GRID=""
#export GRID=""
# for production account, non group-writable, use GRID="-g"
#export GRID="-g"
# for submitting as myself, use GRID="-g --use_gftp"
export GRID="-g --use_gftp"
#export GRID="-g --opportunistic"

if [ -e "/grid/fermiapp/minerva/condor/setup.minerva.condor.sh" ]; then
    source /grid/fermiapp/minerva/condor/setup.minerva.condor.sh
fi

# set up fermi products and SAM, DCAP if at Fermilab
source /grid/fermiapp/products/minerva/etc/setups.sh
source /grid/fermiapp/products/common/etc/setups.sh
setup cpn -z /grid/fermiapp/products/common/db
#setup jobsub_tools
setup jobsub_tools v1_2q 
#    setup jobsub_tools v1_1t -z /grid/fermiapp/products/common/db/


# make sure that didn't screw up any important paths
source setup_beamsim.sh

######################## le010z185i ################################
#export BEAMCONFIG="le100z-200i"
export BEAMCONFIG="le010z185i"
export DOWATER="false"
export WATERCM="3"
export DOIMPWT="false"
#export POT="500000" 
#export NJOBS="50"
export POT="500" 
export NJOBS="2"
export RUN="2"
export OUTDIR="/minerva/data/users/$USER/flux/ftfp_bert2/"
export LOGFILE="${OUTDIR}/g4numi_${BEAMCONFIG}_${RUN}_\${PROCESS}.log"
mkdir -p $OUTDIR
# make it writable for minervaana
chmod 01775 $OUTDIR
echo "Output to $OUTDIR"
export G4NUMIAREA=`pwd`
echo "Will use code and scripts from $G4NUMIAREA"

jobsub $GRID -N $NJOBS -dG4NUMI $OUTDIR -e G4NUMIAREA  -e BEAMCONFIG -e DOWATER -e WATERCM -e POT -e RUN -e DOIMPWT -L $LOGFILE $BEAMSIM/g4numi_job.sh

