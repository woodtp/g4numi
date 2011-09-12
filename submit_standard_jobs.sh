#!/bin/bash
export CONDOR_TMP=/minerva/app/users/condor-tmp/g4numi
export CONDOR_EXEC=/minerva/app/users/condor-exec/g4numi
export X509_USER_PROXY=/scratch/minervadat/grid3/minervadat.proxy.dat
# do first

#export GRID=""
export GRID="-g -opportunistic"
######################## le010z185i ################################
export BEAMCONFIG="le010z185i"
export DOWATER="false"
export WATERCM="3"
export POT="500000" 
export NJOBS="500"
export RUN="1"
export OUTDIR="/minerva/data/flux/g4numi/v1/le010z185i"
export LOGDIR="/minerva/data/flux/g4numi/v1/le010z185i/log"
mkdir -p $OUTDIR
mkdir -p $LOGDIR
export LOGFILE="${LOGDIR}/g4numi_${BEAMCONFIG}_${RUN}_\${PROCESS}.log"
echo $LOGFILE
minerva_jobsub $GRID -e BEAMCONFIG -e DOWATER -e WATERCM -e POT -e RUN -N $NJOBS \
    -dG4NUMI $OUTDIR -L $LOGFILE \
    /grid/fermiapp/minerva/beamsim/g4numi/v1/g4numi_job.sh

####################### le010z-185i ##################################
export BEAMCONFIG="le010z-185i"
export DOWATER="false"
export WATERCM="3"
export POT="500000" 
export NJOBS="500"
export RUN="1"
export OUTDIR="/minerva/data/flux/g4numi/v1/le010z-185i"
export LOGDIR="/minerva/data/flux/g4numi/v1/le010z-185i/log"
mkdir -p $OUTDIR
mkdir -p $LOGDIR
export LOGFILE="${LOGDIR}/g4numi_${BEAMCONFIG}_${RUN}_\${PROCESS}.log"
echo $LOGFILE
minerva_jobsub $GRID -e BEAMCONFIG -e DOWATER -e WATERCM -e POT -e RUN -N $NJOBS \
    -dG4NUMI $OUTDIR -L $LOGFILE \
    /grid/fermiapp/minerva/beamsim/g4numi/v1/g4numi_job.sh


####################### le100z200i #####################################
export BEAMCONFIG="le100z200i"
export DOWATER="false"
export WATERCM="3"
export POT="500000" 
export NJOBS="500"
export RUN="1"
export OUTDIR="/minerva/data/flux/g4numi/v1/le100z200i"
export LOGDIR="/minerva/data/flux/g4numi/v1/le100z200i/log"
mkdir -p $OUTDIR
mkdir -p $LOGDIR
export LOGFILE="${LOGDIR}/g4numi_${BEAMCONFIG}_${RUN}_\${PROCESS}.log"
echo $LOGFILE
minerva_jobsub $GRID -e BEAMCONFIG -e DOWATER -e WATERCM -e POT -e RUN -N $NJOBS \
    -dG4NUMI $OUTDIR -L $LOGFILE \
    /grid/fermiapp/minerva/beamsim/g4numi/v1/g4numi_job.sh


################################### le100z-200i ##########################################
export BEAMCONFIG="le100z-200i"
export DOWATER="false"
export WATERCM="3"
export POT="500000" 
export NJOBS="500"
export RUN="1"
export OUTDIR="/minerva/data/flux/g4numi/v1/le100z-200i"
export LOGDIR="/minerva/data/flux/g4numi/v1/le100z-200i/log"
mkdir -p $OUTDIR
mkdir -p $LOGDIR
export LOGFILE="${LOGDIR}/g4numi_${BEAMCONFIG}_${RUN}_\${PROCESS}.log"
echo $LOGFILE
minerva_jobsub $GRID -e BEAMCONFIG -e DOWATER -e WATERCM -e POT -e RUN -N $NJOBS \
    -dG4NUMI $OUTDIR -L $LOGFILE \
    /grid/fermiapp/minerva/beamsim/g4numi/v1/g4numi_job.sh

################################### le250z200i ###########################################
export BEAMCONFIG="le250z200i"
export DOWATER="false"
export WATERCM="3"
export POT="500000" 
export NJOBS="500"
export RUN="1"
export OUTDIR="/minerva/data/flux/g4numi/v1/le250z200i"
export LOGDIR="/minerva/data/flux/g4numi/v1/le250z200i/log"
mkdir -p $OUTDIR
mkdir -p $LOGDIR
export LOGFILE="${LOGDIR}/g4numi_${BEAMCONFIG}_${RUN}_\${PROCESS}.log"
echo $LOGFILE
minerva_jobsub $GRID -e BEAMCONFIG -e DOWATER -e WATERCM -e POT -e RUN -N $NJOBS \
    -dG4NUMI $OUTDIR -L $LOGFILE \
    /grid/fermiapp/minerva/beamsim/g4numi/v1/g4numi_job.sh

############################## le250z-200i ########################################
export BEAMCONFIG="le250z-200i"
export DOWATER="false"
export WATERCM="3"
export POT="500000" 
export NJOBS="500"
export RUN="1"
export OUTDIR="/minerva/data/flux/g4numi/v1/le250z-200i"
export LOGDIR="/minerva/data/flux/g4numi/v1/le250z-200i/log"
mkdir -p $OUTDIR
mkdir -p $LOGDIR
export LOGFILE="${LOGDIR}/g4numi_${BEAMCONFIG}_${RUN}_\${PROCESS}.log"
echo $LOGFILE
minerva_jobsub $GRID -e BEAMCONFIG -e DOWATER -e WATERCM -e POT -e RUN -N $NJOBS \
    -dG4NUMI $OUTDIR -L $LOGFILE \
    /grid/fermiapp/minerva/beamsim/g4numi/v1/g4numi_job.sh

############################## le010z000i ########################################
export BEAMCONFIG="le010z000i"
export DOWATER="false"
export WATERCM="3"
export POT="500000" 
export NJOBS="500"
export RUN="1"
export OUTDIR="/minerva/data/flux/g4numi/v1/le010z000i"
export LOGDIR="/minerva/data/flux/g4numi/v1/le010z000i/log"
mkdir -p $OUTDIR
mkdir -p $LOGDIR
export LOGFILE="${LOGDIR}/g4numi_${BEAMCONFIG}_${RUN}_\${PROCESS}.log"
echo $LOGFILE
minerva_jobsub $GRID -e BEAMCONFIG -e DOWATER -e WATERCM -e POT -e RUN -N $NJOBS \
    -dG4NUMI $OUTDIR -L $LOGFILE \
    /grid/fermiapp/minerva/beamsim/g4numi/v1/g4numi_job.sh

############################## le150z200i ########################################
export BEAMCONFIG="le150z200i"
export DOWATER="false"
export WATERCM="3"
export POT="500000" 
export NJOBS="500"
export RUN="1"
export OUTDIR="/minerva/data/flux/g4numi/v1/le150z200i"
export LOGDIR="/minerva/data/flux/g4numi/v1/le150z200i/log"
mkdir -p $OUTDIR
mkdir -p $LOGDIR
export LOGFILE="${LOGDIR}/g4numi_${BEAMCONFIG}_${RUN}_\${PROCESS}.log"
echo $LOGFILE
minerva_jobsub $GRID -e BEAMCONFIG -e DOWATER -e WATERCM -e POT -e RUN -N $NJOBS \
    -dG4NUMI $OUTDIR -L $LOGFILE \
    /grid/fermiapp/minerva/beamsim/g4numi/v1/g4numi_job.sh

############################## le150z-200i ########################################
export BEAMCONFIG="le150z-200i"
export DOWATER="false"
export WATERCM="3"
export POT="500000" 
export NJOBS="500"
export RUN="1"
export OUTDIR="/minerva/data/flux/g4numi/v1/le150z-200i"
export LOGDIR="/minerva/data/flux/g4numi/v1/le150z-200i/log"
mkdir -p $OUTDIR
mkdir -p $LOGDIR
export LOGFILE="${LOGDIR}/g4numi_${BEAMCONFIG}_${RUN}_\${PROCESS}.log"
echo $LOGFILE
minerva_jobsub $GRID -e BEAMCONFIG -e DOWATER -e WATERCM -e POT -e RUN -N $NJOBS \
    -dG4NUMI $OUTDIR -L $LOGFILE \
    /grid/fermiapp/minerva/beamsim/g4numi/v1/g4numi_job.sh
