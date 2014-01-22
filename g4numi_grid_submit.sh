#!/bin/bash
#export CONDOR_TMP=/minerva/app/users/condor-tmp/${USER}
#export CONDOR_EXEC=/minerva/app/users/condor-exec/${USER}
#export X509_USER_PROXY=/scratch/minervadat/grid3/minervadat.proxy.dat

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
domain=`uname -n | cut -f 1- -d . | tail -c 9`
if [ -e "/grid/fermiapp/products/minerva/etc/setups.sh" -a $domain  == "fnal.gov" ]; then
    source /grid/fermiapp/products/minerva/etc/setups.sh
    source /grid/fermiapp/products/common/etc/setups.sh
    setup cpn -z /grid/fermiapp/products/common/db
    #setup jobsub_tools -z /grid/fermiapp/products/common/db
    #Switch back to the above line on Monday after Dennis fixes the bug
    setup -t jobsub_tools
#    setup jobsub_tools v1_1t -z /grid/fermiapp/products/common/db/
fi

# make sure that didn't screw up any important paths
source setup_beamsim.sh

######################## le010z185i ################################
#export BEAMCONFIG="le100z-200i"
export BEAMCONFIG="le010z185i"
export DOWATER="false"
export WATERCM="3"
export DOIMPWT="false"
export POT="500000" 
export NJOBS="50"
#export POT="500" 
#export NJOBS="2"
export RUN="2"
export OUTDIR="/minerva/data/users/$USER/flux/ftfp_bert/"
export LOGFILE="${OUTDIR}/g4numi_${BEAMCONFIG}_${RUN}_\${PROCESS}.log"
mkdir -p $OUTDIR
# make it writable for minervaana
chmod 01775 $OUTDIR
echo "Output to $OUTDIR"
export G4NUMIAREA=`pwd`
echo "Will use code and scripts from $G4NUMIAREA"

jobsub $GRID -N $NJOBS -dG4NUMI $OUTDIR -e G4NUMIAREA  -e BEAMCONFIG -e DOWATER -e WATERCM -e POT -e RUN -e DOIMPWT -L $LOGFILE /minerva/app/users/$USER/cmtuser/NumiAna/numisoft/g4numi/g4numi_job.sh 

#
