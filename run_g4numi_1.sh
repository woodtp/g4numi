#!/bin/bash
clear

sm()
{
   export MINOS_SETUP_DIR="/grid/fermiapp/minos/minossoft/setup"
   source /afs/fnal.gov/ups/etc/setups.sh
   #export MINOS_SETUP_DIR="/afs/fnal.gov/files/code/e875/general/minossoft/setup"
   #source /afs/fnal.gov/files/code/e875/general/ups/etc/setups.sh
   source ${MINOS_SETUP_DIR}/setup_minossoft_FNALU.sh $*
   export srt=$SRT_PUBLIC_CONTEXT
}
setup_g4numi_prod()
{
    sm S10-11-11-R2-05 -O
    
    setup -q g77-OpenGL -f Linux+2.6 geant4 v4_9_2_p03

    export G4WORKDIR="/minerva/app/users/mjerkins/geant4"
    echo G4WORKDIR is ${G4WORKDIR}

    setup g4photon
    export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:$G4LIB/plists/Linux-g++"

}
setup_g4numi_prod


./g4numi macros/G4Nevents_le010z185i_${PROCESS}.mac >& /minerva/app/users/condor-tmp/mjerkins/logerr_le010z185i_${PROCESS}.err > /minerva/app/users/condor-tmp/mjerkins/logfile_le010z185i_${PROCESS}.log


