#!/bin/bash
#-*-Shell-Script-*- #

# $Header: /cvs/projects/numi-beam-sim/numi-beam-sim/g4numi/Attic/setup_beamsim.sh,v 1.1.2.13 2015/02/16 00:09:06 laliaga Exp $


setup_beamsim(){
    . "/grid/fermiapp/products/minerva/etc/setups.sh"
    local TOP=${PWD}
#   setup -f Linux+2.6-2.5 gcc v3_4_3

# For using the modified version of geant4 (that contains a class to store hadronic interaction 
# kinematics), you need to add this line:   . /grid/fermiapp/minerva/beamsim/geant4/geant4.9.2.p03_mod/env.sh 
# and comment out the next line (which defines the ups version). 
# Also, you need to use the processor directive to define USEMODGEANT4
    setup -q g77-OpenGL -f Linux+2.6 geant4 v4_9_2_p03
#    . /grid/fermiapp/minerva/beamsim/geant4/geant4.9.2.p03_mod/env.sh 
    setup -q debug -f Linux+2.6-2.5 root v5_30_00

    if [ -e "/grid/fermiapp/minerva/condor/setup.minerva.condor.sh" ]; then
	source /grid/fermiapp/minerva/condor/setup.minerva.condor.sh
    fi


    source /grid/fermiapp/products/minerva/etc/setups.sh
    source /grid/fermiapp/products/common/etc/setups.sh
    setup cpn -z /grid/fermiapp/products/common/db
    setup jobsub_tools v1_2q -z /grid/fermiapp/products/common/db 


    export BEAMSIM="${TOP}"
    export G4NUMIVER="v5"
    setup g4photon
    export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:$G4LIB/plists/Linux-g++"
#    export LD_LIBRARY_PATH="${PWD}/tmpSolution:${LD_LIBRARY_PATH}"


    export G4WORKDIR="${TOP}"
    echo "G4WORKDIR is ${G4WORKDIR}"

}
setup_beamsim

