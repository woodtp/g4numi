#!/bin/bash
#-*-Shell-Script-*- #

# $Header: /cvs/projects/numi-beam-sim/numi-beam-sim/g4numi/Attic/setup_beamsim.sh,v 1.1.2.15 2015/11/28 05:08:58 laliaga Exp $


setup_beamsim(){
    . "/grid/fermiapp/products/minerva/etc/setups.sh"
    local TOP=${PWD}
#   setup -f Linux+2.6-2.5 gcc v3_4_3

# This line sets the modified version of geant4 (that contains a class to store hadronic interaction 
# kinematics):
    . /grid/fermiapp/minerva/beamsim/geant4/geant4.9.2.p03_mod/env.sh 

#    setup -q g77-OpenGL -f Linux+2.6 geant4 v4_9_2_p03

    setup -q debug -f Linux+2.6-2.5 root v5_30_00

    if [ -e "/grid/fermiapp/minerva/condor/setup.minerva.condor.sh" ]; then
	source /grid/fermiapp/minerva/condor/setup.minerva.condor.sh
    fi


    source /grid/fermiapp/products/minerva/etc/setups.sh

    # setup for jobsub client
    # according to the prescription in Mike Kirby's talk
    # minerva doc-10551, Dec 2014    
    source /grid/fermiapp/products/common/etc/setups.sh
    setup jobsub_client
    setup ifdhc

    export IFDH_BASE_URI="http://samweb-minerva.fnal.gov:20004/sam/minerva/api"

    export BEAMSIM="${TOP}"
    export G4NUMIVER="v6"
    setup g4photon
    export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:$G4LIB/plists/Linux-g++"
#    export LD_LIBRARY_PATH="${PWD}/tmpSolution:${LD_LIBRARY_PATH}"


    export G4WORKDIR="${TOP}"
    echo "G4WORKDIR is ${G4WORKDIR}"

}
setup_beamsim

