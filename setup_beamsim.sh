#!/bin/bash
#-*-Shell-Script-*- #

# $Header: /cvs/projects/numi-beam-sim/numi-beam-sim/g4numi/Attic/setup_beamsim.sh,v 1.1.2.16 2017/11/01 17:33:41 bmesserl Exp $


setup_beamsim(){

    local TOP=${PWD}
    export G4WORKDIR="${TOP}"
    echo "G4WORKDIR is ${G4WORKDIR}"

    # geant4
    # We use a modified version "that contains a class to store hadronic 
    #   interaction kinematics. The version is geant4.9.2.p03, which is now
    #   setup with /grid/fermiapp/minerva/beamsim/x86_64/geant4/env.sh
    source /grid/fermiapp/minerva/beamsim/x86_64/geant4/env.sh

    source /grid/fermiapp/minerva/beamsim/x86_64/root/bin/thisroot.sh

    # setup for jobsub client
    # according to the prescription in Mike Kirby's talk
    # minerva doc-10551, Dec 2014    
    source /grid/fermiapp/products/common/etc/setups.sh
    setup jobsub_client
    setup ifdhc
    export IFDH_GRIDFTP_EXTRA="-st 10"
    export IFDH_CP_MAXRETRIES=2
    export IFDH_BASE_URI="http://samweb-minerva.fnal.gov:20004/sam/minerva/api"

    export BEAMSIM="${TOP}"
    export G4NUMIVER="v6"
    setup g4photon
    export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:$G4LIB/plists/Linux-g++"

}
setup_beamsim
