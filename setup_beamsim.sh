#!/bin/bash
#-*-Shell-Script-*- #

# $Header: /cvs/projects/numi-beam-sim/numi-beam-sim/g4numi/Attic/setup_beamsim.sh,v 1.1.2.8 2014/02/13 21:12:03 kordosky Exp $


setup_beamsim(){
    . "/grid/fermiapp/products/minerva/etc/setups.sh"
    local TOP=${PWD}
#   setup -f Linux+2.6-2.5 gcc v3_4_3
# Comment OUT the line below if you are using a local geant4 - g4numi_environment.sh will do this correctly for you
    setup -q g77-OpenGL -f Linux+2.6 geant4 v4_9_2_p03
    setup -f Linux+2.6 -q GCC_3_4_6 root v5_22_00j

    export BEAMSIM="${TOP}"
    export G4NUMIVER="v4"
    setup g4photon
    export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:$G4LIB/plists/Linux-g++"
#    export LD_LIBRARY_PATH="${PWD}/tmpSolution:${LD_LIBRARY_PATH}"


    export G4WORKDIR="${TOP}"
    echo "G4WORKDIR is ${G4WORKDIR}"

}
setup_beamsim

