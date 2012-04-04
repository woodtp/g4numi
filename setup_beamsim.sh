#!/bin/bash
#-*-Shell-Script-*- #

echo "You should run now setup_beamsim -w 1"
setup_beamsim(){
    . "/grid/fermiapp/products/minerva/etc/setups.sh"
    local TOP=${HOME}
    local SETWDIR=0
    while getopts "w" opt
      do  case $opt in
	  w) SETWDIR=1 ;;    
	  h) dohelp=1 ;;      
	  [?]) echo "bad command line option"
	      exit 1;;
      esac
    done
    shift ${OPTIND-1}
    
    if [ "$HOSTNAME" == "if01.fnal.gov" ]; then
	echo 'Setting up MINERVA batch submission using minerva_jobsub'
	source /grid/fermiapp/minerva/condor/setup.minerva.condor.sh
	echo ''
    fi
	setup -q g77-OpenGL -f Linux+2.6 geant4 v4_9_2_p03
	setup -f Linux+2.6 -q GCC_3_4_6 root v5_22_00j
# 	setup python # sets up the local version....
	setup g4photon
	export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:$G4LIB/plists/Linux-g++"
        export LD_LIBRARY_PATH="${PWD}/tmpSolution:${LD_LIBRARY_PATH}"

# people who just want to run g4numi probably do not 
# need to set G4WORKDIR
    if (( SETWDIR==1 )); then
	export G4WORKDIR="${TOP}"
	echo "G4WORKDIR is ${G4WORKDIR}"
    else
	unset G4WORKDIR
    fi
    
    unset OPTSTRING
    unset OPTIND
}
setup_beamsim -w 1
