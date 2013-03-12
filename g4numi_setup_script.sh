if [  -f "/afs/fnal.gov/ups/etc/setups.sh" ]
then
    . "/afs/fnal.gov/ups/etc/setups.sh"
    if ups exist login
    then
        setup login
    fi
fi  

ulimit -c unlimited

sm()
{
   export MINOS_SETUP_DIR="/grid/fermiapp/minos/minossoft/setup"
   source /afs/fnal.gov/ups/etc/setups.sh
   source ${MINOS_SETUP_DIR}/setup_minossoft_FNALU.sh $*
   export srt=$SRT_PUBLIC_CONTEXT
}

smc()
{
   if [ ! -e ".base_release" ]
   then
       echo File .base_release does not exist in current directory; return 0
   fi

   SNAP_RELEASE=`cat .base_release`
   TEST_RELEASE=`pwd`
   echo Setting up local test release ${TEST_RELEASE}
   echo Setting up snapshot release ${SNAP_RELEASE}
   MINOS_SETUP_DIR="/grid/fermiapp/minos/minossoft/setup"
   source /afs/fnal.gov/ups/etc/setups.sh
   source ${MINOS_SETUP_DIR}/setup_minossoft_FNALU.sh ${SNAP_RELEASE} $*
   srt_setup -a

   export BEAM_DATA_DIR="/minerva/data/numi/beamdata/"
   export FLUX_FILES_DIR="/minerva/data/numi/fluxfiles/"	
}

setup_g4numi_prod()
{
    sm -r R2.3 -O

    setup -q GCC_3_4_3-OpenGL -f Linux+2.4 geant4 v4_9_2_p03

    setup g4photon
    export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:$G4LIB/plists/Linux-g++"

}

