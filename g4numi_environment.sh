# these variables to customize based on actual installation location
export G4NUMIDIR=/minerva/app/users/${USER}/cmtuser/Minerva_g4numi/NumiAna/numisoft/g4numi
export G4BASEDIR=$G4NUMIDIR/g4install
export G4INSTALL=$G4INSTALL/geant4.9.2.p03
export CLHEP_BASE_DIR=/afs/fnal.gov/ups/clhep/v2_0_4_5/Linux+2.6-2.5

# relative paths
export G4WORKDIR=$G4NUMIDIR
export G4INCLUDE=$G4INSTALL/include
export G4SYSTEM=Linux-g++
export G4LIB=$G4INSTALL/lib
export G4LEVELGAMMADATA=$G4INSTALL/data/PhotonEvaporation2.0
export G4RADIOACTIVEDATA=$G4INSTALL/data/RadioactiveDecay3.2
export G4LEDATA=$G4INSTALL/data/G4EMLOW6.2
export G4NEUTRONHPDATA=$G4INSTALL/data/G4NDL3.13
export G4ABLADATA=$G4INSTALL/data/G4ABLA3.0
export CLHEP_INCLUDE_DIR=$CLHEP_BASE_DIR/include
export CLHEP_LIB_DIR=$CLHEP_BASE_DIR/lib
export CLHEP_LIB=CLHEP
export G4UI_USE_TCSH=1
export G4LIB_BUILD_SHARED=1
export G4LIB_USE_GRANULAR=1
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:$G4INSTALL/lib/${G4SYSTEM}
export GEANT4_DIR=$G4INSTALL
export GEANT4SOURCE_DIR=$G4BASEDIR
