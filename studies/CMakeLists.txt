include_directories( $ENV{GEANT4_FQ_DIR}/include/Geant4 )
include_directories( $ENV{CLHEP_INC}                    )

cet_make(LIBRARY_NAME  studies
         LIBRARIES     ${G4_LIB_LIST}
         )


install_source()
install_headers()

