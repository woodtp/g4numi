# FOR CLASSES_DEF.XML!!
include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

include_directories( $ENV{CLHEP_INC}                    )
include_directories( $ENV{GEANT4_FQ_DIR}/include/Geant4 )
include_directories( $ENV{ROOT_INC}                     )
include_directories( ${CMAKE_SOURCE_DIR}/include        )

cet_make(LIBRARY_NAME  src
         LIBRARIES     dataProducts
                       ${G4_LIB_LIST}
                       ${CLHEPLIB}
	               ${ROOT_BASIC_LIB_LIST}
         )

install_source()
