# FOR CLASSES_DEF.XML!!
include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

include_directories( $ENV{ROOT_INC}                     )
include_directories( $ENV{GEANT4_FQ_DIR}/include/Geant4 )

cet_make(LIBRARY_NAME   dataProducts
         LIBRARIES      ${G4_LIB_LIST}
                        ${CLHEPLIB}
                        ${ROOT_BASIC_LIB_LIST}
         NO_DICTIONARY
         )

build_dictionary(DICTIONARY_LIBRARIES dataProducts)

install_source()
install_headers()