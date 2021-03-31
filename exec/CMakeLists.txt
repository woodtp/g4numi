include_directories( $ENV{GEANT4_FQ_DIR}/include/Geant4 )
include_directories( $ENV{CLHEP_INC}                    )
include_directories( ${CMAKE_SOURCE_DIR}/include        )

cet_make_exec(g4numi
              SOURCE     g4numi.cc
              LIBRARIES  ${G4_LIB_LIST}
                         src
             )

cet_make_exec(g4na49
              SOURCE     g4na49.cc
              LIBRARIES  ${G4_LIB_LIST}
                         src
             )
