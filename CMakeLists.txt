# use cmake 2.8 or later
cmake_minimum_required (VERSION 3.2)

project(numi-beam-sim-g4numi)

# cetbuildtools contains our cmake modules
SET ( CETBUILDTOOLS_VERSION $ENV{CETBUILDTOOLS_VERSION} )
if (NOT CETBUILDTOOLS_VERSION)
    MESSAGE (FATAL_ERROR "ERROR: setup cetbuildtools to get the cmake modules")
endif()

set(CMAKE_MODULE_PATH $ENV{CANVAS_DIR}/Modules
                      $ENV{CETBUILDTOOLS_DIR}/Modules
                      ${CMAKE_MODULE_PATH})

include(CetCMakeEnv)
cet_cmake_env()

set_install_root()

cet_set_compiler_flags(DIAGS CAUTIOUS
  #WERROR
  NO_UNDEFINED
  EXTRA_FLAGS -Wno-pedantic -Wno-unused-local-typedefs -Wno-narrowing -Wno-deprecated-declarations -Wno-unused-variable -Wno-missing-braces -Wno-undefined-var-template
  EXTRA_CXX_FLAGS
)

cet_enable_asserts()

# issue 5589: if we are not in a debugging build, disable the dispatch of debug messages in messagefacility
string(TOUPPER ${CMAKE_BUILD_TYPE} BTYPE_UC )
if( ${BTYPE_UC} STREQUAL "OPT" OR
    ${BTYPE_UC} STREQUAL "PROF" OR
    ${BTYPE_UC} STREQUAL "RELEASE" OR
    ${BTYPE_UC} STREQUAL "MINSIZEREL" )
  message ( STATUS "Disabling LOG_DEBUG from messagefacility library" )
  add_definitions(-DML_NDEBUG)
endif()

add_definitions( -DCMAKEBUILD )

cet_report_compiler_flags()

# these are minimum required versions, not the actual product versions
# put products in this list that you want to have set up with
# their environment by ups
find_ups_product( geant4         v4_10_4  )
find_ups_product( canvas_root_io v1_05_02 )
find_ups_product( dk2nudata      v01_08_00g )

# GDMLParser in G4Percistency requires access to XERCES_C headers
find_ups_product( xerces_c )

find_ups_root( root v6_00_00 )

include(ArtDictionary)

# collect the G4 libraries
set( ENV{GEANT4_LIB} $ENV{GEANT4_FQ_DIR}/lib )
#MESSAGE( "Geant4 libraries are in $ENV{GEANT4_LIB}" )
cet_find_library( G4ANALYSIS           NAMES G4analysis          PATHS ENV GEANT4_LIB NO_DEFAULT_PATH )
cet_find_library( G4DIGITS_HITS        NAMES G4digits_hits       PATHS ENV GEANT4_LIB NO_DEFAULT_PATH )
cet_find_library( G4ERROR_PROPAGATION  NAMES G4error_propagation PATHS ENV GEANT4_LIB NO_DEFAULT_PATH )
cet_find_library( G4EVENT              NAMES G4event             PATHS ENV GEANT4_LIB NO_DEFAULT_PATH )
cet_find_library( G4FR                 NAMES G4FR                PATHS ENV GEANT4_LIB NO_DEFAULT_PATH )
cet_find_library( G4GEOMETRY           NAMES G4geometry          PATHS ENV GEANT4_LIB NO_DEFAULT_PATH )
cet_find_library( G4GL2PS              NAMES G4gl2ps             PATHS ENV GEANT4_LIB NO_DEFAULT_PATH )
cet_find_library( G4GLOBAL             NAMES G4global            PATHS ENV GEANT4_LIB NO_DEFAULT_PATH )
cet_find_library( G4GMOCREN            NAMES G4GMocren           PATHS ENV GEANT4_LIB NO_DEFAULT_PATH )
cet_find_library( G4GRAPHICS_REPS      NAMES G4graphics_reps     PATHS ENV GEANT4_LIB NO_DEFAULT_PATH )
cet_find_library( G4INTERCOMS          NAMES G4intercoms         PATHS ENV GEANT4_LIB NO_DEFAULT_PATH )
cet_find_library( G4INTERFACES         NAMES G4interfaces        PATHS ENV GEANT4_LIB NO_DEFAULT_PATH )
cet_find_library( G4MATERIALS          NAMES G4materials         PATHS ENV GEANT4_LIB NO_DEFAULT_PATH )
cet_find_library( G4MODELING           NAMES G4modeling          PATHS ENV GEANT4_LIB NO_DEFAULT_PATH )
cet_find_library( G4OPENGL             NAMES G4OpenGL            PATHS ENV GEANT4_LIB NO_DEFAULT_PATH )
cet_find_library( G4PARMODELS          NAMES G4parmodels         PATHS ENV GEANT4_LIB NO_DEFAULT_PATH )
cet_find_library( G4PARTICLES          NAMES G4particles         PATHS ENV GEANT4_LIB NO_DEFAULT_PATH )
cet_find_library( G4PERSISTENCY        NAMES G4persistency       PATHS ENV GEANT4_LIB NO_DEFAULT_PATH )
cet_find_library( G4PHYSICSLISTS       NAMES G4physicslists      PATHS ENV GEANT4_LIB NO_DEFAULT_PATH )
cet_find_library( G4PROCESSES          NAMES G4processes         PATHS ENV GEANT4_LIB NO_DEFAULT_PATH )
cet_find_library( G4RAYTRACER          NAMES G4RayTracer         PATHS ENV GEANT4_LIB NO_DEFAULT_PATH )
cet_find_library( G4READOUT            NAMES G4readout           PATHS ENV GEANT4_LIB NO_DEFAULT_PATH )
cet_find_library( G4RUN                NAMES G4run               PATHS ENV GEANT4_LIB NO_DEFAULT_PATH )
cet_find_library( G4TRACKING           NAMES G4tracking          PATHS ENV GEANT4_LIB NO_DEFAULT_PATH )
cet_find_library( G4TRACK              NAMES G4track             PATHS ENV GEANT4_LIB NO_DEFAULT_PATH )
cet_find_library( G4TREE               NAMES G4Tree              PATHS ENV GEANT4_LIB NO_DEFAULT_PATH )
cet_find_library( G4VISHEPREP          NAMES G4visHepRep         PATHS ENV GEANT4_LIB NO_DEFAULT_PATH )
cet_find_library( G4VIS_MANAGEMENT     NAMES G4vis_management    PATHS ENV GEANT4_LIB NO_DEFAULT_PATH )
cet_find_library( G4VISXXX             NAMES G4visXXX            PATHS ENV GEANT4_LIB NO_DEFAULT_PATH )
cet_find_library( G4VRML               NAMES G4VRML              PATHS ENV GEANT4_LIB NO_DEFAULT_PATH )
cet_find_library( G4ZLIM               NAMES G4zlib              PATHS ENV GEANT4_LIB NO_DEFAULT_PATH )
# PERSISTENCY needs XERCES_C library
cet_find_library( XERCES_C             NAMES xerces-c            PATHS ENV XERCES_C_LIB NO_DEFAULT_PATH )

set( G4_LIB_LIST ${G4ANALYSIS}
                 ${G4DIGITS_HITS}
                 ${G4ERROR_PROPAGATION}
                 ${G4EVENT}
                 ${G4FR}
                 ${G4GEOMETRY}
                 ${G4GL2PS}
                 ${G4GLOBAL}
                 ${G4GMOCREN}
                 ${G4GRAPHICS_REPS}
                 ${G4INTERCOMS}
                 ${G4INTERFACES}
                 ${G4MATERIALS}
                 ${G4MODELING}
                 ${G4OPENGL}
                 ${G4PARMODELS}
                 ${G4PARTICLES}
                 ${G4PERSISTENCY}
                 ${XERCES_C}
                 ${G4PHYSICSLISTS}
                 ${G4PROCESSES}
                 ${G4RAYTRACER}
                 ${G4READOUT}
                 ${G4RUN}
                 ${G4TRACKING}
                 ${G4TRACK}
                 ${G4TREE}
                 ${G4VISHEPREP}
                 ${G4VIS_MANAGEMENT}
                 ${G4VISXXX}
                 ${G4VRML}
                 ${G4ZLIM}
   )

#MESSAGE("Geant4 library list is ${G4_LIB_LIST}" )

cet_find_library( CLHEPLIB NAMES CLHEP PATHS ENV CLHEP_LIB_DIR NO_DEFAULT_PATH )

# ROOT libraries
set( ENV{ROOTLIB} $ENV{ROOTSYS}/lib )
cet_find_library( ROOTCORE       NAMES  Core       PATHS ENV ROOTLIB NO_DEFAULT_PATH )
cet_find_library( ROOTCINT       NAMES  Cling      PATHS ENV ROOTLIB NO_DEFAULT_PATH )
cet_find_library( ROOTRIO        NAMES  RIO        PATHS ENV ROOTLIB NO_DEFAULT_PATH )
cet_find_library( ROOTNET        NAMES  Net        PATHS ENV ROOTLIB NO_DEFAULT_PATH )
cet_find_library( ROOTIMT        NAMES  Imt        PATHS ENV ROOTLIB NO_DEFAULT_PATH )
cet_find_library( ROOTHIST       NAMES  Hist       PATHS ENV ROOTLIB NO_DEFAULT_PATH )
cet_find_library( ROOTGRAF       NAMES  Graf       PATHS ENV ROOTLIB NO_DEFAULT_PATH )
cet_find_library( ROOTGRAF3D     NAMES  Graf3d     PATHS ENV ROOTLIB NO_DEFAULT_PATH )
cet_find_library( ROOTGPAD       NAMES  Gpad       PATHS ENV ROOTLIB NO_DEFAULT_PATH )
cet_find_library( ROOTTREE       NAMES  Tree       PATHS ENV ROOTLIB NO_DEFAULT_PATH )
cet_find_library( ROOTRINT       NAMES  Rint       PATHS ENV ROOTLIB NO_DEFAULT_PATH )
cet_find_library( ROOTPOSTSCRIPT NAMES  Postscript PATHS ENV ROOTLIB NO_DEFAULT_PATH )
cet_find_library( ROOTMATRIX     NAMES  Matrix     PATHS ENV ROOTLIB NO_DEFAULT_PATH )
cet_find_library( ROOTPHYSICS    NAMES  Physics    PATHS ENV ROOTLIB NO_DEFAULT_PATH )
cet_find_library( ROOTMATHCORE   NAMES  MathCore   PATHS ENV ROOTLIB NO_DEFAULT_PATH )
cet_find_library( ROOTTHREAD     NAMES  Thread     PATHS ENV ROOTLIB NO_DEFAULT_PATH )
set(ROOT_BASIC_LIB_LIST ${ROOTCORE}
                        ${ROOTCINT}
                        ${ROOTRIO}
                        ${ROOTNET}
                        ${ROOTIMT}
                        ${ROOTHIST}
                        ${ROOTGRAF}
                        ${ROOTGRAF3D}
                        ${ROOTGPAD}
                        ${ROOTTREE}
                        ${ROOTRINT}
                        ${ROOTPOSTSCRIPT}
                        ${ROOTMATRIX}
                        ${ROOTPHYSICS}
                        ${ROOTMATHCORE}
                        ${ROOTTHREAD}
)
# source

add_subdirectory(dataProducts)
add_subdirectory(exec)
add_subdirectory(include)
add_subdirectory(macros)
#add_subdirectory(rndm) # nothing in this directory to do
add_subdirectory(plot)
add_subdirectory(src)
add_subdirectory(studies)

# ups - table and config files
add_subdirectory(ups)

# packaging utility
include(UseCPack)

install_scripts()

# install root macros into the ups product
file(GLOB root_macros *.C)
install( FILES ${root_macros} DESTINATION ${product}/${version}/macros/root )
