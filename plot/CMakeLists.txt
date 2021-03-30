
cet_make(LIBRARY_NAME  plot
         EXCLUDE       compare_flux.C
         LIBRARIES     ${G4_LIB_LIST}
	               ${ROOT_BASIC_LIB_LIST}
	               
         )

install_source()
install_headers()

file(GLOB root_macros *.C)
install( FILES ${root_macros} DESTINATION ${product}/${version}/macros/root )
