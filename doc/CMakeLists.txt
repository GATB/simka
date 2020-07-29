
find_package (Doxygen)

if (DOXYGEN_FOUND)

    CONFIGURE_FILE (${CMAKE_CURRENT_SOURCE_DIR}/doxygen/gatb-simka.doxyfile  ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile  @ONLY)
    
    ADD_CUSTOM_TARGET (doc-gatb-simka 
      ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
    )
  
endif (DOXYGEN_FOUND)
