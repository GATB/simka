project(simka)

cmake_minimum_required(VERSION 2.6)

################################################################################
# The version number.
################################################################################
SET (gatb-tool_VERSION_MAJOR 1)
SET (gatb-tool_VERSION_MINOR 5)
SET (gatb-tool_VERSION_PATCH 3)

IF (DEFINED MAJOR)
    SET (gatb-tool_VERSION_MAJOR ${MAJOR})
ENDIF()
IF (DEFINED MINOR)
    SET (gatb-tool_VERSION_MINOR ${MINOR})
ENDIF()
IF (DEFINED PATCH)
    SET (gatb-tool_VERSION_PATCH ${PATCH})
ENDIF()

set (gatb-tool-version ${gatb-tool_VERSION_MAJOR}.${gatb-tool_VERSION_MINOR}.${gatb-tool_VERSION_PATCH})

# However, continuous integration has priority over local compilation
IF (DEFINED JENKINS_TAG)
    SET (gatb-tool-version ${JENKINS_TAG})
ENDIF()

################################################################################
# Define cmake modules directory
################################################################################
SET (GATB_CORE_HOME  ${PROJECT_SOURCE_DIR}/thirdparty/gatb-core/gatb-core)
SET (CMAKE_MODULE_PATH ${GATB_CORE_HOME}/cmake)

################################################################################
# THIRD PARTIES
################################################################################

# We don't want to install some GATB-CORE artifacts
 SET (GATB_CORE_EXCLUDE_TOOLS     1)
 SET (GATB_CORE_EXCLUDE_TESTS     1)
 SET (GATB_CORE_EXCLUDE_EXAMPLES  1)

# GATB CORE
include (GatbCore)

################################################################################
# TOOL
################################################################################

# we get compilation definitions from the gatb-core part
add_definitions (${gatb-core-flags})

# we add a new compilation variable
if (PRINTALL)
SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -DPRINTALL" )
endif()

# we give the headers directories from :
#       - from project source
#       - from GATB-CORE source
#       - from simka source
include_directories (include ${gatb-core-includes} ${PROJECT_SOURCE_DIR}/src/core ${PROJECT_SOURCE_DIR}/src/minikc ${PROJECT_SOURCE_DIR}/src)

# we generate one file per template specialization
FOREACH (KSIZE ${gatb-core-klist})
    configure_file (
        ${PROJECT_SOURCE_DIR}/src/core/SimkaAlgorithmTemplate.cpp.in
        ${PROJECT_BINARY_DIR}/src/core/template/SimkaAlgorithmTemplate_${KSIZE}.cpp
    )
ENDFOREACH ()

# we define the files to be compiled
file (GLOB_RECURSE  ProjectFiles  src/core/Simka*  ${PROJECT_BINARY_DIR}/src/core/template/*.cpp)
file (GLOB_RECURSE  SimkaMinFiles  src/simkaMin/MurmurHash3.h src/simkaMin/MurmurHash3.cpp src/simkaMin/*.hpp)#   ${PROJECT_BINARY_DIR}/src/core/template/*.cpp)

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
set(PROJECT_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/bin)

add_executable        (simka  src/SimkaPotara.cpp ${ProjectFiles})
target_link_libraries (simka  ${gatb-core-libraries})

add_executable        (simkaCountProcess  src/minikc/SimkaCountProcess.cpp ${ProjectFiles})
target_link_libraries (simkaCountProcess  ${gatb-core-libraries})
add_executable        (simkaCount  src/SimkaCount.cpp ${ProjectFiles})
target_link_libraries (simkaCount  ${gatb-core-libraries})

add_executable        (simkaMerge  src/SimkaMerge.cpp ${ProjectFiles})
target_link_libraries (simkaMerge  ${gatb-core-libraries})

add_executable        (simkaMinCore  src/simkaMin/SimkaMin.cpp ${SimkaMinFiles})
target_link_libraries (simkaMinCore  ${gatb-core-libraries})

################################################################################
#  DOCUMENTATION GENERATION 
################################################################################
IF (EXISTS "${PROJECT_SOURCE_DIR}/doc")
    ADD_SUBDIRECTORY(doc EXCLUDE_FROM_ALL)
ENDIF()

################################################################################
#  PACKAGING
################################################################################

SET (CPACK_PACKAGE_DESCRIPTION_SUMMARY  "gatb-tool ${PROJECT_NAME}")
SET (CPACK_PACKAGE_VENDOR               "Genscale team (INRIA)")
SET (CPACK_PACKAGE_VERSION_MAJOR        "${gatb-tool_VERSION_MAJOR}")
SET (CPACK_PACKAGE_VERSION_MINOR        "${gatb-tool_VERSION_MINOR}")
SET (CPACK_PACKAGE_VERSION_PATCH        "${gatb-tool_VERSION_PATCH}")
SET (CPACK_PACKAGE_VERSION              "${gatb-tool-version}")

# We chose the kind of archive we want to generate
SET (CPACK_GENERATOR            "TGZ")
SET (CPACK_SOURCE_GENERATOR     "TGZ")

# We ignore unwanted files for the source archive
SET (CPACK_SOURCE_IGNORE_FILES
    "^${PROJECT_SOURCE_DIR}/\\.git/"     ;
    "^${PROJECT_SOURCE_DIR}/\\.gitmodules" ;
    "^${PROJECT_SOURCE_DIR}/\\.gitignore" ;
    "^${PROJECT_SOURCE_DIR}/build/"  ;
    "^${PROJECT_SOURCE_DIR}/dependency/"  ;
    "^${GATB_CORE_HOME}/\\.cproject" ;
    "^${GATB_CORE_HOME}/\\.git/"     ;
    "^${GATB_CORE_HOME}/\\.project"  ;
    "^${GATB_CORE_HOME}/\\.gitignore"
)

# For creating the BINARY package we include the files we want
INSTALL (DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin       DESTINATION .)
INSTALL (DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/example   DESTINATION .)
INSTALL (DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/simkaMin   DESTINATION .)
INSTALL (DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/scripts   DESTINATION . FILES_MATCHING REGEX ".*\\.(py|r)$"  PATTERN "jenkins" EXCLUDE)
INSTALL (FILES     ${CMAKE_CURRENT_SOURCE_DIR}/README.md DESTINATION .)
INSTALL (FILES     ${CMAKE_CURRENT_SOURCE_DIR}/LICENSE   DESTINATION .)

# We include the "bin" tag into binary archive file name
set (CPACK_PACKAGE_FILE_NAME  ${PROJECT_NAME}-${CPACK_PACKAGE_VERSION}-bin-${CMAKE_SYSTEM_NAME})

# To be done at the end.
INCLUDE (CPack)
