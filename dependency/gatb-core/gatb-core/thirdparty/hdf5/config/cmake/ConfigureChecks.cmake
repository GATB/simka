#-----------------------------------------------------------------------------
# Include all the necessary files for macros
#-----------------------------------------------------------------------------
INCLUDE (${CMAKE_ROOT}/Modules/CheckFunctionExists.cmake)
INCLUDE (${CMAKE_ROOT}/Modules/CheckIncludeFile.cmake)
INCLUDE (${CMAKE_ROOT}/Modules/CheckIncludeFileCXX.cmake)
INCLUDE (${CMAKE_ROOT}/Modules/CheckIncludeFiles.cmake)
INCLUDE (${CMAKE_ROOT}/Modules/CheckLibraryExists.cmake)
INCLUDE (${CMAKE_ROOT}/Modules/CheckSymbolExists.cmake)
INCLUDE (${CMAKE_ROOT}/Modules/CheckTypeSize.cmake)
INCLUDE (${CMAKE_ROOT}/Modules/CheckVariableExists.cmake)
INCLUDE (${CMAKE_ROOT}/Modules/CheckFortranFunctionExists.cmake)
INCLUDE (${CMAKE_ROOT}/Modules/TestBigEndian.cmake)
INCLUDE (${CMAKE_ROOT}/Modules/TestForSTDNamespace.cmake)

#-----------------------------------------------------------------------------
# APPLE/Darwin setup
#-----------------------------------------------------------------------------
IF (APPLE)
  LIST(LENGTH CMAKE_OSX_ARCHITECTURES ARCH_LENGTH)
  IF(ARCH_LENGTH GREATER 1)
    set (CMAKE_OSX_ARCHITECTURES "" CACHE STRING "" FORCE)
    message(FATAL_ERROR "Building Universal Binaries on OS X is NOT supported by the HDF5 project. This is"
    "due to technical reasons. The best approach would be build each architecture in separate directories"
    "and use the 'lipo' tool to combine them into a single executable or library. The 'CMAKE_OSX_ARCHITECTURES'"
    "variable has been set to a blank value which will build the default architecture for this system.")
  ENDIF()
  SET (H5_AC_APPLE_UNIVERSAL_BUILD 0)
ENDIF (APPLE)

# Check for Darwin (not just Apple - we also want to catch OpenDarwin)
IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin") 
    SET (H5_HAVE_DARWIN 1) 
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

#-----------------------------------------------------------------------------
# Option to Clear File Buffers before write --enable-clear-file-buffers
#-----------------------------------------------------------------------------
OPTION (HDF5_Enable_Clear_File_Buffers "Securely clear file buffers before writing to file" ON)
IF (HDF5_Enable_Clear_File_Buffers)
  SET (H5_CLEAR_MEMORY 1)
ENDIF (HDF5_Enable_Clear_File_Buffers)
MARK_AS_ADVANCED (HDF5_Enable_Clear_File_Buffers)

#-----------------------------------------------------------------------------
# Option for --enable-strict-format-checks
#-----------------------------------------------------------------------------
OPTION (HDF5_STRICT_FORMAT_CHECKS "Whether to perform strict file format checks" OFF)
IF (HDF5_STRICT_FORMAT_CHECKS)
  SET (H5_STRICT_FORMAT_CHECKS 1)
ENDIF (HDF5_STRICT_FORMAT_CHECKS)
MARK_AS_ADVANCED (HDF5_STRICT_FORMAT_CHECKS)

#-----------------------------------------------------------------------------
# Option for --enable-metadata-trace-file
#-----------------------------------------------------------------------------
OPTION (HDF5_METADATA_TRACE_FILE "Enable metadata trace file collection" OFF)
IF (HDF5_METADATA_TRACE_FILE)
  SET (H5_METADATA_TRACE_FILE 1)
ENDIF (HDF5_METADATA_TRACE_FILE)
MARK_AS_ADVANCED (HDF5_METADATA_TRACE_FILE)

# ----------------------------------------------------------------------
# Decide whether the data accuracy has higher priority during data
# conversions.  If not, some hard conversions will still be prefered even
# though the data may be wrong (for example, some compilers don't
# support denormalized floating values) to maximize speed.
#
OPTION (HDF5_WANT_DATA_ACCURACY "IF data accuracy is guaranteed during data conversions" ON)
IF (HDF5_WANT_DATA_ACCURACY)
  SET (H5_WANT_DATA_ACCURACY 1)
ENDIF(HDF5_WANT_DATA_ACCURACY)
MARK_AS_ADVANCED (HDF5_WANT_DATA_ACCURACY)

# ----------------------------------------------------------------------
# Decide whether the presence of user's exception handling functions is
# checked and data conversion exceptions are returned.  This is mainly
# for the speed optimization of hard conversions.  Soft conversions can
# actually benefit little.
#
OPTION (HDF5_WANT_DCONV_EXCEPTION "exception handling functions is checked during data conversions" ON)
IF (HDF5_WANT_DCONV_EXCEPTION)
  SET (H5_WANT_DCONV_EXCEPTION 1)
ENDIF (HDF5_WANT_DCONV_EXCEPTION)
MARK_AS_ADVANCED (HDF5_WANT_DCONV_EXCEPTION)

# ----------------------------------------------------------------------
# Check if they would like the function stack support compiled in
#
OPTION (HDF5_ENABLE_CODESTACK "Enable the function stack tracing (for developer debugging)." OFF)
IF (HDF5_ENABLE_CODESTACK)
  SET (H5_HAVE_CODESTACK 1)
ENDIF (HDF5_ENABLE_CODESTACK)
MARK_AS_ADVANCED (HDF5_ENABLE_CODESTACK)

OPTION (HDF5_ENABLE_HSIZET "Enable datasets larger than memory" ON)

# ----------------------------------------------------------------------
# Set the flag to indicate that the machine can handle converting
# floating-point to long long values.
# (This flag should be _unset_ for all machines)
#
#  SET (H5_HW_FP_TO_LLONG_NOT_WORKS 0)

# so far we have no check for this
SET(H5_HAVE_TMPFILE 1)

#-----------------------------------------------------------------------------
# This MACRO checks IF the symbol exists in the library and IF it
# does, it appends library to the list.
#-----------------------------------------------------------------------------
SET (LINK_LIBS "")
MACRO (CHECK_LIBRARY_EXISTS_CONCAT LIBRARY SYMBOL VARIABLE)
  CHECK_LIBRARY_EXISTS ("${LIBRARY};${LINK_LIBS}" ${SYMBOL} "" ${VARIABLE})
  IF (${VARIABLE})
    SET (LINK_LIBS ${LINK_LIBS} ${LIBRARY})
  ENDIF (${VARIABLE})
ENDMACRO (CHECK_LIBRARY_EXISTS_CONCAT)

# ----------------------------------------------------------------------
# WINDOWS Hard code Values
# ----------------------------------------------------------------------

SET (WINDOWS)
IF (WIN32)
  IF (MINGW)
    SET (H5_HAVE_MINGW 1)
    SET (WINDOWS 1) # MinGW tries to imitate Windows
    SET (CMAKE_REQUIRED_FLAGS "-DWIN32_LEAN_AND_MEAN=1 -DNOGDI=1")
  ENDIF (MINGW)
  SET (H5_HAVE_WIN32_API 1)
  SET (CMAKE_REQUIRED_LIBRARIES "ws2_32.lib;wsock32.lib")
  IF (NOT UNIX AND NOT CYGWIN AND NOT MINGW)
    SET (WINDOWS 1)
    SET (CMAKE_REQUIRED_FLAGS "/DWIN32_LEAN_AND_MEAN=1 /DNOGDI=1")
    IF (MSVC)
      SET (H5_HAVE_VISUAL_STUDIO 1)
    ENDIF (MSVC)
  ENDIF (NOT UNIX AND NOT CYGWIN AND NOT MINGW)
ENDIF (WIN32)

#IF (WIN32)
#  SET (DEFAULT_STREAM_VFD OFF)
#ELSE (WIN32)
#  SET (DEFAULT_STREAM_VFD ON)
#ENDIF (WIN32)
#OPTION (HDF5_STREAM_VFD "Compile Stream Virtual File Driver support" ${DEFAULT_STREAM_VFD})

# TODO --------------------------------------------------------------------------
# Should the Default Virtual File Driver be compiled?
# This is hard-coded now but option should added to match configure
#
SET (H5_DEFAULT_VFD H5FD_SEC2)

IF (WINDOWS)
  SET (H5_HAVE_WINDOWS 1)
  # ----------------------------------------------------------------------
  # Set the flag to indicate that the machine has window style pathname,
  # that is, "drive-letter:\" (e.g. "C:") or "drive-letter:/" (e.g. "C:/").
  # (This flag should be _unset_ for all machines, except for Windows)
  SET (H5_HAVE_WINDOW_PATH 1)
ENDIF (WINDOWS)

IF (WINDOWS)
  SET (H5_HAVE_STDDEF_H 1)
  SET (H5_HAVE_SYS_STAT_H 1)
  SET (H5_HAVE_SYS_TYPES_H 1)
  SET (H5_HAVE_LIBM 1)
  SET (H5_HAVE_STRDUP 1)
  SET (H5_HAVE_SYSTEM 1)
  SET (H5_HAVE_LONGJMP 1)
  IF (NOT MINGW)
    SET (H5_HAVE_GETHOSTNAME 1)
  ENDIF (NOT MINGW)
  SET (H5_HAVE_GETCONSOLESCREENBUFFERINFO 1)
  SET (H5_HAVE_FUNCTION 1)
  SET (H5_GETTIMEOFDAY_GIVES_TZ 1)
  SET (H5_HAVE_TIMEZONE 1)
  SET (H5_HAVE_GETTIMEOFDAY 1)
  SET (H5_LONE_COLON 0)
  IF (MINGW)
    SET (H5_HAVE_WINSOCK2_H 1)
  ENDIF (MINGW)
  SET (H5_HAVE_LIBWS2_32 1)
  SET (H5_HAVE_LIBWSOCK32 1)
ENDIF (WINDOWS)

#-----------------------------------------------------------------------------
# These tests need to be manually SET for windows since there is currently
# something not quite correct with the actual test implementation. This affects
# the 'dt_arith' test and most likely lots of other code
# ----------------------------------------------------------------------------
SET (H5_FP_TO_ULLONG_RIGHT_MAXIMUM "" CACHE INTERNAL "")

# ----------------------------------------------------------------------
# END of WINDOWS Hard code Values
# ----------------------------------------------------------------------

IF (CYGWIN)
  SET (H5_HAVE_LSEEK64 0)
ENDIF (CYGWIN)

#-----------------------------------------------------------------------------
#  Check for the math library "m"
#-----------------------------------------------------------------------------
IF (NOT WINDOWS)
  CHECK_LIBRARY_EXISTS_CONCAT ("m" ceil     H5_HAVE_LIBM)
  CHECK_LIBRARY_EXISTS_CONCAT ("dl" dlopen     H5_HAVE_LIBDL)
  CHECK_LIBRARY_EXISTS_CONCAT ("ws2_32" WSAStartup  H5_HAVE_LIBWS2_32)
  CHECK_LIBRARY_EXISTS_CONCAT ("wsock32" gethostbyname H5_HAVE_LIBWSOCK32)
ENDIF (NOT WINDOWS)

CHECK_LIBRARY_EXISTS_CONCAT ("ucb"    gethostname  H5_HAVE_LIBUCB)
CHECK_LIBRARY_EXISTS_CONCAT ("socket" connect      H5_HAVE_LIBSOCKET)
CHECK_LIBRARY_EXISTS ("c" gethostbyname "" NOT_NEED_LIBNSL)

IF (NOT NOT_NEED_LIBNSL)
  CHECK_LIBRARY_EXISTS_CONCAT ("nsl"    gethostbyname  H5_HAVE_LIBNSL)
ENDIF (NOT NOT_NEED_LIBNSL)

# For other tests to use the same libraries
SET (CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES} ${LINK_LIBS})

SET (USE_INCLUDES "")
IF (WINDOWS)
  SET (USE_INCLUDES ${USE_INCLUDES} "windows.h")
ENDIF (WINDOWS)

IF (NOT WINDOWS)
  TEST_BIG_ENDIAN(H5_WORDS_BIGENDIAN)
ENDIF (NOT WINDOWS)

# For other specific tests, use this MACRO.
MACRO (HDF5_FUNCTION_TEST OTHER_TEST)
  IF ("H5_${OTHER_TEST}" MATCHES "^H5_${OTHER_TEST}$")
    SET (MACRO_CHECK_FUNCTION_DEFINITIONS "-D${OTHER_TEST} ${CMAKE_REQUIRED_FLAGS}")
    SET (OTHER_TEST_ADD_LIBRARIES)
    IF (CMAKE_REQUIRED_LIBRARIES)
      SET (OTHER_TEST_ADD_LIBRARIES "-DLINK_LIBRARIES:STRING=${CMAKE_REQUIRED_LIBRARIES}")
    ENDIF (CMAKE_REQUIRED_LIBRARIES)

    FOREACH (def ${HDF5_EXTRA_TEST_DEFINITIONS})
      SET (MACRO_CHECK_FUNCTION_DEFINITIONS "${MACRO_CHECK_FUNCTION_DEFINITIONS} -D${def}=${${def}}")
    ENDFOREACH (def)

    FOREACH (def
        HAVE_SYS_TIME_H
        HAVE_UNISTD_H
        HAVE_SYS_TYPES_H
        HAVE_SYS_SOCKET_H
    )
      IF ("${H5_${def}}")
        SET (MACRO_CHECK_FUNCTION_DEFINITIONS "${MACRO_CHECK_FUNCTION_DEFINITIONS} -D${def}")
      ENDIF ("${H5_${def}}")
    ENDFOREACH (def)

    IF (LARGEFILE)
      SET (MACRO_CHECK_FUNCTION_DEFINITIONS
          "${MACRO_CHECK_FUNCTION_DEFINITIONS} -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -D_LARGEFILE_SOURCE"
      )
    ENDIF (LARGEFILE)

    #MESSAGE (STATUS "Performing ${OTHER_TEST}")
    TRY_COMPILE (${OTHER_TEST}
        ${CMAKE_BINARY_DIR}
        ${HDF5_RESOURCES_DIR}/HDF5Tests.c
        CMAKE_FLAGS -DCOMPILE_DEFINITIONS:STRING=${MACRO_CHECK_FUNCTION_DEFINITIONS}
        "${OTHER_TEST_ADD_LIBRARIES}"
        OUTPUT_VARIABLE OUTPUT
    )
    IF (${OTHER_TEST})
      SET (H5_${OTHER_TEST} 1 CACHE INTERNAL "Other test ${FUNCTION}")
      MESSAGE (STATUS "Performing Other Test ${OTHER_TEST} - Success")
    ELSE (${OTHER_TEST})
      MESSAGE (STATUS "Performing Other Test ${OTHER_TEST} - Failed")
      SET (H5_${OTHER_TEST} "" CACHE INTERNAL "Other test ${FUNCTION}")
      FILE (APPEND ${CMAKE_BINARY_DIR}/CMakeFiles/CMakeError.log
          "Performing Other Test ${OTHER_TEST} failed with the following output:\n"
          "${OUTPUT}\n"
      )
    ENDIF (${OTHER_TEST})
  ENDIF ("H5_${OTHER_TEST}" MATCHES "^H5_${OTHER_TEST}$")
ENDMACRO (HDF5_FUNCTION_TEST)

#-----------------------------------------------------------------------------
# Check for these functions before the time headers are checked
#-----------------------------------------------------------------------------
HDF5_FUNCTION_TEST (STDC_HEADERS)

CHECK_FUNCTION_EXISTS (difftime          H5_HAVE_DIFFTIME)
#CHECK_FUNCTION_EXISTS (gettimeofday      H5_HAVE_GETTIMEOFDAY)
#  Since gettimeofday is not defined any where standard, lets look in all the
#  usual places. On MSVC we are just going to use ::clock()
IF (NOT MSVC)
  IF ("H5_HAVE_TIME_GETTIMEOFDAY" MATCHES "^H5_HAVE_TIME_GETTIMEOFDAY$")
    TRY_COMPILE (HAVE_TIME_GETTIMEOFDAY
        ${CMAKE_BINARY_DIR}
        ${HDF5_RESOURCES_DIR}/GetTimeOfDayTest.cpp
        COMPILE_DEFINITIONS -DTRY_TIME_H
        OUTPUT_VARIABLE OUTPUT
    )
    IF (HAVE_TIME_GETTIMEOFDAY STREQUAL "TRUE")
      SET (H5_HAVE_TIME_GETTIMEOFDAY "1" CACHE INTERNAL "H5_HAVE_TIME_GETTIMEOFDAY")
      SET (H5_HAVE_GETTIMEOFDAY "1" CACHE INTERNAL "H5_HAVE_GETTIMEOFDAY")
    ENDIF (HAVE_TIME_GETTIMEOFDAY STREQUAL "TRUE")
  ENDIF ("H5_HAVE_TIME_GETTIMEOFDAY" MATCHES "^H5_HAVE_TIME_GETTIMEOFDAY$")

  IF ("H5_HAVE_SYS_TIME_GETTIMEOFDAY" MATCHES "^H5_HAVE_SYS_TIME_GETTIMEOFDAY$")
    TRY_COMPILE (HAVE_SYS_TIME_GETTIMEOFDAY
        ${CMAKE_BINARY_DIR}
        ${HDF5_RESOURCES_DIR}/GetTimeOfDayTest.cpp
        COMPILE_DEFINITIONS -DTRY_SYS_TIME_H
        OUTPUT_VARIABLE OUTPUT
    )
    IF (HAVE_SYS_TIME_GETTIMEOFDAY STREQUAL "TRUE")
      SET (H5_HAVE_SYS_TIME_GETTIMEOFDAY "1" CACHE INTERNAL "H5_HAVE_SYS_TIME_GETTIMEOFDAY")
      SET (H5_HAVE_GETTIMEOFDAY "1" CACHE INTERNAL "H5_HAVE_GETTIMEOFDAY")
    ENDIF (HAVE_SYS_TIME_GETTIMEOFDAY STREQUAL "TRUE")
  ENDIF ("H5_HAVE_SYS_TIME_GETTIMEOFDAY" MATCHES "^H5_HAVE_SYS_TIME_GETTIMEOFDAY$")

  IF (NOT HAVE_SYS_TIME_GETTIMEOFDAY AND NOT H5_HAVE_GETTIMEOFDAY)
    MESSAGE (STATUS "---------------------------------------------------------------")
    MESSAGE (STATUS "Function 'gettimeofday()' was not found. HDF5 will use its")
    MESSAGE (STATUS "  own implementation.. This can happen on older versions of")
    MESSAGE (STATUS "  MinGW on Windows. Consider upgrading your MinGW installation")
    MESSAGE (STATUS "  to a newer version such as MinGW 3.12")
    MESSAGE (STATUS "---------------------------------------------------------------")
  ENDIF (NOT HAVE_SYS_TIME_GETTIMEOFDAY AND NOT H5_HAVE_GETTIMEOFDAY)
ENDIF (NOT MSVC)

# Find the library containing clock_gettime()
IF (NOT WINDOWS)
  CHECK_FUNCTION_EXISTS(clock_gettime CLOCK_GETTIME_IN_LIBC)
  CHECK_LIBRARY_EXISTS(rt clock_gettime "" CLOCK_GETTIME_IN_LIBRT)
  CHECK_LIBRARY_EXISTS(posix4 clock_gettime "" CLOCK_GETTIME_IN_LIBPOSIX4)
  IF(CLOCK_GETTIME_IN_LIBC)
    SET(H5_HAVE_CLOCK_GETTIME 1)
  ELSEIF(CLOCK_GETTIME_IN_LIBRT)
    SET(H5_HAVE_CLOCK_GETTIME 1)
    LIST(APPEND LINK_LIBS rt)
  ELSEIF(CLOCK_GETTIME_IN_LIBPOSIX4)
    SET(H5_HAVE_CLOCK_GETTIME 1)
    LIST(APPEND LINK_LIBS posix4)
  ENDIF(CLOCK_GETTIME_IN_LIBC)
ENDIF (NOT WINDOWS)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Check IF header file exists and add it to the list.
#-----------------------------------------------------------------------------
MACRO (CHECK_INCLUDE_FILE_CONCAT FILE VARIABLE)
  CHECK_INCLUDE_FILES ("${USE_INCLUDES};${FILE}" ${VARIABLE})
  IF (${VARIABLE})
    SET (USE_INCLUDES ${USE_INCLUDES} ${FILE})
  ENDIF (${VARIABLE})
ENDMACRO (CHECK_INCLUDE_FILE_CONCAT)

#-----------------------------------------------------------------------------
#  Check for the existence of certain header files
#-----------------------------------------------------------------------------
CHECK_INCLUDE_FILE_CONCAT ("sys/resource.h"  H5_HAVE_SYS_RESOURCE_H)
CHECK_INCLUDE_FILE_CONCAT ("sys/time.h"      H5_HAVE_SYS_TIME_H)
CHECK_INCLUDE_FILE_CONCAT ("unistd.h"        H5_HAVE_UNISTD_H)
CHECK_INCLUDE_FILE_CONCAT ("sys/ioctl.h"     H5_HAVE_SYS_IOCTL_H)
CHECK_INCLUDE_FILE_CONCAT ("sys/stat.h"      H5_HAVE_SYS_STAT_H)
CHECK_INCLUDE_FILE_CONCAT ("sys/socket.h"    H5_HAVE_SYS_SOCKET_H)
CHECK_INCLUDE_FILE_CONCAT ("sys/types.h"     H5_HAVE_SYS_TYPES_H)
CHECK_INCLUDE_FILE_CONCAT ("stddef.h"        H5_HAVE_STDDEF_H)
CHECK_INCLUDE_FILE_CONCAT ("setjmp.h"        H5_HAVE_SETJMP_H)
CHECK_INCLUDE_FILE_CONCAT ("features.h"      H5_HAVE_FEATURES_H)
CHECK_INCLUDE_FILE_CONCAT ("dirent.h"        H5_HAVE_DIRENT_H)
CHECK_INCLUDE_FILE_CONCAT ("stdint.h"        H5_HAVE_STDINT_H)

# IF the c compiler found stdint, check the C++ as well. On some systems this
# file will be found by C but not C++, only do this test IF the C++ compiler
# has been initialized (e.g. the project also includes some c++)
IF (H5_HAVE_STDINT_H AND CMAKE_CXX_COMPILER_LOADED)
  CHECK_INCLUDE_FILE_CXX ("stdint.h" H5_HAVE_STDINT_H_CXX)
  IF (NOT H5_HAVE_STDINT_H_CXX)
    SET (H5_HAVE_STDINT_H "" CACHE INTERNAL "Have includes HAVE_STDINT_H")
    SET (USE_INCLUDES ${USE_INCLUDES} "stdint.h")
  ENDIF (NOT H5_HAVE_STDINT_H_CXX)
ENDIF (H5_HAVE_STDINT_H AND CMAKE_CXX_COMPILER_LOADED)

# Darwin
CHECK_INCLUDE_FILE_CONCAT ("mach/mach_time.h" H5_HAVE_MACH_MACH_TIME_H)

# Windows
CHECK_INCLUDE_FILE_CONCAT ("io.h"            H5_HAVE_IO_H)
IF (NOT CYGWIN)
  CHECK_INCLUDE_FILE_CONCAT ("winsock2.h"      H5_HAVE_WINSOCK2_H)
ENDIF (NOT CYGWIN)
CHECK_INCLUDE_FILE_CONCAT ("sys/timeb.h"     H5_HAVE_SYS_TIMEB_H)

IF (CMAKE_SYSTEM_NAME MATCHES "OSF")
  CHECK_INCLUDE_FILE_CONCAT ("sys/sysinfo.h" H5_HAVE_SYS_SYSINFO_H)
  CHECK_INCLUDE_FILE_CONCAT ("sys/proc.h"    H5_HAVE_SYS_PROC_H)
ELSE (CMAKE_SYSTEM_NAME MATCHES "OSF")
  SET (H5_HAVE_SYS_SYSINFO_H "" CACHE INTERNAL "" FORCE)
  SET (H5_HAVE_SYS_PROC_H    "" CACHE INTERNAL "" FORCE)
ENDIF (CMAKE_SYSTEM_NAME MATCHES "OSF")

CHECK_INCLUDE_FILE_CONCAT ("globus/common.h" H5_HAVE_GLOBUS_COMMON_H)
CHECK_INCLUDE_FILE_CONCAT ("pdb.h"           H5_HAVE_PDB_H)
CHECK_INCLUDE_FILE_CONCAT ("pthread.h"       H5_HAVE_PTHREAD_H)
CHECK_INCLUDE_FILE_CONCAT ("srbclient.h"     H5_HAVE_SRBCLIENT_H)
CHECK_INCLUDE_FILE_CONCAT ("string.h"        H5_HAVE_STRING_H)
CHECK_INCLUDE_FILE_CONCAT ("strings.h"       H5_HAVE_STRINGS_H)
CHECK_INCLUDE_FILE_CONCAT ("time.h"          H5_HAVE_TIME_H)
CHECK_INCLUDE_FILE_CONCAT ("stdlib.h"        H5_HAVE_STDLIB_H)
CHECK_INCLUDE_FILE_CONCAT ("memory.h"        H5_HAVE_MEMORY_H)
CHECK_INCLUDE_FILE_CONCAT ("dlfcn.h"         H5_HAVE_DLFCN_H)
CHECK_INCLUDE_FILE_CONCAT ("inttypes.h"      H5_HAVE_INTTYPES_H)
CHECK_INCLUDE_FILE_CONCAT ("netinet/in.h"    H5_HAVE_NETINET_IN_H)

#-----------------------------------------------------------------------------
#  Check for large file support
#-----------------------------------------------------------------------------

# The linux-lfs option is deprecated.
SET (LINUX_LFS 0)

SET (HDF5_EXTRA_C_FLAGS)
SET (HDF5_EXTRA_FLAGS)
IF (NOT WINDOWS)
  # Linux Specific flags
  # This was originally defined as _POSIX_SOURCE which was updated to
  # _POSIX_C_SOURCE=199506L to expose a greater amount of POSIX
  # functionality so clock_gettime and CLOCK_MONOTONIC are defined
  # correctly.
  # POSIX feature information can be found in the gcc manual at:
  # http://www.gnu.org/s/libc/manual/html_node/Feature-Test-Macros.html
  SET (HDF5_EXTRA_C_FLAGS -D_POSIX_C_SOURCE=199506L)
  SET (HDF5_EXTRA_FLAGS -D_BSD_SOURCE)
  
  OPTION (HDF5_ENABLE_LARGE_FILE "Enable support for large (64-bit) files on Linux." ON)
  IF (HDF5_ENABLE_LARGE_FILE)
    SET (msg "Performing TEST_LFS_WORKS")
    TRY_RUN (TEST_LFS_WORKS_RUN   TEST_LFS_WORKS_COMPILE
        ${HDF5_BINARY_DIR}/CMake
        ${HDF5_RESOURCES_DIR}/HDF5Tests.c
        CMAKE_FLAGS -DCOMPILE_DEFINITIONS:STRING=-DTEST_LFS_WORKS
        OUTPUT_VARIABLE OUTPUT
    )
    IF (TEST_LFS_WORKS_COMPILE)
      IF (TEST_LFS_WORKS_RUN  MATCHES 0)
        SET (TEST_LFS_WORKS 1 CACHE INTERNAL ${msg})
        SET (LARGEFILE 1)
        SET (HDF5_EXTRA_FLAGS ${HDF5_EXTRA_FLAGS} -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -D_LARGEFILE_SOURCE)
        MESSAGE (STATUS "${msg}... yes")
      ELSE (TEST_LFS_WORKS_RUN  MATCHES 0)
        SET (TEST_LFS_WORKS "" CACHE INTERNAL ${msg})
        MESSAGE (STATUS "${msg}... no")
        FILE (APPEND ${CMAKE_BINARY_DIR}/CMakeFiles/CMakeError.log
              "Test TEST_LFS_WORKS Run failed with the following output and exit code:\n ${OUTPUT}\n"
        )
      ENDIF (TEST_LFS_WORKS_RUN  MATCHES 0)
    ELSE (TEST_LFS_WORKS_COMPILE )
      SET (TEST_LFS_WORKS "" CACHE INTERNAL ${msg})
      MESSAGE (STATUS "${msg}... no")
      FILE (APPEND ${CMAKE_BINARY_DIR}/CMakeFiles/CMakeError.log
          "Test TEST_LFS_WORKS Compile failed with the following output:\n ${OUTPUT}\n"
      )
    ENDIF (TEST_LFS_WORKS_COMPILE)
  ENDIF (HDF5_ENABLE_LARGE_FILE)
  SET (CMAKE_REQUIRED_DEFINITIONS ${CMAKE_REQUIRED_DEFINITIONS} ${HDF5_EXTRA_FLAGS})
ENDIF (NOT WINDOWS)

ADD_DEFINITIONS (${HDF5_EXTRA_FLAGS})

#-----------------------------------------------------------------------------
# Check for HAVE_OFF64_T functionality
#-----------------------------------------------------------------------------
IF (NOT WINDOWS OR MINGW)
  HDF5_FUNCTION_TEST (HAVE_OFF64_T)
  IF (H5_HAVE_OFF64_T)
    CHECK_FUNCTION_EXISTS (lseek64            H5_HAVE_LSEEK64)
    CHECK_FUNCTION_EXISTS (fseeko64           H5_HAVE_FSEEKO64)
    CHECK_FUNCTION_EXISTS (ftello64           H5_HAVE_FTELLO64)
    CHECK_FUNCTION_EXISTS (ftruncate64        H5_HAVE_FTRUNCATE64)
  ENDIF (H5_HAVE_OFF64_T)

  CHECK_FUNCTION_EXISTS (fseeko               H5_HAVE_FSEEKO)
  CHECK_FUNCTION_EXISTS (ftello               H5_HAVE_FTELLO)

  HDF5_FUNCTION_TEST (HAVE_STAT64_STRUCT)
  IF (HAVE_STAT64_STRUCT)
    CHECK_FUNCTION_EXISTS (fstat64            H5_HAVE_FSTAT64)
    CHECK_FUNCTION_EXISTS (stat64             H5_HAVE_STAT64)
  ENDIF (HAVE_STAT64_STRUCT)
ENDIF (NOT WINDOWS OR MINGW)

#-----------------------------------------------------------------------------
#  Check the size in bytes of all the int and float types
#-----------------------------------------------------------------------------
MACRO (H5_CHECK_TYPE_SIZE type var)
  SET (aType ${type})
  SET (aVar  ${var})
#  MESSAGE (STATUS "Checking size of ${aType} and storing into ${aVar}")
  CHECK_TYPE_SIZE (${aType}   ${aVar})
  IF (NOT ${aVar})
    SET (${aVar} 0 CACHE INTERNAL "SizeOf for ${aType}")
#    MESSAGE (STATUS "Size of ${aType} was NOT Found")
  ENDIF (NOT ${aVar})
ENDMACRO (H5_CHECK_TYPE_SIZE)

H5_CHECK_TYPE_SIZE (char           H5_SIZEOF_CHAR)
H5_CHECK_TYPE_SIZE (short          H5_SIZEOF_SHORT)
H5_CHECK_TYPE_SIZE (int            H5_SIZEOF_INT)
H5_CHECK_TYPE_SIZE (unsigned       H5_SIZEOF_UNSIGNED)
IF (NOT APPLE)
  H5_CHECK_TYPE_SIZE (long         H5_SIZEOF_LONG)
ENDIF (NOT APPLE)
H5_CHECK_TYPE_SIZE ("long long"    H5_SIZEOF_LONG_LONG)
H5_CHECK_TYPE_SIZE (__int64        H5_SIZEOF___INT64)
IF (NOT H5_SIZEOF___INT64)
  SET (H5_SIZEOF___INT64 0)
ENDIF (NOT H5_SIZEOF___INT64)

H5_CHECK_TYPE_SIZE (float          H5_SIZEOF_FLOAT)
H5_CHECK_TYPE_SIZE (double         H5_SIZEOF_DOUBLE)
H5_CHECK_TYPE_SIZE ("long double"  H5_SIZEOF_LONG_DOUBLE)

H5_CHECK_TYPE_SIZE (int8_t         H5_SIZEOF_INT8_T)
H5_CHECK_TYPE_SIZE (uint8_t        H5_SIZEOF_UINT8_T)
H5_CHECK_TYPE_SIZE (int_least8_t   H5_SIZEOF_INT_LEAST8_T)
H5_CHECK_TYPE_SIZE (uint_least8_t  H5_SIZEOF_UINT_LEAST8_T)
H5_CHECK_TYPE_SIZE (int_fast8_t    H5_SIZEOF_INT_FAST8_T)
H5_CHECK_TYPE_SIZE (uint_fast8_t   H5_SIZEOF_UINT_FAST8_T)

H5_CHECK_TYPE_SIZE (int16_t        H5_SIZEOF_INT16_T)
H5_CHECK_TYPE_SIZE (uint16_t       H5_SIZEOF_UINT16_T)
H5_CHECK_TYPE_SIZE (int_least16_t  H5_SIZEOF_INT_LEAST16_T)
H5_CHECK_TYPE_SIZE (uint_least16_t H5_SIZEOF_UINT_LEAST16_T)
H5_CHECK_TYPE_SIZE (int_fast16_t   H5_SIZEOF_INT_FAST16_T)
H5_CHECK_TYPE_SIZE (uint_fast16_t  H5_SIZEOF_UINT_FAST16_T)

H5_CHECK_TYPE_SIZE (int32_t        H5_SIZEOF_INT32_T)
H5_CHECK_TYPE_SIZE (uint32_t       H5_SIZEOF_UINT32_T)
H5_CHECK_TYPE_SIZE (int_least32_t  H5_SIZEOF_INT_LEAST32_T)
H5_CHECK_TYPE_SIZE (uint_least32_t H5_SIZEOF_UINT_LEAST32_T)
H5_CHECK_TYPE_SIZE (int_fast32_t   H5_SIZEOF_INT_FAST32_T)
H5_CHECK_TYPE_SIZE (uint_fast32_t  H5_SIZEOF_UINT_FAST32_T)

H5_CHECK_TYPE_SIZE (int64_t        H5_SIZEOF_INT64_T)
H5_CHECK_TYPE_SIZE (uint64_t       H5_SIZEOF_UINT64_T)
H5_CHECK_TYPE_SIZE (int_least64_t  H5_SIZEOF_INT_LEAST64_T)
H5_CHECK_TYPE_SIZE (uint_least64_t H5_SIZEOF_UINT_LEAST64_T)
H5_CHECK_TYPE_SIZE (int_fast64_t   H5_SIZEOF_INT_FAST64_T)
H5_CHECK_TYPE_SIZE (uint_fast64_t  H5_SIZEOF_UINT_FAST64_T)

IF (NOT APPLE)
  H5_CHECK_TYPE_SIZE (size_t       H5_SIZEOF_SIZE_T)
  H5_CHECK_TYPE_SIZE (ssize_t      H5_SIZEOF_SSIZE_T)
  IF (NOT H5_SIZEOF_SSIZE_T)
    SET (H5_SIZEOF_SSIZE_T 0)
  ENDIF (NOT H5_SIZEOF_SSIZE_T)
  H5_CHECK_TYPE_SIZE (ptrdiff_t    H5_SIZEOF_PTRDIFF_T)
ENDIF (NOT APPLE)

H5_CHECK_TYPE_SIZE (off_t          H5_SIZEOF_OFF_T)
H5_CHECK_TYPE_SIZE (off64_t        H5_SIZEOF_OFF64_T)
IF (NOT H5_SIZEOF_OFF64_T)
  SET (H5_SIZEOF_OFF64_T 0)
ENDIF (NOT H5_SIZEOF_OFF64_T)

IF (NOT WINDOWS)
  #-----------------------------------------------------------------------------
  # Check if the dev_t type is a scalar type
  #-----------------------------------------------------------------------------
  HDF5_FUNCTION_TEST (DEV_T_IS_SCALAR)

  # ----------------------------------------------------------------------
  # Check for MONOTONIC_TIMER support (used in clock_gettime).  This has
  # to be done after any POSIX/BSD defines to ensure that the test gets
  # the correct POSIX level on linux.
  CHECK_VARIABLE_EXISTS (CLOCK_MONOTONIC HAVE_CLOCK_MONOTONIC)

  #-----------------------------------------------------------------------------
  # Check a bunch of time functions
  #-----------------------------------------------------------------------------
  FOREACH (test
      HAVE_TM_GMTOFF
      HAVE___TM_GMTOFF
#      HAVE_TIMEZONE
      HAVE_STRUCT_TIMEZONE
      GETTIMEOFDAY_GIVES_TZ
      TIME_WITH_SYS_TIME
      HAVE_TM_ZONE
      HAVE_STRUCT_TM_TM_ZONE
  )
    HDF5_FUNCTION_TEST (${test})
  ENDFOREACH (test)
  IF (NOT CYGWIN AND NOT MINGW)
      HDF5_FUNCTION_TEST (HAVE_TIMEZONE)
#      HDF5_FUNCTION_TEST (HAVE_STAT_ST_BLOCKS)
  ENDIF (NOT CYGWIN AND NOT MINGW)

  # ----------------------------------------------------------------------
  # Does the struct stat have the st_blocks field?  This field is not Posix.
  #
  HDF5_FUNCTION_TEST (HAVE_STAT_ST_BLOCKS)

  # ----------------------------------------------------------------------
  # How do we figure out the width of a tty in characters?
  #
  CHECK_FUNCTION_EXISTS (ioctl             H5_HAVE_IOCTL)
  HDF5_FUNCTION_TEST (HAVE_STRUCT_VIDEOCONFIG)
  HDF5_FUNCTION_TEST (HAVE_STRUCT_TEXT_INFO)
  CHECK_FUNCTION_EXISTS (_getvideoconfig   H5_HAVE__GETVIDEOCONFIG)
  CHECK_FUNCTION_EXISTS (gettextinfo       H5_HAVE_GETTEXTINFO)
  CHECK_FUNCTION_EXISTS (_scrsize          H5_HAVE__SCRSIZE)
  CHECK_FUNCTION_EXISTS (GetConsoleScreenBufferInfo    H5_HAVE_GETCONSOLESCREENBUFFERINFO)
  CHECK_SYMBOL_EXISTS (TIOCGWINSZ "sys/ioctl.h" H5_HAVE_TIOCGWINSZ)
  CHECK_SYMBOL_EXISTS (TIOCGETD   "sys/ioctl.h" H5_HAVE_TIOCGETD)
ENDIF (NOT WINDOWS)

#-----------------------------------------------------------------------------
# Check for some functions that are used
#
CHECK_FUNCTION_EXISTS (alarm             H5_HAVE_ALARM)
#CHECK_FUNCTION_EXISTS (BSDgettimeofday   H5_HAVE_BSDGETTIMEOFDAY)
CHECK_FUNCTION_EXISTS (fork              H5_HAVE_FORK)
CHECK_FUNCTION_EXISTS (frexpf            H5_HAVE_FREXPF)
CHECK_FUNCTION_EXISTS (frexpl            H5_HAVE_FREXPL)

CHECK_FUNCTION_EXISTS (gethostname       H5_HAVE_GETHOSTNAME)
CHECK_FUNCTION_EXISTS (getpwuid          H5_HAVE_GETPWUID)
CHECK_FUNCTION_EXISTS (getrusage         H5_HAVE_GETRUSAGE)
CHECK_FUNCTION_EXISTS (lstat             H5_HAVE_LSTAT)

CHECK_FUNCTION_EXISTS (rand_r            H5_HAVE_RAND_R)
CHECK_FUNCTION_EXISTS (random            H5_HAVE_RANDOM)
CHECK_FUNCTION_EXISTS (setsysinfo        H5_HAVE_SETSYSINFO)

CHECK_FUNCTION_EXISTS (signal            H5_HAVE_SIGNAL)
CHECK_FUNCTION_EXISTS (longjmp           H5_HAVE_LONGJMP)
CHECK_FUNCTION_EXISTS (setjmp            H5_HAVE_SETJMP)
CHECK_FUNCTION_EXISTS (siglongjmp        H5_HAVE_SIGLONGJMP)
CHECK_FUNCTION_EXISTS (sigsetjmp         H5_HAVE_SIGSETJMP)
CHECK_FUNCTION_EXISTS (sigaction         H5_HAVE_SIGACTION)
CHECK_FUNCTION_EXISTS (sigprocmask       H5_HAVE_SIGPROCMASK)

CHECK_FUNCTION_EXISTS (snprintf          H5_HAVE_SNPRINTF)
CHECK_FUNCTION_EXISTS (srandom           H5_HAVE_SRANDOM)
CHECK_FUNCTION_EXISTS (strdup            H5_HAVE_STRDUP)
CHECK_FUNCTION_EXISTS (symlink           H5_HAVE_SYMLINK)
CHECK_FUNCTION_EXISTS (system            H5_HAVE_SYSTEM)

CHECK_FUNCTION_EXISTS (tmpfile           H5_HAVE_TMPFILE)
CHECK_FUNCTION_EXISTS (vasprintf         H5_HAVE_VASPRINTF)
CHECK_FUNCTION_EXISTS (waitpid           H5_HAVE_WAITPID)

CHECK_FUNCTION_EXISTS (vsnprintf         H5_HAVE_VSNPRINTF)
IF (NOT WINDOWS)
  IF (H5_HAVE_VSNPRINTF)
    HDF5_FUNCTION_TEST (VSNPRINTF_WORKS)
  ENDIF (H5_HAVE_VSNPRINTF)
ENDIF (NOT WINDOWS)

#-----------------------------------------------------------------------------
# sigsetjmp is special; may actually be a macro
#-----------------------------------------------------------------------------
IF (NOT H5_HAVE_SIGSETJMP)
  IF (H5_HAVE_SETJMP_H)
    CHECK_SYMBOL_EXISTS (sigsetjmp "setjmp.h" H5_HAVE_MACRO_SIGSETJMP)
    IF (H5_HAVE_MACRO_SIGSETJMP)
      SET (H5_HAVE_SIGSETJMP 1)
    ENDIF (H5_HAVE_MACRO_SIGSETJMP)
  ENDIF (H5_HAVE_SETJMP_H)
ENDIF (NOT H5_HAVE_SIGSETJMP)

#-----------------------------------------------------------------------------
# Check for Symbols
CHECK_SYMBOL_EXISTS (tzname "time.h" H5_HAVE_DECL_TZNAME)

#-----------------------------------------------------------------------------
# Check a bunch of other functions
#-----------------------------------------------------------------------------
IF (NOT WINDOWS)
  FOREACH (test
      LONE_COLON
      HAVE_ATTRIBUTE
      HAVE_C99_FUNC
      HAVE_FUNCTION
      HAVE_C99_DESIGNATED_INITIALIZER
      SYSTEM_SCOPE_THREADS
      HAVE_SOCKLEN_T
      CXX_HAVE_OFFSETOF
  )
    HDF5_FUNCTION_TEST (${test})
  ENDFOREACH (test)
ENDIF (NOT WINDOWS)

# For other CXX specific tests, use this MACRO.
MACRO (HDF5_CXX_FUNCTION_TEST OTHER_TEST)
  IF ("${OTHER_TEST}" MATCHES "^${OTHER_TEST}$")
    SET (MACRO_CHECK_FUNCTION_DEFINITIONS "-D${OTHER_TEST} ${CMAKE_REQUIRED_FLAGS}")
    SET (OTHER_TEST_ADD_LIBRARIES)
    IF (CMAKE_REQUIRED_LIBRARIES)
      SET (OTHER_TEST_ADD_LIBRARIES "-DLINK_LIBRARIES:STRING=${CMAKE_REQUIRED_LIBRARIES}")
    ENDIF (CMAKE_REQUIRED_LIBRARIES)

    FOREACH (def ${HDF5_EXTRA_TEST_DEFINITIONS})
      SET (MACRO_CHECK_FUNCTION_DEFINITIONS "${MACRO_CHECK_FUNCTION_DEFINITIONS} -D${def}=${${def}}")
    ENDFOREACH (def)

    FOREACH (def
        HAVE_SYS_TIME_H
        HAVE_UNISTD_H
        HAVE_SYS_TYPES_H
        HAVE_SYS_SOCKET_H
    )
      IF ("${H5_${def}}")
        SET (MACRO_CHECK_FUNCTION_DEFINITIONS "${MACRO_CHECK_FUNCTION_DEFINITIONS} -D${def}")
      ENDIF ("${H5_${def}}")
    ENDFOREACH (def)

    IF (LARGEFILE)
      SET (MACRO_CHECK_FUNCTION_DEFINITIONS
          "${MACRO_CHECK_FUNCTION_DEFINITIONS} -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -D_LARGEFILE_SOURCE"
      )
    ENDIF (LARGEFILE)

    #MESSAGE (STATUS "Performing ${OTHER_TEST}")
    TRY_COMPILE (${OTHER_TEST}
        ${CMAKE_BINARY_DIR}
        ${HDF5_RESOURCES_DIR}/HDF5CXXTests.cpp
        CMAKE_FLAGS -DCOMPILE_DEFINITIONS:STRING=${MACRO_CHECK_FUNCTION_DEFINITIONS}
        "${OTHER_TEST_ADD_LIBRARIES}"
        OUTPUT_VARIABLE OUTPUT
    )
    IF ("${OTHER_TEST}" EQUAL 0)
      SET (${OTHER_TEST} 1 CACHE INTERNAL "CXX test ${FUNCTION}")
      MESSAGE (STATUS "Performing CXX Test ${OTHER_TEST} - Success")
    ELSE ("${OTHER_TEST}" EQUAL 0)
      MESSAGE (STATUS "Performing CXX Test ${OTHER_TEST} - Failed")
      SET (${OTHER_TEST} "" CACHE INTERNAL "CXX test ${FUNCTION}")
      FILE (APPEND ${CMAKE_BINARY_DIR}/CMakeFiles/CMakeError.log
          "Performing CXX Test ${OTHER_TEST} failed with the following output:\n"
          "${OUTPUT}\n"
      )
    ENDIF ("${OTHER_TEST}" EQUAL 0)
  ENDIF ("${OTHER_TEST}" MATCHES "^${OTHER_TEST}$")
ENDMACRO (HDF5_CXX_FUNCTION_TEST)

#-----------------------------------------------------------------------------
# Check a bunch of cxx functions
#-----------------------------------------------------------------------------
IF (CMAKE_CXX_COMPILER_LOADED)
  FOREACH (test
      OLD_HEADER_FILENAME
      H5_NO_NAMESPACE
      H5_NO_STD
      BOOL_NOTDEFINED
      NO_STATIC_CAST
  )
    HDF5_CXX_FUNCTION_TEST (${test})
  ENDFOREACH (test)
ENDIF (CMAKE_CXX_COMPILER_LOADED)

#-----------------------------------------------------------------------------
#  Check if Direct I/O driver works
#-----------------------------------------------------------------------------
IF (NOT WINDOWS)
  OPTION (HDF5_ENABLE_DIRECT_VFD "Build the Direct I/O Virtual File Driver" ON)
  IF (HDF5_ENABLE_DIRECT_VFD)
    SET (msg "Performing TEST_DIRECT_VFD_WORKS")
    SET (MACRO_CHECK_FUNCTION_DEFINITIONS "-DTEST_DIRECT_VFD_WORKS -D_GNU_SOURCE ${CMAKE_REQUIRED_FLAGS}")
    TRY_RUN (TEST_DIRECT_VFD_WORKS_RUN   TEST_DIRECT_VFD_WORKS_COMPILE
        ${HDF5_BINARY_DIR}/CMake
        ${HDF5_RESOURCES_DIR}/HDF5Tests.c
        CMAKE_FLAGS -DCOMPILE_DEFINITIONS:STRING=${MACRO_CHECK_FUNCTION_DEFINITIONS}
        OUTPUT_VARIABLE OUTPUT
    )
    IF (TEST_DIRECT_VFD_WORKS_COMPILE)
      IF (TEST_DIRECT_VFD_WORKS_RUN  MATCHES 0)
        HDF5_FUNCTION_TEST (HAVE_DIRECT)
        SET (CMAKE_REQUIRED_DEFINITIONS "${CMAKE_REQUIRED_DEFINITIONS} -D_GNU_SOURCE")
        ADD_DEFINITIONS ("-D_GNU_SOURCE")
      ELSE (TEST_DIRECT_VFD_WORKS_RUN  MATCHES 0)
        SET (TEST_DIRECT_VFD_WORKS "" CACHE INTERNAL ${msg})
        MESSAGE (STATUS "${msg}... no")
        FILE (APPEND ${CMAKE_BINARY_DIR}/CMakeFiles/CMakeError.log
              "Test TEST_DIRECT_VFD_WORKS Run failed with the following output and exit code:\n ${OUTPUT}\n"
        )
      ENDIF (TEST_DIRECT_VFD_WORKS_RUN  MATCHES 0)
    ELSE (TEST_DIRECT_VFD_WORKS_COMPILE )
      SET (TEST_DIRECT_VFD_WORKS "" CACHE INTERNAL ${msg})
      MESSAGE (STATUS "${msg}... no")
      FILE (APPEND ${CMAKE_BINARY_DIR}/CMakeFiles/CMakeError.log
          "Test TEST_DIRECT_VFD_WORKS Compile failed with the following output:\n ${OUTPUT}\n"
      )
    ENDIF (TEST_DIRECT_VFD_WORKS_COMPILE)
  ENDIF (HDF5_ENABLE_DIRECT_VFD)
ENDIF (NOT WINDOWS)

#-----------------------------------------------------------------------------
#  Check for the Stream VFD driver
#-----------------------------------------------------------------------------
IF (HDF5_STREAM_VFD)
  CHECK_INCLUDE_FILE_CONCAT ("netdb.h"       H5_HAVE_NETDB_H)
  CHECK_INCLUDE_FILE_CONCAT ("netinet/tcp.h" H5_HAVE_NETINET_TCP_H)
  CHECK_INCLUDE_FILE_CONCAT ("sys/filio.h"   H5_HAVE_SYS_FILIO_H)
  SET (H5_HAVE_STREAM 1)
ENDIF (HDF5_STREAM_VFD)

#-----------------------------------------------------------------------------
# Check if InitOnceExecuteOnce is available
#-----------------------------------------------------------------------------
IF (WINDOWS)
  IF (NOT HDF5_NO_IOEO_TEST)
  MESSAGE (STATUS "Checking for InitOnceExecuteOnce:")
  IF("${H5_HAVE_IOEO}" MATCHES "^${H5_HAVE_IOEO}$")
    IF (LARGEFILE)
      SET (CMAKE_REQUIRED_DEFINITIONS
          "${CURRENT_TEST_DEFINITIONS} -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -D_LARGEFILE_SOURCE"
      )
    ENDIF (LARGEFILE)
    SET(MACRO_CHECK_FUNCTION_DEFINITIONS 
      "-DHAVE_IOEO ${CMAKE_REQUIRED_FLAGS}")
    IF(CMAKE_REQUIRED_LIBRARIES)
      SET(CHECK_C_SOURCE_COMPILES_ADD_LIBRARIES
        "-DLINK_LIBRARIES:STRING=${CMAKE_REQUIRED_LIBRARIES}")
    ELSE(CMAKE_REQUIRED_LIBRARIES)
      SET(CHECK_C_SOURCE_COMPILES_ADD_LIBRARIES)
    ENDIF(CMAKE_REQUIRED_LIBRARIES)
    IF(CMAKE_REQUIRED_INCLUDES)
      SET(CHECK_C_SOURCE_COMPILES_ADD_INCLUDES
        "-DINCLUDE_DIRECTORIES:STRING=${CMAKE_REQUIRED_INCLUDES}")
    ELSE(CMAKE_REQUIRED_INCLUDES)
      SET(CHECK_C_SOURCE_COMPILES_ADD_INCLUDES)
    ENDIF(CMAKE_REQUIRED_INCLUDES)

    TRY_RUN(HAVE_IOEO_EXITCODE HAVE_IOEO_COMPILED
      ${CMAKE_BINARY_DIR}
      ${HDF5_RESOURCES_DIR}/HDF5Tests.c
      COMPILE_DEFINITIONS ${CMAKE_REQUIRED_DEFINITIONS}
      CMAKE_FLAGS -DCOMPILE_DEFINITIONS:STRING=${MACRO_CHECK_FUNCTION_DEFINITIONS}
      -DCMAKE_SKIP_RPATH:BOOL=${CMAKE_SKIP_RPATH}
      "${CHECK_C_SOURCE_COMPILES_ADD_LIBRARIES}"
      "${CHECK_C_SOURCE_COMPILES_ADD_INCLUDES}"
      COMPILE_OUTPUT_VARIABLE OUTPUT)
    # if it did not compile make the return value fail code of 1
    IF(NOT HAVE_IOEO_COMPILED)
      SET(HAVE_IOEO_EXITCODE 1)
    ENDIF(NOT HAVE_IOEO_COMPILED)
    # if the return value was 0 then it worked
    IF("${HAVE_IOEO_EXITCODE}" EQUAL 0)
      SET(H5_HAVE_IOEO 1 CACHE INTERNAL "Test InitOnceExecuteOnce")
      MESSAGE(STATUS "Performing Test InitOnceExecuteOnce - Success")
      FILE(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeOutput.log 
        "Performing C SOURCE FILE Test InitOnceExecuteOnce succeded with the following output:\n"
        "${OUTPUT}\n"
        "Return value: ${HAVE_IOEO}\n")
    ELSE("${HAVE_IOEO_EXITCODE}" EQUAL 0)
      IF(CMAKE_CROSSCOMPILING AND "${HAVE_IOEO_EXITCODE}" MATCHES  "FAILED_TO_RUN")
        SET(H5_HAVE_IOEO "${HAVE_IOEO_EXITCODE}")
      ELSE(CMAKE_CROSSCOMPILING AND "${HAVE_IOEO_EXITCODE}" MATCHES  "FAILED_TO_RUN")
        SET(H5_HAVE_IOEO "" CACHE INTERNAL "Test InitOnceExecuteOnce")
      ENDIF(CMAKE_CROSSCOMPILING AND "${HAVE_IOEO_EXITCODE}" MATCHES  "FAILED_TO_RUN")

      MESSAGE(STATUS "Performing Test InitOnceExecuteOnce - Failed")
      FILE(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log 
        "Performing InitOnceExecuteOnce Test  failed with the following output:\n"
        "${OUTPUT}\n"
        "Return value: ${HAVE_IOEO_EXITCODE}\n")
    ENDIF("${HAVE_IOEO_EXITCODE}" EQUAL 0)
  ENDIF("${H5_HAVE_IOEO}" MATCHES "^${H5_HAVE_IOEO}$")
  ENDIF (NOT HDF5_NO_IOEO_TEST)
ENDIF (WINDOWS)

#-----------------------------------------------------------------------------
# Option to see if GPFS is available on this filesystem --enable-gpfs
#-----------------------------------------------------------------------------
OPTION (HDF5_ENABLE_GPFS "Enable GPFS hints for the MPI/POSIX file driver" OFF)
IF (HDF5_ENABLE_GPFS)
  CHECK_INCLUDE_FILE_CONCAT ("gpfs.h"        HAVE_GPFS)
  IF (HAVE_GPFS)
    HDF5_FUNCTION_TEST (HAVE_GPFS)  
  ENDIF (HAVE_GPFS)
ENDIF (HDF5_ENABLE_GPFS)
MARK_AS_ADVANCED (HDF5_ENABLE_GPFS)

#-----------------------------------------------------------------------------
# Determine how 'inline' is used
#-----------------------------------------------------------------------------
SET (HDF5_EXTRA_TEST_DEFINITIONS INLINE_TEST_INLINE)
FOREACH (inline_test inline __inline__ __inline)
  SET (INLINE_TEST_INLINE ${inline_test})
  HDF5_FUNCTION_TEST (INLINE_TEST_${inline_test})
ENDFOREACH (inline_test)

SET (HDF5_EXTRA_TEST_DEFINITIONS)
IF (INLINE_TEST___inline__)
  SET (H5_inline __inline__)
ELSE (INLINE_TEST___inline__)
  IF (INLINE_TEST___inline)
    SET (H5_inline __inline)
  ELSE (INLINE_TEST___inline)
    IF (INLINE_TEST_inline)
      SET (H5_inline inline)
    ENDIF (INLINE_TEST_inline)
  ENDIF (INLINE_TEST___inline)
ENDIF (INLINE_TEST___inline__)

#-----------------------------------------------------------------------------
# Check how to print a Long Long integer
#-----------------------------------------------------------------------------
IF (NOT H5_PRINTF_LL_WIDTH OR H5_PRINTF_LL_WIDTH MATCHES "unknown")
  SET (PRINT_LL_FOUND 0)
  MESSAGE (STATUS "Checking for appropriate format for 64 bit long:")
  FOREACH (HDF5_PRINTF_LL l64 l L q I64 ll)
    SET (CURRENT_TEST_DEFINITIONS "-DPRINTF_LL_WIDTH=${HDF5_PRINTF_LL}")
    IF (H5_SIZEOF_LONG_LONG)
      SET (CURRENT_TEST_DEFINITIONS "${CURRENT_TEST_DEFINITIONS} -DHAVE_LONG_LONG")
    ENDIF (H5_SIZEOF_LONG_LONG)
    TRY_RUN (HDF5_PRINTF_LL_TEST_RUN   HDF5_PRINTF_LL_TEST_COMPILE
        ${HDF5_BINARY_DIR}/CMake
        ${HDF5_RESOURCES_DIR}/HDF5Tests.c
        CMAKE_FLAGS -DCOMPILE_DEFINITIONS:STRING=${CURRENT_TEST_DEFINITIONS}
        OUTPUT_VARIABLE OUTPUT
    )
    IF (HDF5_PRINTF_LL_TEST_COMPILE)
      IF (HDF5_PRINTF_LL_TEST_RUN MATCHES 0)
        SET (H5_PRINTF_LL_WIDTH "\"${HDF5_PRINTF_LL}\"" CACHE INTERNAL "Width for printf for type `long long' or `__int64', us. `ll")
        SET (PRINT_LL_FOUND 1)
      ELSE (HDF5_PRINTF_LL_TEST_RUN MATCHES 0)
        MESSAGE ("Width with ${HDF5_PRINTF_LL} failed with result: ${HDF5_PRINTF_LL_TEST_RUN}")
      ENDIF (HDF5_PRINTF_LL_TEST_RUN MATCHES 0)
    ELSE (HDF5_PRINTF_LL_TEST_COMPILE)
      FILE (APPEND ${CMAKE_BINARY_DIR}/CMakeFiles/CMakeError.log
          "Test H5_PRINTF_LL_WIDTH for ${HDF5_PRINTF_LL} failed with the following output:\n ${OUTPUT}\n"
      )
    ENDIF (HDF5_PRINTF_LL_TEST_COMPILE)
  ENDFOREACH (HDF5_PRINTF_LL)

  IF (PRINT_LL_FOUND)
    MESSAGE (STATUS "Checking for apropriate format for 64 bit long: found ${H5_PRINTF_LL_WIDTH}")
  ELSE (PRINT_LL_FOUND)
    MESSAGE (STATUS "Checking for apropriate format for 64 bit long: not found")
    SET (H5_PRINTF_LL_WIDTH "\"unknown\"" CACHE INTERNAL
        "Width for printf for type `long long' or `__int64', us. `ll"
    )
  ENDIF (PRINT_LL_FOUND)
ENDIF (NOT H5_PRINTF_LL_WIDTH OR H5_PRINTF_LL_WIDTH MATCHES "unknown")

# ----------------------------------------------------------------------
# Set the flag to indicate that the machine can handle converting
# denormalized floating-point values.
# (This flag should be set for all machines, except for the Crays, where
# the cache value is set in it's config file)
#
SET (H5_CONVERT_DENORMAL_FLOAT 1)

#-----------------------------------------------------------------------------
#  Are we going to use HSIZE_T
#-----------------------------------------------------------------------------
IF (HDF5_ENABLE_HSIZET)
  SET (H5_HAVE_LARGE_HSIZET 1)
ENDIF (HDF5_ENABLE_HSIZET)

#-----------------------------------------------------------------------------
# Macro to determine the various conversion capabilities
#-----------------------------------------------------------------------------
MACRO (H5ConversionTests TEST msg)
  IF ("${TEST}" MATCHES "^${TEST}$")
   # MESSAGE (STATUS "===> ${TEST}")
    TRY_RUN (${TEST}_RUN   ${TEST}_COMPILE
        ${HDF5_BINARY_DIR}/CMake
        ${HDF5_RESOURCES_DIR}/ConversionTests.c
        CMAKE_FLAGS -DCOMPILE_DEFINITIONS:STRING=-D${TEST}_TEST
        OUTPUT_VARIABLE OUTPUT
    )
    IF (${TEST}_COMPILE)
      IF (${TEST}_RUN  MATCHES 0)
        SET (${TEST} 1 CACHE INTERNAL ${msg})
        MESSAGE (STATUS "${msg}... yes")
      ELSE (${TEST}_RUN  MATCHES 0)
        SET (${TEST} "" CACHE INTERNAL ${msg})
        MESSAGE (STATUS "${msg}... no")
        FILE (APPEND ${CMAKE_BINARY_DIR}/CMakeFiles/CMakeError.log
              "Test ${TEST} Run failed with the following output and exit code:\n ${OUTPUT}\n"
        )
      ENDIF (${TEST}_RUN  MATCHES 0)
    ELSE (${TEST}_COMPILE )
      SET (${TEST} "" CACHE INTERNAL ${msg})
      MESSAGE (STATUS "${msg}... no")
      FILE (APPEND ${CMAKE_BINARY_DIR}/CMakeFiles/CMakeError.log
          "Test ${TEST} Compile failed with the following output:\n ${OUTPUT}\n"
      )
    ENDIF (${TEST}_COMPILE)

  ENDIF("${TEST}" MATCHES "^${TEST}$")
ENDMACRO (H5ConversionTests)

#-----------------------------------------------------------------------------
# Macro to make some of the conversion tests easier to write/read
#-----------------------------------------------------------------------------
MACRO (H5MiscConversionTest  VAR TEST msg)
  IF ("${TEST}" MATCHES "^${TEST}$")
    IF (${VAR})
      SET (${TEST} 1 CACHE INTERNAL ${msg})
      MESSAGE (STATUS "${msg}... yes")
    ELSE (${VAR})
      SET (${TEST} "" CACHE INTERNAL ${msg})
      MESSAGE (STATUS "${msg}... no")
    ENDIF (${VAR})
  ENDIF ("${TEST}" MATCHES "^${TEST}$")
ENDMACRO (H5MiscConversionTest)

#-----------------------------------------------------------------------------
# Check various conversion capabilities
#-----------------------------------------------------------------------------

# -----------------------------------------------------------------------
# Set flag to indicate that the machine can handle conversion from
# long double to integers accurately.  This flag should be set "yes" for
# all machines except all SGIs.  For SGIs, some conversions are
# incorrect and its cache value is set "no" in its config/irix6.x and
# irix5.x.
#
H5MiscConversionTest (H5_SIZEOF_LONG_DOUBLE H5_LDOUBLE_TO_INTEGER_ACCURATE "checking IF converting from long double to integers is accurate")
# -----------------------------------------------------------------------
# Set flag to indicate that the machine can do conversion from
# long double to integers regardless of accuracy.  This flag should be
# set "yes" for all machines except HP-UX 11.00.  For HP-UX 11.00, the
# compiler has 'floating exception' when converting 'long double' to all
# integers except 'unsigned long long'.  Other HP-UX systems are unknown
# yet. (1/8/05 - SLU)
#
H5ConversionTests (H5_LDOUBLE_TO_INTEGER_WORKS "Checking IF converting from long double to integers works")
# -----------------------------------------------------------------------
# Set flag to indicate that the machine can handle conversion from
# integers to long double.  (This flag should be set "yes" for all
# machines except all SGIs, where some conversions are
# incorrect and its cache value is set "no" in its config/irix6.x and
# irix5.x)
#
H5MiscConversionTest (H5_SIZEOF_LONG_DOUBLE H5_INTEGER_TO_LDOUBLE_ACCURATE "checking IF accurately converting from integers to long double")
# ----------------------------------------------------------------------
# Set the flag to indicate that the machine can accurately convert
# 'unsigned long' to 'float' values.
# (This flag should be set for all machines, except for Pathscale compiler
# on Sandia's Linux machine where the compiler interprets 'unsigned long'
# values as negative when the first bit of 'unsigned long' is on during
# the conversion to float.)
#
H5ConversionTests (H5_ULONG_TO_FLOAT_ACCURATE "Checking IF accurately converting unsigned long to float values")
# ----------------------------------------------------------------------
# Set the flag to indicate that the machine can accurately convert
# 'unsigned (long) long' values to 'float' and 'double' values.
# (This flag should be set for all machines, except for the SGIs, where
# the cache value is set in the config/irix6.x config file) and Solaris
# 64-bit machines, where the short program below tests if round-up is
# correctly handled.
#
H5ConversionTests (H5_ULONG_TO_FP_BOTTOM_BIT_ACCURATE "Checking IF accurately converting unsigned long long to floating-point values")
# ----------------------------------------------------------------------
# Set the flag to indicate that the machine can accurately convert
# 'float' or 'double' to 'unsigned long long' values.
# (This flag should be set for all machines, except for PGI compiler
# where round-up happens when the fraction of float-point value is greater
# than 0.5.
#
H5ConversionTests (H5_FP_TO_ULLONG_ACCURATE "Checking IF accurately roundup converting floating-point to unsigned long long values" )
# ----------------------------------------------------------------------
# Set the flag to indicate that the machine can accurately convert
# 'float', 'double' or 'long double' to 'unsigned long long' values.
# (This flag should be set for all machines, except for HP-UX machines
# where the maximal number for unsigned long long is 0x7fffffffffffffff
# during conversion.
#
H5ConversionTests (H5_FP_TO_ULLONG_RIGHT_MAXIMUM "Checking IF right maximum converting floating-point to unsigned long long values" )
# ----------------------------------------------------------------------
# Set the flag to indicate that the machine can accurately convert
# 'long double' to 'unsigned int' values.  (This flag should be set for
# all machines, except for some Intel compilers on some Linux.)
#
H5ConversionTests (H5_LDOUBLE_TO_UINT_ACCURATE "Checking IF correctly converting long double to unsigned int values")
# ----------------------------------------------------------------------
# Set the flag to indicate that the machine can _compile_
# 'unsigned long long' to 'float' and 'double' typecasts.
# (This flag should be set for all machines.)
#
IF (H5_ULLONG_TO_FP_CAST_WORKS MATCHES ^H5_ULLONG_TO_FP_CAST_WORKS$)
  SET (H5_ULLONG_TO_FP_CAST_WORKS 1 CACHE INTERNAL "Checking IF compiling unsigned long long to floating-point typecasts work")
  MESSAGE (STATUS "Checking IF compiling unsigned long long to floating-point typecasts work... yes")
ENDIF (H5_ULLONG_TO_FP_CAST_WORKS MATCHES ^H5_ULLONG_TO_FP_CAST_WORKS$)
# ----------------------------------------------------------------------
# Set the flag to indicate that the machine can _compile_
# 'long long' to 'float' and 'double' typecasts.
# (This flag should be set for all machines.)
#
IF (H5_LLONG_TO_FP_CAST_WORKS MATCHES ^H5_LLONG_TO_FP_CAST_WORKS$)
  SET (H5_LLONG_TO_FP_CAST_WORKS 1 CACHE INTERNAL "Checking IF compiling long long to floating-point typecasts work")
  MESSAGE (STATUS "Checking IF compiling long long to floating-point typecasts work... yes")
ENDIF (H5_LLONG_TO_FP_CAST_WORKS MATCHES ^H5_LLONG_TO_FP_CAST_WORKS$)
# ----------------------------------------------------------------------
# Set the flag to indicate that the machine can convert from
# 'unsigned long long' to 'long double' without precision loss.
# (This flag should be set for all machines, except for FreeBSD(sleipnir)
# where the last 2 bytes of mantissa are lost when compiler tries to do
# the conversion, and Cygwin where compiler doesn't do rounding correctly.)
#
H5ConversionTests (H5_ULLONG_TO_LDOUBLE_PRECISION "Checking IF converting unsigned long long to long double with precision")
# ----------------------------------------------------------------------
# Set the flag to indicate that the machine can handle overflow converting
# all floating-point to all integer types.
# (This flag should be set for all machines, except for Cray X1 where
# floating exception is generated when the floating-point value is greater
# than the maximal integer value).
#
H5ConversionTests (H5_FP_TO_INTEGER_OVERFLOW_WORKS  "Checking IF overflows normally converting floating-point to integer values")
# ----------------------------------------------------------------------
# Set the flag to indicate that the machine is using a special algorithm to convert
# 'long double' to '(unsigned) long' values.  (This flag should only be set for 
# the IBM Power6 Linux.  When the bit sequence of long double is 
# 0x4351ccf385ebc8a0bfcc2a3c3d855620, the converted value of (unsigned)long 
# is 0x004733ce17af227f, not the same as the library's conversion to 0x004733ce17af2282.
# The machine's conversion gets the correct value.  We define the macro and disable
# this kind of test until we figure out what algorithm they use.
#
IF (H5_LDOUBLE_TO_LONG_SPECIAL MATCHES ^H5_LDOUBLE_TO_LONG_SPECIAL$)
  SET (H5_LDOUBLE_TO_LONG_SPECIAL 0 CACHE INTERNAL "Define if your system converts long double to (unsigned) long values with special algorithm")
  MESSAGE (STATUS "Checking IF your system converts long double to (unsigned) long values with special algorithm... no")
ENDIF (H5_LDOUBLE_TO_LONG_SPECIAL MATCHES ^H5_LDOUBLE_TO_LONG_SPECIAL$)
# ----------------------------------------------------------------------
# Set the flag to indicate that the machine is using a special algorithm
# to convert some values of '(unsigned) long' to 'long double' values.  
# (This flag should be off for all machines, except for IBM Power6 Linux, 
# when the bit sequences are 003fff..., 007fff..., 00ffff..., 01ffff..., 
# ..., 7fffff..., the compiler uses a unknown algorithm.  We define a 
# macro and skip the test for now until we know about the algorithm.
#
IF (H5_LONG_TO_LDOUBLE_SPECIAL MATCHES ^H5_LONG_TO_LDOUBLE_SPECIAL$)
  SET (H5_LONG_TO_LDOUBLE_SPECIAL 0 CACHE INTERNAL "Define if your system can convert (unsigned) long to long double values with special algorithm")
  MESSAGE (STATUS "Checking IF your system can convert (unsigned) long to long double values with special algorithm... no")
ENDIF (H5_LONG_TO_LDOUBLE_SPECIAL MATCHES ^H5_LONG_TO_LDOUBLE_SPECIAL$)
# ----------------------------------------------------------------------
# Set the flag to indicate that the machine can accurately convert
# 'long double' to '(unsigned) long long' values.  (This flag should be set for
# all machines, except for Mac OS 10.4 and SGI IRIX64 6.5.  When the bit sequence
# of long double is 0x4351ccf385ebc8a0bfcc2a3c..., the values of (unsigned)long long
# start to go wrong on these two machines.  Adjusting it higher to
# 0x4351ccf385ebc8a0dfcc... or 0x4351ccf385ebc8a0ffcc... will make the converted
# values wildly wrong.  This test detects this wrong behavior and disable the test.
#
H5ConversionTests (H5_LDOUBLE_TO_LLONG_ACCURATE "Checking IF correctly converting long double to (unsigned) long long values")
# ----------------------------------------------------------------------
# Set the flag to indicate that the machine can accurately convert
# '(unsigned) long long' to 'long double' values.  (This flag should be set for
# all machines, except for Mac OS 10.4, when the bit sequences are 003fff...,
# 007fff..., 00ffff..., 01ffff..., ..., 7fffff..., the converted values are twice
# as big as they should be.
#
H5ConversionTests (H5_LLONG_TO_LDOUBLE_CORRECT "Checking IF correctly converting (unsigned) long long to long double values")
# ----------------------------------------------------------------------
# Set the flag to indicate that the machine generates bad code
# for the H5V_log2_gen() routine in src/H5Vprivate.h
# (This flag should be set to no for all machines, except for SGI IRIX64,
# where the cache value is set to yes in it's config file)
#
IF (H5_BAD_LOG2_CODE_GENERATED MATCHES ^H5_BAD_LOG2_CODE_GENERATED$)
  SET (H5_BAD_LOG2_CODE_GENERATED 0 CACHE INTERNAL "Define if your system generates wrong code for log2 routine")
  MESSAGE (STATUS "Checking IF your system generates wrong code for log2 routine... no")
ENDIF (H5_BAD_LOG2_CODE_GENERATED MATCHES ^H5_BAD_LOG2_CODE_GENERATED$)
# ----------------------------------------------------------------------
# Check if pointer alignments are enforced
#
H5ConversionTests (H5_NO_ALIGNMENT_RESTRICTIONS "Checking IF alignment restrictions are strictly enforced")

# Define a macro for Cygwin (on XP only) where the compiler has rounding
#   problem converting from unsigned long long to long double */
IF (CYGWIN)
  SET (H5_CYGWIN_ULLONG_TO_LDOUBLE_ROUND_PROBLEM 1)
ENDIF (CYGWIN)
