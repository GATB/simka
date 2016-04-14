# vfdTest.cmake executes a command and captures the output in a file. Command uses specified VFD.
# Exit status of command can also be compared.

# arguments checking
IF (NOT TEST_PROGRAM)
  MESSAGE (FATAL_ERROR "Require TEST_PROGRAM to be defined")
ENDIF (NOT TEST_PROGRAM)
#IF (NOT TEST_ARGS)
#  MESSAGE (STATUS "Require TEST_ARGS to be defined")
#ENDIF (NOT TEST_ARGS)
#IF (NOT TEST_EXPECT)
#  MESSAGE (STATUS "Require TEST_EXPECT to be defined")
#ENDIF (NOT TEST_EXPECT)
IF (NOT TEST_FOLDER)
  MESSAGE ( FATAL_ERROR "Require TEST_FOLDER to be defined")
ENDIF (NOT TEST_FOLDER)
IF (NOT TEST_VFD)
  MESSAGE (FATAL_ERROR "Require TEST_VFD to be defined")
ENDIF (NOT TEST_VFD)

SET (ERROR_APPEND 1)

MESSAGE (STATUS "USING ${TEST_VFD} ON COMMAND: ${TEST_PROGRAM} ${TEST_ARGS}")

SET (ENV{HDF5_DRIVER} "${TEST_VFD}")
# run the test program, capture the stdout/stderr and the result var
EXECUTE_PROCESS (
    COMMAND ${TEST_PROGRAM} ${TEST_ARGS}
    WORKING_DIRECTORY ${TEST_FOLDER}
    OUTPUT_FILE ${TEST_OUTPUT}_${TEST_VFD}.out
    ERROR_FILE ${TEST_OUTPUT}_${TEST_VFD}.err
    OUTPUT_VARIABLE TEST_ERROR
    ERROR_VARIABLE TEST_ERROR
)

MESSAGE (STATUS "COMMAND Result: ${TEST_RESULT}")

IF (ERROR_APPEND)
  FILE (READ ${TEST_FOLDER}/${TEST_OUTPUT}_${TEST_VFD}.err TEST_STREAM)
  FILE (APPEND ${TEST_FOLDER}/${TEST_OUTPUT}_${TEST_VFD}.out "${TEST_STREAM}") 
ENDIF (ERROR_APPEND)

# if the return value is !=${TEST_EXPECT} bail out
IF (NOT ${TEST_RESULT} STREQUAL ${TEST_EXPECT})
  MESSAGE ( FATAL_ERROR "Failed: Test program ${TEST_PROGRAM} exited != ${TEST_EXPECT}.\n${TEST_ERROR}")
ENDIF (NOT ${TEST_RESULT} STREQUAL ${TEST_EXPECT})

# everything went fine...
MESSAGE ("Passed: The ${TEST_PROGRAM} program used vfd ${TEST_VFD}")

