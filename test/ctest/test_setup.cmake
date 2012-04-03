## Don't edit this file unless you really know what you are doing
# See ctest_hoomd.cmake for complete documentation

# lets create a build name to idenfity all these options in a string. It will look like
# Linux-gcc412-trunk-single-cuda (for a single precision build with cuda)
# Linux-gcc412-hoomd-0.8-double (for a double precision build without cuda and in the hoomd-0.8 branch)
SET (BUILDNAME "${SYSTEM_NAME}-${COMPILER_NAME}-${HOOMD_BRANCH}")

if (ENABLE_STATIC MATCHES "ON")
    SET (BUILDNAME "${BUILDNAME}-static")
else(ENABLE_STATIC MATCHES "ON")
    SET (BUILDNAME "${BUILDNAME}-shared")
endif(ENABLE_STATIC MATCHES "ON")

if (SINGLE_PRECISION MATCHES "ON")
    SET (BUILDNAME "${BUILDNAME}-single")
else (SINGLE_PRECISION MATCHES "ON")
    SET (BUILDNAME "${BUILDNAME}-double")
endif (SINGLE_PRECISION MATCHES "ON")

if (ENABLE_CUDA MATCHES "ON")
    SET (BUILDNAME "${BUILDNAME}-cuda")
endif (ENABLE_CUDA MATCHES "ON")

if (ENABLE_OPENMP MATCHES "ON")
    SET (BUILDNAME "${BUILDNAME}-openmp")
endif (ENABLE_OPENMP MATCHES "ON")

if (NOT BUILD_TYPE MATCHES "Release")
    SET (BUILDNAME "${BUILDNAME}-${BUILD_TYPE}")
endif (NOT BUILD_TYPE MATCHES "Release")

if (ENABLE_COVERAGE)
    SET (COVERAGE_FLAGS "-fprofile-arcs -ftest-coverage")
endif (ENABLE_COVERAGE)

SET (CTEST_COMMAND "ctest -D ${TEST_GROUP} ${IGNORE_TESTS}")
if (MEMORYCHECK_COMMAND)
    set (CTEST_COMMAND "${CTEST_COMMAND}" 
            "ctest -D ${TEST_GROUP}MemCheck -D ${TEST_GROUP}Submit ${IGNORE_TESTS}")
endif (MEMORYCHECK_COMMAND)

SET (CTEST_INITIAL_CACHE "
CMAKE_GENERATOR:INTERNAL=Unix Makefiles
MAKECOMMAND:STRING=/usr/bin/make -i
BUILDNAME:STRING=${BUILDNAME}
SITE:STRING=${SITE_NAME}
CMAKE_BUILD_TYPE:STRING=${BUILD_TYPE}
ENABLE_CUDA:BOOL=${ENABLE_CUDA}
ENABLE_DOXYGEN:BOOL=OFF
ENABLE_OPENMP:BOOL=${ENABLE_OPENMP}
SINGLE_PRECISION:BOOL=${SINGLE_PRECISION}
ENABLE_STATIC:BOOL=${ENABLE_STATIC}
ENABLE_TEST_ALL:BOOL="ON"
CMAKE_C_FLAGS:STRING=${COVERAGE_FLAGS}
CMAKE_CXX_FLAGS:STRING=${COVERAGE_FLAGS}
MEMORYCHECK_COMMAND:FILEPATH=${MEMORYCHECK_COMMAND}
MEMORYCHECK_SUPPRESSIONS_FILE:FILEPATH=${CTEST_CHECKOUT_DIR}/test/unit/combined_valgrind.supp
CUDA_ARCH_LIST:STRING=${CUDA_ARCH_LIST}
")
