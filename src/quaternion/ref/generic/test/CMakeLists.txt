set(SOURCE_FILES_QUATERNION_GENERIC_REF
    intbig.c
    algebra.c
    ideal.c
    dim4.c
    dim2.c
    integers.c
    lattice.c
    lat_ball.c
    finit.c
    normeq.c
    randomized.c
)

# if (NOT GMP_LIBRARY STREQUAL "MINI")
#     list(APPEND SOURCE_FILES_QUATERNION_GENERIC_REF ${PROJECT_SOURCE_DIR}/src/mini-gmp/mini-gmp-extra.c)
# endif()

set(SOURCE_FILES_QUATERNION_GENERIC_REF_TESTS
    ${SOURCE_FILES_QUATERNION_GENERIC_REF}
    ../hnf/hnf_tests.c
    ../lll/lll_verification.c
    ../lll/lll_tests.c
)

# Build tests for each security level
FOREACH(SVARIANT ${SVARIANT_S})
    string(TOUPPER ${SVARIANT} SVARIANT_UPPER)
    string(TOLOWER ${SVARIANT} SVARIANT_LOWER)

    set(TEST_NAME sqisign_test_quaternion_${SVARIANT_LOWER})

    add_executable(${TEST_NAME} ${SOURCE_FILES_QUATERNION_GENERIC_REF_TESTS} test_quaternions.c)
    target_link_libraries(${TEST_NAME} ${LIB_QUATERNION_${SVARIANT_UPPER}} sqisign_common_sys)
    target_include_directories(${TEST_NAME} PRIVATE ../internal_quaternion_headers ${INC_COMMON} ${INC_QUATERNION} ${INC_PUBLIC} ${PROJECT_SOURCE_DIR}/src/mini-gmp)
    target_compile_definitions(${TEST_NAME} PUBLIC SQISIGN_VARIANT=${SVARIANT_LOWER})

    # MSAN and GMP lead to false positives, see
    # https://gmplib.org/list-archives/gmp-bugs/2019-March/004529.html
    if(NOT CMAKE_BUILD_TYPE STREQUAL "MSAN")
        add_test(${TEST_NAME} ${TEST_NAME})
    endif()
ENDFOREACH()

# Note: sqisign_bm_quaternion benchmark target needs to be defined
# Currently this is a placeholder and the actual benchmark executable is not implemented
