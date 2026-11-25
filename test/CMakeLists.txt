# SPDX-License-Identifier: Apache-2.0

FOREACH(SVARIANT ${SVARIANT_S})

    string(TOLOWER ${SVARIANT} SVARIANT_LOWER)
    string(TOUPPER ${SVARIANT} SVARIANT_UPPER)

    add_executable(sqisign_test_kat_${SVARIANT} test_kat.c)
    target_link_libraries(sqisign_test_kat_${SVARIANT} sqisign_${SVARIANT_LOWER}_test)
    target_include_directories(sqisign_test_kat_${SVARIANT} PRIVATE ${PROJECT_SOURCE_DIR}/src/nistapi/${SVARIANT_LOWER} ${INC_PUBLIC} ${INC_INTBIG} ${INC_PRECOMP_${SVARIANT_UPPER}} ${INC_QUATERNION} ${INC_KLPT} ${INC_GF_${SVARIANT_UPPER}} ${INC_EC} ${INC_COMMON} ${INC_ID2ISO} ${INC_PROTOCOLS})
    add_test(NAME sqisign_${SVARIANT}_KAT COMMAND sqisign_test_kat_${SVARIANT})
    set_tests_properties(sqisign_${SVARIANT}_KAT PROPERTIES
        WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/KAT" TIMEOUT 7200)

    add_executable(sqisign_test_parameter_encoding_${SVARIANT} test_parameter_encoding.c)
    target_link_libraries(sqisign_test_parameter_encoding_${SVARIANT} sqisign_${SVARIANT}_test)
    target_include_directories(sqisign_test_parameter_encoding_${SVARIANT} PRIVATE
        ${PROJECT_SOURCE_DIR}/src/nistapi/${SVARIANT_LOWER})
    add_test(NAME sqisign_${SVARIANT}_PARAMETER_ENCODING COMMAND sqisign_test_parameter_encoding_${SVARIANT})

    if (ENABLE_SIGN)

        add_executable(sqisign_test_random_ideal_precision_${SVARIANT} test_random_ideal_precision.c)
        target_link_libraries(sqisign_test_random_ideal_precision_${SVARIANT}
            ${LIB_QUATERNION_${SVARIANT_UPPER}}
            ${LIB_PRECOMP_${SVARIANT_UPPER}}
            sqisign_common_test)
        target_include_directories(sqisign_test_random_ideal_precision_${SVARIANT} PRIVATE
            ${INC_PUBLIC}
            ${INC_COMMON}
            ${INC_QUATERNION}
            ${INC_PRECOMP_${SVARIANT_UPPER}}
            ${PROJECT_SOURCE_DIR}/src/quaternion/ref/generic/internal_quaternion_headers)
        if (SVARIANT STREQUAL "lvl1")
            target_compile_definitions(sqisign_test_random_ideal_precision_${SVARIANT} PRIVATE
                EXPECTED_IBZ_LIMBS=150 PAPER_IBZ_BITS=9555 LINE13_IBZ_BITS=9565 EXPECTED_MIX_BITS=582)
        elseif (SVARIANT STREQUAL "lvl3")
            target_compile_definitions(sqisign_test_random_ideal_precision_${SVARIANT} PRIVATE
                EXPECTED_IBZ_LIMBS=228 PAPER_IBZ_BITS=14572 LINE13_IBZ_BITS=14569 EXPECTED_MIX_BITS=889)
        elseif (SVARIANT STREQUAL "lvl5")
            target_compile_definitions(sqisign_test_random_ideal_precision_${SVARIANT} PRIVATE
                EXPECTED_IBZ_LIMBS=300 PAPER_IBZ_BITS=19145 LINE13_IBZ_BITS=19145 EXPECTED_MIX_BITS=1181)
        endif()
        add_test(NAME sqisign_${SVARIANT}_RANDOM_IDEAL_PRECISION
                 COMMAND sqisign_test_random_ideal_precision_${SVARIANT})
        set_tests_properties(sqisign_${SVARIANT}_RANDOM_IDEAL_PRECISION PROPERTIES TIMEOUT 180)

        add_executable(sqisign_test_scheme_${SVARIANT} test_sqisign.c)
        # Use the deterministic test RNG so --seed reproduces the complete
        # KeyGen/Sign/Verify transcript.
        target_link_libraries(sqisign_test_scheme_${SVARIANT} sqisign_${SVARIANT_LOWER}_test)
        target_include_directories(sqisign_test_scheme_${SVARIANT} PUBLIC ${PROJECT_SOURCE_DIR}/src/common ${INC_PUBLIC} ${PROJECT_SOURCE_DIR}/src/nistapi/${SVARIANT_LOWER})
        if(APPLE)
            # find_uv uses parameter-sized VLAs.  lvl3/lvl5 exceed macOS's
            # default 8 MiB (and shell hard 64 MiB) stack; encode the stack
            # required by the reference self-test in LC_MAIN.
            target_link_options(sqisign_test_scheme_${SVARIANT} PRIVATE "-Wl,-stack_size,0x20000000")
        endif()
        add_test(NAME sqisign_${SVARIANT}_SELFTEST
                 COMMAND sqisign_test_scheme_${SVARIANT} --seed=1 --msglen=32)
        set_tests_properties(sqisign_${SVARIANT}_SELFTEST PROPERTIES TIMEOUT 7200)

    endif()

ENDFOREACH()
