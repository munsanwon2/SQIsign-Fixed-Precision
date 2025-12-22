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
    printer.c
    lll/rationals.c
    lll/l2.c
    lll/lll_verification.c
    lll/lll_applications.c
    lll/rationals.c
    normeq.c
    hnf/ibz_division.c
    hnf/hnf_internal.c
    hnf/hnf.c
    test/random_input_generation.c
    dpe.c
)

# Build separate libraries for each security level with appropriate SQISIGN_VARIANT
FOREACH(SVARIANT ${SVARIANT_S})
    string(TOUPPER ${SVARIANT} SVARIANT_UPPER)
    string(TOLOWER ${SVARIANT} SVARIANT_LOWER)

    set(LIB_NAME sqisign_quaternion_${SVARIANT_LOWER})

    add_library(${LIB_NAME} STATIC ${SOURCE_FILES_QUATERNION_GENERIC_REF})
    target_link_libraries(${LIB_NAME} m)
    target_include_directories(${LIB_NAME} PRIVATE common ${INC_PUBLIC} ${INC_COMMON} ${INC_QUATERNION} internal_quaternion_headers)
    target_compile_options(${LIB_NAME} PRIVATE ${C_OPT_FLAGS})
    target_compile_definitions(${LIB_NAME} PUBLIC SQISIGN_VARIANT=${SVARIANT_LOWER})

    # Update the cache variables to point to the level-specific libraries
    set(LIB_QUATERNION_${SVARIANT_UPPER} ${LIB_NAME} CACHE INTERNAL "LIB")
ENDFOREACH()

# For backward compatibility, create an alias for the generic library (defaults to lvl5)
set(LIB_QUATERNION sqisign_quaternion_lvl5 CACHE INTERNAL "LIB")

add_subdirectory(test)
