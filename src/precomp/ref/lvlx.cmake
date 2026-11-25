set(SOURCE_FILES_PRECOMP_${SVARIANT_UPPER}_REF
    ec_params.c
    e0_basis.c
    hd_splitting_transforms.c
)

if(ENABLE_SIGN)
    set(SOURCE_FILES_PRECOMP_${SVARIANT_UPPER}_REF
        ${SOURCE_FILES_PRECOMP_${SVARIANT_UPPER}_REF}
        torsion_constants.c
        quaternion_data.c
        endomorphism_action.c
    )
endif()

add_library(${LIB_PRECOMP_${SVARIANT_UPPER}} STATIC ${SOURCE_FILES_PRECOMP_${SVARIANT_UPPER}_REF})
target_include_directories(${LIB_PRECOMP_${SVARIANT_UPPER}} PRIVATE ${INC_QUATERNION} ${PROJECT_SOURCE_DIR}/quaternion/ref/generic/include ${INC_PUBLIC} ${PROJECT_SOURCE_DIR}/src/hd/ref/include ${PROJECT_SOURCE_DIR}/src/ec/ref/include ${PROJECT_SOURCE_DIR}/src/ec/ref/${SVARIANT_LOWER}/include ${INC_GF} ${INC_GF_${SVARIANT_UPPER}} ${INC_PRECOMP_${SVARIANT_UPPER}} ${INC_COMMON})
target_compile_options(${LIB_PRECOMP_${SVARIANT_UPPER}} PRIVATE ${C_OPT_FLAGS})
target_compile_definitions(${LIB_PRECOMP_${SVARIANT_UPPER}} PUBLIC SQISIGN_VARIANT=${SVARIANT_LOWER})

find_program(SAGEMATH sage)
set(PRECOMP_SAGE ${CMAKE_COMMAND} -E env
    "PYTHONPATH=${PROJECT_SOURCE_DIR}/scripts/precomp:$ENV{PYTHONPATH}"
    "PYTHONDONTWRITEBYTECODE=1" "PYTHONUNBUFFERED=1" "${SAGEMATH}")
add_custom_target(precomp_${SVARIANT_LOWER}
    DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/sqisign_parameters.txt"
    COMMAND ${PRECOMP_SAGE} "${PROJECT_SOURCE_DIR}/scripts/precomp/precompute_torsion_constants.sage"
    COMMAND ${PRECOMP_SAGE} "${PROJECT_SOURCE_DIR}/scripts/precomp/precompute_quaternion_constants.sage"
    COMMAND ${PRECOMP_SAGE} "${PROJECT_SOURCE_DIR}/scripts/precomp/precompute_sizes.sage"
    COMMAND ${PRECOMP_SAGE} "${PROJECT_SOURCE_DIR}/scripts/precomp/precompute_field_constants.sage"
    COMMAND ${PRECOMP_SAGE} "${PROJECT_SOURCE_DIR}/scripts/precomp/precompute_quaternion_data.sage"
    COMMAND ${PRECOMP_SAGE} "${PROJECT_SOURCE_DIR}/scripts/precomp/precompute_endomorphism_action.sage"
    COMMAND ${PRECOMP_SAGE} "${PROJECT_SOURCE_DIR}/scripts/precomp/ec_params.sage"
    COMMAND ${PRECOMP_SAGE} "${PROJECT_SOURCE_DIR}/scripts/precomp/precompute_hd_splitting.sage"
    COMMAND ${PRECOMP_SAGE} "${PROJECT_SOURCE_DIR}/scripts/precomp/precompute_E0_basis.sage"
    WORKING_DIRECTORY
    "${CMAKE_CURRENT_SOURCE_DIR}"
    VERBATIM
)

set_directory_properties(PROPERTIES CLEAN_NO_CUSTOM true)
set_target_properties(precomp_${SVARIANT_LOWER} PROPERTIES EXCLUDE_FROM_ALL TRUE)

add_dependencies(precomp precomp_${SVARIANT_LOWER})
