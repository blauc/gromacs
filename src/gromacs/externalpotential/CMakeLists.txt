file(GLOB EXTERNALPOTENTIALMODULE_SOURCES *.cpp)
set(LIBGROMACS_SOURCES ${LIBGROMACS_SOURCES} ${EXTERNALPOTENTIALMODULE_SOURCES} PARENT_SCOPE)

gmx_install_headers(
    externalpotential.h
    externalpotentialutil.h
    externalpotentialdata.h
    externalpotentialregistration.h
    )

add_subdirectory(modules)

if (BUILD_TESTING)
    add_subdirectory(tests)
endif()
