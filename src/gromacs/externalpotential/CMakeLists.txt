file(GLOB EXTERNALPOTENTIALMODULE_SOURCES *.cpp modules/*/*.cpp)
set(LIBGROMACS_SOURCES ${LIBGROMACS_SOURCES} ${EXTERNALPOTENTIALMODULE_SOURCES} PARENT_SCOPE)

gmx_install_headers(
    externalpotential.h
    externalpotentialmanager.h
    externalpotentialregistration.h
    )

if (BUILD_TESTING)
    add_subdirectory(tests)
endif()
