cmake_minimum_required(VERSION 2.8)

project(somatic_sniper)

set(CMAKE_INSTALL_PREFIX "/usr")
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake;${CMAKE_SOURCE_DIR}/build-common/cmake")
include(TestHelper)
include(VersionHelper)
include(ProjectHelper)

# NOTE: for sniper we want the exe suffix to be like 0.7.4, not just 0.7
#       which is the default from VersionHelper
# I have diverged Sniper from normal build-common to support 4 part versioning where the last part is an (optional) patch version
#set(EXE_VERSION_SUFFIX ${FULL_VERSION})

configure_file(version/version.h.in version/version.h @ONLY)
include_directories(${PROJECT_BINARY_DIR}/version)

include(BuildSamtools)
include_directories(${Samtools_INCLUDE_DIRS})

add_custom_target(deps ALL)
add_dependencies(deps samtools-lib)

# on linux mint 12, libgtest requires pthreads, so link thread libraries
# everywhere
find_package(Threads)
link_libraries(${CMAKE_THREAD_LIBS_INIT})

if (CMAKE_COMPILER_IS_GNUCC)
    set(CMAKE_C_FLAGS "-Wall")
endif ()

if (CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "-Wall -ansi -pedantic -Wno-long-long -Wno-variadic-macros")
endif ()

# make sure to pick up headers from library dirs
include_directories("src/lib")

# unit tests
enable_testing(true)
add_projects(test/lib)

# main project

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
add_projects(src/lib)
add_projects(src/exe)

# integration testing
add_subdirectory(integration-test)

# generate documentation
add_subdirectory(docs)

include(PackageInfo.cmake)
