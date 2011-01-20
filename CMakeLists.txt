cmake_minimum_required(VERSION 2.6)

project(somatic_sniper)

# TODO: replace with FindSamtools.cmake
set(SAMTOOLS_INC /gsc/pkg/bio/samtools/samtools-0.1.6)
set(SAMTOOLS_LIB /gsc/pkg/bio/samtools/samtools-0.1.6)
include_directories(${SAMTOOLS_INC})
link_directories(${SAMTOOLS_LIB})

# put binaries in ./bin
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)

# perform out of source build in build directory
add_subdirectory(src build)
