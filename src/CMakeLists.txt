cmake_minimum_required(VERSION 2.6)

set(EXECUTABLE_NAME bam-somaticsniper)

set(SOURCES
    sniper_maqcns.c
    sniper_maqcns.h
    sniper_pileup.c
    somatic_sniper.c
    somatic_sniper.h
)

add_executable(${EXECUTABLE_NAME} ${SOURCES})
target_link_libraries(${EXECUTABLE_NAME} z m bam)

install(TARGETS ${EXECUTABLE_NAME} DESTINATION bin/)
