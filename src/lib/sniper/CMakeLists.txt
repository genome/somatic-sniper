cmake_minimum_required(VERSION 2.8)

set(SOURCES
    mean_qualities.c
    mean_qualities.h
    sniper_maqcns.c
    sniper_maqcns.h
    sniper_pileup.c
    somatic_sniper.c
    somatic_sniper.h
)

add_library(sniper ${SOURCES})
target_link_libraries(sniper z m bam)
