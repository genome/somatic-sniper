cmake_minimum_required(VERSION 2.8)

set(SOURCES
    dqstats.c
    dqstats.h
    output_bed.c
    output_bed.h
    output_classic.c
    output_classic.h
    output_format.c
    output_format.h
    output_vcf.c
    output_vcf.h
    sniper_maqcns.c
    sniper_maqcns.h
    sniper_pileup.c
    somatic_sniper.c
    somatic_sniper.h
)

add_library(sniper ${SOURCES})
target_link_libraries(sniper m ${Samtools_LIBRARIES})
