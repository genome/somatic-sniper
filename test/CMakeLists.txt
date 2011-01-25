cmake_minimum_required(VERSION 2.6)

set(TEST_LIBS z m bam)
include_directories(${GTEST_INCLUDE_DIRS})

def_test(MeanQualities mean_qualities.c)
