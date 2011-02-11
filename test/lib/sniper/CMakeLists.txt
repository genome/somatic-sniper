cmake_minimum_required(VERSION 2.8)

set(TEST_LIBS sniper)
include_directories(${GTEST_INCLUDE_DIRS})

def_test(MeanQualities mean_qualities.c)
