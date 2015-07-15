cmake_minimum_required(VERSION 2.8)

set(TEST_LIBS sniper)
include_directories(${GTEST_INCLUDE_DIRS})

set(TEST_SOURCES
    TestAlleleUtil.cpp
    TestDqStats.cpp
    )

add_unit_tests(SniperUnitTests ${TEST_SOURCES})
