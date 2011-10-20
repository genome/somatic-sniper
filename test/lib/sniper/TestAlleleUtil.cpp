#include "sniper/allele_util.h"

#include <gtest/gtest.h>
#include <stdio.h>

#include <string>

using namespace std;

TEST(AlleleUtil, count_alleles) {
    ASSERT_EQ(0, count_alleles(0));
    ASSERT_EQ(1, count_alleles(1));
    ASSERT_EQ(1, count_alleles(2));
    ASSERT_EQ(2, count_alleles(3));
    ASSERT_EQ(1, count_alleles(4));
    ASSERT_EQ(2, count_alleles(5));
    ASSERT_EQ(2, count_alleles(6));
    ASSERT_EQ(3, count_alleles(7));
    ASSERT_EQ(1, count_alleles(8));
}

TEST(AlleleUtil, is_loh) {
    ASSERT_EQ(0, is_loh(1, 1));
    ASSERT_EQ(0, is_loh(3, 2));
    ASSERT_EQ(0, is_loh(4, 2));
    ASSERT_EQ(0, is_loh(2, 4));
    ASSERT_EQ(1, is_loh(2, 3));
}
