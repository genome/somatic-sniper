#include "sniper/allele_util.h"

#include <gtest/gtest.h>
#include <algorithm>
#include <utility>

using namespace std;

namespace {
    enum Allele {
        A = 1,
        C = 2,
        G = 4,
        T = 8
    };

    template<typename T1, typename T2>
    pair<int, int> mkpair(T1 a, T2 b) {
        return make_pair(int(a), int(b));
    }
}
    
TEST(AlleleUtil, count_alleles) {
    ASSERT_EQ(0, count_alleles(0));
    ASSERT_EQ(1, count_alleles(A));
    ASSERT_EQ(1, count_alleles(C));
    ASSERT_EQ(2, count_alleles(A|C));
    ASSERT_EQ(1, count_alleles(G));
    ASSERT_EQ(2, count_alleles(A|G));
    ASSERT_EQ(2, count_alleles(C|G));
    ASSERT_EQ(3, count_alleles(A|C|G));
    ASSERT_EQ(1, count_alleles(T));
    ASSERT_EQ(2, count_alleles(A|T));
    ASSERT_EQ(2, count_alleles(C|T));
    ASSERT_EQ(3, count_alleles(A|C|T));
    ASSERT_EQ(2, count_alleles(G|T));
    ASSERT_EQ(3, count_alleles(A|G|T));
    ASSERT_EQ(3, count_alleles(A|C|T));
    ASSERT_EQ(4, count_alleles(A|C|T|G));
}

TEST(AlleleUtil, is_loh) {
    // 1, 2, 4, and 8 are single alleles, so LOH can't really happen
    for (int i = 0; i < 4; ++i) {
        int value = 1 << i;
        for (int j = 1; j <= 8; ++j)
            ASSERT_EQ(0, is_loh(j, value));
    }

    // these are all the possible ways that LOH can happen with 2/3 alleles
    // (we do not concern ourselves with N until later)
    vector< pair<int,int> > loh_pairs;
    loh_pairs.push_back(mkpair(A, A|C));
    loh_pairs.push_back(mkpair(C, A|C));

    loh_pairs.push_back(mkpair(A, A|G));
    loh_pairs.push_back(mkpair(G, A|G));

    loh_pairs.push_back(mkpair(A, A|T));
    loh_pairs.push_back(mkpair(T, A|T));

    loh_pairs.push_back(mkpair(C, C|G));
    loh_pairs.push_back(mkpair(G, C|G));

    loh_pairs.push_back(mkpair(C, C|T));
    loh_pairs.push_back(mkpair(T, C|T));

    loh_pairs.push_back(mkpair(G, G|T));
    loh_pairs.push_back(mkpair(T, G|T));

    loh_pairs.push_back(mkpair(A, A|C|G));
    loh_pairs.push_back(mkpair(C, A|C|G));
    loh_pairs.push_back(mkpair(G, A|C|G));
    loh_pairs.push_back(mkpair(A|C, A|C|G));
    loh_pairs.push_back(mkpair(A|G, A|C|G));
    loh_pairs.push_back(mkpair(C|G, A|C|G));

    loh_pairs.push_back(mkpair(A, A|C|T));
    loh_pairs.push_back(mkpair(C, A|C|T));
    loh_pairs.push_back(mkpair(T, A|C|T));
    loh_pairs.push_back(mkpair(A|C, A|C|T));
    loh_pairs.push_back(mkpair(A|T, A|C|T));
    loh_pairs.push_back(mkpair(C|T, A|C|T));

    loh_pairs.push_back(mkpair(A, A|G|T));
    loh_pairs.push_back(mkpair(G, A|G|T));
    loh_pairs.push_back(mkpair(T, A|G|T));
    loh_pairs.push_back(mkpair(A|G, A|G|T));
    loh_pairs.push_back(mkpair(A|T, A|G|T));
    loh_pairs.push_back(mkpair(G|T, A|G|T));

    loh_pairs.push_back(mkpair(C, C|G|T));
    loh_pairs.push_back(mkpair(G, C|G|T));
    loh_pairs.push_back(mkpair(T, C|G|T));
    loh_pairs.push_back(mkpair(C|G, C|G|T));
    loh_pairs.push_back(mkpair(C|T, C|G|T));
    loh_pairs.push_back(mkpair(G|T, C|G|T));

    for (int orig = 1; orig < 15; ++orig) {
        for (int mut = 1; mut < 15; ++mut) {
            pair<int,int> p = make_pair(mut, orig);
            int expected = find(loh_pairs.begin(), loh_pairs.end(), p) != loh_pairs.end();
            ASSERT_EQ(expected, is_loh(int(p.first), int(p.second)));
        }
    }

    // deal with N here
    for (int i = 1; i < 15; ++i)
        ASSERT_EQ(1, is_loh(i, A|C|G|T));
}

TEST(AlleleUtil, should_filter_as_loh) {
    // these are all the possible ways that LOH can happen with 2/3 alleles
    // (we do not concern ourselves with N until later)
    int ref_base = A;

    ASSERT_TRUE(should_filter_as_loh(ref_base, A, A|G));
    ASSERT_TRUE(should_filter_as_loh(ref_base, G, A|G));
    ASSERT_TRUE(should_filter_as_loh(ref_base, G, C|G));
    ASSERT_TRUE(should_filter_as_loh(ref_base, C, C|G));
    ASSERT_FALSE(is_loh(A|G, G));
    ASSERT_TRUE(is_loh(G, A|G));
    ASSERT_TRUE(G != ref_base);
    ASSERT_TRUE(should_filter_as_loh(ref_base, A|G, G));

    ASSERT_FALSE(should_filter_as_loh(ref_base, A|G, A));
    ASSERT_FALSE(should_filter_as_loh(ref_base, A, A));
    ASSERT_FALSE(should_filter_as_loh(ref_base, G, A));
    ASSERT_FALSE(should_filter_as_loh(ref_base, G, G));
    ASSERT_FALSE(should_filter_as_loh(ref_base, A, G));
    ASSERT_FALSE(should_filter_as_loh(ref_base, A|G, A|G));
}
