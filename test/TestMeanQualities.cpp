#include "mean_qualities.c"

#include <gtest/gtest.h>
#include <bam.h>
#include <stdio.h>

#include <string>

using namespace std;

namespace {
    string get_temp_value(FILE* fh) {
        string rv;
        char buf[4096];
        int br = 0;
        fseek(fh, 0, SEEK_SET);
        while((br = fread(buf, 1, sizeof(buf), fh)) > 0) {
            buf[br] = 0;
            rv.append(buf);
        }
        return rv; 
    }
}

TEST(MeanQualities, compute) {
    
}

TEST(MeanQualities, print) {
    uint32_t qual[4] = { 1,2,3,4 };
    FILE* fh = tmpfile();
    print_mean_quality_values(fh, 1, qual);
    string result = get_temp_value(fh);
    ASSERT_EQ("1", result);
    fclose(fh);


    fh = tmpfile();
    print_mean_quality_values(fh, 6, qual);
    result = get_temp_value(fh);
    ASSERT_EQ("2,3", result);
    fclose(fh);


    fh = tmpfile();
    print_mean_quality_values(fh, 3, qual);
    result = get_temp_value(fh);
    ASSERT_EQ("1,2", result);
    fclose(fh);

    fh = tmpfile();
    print_mean_quality_values(fh, 0, qual);
    result = get_temp_value(fh);
    ASSERT_EQ("0", result);
    fclose(fh);

}

