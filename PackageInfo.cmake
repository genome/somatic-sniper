cmake_minimum_required(VERSION 2.8)

# .deb packaging
set(ARCH "i686")
if(${CMAKE_SIZEOF_VOID_P} MATCHES 8)
    set(ARCH "x86_64")
endif ()

# The format of the description field is a short summary line followed by a
# longer paragraph indented by a single space on each line
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY
"Tool for detecting somatic mutations between tumor and normal.
 The purpose of this program is to identify single nucleotide
 positions that are different between tumor and normal (or, in
 theory, any two bam files). It takes a tumor bam and a normal
 bam and compares the two, outputting the differences.  It uses
 the genotype likelihood model of MAQ (as implemented in Samtools)
 and then calculates the probability that the tumor and normal
 genotypes are different. This probability is reported as a
 somatic score. The somatic score is the Phred-scaled probability
 (between 0 to 255) that the Tumor and Normal genotypes are not
 different where 0 means there is no probability that the
 genotypes are different and 255 means there is a probability of 
 1-10^(255/-10) that the genotypes are different between tumor and
 normal. This is consistent with how the SAM format reports such
 probabilities.")

set(CPACK_PACKAGE_NAME "somatic-sniper${EXE_VERSION_SUFFIX}")
set(CPACK_PACKAGE_VENDOR "wugc")
set(CPACK_PACKAGE_VERSION ${FULL_VERSION})
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "Dave Larson <dlarson@genome.wustl.edu>")
set(CPACK_SYSTEM_NAME "Linux-${ARCH}")
set(CPACK_TOPLEVEL_TAG "Linux-${ARCH}")
set(CPACK_DEBIAN_PACKAGE_SECTION science)
set(CPACK_DEBIAN_PACKAGE_PRIORITY optional)
set(CPACK_DEBIAN_PACKAGE_REPLACES "")
set(CPACK_DEBIAN_PACKAGE_DEPENDS "libc6 (>= 2.4), libgcc1 (>= 1:4.1.1-21), libstdc++6 (>= 4.2.1-4)")
if (CMAKE_BUILD_TYPE MATCHES package)
    set(CPACK_GENERATOR "DEB")
else(CMAKE_BUILD_TYPE MATCHES package)
    set(CPACK_GENERATOR "TGZ")
endif(CMAKE_BUILD_TYPE MATCHES package)

configure_file(debian/postinst.in debian/postinst @ONLY)
configure_file(debian/prerm.in debian/prerm @ONLY)
set(CPACK_DEBIAN_PACKAGE_CONTROL_EXTRA "debian/postinst;debian/prerm")

include(CPack)
