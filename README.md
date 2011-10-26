# Somatic Sniper

p. The purpose of this program is to identify single nucleotide positions that
are different between tumor and normal (or, in theory, any two bam files).
It takes a tumor bam and a normal bam and compares the two to determine
the differences. Complete documentation is available at the project
[web site](http://gmt.genome.wustl.edu/somatic-sniper/).

## Build Dependencies

* git
* cmake 2.8+ ([cmake.org](http://cmake.org))
* samtools 0.1.6 ([sourceforge download page](http://sourceforge.net/projects/samtools/files/samtools/0.1.6/))

## Build Instructions

### Build Samtools

* Download and extract samtools 0.1.6 to a directory of your choosing. 
* Run make
* Set an environment variable SAMTOOLS_ROOT to point this directory (e.g., export SAMTOOLS_ROOT=`pwd`).

### Clone the Somatic Sniper repository
* Recursively clone the git repository

    git clone --recursive git://github.com/genome/somatic-sniper.git

### Compile Somatic Sniper
* Somatic Sniper does not support in source builds. Create a new build directory, enter it, and run:

    cmake /path/to/somatic-sniper/repo
    make

* Binaries can be found in the bin/ subdirectory of your build directory
