# SomaticSniper

   The purpose of this program is to identify single nucleotide positions that are different between tumor and normal (or in theory, any two bam files). It takes a tumor bam and a normal bam and compares the two to determine the differences. Complete documentation is available at the project [web site](http://gmt.genome.wustl.edu/somatic-sniper/) or in the included [pdf](https://github.com/genome/somatic-sniper/blob/master/docs/sniper_manual.pdf).

## Build Dependencies

* git
* cmake 2.8+ ([cmake.org](http://cmake.org))

## Build Instructions

### Build dependencies

* For APT-based systems (Debian, Ubuntu), install the following packages:

```
sudo apt-get install build-essential git-core cmake zlib1g-dev libncurses-dev
```

* For RPM-based systems (Fedora, CentOS, RHEL), install the following packages instead:

```
sudo yum groupinstall "Development tools" 
sudo yum install zlib-devel ncurses-devel cmake
```

Note that for some RPM based systems (like RHEL), you will need to install cmake 2.8 or greater yourself as the packaged version is much older.

### Clone the SomaticSniper repository

Recursively clone the git repository

```
git clone git://github.com/genome/somatic-sniper.git
```

### Build SomaticSniper

SomaticSniper does not support in-source builds. So create a subdirectory, enter it, build, and run tests:

```
mkdir somatic-sniper/build
cd somatic-sniper/build
cmake ../
make deps
make -j
make test
```

The binary `bam-somaticsniper` can then be found under `somatic-sniper/build/bin`. If you have administrative rights, then run `sudo make install` to install the tool for all users under `/usr/bin`.

## FAQ

### I get lots of compile errors indicating that files are missing. How do I fix this?

SomaticSniper requires that it be linked to an old version of samtools (v0.1.6). This typically happens because you have attempted to link to a newer version. As of version 1.0.3, SomaticSniper will take care of this during the compilation process.

### I get errors from cmake about missing modules. How do I fix this?

As of commit 09ef624e5bb275e0fd62396a14a878711e746cb9 or version 1.0.4, this should no longer be an issue and tarballs from github should function as intended. In earlier versions, SomaticSniper contained a git submodule called build-common. This submodeule contains helper modules for cmake. If you downloaded the source as a tarball from github or forgot to do a recursive clone using git, then you will not have this submodule and will see cmake errors. If you are using git, we recommend you go back and use the --recursive option when cloning the SomaticSniper repository. If you cannot use git, follow the instructions below to remedy the situation.

1. Download the build-common module separately [here](https://github.com/genome/build-common/tarball/master)
2. Extract that tarball and rename the directory it creates to 'build-common'.
3. Replace the empty build-common subdirectory in the sniper directory with directory you just created.
4. Resume following the build instructions.

## User Support

Please first search [Biostar](http://www.biostars.org) and then ask a question there if needed.  We automatically monitor [Biostar](http://www.biostars.org) for questions related to our tools.
