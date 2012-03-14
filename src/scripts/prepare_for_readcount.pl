#!/usr/bin/perl

use strict;
use warnings;

use IO::File;
use Getopt::Long;

my $snp_file;
my $out_file;
my $help;

my $opt_result;

$opt_result = GetOptions(
    'snp-file=s' => \$snp_file,
    'out-file=s' => \$out_file,
    'help' => \$help,
);

unless($opt_result) {
    die help_text();
}

if($help) {
    print STDOUT help_text();
    exit 0;
}

unless($snp_file) {
    warn "You must provide a file to be converted\n";
    die help_text();
}

unless(-s $snp_file) {
    die "Can not find valid bam-somaticsniper output file: $snp_file\n";
}

$out_file ||= $snp_file . '.pos';
my $out_fh = IO::File->new($out_file,"w") or die "Unable to open $out_file for writing: $!\n";
my $snp_fh = IO::File->new($snp_file) or die "Unable to open $snp_file for reading: $1\n";

while (my $snp = $snp_fh->getline) {
    chomp $snp;
    my @fields = split /\t/, $snp;
    print $out_fh join("\t",@fields[0,1,1]),"\n";
}


sub help_text {
    return <<HELP;
prepare_for_readcount - Convert bam-somaticsniper output to an appropriate list of positions for bam-readcount's -l option

SYNOPSIS
prepare_for_readcount [options]

OPTIONS
--snp-file              the input bam-somaticsniper output file (requires v1.0.0 or greater output)
--out-file              snp output file after filter
--help                  this message

DESCRIPTION
This program extracts positions for use with bam-readcount's -l option.

AUTHORS
Dave Larson     Extraction from standard TGI framework and modification specifically for bam-somaticsniper

SUPPORT
For user support please mail genome-dev\@genome.wustl.edu.
HELP
}
