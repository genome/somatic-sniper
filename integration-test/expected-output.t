#!/usr/bin/env perl

use warnings;
use strict;

use Test::More tests => 6;
use File::Basename qw/dirname/;
use File::Temp;

use_ok('Genome') or plan skip_all => "Couldn't 'use Genome', skipping tests";

my $rsb = Genome::Model::Build::ImportedReferenceSequence->get(name => 'NCBI-human-build36');
ok($rsb, 'got reference sequence build');

my $dir = (dirname($0) || ".");
my $exe = "$dir/../bin/bam-somaticsniper";
my $tumor_bam = "$dir/data/tumor.bam";
my $normal_bam = "$dir/data/normal.bam";
my $expected_file = "$dir/data/expected.snps";

my $tmpfh = File::Temp->new;
my $output_file = $tmpfh->filename;

my $ref_fasta = $rsb->fasta_file;
ok(-s $ref_fasta, 'reference fasta file exists');

ok(-s $exe, "sniper executable exists at $dir")
    or die "did you build the project?";

my $rv = system("$exe -f $ref_fasta $tumor_bam $normal_bam $output_file");
is($rv & 127, 0, 'sniper returned success');

my $diff = Genome::Sys->diff_file_vs_file($output_file, $expected_file);
ok(!$diff, 'output matched expected result')
    or diag("diff results:\n" . $diff);
