#!/usr/bin/perl

use strict;
use warnings;

use IO::File;
use Getopt::Long;

my $snp_file;
my $lq_output;
my $min_mapping_quality = 40;
my $min_somatic_score = 40;
my $out_file;
my $help;

my $opt_result;

$opt_result = GetOptions(
    'snp-file=s' => \$snp_file,
    'lq-output=s' => \$lq_output,
    'min-mapping-quality=i' => \$min_mapping_quality,
    'min-somatic-score=i' => \$min_somatic_score,
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
    warn "You must provide a file to be filtered\n";
    die help_text();
}

unless(-s $snp_file) {
    die "Can not find valid bam-somaticsniper output file: $snp_file\n";
}

$out_file ||= $snp_file . '.hc';
my $out_fh = IO::File->new($out_file,"w") or die "Unable to open $out_file for writing: $!\n";
my $lq_out_fh = undef;
if(defined($lq_output)) {
    $lq_out_fh = IO::File->new($lq_output,"w") or die "Unable to open $lq_output for writing: $!\n";
}
my $snp_fh = IO::File->new($snp_file) or die "Unable to open $snp_file for reading: $1\n";

my $sniper_vcf = 0; #attempt to auto-detect whether or not the user is providing VCF or classic input

while (my $snp = $snp_fh->getline) {
    if($snp =~ /^##fileformat=VCF/) {   #not checking version. Note that this also doesn't attempt to verify that this is a sniper output file.
        $sniper_vcf = 1;
    }
    if($snp =~ /^#/) {   #pass through VCF comments etc
        print $out_fh $snp;
        next;
    }

    chomp $snp;
    my ($mean_tumor_mapq,$somatic_score);

    my @fields = split("\t", $snp);
    if ($sniper_vcf) {
        my ($ref, $alts, $format,$tumor_sample) = @fields[3,4,8,10];
        my @tumor_fields = split /:/, $tumor_sample;
        my $index = 0;
        my %format_keys = map { $_ => $tumor_fields[$index++] } split /:/,$format;
        #these are in order ACGT
        my @alleles = ($ref, split /,/, $alts);
        my %gt_alleles = map {$_ => 1} split /\//, $format_keys{GT};
        my @used_alleles;
        for my $allele_index (keys %gt_alleles) {
            push @used_alleles, $alleles[$allele_index];
        }
        @used_alleles = sort @used_alleles;
        my $aindex = 0;
        my %mapq_for_allele = map { $used_alleles[$aindex++] => $_ } split /,/, $format_keys{AMQ};
        delete $mapq_for_allele{$ref};
        $mean_tumor_mapq = join(",", values %mapq_for_allele);  #could be one for each alt. This is kind of inefficient but shouldn't matter too much
        $somatic_score = $format_keys{SSC};
    } else {
        ($mean_tumor_mapq, $somatic_score) = @fields[18,5];
    }
    
    my $pass;
    for my $mapq (split /,/, $mean_tumor_mapq) {
        $pass ||= ($mapq >= $min_mapping_quality);
    }
    $pass &&= ($somatic_score >= $min_somatic_score);
    if($pass) {
        print $out_fh $snp,"\n";
    }
    else {
        print $lq_out_fh $snp,"\n" if($lq_out_fh);
    }
}
$snp_fh->close;
$out_fh->close;
if(defined($lq_output)){
    $lq_out_fh->close;
}


sub help_text {
    return <<HELP;
highconfidence - High confidence filtering for SomaticSniper

SYNOPSIS
highconfidence [options]

OPTIONS
--snp-file              the input bam-somaticsniper output file (requires v1.0.0 or greater output)
--lq-output             this is an optional place to stick sn(p|v)s which have failed to pass this filter
--min-mapping-quality   min mapping quality of the reads supporting the variant in the tumor, default $min_mapping_quality 
--min-somatic-score     minimum somatic score, default $min_somatic_score
--out-file              snp output file after filter
--help                  this message

DESCRIPTION
This program will filter bam-somaticsniper output into a highconfidence set based on mapping quality and somatic score

AUTHORS
Dave Larson     Extraction from standard TGI framework and modification specifically for bam-somaticsniper

SUPPORT
For user support please mail genome-dev\@genome.wustl.edu.
HELP
}
