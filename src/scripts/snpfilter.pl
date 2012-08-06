#!/usr/bin/perl

use strict;
use warnings;

use IO::File;
use Getopt::Long;

my %iub_as_string = (
    A => 'AA',
    C => 'CC',
    G => 'GG',
    T => 'TT',
    M => 'AC',
    K => 'GT',
    Y => 'CT',
    R => 'AG',
    W => 'AT',
    S => 'CG',
    D => 'AGT',
    B => 'CGT',
    H => 'ACT',
    V => 'ACG',
    N => 'ACGT',
);

my $snp_file;
my $lq_output;
my $min_mapping_quality = 40;
my $min_cns_qual = 20;
my $min_read_depth = 3;
my $max_read_depth = 100_000_000;
my $snp_win_size = 10;
my $max_snp_per_win = 2;
my $min_snp_qual = 20;
my $out_file;
my $indel_file;
my $indel_win_size = 10;
my $min_indel_score = 50;
my $tumor_variant_only = 0;
my $include_loh = 0;
my $help;

my $opt_result;

$opt_result = GetOptions(
    'snp-file=s' => \$snp_file,
    'lq-output=s' => \$lq_output,
    'min-mapping-quality=i' => \$min_mapping_quality,
    'min-cns-qual=i' => \$min_cns_qual,
    'min-read-depth=i' => \$min_read_depth,
    'max-read-depth=i' => \$max_read_depth,
    'snp-win-size=i' => \$snp_win_size,
    'max-snp-per-win=i' => \$max_snp_per_win,
    'min-snp-qual=i' => \$min_snp_qual,
    'out-file=s' => \$out_file,
    'indel-file=s' => \$indel_file,
    'indel-win-size=i' => \$indel_win_size,
    'min-indel-score=i' => \$min_indel_score,
    'tumor-variant-only!' => \$tumor_variant_only,
    'include-loh!' => \$include_loh,
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
    die "Can not find valid SAM snp file: $snp_file\n";
}

my %indel_filter;

if ($indel_file) {
    my $indel_fh = IO::File->new($indel_file) or die "Unable to open $indel_file: $!\n";

    while (my $indel = $indel_fh->getline) {
        my ($chr, $pos, $id, $indel_seq, $indel_score) = $indel =~ /^(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+\S+\s+(\S+)\s+/;
        next if $id ne '*' or $indel_seq eq '*/*' or $indel_score < $min_indel_score;
        #map{$indel_filter{$chr, $_}= 1}($pos - $indel_win_size .. $pos + $indel_win_size);
        $indel_filter{$chr, $pos} = 1;
    }
    $indel_fh->close;
}

my @snps = ();
my $last_chr = '';

$out_file ||= $snp_file . '.SNPfilter';
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
    my ($chr, $pos, $id, $ref, $var, $cns_qual, $rd_depth, $map_qual, $snp_qual, $tumor_gt, $normal_var, $somatic_status);
    if ($sniper_vcf) {
        my ($format, @samples);
        ($chr, $pos, $id, $ref, $var, undef , undef, undef, $format, @samples) =split("\t", $snp);
        my @tumor_fields = split /:/, $samples[1];
        my $index=0;
        my %format_keys = map { $_ => $tumor_fields[$index++] } split /:/,$format;
        $cns_qual = $format_keys{GQ};
        $snp_qual = $format_keys{VAQ};
        $map_qual = $format_keys{MQ};
        $rd_depth = $format_keys{DP};
        $tumor_gt = $format_keys{GT};
        $somatic_status = $format_keys{SS};
    } else {
        ($chr, $pos, $ref, $var, $normal_var, undef, $cns_qual, $snp_qual, $map_qual, undef, undef, undef, $rd_depth) = split("\t", $snp);
    }
    my $test = 0;
    for my $range_pos ($pos - $indel_win_size .. $pos + $indel_win_size) {
        if ($indel_filter{$chr, $range_pos}) {  
            $test++;
            last;
        }
    }
    if ($test) {
        if($lq_output){
            print $lq_out_fh $snp,"\n";
        }
        next;
    }
    next if $indel_filter{$chr,$pos};

    my $pass;
    $pass = 1 if $map_qual >= $min_mapping_quality and $rd_depth >= $min_read_depth and $rd_depth <= $max_read_depth;
    $pass = 0 unless $cns_qual >= $min_cns_qual || $snp_qual >= $min_snp_qual;

    if($tumor_variant_only && ( (defined($tumor_gt) && $tumor_gt eq '0/0') || ( !defined($tumor_gt) && $var eq $ref) )) {
        $pass = 0;
    }
    
    if(!$include_loh && ( (defined $somatic_status && $somatic_status eq '3') || (!defined $somatic_status && is_loh($var, $normal_var)))) {
        $pass = 0;
    }

    unless( $pass ) {
        if($lq_output){
            print $lq_out_fh $snp,"\n";
        }
        next;
    }

    if ($chr ne $last_chr) {
        map{$out_fh->print($_->{line}) if $_->{pass}}@snps;
        if(defined($lq_output)){
            map{$lq_out_fh->print($_->{line}) unless $_->{pass}}@snps;
        }
        @snps = ();       #reset
        $last_chr = $chr; #reset
    }

    push @snps, {
        line => $snp."\n",
        pos  => $pos,
        pass => 1,
    };

    if ($#snps == $max_snp_per_win) {
        if ($snps[$#snps]->{pos} - $snps[0]->{pos} < $snp_win_size) {
            map{$_->{pass} = 0}@snps;
        }
        if ($snps[0]->{pass}) {
            $out_fh->print($snps[0]->{line});
        }
        else {
            if(defined($lq_output)){
                $lq_out_fh->print($snps[0]->{line});
            }
        }
        shift @snps; # keep the size of @snps, moving the window snp by snp, check the snp density in a window for all snps.
    }
}
map{$out_fh->print($_->{line}) if $_->{pass}}@snps;
if(defined($lq_output)){
    map{$lq_out_fh->print($_->{line}) unless $_->{pass}}@snps;
}

$snp_fh->close;
$out_fh->close;
if(defined($lq_output)){
    $lq_out_fh->close;
}


sub is_loh {
    my ($tumor, $normal) = @_;
    if($normal =~ /[MKYRWS]/ && (index($iub_as_string{$normal},$tumor) > -1)) {
        return 1;
    }
    else {
        return 0;
    }
}

sub help_text {
    return <<HELP;
snpfilter - Basic filtering for SomaticSniper

SYNOPSIS
snpfilter [options] [file ...]

OPTIONS
--snp-file              the input bam-somaticsniper output file (requires v1.0.0 or greater output)
--lq-output             this is an optional place to stick sn(p|v)s which have failed to pass this filter
--min-mapping-quality   min mapping quality of the reads covering the SNP, default $min_mapping_quality 
--min-cns-qual          minimum consensus quality, default $min_cns_qual
--min-read-depth        minimum read depth to call a SNP, default $min_read_depth
--max-read-depth        maximum read depth to call a SNP, default $max_read_depth
--snp-win-size          window size for filtering dense SNPs, default $snp_win_size
--max-snp-per-win       maximum number of SNPs in a sized window, default $max_snp_per_win
--min-snp-qual          check minimum snp quality if consensus qual is lower than min_cns_qual, default $min_snp_qual
--out-file              snp output file after filter
--indel-file            path of samtools *pileup* format indel file to be used as a filter to screen out snps close to indel
--indel-win-size        window size of indel position in which SNPs should be filtered out, default $indel_win_size
--min-indel-score       minimum samtools indel score, default $min_indel_score
--tumor-variant-only    whether or not to pass homozygous ref calls in the tumor (off by default)
--include-loh           whether or not to pass sites likely to be loss of heterozygosity (on by default)
--help                  this message

DESCRIPTION
This program will filter bam-somaticsniper output with some basic filters inspired by maq.pl SNPfilter

AUTHORS
Dave Larson     Extraction from standard TGI framework and modification specifically for bam-somaticsniper
Feiyu Du        Original code

SUPPORT
For user support please mail genome-dev\@genome.wustl.edu.
HELP
}
