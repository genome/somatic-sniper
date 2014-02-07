#!/usr/bin/env perl

use warnings;
use strict;
use Getopt::Long;
use IO::File;

## Define filtering parameters ##

# Minimums for variant allele fraction, and the number of variant-supporting reads
my $min_var_frac = 0.05;
my $min_var_count = 3;

# Minimum avg relative distance of variant from start/end of read
my $min_read_pos = 0.10;

# Minimum representation of variant allele on each strand
my $min_strandedness = 0.01;

# Maximum difference of mismatch quality sum between var/ref reads (paralog filter)
my $max_mmqs_diff = 50;

# Maximum mismatch quality sum of reference-supporting reads
my $max_var_mmqs = 100;

# Maximum difference of mapping quality between variant and reference reads
my $max_mapqual_diff = 30;

# Maximum difference of average supporting read length between variant and reference reads (paralog filter)
my $max_readlen_diff = 25;

# Minimum average distance to effective 3prime end of read (real end or Q2) for variant-supporting reads
my $min_var_dist_3 = 0.2;

## Parse arguments ##

my ( $output_file, $var_file, $readcount_file, $help );
my $opt_result = GetOptions(
    'var-file=s' => \$var_file,
    'readcount-file=s' => \$readcount_file,
    'output-file=s' => \$output_file,
    'min-read-pos=f' => \$min_read_pos,
    'min-var-frac=f' => \$min_var_frac,
    'min-var-count=f' => \$min_var_count,
    'min-strandedness=f' => \$min_strandedness,
    'max-mmqs-diff=f' => \$max_mmqs_diff,
    'max-mapqual-diff=f' => \$max_mapqual_diff,
    'max-readlen-diff=f' => \$max_readlen_diff,
    'min-var-dist-3=f' => \$min_var_dist_3,
    'max-var-mmqs=f' => \$max_var_mmqs,
    'help' => \$help,
) or die help_text();

# If help text was explicitly requested, print it to STDOUT rather than STDERR
if( $help ) {
    print help_text();
    exit 0;
}

# Make sure input files are properly defined and non-empty
unless( $var_file and $readcount_file ) {
    warn "var-file and readcount-file are required arguments!\n";
    die help_text();
}
die "Variant file not found, or is empty: $var_file\n" unless( -s $var_file );
die "Readcount file not found, or is empty: $readcount_file\n" unless( -s $readcount_file );

# If output file is undefined, use the input filename as a prefix
$output_file = "$var_file.fpfilter" unless( defined $output_file );

## Load the read counts into a hash for quick lookup ##

my %readcounts_by_position = ();
my $rc_fh = IO::File->new( $readcount_file ) or die "Can't open file $readcount_file: $!\n";
while( my $line = $rc_fh->getline ) {
    chomp( $line );
    my ( $chrom, $position ) = split( /\t/, $line );
    $readcounts_by_position{"$chrom\t$position"} = $line;
}
$rc_fh->close;

## Parse the input file and write pass/fail status to output file ##

my $input_fh = new IO::File->new( $var_file ) or die "Can't open file $var_file: $!\n";

# Open the output file for writing, and write a header line with column names
my $output_fh = IO::File->new( $output_file, "w" ) or die "Can't open file $output_file: $!\n";
$output_fh->print( "#CHROM\tPOS\tREF\tVAR\tDEPTH\tRAF\tVAF\tFILTER\tFILTER_DETAILS\n" );

# Initialize all variant fail/pass counters to zero
my %stats = map{($_,0)} qw( num_variants num_fail_pos num_fail_strand num_fail_varcount
    num_fail_varfrac num_fail_mmqs num_fail_var_mmqs num_fail_mapqual num_fail_readlen
    num_fail_dist3 num_no_readcounts num_pass_filter );

my $input_is_vcf = undef;
while( my $line = $input_fh->getline ) {

    # Skip blank lines
    next if( $line =~ m/^\s*$/ );

    # If the first non-blank line in the file is a VCF version tag, then remember that this as a VCF
    unless( defined $input_is_vcf ) {
        $input_is_vcf = 0;
        $input_is_vcf = 1 if( $line =~ /^##fileformat=VCF/ );
    }

    # Skip comment lines
    next if( $line =~ /^#/ );

    chomp( $line );
    my @fields = split( "\t", $line );

    # We'll only need the chromosome name, genomic position, and ref/var alleles from input
    my ( $chrom, $position, $ref, $var );
    if( $input_is_vcf ) {
        # Follow existing convention of alphabetical order to choose from alternate alleles
        ( $chrom, $position, undef, $ref, my $alt ) = @fields[0..4];
        ( $var ) = sort split( /,/, $alt );
    }
    else {
        ( $chrom, $position, $ref, $var ) = @fields;
    }
    $ref = uc( $ref );
    $var = uc( $var );

    # If the variant allele isn't all ACGTs then it likely contains IUPAC codes that need conversion
    if( $var !~ /^[ACGT]+$/ ) {
        $var = join( "", map{iupac_to_base( $ref, $_ )} split( //, $var ));
    }

    # Accumulate the total number of variants, for the final summarized report
    $stats{'num_variants'}++;

    # Proceed only if readcounts are available for this position
    if( $readcounts_by_position{"$chrom\t$position"} ) {
        my $readcounts = $readcounts_by_position{"$chrom\t$position"};
        my $ref_result = read_counts_by_allele( $readcounts, $ref );
        my $var_result = read_counts_by_allele( $readcounts, $var );

        # Proceed only if readcounts are available for reference and variant alleles
        if( $ref_result && $var_result ) {
            # Parse out the bam-readcounts details for each allele. The fields should be in this order:
            # total_reads : num_reads : avg_mapqual : avg_basequal : avg_semq : reads_plus : reads_minus :
            # avg_clip_read_pos : avg_mmqs : reads_q2 : avg_dist_to_q2 : avgRLclipped : avg_eff_3'_dist
            my ( $total_depth, $ref_count, $ref_map_qual, $ref_base_qual, $ref_semq, $ref_plus, $ref_minus, $ref_pos,
                $ref_subs, $ref_mmqs, $ref_q2_reads, $ref_q2_dist, $ref_avg_rl, $ref_dist_3 ) = split( /\t/, $ref_result );
            my ( undef, $var_count, $var_map_qual, $var_base_qual, $var_semq, $var_plus, $var_minus, $var_pos,
                $var_subs, $var_mmqs, $var_q2_reads, $var_q2_dist, $var_avg_rl, $var_dist_3 ) = split( /\t/, $var_result );

            # Set conservative defaults if read positions or mismatch quality sums are not available
            $ref_dist_3 = 0.5 unless( $ref_dist_3 );
            $ref_mmqs = 50 unless( $ref_mmqs );
            $var_mmqs = 0 unless( $var_mmqs );
            my $mismatch_qualsum_diff = $var_mmqs - $ref_mmqs;

            # Determine map qual diff between ref/var reads
            my $mapqual_diff = $ref_map_qual - $var_map_qual;

            # Determine difference in average supporting read length
            my $readlen_diff = $ref_avg_rl - $var_avg_rl;

            # Set max strandedness cutoff
            my $max_strandedness = 1 - $min_strandedness;

            # Set conservative default values for reference and variant strandedness
            my ( $ref_strandedness, $var_strandedness ) = ( 0.5, 0.5 );

            # Determine reference strandedness
            if(( $ref_plus + $ref_minus ) > 0 ) {
                $ref_strandedness = $ref_plus / ( $ref_plus + $ref_minus );
                $ref_strandedness = sprintf( "%.2f", $ref_strandedness );
            }

            # Determine variant strandedness
            if(( $var_plus + $var_minus ) > 0 ) {
                $var_strandedness = $var_plus / ( $var_plus + $var_minus );
                $var_strandedness = sprintf( "%.2f", $var_strandedness );
            }

            # Determine reference/variant allele fractions, and add as columns for output
            my $raf = sprintf( "%.4f", $total_depth ? $ref_count/$total_depth : 0 );
            my $vaf = sprintf( "%.4f", $total_depth ? $var_count/$total_depth : 0 );
            $line = "$chrom\t$position\t$ref\t$var\t$total_depth\t$raf\t$vaf";

            ## We must have non-zero variant read counts to proceed ##
            if( $var_count && ( $var_plus + $var_minus )) {

                ## FAILURE: READ POSITION ##
                if(($var_pos < $min_read_pos)) {
                    $line .= "\tReadPos\t$var_pos < $min_read_pos\n";
                    $stats{'num_fail_pos'}++;
                }

                ## FAILURE: Variant is strand-specific but reference is NOT strand-specific ##
                elsif(($var_strandedness < $min_strandedness || $var_strandedness > $max_strandedness) &&
                    ($ref_strandedness >= $min_strandedness && $ref_strandedness <= $max_strandedness)) {
                    ## Print failure to output file if desired ##
                    $line .= "\tStrandedness\tRef=$ref_strandedness,Var=$var_strandedness,MinMax=[$min_strandedness,$max_strandedness]\n";
                    $stats{'num_fail_strand'}++;
                }

                ## FAILURE: Variant allele count does not meet minimum ##
                elsif($var_count < $min_var_count) {
                    $line .= "\tVarCount\t$var_count < $min_var_count\n";
                    $stats{'num_fail_varcount'}++;
                }

                ## FAILURE: Variant allele fraction does not meet minimum ##
                elsif($vaf < $min_var_frac) {
                    $line .= "\tVarFrac\t$vaf < $min_var_frac\n";
                    $stats{'num_fail_varfrac'}++;
                }

                ## FAILURE: Paralog filter for sites where variant allele mismatch-quality-sum is significantly higher than the reference allele MMQS
                elsif( $mismatch_qualsum_diff > $max_mmqs_diff ) {
                    ## Print failure to output file if desired ##
                    $line .= "\tMismatchQualsum\t$var_mmqs-$ref_mmqs=$mismatch_qualsum_diff > $max_mmqs_diff\n";
                    $stats{'num_fail_mmqs'}++;
                }

                ## FAILURE: Mapping quality difference exceeds allowable maximum ##
                elsif($mapqual_diff > $max_mapqual_diff) {
                    $line .= "\tMapQual\t$ref_map_qual-$var_map_qual=$mapqual_diff > $max_mapqual_diff\n";
                    $stats{'num_fail_mapqual'}++;
                }

                ## FAILURE: Read length difference exceeds allowable maximum ##
                elsif( $readlen_diff > $max_readlen_diff ) {
                    $line .= "\tReadLen\t$ref_avg_rl-$var_avg_rl=$readlen_diff > $max_readlen_diff\n";
                    $stats{'num_fail_readlen'}++;
                }

                ## FAILURE: Avg distance from 3' ends of reads is lower than allowed minimum ##
                elsif( $var_dist_3 < $min_var_dist_3 ) {
                    $line .= "\tVarDist3\t$var_dist_3 < $min_var_dist_3\n";
                    $stats{'num_fail_dist3'}++;
                }

                ## FAILURE: Mismatch quality sum indicative of errors from misalignment ##
                elsif( $var_mmqs > $max_var_mmqs ) {
                    $line .= "\tVarMMQS\t$var_mmqs > $max_var_mmqs\n";
                    $stats{'num_fail_var_mmqs'}++;
                }

                ## SUCCESS: Passes all filters above ##
                else {
                    $line .= "\tPASS\t.\n";
                    $stats{'num_pass_filter'}++;
                }
            }
            else {
                $line .= "\tNoVariantReads\t.\n";
                $stats{'num_no_readcounts'}++;
            }
        }
        else {
            $line = "$chrom\t$position\t$ref\t$var\t0\t0\t0\tNoReadCounts\t.\n";
            $stats{'num_no_readcounts'}++;
        }
    }
    else {
        $line = "$chrom\t$position\t$ref\t$var\t0\t0\t0\tNoReadCounts\t.\n";
        $stats{'num_no_readcounts'}++;
    }

    # Print to output file along with filter status and additional information
    $output_fh->print( $line );
}

$input_fh->close;
$output_fh->close;

## Print filtering stats ##

print $stats{'num_variants'} . " variants\n";
print $stats{'num_no_readcounts'} . " failed to get readcounts for variant allele\n";
print $stats{'num_fail_pos'} . " were near the ends of the supporting reads (position < $min_read_pos)\n";
print $stats{'num_fail_strand'} . " had strandedness < $min_strandedness (most supporting reads are in the same direction)\n";
print $stats{'num_fail_varcount'} . " had var_count < $min_var_count (not enough supporting reads)\n";
print $stats{'num_fail_varfrac'} . " had var_frac < $min_var_frac (low-fraction variants are likely artifacts or from crosstalk between samples in the same lane)\n";
print $stats{'num_fail_mmqs'} . " had mismatch qualsum difference > $max_mmqs_diff (likely a result of paralogous misalignments)\n";
print $stats{'num_fail_var_mmqs'} . " had variant MMQS > $max_var_mmqs (likely a result of paralogous misalignments)\n" if($stats{'num_fail_var_mmqs'});
print $stats{'num_fail_mapqual'} . " had mapping quality difference > $max_mapqual_diff\n";
print $stats{'num_fail_readlen'} . " had read length difference > $max_readlen_diff\n";
print $stats{'num_fail_dist3'} . " had var_distance_to_3' < $min_var_dist_3 (illumina errors are more frequent at the 3' ends of reads)\n";
print $stats{'num_pass_filter'} . " passed the strand filter\n";

## iupac_to_base - Convert IUPAC ambiguity codes to variant bases ##

sub iupac_to_base {
    my ( $allele1, $allele2 ) = @_;

    return( $allele2 ) if( $allele2 eq "A" or $allele2 eq "C" or $allele2 eq "G" or $allele2 eq "T" );

    # Choose the most likely base-pair, or the default for triallelic variants
    if( $allele2 eq "M" ) {
        return( "C" ) if( $allele1 eq "A" );
        return( "A" ) if( $allele1 eq "C" );
        return( "A" );
    } elsif( $allele2 eq "R" ) {
        return( "G" ) if( $allele1 eq "A" );
        return( "A" ) if( $allele1 eq "G" );
        return( "A" );
    } elsif( $allele2 eq "W" ) {
        return( "T" ) if( $allele1 eq "A" );
        return( "A" ) if( $allele1 eq "T" );
        return( "A" );
    } elsif( $allele2 eq "S" ) {
        return( "C" ) if( $allele1 eq "G" );
        return( "G" ) if( $allele1 eq "C" );
        return( "C" );
    } elsif( $allele2 eq "Y" ) {
        return( "C" ) if( $allele1 eq "T" );
        return( "T" ) if( $allele1 eq "C" );
        return( "C" );
    } elsif( $allele2 eq "K" ) {
        return( "G" ) if( $allele1 eq "T" );
        return( "T" ) if( $allele1 eq "G" );
        return( "G" );
    }

    die "Failed to interpret variant allele: $allele2";
}

## read_counts_by_allele - Retrieve relevant read counts for a certain allele ##

sub read_counts_by_allele {
    my ( $line, $allele ) = @_;

    my @lineContents = split(/\t/, $line);
    my $numContents = @lineContents;
    my $total_depth = $lineContents[3];

    for(my $colCounter = 5; $colCounter < $numContents; $colCounter++) {
        my $this_allele = $lineContents[$colCounter];
        my @alleleContents = split(/\:/, $this_allele);
        if($alleleContents[0] eq $allele) {
            my $numAlleleContents = @alleleContents;

            return("") if($numAlleleContents < 8);

            my $return_string = "";
            my $return_sum = 0;
            for(my $printCounter = 1; $printCounter < $numAlleleContents; $printCounter++) {
                $return_sum += $alleleContents[$printCounter];
                $return_string .= "\t" if($return_string);
                $return_string .= $alleleContents[$printCounter];
            }

            return("$total_depth\t$return_string");

        }
    }

    return("");
}

## help_text - Returns usage syntax and documentation ##

sub help_text {
    return <<HELP;

fpfilter.pl - Advanced filtering for variant calls

SYNOPSIS
perl fpfilter.pl [options] [file ...]

OPTIONS
--var-file          List of variants in VCF, or tab-delimited list of "CHR POS REF VAR"
--readcount-file    The output of bam-readcount for the genomic loci in var-file
--output-file       Output file, tab-delimited list of variants with appended columns for filter status
--min-var-frac      Minimum variant allele fraction [$min_var_frac]
--min-var-count     Minimum number of variant-supporting reads [$min_var_count]
--min-read-pos      Minimum avg relative distance of variant from start/end of read [$min_read_pos]
--min-strandedness  Minimum representation of variant allele on each strand [$min_strandedness]
--max-mmqs-diff     Maximum difference of mismatch quality sum between var/ref reads (paralog filter) [$max_mmqs_diff]
--max-var-mmqs      Maximum mismatch quality sum of reference-supporting reads [$max_var_mmqs]
--max-mapqual-diff  Maximum difference of mapping quality between variant and reference reads [$max_mapqual_diff]
--max-readlen-diff  Maximum difference of average supporting read length between var/ref reads (paralog filter) [$max_readlen_diff]
--min-var-dist-3    Minimum avg distance to effective 3' end of read (real end or Q2) for variant-supporting reads [$min_var_dist_3]
--help              Show this message

DESCRIPTION
This program will filter a given list of variants output with a variety of filters as detailed in the VarScan2 paper (http://www.ncbi.nlm.nih.gov/pubmed/22300766). It relies on data generated by the bam-readcount utility (https://github.com/genome/bam-readcount).

This filter was calibrated on 100bp paired-end Illumina reads where we observed ~5% false negative rates. It is likely to be overly stringent for longer reads and may be less effective on shorter reads.

The output file is tab-delimited and contains a column with a filter tag that can be written back into a VCF using vcf-annotate.

AUTHORS
Dan Koboldt     Original code
Dave Larson     Modification specifically for bam-somaticsniper
Cyriac Kandoth  Fixed allele-fraction, improved docs, support for generic VCFs, and for MNPs/indels

SUPPORT
Please visit BioStars.org to see if your question was asked before. Else, post a new question, tag it with appropriate keywords, and the authors will find it.

HELP
}
