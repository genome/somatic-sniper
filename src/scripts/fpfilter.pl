#!/gsc/bin/perl

use warnings;
use strict;

use Getopt::Long;
use IO::File;




## Define filtering parameters ##

my $min_read_pos = 0.10;
my $max_read_pos = 1 - $min_read_pos;
my $min_var_freq = 0.05;
my $min_var_count = 4;

my $min_strandedness = 0.01;
my $max_strandedness = 1 - $min_strandedness;

my $max_mm_qualsum_diff = 50;
my $max_mapqual_diff = 30;
my $max_readlen_diff = 25;
my $min_var_dist_3 = 0.20;
my $max_var_mm_qualsum = 100;


## Parse arguments ##

my $output_basename;
my $verbose = 0;

my $snp_file;
my $readcount_file;
my $help;


my $opt_result;

$opt_result = GetOptions(
    'snp-file=s' => \$snp_file,
    'readcount-file=s' => \$readcount_file,
    'output-basename=s'   => \$output_basename,
    'verbose=s'   => \$verbose,
    'min-read-pos=f' => \$min_read_pos,
    'max-read-pos=f' => \$max_read_pos,
    'min-var-freq=f' => \$min_var_freq,
    'min-var-count=f' => \$min_var_count,
    'min-strandedness=f' => \$min_strandedness,
    'max-strandedness=f' => \$max_strandedness,
    'max-mm-qualsum-diff=f' => \$max_mm_qualsum_diff,
    'max-mapqual-diff=f' => \$max_mapqual_diff,
    'max-readlen-diff=f' => \$max_readlen_diff,
    'min-var-dist-3=f' => \$min_var_dist_3,
    'max_var_mm_qualsum=f' => \$max_var_mm_qualsum,
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

unless($readcount_file) {
    warn "You must provide a bam-readcount output file to use for filtering\n";
    die help_text();
}

unless(-s $readcount_file) {
    die "Can not find valid bam-readcount output file: $readcount_file\n";
}

$output_basename ||= $snp_file;

my %stats = ();
$stats{'num_variants'} = $stats{'num_fail_pos'} = $stats{'num_fail_strand'} = $stats{'num_fail_varcount'} = $stats{'num_fail_varfreq'} = $stats{'num_fail_mmqs'} = $stats{'num_fail_var_mmqs'} = $stats{'num_fail_mapqual'} = $stats{'num_fail_readlen'} = $stats{'num_fail_dist3'} = $stats{'num_pass_filter'} = $stats{'num_no_readcounts'} = 0;

## Load the read counts ##

my %readcounts_by_position = ();

my $rc_input = IO::File->new($readcount_file) or die "Unable to open $readcount_file for reading: $!\n";
my $lineCounter = 0;

while (<$rc_input>)
{
    chomp;
    my $line = $_;
    $lineCounter++;
    my ($chrom, $position) = split(/\t/, $line);
    $readcounts_by_position{"$chrom\t$position"} = $line;
}

close($rc_input);


## Open the output files ##

my $pass_fh = IO::File->new("$output_basename.fp_pass","w") or die "Can't open output file $output_basename.fp_pass: $!\n";
my $fail_fh = IO::File->new("$output_basename.fp_fail","w") or die "Can't open output file $output_basename.fp_fail: $!\n";

## Parse the input file ##


my $input = new IO::File->new($snp_file) or die "Unable to open input file $snp_file: $!\n";
$lineCounter = 0;
my $sniper_vcf = 0;

while (my $line = $input->getline)
{
    if($line =~ /^##fileformat=VCF/) {   #not checking version. Note that this also doesn't attempt to verify that this is a sniper output file.
        $sniper_vcf = 1;
    }
    if($line =~ /^#/) {   #pass through VCF comments etc
        print $pass_fh $line;
        next;
    }

    chomp $line;
    $lineCounter++;
    my @fields = split("\t", $line);
    my ($chrom, $position, $ref, $var);
    if($sniper_vcf) {
        my ($alt, $format, $tumor_sample);
        ($chrom, $position, $ref,$alt,$format, $tumor_sample) = @fields[0,1,3,4,8,10];
        my @tumor_fields = split /:/, $tumor_sample;
        my $index = 0;
        my %format_keys = map { $_ => $tumor_fields[$index++] } split /:/,$format;
        #these are in order ACGT
        my @alleles = ($ref, split /,/, $alt);
        my %gt_alleles = map {$_ => 1} grep { $_ > 0 } split /\//, $format_keys{GT};
        my @used_alleles;
        for my $allele_index (keys %gt_alleles) {
            push @used_alleles, $alleles[$allele_index];
        }
        ($var) = sort @used_alleles; #follow existing convention of fp filter using alphabetical order to choose a single base on triallelic sites
    }
    else {
        ($chrom, $position, $ref, $var) = @fields[0,1,2,3];
    }
    $ref = uc($ref);
    $var = uc($var);

    if(!($var =~ /[ACGT]/)) {
        $var = iupac_to_base($ref, $var);
    }

    $stats{'num_variants'}++;

    if($readcounts_by_position{"$chrom\t$position"})
    {
        my $readcounts = $readcounts_by_position{"$chrom\t$position"};
        my $ref_result = read_counts_by_allele($readcounts, $ref);
        my $var_result = read_counts_by_allele($readcounts, $var);

        if($ref_result && $var_result)
        {
            ## Parse out the bam-readcounts details for each allele. The fields should be: ##
            #num_reads : avg_mapqual : avg_basequal : avg_semq : reads_plus : reads_minus : avg_clip_read_pos : avg_mmqs : reads_q2 : avg_dist_to_q2 : avgRLclipped : avg_eff_3'_dist
            my ($ref_count, $ref_map_qual, $ref_base_qual, $ref_semq, $ref_plus, $ref_minus, $ref_pos, $ref_subs, $ref_mmqs, $ref_q2_reads, $ref_q2_dist, $ref_avg_rl, $ref_dist_3) = split(/\t/, $ref_result);
            my ($var_count, $var_map_qual, $var_base_qual, $var_semq, $var_plus, $var_minus, $var_pos, $var_subs, $var_mmqs, $var_q2_reads, $var_q2_dist, $var_avg_rl, $var_dist_3) = split(/\t/, $var_result);

            my $ref_strandedness = my $var_strandedness = 0.50;
            $ref_dist_3 = 0.5 if(!$ref_dist_3);

            ## Use conservative defaults if we can't get mismatch quality sums ##
            $ref_mmqs = 50 if(!$ref_mmqs);
            $var_mmqs = 0 if(!$var_mmqs);
            my $mismatch_qualsum_diff = $var_mmqs - $ref_mmqs;

            ## Determine map qual diff ##

            my $mapqual_diff = $ref_map_qual - $var_map_qual;


            ## Determine difference in average supporting read length ##

            my $readlen_diff = $ref_avg_rl - $var_avg_rl;


            ## Determine ref strandedness ##

            if(($ref_plus + $ref_minus) > 0) {
                $ref_strandedness = $ref_plus / ($ref_plus + $ref_minus);
                $ref_strandedness = sprintf("%.2f", $ref_strandedness);
            }

            ## Determine var strandedness ##

            if(($var_plus + $var_minus) > 0) {
                $var_strandedness = $var_plus / ($var_plus + $var_minus);
                $var_strandedness = sprintf("%.2f", $var_strandedness);
            }

            if($var_count && ($var_plus + $var_minus))
            {
                ## We must obtain variant read counts to proceed ##

                my $var_freq = $var_count / ($ref_count + $var_count);

                ## FAILURE 1: READ POSITION ##
                if(($var_pos < $min_read_pos)) { # || $var_pos > $max_read_pos)) {
                    print $fail_fh "$line\t$ref_pos\t$var_pos\t$ref_strandedness\t$var_strandedness\tReadPos<$min_read_pos\n";
                    print "$line\t$ref_pos\t$var_pos\t$ref_strandedness\t$var_strandedness\tReadPos<$min_read_pos\n"if ($verbose);
                    $stats{'num_fail_pos'}++;
                }

                ## FAILURE 2: Variant is strand-specific but reference is NOT strand-specific ##
                elsif(($var_strandedness < $min_strandedness || $var_strandedness > $max_strandedness) && ($ref_strandedness >= $min_strandedness && $ref_strandedness <= $max_strandedness)) {
                    ## Print failure to output file if desired ##
                    print $fail_fh "$line\t$ref_pos\t$var_pos\t$ref_strandedness\t$var_strandedness\tStrandedness: Ref=$ref_strandedness Var=$var_strandedness\n";
                    print "$line\t$ref_pos\t$var_pos\t$ref_strandedness\t$var_strandedness\tStrandedness: Ref=$ref_strandedness Var=$var_strandedness\n"if ($verbose);
                    $stats{'num_fail_strand'}++;
                }

                ## FAILURE : Variant allele count does not meet minimum ##
                elsif($var_count < $min_var_count) {
                    print $fail_fh "$line\t$ref_pos\t$var_pos\t$ref_strandedness\t$var_strandedness\tVarCount:$var_count\n";
                    print "$line\t$ref_pos\t$var_pos\t$ref_strandedness\t$var_strandedness\tVarCount:$var_count\n" if ($verbose);
                    $stats{'num_fail_varcount'}++;
                }

                ## FAILURE : Variant allele frequency does not meet minimum ##
                elsif($var_freq < $min_var_freq) {
                    print $fail_fh "$line\t$ref_pos\t$var_pos\t$ref_strandedness\t$var_strandedness\tVarFreq:$var_freq\n";
                    print "$line\t$ref_pos\t$var_pos\t$ref_strandedness\t$var_strandedness\tVarFreq:$var_freq\n" if ($verbose);
                    $stats{'num_fail_varfreq'}++;
                }

                ## FAILURE 3: Paralog filter for sites where variant allele mismatch-quality-sum is significantly higher than reference allele mmqs
                elsif($mismatch_qualsum_diff> $max_mm_qualsum_diff) {
                    ## Print failure to output file if desired ##
                    print $fail_fh "$line\t$ref_pos\t$var_pos\t$ref_strandedness\t$var_strandedness\tMismatchQualsum:$var_mmqs-$ref_mmqs=$mismatch_qualsum_diff\n";
                    print "$line\t$ref_pos\t$var_pos\t$ref_strandedness\t$var_strandedness\tMismatchQualsum:$var_mmqs-$ref_mmqs=$mismatch_qualsum_diff" if ($verbose);
                    $stats{'num_fail_mmqs'}++;
                }

                ## FAILURE 4: Mapping quality difference exceeds allowable maximum ##
                elsif($mapqual_diff > $max_mapqual_diff) {
                    print $fail_fh "$line\t$ref_pos\t$var_pos\t$ref_strandedness\t$var_strandedness\tMapQual:$ref_map_qual-$var_map_qual=$mapqual_diff\n";
                    print "$line\t$ref_pos\t$var_pos\t$ref_strandedness\t$var_strandedness\tMapQual:$ref_map_qual-$var_map_qual=$mapqual_diff" if ($verbose);
                    $stats{'num_fail_mapqual'}++;
                }

                ## FAILURE 5: Read length difference exceeds allowable maximum ##
                elsif($readlen_diff > $max_readlen_diff) {
                    print $fail_fh "$line\t$ref_pos\t$var_pos\t$ref_strandedness\t$var_strandedness\tReadLen:$ref_avg_rl-$var_avg_rl=$readlen_diff\n";
                    print "$line\t$ref_pos\t$var_pos\t$ref_strandedness\t$var_strandedness\tReadLen:$ref_avg_rl-$var_avg_rl=$readlen_diff" if ($verbose);
                    $stats{'num_fail_readlen'}++;
                }

                ## FAILURE 5: Read length difference exceeds allowable maximum ##
                elsif($var_dist_3 < $min_var_dist_3) {
                    print $fail_fh "$line\t$ref_pos\t$var_pos\t$ref_strandedness\t$var_strandedness\tVarDist3:$var_dist_3\n";
                    print "$line\t$ref_pos\t$var_pos\t$ref_strandedness\t$var_strandedness\tVarDist3:$var_dist_3\n" if ($verbose);
                    $stats{'num_fail_dist3'}++;
                }

                elsif($max_var_mm_qualsum && $var_mmqs > $max_var_mm_qualsum) {
                    print $fail_fh "$line\t$ref_pos\t$var_pos\t$ref_strandedness\t$var_strandedness\tVarMMQS: $var_mmqs > $max_var_mm_qualsum\n";
                    print "$line\t$ref_pos\t$var_pos\t$ref_strandedness\t$var_strandedness\tVarMMQS: $var_mmqs > $max_var_mm_qualsum\n" if ($verbose);
                    $stats{'num_fail_var_mmqs'}++;
                }

                ## SUCCESS: Pass Filter ##
                else {
                    $stats{'num_pass_filter'}++;
                    ## Print output, and append strandedness information ##
                    print $pass_fh "$line\n";
                    print "$line\tPASS\n" if($verbose);
                }

            }
        }
        else
        {
            $stats{'num_no_readcounts'}++;
            print $fail_fh "$line\tno_readcounts\n";
        }
    }
    else
    {
        $stats{'num_no_readcounts'}++;
        print $fail_fh "$line\tno_readcounts\n";
    }


}

close($input);

close($pass_fh);
close($fail_fh);

## Print filtering stats ##

print $stats{'num_variants'} . " variants\n";
print $stats{'num_no_readcounts'} . " failed to get readcounts for variant allele\n";
print $stats{'num_fail_pos'} . " had read position < $min_read_pos\n";
print $stats{'num_fail_strand'} . " had strandedness < $min_strandedness\n";
print $stats{'num_fail_varcount'} . " had var_count < $min_var_count\n";
print $stats{'num_fail_varfreq'} . " had var_freq < $min_var_freq\n";

print $stats{'num_fail_mmqs'} . " had mismatch qualsum difference > $max_mm_qualsum_diff\n";
print $stats{'num_fail_var_mmqs'} . " had variant MMQS > $max_var_mm_qualsum\n" if($stats{'num_fail_var_mmqs'});
print $stats{'num_fail_mapqual'} . " had mapping quality difference > $max_mapqual_diff\n";
print $stats{'num_fail_readlen'} . " had read length difference > $max_readlen_diff\n";
print $stats{'num_fail_dist3'} . " had var_distance_to_3' < $min_var_dist_3\n";

print $stats{'num_pass_filter'} . " passed the strand filter\n";


################################################################################

=head3	iupac_to_base

    Convert IUPAC ambiguity codes to variant bases


=cut


sub iupac_to_base {
    (my $allele1, my $allele2) = @_;

    return($allele2) if($allele2 eq "A" || $allele2 eq "C" || $allele2 eq "G" || $allele2 eq "T");

    if($allele2 eq "M") {
        return("C") if($allele1 eq "A");
        return("A") if($allele1 eq "C");
        return("A");    ## Default for triallelic variant
    } elsif($allele2 eq "R") {
        return("G") if($allele1 eq "A");
        return("A") if($allele1 eq "G");
        return("A");     ## Default for triallelic variant
    } elsif($allele2 eq "W") {
        return("T") if($allele1 eq "A");
        return("A") if($allele1 eq "T");
        return("A");    ## Default for triallelic variant
    } elsif($allele2 eq "S") {
        return("C") if($allele1 eq "G");
        return("G") if($allele1 eq "C");
        return("C");    ## Default for triallelic variant
    } elsif($allele2 eq "Y") {
        return("C") if($allele1 eq "T");
        return("T") if($allele1 eq "C");
        return("C");    ## Default for triallelic variant
    } elsif($allele2 eq "K") {
        return("G") if($allele1 eq "T");
        return("T") if($allele1 eq "G");
        return("G");    ## Default for triallelic variant
    }

    return($allele2);
}


################################################################################

=head3	read_counts_by_allele

    Retrieve relevant read counts for a certain allele 


=cut

sub read_counts_by_allele {
    (my $line, my $allele) = @_;

    my @lineContents = split(/\t/, $line);
    my $numContents = @lineContents;

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

            return($return_string);

        }
    }

    return("");
}


sub help_text {
    return <<HELP;
fpfilter - Advanced filtering for SomaticSniper

SYNOPSIS
fpfilter [options] [file ...]

OPTIONS
--snp-file              the input bam-somaticsniper output file (requires v1.0.0 or greater output)
--readcount-file        the output of bam-readcount for the snp-file
--help                  this message

DESCRIPTION
This program will filter bam-somaticsniper output with a variety of filters as detailed in the VarScan2 paper (http://www.ncbi.nlm.nih.gov/pubmed/22300766). It requires the bam-readcount utility (https://github.com/genome/bam-readcount). This is more convenient than the filtering described in the SomaticSniper paper, but more heuristic. Regardless, the principles are similar and we observe ~5% false negative rate with this filter.

This filter was calibrated on 100bp PE Illumina reads. It is likely to be overly stringent for longer reads and may be less effective on shorter reads.

AUTHORS
Dan Koboldt     Original code
Dave Larson     Modification specifically for bam-somaticsniper

SUPPORT
For user support please mail genome-dev\@genome.wustl.edu.
HELP
}
