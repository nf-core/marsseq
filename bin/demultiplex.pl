#!/usr/bin/env perl
use strict;

if ($#ARGV == 0 and $ARGV[0] eq "--version") {
  print "v1.0";
  exit;
}

if ( $#ARGV < 11 ) {
    die(
        "usage: demultiplex.pl \
                amp_batch [str] \
                pool_batcode [str] \
                wells_cells.txt \
                gene_intervals_file.txt \
                spike-seq.txt \
                oligos.txt \
                trimmed.sam [file] \
                scdb_path [str output] \
                well_barcode_min_quality_thresh [int] \
                well_barcode_max_num_of_bad_quality_bp [int] \
                pool_barcode_min_quality_thresh [int] \
                pool_barcode_max_num_of_bad_quality_bp [int] \
                [downsampling rate]\n"
    );
}

my $amp_batch                              = $ARGV[0];
my $pool_barcode                           = $ARGV[1];
my $wells_cells_fn                         = $ARGV[2];
my $gene_intervals_fn                      = $ARGV[3];
my $spike_seq_fn                           = $ARGV[4];
my $oligos_fn                              = $ARGV[5];
my $trimmed_sam                            = $ARGV[6];
my $scdb_path                              = $ARGV[7];
my $well_barcode_min_quality_thresh        = $ARGV[8];
my $well_barcode_max_num_of_bad_quality_bp = $ARGV[9];
my $pool_barcode_min_quality_thresh        = $ARGV[10];
my $pool_barcode_max_num_of_bad_quality_bp = $ARGV[11];

my $downsampling_factor = 1;
if ( $#ARGV == 2 ) {
    $downsampling_factor = $ARGV[12];
}

my $genomic_bin_size            = 100;
my $hamming_thresh_cell_barcode = 1;

my $hamming_thresh_UMI = 1;

my $qual_string = "#$%&'()*+,^./0123456789:;<=>?\@ABCDEFGHIJ";

my %quality_hash = {};
for ( my $qi = 0 ; $qi < length($qual_string) ; $qi++ ) {
    $quality_hash{ substr( $qual_string, $qi, 1 ) } = $qi + 2;
}

# my $seq_batch = "";
# my $config_fn      = $scdb_path . "/config/config.txt";

#output_files
my $amp_batch_outdir = $amp_batch;
if ( $downsampling_factor > 1 ) {
    $amp_batch_outdir = $amp_batch . "_DS" . $downsampling_factor;
}
my $umitab_fn = $scdb_path . "/output/umi.tab/" . $amp_batch_outdir . ".txt";
my $offsetab_fn =
  $scdb_path . "/output/offset.tab/" . $amp_batch_outdir . ".txt";
my $singleofftab_fn =
  $scdb_path . "/output/singleton_offset.tab/" . $amp_batch_outdir . ".txt";

# if ( !-d $scdb_path . "/output/umi.tab/" ) {
#     mkdir $scdb_path . "/output/umi.tab/";
# }
# if ( !-d $scdb_path . "/output/offset.tab/" ) {
#     mkdir $scdb_path . "/output/offset.tab/";
# }
# if ( !-d $scdb_path . "/output/singleton_offset.tab/" ) {
#     mkdir $scdb_path . "/output/singleton_offset.tab/";
# }
# if ( !-d $scdb_path . "/output/QC/read_stats" ) {
#     mkdir $scdb_path . "/output/QC/read_stats";
# }
# if ( !-d $scdb_path . "/output/QC/read_stats_amp_batch" ) {
#     mkdir $scdb_path . "/output/QC/read_stats_amp_batch";
# }
# if ( !-d $scdb_path . "/output/QC/umi_stats" ) {
#     mkdir $scdb_path . "/output/QC/umi_stats";
# }
# if ( !-d $scdb_path . "/output/QC/noffsets_per_umi_distrib" ) {
#     mkdir $scdb_path . "/output/QC/noffsets_per_umi_distrib";
# }
# if ( !-d $scdb_path . "/output/QC/nreads_per_umi_distrib" ) {
#     mkdir $scdb_path . "/output/QC/nreads_per_umi_distrib";
# }
# if ( !-d $scdb_path . "/output/QC/umi_nuc_per_pos" ) {
#     mkdir $scdb_path . "/output/QC/umi_nuc_per_pos";
# }

# if ( !-d $scdb_path . "/_debug/" . $amp_batch_outdir ) {
#     mkdir $scdb_path . "/_debug/" . $amp_batch_outdir;
# }

my $debug_offsets_fn =
  $scdb_path . "/_debug/" . $amp_batch_outdir . "/offsets.txt";

#my $debug_offsets2_fn = $scdb_path."/processed_data/".$amp_batch."/offsets2.txt";
my $debug_umis_fn = $scdb_path . "/_debug/" . $amp_batch_outdir . "/UMIs.txt";

my $read_stats_fn =
  $scdb_path . "/output/QC/read_stats/" . $amp_batch_outdir . ".txt";
my $read_stats_amp_batch_fn =
  $scdb_path . "/output/QC/read_stats_amp_batch/" . $amp_batch_outdir . ".txt";
my $umi_stats_fn =
  $scdb_path . "/output/QC/umi_stats/" . $amp_batch_outdir . ".txt";
my $umi_nuc_per_pos_fn =
  $scdb_path . "/output/QC/umi_nuc_per_pos/" . $amp_batch_outdir . ".txt";
my $noffsets_per_umi_fn =
    $scdb_path
  . "/output/QC/noffsets_per_umi_distrib/"
  . $amp_batch_outdir . ".txt";
my $nreads_per_umi_fn =
    $scdb_path
  . "/output/QC/nreads_per_umi_distrib/"
  . $amp_batch_outdir . ".txt";

###########################################################################
# Read config file
# my %config_hash;
# open( CONFIG_FILE, $config_fn ) || die "ERROR: cannot open file $config_fn.\n";
# while (<CONFIG_FILE>) {
#     chomp;
#     my @row = split("=");
#     $config_hash{ $row[0] } = $row[1];
#     print $row[0] . "\t" . $row[1] . "\n";
# }
# my $well_barcode_min_quality_thresh = 27;
# if ( exists $config_hash{"well_barcode_min_quality_thresh"} ) {
#     $well_barcode_min_quality_thresh =
#       $config_hash{"well_barcode_min_quality_thresh"};
# }
# my $well_barcode_max_num_of_bad_quality_bp = 1;
# if ( exists $config_hash{"well_barcode_max_num_of_bad_quality_bp"} ) {
#     $well_barcode_max_num_of_bad_quality_bp =
#       $config_hash{"well_barcode_max_num_of_bad_quality_bp"};
# }
# my $pool_barcode_min_quality_thresh = 27;
# if ( exists $config_hash{"pool_barcode_min_quality_thresh"} ) {
#     $pool_barcode_min_quality_thresh =
#       $config_hash{"pool_barcode_min_quality_thresh"};
# }
# my $pool_barcode_max_num_of_bad_quality_bp = 1;
# if ( exists $config_hash{"pool_barcode_max_num_of_bad_quality_bp"} ) {
#     $pool_barcode_max_num_of_bad_quality_bp =
#       $config_hash{"pool_barcode_max_num_of_bad_quality_bp"};
# }

my $fdr = 0.25;

# if ( exists $config_hash{"fdr_offset_err"} ) {
#     $fdr = $config_hash{"fdr_offset_err"};
# }

# my $gene_intervals_fn = $scdb_path . "/" . $config_hash{"gene_intervals_file"};
# my $spike_seq_fn      = $scdb_path . "/" . $config_hash{"spike_seq_file"};
# my $oligos_fn         = $scdb_path . "/" . $config_hash{"oligos_file"};
my @oligos = ();
open( OLIGO_FILE, $oligos_fn ) || die "ERROR: cannot open file $oligos_fn.\n";
while (<OLIGO_FILE>) {
    chomp;
    my @line = split("=");
    push( @oligos, $line[0] );
}

print "fdr=$fdr\n";
print "well_barcode_min_quality_thresh=$well_barcode_min_quality_thresh\n";

###################################################################

sub map_to_gene {

    my $chr           = $_[0];
    my $coor          = $_[1];
    my $strand        = $_[2];
    my $gene_hash_ref = $_[3];
    my $bin           = int( $coor / $genomic_bin_size );
    my $key           = $chr . "_" . $strand . "_" . $bin;

    if ( !exists $gene_hash_ref->{$key} ) {
        return ("");
    }
    else {
        my $gene = $gene_hash_ref->{$key};
        return ($gene);
    }
}

##################################################################
sub fdr_thresh {

    my @pvalues = @{ $_[0] };

    my @sorted_pvalues = sort { $a <=> $b } @pvalues;
    my $n              = ( $#sorted_pvalues + 1 );
    for ( my $j = ( $#sorted_pvalues + 1 ) ; $j >= 1 ; $j-- ) {

        if ( $sorted_pvalues[ $j - 1 ] <= ( $fdr * $j / $n ) ) {
            return ( $sorted_pvalues[ $j - 1 ] );
        }
    }
    return (0);
}

sub logfact {
    return gammln( shift(@_) + 1.0 );
}

sub hypergeom {

    # There are m "bad" and n "good" balls in an urn.
    # Pick N of them. The probability of i or more successful selections:
    # (m!n!N!(m+n-N)!)/(i!(n-i)!(m+i-N)!(N-i)!(m+n)!)
    my ( $n, $m, $N, $i ) = @_;

    #   if ($m==0||$n==0){
    #	 return(1);
    #   }

    my $loghyp1 =
      logfact($m) + logfact($n) + logfact($N) + logfact( $m + $n - $N );
    my $loghyp2 =
      logfact($i) +
      logfact( $n - $i ) +
      logfact( $m + $i - $N ) +
      logfact( $N - $i ) +
      logfact( $m + $n );
    return exp( $loghyp1 - $loghyp2 );
}

sub gammln {
    my $xx  = shift;
    my @cof = (
        76.18009172947146,   -86.50532032941677,
        24.01409824083091,   -1.231739572450155,
        0.12086509738661e-2, -0.5395239384953e-5
    );
    my $y   = my $x = $xx;
    my $tmp = $x + 5.5;
    $tmp -= ( $x + .5 ) * log($tmp);
    my $ser = 1.000000000190015;
    for my $j ( 0 .. 5 ) {
        $ser += $cof[$j] / ++$y;
    }
    -$tmp + log( 2.5066282746310005 * $ser / $x );
}

sub hamming {
    my ( $s1, $s2, $max ) = @_;

    my ($dist) = 0;
    for ( my ($i) = 0 ; $i < length($s1) ; $i++ ) {
        if ( substr( $s1, $i, 1 ) ne substr( $s2, $i, 1 ) ) {
            $dist++;
            if ( $dist >= $max ) {
                return ($dist);
            }
        }
    }
    return ($dist);
}

############################################################################################################
#
#  filter_errors
#
#  This function makes filtering decisions for UMI and barcode errors and return p-values for offset errors.
#
#  Parameter:
#  $hash_in_ref -       ref for hash that stores all the input information of a gene
#                   keys: tab delimited (offset, well ID, barcode sequence, UMI sequence)
#                   values: number of reads
#  $hash_out -      ref for hash that  stores all the input+output information of a gene
#                   keys: tab delimited (offset, well ID, barcode sequence, UMI sequence, ...)
#  $cur_gene -      The name of the gene
#  $offsets_p1_ref -    ref for hash that stores the hypergeometic p-values of the "lonely offset" test.
#                   keys: tab delimited (gene,offset)
#                   values: raw p-value
#  $offsets_p2_ref -    ref for hash that stores the hypergeometic p-values of the "read poor offset" test.
#                   keys: tab delimited (gene,offset)
#                   values: raw p-value

sub filter_errors {
    my ( $hash_in_ref, $hash_out, $UMI_hash, $cur_gene, $offsets_p1_ref,
        $offsets_p2_ref, $hamming_thresh_UMI, $hamming_thresh_cell_barcode )
      = @_;
    my %filt_UMI_hash     = ();
    my %filt_barcode_hash = ();

    my (%UMI_off_map);
    my (%UMI_off_list);
    my (%UMI_count);

    my (%barcode_off_map);
    my (%barcode_off_list);
    my (%barcode_count);
    my (%barcode_to_wellid);

    my (%cell_offset_to_umis);
    my (%offset_to_mol);
    my (%mol_to_n_offsets);
    my (%offset_stats_hash);
    my %hash_in = %$hash_in_ref;

    # my %offsets_p1=%{$offsets_p1_ref};
    #my %offsets_p2=%{$offsets_p2_ref};

    my @data                                  = keys %hash_in;
    my ($n_total_offsets_with_single_read)    = 0;
    my ($n_total_offsets_with_multiple_reads) = 0;
    my ($tot);

    for ( my ($i) = 0 ; $i <= $#data ; $i++ ) {
        my ( $offset, $wellid, $barcode, $UMI ) = split( "\t", $data[$i] );
        $barcode_to_wellid{$barcode} = $wellid;
        if ( !exists( $UMI_off_map{$wellid} ) ) {
            $UMI_off_map{$wellid} = {};
        }
        $UMI_off_map{$wellid}->{"$UMI\t$offset"} = 1;
        if ( !exists( $UMI_off_list{$wellid}->{$UMI} ) ) {
            $UMI_off_list{$wellid}->{$UMI} = ();
        }
        push( @{ $UMI_off_list{$wellid}->{$UMI} }, $offset );
        if ( !exists $UMI_count{$wellid} ) {
            $UMI_count{$wellid} = {};
        }
        $UMI_count{$wellid}->{$UMI}++;

        if ( !exists( $barcode_off_map{$UMI} ) ) {
            $barcode_off_map{$UMI} = {};
        }
        $barcode_off_map{$UMI}->{"$barcode\t$offset"} = 1;
        if ( !exists( $barcode_off_list{$UMI}->{$barcode} ) ) {
            $barcode_off_list{$UMI}->{$barcode} = ();
        }
        push( @{ $barcode_off_list{$UMI}->{$barcode} }, $offset );
        if ( !exists $barcode_count{$UMI} ) {
            $barcode_count{$UMI} = {};
        }
        $barcode_count{$UMI}->{$barcode}++;

        my $reads = $hash_in{ $data[$i] };
        if ( !exists( $cell_offset_to_umis{$wellid} ) ) {
            $cell_offset_to_umis{$wellid} = {};
        }
        if ( !exists( $cell_offset_to_umis{$wellid}->{$offset} ) ) {
            $cell_offset_to_umis{$wellid}->{$offset} = {};
        }
        $cell_offset_to_umis{$wellid}->{$offset}->{$UMI} = 1;

        # Need to decide whether to use stats from empty wells
        #	if (exists(%$single_cells->{$wellid})){

        if ( !exists( $offset_to_mol{$offset} ) ) {
            $offset_to_mol{$offset} = {};
        }
        $offset_to_mol{$offset}->{ $wellid . "\t" . $UMI } = $reads;
        if ( $reads == 1 ) {
            $n_total_offsets_with_single_read++;
        }
        if ( $reads > 1 ) {
            $n_total_offsets_with_multiple_reads++;
        }
        $mol_to_n_offsets{ $wellid . "\t" . $UMI }++;

        #	}
        $tot++;
    }

    my @all_wellids = keys %UMI_count;
    for ( my $i = 0 ; $i <= $#all_wellids ; $i++ ) {
        my $wellid = $all_wellids[$i];
        my (@all_UMIs) = keys %{ $UMI_count{$wellid} };

        #	print STDERR "Total #UMIs=".$#all_UMIs."\n";
        my (@UMIs) =
          sort { $UMI_count{$wellid}->{$a} <=> $UMI_count{$wellid}->{$b} }
          @all_UMIs;
        for ( my ($i) = 0 ; $i <= $#UMIs ; $i++ ) {
            my ($UMI)    = $UMIs[$i];
            my ($i_offs) = $UMI_off_list{$wellid}->{$UMI};
            my ($filt)   = 0;
            for ( my ($j) = $i + 1 ; $j <= $#UMIs ; $j++ ) {
                my ($dom_UMI) = $UMIs[$j];
                if ( hamming( $UMI, $dom_UMI, $hamming_thresh_UMI + 1 ) >
                    $hamming_thresh_UMI )
                {
                    next;
                }
                my ($off);
                my ($indep) = 0;
                foreach $off ( @{$i_offs} ) {
                    if ( !exists( $UMI_off_map{$wellid}->{"$dom_UMI\t$off"} ) )
                    {
                        $indep = 1;
                        last;
                    }
                }
                if ( $indep == 0 ) {
                    $filt = 1;
                }
            }
            $filt_UMI_hash{ $wellid . "_" . $UMI } = $filt;
        }
    }

    my @all_UMIs = keys %barcode_count;
    for ( my $i = 0 ; $i <= $#all_UMIs ; $i++ ) {
        my $UMI = $all_UMIs[$i];
        my (@all_barcodes) = keys %{ $barcode_count{$UMI} };
        my (@barcodes) =
          sort { $barcode_count{$UMI}->{$a} <=> $barcode_count{$UMI}->{$b} }
          @all_barcodes;
        for ( my ($i) = 0 ; $i <= $#barcodes ; $i++ ) {
            my ($barcode) = $barcodes[$i];
            my ($i_offs)  = $barcode_off_list{$UMI}->{$barcode};
            my ($filt)    = 0;
            for ( my ($j) = $i + 1 ; $j <= $#barcodes ; $j++ ) {
                my ($dom_barcode) = $barcodes[$j];
                if (
                    hamming(
                        $barcode, $dom_barcode,
                        $hamming_thresh_cell_barcode + 1
                    ) > $hamming_thresh_cell_barcode |
                    $barcode_to_wellid{$barcode} eq
                    $barcode_to_wellid{$dom_barcode}
                  )
                {
                    next;
                }
                my ($off);
                my ($indep) = 0;
                foreach $off ( @{$i_offs} ) {
                    if (
                        !exists(
                            $barcode_off_map{$UMI}->{"$dom_barcode\t$off"}
                        )
                      )
                    {
                        $indep = 1;
                        last;
                    }
                }
                if ( $indep == 0 ) {
                    $filt = 1;
                }
            }
            $filt_barcode_hash{ $barcode . "_" . $UMI } = $filt;

        }

    }

    my ($n_total_lonely_offsets)       = 0;
    my ($n_total_offsets_with_friends) = 0;

    foreach my $mol ( keys %mol_to_n_offsets ) {
        if ( $mol_to_n_offsets{$mol} == 1 ) {
            $n_total_lonely_offsets++;
        }
        else {
            $n_total_offsets_with_friends += $mol_to_n_offsets{$mol};
        }
    }

    foreach my $offset ( keys %offset_to_mol ) {

        my ($n_current_lonely_offsets)       = 0;
        my ($n_current_offsets_with_friends) = 0;

        my ($n_current_offsets_with_single_read)   = 0;
        my ($n_current_offsets_with_multiple_read) = 0;
        my ( @v_n_reads_per_offset,
            @v_n_offsets_per_molecules_participating_in_this_offset, $i );

        foreach my $mol ( keys %{ $offset_to_mol{$offset} } ) {
            $v_n_reads_per_offset[$i] = $offset_to_mol{$offset}{$mol};
            $v_n_offsets_per_molecules_participating_in_this_offset[$i] =
              $mol_to_n_offsets{$mol};

            if ( $mol_to_n_offsets{$mol} == 1 ) {
                $n_current_lonely_offsets++;
            }
            else {
                $n_current_offsets_with_friends++;
            }

            if ( ${ $offset_to_mol{$offset} }{$mol} == 1 ) {
                $n_current_offsets_with_single_read++;
            }
            else {
                $n_current_offsets_with_multiple_read++;
            }

            $i++;
        }

        my ($n_molecules_with_this_offset) =
          scalar keys %{ $offset_to_mol{$offset} };

        my $hypergeom_p1 = hypergeom(
            $n_total_lonely_offsets,
            $n_total_offsets_with_friends,
            $n_current_lonely_offsets + $n_current_offsets_with_friends,
            $n_current_lonely_offsets
        );
        my $hypergeom_p2 = hypergeom(
            $n_total_offsets_with_single_read,
            $n_total_offsets_with_multiple_reads,
            $n_current_offsets_with_single_read +
              $n_current_offsets_with_multiple_read,
            $n_current_offsets_with_single_read
        );

        my $stats =
"$n_total_lonely_offsets\t$n_total_offsets_with_friends\t$n_current_lonely_offsets\t$n_current_offsets_with_friends\t$n_total_offsets_with_single_read\t$n_total_offsets_with_multiple_reads\t$n_current_offsets_with_single_read\t$n_current_offsets_with_multiple_read";

        #	print $cur_gene."\t".$offset."\t".$stats."\n";

        $offset_stats_hash{ $cur_gene . "\t" . $offset } = $stats;
        $offsets_p1_ref->{ $cur_gene . "\t" . $offset }  = $hypergeom_p1;
        $offsets_p2_ref->{ $cur_gene . "\t" . $offset }  = $hypergeom_p2;
    }

    foreach my $s ( keys %hash_in ) {
        my ( $offset, $wellid, $barcode, $UMI ) = split( "\t", $s );

        my $s2 =
            $s . "\t"
          . $hash_in{$s} . "\t"
          . $filt_UMI_hash{ $wellid . "_" . $UMI } . "\t"
          . $filt_barcode_hash{ $barcode . "_" . $UMI } . "\t"
          . $offsets_p1_ref->{ $cur_gene . "\t" . $offset } . "\t"
          . $offsets_p2_ref->{ $cur_gene . "\t" . $offset } . "\t"
          . $offset_stats_hash{ $cur_gene . "\t" . $offset };

        $hash_out->{$s2} = $hash_in{$s};
    }
}

#############################################################################

my $index;

my %hash1;
my %hash2;

my %number_of_cells;

my %UMI_hash;
my %umi_stats;
my %nuc_per_pos_stats;
my @well_list = ();

# read sample_index
# make barcode hash
# read amplification batch, map reads to genes, store in hashes for filtering
# write debug files
# write cell gene matrix.

my %cell_barcode_to_well_id_hash            = ();
my %extended_cell_barcode_to_well_id_hash   = ();
my %extended_cell_barcode_in_amp_batch_hash = ();
my %extended_pool_barcode_hash              = ();
my %column_num;
my %gene_hash;
my %read_counter;
my %oligo_counter;

###########################################################################
#Read gene intervals
my %binned_coordinate_to_gene_hash;
open( GENE_INTERVALS_FILE, $gene_intervals_fn )
  || die "ERROR: cannot open file $gene_intervals_fn.\n";
my $line = <GENE_INTERVALS_FILE>;
chomp $line;
my @column_names = split( "\t", $line );
for ( my $i = 0 ; $i <= $#column_names ; $i++ ) {
    $column_num{ $column_names[$i] } = $i;
}
while (<GENE_INTERVALS_FILE>) {
    chomp;
    my @row       = split("\t");
    my $chr       = $row[ $column_num{"chrom"} ];
    my $start     = $row[ $column_num{"start"} ];
    my $end       = $row[ $column_num{"end"} ];
    my $strand    = $row[ $column_num{"strand"} ];
    my $gene_name = $row[ $column_num{"gene_name"} ];
    $gene_hash{$gene_name} = 1;

    for (
        my $bin = int( $start / $genomic_bin_size ) ;
        $bin <= int( $end / $genomic_bin_size ) ;
        $bin++
      )
    {
        my $key = $chr . "_" . $strand . "_" . $bin;
        $binned_coordinate_to_gene_hash{$key} = $gene_name;
    }
}

my %spike_ins;
open( SPIKE_FILE, $spike_seq_fn )
  || die "ERROR: cannot open file $spike_seq_fn.\n";
<SPIKE_FILE>;
while (<SPIKE_FILE>) {
    my @l = split( "\t", $_ );
    $spike_ins{ $l[0] } = 1;
}

my @genes = sort keys %gene_hash;
for ( my $j = 0 ; $j <= $#genes ; $j++ ) {
    $hash1{ $genes[$j] } = {};
    $hash2{ $genes[$j] } = {};
}

###########################################################################
# Read amplification batch details
# open( AMP_BATCHES_FILE, $amp_batches_fn )
#   || die "ERROR: cannot open file $amp_batches_fn.\n";
# $line = <AMP_BATCHES_FILE>;
# chomp $line;
# @column_names = split( "\t", $line );
# for ( my $i = 0 ; $i <= $#column_names ; $i++ ) {
#     $column_num{ $column_names[$i] } = $i;
# }
# my $pool_barcode;
# my $protocol_version;
# my $seq_batch;

# my $cur_amp_batch = "";
# while (<AMP_BATCHES_FILE>) {
#     chomp;
#     my @row = split("\t");
#     $cur_amp_batch = $row[ $column_num{"Amp_batch_ID"} ];
#     if ( $cur_amp_batch eq $amp_batch ) {
#         $pool_barcode     = $row[ $column_num{"Pool_barcode"} ];
#         $protocol_version = $row[ $column_num{"Protocol_version_ID"} ];
#         $seq_batch        = $row[ $column_num{"Seq_batch_ID"} ];
#         last;
#     }
# }
# if ( $cur_amp_batch eq "" ) {
#     die("ERROR: $amp_batch does not exist!");
# }

###########################################################################
#Read barcodes

open( INDEX_FILE, $wells_cells_fn )
  || die "ERROR: cannot open file $wells_cells_fn.\n";

%column_num = ();
my @column_names = split( "\t", <INDEX_FILE> );
for ( my $i = 0 ; $i <= $#column_names ; $i++ ) {
    $column_num{ $column_names[$i] } = $i;
}

while (<INDEX_FILE>) {
    chomp;
    my @row           = split("\t");
    my $cur_amp_batch = $row[ $column_num{"Amp_batch_ID"} ];

  #  print $cur_amp_batch."\t".$amp_batch."\t".$column_num{"Amp_batch_ID"}."\n";
    if ( $cur_amp_batch eq $amp_batch ) {
        my $sample_index   = $row[ $column_num{"Well_ID"} ];
        my $sample_barcode = $row[ $column_num{"Cell_barcode"} ];
        $number_of_cells{$sample_index} =
          $row[ $column_num{"Number_of_cells"} ];
        push( @well_list, $sample_index );

        $extended_cell_barcode_in_amp_batch_hash{$sample_barcode} = 1;
        $cell_barcode_to_well_id_hash{$sample_barcode}          = $sample_index;
        $extended_cell_barcode_to_well_id_hash{$sample_barcode} = $sample_index;
        $read_counter{$sample_index}                            = {};
        $oligo_counter{$sample_index}                           = {};
    }
}

my @nucs = ( "A", "C", "G", "T" );
foreach my $barcode ( keys %cell_barcode_to_well_id_hash ) {
    for ( my $i = 0 ; $i < length($barcode) ; $i++ ) {
        if (   ( substr( $barcode, $i, 1 ) ne "_" )
            && ( substr( $barcode, $i, 1 ) ne "N" ) )
        {
            for ( my $to_nuci = 0 ; $to_nuci < 4 ; $to_nuci++ ) {
                my $barcode2 = $barcode;
                substr( $barcode2, $i, 1 ) = $nucs[$to_nuci];
                if ( exists $cell_barcode_to_well_id_hash{$barcode2} ) {
                    next;
                }
                if ( !exists $extended_cell_barcode_to_well_id_hash{$barcode2} )
                {
                    $extended_cell_barcode_to_well_id_hash{$barcode2} =
                      $cell_barcode_to_well_id_hash{$barcode};
                    $extended_cell_barcode_in_amp_batch_hash{$barcode2} =
                      $extended_cell_barcode_in_amp_batch_hash{$barcode};
                }
                else {
                    $extended_cell_barcode_to_well_id_hash{$barcode2}   = -1;
                    $extended_cell_barcode_in_amp_batch_hash{$barcode2} = -1;
                }

#		print $barcode."\t".$barcode2."\t".$cell_barcode_to_well_id_hash{$barcode}."\t".$extended_cell_barcode_to_well_id_hash{$barcode}."\t".$extended_cell_barcode_in_amp_batch_hash{$barcode2}."\n";
            }
        }
    }
}

my @nucs = ( "A", "C", "G", "T" );
for ( my $i = 0 ; $i < length($pool_barcode) ; $i++ ) {
    if (   ( substr( $pool_barcode, $i, 1 ) ne "_" )
        && ( substr( $pool_barcode, $i, 1 ) ne "N" ) )
    {
        for ( my $to_nuci = 0 ; $to_nuci < 4 ; $to_nuci++ ) {
            my $barcode2 = $pool_barcode;
            substr( $barcode2, $i, 1 ) = $nucs[$to_nuci];
            $extended_pool_barcode_hash{$barcode2} = 1;
        }
    }
}

#############################################################
# Read Reads
#

# my @fn = <$scdb_path/_trimmed_mapped_reads/$seq_batch/*.sam>;

# if ( $#fn <= 0 ) {
#     die("ERROR: Can't read $trimmed_sam file\n");
# }
# print "reading $#fn sam files\n";

my %umi_status;

my %reads_per_umi;
my %offsets_per_umi;

# Reading mapped reads and filtering reads that are not mapped to gene intervals
my @parsed_line;
my @parsed_header;
my $counter_reads_with_unknown_barcodes = 0;
my $counter_bad_quality_well_barcodes   = 0;
my $counter_bad_quality_pool_barcodes   = 0;
my $counter_bad_quality_UMI             = 0;
my $counter_amp_batch_reads             = 0;
my $counter_seq_batch_reads             = 0;

my $fh = select(STDOUT);
$| = 1;
select($fh);

# foreach my $fn (@fn) {
open( DATA, $trimmed_sam ) || die "ERROR: cannot open file $trimmed_sam.\n";
print "Reading $trimmed_sam ...";
my $ii = 0;
while (<DATA>) {
    chomp;
    my $line = $_;
    if ( $line =~ /^@/ ) {
        next;
    }
    $ii++;

    if ( $ii % $downsampling_factor != 0 ) {
        next;
    }

    @parsed_line = split( "\t", $line );
    my $barcode_info = $parsed_line[0];
    $barcode_info =~ s/^.*barcode=//;
    @parsed_header = split( "-", $barcode_info );

    my $cur_pool_barcode = $parsed_header[ $#parsed_header - 2 ];

    $counter_seq_batch_reads++;
    if ( $extended_pool_barcode_hash{$cur_pool_barcode} == 0 ) {
        next;
    }

    my $cur_pool_quality              = $parsed_header[ $#parsed_header - 5 ];
    my $n_bad_poolbarcode_quality_bps = 0;
    my @cur_poolbarcode_arr           = split( //, $cur_pool_quality );

    for ( my $qi = 0 ; $qi < $#cur_poolbarcode_arr ; $qi++ ) {
        my $bp_qual = ( 0 + $quality_hash{ $cur_poolbarcode_arr[$qi] } );
        if ( ($bp_qual) < $pool_barcode_min_quality_thresh ) {
            $n_bad_poolbarcode_quality_bps++;
        }
    }

    if ( $n_bad_poolbarcode_quality_bps >
        $pool_barcode_max_num_of_bad_quality_bp )
    {
        $counter_bad_quality_pool_barcodes++;
        next;
    }

    $counter_amp_batch_reads++;

    my $cur_wellbarcode_quality = $parsed_header[ $#parsed_header - 4 ];
    my $cur_UMI_quality         = $parsed_header[ $#parsed_header - 3 ];

    my $cell_barcode = $parsed_header[ $#parsed_header - 1 ];
    my $UMI          = $parsed_header[$#parsed_header];
    my $flag         = $parsed_line[1];
    my $mapq         = $parsed_line[4];

    my $strand = 1;
    my $well_id;

    #	my $well_barcode_min_bp_quality=1000;
    my $n_bad_quality_bps   = 0;
    my @cur_wellbarcode_arr = split( //, $cur_wellbarcode_quality );

    for ( my $qi = 0 ; $qi < $#cur_wellbarcode_arr ; $qi++ ) {
        my $bp_qual = ( 0 + $quality_hash{ $cur_wellbarcode_arr[$qi] } );
        if ( ($bp_qual) < $well_barcode_min_quality_thresh ) {
            $n_bad_quality_bps++;
        }
    }

    if ( $n_bad_quality_bps > $well_barcode_max_num_of_bad_quality_bp ) {
        $counter_bad_quality_well_barcodes++;
        next;
    }

    if ( $UMI =~ m/N/ ) {
        $counter_bad_quality_UMI++;
        next;
    }

    if ( $extended_cell_barcode_in_amp_batch_hash{$cell_barcode} < 1 ) {
        $counter_reads_with_unknown_barcodes++;
        next;
    }

    $well_id = $extended_cell_barcode_to_well_id_hash{$cell_barcode};
    if ( $well_id < 0 ) {
        next;
    }

    # At this point we know that this read has a valid pool and well barcode
    $read_counter{$well_id}->{"total"}++;
    my $oligo_est = $parsed_header[ $#parsed_header - 6 ];
    if ( $flag & 4 ) {
        if ( $oligo_est eq "NA" ) {
            $read_counter{$well_id}->{"unmapped"}++;
        }
        else {
            $oligo_counter{$well_id}->{$oligo_est}++;
        }
        next;
    }
## checking multi-mapping
    my $as = -1000;
    my $xs = -1000;
    for ( my $j = 9 ; $j <= $#parsed_line ; $j++ ) {

        if ( $parsed_line[$j] =~ /^AS:i:/ ) {
            ( $as = $parsed_line[$j] ) =~ s/^AS:i://;
        }
        if ( $parsed_line[$j] =~ /^XS:i:/ ) {
            ( $xs = $parsed_line[$j] ) =~ s/^XS:i://;
        }
    }
    if ( ( $xs >= $as ) & ( $as > -1000 ) ) {
        $read_counter{$well_id}->{"non_unique_mapping"}++;
        next;
    }

    if ( $mapq < 30 ) {
        $read_counter{$well_id}->{"lowmapq"}++;
        next;
    }

    if ( $flag & 16 ) {
        $strand = -1;
    }

    my $chr  = $parsed_line[2];
    my $coor = $parsed_line[3];
    my $gene =
      map_to_gene( $chr, $coor, $strand, \%binned_coordinate_to_gene_hash );

    if ( $gene eq "" ) {
        $read_counter{$well_id}->{"mapped_to_nongenic"}++;
        next;
    }

    my $pos = $chr . "_" . $strand . "_" . $coor;
    if ( !exists $UMI_hash{$UMI} ) {
        $UMI_hash{$UMI} = {};
    }
    $UMI_hash{$UMI}->{"$well_id\t$gene\t$pos\t$chr\t$strand\t$coor"} = 1;

    # incrementing read counts;
    $hash1{$gene}->{"$pos\t$well_id\t$cell_barcode\t$UMI"}++;

}
print "Read " . ( int( $ii / $downsampling_factor ) ) . " reads\n";
close(DATA);

# }

print "Finished reading sam files\n";

##############################################################################
#open(TMP,">$debug_offsets2_fn" )|| die "could not open debug file\n";
#my @UMIS=keys %UMI_hash;
#for (my $ii=0;$ii<=$#UMIS;$ii++){
#  my @offs= keys %{$UMI_hash{$UMIS[$ii]}};
#  for (my $jj=0;$jj<=$#offs;$jj++){
#	print TMP $UMIS[$ii]."\t".$offs[$jj]."\n";
#  }
#}
#close(TMP);
##############################################################################

# Now we will filter reads with UMI or barcode single bp error.
print "Filtering single bp erros\n";

open( OUT_DEBUG_OFFSETS, ">$debug_offsets_fn" )
  || die "ERROR: cannot open file $debug_offsets_fn.\n";
print OUT_DEBUG_OFFSETS
"gene_name\toffset\twell_id\tbarcode\tUMI\treads\tfiltered_UMI_err\tfiltered_barcode_err\tp_lonely_offset\tp_readpoor_offset\tn_total_lonely_offsets\tn_total_offsets_with_friends\tn_current_lonely_offsets\tn_current_offsets_with_friends\tn_total_offsets_with_single_read\tn_total_offsets_with_multiple_reads\tn_current_offsets_with_single_read\tn_current_offsets_with_multiple_read\n";

my @genes = sort keys %gene_hash;
my %offsets_p1;
my %offsets_p2;
for ( my $j = 0 ; $j <= $#genes ; $j++ ) {
    my $nrows = keys %{ $hash1{ $genes[$j] } };

    #print $genes[$j]."\t$j\t$nrows\n";
    filter_errors(
        \%{ $hash1{ $genes[$j] } }, \%{ $hash2{ $genes[$j] } },
        \%UMI_hash,                 $genes[$j],
        \%offsets_p1,               \%offsets_p2,
        $hamming_thresh_UMI,        $hamming_thresh_cell_barcode
    );
    my @data = keys %{ $hash2{ $genes[$j] } };
    for ( my $i = 0 ; $i <= $#data ; $i++ ) {
        print OUT_DEBUG_OFFSETS $genes[$j] . "\t" . $data[$i] . "\n";
    }
}

close(OUT_DEBUG_OFFSETS);

my @p1_vec = values %offsets_p1;
my @p2_vec = values %offsets_p2;
print $#p1_vec. "\n";
for ( my $a = 0 ; $a <= $#p1_vec ; $a++ ) {
    if ( $p1_vec[$a] eq "" ) {
        print $a;
    }
}
my $fdr_thresh1 = fdr_thresh( \@p1_vec );
my $fdr_thresh2 = fdr_thresh( \@p2_vec );
print "FDR_THRESHES:\t$fdr_thresh1\t$fdr_thresh2\n";

##############################################################################

# Writing umitab, offsetab, singleofftab
print "Writing output files\n";
open( UMITAB_FILE, ">$umitab_fn" )
  || die "ERROR: cannot open file $umitab_fn.\n";
open( OFFSETAB_FILE, ">$offsetab_fn" )
  || die "ERROR: cannot open file $offsetab_fn.\n";
open( SINGLEOFFTAB_FILE, ">$singleofftab_fn" )
  || die "ERROR: cannot open file $singleofftab_fn.\n";
print UMITAB_FILE $well_list[0];
print OFFSETAB_FILE $well_list[0];
print SINGLEOFFTAB_FILE $well_list[0];

for ( my $i = 1 ; $i <= $#well_list ; $i++ ) {
    print UMITAB_FILE "\t" . $well_list[$i];
    print OFFSETAB_FILE "\t" . $well_list[$i];
    print SINGLEOFFTAB_FILE "\t" . $well_list[$i];
}
print UMITAB_FILE "\n";
print OFFSETAB_FILE "\n";
print SINGLEOFFTAB_FILE "\n";

print $#genes. " genes\n";
for ( my $j = 0 ; $j <= $#genes ; $j++ ) {

    my $gene          = $genes[$j];
    my $spike_or_gene = "gene";
    if ( exists $spike_ins{$gene} ) {
        $spike_or_gene = "spike";
    }

    my @data = keys %{ $hash2{$gene} };

    my %well_to_UMIs    = {};
    my %well_to_UMIs_ok = {};

    my %wellgene_to_noffsets           = {};
    my %wellgene_to_nsingleton_offsets = {};
    for ( my $i = 0 ; $i <= $#data ; $i++ ) {

        # Iterating over all IVT products of the gene
        my ( $off, $wellid, $barcode, $UMI, $nreads, $filt1, $filt2, $p1, $p2 )
          = split( "\t", $data[$i] );
        my $well_gene     = $wellid . "\t" . $gene;
        my $well_gene_umi = $wellid . "\t" . $gene . "\t" . $UMI;
        if ( !exists $well_to_UMIs{$wellid} ) {
            $well_to_UMIs{$wellid} = {};
        }
        if ( !exists $well_to_UMIs_ok{$wellid} ) {
            $well_to_UMIs_ok{$wellid} = {};
        }
        $well_to_UMIs{$wellid}->{$UMI} = 1;

        # Offset filtering: filtering UMI if all its IVT products are filtered
        if ( $p1 < $fdr_thresh1 ) {
            $read_counter{$wellid}->{ $spike_or_gene . "\tfilt3" } += $nreads;
            if ( !exists $umi_status{$well_gene_umi} ) {
                $umi_status{$well_gene_umi} = -3;
            }
        }
        elsif ( $p2 < $fdr_thresh2 ) {
            $read_counter{$wellid}->{ $spike_or_gene . "\tfilt4" } += $nreads;
            if ( !exists $umi_status{$well_gene_umi} ) {
                $umi_status{$well_gene_umi} = -4;
            }
        }

        # UMI filtering
        elsif ( $filt1 == 1 ) {
            $read_counter{$wellid}->{ $spike_or_gene . "\tfilt1" } += $nreads;
            if ( $umi_status{$well_gene_umi} >= 0 ) {
                $umi_status{$well_gene_umi} = -1;
            }
        }

        # Well barcode filtering
        elsif ( $filt2 == 1 ) {
            $read_counter{$wellid}->{ $spike_or_gene . "\tfilt2" } += $nreads;
            if ( $umi_status{$well_gene_umi} >= 0 ) {
                $umi_status{$well_gene_umi} = -2;
            }
        }
        else {
            ##This (gene,well,UMI,OFFSET) is not going to be filtered
            if ( $spike_or_gene eq "spike" ) {
                $read_counter{$wellid}->{"spike_mapped"} += $nreads;
            }
            else {
                $read_counter{$wellid}->{"gene_mapped"} += $nreads;
            }

            $well_to_UMIs_ok{$wellid}->{$UMI} = 1;
            $wellgene_to_noffsets{$well_gene}++;
            $umi_status{$well_gene_umi} = 1;
            my @umiarr = split( "", $UMI );
            for ( my $ni = 0 ; $ni <= $#umiarr ; $ni++ ) {
                $nuc_per_pos_stats{ $ni . "\t" . $umiarr[$ni] }++;
            }

            $reads_per_umi{$well_gene_umi} += $nreads;
            $offsets_per_umi{$well_gene_umi}++;

            if ( $offsets_per_umi{$well_gene_umi} == 1 ) {
                $wellgene_to_nsingleton_offsets{$well_gene}++;
            }
            elsif ( $offsets_per_umi{$well_gene_umi} == 2 ) {
                $wellgene_to_nsingleton_offsets{$well_gene}--;
            }
        }
    }

    print UMITAB_FILE $gene;
    print OFFSETAB_FILE $gene;
    print SINGLEOFFTAB_FILE $gene;
    for ( my $i = 0 ; $i <= $#well_list ; $i++ ) {
        my $wellid  = $well_list[$i];
        my @UMIs    = keys %{ $well_to_UMIs{$wellid} };
        my @UMIs_ok = keys %{ $well_to_UMIs_ok{$wellid} };

        foreach my $UMI (@UMIs) {
            $umi_stats{ $wellid . "\t"
                  . $spike_or_gene . "\t"
                  . $umi_status{ $wellid . "\t" . $gene . "\t" . $UMI } }++;
        }

        my $number_of_UMIs = 1 + $#UMIs_ok;
        my $number_of_offsets =
          0 + $wellgene_to_noffsets{ $well_list[$i] . "\t" . $gene };
        my $number_of_singleton_offsets =
          0 + $wellgene_to_nsingleton_offsets{ $well_list[$i] . "\t" . $gene };

        print UMITAB_FILE "\t" . $number_of_UMIs;
        print OFFSETAB_FILE "\t" . $number_of_offsets;
        print SINGLEOFFTAB_FILE "\t" . $number_of_singleton_offsets;
    }
    print UMITAB_FILE "\n";
    print OFFSETAB_FILE "\n";
    print SINGLEOFFTAB_FILE "\n";
}
close(UMITAB_FILE);
close(OFFSETAB_FILE);
close(SINGLEOFFTAB_FILE);

# Writing stats

open( OUT_READ_STATS_AMP_BATCH, ">$read_stats_amp_batch_fn" )
  || die "ERROR: cannot open file $read_stats_amp_batch_fn.\n";
print OUT_READ_STATS_AMP_BATCH
"nreads_with_unknown_cell_barcode\tnreads_bad_quality_well_barcode\tnreads_bad_quality_pool_barcode\tnreads_bad_quality_UMI\tnreads_amp_batch\tnreads_seq_batch\n";

print OUT_READ_STATS_AMP_BATCH
"$counter_reads_with_unknown_barcodes\t$counter_bad_quality_pool_barcodes\t$counter_bad_quality_well_barcodes\t$counter_bad_quality_UMI\t$counter_amp_batch_reads\t$counter_seq_batch_reads\n";
close(OUT_READ_STATS_AMP_BATCH);

my $oligos_header = "";
for ( my $oi = 0 ; $oi <= $#oligos ; $oi++ ) {
    $oligos_header = $oligos_header . $oligos[$oi] . "\t";
}

open( OUT_READ_STATS_PER_CELL, ">$read_stats_fn" )
  || die "ERROR: cannot open file $read_stats_fn.\n";
print OUT_READ_STATS_PER_CELL
"well_id\tunmapped\tnon_unique_mapping\tlow_mapq\tmapped_to_nongenic\tspike_UMI_err\tspike_barcode_err\tspike_lonely_offset_err\tspike_readpoor_offset_err\tspike_mapped\tgene_UMI_err\tgene_barcode_err\tgene_lonely_offset_err\tgene_readpoor_offset_err\tgene_mapped\t"
  . $oligos_header
  . "total\n";

for ( my $j = 0 ; $j <= $#well_list ; $j++ ) {
    my $well_id = $well_list[$j];

    print OUT_READ_STATS_PER_CELL "$well_id\t";
    print OUT_READ_STATS_PER_CELL ( 0 + $read_counter{$well_id}->{"unmapped"} )
      . "\t";
    print OUT_READ_STATS_PER_CELL (
        0 + $read_counter{$well_id}->{"non_unique_mapping"} )
      . "\t";
    print OUT_READ_STATS_PER_CELL ( 0 + $read_counter{$well_id}->{"lowmapq"} )
      . "\t";
    print OUT_READ_STATS_PER_CELL (
        0 + $read_counter{$well_id}->{"mapped_to_nongenic"} )
      . "\t";
    print OUT_READ_STATS_PER_CELL (
        0 + $read_counter{$well_id}->{"spike\tfilt1"} )
      . "\t";
    print OUT_READ_STATS_PER_CELL (
        0 + $read_counter{$well_id}->{"spike\tfilt2"} )
      . "\t";
    print OUT_READ_STATS_PER_CELL (
        0 + $read_counter{$well_id}->{"spike\tfilt3"} )
      . "\t";
    print OUT_READ_STATS_PER_CELL (
        0 + $read_counter{$well_id}->{"spike\tfilt4"} )
      . "\t";
    print OUT_READ_STATS_PER_CELL (
        0 + $read_counter{$well_id}->{"spike_mapped"} )
      . "\t";
    print OUT_READ_STATS_PER_CELL (
        0 + $read_counter{$well_id}->{"gene\tfilt1"} )
      . "\t";
    print OUT_READ_STATS_PER_CELL (
        0 + $read_counter{$well_id}->{"gene\tfilt2"} )
      . "\t";
    print OUT_READ_STATS_PER_CELL (
        0 + $read_counter{$well_id}->{"gene\tfilt3"} )
      . "\t";
    print OUT_READ_STATS_PER_CELL (
        0 + $read_counter{$well_id}->{"gene\tfilt4"} )
      . "\t";
    print OUT_READ_STATS_PER_CELL (
        0 + $read_counter{$well_id}->{"gene_mapped"} )
      . "\t";

    for ( my $oi = 0 ; $oi <= $#oligos ; $oi++ ) {
        print OUT_READ_STATS_PER_CELL (
            0 + $oligo_counter{$well_id}->{ $oligos[$oi] } )
          . "\t";
    }
    print OUT_READ_STATS_PER_CELL ( 0 + $read_counter{$well_id}->{"total"} )
      . "\n";
}

close(OUT_READ_STATS_PER_CELL);

open( OUT_UMI_STATS_PER_CELL, ">$umi_stats_fn" )
  || die "ERROR: cannot open file $umi_stats_fn.\n";
print OUT_UMI_STATS_PER_CELL
"well_id\tgene_UMI_err\tgene_barcode_err\tgene_lonely_offset_err\tgene_readpoor_offset_err\tgene_ok\tspike_UMI_err\tspike_barcode_err\tspike_lonely_offset_err\tspike_readpoor_offset_err\tspike_ok\n";
for ( my $j = 0 ; $j <= $#well_list ; $j++ ) {
    my $well_id = $well_list[$j];
    print OUT_UMI_STATS_PER_CELL "$well_id\t"
      . ( 0 + $umi_stats{"$well_id\tgene\t-1"} ) . "\t"
      . ( 0 + $umi_stats{"$well_id\tgene\t-2"} ) . "\t"
      . ( 0 + $umi_stats{"$well_id\tgene\t-3"} ) . "\t"
      . ( 0 + $umi_stats{"$well_id\tgene\t-4"} ) . "\t"
      . ( 0 + $umi_stats{"$well_id\tgene\t1"} ) . "\t"
      . ( 0 + $umi_stats{"$well_id\tspike\t-1"} ) . "\t"
      . ( 0 + $umi_stats{"$well_id\tspike\t-2"} ) . "\t"
      . ( 0 + $umi_stats{"$well_id\tspike\t-3"} ) . "\t"
      . ( 0 + $umi_stats{"$well_id\tspike\t-4"} ) . "\t"
      . ( 0 + $umi_stats{"$well_id\tspike\t1"} ) . "\n";
}

close(OUT_UMI_STATS_PER_CELL);

open( OUT_NUC_PER_POS_STATS, ">$umi_nuc_per_pos_fn" )
  || die "ERROR: cannot open file $umi_nuc_per_pos_fn.\n";
my $UMI_size = ( keys %nuc_per_pos_stats ) / 4;
print OUT_NUC_PER_POS_STATS "position\tA\tC\tG\tT\n";
for ( my $ni = 0 ; $ni < $UMI_size ; $ni++ ) {
    print OUT_NUC_PER_POS_STATS $ni . "\t"
      . $nuc_per_pos_stats{ $ni . "\tA" } . "\t"
      . $nuc_per_pos_stats{ $ni . "\tC" } . "\t"
      . $nuc_per_pos_stats{ $ni . "\tG" } . "\t"
      . $nuc_per_pos_stats{ $ni . "\tT" } . "\n";
}
close(OUT_NUC_PER_POS_STATS);

my @umis                  = keys %umi_status;
my %noffsets_per_umi_hash = {};
my %nreads_per_umi_hash   = {};
open( OUT_DEBUG_UMIS, ">$debug_umis_fn" )
  || die "ERROR: annot open file $debug_umis_fn.\n";
print OUT_DEBUG_UMIS
  "well_id\tgene\tumi\tis_spike\tstatus\tnreads\tn_offsets\n";
for ( my $i = 0 ; $i <= $#umis ; $i++ ) {
    my $is_spike = 0;
    my @arr      = split( "\t", $umis[$i] );
    if ( exists $spike_ins{ $arr[1] } ) {
        $is_spike = 1;
    }
    my $noffsets_per_cur_umi = $offsets_per_umi{ $umis[$i] };
    my $nreads_per_cur_umi   = $reads_per_umi{ $umis[$i] };

    my $well_id = $arr[0];

    if ( $umi_status{ $umis[$i] } == 1 ) {
        if ( $noffsets_per_cur_umi > 20 ) {
            $noffsets_per_cur_umi = 20;
        }
        $noffsets_per_umi_hash{ $is_spike . "\t"
              . $number_of_cells{$well_id} . "\t"
              . $noffsets_per_cur_umi }++;

        if ( $nreads_per_cur_umi > 100 ) {
            $nreads_per_cur_umi = 100;
        }
        $nreads_per_umi_hash{ $is_spike . "\t"
              . $number_of_cells{$well_id} . "\t"
              . $nreads_per_cur_umi }++;
    }
    print OUT_DEBUG_UMIS $umis[$i] . "\t"
      . $is_spike . "\t"
      . $umi_status{ $umis[$i] } . "\t"
      . $reads_per_umi{ $umis[$i] } . "\t"
      . $offsets_per_umi{ $umis[$i] } . "\n";
}
close(OUT_DEBUG_UMIS);

open( OFFSETS_PER_UMI, ">$noffsets_per_umi_fn" )
  || die "ERROR: cannot open file $noffsets_per_umi_fn.\n";
print OFFSETS_PER_UMI
  "noffsets_per_umi\tspikes\tspike_neg_ctrl\tgenes\tgenes_neg_ctrl\n";
for ( my $i = 0 ; $i <= 20 ; $i++ ) {
    print OFFSETS_PER_UMI $i . "\t"
      . ( 0 + $noffsets_per_umi_hash{"1\t1\t$i"} ) . "\t"
      . ( 0 + $noffsets_per_umi_hash{"1\t0\t$i"} ) . "\t"
      . ( 0 + $noffsets_per_umi_hash{"0\t1\t$i"} ) . "\t"
      . ( 0 + $noffsets_per_umi_hash{"0\t0\t$i"} ) . "\n";
}
close(OFFSETS_PER_UMI);

open( READS_PER_UMI, ">$nreads_per_umi_fn" )
  || die "ERROR: cannot open file $nreads_per_umi_fn.\n";
print READS_PER_UMI
  "noffsets_per_umi\tspikes\tspike_neg_ctrl\tgenes\tgenes_neg_ctrl\n";
for ( my $i = 0 ; $i <= 100 ; $i++ ) {
    print READS_PER_UMI $i . "\t"
      . ( 0 + $nreads_per_umi_hash{"1\t1\t$i"} ) . "\t"
      . ( 0 + $nreads_per_umi_hash{"1\t0\t$i"} ) . "\t"
      . ( 0 + $nreads_per_umi_hash{"0\t1\t$i"} ) . "\t"
      . ( 0 + $nreads_per_umi_hash{"0\t0\t$i"} ) . "\n";
}
close(READS_PER_UMI);

# filt1 UMI error
# filt2 BARCODE error
# filt3 lonely offset error
# filt4 read-poor offset error
