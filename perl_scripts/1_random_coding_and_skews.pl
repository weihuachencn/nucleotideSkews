#!/usr/bin/perl -w
use strict;
use POSIX;
use Data::Dumper;
use List::Util qw(shuffle);

## *********************************************
## ** version history **
## *********************************************
my $ver = '1.1';
my $last_modified = 'Oct 4, 2013';

$ver = '1.1';
$last_modified = 'Oct 6, 2013';
## -- also calculate average aa-cost --

## *********************************************
## ** GET opts **
## *********************************************
use Getopt::Long;
my %opts=();
GetOptions(\%opts,"i:s","o:s", "start:f", "end:f", "incr:s", "debug" );

if ( !$opts{i} or !$opts{o} ){
    print "--------------------------------------------------------------------------------------------------
    \t\tversion : $ver by Weihua Chen; last modified : $last_modified
--------------------------------------------------------------------------------------------------
    USAGE: perl $0
        -i sequence length (bp), 2000000, for example
        -o output file
        
      [ optional ]
        -start start GC content, default = 20
        -end end GC content, default = 80
        -incr increasing step, default = 0.1
        -debug, if true, stop at 100 iterations 
--------------------------------------------------------------------------------------------------\n";
    exit;
}

## --- get user input parameters ---
if( defined $opts{debug} ){
    print STDERR "\t=============================================\n\t\truning in debug mode ... \n\t\t will stop at 100 iterations ... \n\t=============================================\n\n";
}


## -- calculate 1st, 2nd and 4th fold generate sites --
my $hrCodes = &codonStats(); ## -- hash{ $codon }{ $pos } = 1,2,4; 0 if the stop codon

## -- load cost of amino acids --
my %hCodon2cost = (); 
open AACOSTS, "cost_of_amino_acid_synthesis.txt" or die;
while(<AACOSTS>){
    next if( /^#/ );
    my ( undef, $aa, $codons, undef, undef, undef, $cost ) = split(/\t/, $_);
    foreach my $codon ( split(/[, ]+/, $codons )){
        $hCodon2cost{ $codon } = $cost;
    }
}
close AACOSTS;

my $len = $opts{i};
open OUT, ">$opts{o}" or die;
print OUT join("\t", "GC skewAT1 skewGC1 skewAT2 skewGC2 skewAT4 skewGC4 costPerAA"), "\n";

my $counts = 0;
for(my $gc = $opts{start}; $gc <=$opts{end}; $gc += $opts{incr}){
    my ( $ats, $gcs ) = ( $len * ( 100 - $gc) / 100, $len * $gc / 100 ); ## -- prepare twice as many as needed in case many of them are stop codons --
    my @arr = ( split(//, "AT" x $ats), split(//, "GC" x $gcs) );
    @arr = shuffle @arr;
    
    print "\t", $counts + 1, "; # of letters : ", scalar @arr, "\n";
    
    my %hValidCodonStats = (); ## -- $hash{ $codon } = $counts --
    my @aValidCodons = ();
    my $valid_codons = $len / 3;
    my $j = 0;
    for( my $i = 0; $i < $valid_codons; ){
        my $codon = join("", @arr[ $j .. $j + 2 ] );
        if( exists $$hrCodes{ $codon } ){
            push @aValidCodons, $codon;
            $hValidCodonStats{ $codon } ++;
            $i ++;
        }
        
        $j += 3;
    }
    
    ## -- calculate costs per aa --
    my $valid_counts = 0;
    my $total_costs = 0;
    while( my ( $codon, $counts ) = each %hValidCodonStats ){
        next if( !exists $hCodon2cost{ $codon } );
        $total_costs += $hCodon2cost{ $codon } * $counts;
        $valid_counts += $counts;
    }
    
    print OUT join("\t", $gc, &ntFrequencyOfDifferentDSites( \@aValidCodons ), $total_costs / $valid_counts ), "\n";
    
    ## -- if debug mode --
    $counts ++;
    if( defined $opts{debug} and $counts >= 100 ){
        close OUT;
        exit;
    }
}
close OUT;


## -- input : arr of strings, [$codon, $codon];
## -- output: skews [  at1 gc1 at2 gc2 at3 gc3 ]
sub ntFrequencyOfDifferentDSites{
    my ( $ar ) = @_;
    
    ## -- note: if $arAllowedAminoAcids is undefined, all amino-acids will be used --
    my %hFirst = ();
    my %hSecond = ();
    my %hFourth = ();
    foreach my $codon ( @{$ar} ){
        my @aBases = split(//, uc $codon);
        for( my $i = 0; $i < @aBases; $i ++ ){
            my $base = $aBases[ $i ];
            my $degeneracy = $$hrCodes{ $codon }{ $i + 1 };
            
            if( $degeneracy == 1 ){
                $hFirst{ $base } ++;
            } elsif( $degeneracy == 2 ){
                $hSecond{ $base } ++;
            } elsif( $degeneracy == 4 ){
                $hFourth{ $base } ++;
            } 
        }
    }
    
    return( &skewCalc( \%hFirst ), &skewCalc( \%hSecond ), &skewCalc( \%hFourth ) );
}

## -- input : hash ref of letters $hash{ A|T|G|C } = $count --
## -- output : [ atskew, gcskew ] 
sub skewCalc{
    my ( $hr ) = @_;
    
    my $a = exists $$hr{ "A" } ? $$hr{ "A" } : 0;
    my $t = exists $$hr{ "T" } ? $$hr{ "T" } : 0;
    my $g = exists $$hr{ "G" } ? $$hr{ "G" } : 0;
    my $c = exists $$hr{ "C" } ? $$hr{ "C" } : 0;
        
    return( ($a - $t) / ( $a + $t ), ( $g - $c ) / ( $g + $c ) );
}


### ---- Jan 07, 2013 : results manually checked ----
## -- according to the webpage at : http://www.bio.ku.dk/nuf/resources/scitab/misc/bactcode.html
sub codonStats{
    
    my %hash = ();
    my @aB1  = split( //, 'TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG' );
    my @aB2  = split( //, 'TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG' );
    my @aB3  = split( //, 'TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG' );
    my @aA   = split( //, 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG' );
    
    ## -- codons --
    my %hCodons = (); ## $hCodons{ $codon } = $AA
    for ( my $i = 0; $i < @aB1; $i ++  ){
        my $aa = $aA[ $i ];
        $aa = 'stop' . $i if( $aa eq '*' );
        
        $hCodons{ $aB1[ $i ] . $aB2[ $i ] . $aB3[ $i ] } = $aa;
    }
    
    #use Data::Dumper;
    #print Dumper %hCodons;
    #print "\n======================================\n";
    
    my @aNT = ('A', 'T', 'G', 'C');
    ## -- degenate count --
    while ( my ( $codon, $aa ) = each %hCodons ){
        next if( $aa =~ /stop/ );
        my @a = split( // , $codon );
        ## -- letter one --
        my $h = 0;
        foreach my $nt ( @aNT ){
            my $key = $nt . $a[1] . $a[2];
            $h ++ if( exists $hCodons{ $key } and $hCodons{ $key } eq $aa );
        }
        $hash{ $codon }{ 1 } = $h;
        
        ## -- letter two --
        $h = 0;
        foreach my $nt ( @aNT ){
            my $key = $a[0] . $nt . $a[2];
            $h ++ if( exists $hCodons{ $key } and $hCodons{ $key } eq $aa );
        }
        $hash{ $codon }{ 2 } = $h;
        
        ## -- letter three --
        $h = 0;
        foreach my $nt ( @aNT ){
            my $key = $a[0] . $a[1] . $nt ;
            $h ++ if( exists $hCodons{ $key } and $hCodons{ $key } eq $aa );
        }
        $hash{ $codon }{ 3 } = $h;
    }
    
    return \%hash;
}