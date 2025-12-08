#!/usr/bin/env perl

use strict;
use warnings;

#Filter somatic variants CSV for potentially interesting variants

my $header = <>;
print $header;

#Prepare column name hash
my %column_names = ();
my %column_index = ();

my @header = split /\t/, $header;

for my $i (0 .. $#header) {
    $column_names{$header[$i]} = 1;
    $column_index{$header[$i]} = $i;
}

#Subroutine to extract specific field by header/column name
sub getField {
    my ($row_array_ref, $column_name) = @_;

    if (scalar @$row_array_ref != scalar @header) {
        die "Data row has a different number of fields (", scalar @$row_array_ref,") compared to the header (", scalar @header, ")!";
    }

    if (exists $column_names{$column_name}) {
        $row_array_ref->[$column_index{$column_name}];
    } else {
        ".";
    }
}

#Go through all data rows
while (<>) {
    chomp;
    my @F = split /\t/;

    my $pass = 0;

    #Extract relevant annotations
    my $clnsig = getField(\@F, "CLNSIG");
    my $func = getField(\@F, "Func.refGene");
    my $exonic_func = getField(\@F, "ExonicFunc.refGene");
    my $dbscsnv_ada = getField(\@F, "dbscSNV_ADA_SCORE");
    my $dbscsnv_rf = getField(\@F, "dbscSNV_RF_SCORE");
    my $cadd_phred = getField(\@F, "CADD_phred");

    #ClinVar pathogenic
    my $pathogenic = $clnsig =~ m/athogenic?!=i/;
    $pass = 1 if ($pathogenic);

    #Stop, frameshift or splicing
    #Add Startloss? Why not, somatic variants would be a random mess, in any case.
    #Can actually probably just rely on CADD>10 for these special ones - Naah
    $pass = 1 if ($exonic_func =~ m/(stop|start|frameshift|splicing|nonsynonymous)/);

    #dbscSNV_ADA_SCORE or dbscSNV_RF_SCORE >= 0.9 - changed to 0.6 according to author's suggestion
    $pass = 1 if (($dbscsnv_ada ne "." and $dbscsnv_ada >= 0.6) or ($dbscsnv_rf ne "." and $dbscsnv_rf >= 0.6));

    print "$_\n" if ($pass or $. == 1);
}

