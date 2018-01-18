#!/usr/bin/env perl

# Script that is run on each sample when testing the generic pipeline.
# It counts the number of reads in the first reads input file
# and writes the count to a file.
# Simply does this by counting lines and dividing by 4.
# If line count is not a multiple of 4 then dies.

use strict;
use warnings;

print "@ARGV\n";
my ($outdir, $sample, $pool, $isolate,
    $seqrep_id, $seqrep_number, $reads1, $reads2) = @ARGV;

-e $reads1 or die "File not found: $reads1";
-e $reads2 or die "File not found: $reads2";
my $count = `zcat $reads1 | wc -l`;
chomp $count;
print "outdir: $outdir\n";

if ($count %4 == 0) {
    open F, ">$outdir/count.txt" or die $!;
    print F "$reads1\t$count\n";
    close F or die $!;
}
else {
    die "Bad line count for file $reads1: $count";
}

