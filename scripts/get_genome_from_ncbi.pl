#!/usr/bin/env perl
use strict;
use warnings;

@ARGV == 2 or die "usage: $0 <accesion> <outfile>";
my $accession = $ARGV[0];
my $outfile = $ARGV[1];

print "Download accession $accession\n";

my $url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&rettype=fasta&retmode=text&id=$accession";
system("wget -nv -O $outfile.$$ '$url'") and die "Error downloading $url";

# They put an empty line at the end! Remove it.
# And remove everything after first whitespace in header
open F_IN, "$outfile.$$" or die $!;
open F_OUT, ">$outfile" or die $!;

while (<F_IN>) {
    chomp;
    next unless /.+/;
    if (/^(>\S+)\s+.*/) {
        print F_OUT "$1\n";
    }
    else {
        print F_OUT "$_\n";
    }
}

close F_IN or die $!;
close F_OUT or die $!;
unlink "$outfile.$$" or die $!;
