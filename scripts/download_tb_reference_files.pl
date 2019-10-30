#!/usr/bin/env perl
use strict;
use warnings;
use File::Basename;
use File::Spec;
use Cwd;


@ARGV == 1 or die "usage: $0 <outdir>";
my $outdir = $ARGV[0];

unless (-d $outdir) {
    mkdir $outdir or die $!;
}
chdir $outdir or die $!;


my @ntm_accessions = qw/CU458896 JALN00000000 CVQQ00000000 CCAW00000000
ACFI00000000 BCSX00000000 BCQY00000000 AFVW00000000 CCBB00000000 BBHD00000000
CCAY00000000 BBFT00000000 ALQB00000000 JAGZ00000000 CP002385 CP012150
MCHX00000000 CP011883 ARBU00000000 LFOF00000000 LDPO00000000 BBGO00000000
MIGZ00000000 FJVO00000000 CP003322 CP006835 CTEE00000000 AL450380 LAWX00000000
CP019882 JXST00000000 CCBF00000000 CVTB00000000 ANPL00000000
CYSI00000000 JMDW00000000 CWKH00000000 BCTA00000000 JYNU00000000 BBHE00000000
ADNV00000000 CP014475 BCND00000000 BBHB00000000 BBGS00000000 BBHF00000000
MAFR00000000 LDCO00000000 CBMO000000000 JTJW00000000 CBMJ00000000
LN831039 AP018165 MLQM00000000 AGVE00000000 NC_008611
ALQA00000000 CP000511 CP003347 BDDI00000000/;


print "------------------- download genomes --------------------\n";
my %urls = (
    'human.fa.gz' => 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.26_GRCh38.p11/GCA_000001405.26_GRCh38.p11_genomic.fna.gz',
    'NC_001802.1.fa' => 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&rettype=fasta&retmode=text&id=NC_001802.1', # HIV-1
    'NC_007605.1.fa' => 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&rettype=fasta&retmode=text&id=NC_007605.1', # EBV wild type
    'NC_009334.1.fa' => 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&rettype=fasta&retmode=text&id=NC_009334.1', # EBV type 2
    'NC_000962.1.fa' => 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&rettype=fasta&retmode=text&id=NC_000962.1', # TB reference genome v1
    'NC_000962.2.fa' => 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&rettype=fasta&retmode=text&id=NC_000962.2', # TB reference genome v2
    'NC_000962.3.fa' => 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&rettype=fasta&retmode=text&id=NC_000962.3', # TB reference genome v3
    'Saliva.tar.bz2' => 'http://downloads.ihmpdcc.org/data/HMBSA/Saliva.tar.bz2',
    'Throat.tar.bz2' => 'http://downloads.ihmpdcc.org/data/HMBSA/Throat.tar.bz2',
    'Tongue_dorsum.tar.bz2' => 'http://downloads.ihmpdcc.org/data/HMBSA/Tongue_dorsum.tar.bz2',
    'Buccal_mucosa.tar.bz2' => 'http://downloads.ihmpdcc.org/data/HMBSA/Buccal_mucosa.tar.bz2',
    'Palatine_Tonsils.tar.bz2' => 'http://downloads.ihmpdcc.org/data/HMBSA/Palatine_Tonsils.tar.bz2',
);

my @ntm_fastas;
for my $accession (@ntm_accessions) {
    if ($accession =~ /^[A-Z]{4}[0-1]{8}$/) {
        $urls{"$accession.fa.gz"} = "enaDataGet $accession";
        push(@ntm_fastas, "$accession.fa.gz");
    }
    else {
        $urls{"$accession.fa"} = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&rettype=fasta&retmode=text&id=' . $accession;
        push(@ntm_fastas, "$accession.fa");
    }
}

for (keys %urls) {
    download_file($_, $urls{$_});
    if (/\.tar.bz2$/) {
        get_scaffolds_from_tar_bz2($_);
        unlink $_ or die $!;
    }
}


my %genomes = (
    Human => {
        fastas => ['human.fa.gz'],
        is_contam => 1,
    },
    Virus => {
        fastas => ['NC_001802.1.fa', 'NC_007605.1.fa', 'NC_009334.1.fa'],
        is_contam => 1,
    },
    Bacteria => {
        fastas => ['Saliva.tar.bz2.fa', 'Throat.tar.bz2.fa', 'Tongue_dorsum.tar.bz2.fa', 'Buccal_mucosa.tar.bz2.fa', 'Palatine_Tonsils.tar.bz2.fa'],
       is_contam => 1,
    },
    TB => {
        fastas => ['NC_000962.3.fa'],
        is_contam => 0,
    },
    NTM => {
        fastas => \@ntm_fastas,
        is_contam => 0,
    },
);


my %fastas_to_keep = ('NC_000962.1.fa' => 1, 'NC_000962.2.fa' => 1, 'NC_000962.3.fa' => 1);

print "------------------- parse fasta files --------------------\n";
my $remove_contam_fa = 'remove_contam.fa.gz';
open my $metadata_fh, ">remove_contam.tsv" or die $!;
# Pipe through seqtk seq to make all line lengths 60 characters.
# Otherwise bwa index and samtools faidx will break later.
open my $contam_fa_fh, "| seqtk seq -l 60 | gzip -c -9 >$remove_contam_fa" or die $!;

for my $group (keys %genomes) {
    my $is_contam = $genomes{$group}{is_contam};
    for my $fasta (@{$genomes{$group}{fastas}}) {
        print "Start $group, $fasta, $is_contam\n";
        format_fasta_and_update_metadata_tsv($group, $is_contam, $fasta, $metadata_fh, $contam_fa_fh);
        unless (defined $fastas_to_keep{$fasta}) {
            unlink $fasta or die $!;
        }
        print " ... end $group, $fasta\n";
    }
}

close $metadata_fh or die $!;
close $contam_fa_fh or die $!;


sub format_fasta_and_update_metadata_tsv {
    my $group = shift;
    my $is_contam = shift;
    my $infile = shift;
    my $tsv_fh = shift;
    my $fasta_fh = shift;

    my $in_fh;
    if ($infile =~ /\.gz$/) {
        open($in_fh, "gunzip -c $infile |") or die $!;
    }
    elsif ($infile =~ /\.bz2$/) {
        open($in_fh, "bunzip2 -c $infile |") or die $!;
    }
    else {
        open($in_fh, $infile) or die $!;
    }


    while (<$in_fh>){
        chomp;
        next unless /.+/;
        if (/^>(\S+)\s*.*/) {
            print $fasta_fh ">$1\n";
            print $tsv_fh "$group\t$is_contam\t$1\n";
        }
        else {
            print $fasta_fh "$_\n";
        }
    }

    close $in_fh or die $!;
}


sub get_scaffolds_from_tar_bz2 {
    my $infile = shift;
    my @files = `tar --list -f $infile`;
    chomp @files;
    @files = grep(/\.scaffolds.fa.bz2$/, @files);
    @files == 1 or die "Error getting scaffolds file from $infile";
    system_call("tar -O -x -f $infile $files[0] | bunzip2 -c > $infile.fa");
    print "Extracted scaffolds file: $files[0]\n";
}


sub download_file {
    my $outname = shift;
    my $url = shift;
    if (-e $outname and -s $outname == 0) {
        unlink $outname;
    }

    if (-e $outname and -s $outname) {
        print "Skipping $outname,$url - already downloaded\n";
    }
    elsif ($url =~ /^enaDataGet (.*)$/) {
        mkdir "tmp.$$";
        system_call("enaDataGet -d tmp.$$ -f fasta $1");
        system_call("mv tmp.$$/* $outname");
        system_call("rm -r tmp.$$");
    }
    else {
        system_call("wget -nv -O $outname '$url'");
    }
}


sub system_call {
    my $cmd = shift;
    print "$cmd\n";
    if (system($cmd)) {
        die "Error running: $cmd";
    }
}
