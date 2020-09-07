#!/usr/bin/perl

use strict;
use warnings;
use Bio::DB::Fasta;

#Script by Dr. Andrew Toseland (A.Toseland@uea.ac.uk).
# Starting from the SINA-aligned, Silva Eukaryote SSU database
# and a list of desired taxonomic groups create:
#   1) A FASTA file containing all sequences from the desired groups
#   2) A list of all NCBI tax_ids for all species selected
#   3) A seq_info.csv file containing details of seq_names, tax_ids etc.
#
# Take from the command line the location of the Silva fasta file and a
# list of required taxonomic groups, one name per line.
# Will also need location of the NCBI names.dmp file to look up the taxa ID

my ($silva_fasta, $ncbi_names, $taxa_list) = @ARGV;
my $silva_db = Bio::DB::Fasta->new($silva_fasta);

# Parse taxon list
my @taxa = &tax_to_list($taxa_list);


# Setup temporary output files
my $tmp_file = './tmp_ids.txt';
my $seq_file = './tmp_seqs.fasta';
my $tax_file = './tmp_tax_list.txt';


# Loop through taxon list, get list of all matching sequence Ids and store in a temp file
foreach (@taxa) {
  my $t = $_;
  my $command = "grep \"$t\" $silva_fasta \| grep -vi \"uncultured\" \|
    grep -vi \"environmental\" \| grep -vi \"unclassified\" \| grep -vi \"unidentified\" \|
    grep -vi \"clone\" \| grep -v \"sp\\.\" \| grep -vi \"symbiont\" \| grep -vi \"metagenome\" \|
    grep -v \"var\\.\" \| grep -v \"cf\\.\" \| sed 's\/>\\(.\*\\)\/\\1\/' >> $tmp_file";
  system("$command");
} # foreach


# Parse silva list
my %id_hash = &parse_silva_list($tmp_file);


open (OUT, ">>$seq_file") or die "Can't open $seq_file\n"; open (TAX, ">>$tax_file") or die "Can't open $taxa_list\n";
print "\"seqname\",\"accession\",\"tax_id\",\"tax_name\",\"is_type\"\n";

# Loop through list of seq IDs, get the species name and lookup the
# NCBI tax ID
# Example Silva ID line:
# AADB02114877.157.2010 Eukaryota;Opisthokonta;Holozoa;Metazoa;Animalia;Craniata;Mammalia;Homo
sapiens (human)

my %done_hash;

foreach my $key (sort keys %id_hash)
{
  if ($id_hash{$key}{str} eq 'N')
  {
    my @levels = split(/;/,$id_hash{$key}{details});
    my $species = $levels[scalar(@levels)-1];
    if ($species =~ m/^(.*)\s+\(.*\)/)
    {
      $species = $1;
    }
    if (!exists $done_hash{$species})
    {
      my $lookup = `grep \"$species\" $ncbi_names`;
      chomp $lookup;
      # Print seq_info.csv, has this format:
      # "seqname","accession","tax_id","tax_name","is_type"
      # "AADB02002333","AADB02002333.14402.16240","9606","Homo sapiens",0
      if ($lookup =~ m/^(\d+)\s+.*/)
      {
        my $tax_id = $1;
        print TAX "$tax_id\n";
        print "\"$key\",\"$key\",\"$tax_id\",\"$species\",";
        if ($id_hash{$key}{str} eq "Y")
        {
          print "1\n";
        }
        else
        {
          print "0\n";
        }
      }
      # Get sequence
      my $seq = $silva_db->seq($key);
      print OUT ">$key\n";
      print OUT "$seq\n";
      $done_hash{$species} = $species;
    }
  }
}

# Tidy up
system("rm $tmp_file");
close OUT or die "Can't close $seq_file\n";
close TAX or die "Can't close $tax_file\n";

system("sort -u $tax_file > tmp_tax_list_uniq.txt");
system("rm $tax_file");



# Example classification:
# Eukaryota;Opisthokonta;Nucletmycea;Fungi;Dikarya;Ascomycota;Saccharomycotina;Saccharomycetes;Saccharomycetidae;Saccharomycetales;Saccharomycetaceae;Saccharomyces;Saccharomyces
# cerevisiae (baker's yeast)
sub parse_silva_list
{
  my $file = shift;
  my %hash;
  open (IN, "<$file") or die "Can't open $file\n";
  while (<IN>)
  {
    my $line = $_;
    chomp $line;
    # Need to check here for strain variants
    if ($line =~ m/^(\S+)\s+(.*)/)
    {
      my ($id, $rest) = ($1,$2);
      $hash{$id}{details} = $rest;
      my @class = split(/;/,$rest);
      # Does the last element contain 3 words
      my $species = $class[scalar(@class)-1];
      if ($species =~ m/^(\S+)\s+(\S+)\s+(.*)/)
      {
        if ($3 =~ m/\(.*\)/)
        {
          $hash{$id}{str} = 'N';
        }
        else
        {
          $hash{$id}{str} = 'Y';
        }
      }
      else
      {
        $hash{$id}{str} = 'N';
      }
    }
  }
  close IN or die "Can't close $file\n";
  return %hash;
}



# Return list of required taxa
sub tax_to_list
{
  my $file = shift;
  my @list;
  open (IN, "<$file") or die "Can't open $file\n";
  while (<IN>)
  {
    my $tax = $_;
    chomp $tax;
    push (@list, $tax);
  }
  close IN or die "Can't close $file\n";
  return @list;
}
