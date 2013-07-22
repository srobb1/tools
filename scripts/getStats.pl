#!/usr/bin/perl -w
use strict;
use Getopt::Long;
if ( !defined @ARGV ) {
  &getHelp();
}
my ($flagstats,$bwas,$characterizer,$compare,$reference);
GetOptions(
  'f|flagstat:s'        => \$flagstats,
  'b|bwa:s'		=> \$bwas,
  'r|reference:s'	=> \$reference,	
  'z|characterizer:s'	=> \$characterizer,	
  'c|compare:s'		=> \$compare,
  'h|help' 		=> \&getHelp,
);

sub getHelp {
  print " 
usage:
$0 [-f flagstat_out_dir][-b bwa_stderr_dir][-r reference_file][-z characterizer_out_file][-c compare_out_file][-h] 

options:
  'f|flagstat:s'        flagstat_dir,
  'b|bwa:s'		bwa_dir,
  'r|reference:s'	reference_file,	
  'z|characterizer:s'	characterizer_file,	
  'c|compare:s'		compare_file
  'h|help' 		getHelp,

";
  exit 1;
}

my @args = ($flagstats,$bwas,$characterizer,$compare);
foreach my $arg(@args){
  getHelp() if !defined $arg;
}

my %stats;

my @flagstats = <$flagstats/*>;
foreach my $flagstat (@flagstats){
  $flagstat =~ /((HEG4_\w|NB_\w|RIL\d+_\d)_[ATGC]{6}_FC\d+L\d+)/;
  my $id = $1;
  my $sample = $2;
  open FLAGSTAT, $flagstat or die "Can't open flagstat file $flagstat for reading:$!\n";
  while(my $line = <FLAGSTAT>){
    chomp $line;
    ## 21946589 + 0 mapped (96.65%:-nan%)
    next unless $line =~/mapped/;
    $line =~/(\d+)\s+\+/;
    $stats{$sample}{$id}{mapped}=$1;
    last;
  }
  close FLAGSTAT;
}

my @bwas = <$bwas/*>;
foreach my $bwa (@bwas){
  $bwa =~ /((HEG4_\w|NB_\w|RIL\d+_\d)_[ATGC]{6}_FC\d+L\d+)/;
  my $id = $1;
  my $sample = $2;
  open BWA, $bwa or die "Can't open BWA file $bwa for reading:$!\n";
  while(my $line = <BWA>){
    chomp $line;
    ## [infer_isize] inferred external isize from 185720 pairs: 167.980 +/- 17.870
    next unless $line =~/\[infer_isize\].*inferred external isize from /;
    $line =~/inferred external isize from \d+ pairs:\s+(.+)$/;
    $stats{$sample}{$id}{insert_size}=$1;
    last;
  }
  close BWA;
}

open CHARACTERIZER, $characterizer or die "Can't open CHARACTERIZER file $characterizer for reading:$!\n";
while(my $line = <CHARACTERIZER>){
  chomp $line;
  next if $line =~ /strain/;
  ## RIL9_0  mping   TTA     Chr8:26035270..26035272 +       2.5     4       heterozygous
  my @line=split /\t/,$line;
  my $class;
  if ($line[7] =~ /homo/){
   $class = 'homo';
  }
  if ($line[7] =~ /het/){
   $class = 'het';
  }
  if ($line[7] =~ /som/){
   $class = 'som';
  }
  $stats{$line[0]}{all}{class}{$class}++;
}
close CHARACTERIZER;

open REFERENCE, $reference or die "Can't open REFERENCE file $reference for reading:$!\n";
my %ref;
while(my $line = <REFERENCE>){
  chomp $line;
  next if $line =~ /strain/;
  ## RIL45_0 mping   Chr10:21716391..21716819        10      6 
  my @line=split /\t/,$line;
  my $strain = $line[0];
  my $pos = $line[2];
  my $total_reads = $line[3] + $line[4];
  if ( !exists $ref{$strain}{$pos} and $total_reads > 0 ){
    $stats{$strain}{all}{sharedWithNB}++ if !exists $ref{$strain}{$pos} and $total_reads > 0;
    $ref{$strain}{$pos}++;
  }  
}
close REFERENCE;

open COMPARE, $compare or die "Can't open COMPARE file $compare for reading:$!\n";
my $found = 0;
while(my $line = <COMPARE>){
  chomp $line;
  ## -----count of how many times this combination of strains share this same set of classifications
  ## RIL15_0 homozygous      69 
  last if ($found and $line =~ /\-/);
  if ($line eq '-----count of how many times this combination of strains share this same set of classifications'){
    $found = 1;
    next;
  }
  if ($found){
    next if $line =~ /^\s*$/;
    next if $line =~ /strain/;
    my @line=split /\s+/,$line;
    ## dont count insertion found in NB and HEG4 as shared with HEG4
    if ($line =~ /,/ and $line =~/HEG4/ and $line !~ /NB/){
      my @strains = split /,/ , $line[0];
      foreach my $strain (@strains){
        next if $strain =~ /HEG4/;
        next if $strain =~ /NB/;
        $stats{$strain}{all}{sharedWithHEG4}+=$line[2];
      }
    }
    if ($line =~ /,/ and $line =~/NB/ and $line !~ /HEG4/ ){
      my @strains = split /,/ , $line[0];
      foreach my $strain (@strains){
        next if $strain =~ /HEG4/;
        next if $strain =~ /NB/;
        $stats{$strain}{all}{sharedWithNB}+=$line[2];
      }
    }
    next if $line =~ /,/;
    my $class;
    if ($line[1] =~ /homo/){
      $class = 'homo';
    }
    if ($line[1] =~ /het/){
      $class = 'het';
    }
    if ($line[1] =~ /som/){
     $class = 'som';
    } 
    $stats{$line[0]}{all}{uniq}{$class}+=$line[2];
  }
}
close COMPARE;


### print stats
print join ("\t", qw(sample coverage total homo_total het_total som_total uniq_total uniq_homo uniq_het uniq_som sharedNB_P sharedHEG4_P)) ,"\n";
my @i_size;
foreach my $sample (keys %stats){
  next if !exists $stats{$sample}{all}{class};
  my $mapped;
  my $homo_total = defined $stats{$sample}{all}{class}{homo} ? $stats{$sample}{all}{class}{homo} : 0;
  my $het_total = defined $stats{$sample}{all}{class}{het} ? $stats{$sample}{all}{class}{het}  : 0;
  my $som_total = defined $stats{$sample}{all}{class}{som} ? $stats{$sample}{all}{class}{som}  : 0;
  my $total = $homo_total + $het_total + $som_total ;
  my $homo_uniq = defined $stats{$sample}{all}{uniq}{homo} ?  $stats{$sample}{all}{uniq}{homo}:0 ;
  my $het_uniq = defined $stats{$sample}{all}{uniq}{het} ? $stats{$sample}{all}{uniq}{het} : 0;
  my $som_uniq = defined $stats{$sample}{all}{uniq}{som} ? $stats{$sample}{all}{uniq}{som}:0 ;
  my $uniq_total = $homo_uniq + $het_uniq + $som_uniq;
  my $sharedWithNB = defined $stats{$sample}{all}{sharedWithNB} ?$stats{$sample}{all}{sharedWithNB}:0 ;
  my $sharedWithHEG4 = defined $stats{$sample}{all}{sharedWithHEG4} ?$stats{$sample}{all}{sharedWithHEG4}:0 ;
  foreach my $id (keys %{$stats{$sample}}){ 
    next if $id eq 'all';
    $mapped += $stats{$sample}{$id}{mapped};
    my $i_size = $stats{$sample}{$id}{insert_size};
    push @i_size , [$sample,$id,$i_size];
  }
  my $coverage = (($mapped * 100) / 373245519) ;
  print join ("\t",$sample,$coverage,$total,$homo_total,$het_total,$som_total,$uniq_total,$homo_uniq,$het_uniq,$som_uniq,$sharedWithNB,$sharedWithHEG4) ,"\n";
}

print "\n\n-insert_size-\n";
print "sample\tflowcell\tinsert_size\n";
foreach my $i_size (sort { $$a[0] cmp $$b[0]} @i_size){
  print join ("\t",@$i_size),"\n";
}
