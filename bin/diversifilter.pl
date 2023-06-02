#!/usr/bin/perl

# this is a simple tool for filtering the entropy file output from diversitools
# the user specifies the p-value cut-off for the nucleotide significance
# calculated in diversitools and a p-value cut-off for the fisher's exact test 
# for the strand bias calculated in a similar way to v-phaser
# 
# Entropy gets recalculated, order of nucleotides gets amended, Tv and Ts gets recalculated
# TO DO:
# NonSyn and Syn CANNOT be recalculated at this time. With Orf info, they could be recalculated in a dirty way
# Recalculate the AA table
#
# example command: ./diversifilter -in out -pQ 0.05 -pS 0.05 -stub out_filtered

use strict;
use Getopt::Long; 
use Math::CDF;
use List::Util 'sum'; # for counting the values in the inserfreq and delfreq hash
use Text::NSP::Measures::2D::Fisher::left;

# global variables
my $stub="out_filtered";
my $pQ=1;
my $pS=1;
my ($in, $help);

&GetOptions(
        'in:s'  => \$in,#the stub for the input file
        'pQ:s'  => \$pQ, #the quality p-value used as a cut-off
        'pS:s'  => \$pS, #the p-value to use as a cut-off for filtering strand-bias  
        "stub:s" => \$stub,
           );

if (($help)&&!($help)||!($in)){
  print "Usage: perl diversifilter.pl -in out -pQ 0.05 -pS 0.05 -stub out_filtered \n";
  print "Usage for mac compiled version: diversifilter_macosx -in out -pQ 0.05 -pS 0.05 -stub out_filtered \n";
  print "Usage for linux compiled version: diversifilter_linux -in out -pQ 0.05 -pS 0.05 -stub out_filtered \n";
  print " -in <txt>  - the stub of a diversitool output\n";
  print " -pQ <txt>  - the quality cut-off for a base (default 0.05)\n";
  print " -pS <txt> - the p-value cut-off for strand-bias (default 0.05)\n";
  print " -stub <txt> - the output in text-tab delimited\n";
  print " -help        - Get this help\n";
  exit();
}

my (%table,%strandbias);
open (ENTROPY,"<$in\_entropy.txt")||die "Can't open $in\_entropy.txt\n";
my $firstline = <ENTROPY>;
chomp($firstline);
my @headers=split(/\t/,$firstline); 
while(<ENTROPY>){
  chomp($_);
  my @elements=split(/\t/,$_);
  my (@order,%strandinfo);
  for (my $i=3; $i<scalar(@headers); $i++){
    $table{$elements[0]}{$elements[1]}{$elements[2]}{$headers[$i]}=$elements[$i];#$table{Sample}{Chr}{Position}
    # print "$headers[$i]\t$elements[$i]\n";
    if ($headers[$i]=~/OrderOfNucs/){
      @order=split(//,$elements[$i]);
    }
    if ($headers[$i]=~/strandbias/){
      if ($elements[$i]=~/A(\d+):(\d+);C(\d+):(\d+);T(\d+):(\d+);G(\d+):(\d+)/){
        $strandinfo{"A"}{"plus"}=$1;
        $strandinfo{"A"}{"neg"}=$2;
        $strandinfo{"C"}{"plus"}=$3;
        $strandinfo{"C"}{"neg"}=$4;
        $strandinfo{"T"}{"plus"}=$5;
        $strandinfo{"T"}{"neg"}=$6;
        $strandinfo{"G"}{"plus"}=$7;
        $strandinfo{"G"}{"neg"}=$8;
      }
    }
  }
  $strandbias{$elements[0]}{$elements[1]}{$elements[2]}{$order[0]}=0;
  for (my $i=1; $i<scalar(@order); $i++){
    #print "$order[0] $strandinfo{$order[0]}{plus} $strandinfo{$order[0]}{neg}\t";
    #print "$order[$i] $strandinfo{$order[$i]}{plus} $strandinfo{$order[$i]}{neg}\n";
    #          TopNuc   OtherNuc
    #  Pos    n11      n12 | n1p
    #  Neg    n21      n22 | n2p
    #           --------------
    #           np1      np2   npp

    my $npp = $strandinfo{$order[0]}{plus}+$strandinfo{$order[0]}{neg}+$strandinfo{$order[$i]}{plus}+$strandinfo{$order[$i]}{neg}; 
    my $n1p = $strandinfo{$order[0]}{plus}+$strandinfo{$order[$i]}{plus}; 
    my $np1 = $strandinfo{$order[0]}{plus}+$strandinfo{$order[0]}{neg};  
    my $n11 = $strandinfo{$order[0]}{plus};
 
    my $left_value = calculateStatistic( n11=>$n11,
                                         n1p=>$n1p,
                                         np1=>$np1,
                                         npp=>$npp);
   #print "Significance os strand bias $left_value\n"; 
   
   $strandbias{$elements[0]}{$elements[1]}{$elements[2]}{$order[$i]}=$left_value;
  }
}

# produces a table with Sample\tChr\tPosition\tRefBase\tCoverage\t
# AvQual\tAcnt\tApval\tCcnt\tCpval\tTcnt\tTpval\tGcnt\tGpval\tentropy(base 2)\tNonRefCnt
# CntTv\tCntTs\tCntNonSyn\tCntSyn\tOrderATCG\tinscnt\tdelcnt\tstrandbias\n

# Sample	Chr	Position	RefBase	Coverage	AvQual	Acnt	Apval	Ccnt	Cpval	Tcnt	Tpval	Gcnt	Gpval	entropy(base e)	NonRefCnt	CntTv	CntTs	OrderOfNucs	strandbias	ins_cnt	ins_mode	del_cnt	del_mode

my @nuc=qw/A C T G/;
open(OUT,">$stub\_filter.txt")||die "Can't open $stub\_filter.txt\n";
print OUT join("\t",@headers,"\n");
for my $sample (keys %table){
  for my $chr (keys %{$table{$sample}}){
    for my $position (sort {$a <=> $b} keys %{$table{$sample}{$chr}}){
      
      my ($Acnt,$Ccnt,$Tcnt,$Gcnt);
      if ($table{$sample}{$chr}{$position}{"Apval"}<=$pQ && $strandbias{$sample}{$chr}{$position}{"A"}<=$pS){
        $Acnt=$table{$sample}{$chr}{$position}{"Acnt"};
      }else{
        $Acnt=0;
      }
      if ($table{$sample}{$chr}{$position}{"Cpval"}<=$pQ && $strandbias{$sample}{$chr}{$position}{"C"}<=$pS){
        $Ccnt=$table{$sample}{$chr}{$position}{"Ccnt"};
      }else{
        $Ccnt=0;
      }
      if ($table{$sample}{$chr}{$position}{"Gpval"}<=$pQ && $strandbias{$sample}{$chr}{$position}{"G"}<=$pS){
        $Gcnt=$table{$sample}{$chr}{$position}{"Gcnt"};
      }else{
        $Gcnt=0;
      }
      if ($table{$sample}{$chr}{$position}{"Tpval"}<=$pQ  && $strandbias{$sample}{$chr}{$position}{"T"}<=$pS){
        $Tcnt=$table{$sample}{$chr}{$position}{"Tcnt"};
      }else{
        $Tcnt=0;
      }
      # recalculate order
       my %NucHash = ('A' => $Acnt, 'C' => $Ccnt, 'T' => $Tcnt, 'G' => $Gcnt);
      my $NucOrder="";
      foreach my $nuc (sort {$NucHash{$b} <=> $NucHash{$a}} keys %NucHash){
        if ($NucHash{$nuc}>0){
          $NucOrder=$NucOrder.$nuc;
        }
      }
      # recalculate entropy
      my $shannon=0;
      my $nucnt=$Acnt+$Ccnt+$Gcnt+$Tcnt;
      if ($nucnt>0){
#      if ($table{$sample}{$chr}{$position}{"Coverage"}>0){
#        foreach my $nuc (@nuc){
          
          # not checked yet, looks wrong
          #my $p = $nucnt / $table{$sample}{$chr}{$position}{"Coverage"};
#          my $p = $nuc / $nucnt;
          # changed nucnt to nuc
          if($Acnt > 0){
            my $p = $Acnt/$nucnt;
            $shannon += -$p*log($p);#natural log, i.e. log base e
          }
           if($Ccnt > 0){
            my $p = $Ccnt/$nucnt;
            $shannon += -$p*log($p);#natural log, i.e. log base e
          }
          if($Gcnt > 0){
            my $p = $Gcnt/$nucnt;
            $shannon += -$p*log($p);#natural log, i.e. log base e
          }
          if($Tcnt > 0){
            my $p = $Tcnt/$nucnt;
            $shannon += -$p*log($p);#natural log, i.e. log base e
          }
         
#        }   
      }else{
        $shannon="<NA>";
      }   
      #recalculate CntTv and CntTS
      my ($Ts,$Tv) = cntTsTv($table{$sample}{$chr}{$position}{"RefBase"},$Acnt,$Tcnt,$Gcnt,$Ccnt);
      
      print OUT "$sample\t$chr\t$position\t".$table{$sample}{$chr}{$position}{"RefBase"}."\t";
      print OUT $table{$sample}{$chr}{$position}{"Coverage"}."\t";
      print OUT $table{$sample}{$chr}{$position}{"AvQual"}."\t";
      print OUT $Acnt."\t".$table{$sample}{$chr}{$position}{"Apval"}."\t";
      print OUT $Ccnt."\t".$table{$sample}{$chr}{$position}{"Cpval"}."\t";
      print OUT $Tcnt."\t".$table{$sample}{$chr}{$position}{"Tpval"}."\t";
      print OUT $Gcnt."\t".$table{$sample}{$chr}{$position}{"Gpval"}."\t";
      print OUT $shannon."\t";
      print OUT $table{$sample}{$chr}{$position}{"NonRefCnt"}."\t";
      print OUT $Tv."\t".$Ts."\t";
      print OUT $NucOrder."\t";
      print OUT $table{$sample}{$chr}{$position}{"strandbias"}."\t";
      print OUT $table{$sample}{$chr}{$position}{"ins_cnt"}."\t";
      print OUT $table{$sample}{$chr}{$position}{"ins_mode"}."\t";
      print OUT $table{$sample}{$chr}{$position}{"del_cnt"}."\t";
      print OUT $table{$sample}{$chr}{$position}{"del_mode"}."\n";
    }
  }
}


sub cntTsTv{
  #count transitions A<=>G (purine to purine) and C<=>T *pyrimidine to pyrimidines
  #count transversions A<=>C, A<=>T, G<=>T, G<=>C 
  my ($refbase, $Acnt, $Tcnt, $Gcnt, $Ccnt)= @_;
  my $Ts=0;
  my $Tv=0;
  if ($refbase=~/A/){
    $Ts=$Gcnt;
    $Tv=$Ccnt+$Tcnt;
  }
  if ($refbase=~/T/){
    $Ts=$Ccnt;
    $Tv=$Acnt+$Gcnt;
  }
  if ($refbase=~/G/){
    $Ts=$Acnt;
    $Tv=$Tcnt+$Ccnt;
  }
  if ($refbase=~/C/){
    $Ts=$Tcnt;
    $Tv=$Acnt+$Gcnt;
  }
  #print "$refbase A:$Acnt T:$Tcnt G:$Gcnt C:$Ccnt Transition:$Ts Transversion:$Tv\n";
  return($Ts,$Tv);
}





