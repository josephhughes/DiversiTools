#!/usr/bin/perl

# use this to obtain the consensus sequence from a bam file
# this is a simple consensus program that just takes the majority at a site 
# irrespecitve of depth
# it will include insertions and deletions if their count is greater than the 
# 50% of the reads in the previous site
# the most frequently found insertion will be used
# 
# perl wysiwyg_consensus.pl -bam SamTestFiles/S1_refHPAI_cons_stampy.bam -ref SamTestFiles/refHPAI_cons.fa -stub out
# mismatches and insertions relative to reference marked in lower-case
# deletions are marked with hyphen (-)
# missing reads (no coverage) is assigned N
# perl wysiwyg_consensus.pl -bam 6146_sort.bam -ref Entero_68_AB601883.fasta -stub entero

use strict;
use Getopt::Long; 
use Bio::DB::Sam;
use Math::CDF;
use Bio::SeqIO;
use List::Util 'sum'; # for counting the values in the inserfreq and delfreq hash

# global variables
my ($bam, $ref,$help, $out,%basefreq,%cumulqual,%refbase,%delfreq,%insfreq);
my $stub="output";

&GetOptions(
	    'bam:s'  => \$bam,#bam file (binary sam)
	    'ref:s'  => \$ref, #reference file in fasta  
        "stub:s" => \$stub,
           );

if (($help)&&!($help)||!($bam)||!($ref)){
 print "Usage : perl wysiwyg_consensus.pl -bam input.bam -ref ref.fa -stub stub \n";
 print " -bam <txt>  - the input a file in bam format\n";
 print " -ref <txt>  - the reference fasta file\n";
 print " -stub <txt> - the output in text-tab delimited\n";
 print " -help        - Get this help\n";
 exit();
}

# high level API
my $sam = Bio::DB::Sam->new(-bam  => $bam,
							 -fasta=> $ref);
my @targets    = $sam->seq_ids;
my %targetHash = map { $_ => 1 } @targets;
my $ins_cnt=0;
my $nocigs_cnt=0;
my $Ncnt=0;
# get a hash of the reference sequences with gene, site, refbase
my %refseq;
my $seq_in  = Bio::SeqIO->new(-format => 'fasta',-file   => $ref);
while( my $seq = $seq_in->next_seq() ) {
  my $id=$seq->display_id();
  if ($targetHash{$id}){
  my $seq_str=$seq->seq();
  my @refbases=split(//,$seq_str);
  for (my $i=0; $i<scalar(@refbases); $i++){
    my $site=1+$i;
    $refseq{$id}{$site}=$refbases[$i];
  }
  }else{
    print "$id is not in the bam file\n";
  }
}



my $increment;
foreach my $target (@targets){
if (!$refseq{$target}){# check that reference is present in the reference file
  print "$target does not exist in the reference input file\n";
}else{
    print "\n$target\n";
    $increment=0;
    my @alignments = $sam->get_features_by_location(-seq_id => $target);
    my $refobj = $sam->segment(-seq_id => $target);
    my $refseq = $refobj->dna;
    my $total = scalar(@alignments);
    my $progress=0;
    for my $a (@alignments) {
      $progress++;
      progress_bar( $progress, $total, 25, '=' );
      # where does the alignment start in the reference sequence
      my $seqid  = $a->seq_id;
      my $name   = $a->display_name;
      my $start  = $a->start; #position in the reference
      my $end    = $a->end;
      my $strand = $a->strand;
      my $cigar  = $a->cigar_str;
      my $paired = $a->get_tag_values('PAIRED');
      # where does the alignment start in the query sequence
      my $query_start = $a->query->start;     
      my $query_end   = $a->query->end;
      my $ref_dna   = $a->dna;        # reference sequence bases
      my $query_dna = $a->query->dna; # query sequence bases
      my @scores    = $a->qscore;     # per-base quality scores
      my $match_qual= $a->qual;       # quality of the match
      my $mismatches=0;
      # print "Ref $start Cigar $cigar $ref_dna $query_dna\n";
      # parse the cigar string and modify the query_dna accordingly
      #print "$name $cigar\n";
      #print "$query_dna\n$ref_dna\n";
      my @cigars=split(/(\d+[MIDNSHPX])/,$cigar);# might want to add = as an option for identical match to reference
      # need a $site variable for the site in the reference and a $readpos for the position in the read
      # these variables depend on the cigar information
      my $site=$start;
      my $refpos=0;
      my $readpos=0;
      for (my $k=0; $k<scalar(@cigars);$k++){
        if ($cigars[$k]){
          # listing all the different cigar type
          if ($cigars[$k]=~/(\d+)[M=X]$/){#alignment match (can be a sequence match or mismatch)
            # we want to deal with all of these
            my $match_len=$1;
            #need to trim $query_dna and $ref_dna for the region of the match
            my $sub_query_dna = substr($query_dna,$readpos,$match_len);
            my $sub_ref_dna=substr($ref_dna,$refpos,$match_len);
            # re-using code from btcutils
            my @bases=split(//,$sub_query_dna);
            my @refbases=split(//,$sub_ref_dna);
            my $cumP=0;
            my $matches=0;
            #print "$cigars[$k] $site $readpos $refpos Strand $strand\n";
            for (my $i=0; $i<$match_len; $i++){
              # chromosome site nuc 
              #print "For loop $site $readpos $i\n";
              #print "$site\t$readpos\n";
              $basefreq{$target}{$site}{$strand}{$bases[$i]}++;
              #print "$bam\t$target\t$site\t$bases[$i]\t$basefreq{$target}{$site}{$strand}{$bases[$i]}\n"; 
              $refbase{$target}{$site}=$refbases[$i];
              my $P = 10**(-$scores[$i]/10);
              if ($cumulqual{$target}{$site}){
                $cumulqual{$target}{$site}=$cumulqual{$bam}{$target}{$site}+$P;
              }else{
                $cumulqual{$target}{$site}=$P;
              }
                $site=$site+1;
                $readpos=$readpos+1;
                $refpos=$refpos+1;
                #print "$name $k\t$cigars[$k] $site $readpos $refpos => match\n";
              }
            #print "$cigars[$k] $site $readpos $refpos => match\n";
            #print "$sub_query_dna\n$sub_ref_dna\n";
            
          }if ($cigars[$k]=~/(\d+)D$/){# deletion relative to reference
            #we want to gather the location and length distribution of these assigning the location to the reference position before the deletion
            my $del_len=$1;
            my $sub_ref_dna = substr($ref_dna,$refpos,$del_len);
            $delfreq{$target}{$site}{$sub_ref_dna}++;
            $refpos=$refpos+$del_len;
            $site=$site+$del_len;
            #print "$cigars[$k] $site $readpos $refpos $sub_ref_dna => deletion\n";
          }if ($cigars[$k]=~/(\d+)I$/){# insertion in the read relative to the reference
            #we want to gather the location and length distribution of these
            my $ins_len=$1;
            my $sub_query_dna = substr($query_dna,$readpos,$ins_len);
            $insfreq{$target}{$site-1}{$sub_query_dna}++;
            # $insfreq{$target}{$site}{$sub_query_dna}++;
            $readpos=$readpos+$ins_len;
            #print "$cigars[$k] $site $readpos $refpos $sub_query_dna => insertion\n";
          }if ($cigars[$k]=~/(\d+)N$/){# skipped region from the reference
            #we want to ignore this or maybe treat this in a similar way to an insertion
            $refpos=$refpos+$1;
            $site=$site+$1;
          }if ($cigars[$k]=~/(\d+)S$/){# soft clipping (clipped sequences present in SEQ)
            # we want to skip this region in the read
            $readpos=$readpos+$1;
            #print "$k\t$cigars[$k] \n";
          }if ($cigars[$k]=~/(\d+)H$/){# hard clipping (clipped sequences NOT present in SEQ)
            # we want to ignore this
            #print "$k\t$cigars[$k] $1 \n";
          }if ($cigars[$k]=~/(\d+)P$/){# padding (silent deletion from padded reference)
            # we can ignore this
            
            #print "$k\t$cigars[$k] $1\n";
#           }if ($cigars[$k]=~/(\d+)=$/){# sequence match
#             #we can deal with these the same way as M
#             # we need to count these sites
#             
#             print "$k\t$cigars[$k] $1\n";
#           }if ($cigars[$k]=~/(\d+)X$/){# sequence mismatch
#             # we need to count these sites
#             print "$k\t$cigars[$k] $1\n";
          }
        }
      }
    }
  }
}

my %shannon;#Shannon-Wiener Diversity Index (also known as Shannon's diversity index, the Shannon-Wiener index, the Shannon-Weaver index and the Shannon entropy)

my @nuc=qw/A C T G/;
#open (OUT, ">$stub\_entropy.txt")||die "can't open $stub\_entropy.txt\n";
# add the reference base in the table
# Nucleotide Table:
# Sample\tChr\tPosition\tRefBase\tCoverage\tAvQual\tAcnt\tApval\tCcnt\tCpval\tTcnt\tTpval\tGcnt\tGpval\tentropy(base e)\tNonRefCnt
# CntTv\tCntTs\tOrderATCG
open (CONS,">$stub\_cons.fa")||die "can't open $stub\_cons.fa\n";

print "Test this\n";
#print OUT "Sample\tChr\tPosition\tRefBase\tCoverage\tAvQual\tAcnt\tApval\tCcnt\tCpval\tTcnt\tTpval\tGcnt\tGpval\tentropy(base e)\tNonRefCnt\tCntTv\tCntTs\tOrderOfNucs\tstrandbias\tins_cnt\tins_mode\tdel_cnt\tdel_mode\n";
foreach my $gene (keys %refseq){
  my $nbsites=0;
  my $sumentropy=0;
  my $rawsumentropy=0;
  print CONS ">$stub\_$gene\n";
  my $gene_len=(sort {$b<=>$a} keys %{$refseq{$gene}})[0];
  print "Test".$gene_len;
  #foreach my $site (sort {$a<=>$b} keys %{$refseq{$gene}}){
  for (my $site=1; $site<=$gene_len; $site++){
    #print "$site $gene_len\n";
    if (keys %{$basefreq{$gene}{$site}}){# if there is information in the bam file about this site
      # strandbias info
      my ($Aplus,$Aneg,$Cplus,$Cneg,$Tplus,$Tneg,$Gplus,$Gneg);
      if ($basefreq{$gene}{$site}{1}{"A"}){ $Aplus=$basefreq{$gene}{$site}{1}{"A"}}else{$Aplus=0}
      if ($basefreq{$gene}{$site}{-1}{"A"}){ $Aneg=$basefreq{$gene}{$site}{-1}{"A"}}else{$Aneg=0}
      if ($basefreq{$gene}{$site}{1}{"C"}){ $Cplus=$basefreq{$gene}{$site}{1}{"C"}}else{$Cplus=0}
      if ($basefreq{$gene}{$site}{-1}{"C"}){ $Cneg=$basefreq{$gene}{$site}{-1}{"C"}}else{$Cneg=0}
      if ($basefreq{$gene}{$site}{1}{"T"}){ $Tplus=$basefreq{$gene}{$site}{1}{"T"}}else{$Tplus=0}
      if ($basefreq{$gene}{$site}{-1}{"T"}){ $Tneg=$basefreq{$gene}{$site}{-1}{"T"}}else{$Tneg=0}
      if ($basefreq{$gene}{$site}{1}{"G"}){ $Gplus=$basefreq{$gene}{$site}{1}{"G"}}else{$Gplus=0}
      if ($basefreq{$gene}{$site}{-1}{"G"}){ $Gneg=$basefreq{$gene}{$site}{-1}{"G"}}else{$Gneg=0}
      my $strandbias="A$Aplus".":$Aneg;"."C$Cplus".":$Cneg;"."T$Tplus".":$Tneg;"."G$Gplus".":$Gneg";

      my ($cntA,$cntC,$cntT,$cntG);
#       if ($basefreq{$bam}{$gene}{$site}{1}{"A"} || $basefreq{$bam}{$gene}{$site}{-1}{"A"}){ $cntA = $basefreq{$bam}{$gene}{$site}{1}{"A"}+$basefreq{$bam}{$gene}{$site}{-1}{"A"}}else{$cntA=0}
#       if ($basefreq{$bam}{$gene}{$site}{1}{"C"} || $basefreq{$bam}{$gene}{$site}{-1}{"C"}){ $cntC = $basefreq{$bam}{$gene}{$site}{1}{"C"}+$basefreq{$bam}{$gene}{$site}{-1}{"C"}}else{$cntC=0}
#       if ($basefreq{$bam}{$gene}{$site}{1}{"T"} || $basefreq{$bam}{$gene}{$site}{-1}{"T"}){ $cntT = $basefreq{$bam}{$gene}{$site}{1}{"T"}+$basefreq{$bam}{$gene}{$site}{-1}{"T"}}else{$cntT=0}
#       if ($basefreq{$bam}{$gene}{$site}{1}{"G"} || $basefreq{$bam}{$gene}{$site}{-1}{"G"}){ $cntG = $basefreq{$bam}{$gene}{$site}{1}{"G"}+$basefreq{$bam}{$gene}{$site}{-1}{"G"}}else{$cntG=0}
      $cntA=$Aplus+$Aneg;
      $cntC=$Cplus+$Cneg;
      $cntT=$Tplus+$Tneg;
      $cntG=$Gplus+$Gneg;
      my $coverage = $cntA + $cntT + $cntC + $cntG;
      my $average_p=$cumulqual{$gene}{$site}/$coverage;
      my $p = $average_p/3;
#       my %prob;
#       $prob{"A"} = 1 - (&Math::CDF::pbinom(($cntA-1), $coverage, $p));# need to double check with Richard about the -1
#       $prob{"C"} = 1 - (&Math::CDF::pbinom(($cntC-1), $coverage, $p));
#       $prob{"T"} = 1 - (&Math::CDF::pbinom(($cntT-1), $coverage, $p));
#       $prob{"G"} = 1 - (&Math::CDF::pbinom(($cntG-1), $coverage, $p));
#       if ($coverage>0){
#         foreach my $nuc (@nuc){
#           my $nucnt=$basefreq{$gene}{$site}{1}{$nuc}+$basefreq{$gene}{$site}{-1}{$nuc};
#           my $p = $nucnt / $coverage;
#           if($nucnt > 0){
#             $shannon{$gene}{$site} += -$p*log($p);#natural log, i.e. log base e
#           }
#         }   
#       }else{
#         $shannon{$gene}{$site}="<NA>";
#       }   
      $nbsites++;
      $sumentropy=$sumentropy+$shannon{$gene}{$site};
      my $refbase=$refbase{$gene}{$site};
      my $nonrefcnt=$coverage-($basefreq{$gene}{$site}{1}{$refbase})-($basefreq{$gene}{$site}{-1}{$refbase});
      my ($Ts,$Tv) = cntTsTv($refbase,$cntA,$cntT,$cntG,$cntC);
      my $NucOrderPlus="";
      my $NucOrderNeg="";
      foreach my $nuc (sort { $basefreq{$gene}{$site}{1}{$b} <=> $basefreq{$gene}{$site}{1}{$a} } keys %{$basefreq{$gene}{$site}{1}}) {
        $NucOrderPlus=$NucOrderPlus.$nuc;
      }
      foreach my $nuc (sort { $basefreq{$gene}{$site}{-1}{$b} <=> $basefreq{$gene}{$site}{-1}{$a} } keys %{$basefreq{$gene}{$site}{-1}}) {
        $NucOrderNeg=$NucOrderNeg.$nuc;
      }
      my %NucHash = ('A' => $cntA, 'C' => $cntC, 'T' => $cntT, 'G' => $cntG);
      my $NucOrder="";
      foreach my $nuc (sort {$NucHash{$b} <=> $NucHash{$a}} keys %NucHash){
        if ($NucHash{$nuc}>0){
          $NucOrder=$NucOrder.$nuc;
        }
      }
      # check this $topNuc for the case where there is a deletion
      
      my $topNuc=(sort {$NucHash{$b} <=> $NucHash{$a}} keys %NucHash)[0];
      
      #print "$site $refbase $NucOrder\n";	
      # get the cnt of insertion and mode and for deletions
      my ($del_cnt,$ins_cnt,$freq_ins,$freq_del);
      if (keys %{$delfreq{$gene}{$site}}){ $del_cnt= sum values %{$delfreq{$gene}{$site}}}else{$del_cnt="<NA>"}
      if (keys %{$insfreq{$gene}{$site}}){ $ins_cnt= sum values %{$insfreq{$gene}{$site}}}else{$ins_cnt="<NA>"}
      # most frequently found insertion and deletion
      #my $freq_ins = &largest_value_mem(%{$insfreq{$gene}{$site}});
      #my $freq_del = &largest_value_mem(%{$delfreq{$gene}{$site}});
     # need to sort in descending order so that the most frequent is first 
      if (keys %{$insfreq{$gene}{$site}}){$freq_ins = (sort {$insfreq{$gene}{$site}{$b} <=> $insfreq{$gene}{$site}{$a}} keys %{$insfreq{$gene}{$site}})[0]}else{$freq_ins="<NA>"};
      if (keys %{$delfreq{$gene}{$site}}){$freq_del = (sort {$delfreq{$gene}{$site}{$b} <=> $delfreq{$gene}{$site}{$a}} keys %{$delfreq{$gene}{$site}})[0]}else{$freq_del ="<NA>"};
      #print OUT "$bam\t$gene\t$site\t".uc($refbase)."\t$coverage\t$average_p\t$cntA\t".$prob{"A"}."\t$cntC\t".$prob{"C"}."\t$cntT\t".$prob{"T"}."\t$cntG\t".$prob{"G"}."\t";
      #print OUT "$shannon{$gene}{$site}\t$nonrefcnt\t$Ts\t$Tv\t$NucOrder\t$strandbias\t$ins_cnt\t$freq_ins\t$del_cnt\t$freq_del\n";
      if (uc($refbase) eq uc($topNuc)){ 
      #}elsif (uc($refbase) eq uc($topNuc)){
        print CONS uc($refbase);
      }else{
        print CONS lc($topNuc);
      }if ($ins_cnt>$coverage/2){
        print CONS lc($freq_ins);
      }if ($del_cnt>$coverage){
        for (my $del_site=0; $del_site<length($freq_del); $del_site++){
          print CONS "-";
        }
        $site=$site+length($freq_del);
      }

    }else{#there is no coverage for that site in the bam
      # this should be a deletion as there is no ATGC at this site 
      # check that previous site does not have del or insertion associated
        #print "No coverage $site\n";
        #print OUT "$bam\t$gene\t$site\t".uc($refseq{$gene}{$site})."\t0\t<NA>\t<NA>\t<NA>\t<NA>\t<NA>\t<NA>\t<NA>\t<NA>\t<NA>\t";
        #print OUT "<NA>\t<NA>\t<NA>\t<NA>\t<NA>\t<NA>\t<NA>\t<NA>\t<NA>\t<NA>\n";
        if (keys %{$insfreq{$gene}{$site}}){
          print "Warning: There are insertions in $gene at site $site\n";
        }
        if (keys %{$insfreq{$gene}{$site}}){
          print "Warning: There are deletions in $gene at site $site\n";
        }else{
          #print CONS "N";
          print CONS "-";
        }
    }
  }
  print LOG "Gene $gene Average entropy = ".$sumentropy/$nbsites."\n";
  print CONS "\n";
}
sub progress_bar {
    my ( $got, $total, $width, $char ) = @_;
    if ($got/$total>=$increment){
    $width ||= 25;
    $char  ||= '=';
    my $num_width = length $total;
    local $| = 1;
    $increment=$increment+0.05;
#     printf "|%-${width}s| Processed %${num_width}s reads of %s (%.2f%%)\r", 
#         $char x (($width-1)*$got/$total). '>', $got, $total, 100*$got/+$total;
    printf "|%-${width}s| Processed %.2f%% of reads out of %s\r", 
        $char x (($width+1)*$got/$total). '>', 100*$increment, $total;

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
