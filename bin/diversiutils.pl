#!/usr/bin/perl

# diversiutils.pl is a re-write of btcutils.pl
# with strand bias and insertion deletions taken into account, median of insertion sites
# use this to obtain the number of each base at each site from a bam file
# you need to give the bam file and the reference fasta
# it also calculate the average per site entropy for the sample
# 
# produces a table with Sample\tChr\tPosition\tRefBase\tCoverage\t
# AvQual\tAcnt\tApval\tCcnt\tCpval\tTcnt\tTpval\tGcnt\tGpval\tentropy(base e)\tentropy(base 2)\tsimpson\tNonRefCnt
# CntTv\tCntTs\tCodonPos\tCntNonSyn\tCntSyn\tOrderATCG\tinscnt\tdelcnt\tstrandbias\n
# example command: ./diversiutils.pl -bam test/S3_refHPAI_cons_stampy.bam -ref test/refHPAI_cons.fa -stub out
# perl diversiutils.pl -bam test/S3_refHPAI_cons_stampy.bam -ref test/refHPAI_cons.fa -stub out
# option to include an orf file which will provide information about the amino acids
# ./diversiutils.pl -bam test/S3_refHPAI_cons_stampy.bam -ref test/refHPAI_cons.fa -orf test/Coding.regions.HPAI.txt
# counts the Ns in the  reads and Ns are taken into account in the coverage
# CDS can be on the reverse stand of the reference.

# Modified 06-02-2017 : check NonRefCnt, changed refbase to upper case
# Change reference to uppercase and convert IUPAC to N: this needs to be done for the input reference and the reference in the bam
# This affects the CntRef value in the table, N could code for A,C,T,G so most bases will be matching to ref 

use strict;
use Getopt::Long; 
use Bio::DB::Sam;
use Math::CDF;
use Bio::SeqIO;
use List::Util 'sum'; # for counting the values in the inserfreq and delfreq hash

# global variables
my ($bam, $ref, $orfs,$help);
my (%basefreq,%refbase,%cumulqual,%readinfo,%delfreq,%insfreq,%readmis);# data storing 
my (%aafreq,%aaorder,%stopfreq,%codonfreq);# storing of AA info
my $stub="output";

&GetOptions(
        'bam:s'  => \$bam,#bam file (binary sam)
        'ref:s'  => \$ref, #reference file in fasta  
        'orfs:s' => \$orfs, #start stop position for each gene labelled the same way as the ref file, keep in mind that a gene may code for multiple proteins
        "stub:s" => \$stub,
           );

if (($help)&&!($help)||!($bam)||!($ref)){
  print " _____  _                    _ _______          _           \n";
  print "|  __ \\(_)                  (_)__   __|        | |          \n";
  print "| |  | |___   _____ _ __ ___ _   | | ___   ___ | |___       \n";
  print "| |  | | \\ \\ / / _ \\ '__/ __| |  | |/ _ \\ / _ \\| / __|  \n";
  print "| |__| | |\\ V /  __/ |  \\__ \\ |  | | (_) | (_) | \\__ \\ \n";
  print "|_____/|_| \\_/ \\___|_|  |___/_|  |_|\\___/ \\___/|_|___/  \n\n";
  print "Usage: perl diversiutils.pl -bam input.bam -ref ref.fa -out stub \n";
  print "Usage for the mac compiled version: diversiutils_macosx -bam input.bam -ref ref.fa -out stub \n";
  print "Usage for the linux compiled version: diversiutils_linux -bam input.bam -ref ref.fa -out stub \n";
  print " -bam <txt>  - the input a file in bam format\n";
  print " -ref <txt>  - the reference fasta file\n";
  print " -orfs <txt> - text tab-delimited file with Protein,Beginning,End,Reference(Chr) of the coding sequence [optional: only if you want information about dN/dS and aa frequencies etc...]\n";
  print " -stub <txt> - the output in text-tab delimited\n";
  print " -help        - Get this help\n";
  exit();
}

open(LOG,">$stub\_log.txt")||die "Can't open output $stub\_log.txt\n";

# remove indexed reference fasta file
my $indexref="$ref\.fai";
if ($indexref){
  system("rm $indexref");
}
# high level API
my $sam = Bio::DB::Sam->new(-bam  => $bam, -fasta=> $ref);
my @targets    = $sam->seq_ids;
my %targetHash = map { $_ => 1 } @targets;
my $ins_cnt=0;

# get a hash of the reference sequences with gene, site, refbase
# check whether the references are in the bam file
my %refseq;
my $seq_in  = Bio::SeqIO->new(-format => 'fasta',-file   => $ref);
while( my $seq = $seq_in->next_seq() ) {
  my $id=$seq->display_id();
  if ($targetHash{$id}){
    my $seq_str=$seq->seq();
    my @refbases=split(//,$seq_str);
    for (my $i=0; $i<scalar(@refbases); $i++){
      my $site=1+$i;
      if ($refbases[$i]=~/A|C|T|G|N/i){
        $refseq{$id}{$site}=uc($refbases[$i]);
      }else{
        $refseq{$id}{$site}="N";
      }
    }
  }else{
    print "$id is not in the bam file\n";
  }
}

# if an orfs file is specified, read in the coding regions and load the IUPAC info and codon ambiguity code
my (%codreg, %IUPAC, %c2p, @proteins);
if ($orfs){
    open (CODING,"<$orfs")||die "Can't open $orfs\n";
    my $firstLine = <CODING>; 
    if ($firstLine !~/Protein\tBeg\tEnd\tReference/){
      die "Error: The input file $orfs does not have the proper HEADER:\nProtein  Beg  End  Reference\n";
    }else{
      while(<CODING>){
        chomp($_);
        my @elements=split(/\t/,$_);
        push(@proteins,$elements[0]);
        $codreg{$elements[3]}{$elements[0]}{"Beg"}=$elements[1];#$codreg{Chr name}{ProteinName}{"Beg"}
        $codreg{$elements[3]}{$elements[0]}{"End"}=$elements[2];
        if ($elements[1]<$elements[2]){
          my $codlength=$elements[2]-$elements[1]+1;
          my $mod = $codlength % 3;
          if ($mod ne 0){
            die "Error: The coding sequence for $elements[0] is not a multiple of 3\n";
          }
        }elsif ($elements[2]<$elements[1]){
          my $codlength=$elements[2]-$elements[1]-1;
          my $mod = $codlength % 3;
          if ($mod ne 0){
            die "Error: The coding sequence for $elements[0] is not a multiple of 3\n";
          }
        }
      }
    }
    # from PAL2NAL
    #load translation table if start and stop specified
    # IUPAC
    # W  weak  A      T  
    # S  strong    C  G  
    # M  amino  A  C    
    # K  keto      G  T
    # R  purine  A    G  
    # Y  pyrimidine    C    T
    # B  not A (B comes after A)    C  G  T  
    # D  not C (D comes after C)  A    G  T
    # H  not G (H comes after G)  A  C    T
    # V  not T (V comes after T and U)  A  C  G  
    # N or -  any base (not a gap)  A  C  G  T  
    $IUPAC{"W"} = "T A";
    $IUPAC{"S"} = "C G";
    $IUPAC{"M"} = "A C";
    $IUPAC{"K"} = "G T";
    $IUPAC{"R"} = "A G";
    $IUPAC{"Y"} = "C T";
    $IUPAC{"B"} = "C G T";
    $IUPAC{"D"} = "A G T";
    $IUPAC{"H"} = "A C T";
    $IUPAC{"V"} = "A C G";
    $IUPAC{"N"} = "A C G T";
    # Ambiguous Amino Acids  3-Letter  1-Letter
    # Asparagine or aspartic acid  Asx  B
    # Glutamine or glutamic acid  Glx  Z
    # Leucine or Isoleucine  Xle  J
    # Unspecified or unknown amino acid  Xaa  X
    #--------------------------#
    $c2p{$_} = "Z" for qw(SAA);
    $c2p{$_} = "J" for qw(MTT MTA);
    $c2p{$_} = "B" for qw(AAT AAC GAC GAT RAT RAC);
    $c2p{$_} = "L" for qw(CTA CTT CTG CTC CTN CTK TTA TTG TTR YTA YTG);
    $c2p{$_} = "R" for qw(CGA CGT CGC CGG CGN AGG AGA AGR MGA MGG);
    $c2p{$_} = "S" for qw(TCA TCG TCT TCC TCN AGT AGC AGY);
    $c2p{$_} = "A" for qw(GCC GCT GCA GCG GCN);
    $c2p{$_} = "G" for qw(GGC GGT GGA GGG GGN);
    $c2p{$_} = "P" for qw(CCA CCT CCG CCC CCN);
    $c2p{$_} = "T" for qw(ACA ACG ACC ACT ACN);
    $c2p{$_} = "V" for qw(GTA GTC GTG GTT GTN);
    $c2p{$_} = "I" for qw(ATT ATC ATY ATA ATW);
    #$c2p{$_} = "_" for qw(TAA TAG TAR TGA);
    $c2p{$_} = "*" for qw(TAA TAG TAR TGA);
    $c2p{$_} = "C" for qw(TGT TGC TGY);
    $c2p{$_} = "D" for qw(GAT GAC GAY);
    $c2p{$_} = "E" for qw(GAA GAG GAR);
    $c2p{$_} = "F" for qw(TTT TTC TTY);
    $c2p{$_} = "H" for qw(CAT CAC CAY);
    $c2p{$_} = "K" for qw(AAA AAG AAR);
    $c2p{$_} = "N" for qw(AAT AAC AAY);
    $c2p{$_} = "Q" for qw(CAA CAG CAR);
    $c2p{$_} = "Y" for qw(TAT TAC TAY);
    $c2p{$_} = "M" for qw(ATG);
    $c2p{$_} = "W" for qw(TGG);
    $c2p{$_} = "X" for qw(... AAN AGN ANA ANC ANG ANN ANT ATN CAN CNA CNC CNG CNN CNT GAN GNA GNC GNG GNN GNT NAA NAC NAG NAN NAT NCA NCC NCG NCN NCT NGA NGC NGG NGN NGT NNA NNC NNG NNT NTA NTC NTG NTN NTT TAN TGN TNA TNC TNG TNN TNT TTN);#TYA CAM
    foreach my $target (@targets){
      if (!keys %{$codreg{$target}}){
        my @orfIDs=keys %codreg;
        die "Error: The reference identifer $target in the bam file does not match the reference identifier in your $orfs:\n @orfIDs\n";
      }
    }
}

# start parsing the information from the bam file
my $increment; #use increment for the progress bar
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
      # in the reference sequence substitute all non-ATGC for N
      ($ref_dna=uc($ref_dna))=~s/[^ACTG]/N/g;
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
              $basefreq{$target}{$site}{$strand}{uc($bases[$i])}++;
              #print "$bam\t$target\t$site\t$bases[$i]\t$basefreq{$bam}{$target}{$site}{$bases[$i]}\n"; 
              $refbase{$target}{$site}=$refbases[$i];
              my $P = 10**(-$scores[$i]/10);
              if ($cumulqual{$target}{$site}){
                $cumulqual{$target}{$site}=$cumulqual{$bam}{$target}{$site}+$P;
              }else{
                $cumulqual{$target}{$site}=$P;
              }
              # create a hash for the read information (position of mismatches relative to the start position of a read and motifs upstream and downstream of a mismatch)
              $readinfo{$readpos+1}{"AvQual"}=$readinfo{$readpos+1}{"AvQual"}+$P;
              if ($refbases[$i]=~/$bases[$i]|N/i){ # adding the case where there is an N in the reference 
                $readinfo{$readpos+1}{"CntRef"}++;
                $readinfo{$readpos+1}{"AvQualRef"}=$readinfo{$readpos+1}{"AvQualRef"}+$P;
              }elsif ($refbases[$i]!~/$bases[$i]/i){ 
                $mismatches++;
                $readinfo{$readpos+1}{"CntNonRef"}++;
                $readinfo{$readpos+1}{"AvQualNonRef"}=$readinfo{$readpos+1}{"AvQualNonRef"}+$P;
                # motif on either side of mismatch
                # not interested in this at this stage
#                 my ($motif_start,$motif_end,$endfromref,$startfromref);
#                 if (length($query_dna)<=($i+10)){
#                   $motif_end=length($query_dna);
#                   $endfromref=$i+10+1-length($query_dna);
#                 }elsif (length($query_dna)>($i+10)){
#                   $motif_end=$i+10;
#                 }if ($i-10<0){
#                   $motif_start=0;
#                   $startfromref=10-$i;
#                 }elsif($i-10>=0){
#                   $motif_start=$i-10;
#                 }if ($start-1+$i-10>=0 && length($refseq)>=$i+10){
#                   my $motif = substr ($query_dna,$motif_start,$motif_end+1-$motif_start);
#                   my $begmotif = lc(substr ($refseq,($start-1+$i-10),$startfromref));
#                   my $endmotif = lc(substr ($refseq,(($start-1+$i-10+21)-$endfromref),$endfromref));
#                   $motif = $begmotif.$motif.$endmotif;
#                   # print "Motif $begmotif$motif$endmotif\n";
#                   my $refmotif = substr ($refseq,($start-1+$i-10),21);
#                   # print "Refer $refmotif\n\n";
#                   if (length($motif)<21){
#                     #the motif is too close to the beginning or end of the reference sequence
#                     #print "Motifend $motif_end EndFromRef $endfromref\t StartMotif $motif_start StartFromRef $startfromref\tlength ".length($query_dna)." position $i Length of reference ".length($refseq)." site $site\n";
#                   }elsif (length($motif)==21){
#                     print MOTIF ">$name\_$mismatches\n$motif\n";
#                   }
#                }
              }elsif ($bases[$i]=~/N/i){
                $readinfo{$readpos+1}{"CntNs"}++;
              }
              if ($orfs){
                #$codreg{Chr name}{ProteinName}{"Beg"}
                foreach my $prot (keys %{$codreg{$target}}){
                  #check here that Beg<End
                  if ($codreg{$target}{$prot}{"Beg"}<$codreg{$target}{$prot}{"End"}){
					  my $noUTR=($site-$codreg{$target}{$prot}{"Beg"})+1;
					  my ($aasite,$mod);
						if ($noUTR==1){
						  $mod=1;
						}elsif ($noUTR>1){
						  $mod = $noUTR % 3;
						}
						# check that the site is in a coding region, that it is the first codon position of the coding region and that the read is long enough for final codon
						if ($site>=$codreg{$target}{$prot}{"Beg"} && $site<=$codreg{$target}{$prot}{"End"} && $mod==1 && $i<(scalar(@bases)-2)){
						  #print "$prot\t$target\tCoding region for $prot Site $site and Modular $mod\n";
						  my $rcodon=$refbases[$i].$refbases[$i+1].$refbases[$i+2];
						  my $qcodon=$bases[$i].$bases[$i+1].$bases[$i+2];
						  #ignore AA if 
						  my $raa= $c2p{uc($rcodon)};
						  my $qaa= $c2p{uc($qcodon)};
						  if (!$qaa){
						    print "Error: No translation for codon $qcodon\n"
						  }
						  #print "$rcodon\t$qcodon\n";
						  #nucsite aasite $rcodon $raa $qcodon $qaa $codonposmis
						  if ($noUTR==1){
							$aasite=1;
						  }else{
							$aasite = ($noUTR+3-$mod)/3;
							#print "No UTR $noUTR modular $mod $aasite\n";
						  }
						  # the position of the mutation
						  # Sample\tChr\tAAPosition\tRefAA\tRefCodon\tCntNonSyn\tCntSyn\tTopAA\tTopAAcnt\tSndAA\tSndAAcnt\tTrdAA\tTrdAAcnt\tAAcoverage\tNbStop
						  $aafreq{$target}{$prot}{$aasite}{"AAcoverage"}++;#this will provide the coverage
						  #print $aasite."\t".$aafreq{$target}{$prot}{$aasite}{"AAcoverage"}."\t".$aafreq{$target}{$prot}{$aasite}{"RefAA"}."\t$target\t$prot\t$rcodon\t$raa\t$qcodon\t$qaa\n";
						  $aaorder{$target}{$prot}{$aasite}{$qaa}++;
						  $codonfreq{$target}{$prot}{$aasite}{$qcodon}++;
						  $aafreq{$target}{$prot}{$aasite}{"RefAA"}=$raa;
						  #print $aafreq{$target}{$prot}{$aasite}{"RefAA"}."\t>$target<\t>$prot<>$aasite<\n";
						  $aafreq{$target}{$prot}{$aasite}{"RefCodon"}=$rcodon;
						  $aafreq{$target}{$prot}{$aasite}{"RefSite"}=$site;
						  my $aamut=$raa.$aasite.$qaa;
						  if (uc($rcodon) ne uc($qcodon)){
						  # the syn and non-syn counts only refer to codons where there are changes
						  # i.e. when the codon is the same as the reference, this is not counted as being a non-synonymous count
						  # the count for syn and non-syn is the number of changes relative to the reference codon, this can range from 1-3
							if(uc($raa) eq uc($qaa)){
							  my $mm=mismatch_count($rcodon,$qcodon);
							  $aafreq{$target}{$prot}{$aasite}{"syn"}+=$mm;
							}if(uc($raa) ne uc($qaa)){
							  my $mm=mismatch_count($rcodon,$qcodon);
							  $aafreq{$target}{$prot}{$aasite}{"nonsyn"}+=$mm;
							}if(uc($qaa)=~/\*/){
							  $stopfreq{$target}{$prot}{$aasite}{$strand}++;
							  $aafreq{$target}{$prot}{$aasite}{"stop"}++;
							}if($refbases[$i] ne $bases[$i]){
							  $aafreq{$target}{$prot}{$aasite}{"firstcodonpos"}++;
							}if($refbases[$i+1] ne $bases[$i+1]){
							  $aafreq{$target}{$prot}{$aasite}{"secondcodonpos"}++;
							}if($refbases[$i+2] ne $bases[$i+2]){
							  $aafreq{$target}{$prot}{$aasite}{"thirdcodonpos"}++;
							}
						  }
						}
                    }if ($codreg{$target}{$prot}{"Beg"}>$codreg{$target}{$prot}{"End"}){# CDS on the reverse strand
                      
                      my $noUTR=($site-$codreg{$target}{$prot}{"End"})+1;
					  my ($aasite,$mod);
						if ($noUTR==1){
						  $mod=1;
						}elsif ($noUTR>1){
						  $mod = $noUTR % 3;
						}
						# check that the site is in a coding region, that it is the first codon position of the coding region and that the read is long enough for final codon
						if ($site>=$codreg{$target}{$prot}{"End"} && $site<=$codreg{$target}{$prot}{"Beg"} && $mod==1 && $i<(scalar(@bases)-2)){
						  #print "$target Proteins $prot in the negative strand ".$codreg{$target}{$prot}{"Beg"}." - ".$codreg{$target}{$prot}{"End"}." have not been fully tested Site $site Position in read $i \n";
						  my $rcodon=&revcomp($refbases[$i].$refbases[$i+1].$refbases[$i+2]);
						  my $qcodon=&revcomp($bases[$i].$bases[$i+1].$bases[$i+2]);
						  my $raa= $c2p{uc($rcodon)};
						  my $qaa= $c2p{uc($qcodon)};
						  if (!$qaa){
						    print "Error: No translation for codon $qcodon\n"
						  }						  
						  #print "$rcodon\t$qcodon\n";
						  #nucsite aasite $rcodon $raa $qcodon $qaa $codonposmis
						  if ($noUTR==1){#this corresponds to the last amino acid in the sequence
							$aasite=($codreg{$target}{$prot}{"Beg"}-$codreg{$target}{$prot}{"End"}-1)/3;
							#print "Last site $aasite\n";
						  }else{# if it is not the last amino acid then we calculate the aa position based on the distance of the site from the beginning of the protein
							$aasite = ($codreg{$target}{$prot}{"Beg"}-$site+1) / 3;
							#print "No UTR $noUTR modular $mod Nucsite ".$site+2." AAsite $aasite\n";
						  }
						  $aafreq{$target}{$prot}{$aasite}{"AAcoverage"}++;#this will provide the coverage
						  $aaorder{$target}{$prot}{$aasite}{$qaa}++;
						  $codonfreq{$target}{$prot}{$aasite}{$qcodon}++;
						  $aafreq{$target}{$prot}{$aasite}{"RefAA"}=$raa;
						  $aafreq{$target}{$prot}{$aasite}{"RefCodon"}=$rcodon;
						  $aafreq{$target}{$prot}{$aasite}{"RefSite"}=$site+2;# plus 2 because we want the first position of the codon and the CDS is in reverse direction
						  my $aamut=$raa.$aasite.$qaa;
						  if (uc($rcodon) ne uc($qcodon)){
						  # the syn and non-syn counts only refer to codons where there are changes
						  # i.e. when the codon is the same as the reference, this is not counted as being a non-synonymous count
						  # the count for syn and non-syn is the number of changes relative to the reference codon, this can range from 1-3
							if(uc($raa) eq uc($qaa)){
							  my $mm=mismatch_count($rcodon,$qcodon);
							  $aafreq{$target}{$prot}{$aasite}{"syn"}+=$mm;
							}if(uc($raa) ne uc($qaa)){
							  my $mm=mismatch_count($rcodon,$qcodon);
							  $aafreq{$target}{$prot}{$aasite}{"nonsyn"}+=$mm;
							}if(uc($qaa)=~/\*/){
							  $stopfreq{$target}{$prot}{$aasite}{$strand}++;
							  $aafreq{$target}{$prot}{$aasite}{"stop"}++;
							}if($refbases[$i] ne $bases[$i]){
							  $aafreq{$target}{$prot}{$aasite}{"firstcodonpos"}++;
							}if($refbases[$i+1] ne $bases[$i+1]){
							  $aafreq{$target}{$prot}{$aasite}{"secondcodonpos"}++;
							}if($refbases[$i+2] ne $bases[$i+2]){
							  $aafreq{$target}{$prot}{$aasite}{"thirdcodonpos"}++;
							}
						  }
						}
                    }
                  }#for each protein loop
                }# end of orfs IF statement
                $site=$site+1;
                $readpos=$readpos+1;
                $refpos=$refpos+1;
                #print "$name $k\t$cigars[$k] $site $readpos $refpos => match\n";
              }
            $readmis{$mismatches}++;
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
            $insfreq{$target}{$site}{$sub_query_dna}++;
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
print "\n";
print LOG "$bam:\n";
# print out the number of reads with 1 , 2, 3, 4, etc.. mismatches
print LOG "Frequency of mismatches per read (mismatches: number of reads)\n";
for my $mismatchcnt (sort {$a<=>$b} keys %readmis){
  print LOG "$mismatchcnt: $readmis{$mismatchcnt}\n";
} 
my %shannon;#Shannon-Wiener Diversity Index (also known as Shannon's diversity index, the Shannon-Wiener index, the Shannon-Weaver index and the Shannon entropy)
my @nuc=qw/A C T G/;
open (OUT, ">$stub\_entropy.txt")||die "can't open $stub\_entropy.txt\n";
# add the reference base in the table
# Nucleotide Table:
# Sample\tChr\tPosition\tRefBase\tCoverage\tAvQual\tAcnt\tApval\tCcnt\tCpval\tTcnt\tTpval\tGcnt\tGpval\tentropy(base e)\tNonRefCnt
# CntTv\tCntTs\tOrderATCG\tstrandbias\tins_cnt\tins_mode\tdel_cnt\tdel_mode\n
# strandbias as A3:5C5:6G4:5T5:7

print OUT "Sample\tChr\tPosition\tRefBase\tCoverage\tAvQual\tAcnt\tApval\tCcnt\tCpval\tTcnt\tTpval\tGcnt\tGpval\tentropy(base e)\tNonRefCnt\tCntTv\tCntTs\tOrderOfNucs\tstrandbias\tins_cnt\tins_mode\tdel_cnt\tdel_mode\tNcnt\n";
foreach my $gene (keys %refseq){
  my $nbsites=0;
  my $sumentropy=0;
  my $rawsumentropy=0;
  foreach my $site (sort {$a<=>$b} keys %{$refseq{$gene}}){
    #print "$site\n";
    if (keys %{$basefreq{$gene}{$site}}){# if there is information in the bam file about this site
      # strandbias info
      my ($Aplus,$Aneg,$Cplus,$Cneg,$Tplus,$Tneg,$Gplus,$Gneg,$Nplus,$Nneg);
      if ($basefreq{$gene}{$site}{1}{"A"}){ $Aplus=$basefreq{$gene}{$site}{1}{"A"}}else{$Aplus=0}
      if ($basefreq{$gene}{$site}{-1}{"A"}){ $Aneg=$basefreq{$gene}{$site}{-1}{"A"}}else{$Aneg=0}
      if ($basefreq{$gene}{$site}{1}{"C"}){ $Cplus=$basefreq{$gene}{$site}{1}{"C"}}else{$Cplus=0}
      if ($basefreq{$gene}{$site}{-1}{"C"}){ $Cneg=$basefreq{$gene}{$site}{-1}{"C"}}else{$Cneg=0}
      if ($basefreq{$gene}{$site}{1}{"T"}){ $Tplus=$basefreq{$gene}{$site}{1}{"T"}}else{$Tplus=0}
      if ($basefreq{$gene}{$site}{-1}{"T"}){ $Tneg=$basefreq{$gene}{$site}{-1}{"T"}}else{$Tneg=0}
      if ($basefreq{$gene}{$site}{1}{"G"}){ $Gplus=$basefreq{$gene}{$site}{1}{"G"}}else{$Gplus=0}
      if ($basefreq{$gene}{$site}{-1}{"G"}){ $Gneg=$basefreq{$gene}{$site}{-1}{"G"}}else{$Gneg=0}
      if ($basefreq{$gene}{$site}{1}{"N"}){ $Nplus=$basefreq{$gene}{$site}{1}{"N"}}else{$Nplus=0}
      if ($basefreq{$gene}{$site}{-1}{"N"}){ $Nneg=$basefreq{$gene}{$site}{-1}{"N"}}else{$Nneg=0}

      my $strandbias="A$Aplus".":$Aneg;"."C$Cplus".":$Cneg;"."T$Tplus".":$Tneg;"."G$Gplus".":$Gneg";

      my ($cntA,$cntC,$cntT,$cntG,$cntN);
      $cntA=$Aplus+$Aneg;
      $cntC=$Cplus+$Cneg;
      $cntT=$Tplus+$Tneg;
      $cntG=$Gplus+$Gneg;
      $cntN=$Nplus+$Nneg;
      my $coverage = $cntA + $cntT + $cntC + $cntG + $cntN;
      my $average_p=$cumulqual{$gene}{$site}/$coverage;
      my $p = $average_p/3;
      my %prob;
      $prob{"A"} = 1 - (&Math::CDF::pbinom(($cntA-1), $coverage, $p));# need to double check with Richard about the -1
      $prob{"C"} = 1 - (&Math::CDF::pbinom(($cntC-1), $coverage, $p));
      $prob{"T"} = 1 - (&Math::CDF::pbinom(($cntT-1), $coverage, $p));
      $prob{"G"} = 1 - (&Math::CDF::pbinom(($cntG-1), $coverage, $p));
      if ($coverage>0){
        foreach my $nuc (@nuc){
          my $nucnt=$basefreq{$gene}{$site}{1}{$nuc}+$basefreq{$gene}{$site}{-1}{$nuc};
          my $p = $nucnt / $coverage;
          if($nucnt > 0){
            $shannon{$gene}{$site} += -$p*log($p);#natural log, i.e. log base e
          }
        }   
      }else{
        $shannon{$gene}{$site}="<NA>";
      }   
      $nbsites++;
      $sumentropy=$sumentropy+$shannon{$gene}{$site};
      my $refbase=uc($refbase{$gene}{$site});  # modified to upper case 2017-02-06
      my $nonrefcnt;
      if ($refbase=~/N/i){
        $nonrefcnt=$coverage-($basefreq{$gene}{$site}{1}{"N"})-($basefreq{$gene}{$site}{-1}{"N"}); # dealing with issues #11
      }else{
        $nonrefcnt=$coverage-($basefreq{$gene}{$site}{1}{$refbase})-($basefreq{$gene}{$site}{-1}{$refbase}); 
      }
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

      #print "$site $refbase $NucOrder\n";	
      # get the cnt of insertion and mode and for deletions
      my ($del_cnt,$ins_cnt,$freq_ins,$freq_del);

      if (keys %{$delfreq{$gene}{$site}}){ $del_cnt= sum values %{$delfreq{$gene}{$site}}}else{$del_cnt="<NA>"}
      if (keys %{$insfreq{$gene}{$site}}){ $ins_cnt= sum values %{$insfreq{$gene}{$site}}}else{$ins_cnt="<NA>"}
      # most frequently found insertion and deletion
      if (keys %{$insfreq{$gene}{$site}}){$freq_ins = (sort {$insfreq{$gene}{$site}{$b} <=> $insfreq{$gene}{$site}{$a}} keys %{$insfreq{$gene}{$site}})[0]}else{$freq_ins="<NA>"};
      if (keys %{$delfreq{$gene}{$site}}){$freq_del = (sort {$delfreq{$gene}{$site}{$b} <=> $delfreq{$gene}{$site}{$a}} keys %{$delfreq{$gene}{$site}})[0]}else{$freq_del ="<NA>"};
      print OUT "$bam\t$gene\t$site\t".uc($refbase)."\t$coverage\t$average_p\t$cntA\t".$prob{"A"}."\t$cntC\t".$prob{"C"}."\t$cntT\t".$prob{"T"}."\t$cntG\t".$prob{"G"}."\t";
      print OUT "$shannon{$gene}{$site}\t$nonrefcnt\t$Ts\t$Tv\t$NucOrder\t$strandbias\t$ins_cnt\t$freq_ins\t$del_cnt\t$freq_del\t$cntN\n";
    }else{#there is no coverage for that site in the bam
        #print "No coverage $site\n";
        my ($del_cnt,$ins_cnt,$freq_ins,$freq_del);# this is for printing the insertion or deletion when there is zero coverage at a site but 100% deletion or 100% insertion or 100% Ns
        if (keys %{$delfreq{$gene}{$site}}){ $del_cnt= sum values %{$delfreq{$gene}{$site}}}else{$del_cnt="<NA>"}
        if (keys %{$insfreq{$gene}{$site}}){ $ins_cnt= sum values %{$insfreq{$gene}{$site}}}else{$ins_cnt="<NA>"}
        # most frequently found insertion and deletion
        if (keys %{$insfreq{$gene}{$site}}){$freq_ins = (sort {$insfreq{$gene}{$site}{$b} <=> $insfreq{$gene}{$site}{$a}} keys %{$insfreq{$gene}{$site}})[0]}else{$freq_ins="<NA>"};
        if (keys %{$delfreq{$gene}{$site}}){$freq_del = (sort {$delfreq{$gene}{$site}{$b} <=> $delfreq{$gene}{$site}{$a}} keys %{$delfreq{$gene}{$site}})[0]}else{$freq_del ="<NA>"};
        my ($Nplus,$Nneg,$cntN);
        if ($basefreq{$gene}{$site}{1}{"N"}){ $Nplus=$basefreq{$gene}{$site}{1}{"N"}}else{$Nplus=0}
        if ($basefreq{$gene}{$site}{-1}{"N"}){ $Nneg=$basefreq{$gene}{$site}{-1}{"N"}}else{$Nneg=0}
        $cntN=$Nplus+$Nneg;
        print OUT "$bam\t$gene\t$site\t".uc($refseq{$gene}{$site})."\t0\t<NA>\t<NA>\t<NA>\t<NA>\t<NA>\t<NA>\t<NA>\t<NA>\t<NA>\t";
        print OUT "<NA>\t<NA>\t<NA>\t<NA>\t<NA>\t<NA>";
        if ($ins_cnt){print OUT "\t$ins_cnt";}else{print OUT "\t<NA>";}
        if ($freq_ins){print OUT "\t$freq_ins";}else{print OUT "\t<NA>";}
        if ($del_cnt){print OUT "\t$del_cnt";}else{print OUT "\t<NA>";}
        if ($freq_del){print OUT "\t$freq_del";}else{print OUT "\t<NA>";}
        if ($cntN){print OUT "\t$cntN";}else{print OUT "\t<NA>";}
        print OUT "\n";
    }
  }
  if ($nbsites>0){
    print LOG "Gene $gene Average entropy = ".$sumentropy/$nbsites."\n";
  }else{
    print LOG "Gene $gene Average entropy = Not Available\n";
  }
}
# Read mismatch table:
# ReadPos\tCntNonRef\tCntRef\tTotalCnt\tFreq\tAvQual (Freq=NonRef/TotalCnt)
open (READ, ">$stub\_read.txt")||die "can't open $stub\_read.txt\n";
print READ "ReadPos\tCntRef\tCntNonRef\tTotalCnt\tFreq\tAvQual\tAvQualRef\tAvQualNonRef\n";
foreach my $readpos (sort {$a<=>$b} keys %readinfo){
  print READ "$readpos\t";
  my $total = $readinfo{$readpos}{"CntNonRef"}+$readinfo{$readpos}{"CntRef"};
  if ($total>0){
    if ($readinfo{$readpos}{"CntRef"}=~/.+/){ print READ $readinfo{$readpos}{"CntRef"}."\t";}else{ print READ "0\t"}
    if ($readinfo{$readpos}{"CntNonRef"}=~/.+/){ print READ $readinfo{$readpos}{"CntNonRef"}."\t";}else{ print READ "0\t"}
    print READ $total."\t".$readinfo{$readpos}{"CntNonRef"}/$total."\t";
    if ($readinfo{$readpos}{"AvQual"}=~/.+/ & $total>0){ print READ $readinfo{$readpos}{"AvQual"}/$total."\t";}else{print READ "0\t"}
    if ($readinfo{$readpos}{"CntRef"}=~/.+/){print READ $readinfo{$readpos}{"AvQualRef"}/$readinfo{$readpos}{"CntRef"}."\t";}else{print READ "0\t";}
    if ($readinfo{$readpos}{"CntNonRef"}=~/.+/){print READ $readinfo{$readpos}{"AvQualNonRef"}/$readinfo{$readpos}{"CntNonRef"}."\n";}else{print READ "0\n"}
  }else{
    print READ "<NA>\t<NA>\t<NA>\t<NA>\t<NA>\t<NA>\t<NA>\n";  
  }
}
# loop for the aa mutations and position of mismatches in the codon
if ($orfs){
  open (AA, ">$stub\_AA.txt")||die "can't open $stub\_AA.txt\n";
  # the position of the mutation
  print AA "Sample\tChr\tProtein\tAAPosition\tRefAA\tRefSite\tRefCodon\tFstCodonPos\tSndCodonPos\tTrdCodonPos\tCntNonSyn\tCntSyn\tNbStop\tTopAA\tTopAAcnt\tSndAA\tSndAAcnt\tTrdAA\tTrdAAcnt\tTopCodon\tTopCodoncnt\tSndCodon\tSndCodoncnt\tTrdCodon\tTrdCodoncnt\tAAcoverage\n";
    foreach my $target (keys %codreg){
      foreach my $prot (@proteins){  #to provide the output in the same order as the coding region input table
        # dealing with lack of coverage within regions of the protein
        # #$codreg{Chr name}{ProteinName}{"Beg"}
        my $beg=$codreg{$target}{$prot}{"Beg"};
        my $end=$codreg{$target}{$prot}{"End"};
        my $aacnt=abs(($end-$beg)/3);
        #foreach my $aasite (sort {$a<=>$b} keys %{$aafreq{$target}{$prot}}){
        for (my $aasite=1; $aasite<=$aacnt; $aasite++){
         if (defined $aafreq{$target}{$prot}{$aasite}{"RefAA"}){
          print AA "$bam\t$target\t$prot\t$aasite\t";
          
          print AA $aafreq{$target}{$prot}{$aasite}{"RefAA"}."\t";
          print AA $aafreq{$target}{$prot}{$aasite}{"RefSite"}."\t";
          print AA $aafreq{$target}{$prot}{$aasite}{"RefCodon"}."\t";
          if ($aafreq{$target}{$prot}{$aasite}{"firstcodonpos"}=~/\d+/){ print AA $aafreq{$target}{$prot}{$aasite}{"firstcodonpos"}."\t";}else{ print AA "<NA>\t"}
          if ($aafreq{$target}{$prot}{$aasite}{"secondcodonpos"}=~/\d+/){ print AA $aafreq{$target}{$prot}{$aasite}{"secondcodonpos"}."\t";}else{ print AA "<NA>\t"}
          if ($aafreq{$target}{$prot}{$aasite}{"thirdcodonpos"}=~/\d+/){ print AA $aafreq{$target}{$prot}{$aasite}{"thirdcodonpos"}."\t";}else{ print AA "<NA>\t"}
          if ($aafreq{$target}{$prot}{$aasite}{"nonsyn"}=~/\d+/){ print AA $aafreq{$target}{$prot}{$aasite}{"nonsyn"}."\t";}else{ print AA "<NA>\t"}
          if ($aafreq{$target}{$prot}{$aasite}{"syn"}=~/\d+/){ print AA $aafreq{$target}{$prot}{$aasite}{"syn"}."\t";}else{ print AA "<NA>\t"}
          my ($stop_plus,$stop_neg);
          if ($stopfreq{$target}{$prot}{$aasite}{1}){ $stop_plus = $stopfreq{$target}{$prot}{$aasite}{1}}else{$stop_plus="0"}
          if ($stopfreq{$target}{$prot}{$aasite}{-1}){ $stop_neg = $stopfreq{$target}{$prot}{$aasite}{-1}}else{$stop_neg="0"}
          
          #print AA $stopfreq{$target}{$prot}{$aasite}{1}.":".$stopfreq{$target}{$prot}{$aasite}{-1}."(".$aafreq{$target}{$prot}{$aasite}{"stop"}.")"."\t";
          print AA $stop_plus.";".$stop_neg."\t";
          my $topAAs=0;
          foreach my $aa (sort { $aaorder{$target}{$prot}{$aasite}{$b} <=> $aaorder{$target}{$prot}{$aasite}{$a} } keys %{$aaorder{$target}{$prot}{$aasite}}) {
            if ($topAAs<3){# provide the top three AAs and their counts
              print AA "$aa\t$aaorder{$target}{$prot}{$aasite}{$aa}\t";
              $topAAs++;
            }
          }          
          while ($topAAs<3){# if there are less than three different AAs
            print AA "<NA>\t<NA>\t";
            $topAAs++;
          }
          my $topCodons=0;
          foreach my $codon (sort { $codonfreq{$target}{$prot}{$aasite}{$b} <=> $codonfreq{$target}{$prot}{$aasite}{$a} } keys %{$codonfreq{$target}{$prot}{$aasite}}) {
            if ($topCodons<3){# provide the top three AAs and their counts
              print AA "$codon\t$codonfreq{$target}{$prot}{$aasite}{$codon}\t";
              $topCodons++;
            }
          }
          while ($topCodons<3){# if there are less than three different codons
            print AA "<NA>\t<NA>\t";
            $topCodons++;
          }
          print AA $aafreq{$target}{$prot}{$aasite}{"AAcoverage"}."\n";
        }else{# close if statement
          print AA "$bam\t$target\t$prot\t$aasite\t<NA>\t".$aafreq{$target}{$prot}{$aasite}{"RefSite"};
          print AA "\t<NA>\t<NA>\t<NA>\t<NA>\t<NA>\t<NA>\t<NA>\t<NA>\t<NA>\t<NA>\t<NA>\t<NA>\t<NA>";
          print AA "\t<NA>\t<NA>\t<NA>\t<NA>\t<NA>\t<NA>\t".$aafreq{$target}{$prot}{$aasite}{"AAcoverage"}."\n";
        }
        }
      }
    }
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

# count mismatches between two strings
sub mismatch_count($$) { length( $_[ 0 ] ) - ( ( $_[ 0 ] ^ $_[ 1 ] ) =~ tr[\0][\0] ) }

# fastest way to get the key with the largest value but uses more memory 
sub largest_value_mem (\%) {
    my $hash   = shift;
    my ($key, @keys) = keys   %$hash;
    my ($big, @vals) = values %$hash;

    for (0 .. $#keys) {
        if ($vals[$_] > $big) {
            $big = $vals[$_];
            $key = $keys[$_];
        }
    }
    $key
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

sub revcomp {
        my $dna = shift;

        # reverse the DNA sequence
        my $revcomp = reverse($dna);

        # complement the reversed DNA sequence
        $revcomp =~ tr/ATGCatgcNn/TACGtacgNn/;
        return $revcomp;
}