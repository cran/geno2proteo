#!/usr/bin/perl -w

use strict;

my $len = $#ARGV;

my $inputFiles = $ARGV[0]; #the text file containing the exons and their DNA seq

my $geneticCodesFile = $ARGV[1];

my $outputFile = $ARGV[2]; #the output file


## first get the genetic codes

 my %geneticCodes;
 
 open DOC, $geneticCodesFile or die "Cannot open the file$geneticCodesFile\n";
 my $num=1;
 while(<DOC>) {
   chomp;
   my @items = split/[\t ]+/;
   #$items[1] = 9 if($items[1] eq '*');
   $geneticCodes{$items[0]} = $items[1];
   #print "$num ($items[0]) ($items[1])\n";
   ++$num;
 }
 close DOC;

 # then get the CDS
 open DOC, $inputFiles or die "Cannot open the file $inputFiles\n";
 
 open outDOC, '>'.$outputFile or die "Cannot open the file $outputFile\n";
 
 my $prevTrans=' ';
 my $currStrand=' ';
 my @exonsInA='';
 my @exonsSeqInA='';
 my @frameInA='';
 my $numExons=0;
 
 while(<DOC>) {
   chomp;
   
   my @items = split/[\t ]+/;
   my $lastItem = scalar @items;
   $lastItem -= 1;
   $items[$lastItem] = uc $items[$lastItem]; # the DNA sequence, in either strand, in upper case
   my @items11 = @items[4,5]; # the 1st is the transcript id, and the 2nd is the rank of exon in the transcript.
   
   #print "3=$items[3], 4=$items[4], $prevTrans, $currStrand\n";
  #  print "$_\n";
   # my $iii=<STDIN>;
   if($prevTrans eq ' ') {
      $exonsInA[$numExons]=$_;
      $exonsSeqInA[$numExons]=$items[$lastItem];
      $frameInA[$numExons] = $items[6];
      $currStrand=$items[3];
      $prevTrans=$items11[0];
      ++$numExons;
   } elsif($prevTrans ne $items11[0] ) { # finish a transcript
      # deal with the previous transcript first
     if($currStrand eq '-') { # reverse the order for the - strand.
       @exonsInA = reverse @exonsInA;
       @exonsSeqInA= reverse @exonsSeqInA;
       @frameInA = reverse @frameInA;
      # for(my $i=0; $i<$numExons; ++$i) { # don't need do this, because it was done already.
         #print "before: $exonsSeqInA[$i]\n";
       #  $exonsSeqInA[$i] = reverseComplement($exonsSeqInA[$i]);
         #print "after: $exonsSeqInA[$i]\n";
         #my $iii=<STDIN>;
       #}
     }
     
     #for(my $exonInd=0; $exonInd < $numExons; ++$exonInd) {
     # print "$exonInd: $exonsInA[$exonInd]; $exonsSeqInA[$exonInd]; $frameInA[$exonInd]\n";
     #}
     #my $iii=<STDIN>;
    
     my $exonLen=0;
     my $remainN=0;
     my $exonInd=0;
     my $frameNow = $frameInA[$exonInd];
     my $exonNow=$exonsSeqInA[$exonInd];
     my $addNs = '';
     if($frameNow == 1) { #if the 1st exon is not in frame, add some additional Ns and one AA as X as well.
        $addNs = 'NN';
     }
     if($frameNow == 2) { #if the 1st exon is not in frame, add some additional Ns and one AA as X as well.
        $addNs = 'N';
     }
     $exonNow = $addNs.$exonNow;
     $exonLen = length($exonNow);
     $remainN = $exonLen % 3;
     my $exonSeqTr;
     my $exonPep;
     if($exonLen == $remainN ) {
       $exonPep = ' ';
     } else {
       $exonSeqTr = substr($exonNow, 0, $exonLen-$remainN); # get the major DNAseq
       $exonPep = trans2pep($exonSeqTr); # translate
     }
     
     my $remainDNAs = substr($exonNow, $exonLen-$remainN, $remainN); # the remain dna seq
     
     my $prevRemainDNA=' ';
     if($frameNow !=0) {
       $prevRemainDNA= $addNs;
     }
 
     $exonsInA[$exonInd] .= "\t$prevRemainDNA\t$remainN\t$remainDNAs\t$exonPep";
     
     $prevRemainDNA=$remainDNAs;
     
     for($exonInd=1; $exonInd < $numExons; ++$exonInd) {
       $exonNow=$exonsSeqInA[$exonInd];
       $exonNow = $remainDNAs.$exonNow;  # add the previous remain one if there are some
       $exonLen = length($exonNow);
       $remainN = $exonLen % 3;
       $exonSeqTr = substr($exonNow, 0, $exonLen-$remainN);
    
       $exonPep = trans2pep($exonSeqTr);#translate
       $remainDNAs = substr($exonNow, $exonLen-$remainN, $remainN); #get the remain one if there are some
       $exonsInA[$exonInd] .=  "\t$prevRemainDNA\t$remainN\t$remainDNAs\t$exonPep";
       
       $frameNow = $frameInA[$exonInd];
       my $len = length($prevRemainDNA);
       if( $frameNow==0 && $len !=0 ) {
         print "The frame stated in the annotation file may not be correct, frameNow=$frameNow, len=$len!!\n";
         print "$exonInd, frame=$frameNow, len=$len\n";
         print "previous: $exonsInA[$exonInd-1]\n";
         print "current: $exonsInA[$exonInd]\n";
         #my $iii=<STDIN>;
       }
       if($frameNow>0 && $len != 3- $frameNow) {
         print "The frame stated in the annotation file may not be correct, frameNow=$frameNow, len=$len!!!\n";
         print "$exonInd, frame=$frameNow, len=$len\n";
         print "previous: $exonsInA[$exonInd-1]\n";
         print "current: $exonsInA[$exonInd]\n";
       }
       $prevRemainDNA=$remainDNAs;
       
       #print "$exonInd: $exonsInA[$exonInd]\n";
       #my $iii=<STDIN>;
       
     }
      # write the results into a file
     for( my $i=0; $i < $numExons; ++$i ) {
       #print "$i $exonsInA[$i]\n";
       print outDOC "$exonsInA[$i]\n";
       # my $iii=<STDIN>;
     }
     
     # then get the current exon.
     $numExons=0;
     @exonsInA='';
     @exonsSeqInA='';
     @frameInA='';
     
     $currStrand=$items[3];
     $prevTrans=$items11[0];
     $exonsInA[$numExons]=$_;
     $exonsSeqInA[$numExons]=$items[$lastItem];
     $frameInA[$numExons] = $items[6];
     ++$numExons;
     
     
   } else { # the same transcript as in the last line
      $exonsInA[$numExons]=$_;
      $exonsSeqInA[$numExons]=$items[$lastItem];
      $frameInA[$numExons] = $items[6];
      ++$numExons;
   
   }
   
 }
 
 #the last trascript in the file
 
    if($currStrand eq '-') { # reverse the order for the - strand.
       @exonsInA = reverse @exonsInA;
       @exonsSeqInA= reverse @exonsSeqInA;
       @frameInA = reverse @frameInA;
      # for(my $i=0; $i<$numExons; ++$i) { # don't need do this, because it was done already.
         #print "before: $exonsSeqInA[$i]\n";
       #  $exonsSeqInA[$i] = reverseComplement($exonsSeqInA[$i]);
         #print "after: $exonsSeqInA[$i]\n";
         #my $iii=<STDIN>;
       #}
     }
     
 
 my $exonLen=0;
     my $remainN=0;
     my $exonInd=0;
     my $frameNow = $frameInA[$exonInd];
     my $exonNow=$exonsSeqInA[$exonInd];
     my $addNs = '';
     if($frameNow == 1) { #if the 1st exon is not in frame, add some additional Ns and one AA as well.
        $addNs = 'NN';
     }
     if($frameNow == 2) { #if the 1st exon is not in frame, add some additional Ns and one AA as well.
        $addNs = 'N';
     }
     $exonNow = $addNs.$exonNow;
     $exonLen = length($exonNow);
     $remainN = $exonLen % 3;
     my $exonSeqTr;
     my $exonPep;
     if($exonLen == $remainN ) {
       $exonPep = ' ';
     } else {
       $exonSeqTr = substr($exonNow, 0, $exonLen-$remainN); # get the major DNAseq
       $exonPep = trans2pep($exonSeqTr); # translate
     }
     
     my $remainDNAs = substr($exonNow, $exonLen-$remainN, $remainN); # the remain dna seq
     
     my $prevRemainDNA=' ';
     if($frameNow !=0) {
       $prevRemainDNA= $addNs;
     }
 
     $exonsInA[$exonInd] .= "\t$prevRemainDNA\t$remainN\t$remainDNAs\t$exonPep";
     
     $prevRemainDNA=$remainDNAs;
     
     for($exonInd=1; $exonInd < $numExons; ++$exonInd) {
       $frameNow = $frameInA[$exonInd];
       my $len = length($prevRemainDNA);
       if( $frameNow==0 && $len !=0 ) {
         print "The frame stated in the annotation file may not be correct, frameNow=$frameNow, len=$len!!\n";
         print "$exonInd, frame=$frameNow, len=$len\n";
         print "previous: $exonsInA[$exonInd-1]\n";
         print "current: $exonsInA[$exonInd]\n";
       }
       if($frameNow>0 && $len != 3- $frameNow) {
         print "The frame stated in the annotation file may not be correct, frameNow=$frameNow, len=$len!!!\n";
         print "$exonInd, frame=$frameNow, len=$len\n";
         print "previous: $exonsInA[$exonInd-1]\n";
         print "current: $exonsInA[$exonInd]\n";
       }
       $exonNow=$exonsSeqInA[$exonInd];
       $exonNow = $remainDNAs.$exonNow;  # add the previous remain one if there are some
       $exonLen = length($exonNow);
       $remainN = $exonLen % 3;
       $exonSeqTr = substr($exonNow, 0, $exonLen-$remainN);
    
       $exonPep = trans2pep($exonSeqTr);#translate
       $remainDNAs = substr($exonNow, $exonLen-$remainN, $remainN); #get the remain one if there are some
       $exonsInA[$exonInd] .=  "\t$prevRemainDNA\t$remainN\t$remainDNAs\t$exonPep";
       $prevRemainDNA=$remainDNAs;
       
       #print "$exonInd: $exonsInA[$exonInd]\n";
       #my $iii=<STDIN>;
       
     }
      # write the results into a file
     for( my $i=0; $i < $numExons; ++$i ) {
       #print "$i $exonsInA[$i]\n";
       print outDOC "$exonsInA[$i]\n";
       # my $iii=<STDIN>;
     }
 
 close DOC;
 close outDOC;
 
 #print "Finished!\n";
 
 sub reverseComplement {
   my $str=$_[0];
   $str =~ s/A/X/g;
   $str =~ s/T/A/g;
   $str =~ s/X/T/g;
   $str =~ s/G/X/g;
   $str =~ s/C/G/g;
   $str =~ s/X/C/g;
   
   reverse($str);
 
 }
 
 sub trans2pep {
   my $str=$_[0];
   my $len = length($str);
   my $pep = ' ';
   #print"$len, $str,\n";
   
   if($len >= 3) {
     $len /= 3;
     $pep = '';
     for(my $i=0; $i<$len; ++$i) {
      my $strNow = substr($str, $i*3, 3);
      if(!defined($geneticCodes{$strNow}) ) {
        #if(substr($strNow, 0, 1) eq 'N') {
          $pep .= 'X';
        #} else {
         # print "no defined genetic code:: $i, $str, ($strNow),\n";
        #}
        # my $iii=<STDIN>;
      } else {
        $pep .= $geneticCodes{$strNow};
      }
    }
   }
   
   $pep;
 }

 
