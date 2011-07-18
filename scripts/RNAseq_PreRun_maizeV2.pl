#!/usr/bin/perl

##      RNAseq_PreRun_maizeV2.pl
## 
##      Copyright 2011 Lin Wang <lw374@cornell.edu> -- Brutnell Lab  
## 
##      Version 1.5 -- 06/15/2011
## 
##      This program is free software; you can redistribute it and/or modify
##      it under the terms of the GNU General Public License as published by
##      the Free Software Foundation; either version 2 of the License, or
##      (at your option) any later version.
## 
##      his program is distributed in the hope that it will be useful,
##      but WITHOUT ANY WARRANTY; without even the implied warranty of
##      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##      GNU General Public License for more details.


						                              
##      This script uses BWA to index the rice genome and pre-processes all the    
##      annotation files needed for the RNA-seq script to run. This needs to be done   	
##      only once if genome and annotation remain the same for the analysis. 	
##      This script was written based on a script originally developed by Lalit		
##      Ponnala for the study of maize transcriptome (PMID:21037569) 							
  										



use warnings;
use strict;
use Getopt::Long;
use File::Basename;  			

## required BioPerl module
use Bio::Seq;    
use Bio::SeqIO;


#--scalar
my $genome_name;	    	## genome sequence file in FASTA format
my $genome_path;			## path to where genome sequence is stored
my $genome_file_type;		## file type of the genome   	
my $GFF3_annotaiton_input;	## genome annotation file in GFF3 format
my $GFF3_path;				## path to where annotation and parsed annotation files is stored
my $output_prefix;			## This will be used as prefix for all output. If this is not given, name of the genome will be used as prefix.
my $output_prefix1;			## Name of the genome, will be used as prefix if default name is not given.
my $output_path;			## Output path. If not given, annotation folder/path will be used
my $help;					## help message for this script
my $tic;					## track start of running time
my $toc;					## track end of running time
my $seq_read_length;		## lenghth of your RNA-seq reads (normally 40 bps, but if index is substracted, it will be 35)
my $chr;					## name of chromosomes in GFF annotation file
my $dis;					## distance used to calculate non-coding region coordinate


#--hash
my %E=(); 			## exon record 
my %T=(); 			## transcript safty check
my %G=(); 			## gene record information
my %TE=(); 			## transcript and corresponding exons
my %GT=();			## gene and corresponding transcripts
my %exonlen=();   	## exon length
my %exonstart=(); 	## exon start postion
my %exonstop=();  	## exon stop postion
my %tstrand=();		## strand information
my %tchr=();		## hash of hash store transcript on each chromosome

# Add support for explcit BWA path
my $bwa_app = '/usr/local3/bin/bwa';


######  adjustable parameter ######
my $minimal_junc_anchor_length = 10;	## minimal exon junction anchor length that is required to be mapped (This is something that need to be think about......)


######------- Get options  ---------------######

GetOptions
( 					
  'g:s'  => \$genome_name,           
  'a:s'  => \$GFF3_annotaiton_input, 
  'p:s'  => \$output_prefix,         
  'o:s'  => \$output_path,           
  'l:i'  => \$seq_read_length,
  'd:i'  => \$dis,         
  'help' => \$help,
);
		  
######------- parsing options ---------------######

if($help||!$genome_name||!$GFF3_annotaiton_input||!$seq_read_length||!$dis){ &usage; exit(0);}	## show help message if input are not in correct format, or help option is used 
($output_prefix1, $genome_path, $genome_file_type) = fileparse($genome_name, qr/\.[^.]*/);	## parse genome input file name
if (!$output_prefix) { $output_prefix = $output_prefix1; }									## if output_prefix is not give, use the genome name as the prefix instead
if (!$output_path) {$GFF3_path = dirname($GFF3_annotaiton_input).'/';}						## extract path of annotation file as output folder if specific output path is not given
else {$GFF3_path = $output_path;}															## if it is given, then use the given path/folder														
if (!(-d "$GFF3_path")) {system("mkdir $GFF3_path");}										## make folder if it's not already there
my $half_exon_junc_length = $seq_read_length - $minimal_junc_anchor_length;					## the calculated length of exon junctions to extract 
																							

######------- Indexing Genome with BWA ---------------######

print "\n---- Checking BWA index files... ----\n";				                        		## print message of starting BWA indexing
my @BWA_filetype = ('.amb', '.ann', '.bwt', '.pac', '.rbwt', '.rsa', '.sa', '.rpac');		        ## define BWA index file types with an array
my $BWA_CHECK = 0;										        									## initate a check 
foreach (@BWA_filetype) { if (-e $genome_path.$output_prefix1.$genome_file_type.$_) {++$BWA_CHECK;}}## count BWA files, add 1 to the checksum each time there is a file type match

#### If there appears to be BWA index files already, ask whether to skip the indexing process
if ($BWA_CHECK==8)
{
  print "It seems that you already have indexed the genome with BWA, do you want to index the genome again?\n";
  print "press \"Y/y\" to index the genome or press any other key to skip this step.\n"; 
  my $Answer = <STDIN>;
  chomp ($Answer);										
  if ($Answer eq "Y"|| $Answer eq "y")                                                         
    { 
      print "\n---- Indexing the genome with BWA... ----\n"; 
      system("$bwa_app index -a is $genome_name");
    } 		
  else {print "Skipping to the next step....\n"; }	                                        
}			  
#### if BWA index files are not complete in the genome folder, index the genome			
else { print "\n---- Indexing the genome with BWA... ----\n"; system("$bwa_app index -a is $genome_name");}							
				
			    
######---------- Start parsing the annotation ---------------######

print "\n\n---- Parsing the GFF annotation file... ----\n";		
$tic=time;

if ($GFF3_path !~ /\/$/)	{$GFF3_path = $GFF3_path.'/';}							## if the path variable does not end with '/', add it
my $exon_coordinate			= $GFF3_path.$output_prefix.'_exon.csv'; 				## exon coordinate file
my $transcript_coordinate	= $GFF3_path.$output_prefix.'_transcript.csv'; 			## transcript coordinate file
my $gene_coordinate			= $GFF3_path.$output_prefix.'_gene.csv'; 				## gene coordinate file
my $trans_exon_junction		= $GFF3_path.$output_prefix.'_transcript_exon_junc.csv';## transcript-exon mapping file
my $gene_trans				= $GFF3_path.$output_prefix.'_gene_transcript.csv';		## gene-transcript mapping file
my $exon_juc_cord			= $GFF3_path.$output_prefix.'_exon_junctions.csv';		## exon-junction coordinate file
my $exon_juc_seq			= $GFF3_path.$output_prefix.'_exon_junctions.fa';		## exon-junction sequence file
my $nonecoding_coordinate	= $GFF3_path.$output_prefix.'_non_coding.csv';			## non-coding file

			  
### creating files for outputs, otherwise quit ####
open(GENERECORD,">$gene_coordinate") or die("Cannot create $gene_coordinate"); 						#### $gene_coordinate = gene coordinate file
open(TRANSCRIPTRECORD,">$transcript_coordinate") or die("Cannot create $transcript_coordinate"); 	#### $transcript_coordinate = transcript coordinate file
open(EXONRECORD,">$exon_coordinate") or die("Cannot create $exon_coordinate"); 						#### $exon_coordinate = exon coordinate file
open(GTFILE,">$gene_trans") or die("Cannot open $gene_trans");										#### $gene_trans = gene and transcripts association file
		  
### open GFF annotation files, otherwise quit ####
open(GFFFILE,$GFF3_annotaiton_input) or die("Cannot open $GFF3_annotaiton_input");


######---------- Start looping through the gff annotation files and parsing the annotation ---------------######

while (<GFFFILE>) 
{					
  next if $_ =~ m/^\#/;										## skip comment lines
  chomp; 		  
  my @elems=split("\t",$_); 		        					## split each line by \t and store strings into @elems
  #if    	($elems[0] !~ /^chr/) {$chr='chr'.$elems[0]}; 	## The first elements contain chromosome info. If there is no "chr", add "chr"
  #if 	($elems[0] =~ /^chr/) {$chr=$elems[0]};			## if there is already "chr" in name, just store it in $chr without change
  $chr = $elems[0];
  next if 	($chr eq 'UNKNOWN'); 	   ## ignore the chromosomes not used for RPKM calculation





#######################################		
###		parsing "exon" records		###
#######################################
				
if ($elems[2] eq 'exon')
 {	  
	$elems[8]=~/Parent=(\S+);Name=(\S+)/; 	## matching parent transcript name and exon name and then stored in $1 and $2, any non-white space character
	my $tr=$1; 					## store match $1 as transcript name 
	my $exon=$2; 				
	
	if (exists $E{$exon}) { if ($E{$exon} ne $chr."\t".$elems[3]."\t".$elems[4]."\t".$elems[6]) { die("exon $exon found with different coordinates!\n"); } }
    $E{$exon}=$chr."\t".$elems[3]."\t".$elems[4]."\t".$elems[6]; 
	if (exists $TE{$tr}) { $TE{$tr}=$TE{$tr}.','.$exon; } else { $TE{$tr}=$exon; }
    print EXONRECORD $exon."\t".$chr."\t".$elems[3]."\t".$elems[4]."\t".$elems[6]."\n";
				  		  
			  
	#### calclulate exon length hash... if length is smaller than 0, report error and quit. This is a safty check
	$exonlen{$exon}=$elems[4]-$elems[3]+1; 
	if ($exonlen{$exon}<=0) { die("Error: exon length is smaller than zero for exon $exon\n"); }
				  
	#### store exon start position in %exonstart hash, key as the exon name
	$exonstart{$exon}=$elems[3]; 
				  
	#### store exon end position in %exonstop hash, key as the exon name
	$exonstop{$exon}=$elems[4];
 }
	


##############################
### parsing "mRNA" records ###
##############################
				
if ($elems[2] eq 'mRNA') 
 {
	$elems[8]=~/^ID=(\S+);Parent=(\S+);Name=(\S+);/; if ($1 ne $3) { die("ID not same as Name for transcript: $3\n"); }
      my $gene=$2; my $tr=$3; if (exists $T{$tr}) { die("transcript $tr already processed!\n"); }
      $T{$tr}=1; if (exists $GT{$gene}) { $GT{$gene}=$GT{$gene}.','.$tr; } else { $GT{$gene}=$tr; }
      print TRANSCRIPTRECORD $tr."\t".$chr."\t".$elems[3]."\t".$elems[4]."\t".$elems[6]."\n";
      $tstrand{$tr}=$elems[6]; $tchr{$chr}{$tr}=1;
    }
	
	
				
##############################
### parsing "gene" records ###
##############################

if ($elems[2] eq 'gene') 
 {
      $elems[8]=~/^ID=(\S+);Name=(\S+);/; if ($1 ne $2) { die("ID not same as Name for gene: $2\n"); }
      my $gene=$2; if (exists $G{$gene}) { die("gene $gene already processed!\n"); }
      $G{$gene}=1;
      print GENERECORD $gene."\t".$chr."\t".$elems[3]."\t".$elems[4]."\t".$elems[6]."\n";


 }
 
 }

## output gene-transcript mapping file 

foreach (keys %GT) { print GTFILE $_."\t".$GT{$_}."\n"; }


close(GFFFILE);
close(EXONRECORD);
close(TRANSCRIPTRECORD);
close(GENERECORD);
close(GTFILE);

## clean-up unwanted hashes
undef %E; 
undef %T; 
undef %G; 
undef %GT;


$toc=time;																			
print "\n\nTime used to process annotation is ".($toc-$tic)." seconds.\n\n";


######----------------------------- Generating junction files ------------------------------------######

#### prepare following output files 
## $trans_exon_junction = transcript-exon-junction mapping file
## $exon_juc_cord = exon-junction coordinate file
## $exon_juc_seq = exon-junction sequence file

open(EJSFILE,">$exon_juc_seq") or die ("Cannot create $exon_juc_seq");
open(EJCFILE,">$exon_juc_cord") or die ("Cannot create $exon_juc_cord");
open(TEJFILE,">$trans_exon_junction") or die ("Cannot create $trans_exon_junction");
			  
			  
my $jcount=0; 
my %TJ=(); 
my %J=();  ## splice junctions
my %L=();  ## length of chromosomes
			  
		  
print "\n---- Extracting junction sequences... ----\n\n";

$tic=time;

### processing genome files in fasta format
my $genome_obj = Bio::SeqIO->new(-file=>$genome_name,-format=>'fasta');

while (my $chromosome_obj=$genome_obj->next_seq) 
 {
	my $chr_ID=$chromosome_obj->id;  			### read in chromosome name
	#$chr_ID = 'chr'.$chr_ID;
	my $chr_seq=$chromosome_obj->seq();			### read in chromosome sequence.. this can be memory intensive	
	$L{$chr_ID}=$chromosome_obj->length();		### store length of chromosome in %L, key as chromosome name
								
	foreach (keys %{$tchr{$chr_ID}}) 	#### loop through each transcript on each chromosome, $_ as the $tr (transcript name)
	  {	
		my @exons=split(",",$TE{$_});	#### put exons of each transcript in an array @exons
		
		if (scalar(@exons)>1) 			#### if there is more than one exon for the transcript
				{		
					my @exonstarts=(); 	#### declare an array for exon start positions
					my @exonstops=();	#### declare an array for exon stop positions
										
					##### loop through all exons for each transcript and store start and stop positions into two array
					for (my $j=0; $j<scalar(@exons); ++$j) 
						{ push(@exonstarts,$exonstart{$exons[$j]});     #### save the exon start position in hash %exonstart
						  push(@exonstops,$exonstop{$exons[$j]}); 		#### save the exon stop position in hash %exonstop
						}
						
					#### if it is the plus strand, sort exon starts in increasing order of starts	
				    if ($tstrand{$_} eq '+')  
						{ 
						  my($sa,$I)=numsortascend(@exonstarts); 		## calling subroutine and return two references
						  @exonstarts=@$sa; 							## de-reference and get sorted exon start positions
						  @exonstops=putarrayinorder(\@exonstops,$I);	## put exon stops in the same order as the sorted exon starts
						  @exons=putarrayinorder(\@exons,$I);			## put exon names in the same sorted order
						}
						
					#### if it is the minus strand, sort exon stops in increasing order of starts								
					elsif ($tstrand{$_} eq '-') 
					    { 
						  my ($sa,$I)=numsortascend(@exonstops); 
						  @exonstops=@$sa; 
						  @exonstarts=putarrayinorder(\@exonstarts,$I); 
						  @exons=putarrayinorder(\@exons,$I);        
						}
				   
				   else { die("Cannot identify strand for transcript: $_\n"); }
					
					
					
                for (my $j=1; $j<scalar(@exons); $j++){ 
				  for (my $i=$j-1; $i<scalar(@exons); $i++){
				     if ($exonstarts[$j]>$exonstops[$j-1]){ 
						my $jkey=$chr_ID."\t".$tstrand{$_}."\t".$exonstops[$j-1]."\t".$exonstarts[$i];	  
				        if (not exists $J{$jkey}){ 
								$jcount++; 
								$J{$jkey}='junc_'.$jcount;						  
								### generating junction sequences
								my $juncseq=substr($chr_seq,$exonstops[$j-1]-$half_exon_junc_length,$half_exon_junc_length).substr($chr_seq,$exonstarts[$i]-1,$half_exon_junc_length);
								### check for empty junciton sequences
								if ($juncseq=~/^\s*$/) { die("empty junction sequence on transcript: $_\n"); }
								### check the length of junction sequences
								if (length($juncseq)!=(2*$half_exon_junc_length)) { die("junction sequence is not of correct length on transcript: $_\n"); }
								### output junction sequence files
								print EJSFILE '>'.$J{$jkey}."\n".$juncseq."\n";
							 }
						  ### output junction sequence cordination file
						  print EJCFILE $J{$jkey}."\t".$exons[$j-1]."\t".$exons[$i]."\t".$jkey."\n";
						  
						  if (exists $TJ{$_}) { $TJ{$_}=$TJ{$_}.','.$J{$jkey}; } 
						  else { $TJ{$_}=$J{$jkey}; }
						}
					}
					}
				  }
				  
				  ### calculate exon length for each transcripts/mRNA
				  my $tlen=0; 
				  foreach my $exon_name (@exons) { $tlen+=$exonlen{$exon_name};}
				  
				  ### output transcript-exon-junction mapping file
                  if (exists $TJ{$_}) {print TEJFILE $_."\t".$tlen."\t".join(",",@exons)."\t".$TJ{$_}."\n";}
				  else {print TEJFILE $_."\t".$tlen."\t".join(",",@exons)."\n"};
				
				}
 }

close(TEJFILE);
close(EJSFILE);
close(EJCFILE);
print "number of coordinate-level-redundancy-removed junctions whose sequences have been written out = ".scalar(keys %J)."\n";

$toc=time;																			
print "\nTime used to extract junction sequences is ".($toc-$tic)." seconds.\n\n";
		  


#### clean up
undef %TE; undef %TJ; undef %J;
undef %exonlen; undef %exonstart; undef %exonstop=(); undef %tstrand; undef %tchr;


######----------------------------- Indexing junction files ------------------------------------######

print "\n---- Indexing the junction file... ----\n\n";			
		
system("$bwa_app index -a is $exon_juc_seq");	




######-------------------------- Processing non-coding regions ---------------------------------######

print "\n---- Extracting noncoding region coordinates.... ----\n";
$tic=time;
			  
my %gstarts=(); 
my %gstops=();
my %GG=();			  
### open gene corrdination file

open(GENEFILE,$gene_coordinate) or die("Cannot open $gene_coordinate");
			  
while (<GENEFILE>) 
			  {
				my $gene_line=$_; 
				chomp($gene_line); 
				next if ($gene_line=~/^\s*$/);
				my @gelems=split("\t",$gene_line);
				if (exists $GG{$gelems[0]}) {next;}
				push(@{$gstarts{$gelems[1]}},$gelems[2]); 
				push(@{$gstops{$gelems[1]}},$gelems[3]);
				$GG{$gelems[0]}=1;
			  }
			  close(GENEFILE);
			  
open(OFILE,">$nonecoding_coordinate") or die("Cannot create $nonecoding_coordinate\n");
			  
foreach my $chr (keys %gstarts) 
	{
		my @starts=@{$gstarts{$chr}}; 
		my @stops=@{$gstops{$chr}};
		my ($sa,$I)=numsortascend(@starts); 
		my @sortedstarts=@$sa; 
		my @sortedstops=putarrayinorder(\@stops,$I);
		
		if (($sortedstarts[0]-$dis)>$half_exon_junc_length) { print OFILE $chr."\t"."1"."\t".($sortedstarts[0]-$dis)."\n"; }
		
		my $maxsortedstop=$sortedstops[0];
				
		for (my $j=1; $j<scalar(@sortedstarts); $j++) 
				{
				  $maxsortedstop=max($maxsortedstop,$sortedstops[$j-1]);
				  if (($sortedstarts[$j]-$maxsortedstop)>(2*$dis+$half_exon_junc_length)) { print OFILE $chr."\t".($maxsortedstop+$dis)."\t".($sortedstarts[$j]-$dis)."\n"; }
				 }
			if (($L{$chr}-$maxsortedstop)>($dis+$half_exon_junc_length)) {print OFILE $chr."\t".($maxsortedstop+$dis)."\t".$L{$chr}."\n"; }
     }
close(OFILE);

$toc=time;																			
print "\nTime used to extract noncoding region coordinates is ".($toc-$tic)." seconds.\n\n";	
print "---- ALL DONE!!! ----\n\n";

exit();



###=====================================================================================================================###
###============================================  Sub-rountines from here on ============================================###
###=====================================================================================================================###


## ------------------------------------------------------------------------------------
## Give out usage information
## ------------------------------------------------------------------------------------
sub usage 
{
    print'    
    Program: RNAseq_PreRun.pl
    Version: 1.5 (06/15/2011)
    Contact: Lin Wang <lw374@cornell.edu>,  BIT Cornell, Brutnell Lab.

    This script requires BWA (alignment via Burrows-Wheeler transformation) to run, make sure it has been installed.
     
    The options are:
    
    -g 	  : input genome sequence file, must be in FASTA format, path can be included
    -a 	  : GFF3 annotation file for the genome sequence.     
    -o	  : path/folder to save parsed output files.. if not given, your annotation folder will be used as output folder
    -p 	  : prefix that will be used for all your output files... If not given, the genome name will be used as prefix.
    -l	  : length of your illumina sequencing reads (e.g. 35), important for extracting junction sequences
    -d	  : distance used to calculate None-coding region coordinate, e.g. 1kb = 1000, 5kb = 5000
    -help : print this help message'."\n\n";
}



## ------------------------------------------------------------------------------------
## MAX: Find maximum of 2 input elements
## USAGE: $maxelem=max($elem1,$elem2);
## INPUT: $elem1,$elem2 = two numbers
## ------------------------------------------------------------------------------------
sub max {
  if (scalar(@_)!=2) { die("There must be 2 numbers passed to this function!"); }
  my ($Element_one,$Element_two)=@_;
  if ($Element_one>=$Element_two) { return $Element_one; }
  else { return $Element_two; }
}
## ------------------------------------------------------------------------------------


## -------------------------------------------------------
## Numeric sort array in ascending order and return order
## Array can contain the same number more than once
## USAGE:
## ($sa,$I)=numsortascend(@arr);
## @arr = array of numbers (int or float)
## $sa = reference to sorted array (retrieve as @$sa)
## $I = reference to sorted indices (retrieve as @$I)
## -------------------------------------------------------
sub numsortascend {
  my @array=@_;  ## input is an array
  my %ind=();    ## set up a hash
  
  for (my $j=0; $j<scalar(@array); $j++)   ## loop through the array 
  {
    #if ($array[$j] !~ /^[+-]?\d+\.?\d*$/ ) { die("this element of the input array does not look like a number: $array[$j]\nQuitting...\n"); }
    push(@{$ind{$array[$j]}},$j); 
  }

  my @sortedarray = sort {$a <=> $b} @array;
  my @sortedind=();
  
  for (my $j=0; $j<scalar(@sortedarray); $j++) 
  {
    push(@sortedind,shift(@{$ind{$sortedarray[$j]}} ));
  }
  
  return (\@sortedarray,\@sortedind);
}
## -------------------------------------------------------


## -------------------------------------------------------
## Numeric sort array in descending order and return order
## Array can contain the same number more than once
## USAGE:
## ($sa,$I)=numsortdescend(@arr);
## @arr = array of numbers (int or float)
## $sa = reference to sorted array (retrieve as @$sa)
## $I = reference to sorted indices (retrieve as @$I)
## -------------------------------------------------------
sub numsortdescend {
  my @array=@_;

  my %ind=();
  for (my $j=0; $j<scalar(@array); $j++) {
    #if ($array[$j] !~ /^[+-]?\d+\.?\d*$/ ) { die("this element of the input array does not look like a number: $array[$j]\nQuitting...\n"); }
    push(@{$ind{$array[$j]}},$j);
  }

  my @sortedarray = sort {$b <=> $a} @array;
  my @sortedind=();
  for (my $j=0; $j<scalar(@sortedarray); $j++) {
    push(@sortedind,shift(@{$ind{$sortedarray[$j]}} ));
  }
  return (\@sortedarray,\@sortedind);
}
# -------------------------------------------------------


## ------------------------------------------------
## Re-arrange elements of array in specified order
## NOTE:
## array can contain either numbers or strings
## USAGE:
## @newarray = putarrayinorder($ref_arr,$ref_I)
## $ref_arr = reference to the original array
## $ref_I = reference to the array of indices
## @newarray = re-arranged array
## ------------------------------------------------
sub putarrayinorder {
  my ($ref_arr,$ref_I)=@_;
  my @array=@$ref_arr; my @I=@$ref_I;
  if (scalar(@array)!=scalar(@I)) { die("arrays of elements and indices have different sizes!"); }
  my @newarray=();
  for(my $j=0; $j<scalar(@I); $j++) { push(@newarray,$array[$I[$j]]); }
  return @newarray;
}
# ------------------------------------------------


