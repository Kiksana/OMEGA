#!/usr/bin/env perl
my $command="$0 ";
my $i=0;
while (defined($ARGV[$i])){
  $command .= "$ARGV[$i] ";
  $i++;
}
use Getopt::Std;
use Cwd;
use strict;
# Default parameters
*LOG=*STDERR;
my $fh;
my $verbose=0;
my $minimum_readCount= 10;
my $clade = 'strain';
my $minimum_relative_abundance = "0.01";
my $normfactor = 1000000;
my $normalization = 1;
my $observed_min_rel_abundance=9999;
my $observed_max_rel_abundance=0;
my $lines=0;
my $hits=0;
my $aveReadLen;
my $total_reads=0;
my $rejected_reads=0;
my $uniq_ratio;
my $nm_over_nuc_ratio;
my $depth;
my $expected_coverage;
my $new_abundance;
my $new_abundance_subtractor = 0.0000000001;
# minimum ratio of uniqReadCount/ReadCount
my $ratio_cutoff="0.005";

# maximum ratio of nm/nucleotides i.e. fraction of errors or mismatches compared to number of basepairs
my $max_nm_nuc_ratio="0.01";

my $min_coverage=0;
my $single_end_reads=0;
my $Nreads;
#
# A flag to determine if failed hits are printed to file defined by option -f
#
my $failed_flag=0;
#
# Process command line
#
getopts('hi:o:vl:n:a:c:Hr:sm:g:f:F')||Usage();
#
# Usage
#
if (defined($Getopt::Std::opt_h)||defined($Getopt::Std::opt_h)){
  # Print help message
  Usage();
}

sub Usage {
  print ("Usage: $0 [-h] [-i name] [-o name] [-l name] [-n number] [-c name] [-a number] [-v] [-r number] [-m number] [-g number] [-H] [-s] [-f name] [-F number] \n");
  print ("Description:\n");
  print ("$0 - OMEGA_module2.pl - clade classification of the OMEGA_module3.pl annotation files. Filtering by criteria selected options\n");
  print ("\n");
  print ("Options:\n");
  print ("  -h  : display this message\n");
  print ("  -i  : input file name [STDIN]\n");
  print ("  -o  : output file name [STDOUT]\n");
  print ("  -s  : read count are treated as single-end reads [off]\n");
  print ("  -n  : minimum number of mapped reads [$minimum_readCount]\n");
  print ("  -a  : minimum abundance % abundance [$minimum_relative_abundance]\n");
  print ("  -r  : minimum uniqueReadCount/ReadCount ratio [$ratio_cutoff]\n");
  print ("  -g  : minimum coverage [$min_coverage]\n");
  print ("  -m  : maximum nucleotide error fraction [$max_nm_nuc_ratio]\n");
  print ("  -c  : collapse at clade level (strain, species, genus, family, order, class, phylum, superfamily) - default [$clade]\n");
  print ("  -f  : output File with failed hits\n");
  print ("  -F  : size relative abundance normalization by factor of 1000000\n");
  print ("  -H  : print Header [off]\n");
  print ("  -l  : logfile [STDERR]\n");
  print ("  -v  : Verbose [off]\n\n");
  print ("       Examples:\n");
  print ("       # At strain -r 0.005 should be used\n");
  print ("       # At species -r 0.005 should be used and option -v prints lowest accepted relative abundance which is used in a second round with -r 0\n");
  print ("       # At all other levels i.e. genus phylum etc -r 0 should be used i.e. uniq read count ratio is not used\n");
  print ("       # option -n 10 sets a lower readCount cutoff. Two strains belonging to the same species with each say 5 reads mapped to each,\n");
  print ("       # results in none of them being selected at strain level, but at species level there are 10 reads and the species will be accepted as a true hit.\n");
  print ("       coverage= 'covered_positions'/'size'\n");
  print ("       paired-end reads: abundance= (100*read_count)/(size*2)\n");
  print ("       single-end reads: abundance= (100*read_count)/(size), specify -s for single-end reads\n");
  print ("       uniqReadCount ratio (option r) = uniqReadCount/readCount\n");
  print ("       misMatch ratio (option m) = misMatches/nucleotides\n");
  print ("       normalization = ('normfactor'/'tot_reads')*('relative_abundance')\n");
  print ("\n");
 exit;
} # Usage

#
# Open input
#
if (not defined($Getopt::Std::opt_i)){
  # Read from standard input
  *INP = *STDIN;
} 
else{
  # Read from file
  if (($Getopt::Std::opt_i=~/\.gz$/) || ($Getopt::Std::opt_i=~/\.Z$/)){
    open(INP,"gunzip -c $Getopt::Std::opt_i |") || die ("can't open file $Getopt::Std::opt_i: $!");
  }
  else{
    open(INP,"<$Getopt::Std::opt_i") || die ("can't open file $Getopt::Std::opt_i: $!");
  }
}
#
# If not file name is given, use standard output
#
if (not defined($Getopt::Std::opt_o)){
  # Output goes to std output
  *OUT = *STDOUT;
} else {
  # Open file to write to
  open(OUT, ">$Getopt::Std::opt_o") || die ("can't open file $Getopt::Std::opt_o: $!");
}
if (defined($Getopt::Std::opt_l)){
    open(LOG,">$Getopt::Std::opt_l");
}
if (defined($Getopt::Std::opt_v)){
    $verbose=1;
}
if (defined($Getopt::Std::opt_n)){
    $minimum_readCount=$Getopt::Std::opt_n;
}
if (defined($Getopt::Std::opt_c)){
    $clade=$Getopt::Std::opt_c;
}
if (defined($Getopt::Std::opt_a)){
    $minimum_relative_abundance=$Getopt::Std::opt_a;
}
if (defined($Getopt::Std::opt_r)){
    $ratio_cutoff=$Getopt::Std::opt_r;
}
if (defined($Getopt::Std::opt_g)){
    $min_coverage=$Getopt::Std::opt_g;
}
if (defined($Getopt::Std::opt_m)){
    $max_nm_nuc_ratio=$Getopt::Std::opt_m;
}
if (defined($Getopt::Std::opt_s)){
    $single_end_reads=1;
}

if (defined($Getopt::Std::opt_f)){
    open(FAILED,">$Getopt::Std::opt_f");
    $failed_flag=1;
}
###############################################################################
# Main
#
###############################################################################
my $datestring = localtime();
my $thisDir=cwd();
if ($verbose){
    print LOG "## Local date and time $datestring - Start program\n";
    print LOG "# $command\n";
    print LOG "# working dir: $thisDir\n\n";
}
my $tot_readCount;
my $sample_acc;
my @w=();
my $db;
my $strain;
my $size;
my $nucleotides;
my $covered_positions;
my $readCount;
my $readCountUniq;
my $nm;

my $superfamily;
my $superfamilyTax;
my $phylum;
my $phylumTax;
my $class;
my $classTax;
my $order;
my $orderTax;
my $family;
my $familyTax;
my $genus;
my $genusTax;
my $species;
my $speciesTax;
my $strain;
my $strainTax;

my $coverage;
my $depth;
my $relative_abundance;
my $desc;


my %collapse_db=();
# store database names
my %d=();
my @databases;

my $taxid;
my $col;
my $lastCol;
my $taxid_column;
my @columns= qw(superfamily superfamilyTax phylum phylumTax class classTax order orderTax family familyTax genus genusTax species speciesTax strain strainTax);
if ($clade eq 'superfamily'){
    $col=11;
    $taxid_column=12;
    $lastCol=14;
}
elsif ($clade eq 'phylum'){
    $col=13;
    $taxid_column=14;
    $lastCol=16;
}
elsif ($clade eq 'class'){
    $col=15;
    $taxid_column=16;
    $lastCol=18;
}
elsif ($clade eq 'order'){
    $col=17;
    $taxid_column=18;
    $lastCol=20;
}
elsif ($clade eq 'family'){
    $col=19;
    $taxid_column=20;
    $lastCol=22;
}
elsif ($clade eq 'genus'){
    $col=21;
    $taxid_column=22;
    $lastCol=24;
}
elsif ($clade eq 'species'){
    $col=23;
    $taxid_column=24;
    $lastCol=26;
}
elsif ($clade eq 'strain'){
    $col=3;
    $taxid_column=26;
    $lastCol=28;
}
else{
    print LOG "unknown value for option -c\n";
    exit;
}
my $name;
if (defined($Getopt::Std::opt_H)){
    print OUT "# study_sample\ttot_readCount\tName\tS_Abundance (%)\tR_Abundance (%)\tSize (bp)\tSeq_count\tNucleotides\tCovered_positions\tCoverage\tDepth\tReads\tReads_uniq\tEdit_distance\tDescription";
    if ($failed_flag){
	print FAILED "# study_sample\ttot_readCount\tName\tS_Abundance (%)\tR_Abundance (%)\tSize (bp)\tSeq_count\tNucleotides\tCovered_positions\tCoverage\tDepth\tReads\tReads_uniq\tEdit_distance\tDescription";
    }

    my $end = $lastCol - 13;
    for (my $i=0;$i<=$end;$i++){
	print OUT "\t$columns[$i]";
	if ($failed_flag){
	    print FAILED "\t$columns[$i]";
	}
    }
    print OUT "\n";
    if ($failed_flag){
	print FAILED "\n";	
    }

}

while (defined($_=<INP>)){
    chomp;
    @w=split(/\t/);
    $tot_readCount = $w[1];
    

    if ($w[0] ne $sample_acc and $sample_acc ne '') {
       
	if (defined($Getopt::Std::opt_F)){
	    $normalization = $normfactor/$tot_readCount;
	}
    }
    if ($clade != 'species' or $clade != 'strain') {
	&collapse($sample_acc);
    } else {

	$size=$w[4];
	$readCount = $w[7];
	$collapse_db{$taxid}{readCount} += $readCount;
	$collapse_db{$taxid}{size} += $size;
	if ($single_end_reads){
	    $collapse_db{$taxid}{relative_abundance} += ($readCount*100)/$size;
	}
	else{
	    $collapse_db{$taxid}{relative_abundance} += ($readCount*100)/($size*2);
	}
	foreach my $key (sort {$collapse_db{$b}{relative_abundance} <=> $collapse_db{$a}{relative_abundance}} keys %collapse_db){
	
	    if ($collapse_db{$taxid}{relative_abundance} > $observed_max_rel_abundance){
		$observed_max_rel_abundance = $collapse_db{$key}{relative_abundance};
	    }
	    if ($collapse_db{$taxid}{relative_abundance} < $observed_min_rel_abundance){
		$observed_min_rel_abundance = $collapse_db{$key}{relative_abundance};
	    }
	}
    }
	
    %collapse_db = ();

    if ($clade != 'species' or $clade != 'strain') {
	$Nreads = $w[1];    
	$sample_acc = $w[0];
	$db=$w[2];
	if (! exists($d{$db})){
	    print LOG "# databases: $db\n" if ($verbose);
	    push(@databases,$db);
	    $d{$db}=$db;
	}
	$strain=$w[3];
	$size=$w[4];
	$nucleotides=$w[5];
	$covered_positions=$w[6];
	$readCount = $w[7];
	$total_reads += $readCount;
	$readCountUniq = $w[8];
	$nm=$w[9];
	$desc=$w[10];
	for (my $i=11;$i<=26;$i++){
	    if ($w[$i] eq 'Unknown'){
		$w[$i] = '-';
	    }
	}
    
	#
	# Taxonomy
	#
    
	$superfamily= $w[11];
	$superfamilyTax= $w[12];

	$phylum=$w[13];
	$phylumTax=$w[14];

	$class=$w[15];
	$classTax=$w[16];

	$order=$w[17];
	$orderTax=$w[18];

	$family=$w[19];
	$familyTax=$w[20];
    
	$genus=$w[21];
	$genusTax=$w[22];

	$species=$w[23];
	$speciesTax=$w[24];

	$strain=$w[25];
	$strainTax=$w[26];

	$coverage=$covered_positions/$size;
	$depth=$nucleotides/$size;

    
	$name=$w[$col];
	$taxid=$w[$taxid_column];

	if ($clade eq 'strain'){
	    $taxid=$name;
	    $collapse_db{$taxid}{desc}=$desc;
	}
	if ($taxid eq '-'){
	    $name="Unclassified";
	}

	if ($clade ne  'strain'){
	    $collapse_db{$taxid}{desc}="collapsed at $clade level";
	}

	$collapse_db{$taxid}{seqCount}++;
	$collapse_db{$taxid}{taxid} = $taxid;
	$collapse_db{$taxid}{name} = $name;
        $collapse_db{$taxid}{total_readCount} = $tot_readCount;
	$collapse_db{$taxid}{size} += $size;
	$collapse_db{$taxid}{nucleotides} += $nucleotides;
	$collapse_db{$taxid}{covered_positions} += $covered_positions;
	$collapse_db{$taxid}{readCount} += $readCount;
	$collapse_db{$taxid}{readCountUniq} += $readCountUniq;
	$collapse_db{$taxid}{nm} += $nm;

	if ($single_end_reads){
	    $collapse_db{$taxid}{relative_abundance} += ($readCount*100)/$size;
	}
	else{
	    $collapse_db{$taxid}{relative_abundance} += ($readCount*100)/($size*2);
	}
	
	$collapse_db{$taxid}{superfamily} = $superfamily;
	$collapse_db{$taxid}{superfamilyTax} = $superfamilyTax;
	
	$collapse_db{$taxid}{phylum} = $phylum;
	$collapse_db{$taxid}{phylumTax} = $phylumTax;
	
	$collapse_db{$taxid}{class} = $class;
	$collapse_db{$taxid}{classTax} = $classTax;
	
	$collapse_db{$taxid}{order} = $order;
	$collapse_db{$taxid}{orderTax} = $orderTax;
	
	$collapse_db{$taxid}{family} = $family;
	$collapse_db{$taxid}{familyTax} = $familyTax;
	
	$collapse_db{$taxid}{genus} = $genus;
	$collapse_db{$taxid}{genusTax} = $genusTax;
	
	$collapse_db{$taxid}{species} = $species;
	$collapse_db{$taxid}{speciesTax} = $speciesTax;
	
	$collapse_db{$taxid}{strain} = $strain;
	$collapse_db{$taxid}{strainTax} = $strainTax;
    }
}

if (defined($Getopt::Std::opt_i)){
    close(INP);
}

if ($clade == 'species' or $clade == 'strain') {
    $new_abundance = $observed_min_rel_abundance - $new_abundance_subtractor;
    $minimum_relative_abundance = $new_abundance;
    $ratio_cutoff = 0;
    @w = ();
    %collapse_db = ();

    open(INP,"<$Getopt::Std::opt_i") || die ("can't open file $Getopt::Std::opt_i: $!");
    
    while (defined($_=<INP>)){
	chomp;
	@w=split(/\t/);
	$tot_readCount = $w[1];
    

	if ($w[0] ne $sample_acc and $sample_acc ne '') {
       
	    if (defined($Getopt::Std::opt_F)){
		$normalization = $normfactor/$tot_readCount;
	    }

	    &collapse($sample_acc);
	    %collapse_db = ();
	}

	$Nreads = $w[1];    
	$sample_acc = $w[0];
	$db=$w[2];
	if (! exists($d{$db})){
	    print LOG "# databases: $db\n" if ($verbose);
	    push(@databases,$db);
	    $d{$db}=$db;
	}
	$strain=$w[3];
	$size=$w[4];
	$nucleotides=$w[5];
	$covered_positions=$w[6];
	$readCount = $w[7];
	$total_reads += $readCount;
	$readCountUniq = $w[8];
	$nm=$w[9];
	$desc=$w[10];
	for (my $i=11;$i<=26;$i++){
	    if ($w[$i] eq 'Unknown'){
		$w[$i] = '-';
	    }
	}
    
	#
	# Taxonomy
	#
    
	$superfamily= $w[11];
	$superfamilyTax= $w[12];

	$phylum=$w[13];
	$phylumTax=$w[14];

	$class=$w[15];
	$classTax=$w[16];

	$order=$w[17];
	$orderTax=$w[18];

	$family=$w[19];
	$familyTax=$w[20];
    
	$genus=$w[21];
	$genusTax=$w[22];

	$species=$w[23];
	$speciesTax=$w[24];

	$strain=$w[25];
	$strainTax=$w[26];

	$coverage=$covered_positions/$size;
	$depth=$nucleotides/$size;

    
	$name=$w[$col];
	$taxid=$w[$taxid_column];

	if ($clade eq 'strain'){
	    $taxid=$name;
	    $collapse_db{$taxid}{desc}=$desc;
	}
	if ($taxid eq '-'){
	    $name="Unclassified";
	}

	if ($clade ne 'strain'){
	    $collapse_db{$taxid}{desc}="collapsed at $clade level";
	}

	$collapse_db{$taxid}{seqCount}++;
	$collapse_db{$taxid}{taxid} = $taxid;
	$collapse_db{$taxid}{name} = $name;
        $collapse_db{$taxid}{total_readCount} = $tot_readCount;
	$collapse_db{$taxid}{size} += $size;
	$collapse_db{$taxid}{nucleotides} += $nucleotides;
	$collapse_db{$taxid}{covered_positions} += $covered_positions;
	$collapse_db{$taxid}{readCount} += $readCount;
	$collapse_db{$taxid}{readCountUniq} += $readCountUniq;
	$collapse_db{$taxid}{nm} += $nm;

	if ($single_end_reads){
	    $collapse_db{$taxid}{relative_abundance} += ($readCount*100)/$size;
	}
	else{
	    $collapse_db{$taxid}{relative_abundance} += ($readCount*100)/($size*2);
	}
	
	$collapse_db{$taxid}{superfamily} = $superfamily;
	$collapse_db{$taxid}{superfamilyTax} = $superfamilyTax;
	
	$collapse_db{$taxid}{phylum} = $phylum;
	$collapse_db{$taxid}{phylumTax} = $phylumTax;
	
	$collapse_db{$taxid}{class} = $class;
	$collapse_db{$taxid}{classTax} = $classTax;
	
	$collapse_db{$taxid}{order} = $order;
	$collapse_db{$taxid}{orderTax} = $orderTax;
	
	$collapse_db{$taxid}{family} = $family;
	$collapse_db{$taxid}{familyTax} = $familyTax;
	
	$collapse_db{$taxid}{genus} = $genus;
	$collapse_db{$taxid}{genusTax} = $genusTax;
	
	$collapse_db{$taxid}{species} = $species;
	$collapse_db{$taxid}{speciesTax} = $speciesTax;
	
	$collapse_db{$taxid}{strain} = $strain;
	$collapse_db{$taxid}{strainTax} = $strainTax;
    }
}

if ($verbose){
    my $perc=$rejected_reads*100/$total_reads;
    print LOG "# Number of hits processed:\t$lines\n";
    print LOG "# Number of hits accepted:\t$hits\n";
    printf LOG "# Total_reads= $total_reads\trejected_reads= $rejected_reads\t % rejected= %.3f\n",$perc;
    print LOG ("# Maximum accepted abundance value=\t$observed_max_rel_abundance\n");
    print LOG ("# Minimum accepted abundance value=\t$new_abundance\n");
}

if (defined($Getopt::Std::opt_i)){
    close(INP);
}
if (defined($Getopt::Std::opt_o)){
    close(OUT);
}
if (defined($Getopt::Std::opt_l)){
    close(LOG);
}

sub collapse {
    foreach my $key (sort {$collapse_db{$b}{relative_abundance} <=> $collapse_db{$a}{relative_abundance}} keys %collapse_db){
	if (! exists($collapse_db{$key}{size})){
	    next;
	}
	$lines++;

	if ($collapse_db{$key}{relative_abundance} < $minimum_relative_abundance){
	    $rejected_reads += $collapse_db{$key}{readCount};
	    if ($failed_flag){
		$fh=*FAILED;
		&printout_collapse($fh,$key);
	    }
	    next;
	}

	#
	# store some min and max observed relative abundances to be written to logfile only
	#
	if ($collapse_db{$key}{readCount} < $minimum_readCount){
	    $rejected_reads += $collapse_db{$key}{readCount};
	    if ($failed_flag){
		$fh=*FAILED;
		&printout_collapse($fh,$key);
	    }
	    next;
	}

	$uniq_ratio=$collapse_db{$key}{readCountUniq}/$collapse_db{$key}{readCount};
	if ($uniq_ratio < $ratio_cutoff){
	    $rejected_reads += $collapse_db{$key}{readCount};
	    if ($failed_flag){
		$fh=*FAILED;
		&printout_collapse($fh,$key);
	    }
	    next;
	}
	
	$nm_over_nuc_ratio=$collapse_db{$key}{nm}/$collapse_db{$key}{nucleotides};
	if ($nm_over_nuc_ratio > $max_nm_nuc_ratio){
	    $rejected_reads += $collapse_db{$key}{readCount};
	    if ($failed_flag){
		$fh=*FAILED;
		&printout_collapse($fh,$key);
	    }
	    next;
	}

    
	$depth=$collapse_db{$key}{nucleotides}/$collapse_db{$key}{size};
	$expected_coverage=1-exp(-$depth);
	my $coverage = $collapse_db{$key}{covered_positions}/$collapse_db{$key}{size};
	
	if ($coverage < $min_coverage){
	    $rejected_reads += $collapse_db{$key}{readCount};
	    if ($failed_flag){
		$fh=*FAILED;
		&printout_collapse($fh,$key);
	    }
	    next;
	}
	$fh=*OUT;
	if ($clade != 'species' or $clade != 'strain') {
	    if ($collapse_db{$key}{relative_abundance} > $observed_max_rel_abundance){
		$observed_max_rel_abundance = $collapse_db{$key}{relative_abundance};
	    }
	    if ($collapse_db{$key}{relative_abundance} < $observed_min_rel_abundance){
		$observed_min_rel_abundance = $collapse_db{$key}{relative_abundance};
	    }
	}
	$hits++;
	&printout_collapse($fh,$key);
    }
}

sub printout_collapse{
    my ($fh, $key)=@_;
    
    my $normalized_relative_abundance = $collapse_db{$key}{relative_abundance} * $normalization;
    printf $fh "%s", $sample_acc;
    printf $fh "\t%d",$collapse_db{$key}{total_readCount};
    printf $fh "\t%s",$collapse_db{$key}{name};
    printf $fh "\t%.3f",$normalized_relative_abundance;
    my $value = $collapse_db{$key}{readCount}*100/$Nreads;
    printf $fh "\t%.3f",$value;

    printf $fh "\t%d",$collapse_db{$key}{size};
    printf $fh "\t%d",$collapse_db{$key}{seqCount};
    printf $fh "\t%d",$collapse_db{$key}{nucleotides};
    
    
    printf $fh "\t%d",$collapse_db{$key}{covered_positions};

    $coverage=$collapse_db{$key}{covered_positions}/$collapse_db{$key}{size};
    printf $fh "\t%.3f",$coverage;
    
    my $depth=$collapse_db{$key}{nucleotides}/$collapse_db{$key}{size};
    printf $fh "\t%.3f",$depth;
    
    printf $fh "\t%d",$collapse_db{$key}{readCount};
    printf $fh "\t%d",$collapse_db{$key}{readCountUniq};
    printf $fh "\t%d",$collapse_db{$key}{nm};
    printf $fh "\t%s",$collapse_db{$key}{desc};
    
    my $column=14;
    if ($column <= $lastCol){
	printf $fh "\t%s", $collapse_db{$key}{superfamily};
	printf $fh "\t%s", $collapse_db{$key}{superfamilyTax};
	$column +=2;
    }
    if ($column <= $lastCol){
	printf $fh "\t%s", $collapse_db{$key}{phylum};
	printf $fh "\t%s", $collapse_db{$key}{phylumTax};
	$column +=2;
    }
    if ($column <= $lastCol){
	printf $fh "\t%s", $collapse_db{$key}{class};
	printf $fh "\t%s", $collapse_db{$key}{classTax};
	$column +=2;
    }
    if ($column <= $lastCol){
	printf $fh "\t%s", $collapse_db{$key}{order};
	printf $fh "\t%s", $collapse_db{$key}{orderTax};
	$column +=2;
    }
    if ($column <= $lastCol){
	printf $fh "\t%s", $collapse_db{$key}{family};
	printf $fh "\t%s", $collapse_db{$key}{familyTax};
	$column +=2;
    }
    if ($column <= $lastCol){
	printf $fh "\t%s", $collapse_db{$key}{genus};
	printf $fh "\t%s", $collapse_db{$key}{genusTax};
	$column +=2;
    }
    if ($column <= $lastCol){
	printf $fh "\t%s", $collapse_db{$key}{species};
	printf $fh "\t%s", $collapse_db{$key}{speciesTax};
	$column +=2;
    }
    if ($column <= $lastCol){
	printf $fh "\t%s", $collapse_db{$key}{strain};
	printf $fh "\t%s", $collapse_db{$key}{strainTax};
    }
    print $fh "\n"
}
