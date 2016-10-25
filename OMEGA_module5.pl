#!/usr/bin/env perl
my $command="$0 ";
my $i=0;
while (defined($ARGV[$i])){
    $command .= "$ARGV[$i] ";
    $i++;
}

use lib '/cm/local/apps/environment-modules/3.2.10/init/';
use Getopt::Std;
use Cwd;
use perl;
use strict;

module("load tools R/3.2.5");

# Default parameters
*LOG=*STDERR;
my $verbose = 0;
my $reRun = 0;
my $out_R = 0;
my $out_r = 0;
my %set;
my %set1;
my %set2;
my %stored;
my %richness;
my @stored_study_sample1;
my @stored_study_sample2;
my $num_samples_set1 = 0;
my $num_samples_set2 = 0;
my %total_reads;
my @total_reads;
my @nn;
my $percent_present_in_samples = "0.0";
my $condition_matrix = "condition.txt";
my $R_dir = "OMEGA_R";
my $R_file = "abundance_mat.txt";
my $r_file = "counts_mat.txt";
my $divers_file_A = "alpha_diversity_A.txt";
my $divers_file_B = "alpha_diversity_B.txt";
my $batch_R = "batch.R";
my $prog_R = "R";
#
# Process command line
#
getopts('hi:vl:A:B:RHf:d:rOob')||Usage();
#
# Usage
#
if (defined($Getopt::Std::opt_h)||defined($Getopt::Std::opt_h)){
  # Print help message
    Usage();
}

sub Usage {
    print ("Usage: $0 [-h] [-i name] [-A name] [-B name] [-R] [-r] ... [-v]\n");
    print ("Description:\n");
    print ("$0 - pipeline for autOmated MEtaGenomic Analysis (OMEGA) - module 5: Construction of tab-separated matrices and comparative analysis between two cohorts via \"DESeq2\" R package.\n");
    print ("\n");
    print ("Options:\n");
    print ("  -h  : display this message\n");
    print ("  -i  : input file name\n");
    print ("  -A  : list of study_sample accessions separated by new line\n");
    print ("  -B  : list of study_sample accessions separated by new line\n");
    print ("  -H  : print header line [OFF]\n");
    print ("  -f  : clade minimum occurance in either set of samples % occurance [$percent_present_in_samples]\n");
    print ("  -R  : R matrix-ready tab separated file with normalized abundances [OFF]\n");
    print ("  -r  : R matrix-ready tab separated file with read counts [OFF]\n");
    print ("  -O  : print to STDOUT tab-separated normalized abundances [OFF]\n");
    print ("  -o  : print to STDOUT tab-separated read counts [OFF]\n");
    print ("  -b  : re-run the module with customized batch.R script found in directory 'OMEGA_R'\n");
    print ("  -l  : logfile [STDERR]\n");
    print ("  -v  : Verbose [OFF]\n");
    print ("Dependencies: R 3.2.5 or later. The following R packages have to be installed: Bioconductor, DESeq2.\n");
    print ("\n");
    exit;
} 
# Usage
#
# Open input
#

if (not defined($Getopt::Std::opt_i)){
  # Read from standard input
    print "You must specify a filename\n";
    exit;
} 
else{
    open(INP,"<", $Getopt::Std::opt_i) || die ("can't open file $Getopt::Std::opt_i: $!");
}
#
# If not file name is given, use standard output
#

if (defined($Getopt::Std::opt_A)){
    open(INPA, "<", $Getopt::Std::opt_A) || die ("can't open file $Getopt::Std::opt_A: $!");
} else {
    print "No chosen accession list selected for comparison [-A]\n";
    exit;
}

if (defined($Getopt::Std::opt_B)){
    open(INPB, "<", $Getopt::Std::opt_B) || die ("can't open file $Getopt::Std::opt_B: $!");
} else {
    print "No chosen accession list selected for comparison [-B]\n";
    exit;
}

if (defined($Getopt::Std::opt_f)){
    $percent_present_in_samples = $Getopt::Std::opt_f;
}

if (defined($Getopt::Std::opt_R)){
    system("mkdir -p $R_dir");
    open(R, ">", "$R_dir/$R_file") || die ("can't open file $Getopt::Std::opt_R: $!\n");
    open(CONDITION, ">", "$R_dir/$condition_matrix") || die "$!\n";
    open(DIVERSA, ">", "$R_dir/$divers_file_A") || die "$!\n";
    open(DIVERSB, ">", "$R_dir/$divers_file_B") || die "$!\n";
}

if (defined($Getopt::Std::opt_r)){
    system("mkdir -p $R_dir");
    open(READS, ">", "$R_dir/$r_file") || die ("can't open file $Getopt::Std::opt_R: $!\n");
    open(CONDITION, ">", "$R_dir/$condition_matrix") || die "$!\n";
    open(DIVERSA, ">", "$R_dir/$divers_file_A") || die "$!\n";
    open(DIVERSB, ">", "$R_dir/$divers_file_B") || die "$!\n";
}

if (defined($Getopt::Std::opt_v)){
    $reRun = 1;
}

if (defined($Getopt::Std::opt_l)){
    open(LOG,">", "$R_dir/$Getopt::Std::opt_l");
}

if (defined($Getopt::Std::opt_v)){
    $verbose = 1;
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
    print LOG "# working dir: $R_dir\n\n";
}

while(defined(my $line = <INPA>)){
    chomp $line;
    $set1{$line} = {};
}
close INPA;

while(defined(my $line = <INPB>)){
    chomp $line;
    $set2{$line} = {};
}
close INPB;

while (defined(my $line = <INP>)) {
    chomp $line;
  next if $line =~ m/^#/;

    my @tmp = split("\t", $line);
    my $study_sample = $tmp[0];
    my $clade = $tmp[2];
    my $tot_readCount = $tmp[1];
    my $abundance = $tmp[3];
    my $reads = $tmp[11];

    if (exists $set1{$study_sample}) {
	$set{$clade}{$study_sample}{study_sample} = $study_sample;
	$set{$clade}{$study_sample}{clade} = $clade;
	$total_reads{$study_sample}{tot_readCount} = $tot_readCount;
	$total_reads{$study_sample}{nn} += $reads;
	$set{$clade}{$study_sample}{abundance} += $abundance;
	$set{$clade}{$study_sample}{reads} += $reads;
	if (! exists($stored{$study_sample})) {

	    push(@stored_study_sample1,$study_sample);
	    $num_samples_set1++;
	    $stored{$study_sample} = {};
	}
    }
    elsif (exists $set2{$study_sample}) {
	$set{$clade}{$study_sample}{study_sample} = $study_sample;
	$set{$clade}{$study_sample}{clade} = $clade;
	$total_reads{$study_sample}{tot_readCount} = $tot_readCount;
	$total_reads{$study_sample}{nn} += $reads;
	$set{$clade}{$study_sample}{abundance} += $abundance;
	$set{$clade}{$study_sample}{reads} += $reads;
	if (! exists($stored{$study_sample})) {

	    push(@stored_study_sample2,$study_sample);
	    $num_samples_set2++;
	    $stored{$study_sample} = {};
	}
    }
}

my $i = 0;
my $percentage_threshold_A = $percent_present_in_samples * $num_samples_set1;
my $percentage_threshold_B = $percent_present_in_samples * $num_samples_set2;

if (defined($Getopt::Std::opt_R)) {
    my $i = 0;
    my $y = 0;
    foreach my $sample_A (@stored_study_sample1) {
	$i++;
	if ($i == 1) {
	    print R "# clade\t$sample_A\t";
	    print DIVERSA "\t$sample_A\t";
	}
	elsif ($i < (scalar(@stored_study_sample1))) {
	    print R "$sample_A\t";
	    print DIVERSA "$sample_A\t";
	} else {
	    print R "$sample_A";
	    print DIVERSA "$sample_A";
	}
    }
    foreach my $sample_B (@stored_study_sample2) {
	$y++;
	if ($y == 1) {
	    print R "\t$sample_B\t";
	    print DIVERSB "\t$sample_B\t";
	}
	elsif ($y < (scalar(@stored_study_sample2))) {
	    print R "$sample_B\t";
	    print DIVERSB "$sample_B\t";
	}
	else {
	    print R "$sample_B";
	    print DIVERSB "$sample_B";
	}
    }
    print R "\n";
    print DIVERSA "\n";
    print DIVERSB "\n";
}

if (defined($Getopt::Std::opt_r)) {
    my $i = 0;
    my $y = 0;
    foreach my $sample_A (@stored_study_sample1) {
	$i++;
	if ($i == 1) {
	    print READS "# clade\t$sample_A\t";
	    print DIVERSA "\t$sample_A\t" if not defined $Getopt::Std::opt_R;
	}
	elsif ($i < (scalar(@stored_study_sample1))) {
	    print READS "$sample_A\t";
	    print DIVERSA "$sample_A\t" if not defined $Getopt::Std::opt_R;
	}
	else {
	    print READS "$sample_A";
	    print DIVERSA "$sample_A" if not defined $Getopt::Std::opt_R;
	}
    }
    foreach my $sample_B (@stored_study_sample2) {
	$y++;
	if ($y == 1) {
	    print READS "\t$sample_B\t";
	    print DIVERSB "\t$sample_B\t" if not defined $Getopt::Std::opt_R;
	}
	elsif ($y < (scalar(@stored_study_sample2))) {
	    print READS "$sample_B\t";
	    print DIVERSB "$sample_B\t" if not defined $Getopt::Std::opt_R;
	}
	else {
	    print READS "$sample_B";
	    print DIVERSB "$sample_B" if not defined $Getopt::Std::opt_R;
	}
    }
    print READS "\n";
    print DIVERSA "\n" if not defined $Getopt::Std::opt_R;
    print DIVERSB "\n" if not defined $Getopt::Std::opt_R;
}

if (defined($Getopt::Std::opt_R) or defined($Getopt::Std::opt_r)) {
    print CONDITION "# sample\tcondition\n";
    foreach my $sample_A (@stored_study_sample1) {
	print CONDITION "$sample_A\tA\n";
    }
    foreach my $sample_B (@stored_study_sample2) {
	print CONDITION "$sample_B\tB\n";
    }
}

if (defined($Getopt::Std::opt_H)) {
    if (defined ($Getopt::Std::opt_o) or defined ($Getopt::Std::opt_O)) {
	print "# clade\t".join("\t", @stored_study_sample1)."\t".join("\t", @stored_study_sample2)."\n";
    }
}

foreach my $clade (sort keys %set) {
    my @samples_1 = ();
    my @samples_2 = ();
    my @samples_1_reads = ();
    my @samples_2_reads = ();
    my $sample_1;
    my $sample_2;
    my $sample_1_reads;
    my $sample_2_reads;
    my $present_A = 0;
    my $present_B = 0;

    foreach my $sample_A (@stored_study_sample1) {

	if (exists $set{$clade}{$sample_A}{abundance}) {
	    $sample_1 = $set{$clade}{$sample_A}{abundance};
	    $sample_1_reads = $set{$clade}{$sample_A}{reads};
	    $richness{$sample_A}{richness} += 1;
	    $present_A++;
	} else {
	    $sample_1 = 0;
	    $sample_1_reads = 0;
	}
	push(@samples_1, $sample_1);
	push(@samples_1_reads, $sample_1_reads);
    }

    foreach my $sample_B (@stored_study_sample2) {

	if (exists $set{$clade}{$sample_B}{abundance}) {
	    $sample_2 = $set{$clade}{$sample_B}{abundance};
	    $sample_2_reads = $set{$clade}{$sample_B}{reads};
	    $richness{$sample_B}{richness} += 1;
	    $present_B++;
	} else {
	    $sample_2 = 0;
	    $sample_2_reads = 0;
	}
	push(@samples_2, $sample_2);
	push(@samples_2_reads, $sample_2_reads);
    }

    if ($present_A > $percentage_threshold_A or $present_B > $percentage_threshold_B) {
	if (defined($Getopt::Std::opt_R)) {
	    print R "$clade\t".join("\t", @samples_1)."\t".join("\t", @samples_2)."\n";
	}

	if (defined($Getopt::Std::opt_r)) {
	    print READS "$clade\t".join("\t", @samples_1_reads)."\t".join("\t", @samples_2_reads)."\n";
	}

	if (defined($Getopt::Std::opt_O)) {
	    print "$clade\t".join("\t", @samples_1)."\t".join("\t", @samples_2)."\n";
	}
	if (defined($Getopt::Std::opt_o)) {
	    print "$clade\t".join("\t", @samples_1_reads)."\t".join("\t", @samples_2_reads)."\n";
	}
    }
}

if (defined($Getopt::Std::opt_r) or defined($Getopt::Std::opt_R)) {
    print DIVERSA "# alpha_diversity\t";
    print DIVERSB "# alpha_diversity\t";
    my $c = 0;
    foreach my $el (@stored_study_sample1) {
	$c++;
	if ($c < scalar(@stored_study_sample1)) {
	    print DIVERSA "$richness{$el}{richness}\t";
	}
	else {
	    print DIVERSA "$richness{$el}{richness}"; 
	}
    }
    $c = 0;
    foreach my $el (@stored_study_sample2) {
	$c++;
	if ($c < scalar(@stored_study_sample2)) {
	    print DIVERSB "$richness{$el}{richness}\t";
	}
	else {
	    print DIVERSB "$richness{$el}{richness}"; 
	} 
    }
}
print LOG "# Doing batch.R\n" if defined $Getopt::Std::opt_l;

#sorting the total number of reads for each sample according their position in the matrix

foreach my $sample (@stored_study_sample1, @stored_study_sample2) {
    my $tot_readCounts;
    my $nn;
    if (exists $total_reads{$sample}{tot_readCount}) {
	$tot_readCounts = $total_reads{$sample}{tot_readCount};
	push(@total_reads, $tot_readCounts);
    }
    if (exists $total_reads{$sample}{nn}) {
	$nn = $total_reads{$sample}{nn};
	push(@nn, $nn);
    }
}

#print the nn

if (defined($Getopt::Std::opt_r)) {
    print READS "nn\t";
    my $difference;
    for (my $i = 0; $i < scalar(@total_reads); $i++) {
	$difference = $total_reads[$i] - $nn[$i];
	print READS "$difference\t" if $i < scalar(@total_reads);
	print READS "$difference" if $i == scalar(@total_reads);
    }
}

if (defined($Getopt::Std::opt_o)) {
    my $difference;
    print "nn\t";
    for (my $i = 0; $i < scalar(@total_reads); $i++) {
	$difference = $total_reads[$i] - $nn[$i];
	print "$difference\t" if $i < scalar(@total_reads);
	print "$difference" if $i == scalar(@total_reads);
    }
}

if (defined($Getopt::Std::opt_i)) {
    close INP;
}

if (defined($Getopt::Std::opt_R)) {
    close R;
    close CONDITION;
    close DIVERS;
}

if (defined($Getopt::Std::opt_r)) {
    close READS;
    close CONDITION;
    close DIVERS;
}

if (defined($Getopt::Std::opt_l)){
    close LOG;
}
=pod
    print join("\n", @nn);
print "\n\n\n";
print join("\n", @total_reads);
=cut

if ($reRun == 0) {
    system("cp /services/tools/omega/1.1/$batch_R $R_dir");
}

system("$prog_R CMD BATCH $R_dir/$batch_R");
