#!/usr/bin/env perl
###############################################################
#
# Global variables
#
###############################################################
my $command="$0 ";
my $i=0;

while (defined($ARGV[$i])){
    $command .= "$ARGV[$i] ";
    $i++;
}

my $ppn = 14;
my $sequencing_mode = 2;
my $study_accession;
my @accession_list;
my $download_txt;
my $fastq_paired_1;
my $fastq_paired_2;
my $fastq_single;
my $download_paired_1_cmd;
my $download_paired_2_cmd;
my $download_single_cmd;
my $sample_acc;
my $read_acc;
my $F_MGlist_paired = "F_list";
my $R_MGlist_paired = "R_list";
my $Paired_MGlists = "PairedListofLists";
my $Single_MGlist = "SingleList";
my $Single_MGlists = "SingleListofLists";
my $fastq_dir = "fastq";
my $MGmapper_dir = "MGmapper";
my $down_stat = "FastqDownload.stat";
my $walltime = 14;
my $memory_usage = 60;
my $workingDir = "ENA";
my $prog_wget = "wget";
my $gb = 'gb';
my $parameterFile;
my $bestModeparam;
my $fullModeparam;
my $verbose=0;

use lib '/cm/local/apps/environment-modules/3.2.10/init/';
use Getopt::Std;
use Cwd;
use strict;
use perl;

module("load tools moab torque perl/5.20.2 mgmapper/2.4");

# Default parameters
#
# Process command line
#
getopts('hi:vl:d:p:m:L:w:e:P:C:F:')||Usage();
#
# Usage
#
if (defined($Getopt::Std::opt_h)){
  # Print help message
    Usage();
}

sub Usage {
    print ("Usage: $0 [-h] [-i name] [-L name] [-d name] [-m number] [-p number] ... [-v]\n");
    print ("Description:\n");
    print ("$0 - pipeline for autOmated MEtaGenomic Analysis (OMEGA) - module 1: Download ENA study accession files and process them through MGmapper.\n");
    print ("\n");
    print ("Options:\n");
    print ("  -h  : display this message\n");
    print ("  -i  : input study accession [STDIN]\n");
    print ("  -L  : input list of tab-separated study accessions\n");
    print ("  -d  : output directory [ENA]\n");
    print ("  -P  : Absolute path to parameter file with custom databases to be used for the alignment [/services/tools/mgmapper/2.4/databases.txt]\n");
    print ("  -C  : Only if -P option is defined: map reads in Best mode against these databases (comma-separated numbers - the order matters)\n");
    print ("  -F  : Only if -P option is defined: map reads in Full mode against these databases (comma-separated numbers)\n");
    print ("  -p  : processors per node [$ppn]\n");
    print ("  -m  : sequencing mode: 1 = single, 2 = paired, 3 = both [paired]\n");
    print ("  -w  : walltime (days) for each submitted MGmapper job [$walltime]\n");
    print ("  -e  : memory (Gb) for each submitted MGmapper job [$memory_usage]\n");
    print ("  -l  : logfile [STDERR]\n");
    print ("  -v  : Verbose [off]\n");
    print ("\n");
    exit;
} # Usage

#
# Open input
#

if (defined($Getopt::Std::opt_i)){
  # Read from standard input
    $study_accession = $Getopt::Std::opt_i;
} elsif (defined($Getopt::Std::opt_L)) {
    open (FILE, "<", $Getopt::Std::opt_L) or die "Can't open $Getopt::Std::opt_L. $!\n";
    while (defined(my $line = <FILE>)) {
	chomp $line;
	my @tmp_array = split("\t", $line);
	push(@accession_list, @tmp_array);
    }
    close FILE;
} else {
    print "Accession not found.\n";
    exit;
}

if (defined($Getopt::Std::opt_v)){
    $verbose=1;
}

if (defined($Getopt::Std::opt_d)) {
    $workingDir = "$Getopt::Std::opt_d";
}

if (! -d $workingDir) {
    system("mkdir -p $workingDir");
    print LOG "## Directory $workingDir was created\n" if ($verbose);
}

if (defined($Getopt::Std::opt_p)) {
    $ppn = $Getopt::Std::opt_p;
}

if (defined($Getopt::Std::opt_m)) {
    $sequencing_mode = $Getopt::Std::opt_m;
}

if (defined($Getopt::Std::opt_w)) {
    $walltime = $Getopt::Std::opt_w;
}

if (defined($Getopt::Std::opt_e)) {
    $memory_usage = $Getopt::Std::opt_e;
}

if (defined($Getopt::Std::opt_P)) {
    $parameterFile = $Getopt::Std::opt_P;
    if (! -e $parameterFile){
    print LOG "can't open file $parameterFile\n";
    print LOG "Done!\n";
    die;
    }
    if (defined($Getopt::Std::opt_C)) {
	$bestModeparam = $Getopt::Std::opt_C;
    }
    if (defined($Getopt::Std::opt_F)) {
	$fullModeparam = $Getopt::Std::opt_F;
    }
}

if (defined($Getopt::Std::opt_l)){
    open(LOG,">", "$workingDir/$Getopt::Std::opt_l") or die "Can't open $workingDir/$study_accession/$Getopt::Std::opt_l. $!\n";
} else {
    *LOG=*STDERR;
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
    print LOG "# working dir: $workingDir\n\n";
}
############################################################################################################
#
# Download the text file from study accession number as input and extract fastq/fastq.gz files from the list
# Running MGmapper
#
############################################################################################################
if (defined($Getopt::Std::opt_i)) {
    print LOG "# Downloading $study_accession list with file addresses.\n" if ($verbose);
    acc_list_download($study_accession);
    print LOG "# Downloading fastq/fastq.gz files from $study_accession list and creating input lists for MGmapper.\n" if ($verbose);
    remove_duplicates($study_accession);
    acc_fastq_download($study_accession);
    print LOG "# Doing: MGmapper\n" if ($verbose);
    create_MGlist($study_accession);
    MGmapper_cmd($study_accession);
} elsif (defined($Getopt::Std::opt_L)) {
    foreach my $el (@accession_list) {
	print LOG "# Downloading $el list with file addresses\n" if ($verbose);
	acc_list_download($el);
	print LOG "# Downloading fastq/fastq.gz files from $el list and creating input lists for MGmapper.\n" if ($verbose);
	remove_duplicates($el);
	acc_fastq_download($el) if (-d "$workingDir/$el");
	print LOG "# Doing: MGmapper\n" if ($verbose);
	create_MGlist($el);
	MGmapper_cmd($el);
    }
}
################################################################################
#
# Subroutinies
#
################################################################################
sub acc_list_download {
    my ($accession) = @_;
    if (-d "$workingDir/$accession" and $accession ne "") {
	print "This study accession already exists. Updating list.\n";
    } else {
	system("mkdir -p $workingDir/$accession");
	print LOG "##  Directory $workingDir/$accession was created.\n" if ($verbose);
    }
    $download_txt .= "$prog_wget -q -nc -O $workingDir/$accession/$accession ";
    $download_txt .= '"http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=';
    $download_txt .= $accession;
    $download_txt .= '&result=read_run&fields=study_accession,secondary_study_accession,sample_accession,secondary_sample_accession,experiment_accession,run_accession,tax_id,scientific_name,instrument_model,library_layout,fastq_ftp,fastq_galaxy,submitted_ftp,submitted_galaxy,cram_index_ftp,cram_index_galaxy&download=txt"';
    print LOG "# Doing: $download_txt\n" if ($verbose);
    system($download_txt);
    if ((-s "$workingDir/$accession/$accession") <= 300) {
	system ("rm -r $workingDir/$accession") if (-e "$workingDir/$accession") ;
	print "Not a valid accession number!\n";
    }
}
sub remove_duplicates {
    my ($accession) = @_;
    system("find $workingDir/$accession -type f -name $F_MGlist_paired -delete -print &>/dev/null");
    system("find $workingDir/$accession -type f -name $R_MGlist_paired -delete -print &>/dev/null");
    system("find $workingDir/$accession -type f -name $Single_MGlist -delete -print &>/dev/null");
}

sub acc_fastq_download {
    my ($accession) = @_;
    open (FILE, "<", "$workingDir/$accession/$accession") or die "Can't open $workingDir/$accession/$accession. $!\n";
    open (OUT, ">", "$workingDir/$accession/$down_stat") or die "Can't open $workingDir/$accession/$down_stat. $!\n";
    my $downloaded = 0;
    my $failed = 0;
    my @tmp_failed;
    while (defined(my $line = <FILE>)) {
	chomp $line;
	if ($line =~ m/PAIRED/) {
	    if ($sequencing_mode == 2 or $sequencing_mode == 3) {
		if ($line =~ m/;(ftp\S+_1\.fastq\.gz);(ftp\S+_2\.fastq\.gz)/ or $line =~ m/(ftp\S+_1\.fastq\.gz);(ftp\S+_2\.fastq\.gz)/) {
		    my @tmp_sample_acc = split("\t", $line);
		    $sample_acc = $tmp_sample_acc[2];
		    $read_acc = $tmp_sample_acc[5];
		    my $i = "_1.fastq.gz";
		    my $y = "_2.fastq.gz";
		    system("mkdir -p $workingDir/$accession/$sample_acc");
		    system("mkdir -p $workingDir/$accession/$sample_acc/$fastq_dir");
		    print LOG "## Directory  $workingDir/$accession/$sample_acc/$fastq_dir was created.\n" if ($verbose) and (! -d "$workingDir/$accession/$sample_acc/$fastq_dir");
		    $download_paired_1_cmd = "$prog_wget -q -nc -o $workingDir/$accession/$sample_acc/wget.log -O $workingDir/$accession/$sample_acc/$fastq_dir/$read_acc$i $1";
		    $download_paired_2_cmd = "$prog_wget -q -nc -o $workingDir/$accession/$sample_acc/wget.log -O $workingDir/$accession/$sample_acc/$fastq_dir/$read_acc$y $2";   
		    system("echo `pwd`/$workingDir/$accession/$sample_acc/$fastq_dir/$read_acc$i >> $workingDir/$accession/$sample_acc/$F_MGlist_paired");
		    system("echo `pwd`/$workingDir/$accession/$sample_acc/$fastq_dir/$read_acc$y >> $workingDir/$accession/$sample_acc/$R_MGlist_paired");
		    print LOG "# Doing: $download_paired_1_cmd\n# Doing: $download_paired_2_cmd\n" if ($verbose);
		    system($download_paired_1_cmd);
		    system($download_paired_2_cmd);
		    $downloaded += 1;
		} elsif ($line =~ m/;(ftp\S+_1\.fastq);(ftp\S+_2\.fastq)/ or $line =~ m/(ftp\S+_1\.fastq);(ftp\S+_2\.fastq)/) {
		    my @tmp_sample_acc = split("\t", $line);
		    $sample_acc = $tmp_sample_acc[2];
		    $read_acc = $tmp_sample_acc[5];
		    my $i = "_1.fastq";
		    my $y = "_2.fastq";
		    system("mkdir -p $workingDir/$accession/$sample_acc");
		    system("mkdir -p $workingDir/$accession/$sample_acc/$fastq_dir");
		    print LOG "## Directory  $workingDir/$accession/$sample_acc/$fastq_dir was created.\n" if ($verbose) and (! -d "$workingDir/$accession/$sample_acc/$fastq_dir");
		    $download_paired_1_cmd = "$prog_wget -q -nc -o $workingDir/$accession/$sample_acc/wget.log -O $workingDir/$accession/$sample_acc/$fastq_dir/$read_acc$i $1";
		    $download_paired_2_cmd = "$prog_wget -q -nc -o $workingDir/$accession/$sample_acc/wget.log -O $workingDir/$accession/$sample_acc/$fastq_dir/$read_acc$y $2";
		    system("echo `pwd`/$workingDir/$accession/$sample_acc/$fastq_dir/$read_acc$i >> $workingDir/$accession/$sample_acc/$F_MGlist_paired");
		    system("echo `pwd`/$workingDir/$accession/$sample_acc/$fastq_dir/$read_acc$y >> $workingDir/$accession/$sample_acc/$R_MGlist_paired");
		    print LOG "# Doing: $download_paired_1_cmd\n# Doing: $download_paired_2_cmd\n" if ($verbose);
		    system($download_paired_1_cmd);
		    system($download_paired_2_cmd);
		    $downloaded += 1;
		} else {
		    $failed += 1;
		    my @tmp_sample_acc = split("\t", $line);
		    my $failed_sample  = "Inconsistent data in: Study_acc $tmp_sample_acc[1]\tSample_acc $tmp_sample_acc[2]\tRead_acc $tmp_sample_acc[5]";
		    push(@tmp_failed, $failed_sample);
		}
	    }
	} elsif ($line =~ m/SINGLE/) {
	    if ($sequencing_mode == 1 or $sequencing_mode == 3) {
		if ($line =~ m/(ftp\S+fastq\.gz)/) {
		    my @tmp_sample_acc = split("\t", $line);
		    $sample_acc = $tmp_sample_acc[2];
		    $read_acc = $tmp_sample_acc[5];
		    my $i = ".fastq.gz";
		    system("mkdir -p $workingDir/$accession/$sample_acc");
		    system("mkdir -p $workingDir/$accession/$sample_acc/$fastq_dir");
		    print LOG "## Directory $workingDir/$accession/$sample_acc/$fastq_dir was created.\n" if ($verbose) and (! -d "$workingDir/$accession/$sample_acc/$fastq_dir");
		    $download_single_cmd = "$prog_wget -q -nc -o $workingDir/$accession/$sample_acc/wget.log -O $workingDir/$accession/$sample_acc/$fastq_dir/$read_acc$i $1";
		    system("echo `pwd`/$workingDir/$accession/$sample_acc/$fastq_dir/$read_acc$i >> $workingDir/$accession/$sample_acc/$Single_MGlist");
		    print LOG "# Doing: $download_single_cmd\n" if ($verbose);
		    system($download_single_cmd);
		    $downloaded += 1;
		} elsif ($line =~ m/(ftp\S+fastq)/) {
		    my @tmp_sample_acc = split("\t", $line);
		    $sample_acc = $tmp_sample_acc[2];
		    $read_acc = $tmp_sample_acc[5];
		    my $i = ".fastq";
		    system("mkdir -p $workingDir/$accession/$sample_acc");
		    system("mkdir -p $workingDir/$accession/$sample_acc/$fastq_dir");
		    print LOG "## Directory  $workingDir/$accession/$sample_acc/$fastq_dir was created.\n" if ($verbose) and (! -d "$workingDir/$accession/$sample_acc/$fastq_dir");
		    $download_single_cmd = "$prog_wget -q -nc -o $workingDir/$accession/$sample_acc/wget.log -O $workingDir/$accession/$sample_acc/$fastq_dir/$read_acc$i $1";
		    system("echo `pwd`/$workingDir/$accession/$sample_acc/$fastq_dir/$read_acc$i >> $workingDir/$accession/$sample_acc/$Single_MGlist");
		    print LOG "# Doing: $download_single_cmd\n" if ($verbose);
		    system($download_single_cmd);
		    $downloaded += 1;
		} else {
		    $failed += 1;
		    $failed += 1;
		    my @tmp_sample_acc = split("\t", $line);
		    my $failed_sample  = "Inconsistent data in: Study_acc $tmp_sample_acc[1]\tSample_acc $tmp_sample_acc[2]\tRead_acc $tmp_sample_acc[5]";
		    push(@tmp_failed, $failed_sample);
		}
	    }
	}
    }
    close FILE;
    print OUT "Downloaded read accessions: $downloaded\n";
    print OUT "Failed to download: $failed\n";
    foreach my $x (@tmp_failed) {
	print OUT "$x\n";
    }
    close OUT;
}

sub create_MGlist {
    my ($accession) = @_;
    
    if ($sequencing_mode == 2 or $sequencing_mode == 3) {
	system("find `pwd`/$workingDir -type f -name $F_MGlist_paired | sort > $workingDir/$accession/tmp_paired_1");
	system("find `pwd`/$workingDir -type f -name $R_MGlist_paired | sort > $workingDir/$accession/tmp_paired_2");
	system("paste -d '\t' $workingDir/$accession/tmp_paired_1 $workingDir/$accession/tmp_paired_2 > $workingDir/$accession/$Paired_MGlists");
    } elsif ($sequencing_mode == 1 or $sequencing_mode == 3) {
	system("find `pwd`/$workingDir -type f -name $Single_MGlist > $workingDir/$accession/$Single_MGlists");
    }
    system("rm -f $workingDir/$accession/tmp_paired_1 $workingDir/$accession/tmp_paired_2");
}

sub MGmapper_cmd {
    my ($accession) = @_;
    if ($sequencing_mode == 2) {
	open (LIST ,"<", "$workingDir/$accession/$Paired_MGlists") or die "Can't open $workingDir/$accession/$Paired_MGlists. $!\n";
	while (defined(my $line = <LIST>)) {
	    chomp $line;
	    $line =~ s/\.\///g;
	    my $tmp_sample_acc = $1 if $line =~ m/\/$accession\/(\S+)\//;
	    my @tmp_list = split("\t", $line);
	    my $re="$tmp_sample_acc".".e";
	    my $ro="$tmp_sample_acc".".o";
	    my $cmd;
	    if (not defined($Getopt::Std::opt_P)) {
		$cmd = "xqsub -V -ro $ro -re $re -N $tmp_sample_acc -d $workingDir -l nodes=1:ppn=$ppn,mem=$memory_usage$gb,walltime=$walltime:00:00:00 -de MGmapper_PE.pl -I $tmp_list[0] -J $tmp_list[1] -d $accession/$tmp_sample_acc/$MGmapper_dir -C 1,5,3,4,2,11,8,9,16,7,10,19 -F 6,12 -c $ppn -M";
	    } else {
		$cmd = "xqsub -V -ro $ro -re $re -N $tmp_sample_acc -d $workingDir -l nodes=1:ppn=$ppn,mem=$memory_usage$gb,walltime=$walltime:00:00:00 -de MGmapper_PE.pl -I $tmp_list[0] -J $tmp_list[1] -d $accession/$tmp_sample_acc/$MGmapper_dir -c $ppn -M -P $parameterFile";
		if (defined($Getopt::Std::opt_C)) {
		    $cmd .= " -C $bestModeparam";
		}
		if (defined($Getopt::Std::opt_F)) {
		    $cmd .= " -F $fullModeparam";
		}
	    }
	    if (-d "$workingDir/$accession/$tmp_sample_acc/$MGmapper_dir") {
		print LOG "# Directory $MGmapper_dir already exists\n" if ($verbose);
	    } else {
		print LOG "# $cmd\n" if ($verbose);
		system("$cmd");
	    }
	}
	close LIST;
    } elsif ($sequencing_mode == 1) {
	open (LIST ,"<", "$workingDir/$accession/$Single_MGlists") or die "Can't open $workingDir/$accession/$Single_MGlists. $!\n";
	while (defined(my $line = <LIST>)) {
	    chomp $line;
	    $line =~ s/\.\///g;
	    my $tmp_sample_acc = $1 if $line =~ m/\/$accession\/(\S+)\//;
	    my $re="$tmp_sample_acc".".e";
	    my $ro="$tmp_sample_acc".".o";
	    my $cmd;
	    if (not defined($Getopt::Std::opt_P)) {
		$cmd = "xqsub -V -ro $ro -re $re -N $tmp_sample_acc -d $workingDir -l nodes=1:ppn=$ppn,mem=$memory_usage$gb,walltime=$walltime:00:00:00 -de MGmapper_SE.pl -f $accession/$Single_MGlists -d $accession/$tmp_sample_acc/$MGmapper_dir -C 1,5,3,4,2,11,8,9,16,7,10,19 -F 6,12 -c $ppn -M";
	    } else {
		$cmd = "xqsub -V -ro $ro -re $re -N $tmp_sample_acc -d $workingDir -l nodes=1:ppn=$ppn,mem=$memory_usage$gb,walltime=$walltime:00:00:00 -de MGmapper_SE.pl -f $accession/$Single_MGlists -d $accession/$tmp_sample_acc/$MGmapper_dir -c $ppn -M -P $parameterFile";
		if (defined($Getopt::Std::opt_C)) {
		    $cmd .= " -C $bestModeparam";
		}
		if (defined($Getopt::Std::opt_F)) {
		    $cmd .= " -F $fullModeparam";
		}
	    }
	    if (-d "$workingDir/$accession/$tmp_sample_acc/$MGmapper_dir") {
		print LOG "# Directory $MGmapper_dir already exists\n" if ($verbose);
	    } else {
		print LOG "# $cmd\n" if ($verbose);
		system("$cmd");
	    }
	}
	close LIST;
    } elsif ($sequencing_mode == 3) {
	open (LIST ,"<", "$workingDir/$accession/$Paired_MGlists") or die "Can't open $workingDir/$accession/$Paired_MGlists. $!\n";
	while (defined(my $line = <LIST>)) {
	    chomp $line;
	    $line =~ s/\.\///g;
	    my $tmp_sample_acc = $1 if $line =~ m/\/$accession\/(\S+)\//;
	    my @tmp_list = split("\t", $line);
	    my $re="$tmp_sample_acc".".e";
	    my $ro="$tmp_sample_acc".".o";
	    my $cmd;
	    if (not defined($Getopt::Std::opt_P)) {
		$cmd = "xqsub -V -ro $ro -re $re -N $tmp_sample_acc -d $workingDir -l nodes=1:ppn=$ppn,mem=$memory_usage$gb,walltime=$walltime:00:00:00 -de MGmapper_PE.pl -I $tmp_list[0] -J $tmp_list[1] -d $accession/$tmp_sample_acc/$MGmapper_dir -C 1,5,3,4,2,11,8,9,16,7,10,19 -F 6,12 -c $ppn -M";
	    } else {
		$cmd = $cmd = "xqsub -V -ro $ro -re $re -N $tmp_sample_acc -d $workingDir -l nodes=1:ppn=$ppn,mem=$memory_usage$gb,walltime=$walltime:00:00:00 -de MGmapper_PE.pl -I $tmp_list[0] -J $tmp_list[1] -d $accession/$tmp_sample_acc/$MGmapper_dir -c $ppn -M -P $parameterFile";
		if (defined($Getopt::Std::opt_C)) {
		    $cmd .= " -C $bestModeparam";
		}
		if (defined($Getopt::Std::opt_F)) {
		    $cmd .= " -F $fullModeparam";
		}
	    }
	    if (-d "$workingDir/$accession/$tmp_sample_acc/$MGmapper_dir") {
		print LOG "# Directory $MGmapper_dir already exists\n" if ($verbose);
	    } else {
		print LOG "# $cmd\n" if ($verbose);
		system("$cmd");
	    }
	}
	close LIST;
	open (LIST ,"<", "$workingDir/$accession/$Single_MGlists") or die "Can't open $workingDir/$accession/$Single_MGlists. $!\n";
	while (defined(my $line = <LIST>)) {
	    chomp $line;
	    $line =~ s/\.\///g;
	    my $tmp_sample_acc = $1 if $line =~ m/\/$accession\/(\S+)\//;
	    my $re="$tmp_sample_acc".".e";
	    my $ro="$tmp_sample_acc".".o";
	    my $cmd;
	    if (not defined($Getopt::Std::opt_P)) {
		$cmd = "xqsub -V -ro $ro -re $re -N $tmp_sample_acc -d $workingDir -l nodes=1:ppn=$ppn,mem=$memory_usage$gb,walltime=$walltime:00:00:00 -de MGmapper_SE.pl -f $accession/$Single_MGlists -d $accession/$tmp_sample_acc/$MGmapper_dir -C 1,5,3,4,2,11,8,9,16,7,10,19 -F 6,12 -c $ppn -M";
	    } else {
		$cmd = "xqsub -V -ro $ro -re $re -N $tmp_sample_acc -d $workingDir -l nodes=1:ppn=$ppn,mem=$memory_usage$gb,walltime=$walltime:00:00:00 -de MGmapper_SE.pl -f $accession/$Single_MGlists -d $accession/$tmp_sample_acc/$MGmapper_dir -c $ppn -M -P $parameterFile";
		if (defined($Getopt::Std::opt_C)) {
		    $cmd .= " -C $bestModeparam";
		}
		if (defined($Getopt::Std::opt_F)) {
		    $cmd .= " -F $fullModeparam";
		}
	    }
	    if (-d "$workingDir/$accession/$tmp_sample_acc/$MGmapper_dir") {
		print LOG "# Directory $MGmapper_dir already exists\n" if ($verbose);
	    } else {
		print LOG "# $cmd\n" if ($verbose);
		system("$cmd");
	    }
	}
	close LIST;
    }
}
