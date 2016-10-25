#!/usr/bin/env perl
my $command="$0 ";
my $i=0;
###################################### 
# WORKING ENVIRONMENT
######################################
my $workingDir;
my $mysql_dir = "MySQL";
######################################
# MYSQL FILENAMES
######################################
my $cfg = "$ENV{HOME}/.my.cnf";
my $database;
my $private = '_private';
open(ACCESS_INFO, "<", $cfg) or die "$cfg file not found\n";
while (defined(my $line = <ACCESS_INFO>)){
    chomp $line;
    if ($line =~ m/^user\s+=\s+(\S+)/) {
	$database = "$1$private";
    }
}
close ACCESS_INFO;
my $study_tabFile = "omega_study.tab";
my $sample_tabFile = "omega_study_sample.tab";
my $study_sample_db = "omega_study_sample_databasesRef.tab";
my $db_ref = "omega_databasesRef.tab";
my $seq_data = "omega_seq_data.tab";
my $taxonomy_file = "omega_taxonomy.tab";
my $runInfo_file = "omega_runInfo.tab";
my $study_description = "NULL";
######################################
# MYSQL TABLES
######################################
my $table_study = "omega_study";
my $table_study_sample = "omega_study_sample";
my $table_ssdb = "omega_study_sample_databasesRef";
my $table_dbRef = "omega_databasesRef";
my $table_seq_data = "omega_seq_data";
my $table_tax = "omega_taxonomy";
my $table_runInfo = "omega_run_info";
###################################### 
# MGmapper dirs and files
######################################
my $MGmapper_dir = "MGmapper";
my $misc_dir = "misc";
my $abundance_file = "abundance.databases.txt";
my $db_count_file = "dbCount.tab";
my $run_file = "run.info";
my $prefix_db_filename = "stat.";
my $suffix_db_filename = ".annot";
###################################### 
# GLOBAL VARIABLES
######################################
my $study_accession;
my @accession_list;
my $sample_acc;
my $sth;
#####################################
# total entries for each table
#####################################
my $tot_table_study = 0;
my $tot_table_sample = 0;
my $tot_table_study_dbref = 0;
my $tot_db_ref = 0;
my $tot_seq_data = 0;
my $tot_tax = 0;
my $tot_run_info = 0;

while (defined($ARGV[$i])){
  $command .= "$ARGV[$i] ";
  $i++;
}

use lib '/cm/local/apps/environment-modules/3.2.10/init/';
use perl;
use Getopt::Std;
use Cwd;
use MysqlTunnel;
use DBI;
use strict;


# Default parameters
*LOG=*STDERR;
my $verbose=0;
#
# Process command line
#
getopts('hi:vl:d:L:D:m:')||Usage();
#
# Usage
#
if (defined($Getopt::Std::opt_h)||defined($Getopt::Std::opt_h)){
  # Print help message
  Usage();
}

sub Usage {
  print ("Usage: $0 [-h] [-i name] [-L name] [-d name] ... [-v]\n");
  print ("Description:\n");
  print ("$0 - pipeline for autOmated MEtaGenomic Analysis (OMEGA) - module 2: Create MySQL tables with 'omega_' prefix. Sort and insert MGmapper .annot data into MySQL database.\n");
  print ("\n");
  print ("Options:\n");
  print ("  -h  : display this message\n");
  print ("  -i  : input study accession\n");
  print ("  -L  : input list of tab separated study accessions\n");
  print ("  -d  : working directory*\n");
  print ("  -m  : mysql output directory [$mysql_dir]\n");
  print ("  -D  : MySQL database [$database]\n");
  print ("  -l  : logfile [STDERR]\n");
  print ("  -v  : Verbose [off]\n");
  print ("* Working directory must be the path to the directory where the study accessions are found [i.e. $ENV{HOME}/ENA]\n");
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

  open (FILE, "<", $Getopt::Std::opt_L) or die "$!\n";
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

if (defined($Getopt::Std::opt_D)) {
    $database = $Getopt::Std::opt_D;
} else {
    print "MySQL database not specified\n";
    exit;
}

if (defined($Getopt::Std::opt_l)){
    open(LOG,">", $Getopt::Std::opt_l);
} else {
    *LOG=*STDERR;
}

if (defined($Getopt::Std::opt_m)) {
    $mysql_dir = $Getopt::Std::opt_m;
}

if (defined($Getopt::Std::opt_d)) {
    $workingDir = "$Getopt::Std::opt_d";
} else {
    print "Working directory not found\n";
    exit;
}

if (! -d $workingDir) {
  system("mkdir -p $workingDir");
  print LOG "## Directory $workingDir was created\n" if ($verbose);
}

###############################################################################
#
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

################################################################################
#
# Creating the input lists for the MySQL database
#
################################################################################

if (defined($Getopt::Std::opt_i)) {
  if (-d "$workingDir/$study_accession") {
    if (! -d "$workingDir/$study_accession/$mysql_dir") {
      system("mkdir -p $workingDir/$study_accession/$mysql_dir");
      print LOG "## Directory $mysql_dir was created\n" if ($verbose);
    }
    createFiles($study_accession);
  } else {
    print LOG "Study accession $study_accession do not exist\n" if ($verbose);
  }
  
} elsif (defined($Getopt::Std::opt_L)) {

  foreach my $acc (@accession_list) {
    if (-d "$workingDir/$acc") {
      if (! -d "$workingDir/$acc/$mysql_dir") {
        system("mkdir -p $workingDir/$acc/$mysql_dir");
        print LOG "## Directory $mysql_dir was created\n" if ($verbose);
      }
      createFiles($acc);
    } else {
      print LOG "Study accession $acc do not exist\n" if ($verbose);
    }
  }
}

################################################################################
#
# Connect to MySQL database
#
################################################################################

my $tunnel = new MysqlTunnel;
my ($dsn, $user, $password) = $tunnel->getparm($database);
my $dbh = DBI->connect($dsn, $user, $password, { RaiseError => 1, AutoCommit => 0 });

my $query = "SELECT VERSION()";
my $statement = $dbh->prepare($query);
$statement->execute();
my @status = $statement->fetchrow_array();
print LOG "MySQL Version: @status\n" if ($verbose);

################################################################################
#
# Create tables if not exist
#
################################################################################

my $sql_table_study = "CREATE TABLE IF NOT EXISTS $table_study (study_acc VARCHAR(20), description TEXT, PRIMARY KEY (study_acc));";
my $sql_table_sample = "CREATE TABLE IF NOT EXISTS $table_study_sample (study_acc VARCHAR(100), sample_acc VARCHAR(100), read_count INT(20), PRIMARY KEY (study_acc,sample_acc));";
my $sql_table_ssdb = "CREATE TABLE IF NOT EXISTS $table_ssdb (study_acc VARCHAR(100), sample_acc VARCHAR(100), db_id VARCHAR(100), db_version VARCHAR(100), PRIMARY KEY (study_acc, sample_acc, db_id));";
my $sql_table_dbRef = "CREATE TABLE IF NOT EXISTS $table_dbRef (db_id VARCHAR(100), db_version VARCHAR(100), PRIMARY KEY (db_id));";
my $sql_table_seq_data = "CREATE TABLE IF NOT EXISTS $table_seq_data (study_acc VARCHAR(100), sample_acc VARCHAR(100), db_id VARCHAR(100), seq_id VARCHAR(100), size INT(20), nucleotide_count INT(20), covered_pos INT(20), read_count INT(20), read_count_unique INT(20), edit_distance INT, description VARCHAR(255), PRIMARY KEY (study_acc, sample_acc, db_id, seq_id));";
my $sql_table_tax = "CREATE TABLE IF NOT EXISTS $table_tax (study_acc VARCHAR(100), sample_acc VARCHAR(100), db_id VARCHAR(100), seq_id VARCHAR(100), superfamily VARCHAR(100), superfamilyTax INT, phylum VARCHAR(100), phylumTax INT, class VARCHAR(100), classTax INT, `order` VARCHAR(100), orderTax INT, family VARCHAR(100), familyTax INT, genus VARCHAR(100), genusTax INT, species VARCHAR(100), speciesTax INT, strain VARCHAR(100), strainTax INT, PRIMARY KEY (study_acc, sample_acc, db_id, seq_id));";
my $sql_table_run_info = "CREATE TABLE IF NOT EXISTS $table_runInfo (study_acc VARCHAR(20), sample_acc VARCHAR(100), seq_mode VARCHAR(20), MGmapper_ver VARCHAR(20), date VARCHAR(100), work_dir TEXT, command TEXT, PRIMARY KEY (study_acc, sample_acc));";
################################################################################
#
# Open each table and insert elements from files
#
################################################################################
if (defined($Getopt::Std::opt_i)){
  
  $study_accession = $Getopt::Std::opt_i;
  SortandInsert($study_accession);
 
} elsif (defined($Getopt::Std::opt_L)) {

    foreach my $acc (@accession_list) {
  
      SortandInsert($acc);
  
    }
  }
  
$statement->finish();
$dbh->disconnect;
$tunnel->close;

print LOG "Done!\n" if ($verbose);
################################################################################
#
# Subroutines
#
################################################################################


sub createFiles {

  my ($accession) = @_;
  open(STUDY, ">", "$workingDir/$study_accession/$mysql_dir/$study_tabFile") or die "$!\n";
  $study_description =~ s/'/ /g;
  print STUDY "$study_accession\t$study_description";
  close STUDY;

  open(FILE, "<", "$workingDir/$study_accession/$study_accession") or die "$!\n";
  open(SAMPLE, ">", "$workingDir/$study_accession/$mysql_dir/$sample_tabFile") or die "$!\n";
  open(SSDB, ">", "$workingDir/$study_accession/$mysql_dir/$study_sample_db") or die "$!\n";
  open(DBREF, ">", "$workingDir/$study_accession/$mysql_dir/$db_ref") or die "$!\n";
  open(SEQDATA, ">", "$workingDir/$study_accession/$mysql_dir/$seq_data") or die "$!\n";
  open(TAX, ">", "$workingDir/$study_accession/$mysql_dir/$taxonomy_file") or die "$!\n";
  open(INFO, ">", "$workingDir/$study_accession/$mysql_dir/$runInfo_file") or die "$!\n";
  my $line = <FILE>;
  while (defined($line = <FILE>)) {
    chomp $line;
    my $count = '';
    my @tmp_array = split("\t", $line);
    $sample_acc = $tmp_array[2];
    if (-e "$workingDir/$study_accession/$sample_acc/$MGmapper_dir/$abundance_file" and ! -z "$workingDir/$study_accession/$sample_acc/$MGmapper_dir/$abundance_file") {

      open(COUNTS, "<", "$workingDir/$study_accession/$sample_acc/$MGmapper_dir/$abundance_file");
      my $count_line = <COUNTS>;
      chomp $count_line;
      my @tmp_counts = split("\t", $count_line);
      $count = $tmp_counts[3];
      close COUNTS;
      
    } 
    else {
      print LOG "$abundance_file file for sample $sample_acc do not exist or has zero size\n" if ($verbose);
    }

    if (-e "$workingDir/$study_accession/$sample_acc/$MGmapper_dir/$misc_dir/$run_file" and ! -z "$workingDir/$study_accession/$sample_acc/$MGmapper_dir/$misc_dir/$run_file") {
      open(RUN, "<", "$workingDir/$study_accession/$sample_acc/$MGmapper_dir/$misc_dir/$run_file");
      my ($sequencing_mode, $MGversion, $command, $workdir, $date) = ('','','','','');
      while (defined(my $runline = <RUN>)) {
        chomp $runline;
        my @tmp_run = split("\t", $runline);
        if ($runline =~ m/^Version\s+(\S+)_(\S+)/) {
          $sequencing_mode = $1;
          $MGversion = $2;
        }
        $command = $tmp_run[1] if $runline =~ m/^Command/;
	$workdir = $tmp_run[1] if $runline =~ m/^WorkDir/;
	$date = $tmp_run[1] if $runline =~ m/^Date/;
      }
      print INFO "$study_accession\t$sample_acc\t$sequencing_mode\t$MGversion\t$command\t$workdir\t$date\n";
    } 
    else {
      print LOG "$run_file file for sample $sample_acc do not exist or has zero size\n" if ($verbose);
    }      
    
    if (-e "$workingDir/$study_accession/$sample_acc/$MGmapper_dir/$misc_dir/$db_count_file" and ! -z "$workingDir/$study_accession/$sample_acc/$MGmapper_dir/$misc_dir/$db_count_file"){
      open(DBCOUNT, "<", "$workingDir/$study_accession/$sample_acc/$MGmapper_dir/$misc_dir/$db_count_file");
      my $dbname = '';
      my $dbversion = '';
      my $dbcount = <DBCOUNT>;
      while (defined($dbcount = <DBCOUNT>)) {
        
        chomp $dbcount;
        my @tmp_dbcount = split("\t", $dbcount);
        $dbname = $tmp_dbcount[0];
	$dbversion = $tmp_dbcount[1];
        my $path_to_db_file_name = "$workingDir/$study_accession/$sample_acc/$MGmapper_dir/$misc_dir/$prefix_db_filename$dbname$suffix_db_filename";
        
        if (-e $path_to_db_file_name and ! -z $path_to_db_file_name) {
          print SSDB "$study_accession\t$sample_acc\t$dbname\t$dbversion\n" if $tmp_dbcount[2] != 0;
          open(IN, "<", $path_to_db_file_name);
          
          while (defined(my $annot_line = <IN>)) {
            chomp $annot_line;
	    $annot_line =~ s/'/ /g; # MySQL do not allow single quotes for input data so I had to escape it
            my @tmp_annot_data = split("\t", $annot_line);
            print SEQDATA "$study_accession\t$sample_acc\t$dbname\t$tmp_annot_data[1]\t$tmp_annot_data[2]\t$tmp_annot_data[3]\t$tmp_annot_data[4]\t$tmp_annot_data[5]\t$tmp_annot_data[6]\t$tmp_annot_data[7]\t$tmp_annot_data[8]\n";
            print TAX "$study_accession\t$sample_acc\t$dbname\t$tmp_annot_data[1]\t$tmp_annot_data[9]\t$tmp_annot_data[10]\t$tmp_annot_data[11]\t$tmp_annot_data[12]\t$tmp_annot_data[13]\t$tmp_annot_data[14]\t$tmp_annot_data[15]\t$tmp_annot_data[16]\t$tmp_annot_data[17]\t$tmp_annot_data[18]\t$tmp_annot_data[19]\t$tmp_annot_data[20]\t$tmp_annot_data[21]\t$tmp_annot_data[22]\t$tmp_annot_data[23]\t$tmp_annot_data[24]\n";
          }
	  close IN;
        }
        print DBREF "$dbname\t$dbversion\n";
        
      }
      close DBREF;
      close DBCOUNT;
    } else {
      print LOG "$db_count_file for sample $sample_acc do not exist or has zero size\n" if ($verbose);
    }

    print SAMPLE "$study_accession\t$sample_acc\t$count\n" if (-d "$workingDir/$study_accession/$sample_acc" and ! -z "$workingDir/$study_accession/$sample_acc");
  }
  close FILE;
  close SAMPLE;
  close SSDB;
  close SEQDATA;
  close TAX;
  close INFO;
}

sub SortandInsert {

  my ($study_accesion) = @_;

  open(STUDY, "<", "$workingDir/$study_accession/$mysql_dir/$study_tabFile") or die "$!\n";
  $sth = $dbh -> prepare($sql_table_study) or die ("Error in SQL: '$sql_table_study'\n");
  $sth->execute() or die ("Error in SQL\n");
  while (defined(my $line = <STUDY>)) {
    chomp $line;
    my @tmp_words = split("\t", $line);
    my $sql_study_acc = $tmp_words[0];
    my $sql_description = $tmp_words[1];
    $sql_table_study = "REPLACE INTO $table_study (study_acc, description) VALUES ('$sql_study_acc', '$sql_description')";
    $sth = $dbh -> prepare($sql_table_study) or die ("Error in SQL: '$sql_table_study'\n");
    $sth->execute() or die ("Error in SQL\n");
    $tot_table_study += 1;
  }
  close STUDY;
  print LOG "$tot_table_study entrie(s) inserted into MySQL table $table_study\n" if ($verbose);

  open(SAMPLE, "<", "$workingDir/$study_accession/$mysql_dir/$sample_tabFile") or die "$!\n";
  $sth = $dbh -> prepare($sql_table_sample) or die ("Error in SQL: '$sql_table_sample'\n");
  $sth->execute() or die ("Error in SQL\n");
  while (defined(my $line = <SAMPLE>)) {
    chomp $line;
    my @tmp_words = split("\t", $line);
    my $sql_study_acc = $tmp_words[0];
    my $sql_sample_acc = $tmp_words[1];
    my $sql_readCount = $tmp_words[2];
    $sql_table_sample = "REPLACE INTO $table_study_sample (study_acc, sample_acc, read_count) VALUES ('$sql_study_acc', '$sql_sample_acc', '$sql_readCount')";
    $sth = $dbh -> prepare($sql_table_sample) or die ("Error in SQL: '$sql_table_sample'\n");
    $sth->execute() or die ("Error in SQL\n");
    $tot_table_sample += 1;
  }
  close SAMPLE;
  print LOG "$tot_table_sample entrie(s) inserted into MySQL table $table_study_sample\n" if ($verbose);

  open(SSDB, "<", "$workingDir/$study_accession/$mysql_dir/$study_sample_db") or die "$!\n";
  $sth = $dbh -> prepare($sql_table_ssdb) or die ("Error in SQL: '$sql_table_ssdb'\n");
  $sth->execute() or die ("Error in SQL\n");
  while (defined(my $line = <SSDB>)) {
    chomp $line;
    my @tmp_words = split("\t", $line);
    my $sql_study_acc = $tmp_words[0];
    my $sql_sample_acc = $tmp_words[1];
    my $sql_db_id = $tmp_words[2];
    $sql_table_ssdb = "REPLACE INTO $table_ssdb (study_acc, sample_acc, db_id) VALUES ('$sql_study_acc', '$sql_sample_acc', '$sql_db_id')";
    $sth = $dbh -> prepare($sql_table_ssdb) or die ("Error in SQL: '$sql_table_ssdb'\n");
    $sth->execute() or die ("Error in SQL\n");
    $tot_table_study_dbref += 1;
  }
  close SSDB;
  print LOG "$tot_table_study_dbref entrie(s) inserted into MySQL table $table_ssdb\n" if ($verbose);

  open(DBREF, "<", "$workingDir/$study_accession/$mysql_dir/$db_ref") or die "$!\n";
  $sth = $dbh -> prepare($sql_table_dbRef) or die ("Error in SQL: '$sql_table_dbRef'\n");
  $sth->execute() or die ("Error in SQL\n");
  while (defined(my $line = <DBREF>)) {
    chomp $line;
    my @tmp_words = split("\t", $line);
    my $sql_db_id = $tmp_words[0];
    $sql_table_dbRef = "REPLACE INTO $table_dbRef (db_id) VALUES ('$sql_db_id')";
    $sth = $dbh -> prepare($sql_table_dbRef) or die ("Error in SQL: '$sql_table_dbRef'\n");
    $sth->execute() or die ("Error in SQL\n");
    $tot_db_ref += 1;
  }
  close DBREF;
  print LOG "$tot_db_ref entrie(s) inserted into MySQL table $table_dbRef\n" if ($verbose);

  open(SEQDATA, "<", "$workingDir/$study_accession/$mysql_dir/$seq_data") or die "$!\n";
  $sth = $dbh -> prepare($sql_table_seq_data) or die ("Error in SQL: '$sql_table_seq_data'\n");
  $sth->execute() or die ("Error in SQL\n");
  while (defined(my $line = <SEQDATA>)) {
    chomp $line;
    my @tmp_words = split("\t", $line);
    my ($sql_study_acc, $sql_sample_acc, $sql_db_id, $sql_seq_id, $sql_size, $sql_nuc_count, $sql_covered_pos, $sql_read_count, $sql_uniq_read, $sql_distance, $sql_description) = ($tmp_words[0], $tmp_words[1], $tmp_words[2], $tmp_words[3], $tmp_words[4], $tmp_words[5], $tmp_words[6], $tmp_words[7], $tmp_words[8], $tmp_words[9], $tmp_words[10]);
    $sql_table_seq_data = "REPLACE INTO $table_seq_data (study_acc, sample_acc, db_id, seq_id, size, nucleotide_count, covered_pos, read_count, read_count_unique, edit_distance, description) VALUES ('$sql_study_acc', '$sql_sample_acc', '$sql_db_id', '$sql_seq_id', '$sql_size', '$sql_nuc_count', '$sql_covered_pos', '$sql_read_count', '$sql_uniq_read', '$sql_distance', '$sql_description')";
    $sth = $dbh -> prepare($sql_table_seq_data) or die ("Error in SQL: '$sql_table_seq_data'\n");
    $sth->execute() or die ("Error in SQL\n");
    $tot_seq_data += 1;
  }
  close SEQDATA;
  print LOG "$tot_seq_data entrie(s) inserted into MySQL table $table_seq_data\n" if ($verbose);

  open(TAX, "<", "$workingDir/$study_accession/$mysql_dir/$taxonomy_file") or die "$!\n";
  $sth = $dbh -> prepare($sql_table_tax) or die ("Error in SQL: '$sql_table_tax'\n");
  $sth->execute() or die ("Error in SQL\n");
  while (defined(my $line = <TAX>)) {
    chomp $line;
    my @tmp_words = split("\t", $line);
    my ($sql_study_acc, $sql_sample_acc, $sql_db_id, $sql_seq_id, $sql_superfamily, $sql_superfamilyTax, $sql_phylum, $sql_phylumTax, $sql_class, $sql_classTax, $sql_order, $sql_orderTax, $sql_family, $sql_familyTax, $sql_genus, $sql_genusTax, $sql_species, $sql_speciesTax, $sql_strain, $sql_strainTax) = ($tmp_words[0], $tmp_words[1], $tmp_words[2], $tmp_words[3], $tmp_words[4], $tmp_words[5], $tmp_words[6], $tmp_words[7], $tmp_words[8], $tmp_words[9], $tmp_words[10], $tmp_words[11], $tmp_words[12], $tmp_words[13], $tmp_words[14], $tmp_words[15], $tmp_words[16], $tmp_words[17], $tmp_words[18], $tmp_words[19]);
    $sql_table_tax = "REPLACE INTO $table_tax (study_acc, sample_acc, db_id, seq_id, superfamily, superfamilyTax, phylum, phylumTax, class, classTax, `order`, orderTax, family, familyTax, genus, genusTax, species, speciesTax, strain, strainTax) VALUES ('$sql_study_acc', '$sql_sample_acc', '$sql_db_id', '$sql_seq_id', '$sql_superfamily', '$sql_superfamilyTax', '$sql_phylum', '$sql_phylumTax', '$sql_class', '$sql_classTax', '$sql_order', '$sql_orderTax', '$sql_family', '$sql_familyTax', '$sql_genus', '$sql_genusTax', '$sql_species', '$sql_speciesTax', '$sql_strain', '$sql_strainTax')";
    $sth = $dbh -> prepare($sql_table_tax) or die ("Error in SQL: '$sql_table_tax'\n");
    $sth->execute() or die ("Error in SQL\n");
    $tot_tax += 1;
  }
  close TAX;
  print LOG "$tot_tax entrie(s) inserted into MySQL table $table_tax\n" if ($verbose);

  open(RUNINFO, "<", "$workingDir/$study_accession/$mysql_dir/$runInfo_file") or die "$!\n";
  $sth = $dbh -> prepare($sql_table_run_info) or die ("Error in SQL: '$sql_table_run_info'\n");
  $sth->execute() or die ("Error in SQL\n");
  while (defined(my $line = <RUNINFO>)) {
    chomp $line;
    my @tmp_words = split("\t", $line);
    my ($sql_study_acc, $sql_sample_acc, $sql_seq_mode, $sql_MGmapper_ver, $sql_date, $sql_work_dir, $sql_command) = ($tmp_words[0], $tmp_words[1], $tmp_words[2], $tmp_words[3], $tmp_words[4], $tmp_words[5], $tmp_words[6]);
    $sql_table_run_info = "REPLACE INTO $table_runInfo (study_acc, sample_acc, seq_mode, MGmapper_ver, date, work_dir, command) VALUES ('$sql_study_acc', '$sql_sample_acc', '$sql_seq_mode', '$sql_MGmapper_ver', '$sql_date', '$sql_work_dir', '$sql_command')";
    $sth = $dbh -> prepare($sql_table_run_info) or die ("Error in SQL: '$sql_table_run_info'\n");
    $sth->execute() or die ("Error in SQL\n");
    $tot_run_info += 1;
  }
  close RUNINFO;
  print LOG "$tot_run_info entrie(s) inserted into MySQL table $table_runInfo\n" if ($verbose);

}

