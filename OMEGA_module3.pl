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
use MysqlTunnel;
use DBI;
#use DBD::mysql;
use perl;
use strict;

################################################################################
# Global variables
################################################################################
my $cfg = "$ENV{HOME}/.my.cnf";
my $database;
my $private = "_private";
open(ACCESS_INFO, "<", $cfg) or die "$cfg file notfound\n";
while (defined(my $line = <ACCESS_INFO>)){
    chomp $line;
    if ($line =~ m/^user\s+=\s+(\S+)/) {
        $database = "$1$private";
    }
}
close ACCESS_INFO;
my $study_accession;
my $db_id = "ALL";
my $sample_accession;
my  $cmd = <<TXT;
SELECT CONCAT_WS('_', a.study_acc, a.sample_acc), c.read_count, a.db_id, a.seq_id, a.size, a.nucleotide_count, 
       a.covered_pos, a.read_count, a.read_count_unique, a.edit_distance, a.description, 
       b.superfamily, b.superfamilyTax, b.phylum, b.phylumTax, b.class, b.classTax, 
       b.order, b.orderTax, b.family, b.familyTax, b.genus, b.genusTax, b.species, 
       b.speciesTax, b.strain, b.strainTax FROM omega_seq_data AS a, omega_taxonomy AS b, omega_study_sample AS c 
 WHERE a.study_acc=b.study_acc AND a.sample_acc=b.sample_acc AND a.db_id=b.db_id AND 
        a.seq_id=b.seq_id AND a.study_acc=c.study_acc AND a.sample_acc=c.sample_acc
TXT

# Default parameters
*LOG=*STDERR;
my $verbose=0;
#
# Process command line
#
getopts('hi:vl:d:s:D:o:SL:')||Usage();
#
# Usage
#
if (defined($Getopt::Std::opt_h)||defined($Getopt::Std::opt_h)){
  # Print help message
  Usage();
}

sub Usage {
  print ("Usage: $0 [-h] [-i 'study accession'] [-s 'sample accession'] [-d 'database ID'] [-D mysql database name] [-S retrieve db_id's] [-l log file] [-v verbose]\n");
  print ("Description:\n");
  print ("$0 - pipeline for autOmated MEtaGenomic Analysis (OMEGA) - module 3: select and retrieve information requests from MySQL database.\n");
  print ("\n");
  print ("Options:\n");
  print ("  -h  : display this message*\n");
  print ("  -i  : input study accession\n");
  print ("  -o  : output filename [STDOUT]\n");
  print ("  -s  : input sample accession* [ALL]\n");
  print ("  -L  : input list of sample accessions separated by new line\n");
  print ("  -d  : input annotation database ID* [$db_id]\n");
  print ("  -D  : mysql database name [$database]\n");
  print ("  -S  : retrieve annotation database id's from MySQL databse and exits [OFF]\n");
  print ("  -l  : logfile [STDERR]\n");
  print ("  -v  : Verbose [off]\n");
  print ("*Note: If more than one study accession, sample accession or database ID is selected, they must be placed within single quotes and separated by space ('db1 db2 db3')\n");
  print ("\n");
 exit;
} # Usage

#
# Open input
#
if (not defined($Getopt::Std::opt_D)) {
    print "MySQL database not found\n";
    exit;
} else {
    $database = $Getopt::Std::opt_D;
}

if (defined($Getopt::Std::opt_S)){
    my $tunnel = new MysqlTunnel;
    my ($dsn, $user, $password) = $tunnel->getparm($database);
    my $dbh = DBI->connect($dsn, $user, $password, { RaiseError => 1, AutoCommit => 0 });
    my $db_retrieve = "SELECT db_id FROM omega_databasesRef";
    print LOG "Doing: $db_retrieve\n" if ($verbose);
    my $fetch = $dbh->prepare("$db_retrieve");
    $fetch->execute();
    my @retrieved;
    while (my @row_array = $fetch->fetchrow_array()){
	print "@row_array". "\n";
    }
    $dbh->disconnect();
    $tunnel->close;
    exit;
}

if (defined($Getopt::Std::opt_i)){
  # Read from standard input
  $study_accession = $Getopt::Std::opt_i;
  my @tmp = split(/\s+/, $study_accession);
  my $i = 0;
  $cmd .= " AND (";
  my $size = scalar(@tmp);
  foreach my $el (@tmp) {
      $i++;
      $cmd .= " a.study_acc = '$el'";
      if (($size > 1) and ($i < $size)){
      	$cmd .= " OR";
      }
  }
  $cmd .= " )";
  } else {
  print "Accession not found.\n";
  exit;
}

if (defined($Getopt::Std::opt_s)){

  $sample_accession = $Getopt::Std::opt_s;
  my @tmp = split(/\s+/, $sample_accession);    
  my $i = 0;
  $cmd .= " AND (";
  my $size = scalar(@tmp);
  foreach my $el (@tmp) {
     $i++;
     $cmd .= " a.sample_acc = '$el'";
     if (($size > 1) and ($i<$size)) {
         $cmd .= " OR";
     }
  }
  $cmd .= " )";

} elsif (defined($Getopt::Std::opt_L)){
	open(IN, "<", $Getopt::Std::opt_L) or die "$!\n";
	my @tmp;
	my $i = 0;
	$cmd .= " AND (";
	while(defined(my $line = <IN>)) {
		chomp $line;
		$cmd .= " a.sample_acc = '$line'";
		if (not eof) {
			$cmd .= " OR";
		}
	}
	$cmd .= " )";
	close IN;
}

if (defined($Getopt::Std::opt_d)){

  $db_id = $Getopt::Std::opt_d;
  my @tmp = split(/\s+/, $db_id);
  my $i = 0;
  $cmd .= " AND (";
  my $size = scalar(@tmp);
  foreach my $el (@tmp) {
     $i++;
     $cmd .= " a.db_id = '$el'";
     if (($size > 1) and ($i<$size)) {
         $cmd .= " OR";
     }
  }
  $cmd .= " )";
}

if (defined($Getopt::Std::opt_v)){
  $verbose=1;
}

if (defined($Getopt::Std::opt_l)){
  open(LOG,">", $Getopt::Std::opt_l);
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
    print LOG "# working dir: $thisDir\n\n";
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
print LOG "MySQL version: @status\n\n" if $verbose;
$statement->finish();

print LOG "Doing: $cmd\n\n" if ($verbose);
my $fetch = $dbh->prepare ("$cmd");
$fetch->execute();

my $tot_entries = 0;
if (defined($Getopt::Std::opt_o)) {
  open(OUT, ">", $Getopt::Std::opt_o) or die "Can not open file for output: $!\n";
  my @all_rows;
  while (my @row_array = $fetch->fetchrow_array()) {
    print OUT join("\t", @row_array). "\n";
    $tot_entries += 1;
  }

  close OUT;
} else {
    my @all_rows;
    while (my @row_array = $fetch->fetchrow_array()) {
      print join("\t", @row_array). "\n";
      $tot_entries += 1;
    }
  }

$dbh->disconnect();
$tunnel->close;

print LOG "$tot_entries entries were retrieved from MySQL database $database\n\n" if ($verbose);
print LOG "Done!\n" if ($verbose);
