#!/usr/bin/env perl
use strict;

#---------------------------------------------------------------------------
# pQuest.pl
# A software/script to prepare the folder structure for xQuest analysis 
# and to extract features using pseudoSH.pl
# Execute pQuest.pl -help to display information and usage options.
#---------------------------------------------------------------------------

#---------------------------------------------------------------------------
# Licence
# This software is licensed under the Apache License, Version 2.0.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
# This software and associated documentation is provided "AS IS",
# without warranties or conditions of any kind.
# See the License for the specific language governing permissions and
# limitations under the License.
#---------------------------------------------------------------------------

use File::Basename;
use File::Copy "cp";
use Getopt::Long;
use File::Spec;
use Cwd;
use Data::Dumper;

# Include modules dir as lib that is relative to the Script path
##########################################################
use FindBin;                           ##
use lib "$FindBin::Bin/../modules";    ##
##########################################################
#use Environment;
#my $env = Environment->new;
use Cwd;
#my $dev = $env->get_path('xquest_stable');
my $dev = getcwd;
my ( $verbose, $help, $listfile, $pseudoSH, $mzxml, $suffix, $link );
### Standard Parameters
$listfile = "files";
$pseudoSH = 1;
my $paramdef   = "param.def";
my $xquestdef  = "xquest.def";
my $xmmdef     = "xmm.def";
my $cwd        = getcwd;
my $path2NAS   = $cwd."/";
my $OS         = $^O;            # get the current OS
my @mzXMLfiles = ();
my $sep        = ".";
GetOptions(
			'list=s'      => \$listfile,
			'path=s'      => \$path2NAS,
			'xquestdef=s' => \$xquestdef,
			'xmmdef=s'    => \$xmmdef,
			'pseudosh=i'  => \$pseudoSH,
			'verbose'     => \$verbose,
			'help'        => \$help,
			'mzxml=s'     => \$mzxml,
			'sep=s'       => \$sep,
			'link=s'		=> \$link,
);
my $bin        = dirname(__FILE__); #File::Spec->catfile( $dev, "bin" ); # Old version was not as versatible
my $version    = "1.3";
my $scriptinfo = {};
$scriptinfo->{'version'} = $version;
$scriptinfo->{'author'}  = "Thomas Walzthoeni (1)";
my $affiliation = "(1) ETH Zurich, Institute of Molecular Systems Biology, Wolfgang Pauli-Str. 16, CH-8093 Zurich";
$scriptinfo->{'affi'} = $affiliation;
my $mailto = "walzthoeni\@imsb.biol.ethz.ch";
$scriptinfo->{'mailto'} = $mailto;
$scriptinfo->{'clog'}->{'1.2'} = "Windows functionality added.\n";
$scriptinfo->{'clog'}->{'1.3'} = "Linking of folders functionality added.\n";
unless ($listfile)
{
	print "No file with list of the basenames of the mzXML files defined use -list option\n";
	&usage();
}
&usage() if $help;

unless ( -e $xquestdef )
{
	die "Error: $xquestdef does not exist $!";
}

unless ( -e $xmmdef )
{
	die "Error: $xmmdef does not exist $!";
}

open INPUT, "<$listfile" or die $!;
while ( my $line = <INPUT> )
{
	chomp($line);
	my @results = split( " ", $line );
	if ($line)
	{
		push @mzXMLfiles, $results[0];
	}
}
close(INPUT);

my $mzsuffix;
foreach my $mzxmlfile (@mzXMLfiles)
{
	my $centroidmzxmlfile;
	my $mzxmlfilename;
	### generate the mzXML filename
	$mzxmlfilename = $mzxmlfile . $sep . "mzXML";
	$centroidmzxmlfile = join "", $path2NAS, $mzxmlfilename;
	if ( $mzxml eq "p" )
	{
		$mzxmlfilename = $mzxmlfile . $sep . "p.mzXML";
		$centroidmzxmlfile = join "", $path2NAS, $mzxmlfilename;
	}
	if ( $mzxml eq "c" )
	{
		$mzxmlfilename = $mzxmlfile . $sep . "c.mzXML";
		$centroidmzxmlfile = join "", $path2NAS, $mzxmlfilename;
	}
	my ( $status, $msg ) = check_file($centroidmzxmlfile);
	if ( $status != 1 )
	{
		die $msg;
	}
	print "Preparing directory structure for $mzxmlfilename\n";
	my $basename = basename($mzxmlfile);
	if ( -e $basename )
	{
		print "Directory exists, should I remove the directory (yes/no)?";
		my $a = <STDIN>;    # Get input
		chop($a);
		if ( $a eq "yes" )
		{
			my $cmd = "rm -rfv $basename";
			print $cmd. "\n";
			system($cmd);
		} else
		{
			print "->Exit\n";
			next;
		}
	}
	
	## Create the directory
	mkdir $basename or die $!;
	mkdir $basename . "/MASTER_RUN" or die $!;
	## Copy the deffiles
	cp( $xmmdef,    $basename );
	cp( $xquestdef, $basename );
	### Check which OS is used:
	if ( $OS =~ m/MSWin/ )
	{
		print "Copy $centroidmzxmlfile to $basename\n";
		## mzxmlfile needs to be copied to the directory
		cp( $centroidmzxmlfile, $basename );
		
		if ($link)
		{
		my $cmd = "cp $link $basename/$basename";
		print "Copy Folder $link to subfolder: $cmd\n";
		mkdir $basename . "/$basename" or die $!;
		cp( $link, $basename."/".$basename );				
		}
		
		
	}
	if ( $OS eq "linux" )
	{
		print "Generate symlink:  ln -s $centroidmzxmlfile $basename/$mzxmlfilename\n";
		my $cmd = "ln -s $centroidmzxmlfile $basename/$mzxmlfilename";
		print("-->$cmd\n");
		system($cmd);
		
		if ($link)
		{
		my $cmd = "ln -s $link $basename/$basename";
		print "Link command: $cmd\n";
		system($cmd);
		}
		
		
		
	}
	if ( $pseudoSH == 1 )
	{
		my $runcmd = "$bin/pseudoSH.pl -mz $centroidmzxmlfile -out $basename/MASTER_RUN/MASTER_RUN.txt";
		print( "Execute pseudSH:" . $runcmd . "\n" );
		my $result = system($runcmd);
		if ($result)
		{
			die "pseudoSH execution returned an error: error code: $result\n";
		}
	}
}


#===  FUNCTION  ================================================================
#  NAME:        check_file
#  PURPOSE:     Check file for existence and readability
#  DESCRIPTION: Check file for existence and readability
#  PARAMETERS:  $filename (including path if script is not executed in the folder)
#  RETURNS:     void
#===============================================================================
sub check_file
{
	my $filename = shift;
	my $verbose  = shift;
	my $msg;
	my $status;

	# -e is for exists, -r for readable
	unless ( ( -e $filename ) && ( -r $filename ) )
	{
		$msg    = "Error: Cannot find/read the file $filename.\n";
		$status = -1;
	} else
	{
		$msg    = "Checking file $filename ok.\n";
		$status = 1;
	}
	return ( $status, $msg );
}

sub usage()
{
print "
	SOFTWARE: ", basename($0), " version $version
	
	AUTHOR: Thomas Walzthoeni

	INFORMATION: A software/script to prepare the folder structure for xQuest analysis 
	and to extract features using pseudoSH.pl.
 	
 	USAGE: ", basename($0), " -Option [Parameter]

	REQUIRED OPTIONS: Values in [] are the default values that are used if the option is not provided
	-list [files] a file with basenames of the mzXML input files (one per line). e.g.: File AXX-YYYYY.mzXML --> use: AXX-YYYYY 
	-path [./] path to the mzXML files, include ending \"/\" e.g. /path/to/mzxmls/
	On windows, the mzXML files are then copied to the individual search directories, on linux a softlink is created.
	-pseudosh [1] extracts a feature map with the program pseudoSH.pl (map is stored in the MASTER_RUN folders), to omit this option use -pseudosh 0
	
	-xquestdef [xquest.def] xquest definition filename, must be present in the directory
	-xmmdef [xmm.def] xmm definition filename, must be present in the directory
	-sep[.] delimiter for mzXML files: for e.g. AXX-YYYYY_c.mzXML use -sep _c. as delimiter (and AXX-YYYYY as the basename in the list file)

	OTHER OPTIONS:
	-mzxml []: define p or c to use .p. or .c. as delimiter for the filenames: eg. if you have filenames like: AXX-YYYYY.c.mzXML
 	-help print this help
 	-link [path] Link (on Linux) or copy (on Windows) a folder into the subfolders

	EXAMPLE: 
	",basename($0)," -list files -path /path/to/mzxmls/
 		
	";	

	exit;
}
