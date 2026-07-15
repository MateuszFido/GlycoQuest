#!/usr/bin/env perl
use strict;
#---------------------------------------------------------------------------
# mergexml.pl
# A software/script to to merge xQuest XML files.
# Execute mergexml.pl -help to display information and usage options.
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

# path to the libs on brutus
use FindBin;
#use lib "$FindBin::Bin/../../perl5";
use File::Basename;
use Getopt::Long;
use File::Spec;
use Cwd;
use Data::Dumper;
use XML::TreeBuilder;
use XML::Element;
use File::Copy;

#---------------------------------------------------------------------------
#  Variables
#---------------------------------------------------------------------------
my ( $filesin, $outfilename, $help, $version, $list, $verbose, $resdir, $nospecxml, $xquestdef, $xmmdef );

#---------------------------------------------------------------------------
#  Default values
#---------------------------------------------------------------------------
$version     = "1.2";
$outfilename = "merged_xquest.xml";
$resdir      = "xmlmergeroutput";
$xquestdef   = "xquest.def";
$xmmdef   = "xmm.def";
$list        = "resultdirectories_fullpath";

#---------------------------------------------------------------------------
#  Options that are passed by arguments
#---------------------------------------------------------------------------
GetOptions(
			'files=s'   => \$filesin,                   ## String of files "file1 file2"
			'out=s'     => \$outfilename,
			'list=s'    => \$list,
			'help'      => \$help,
			'v'         => \$verbose,
			'resdir=s'  => \$resdir,
			'nospecxml' => \$nospecxml,
			'xquestdef' => \$xquestdef,
			'xmmdef' => \$xmmdef,
);
if ($verbose) { $verbose = 1 }

#---------------------------------------------------------------------------
#  Show usage
#---------------------------------------------------------------------------
&usage() unless ( $filesin || $list );
## Undefine $list if filenames are provided
if ($filesin)
{
	$list = undef;
}
&usage() if $help;

#---------------------------------------------------------------------------
#  Check files from input string $filesin or $list
#---------------------------------------------------------------------------
my @files;
my @specfiles;
my $dir = getcwd;
if ($filesin)
{
	my @filenames = split( / /, $filesin );
	foreach my $file (@filenames)
	{
		my $filename = File::Spec->catfile( $dir, $file );
		check_file( $filename, $verbose );
		push @files, $filename;
	}
}
if ($list)
{
	check_file( $list, $verbose );
	### make a resultdirectory
	my $resultdir = File::Spec->catfile( $dir, $resdir );
	_create_dir( $resultdir, $verbose, "xml merger resultdirectory" );
	my $filenamelist = File::Spec->catfile( $dir, $list );
	my @directories = read_file( $filenamelist, $verbose );
	foreach my $listdir (@directories)
	{
		my $filename = File::Spec->catfile( $dir, $listdir, "xquest.xml" );
		check_file( $filename, $verbose );
		push @files, $filename;
		unless ($nospecxml)
		{
			my $bnlist = basename($listdir);
			my $filename = File::Spec->catfile( $dir, $listdir, "$bnlist.spec.xml" );
			#check_file( $filename, $verbose );
			push @specfiles, $filename;
		}
	}
}else{
my $resultdir = File::Spec->catfile( $dir, $resdir );
_create_dir( $resultdir, $verbose, "xml merger resultdirectory" );
}
unless (@files)
{
	print "Error: no files found\n";
	exit;
}

#---------------------------------------------------------------------------
#  Parse all xml files into one hash, index with filname
#---------------------------------------------------------------------------
my %xmlhash;
foreach my $file (@files)
{
	($verbose) && print "Parsing file $file\n";
	my $tree = XML::TreeBuilder->new;    # empty tree
	$tree->parse_file($file);
	$xmlhash{$file} = $tree;
}

#print ($xmlhash{"xquest.xml"}->dump);
#---------------------------------------------------------------------------
#  Put the root elements attributes and all spectrum_search elements into one hash
#---------------------------------------------------------------------------
my %spectrumelements;
my %rootelements;
foreach my $file (@files)
{
	my $root = $xmlhash{$file}->find_by_tag_name('xquest_results');
	$rootelements{$file} = $root;
}

#---------------------------------------------------------------------------
# Save the merged XML
#---------------------------------------------------------------------------
#my $root = XML::Element->new('xquest_results');
my $xml_pi = XML::Element->new( '~pi', text => 'xml version="1.0"' );
## Create a new element
my $root = XML::Element->new('xquest_merger');
$root->attr( 'version', $version );
my $filesstring = join( ",", @files );
$root->attr( 'inputfiles', $filesstring );
foreach my $mergedfileelements ( keys %xmlhash )
{
## push the element into the xquest merger element
	$root->push_content( $xmlhash{$mergedfileelements} );
}
my $resultxml;
$resultxml .= $xml_pi->as_XML;
$resultxml .= $root->as_XML();
$outfilename = File::Spec->catfile( $dir, $resdir, $outfilename );
print $verbose && "Saving resultxml to file $outfilename\n";
save_to_file( $outfilename, $resultxml );
### save spectra
unless ($nospecxml)
{

	foreach my $specfile (@specfiles)
	{
		$outfilename = File::Spec->catfile( $dir, $resdir );
		$verbose && print "Saving Specfile $specfile to resultfolder $outfilename\n";
		unless (-e $specfile){
		warn "Warning: Cannot find spectrum xml file $specfile, won't copy file to resultdirectory.\n";	
		}else{
		copy( $specfile, $outfilename ) or die "Copy failed: $!";
	
		}
		
	}
}
print "\nMerged ", scalar(@files), " files succesfully!\n";

#---------------------------------------------------------------------------
#  Copy the xquest definition files and the database files into the mergedfolder
#---------------------------------------------------------------------------
$verbose && print "Copy definition files to $outfilename\n";
copy( $xquestdef, $outfilename ) or warn "Copy failed: $!";
copy( $xmmdef, $outfilename ) or warn "Copy failed: $!";

## read xquestdef
my $PARAMS = {};
open DEF, "<$xquestdef" or die "cannot open xquest definition file $xquestdef $!";
while ( my $line = <DEF> )
{
	if ( $line ne "" )
	{
		my @results = split( ' ', $line );
		$PARAMS->{ $results[0] } = $results[1];
	}
}
close(DEF);
##
my $targetdb = $PARAMS->{'database'};
my $decoydb  = $PARAMS->{'database_dc'};
unless ( -e $targetdb )
{
	$verbose && warn "Cannot find the target database: $targetdb\n";
} else
{
	$verbose && print "Copy target db to $outfilename\n";
	copy( $targetdb, $outfilename ) or warn "Copy failed: $!";
}
unless ( -e $decoydb )
{
	$verbose && warn "Cannot find the decoy database.\n";
} else
{
	$verbose && print "Copy decoy db to $outfilename\n";
	copy( $decoydb, $outfilename ) or warn "Copy failed: $!";
}

#---------------------------------------------------------------------------
#  Delete all XML tree objects, and HTML elements
#---------------------------------------------------------------------------
foreach my $file (@files)
{
	$xmlhash{$file}->delete;
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
	if ($verbose)
	{
		print "Checking file: " . $filename . "...";
	}

	# -e is for exists, -r for readable
	unless ( ( -e $filename ) && ( -r $filename ) )
	{
		print "Error: Cannot find/read the file $filename.\n";
		exit;
	} else
	{
		if ($verbose) { print " ok.\n" }
	}
	return;
}

#===  FUNCTION  ================================================================
#  NAME:        read_file
#  PURPOSE:     read_file
#  DESCRIPTION: Read file line by line and put results into an array
#  PARAMETERS:  $filename (including path if script is not executed in the folder)
#  RETURNS:     array
#===============================================================================
sub read_file
{
	my $filename = shift;
	my $verbose  = shift;
	my @array;
	if ($verbose) { print "Reading from file $filename\n" }
	open FILE, $filename or die $!;
	while ( my $line = <FILE> )
	{
		if ($verbose) { print "Reading line $line" }
		chomp($line);
		push( @array, $line );
	}
	($verbose) && print "\n";
	close FILE;
	return @array;
}

#===  FUNCTION  ================================================================
#  NAME:        save_to_file
#  PURPOSE:     save_to_file
#  DESCRIPTION: Save a string to a file
#  PARAMETERS:  $filename, $string, $append (1 if string should be appended)
#  RETURNS:     1 if done correctly
#===============================================================================
sub save_to_file
{
	my $filename = shift;
	my $text     = shift;
	my $append   = shift;
	if ($append)
	{
		open MYOUTFILE, ">>", "$filename"
		  or die "cannot open file $filename $!";
	} else
	{
		open MYOUTFILE, ">", "$filename"
		  or die "cannot open file $filename $!";
	}
	print MYOUTFILE $text;
	close MYOUTFILE;
	return 1;
}

sub _create_dir
{
	my $dirpath = shift;
	my $verbose = shift;
	my $msg     = shift;
	unless ( -e $dirpath )
	{
		$verbose && print "Create $msg directory: $dirpath\n";
		mkdir($dirpath);
	} else
	{
		$verbose && print "$msg: $dirpath exists\n";
	}
	return;
}

sub usage()
{
	print "
	SOFTWARE: ", basename($0), " version $version
	
	AUTHOR: Thomas Walzthoeni

	INFORMATION: A software/script to to merge xQuest XML files (merging of individual search results into one result xml file).
 	Creates a folder, with the merged.xml file, the spectrum xml files, the sequence databases and the xquest.def file. 
	These files are used by the xQuest/xProphet viewer to view and browse the results.
 	
 	USAGE: ", basename($0), " -Option [Parameter]

	REQUIRED OPTIONS [defaults]:
	-files \"filename1.xml filename2.xml\" or -list [resultdirectories_fullpath]
	
	Info: -list [resultdirectories_fullpath] a file with the path to the folders (one per line) of the searches where a xquest.xml must be present.
	This file is availaible in the rootdirectory if runXquest.pl was used to start the search.
	Using this option also tries to find the spectrum xml files of the corresponding search, which is neccessary to view the spectra in the viewer.
	If the -files option is used, filenames/paths of/to the xquest.xml files that should be merged must be provided (in double quotes separated by whitespace \" ... \").
 	Please note that the spectrum xml files must be copied manually to the resultfolder if this option is used.

	-out [merged_xquest.xml] outputfile of the result xml file. 	
 	-xquestdef [$xquestdef] name of xquest definition file
 	-xmmdef [$xmmdef] name of xmm definition file
 	-resdir [$resdir] resultdirectory name

 	OTHER OPTIONS:
 	-nospecxml don't copy specxml files
 	-v verbose, print progress information (recommended)
 	-help print this help
	
	EXAMPLE:
	",basename($0)," -v
 	
";
	exit;
}
