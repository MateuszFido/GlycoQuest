#!/usr/bin/env perl
use strict;
use File::Basename;
use Getopt::Long;
use File::Spec;
use Cwd;
use Data::Dumper;
use FindBin;
#use lib "$FindBin::Bin/../../perl5";
use XML::TreeBuilder;
use XML::Element;

#---------------------------------------------------------------------------
# SOFTWARE: annotatexml.pl
# A software to annotate xquest.xml files
# Execute annotatexml.pl -help to display information and usage options
#---------------------------------------------------------------------------

#---------------------------------------------------------------------------
# This software is licensed under the Apache License, Version 2.0.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
# This software and associated documentation is provided "AS IS",
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND.
# See the License for the specific language governing permissions and
# limitations under the License.
#---------------------------------------------------------------------------

my ( $xmlfile, $annotationfile, $verbose, $outfile, $help, $annotatenative );
$xmlfile = "xquest.xml";
$outfile = "annotated_xquest.xml";
my $version = "1.0";
my %PARAMS  = (        ## the options that are passed as arguments ovverride the standard values if they are set
	## common settings
	'xmlfile'        => $xmlfile,
	'annotationfile' => $annotationfile,
	'native'         => $annotatenative,
	'out=s'          => $outfile,
	'v'              => $verbose,
	'help'           => $help,
);

GetOptions( \%PARAMS, 'xmlfile=s', 'annotationfile=s', 'v', 'nospec', 'help', 'out', 'native' );
my $PARAMS = \%PARAMS;
if ( $PARAMS->{'help'} )
{
	&usage;
}
if ( $PARAMS->{'v'} )
{
	print "Input Parameters:\n";
	print_params($PARAMS);
}
unless ( $PARAMS->{'xmlfile'} && ( $PARAMS->{'annotationfile'} || $PARAMS->{'native'} ) )
{
	print "Cannot find xquest xml file or annotationfile or -native options\n";
	&usage;
}
if ( $PARAMS->{'help'} )
{
	&usage;
}
check_file( $PARAMS->{'xmlfile'}, $PARAMS->{'v'} );
my $conditionshash = {};
my @annotationarray;
if ( $PARAMS->{'annotationfile'} )
{
	check_file( $PARAMS->{'annotationfile'}, $PARAMS->{'v'} );
### Read the Annotationfile
	@annotationarray = read_file( $PARAMS->{'annotationfile'}, $PARAMS->{'v'} );
### Create a hash with the condition as value and the annotation as key
	my $i = 0;
	foreach my $line (@annotationarray)
	{
		my @array = split( /\t/, $line );
		### generate the condition
		my $string = "$array[0]::$array[1]::$array[2]";
		$conditionshash->{ $array[0] } = $array[1];
		$i++;
	}
}
### parse the XML file
my $tree = XML::TreeBuilder->new();
$PARAMS->{'v'} && print "Parsing $PARAMS->{'xmlfile'} ... ";
$tree->parse_file( $PARAMS->{'xmlfile'} );
$PARAMS->{'v'} && print "done.\n";
### Delete the spectra if -nospec is selected
if ( $PARAMS->{'nospec'} )
{
	$PARAMS->{'v'} && print "-nospec option selected, will remove all spectrum elements\n";
	my @spectra = map { $_ } $tree->find('spectrum');
	foreach my $spec (@spectra)
	{
		$spec->delete_content();
	}
	print "Deleted ", scalar(@spectra), " spectrum Elements\n";
}

# find all xquest_results objects (if merged) in the xml $tree and save them in the @spectra array
my @xqresults            = map { $_ } $tree->find('xquest_results');
my $annotationcounter    = 0;
my $notannotationcounter = 0;
my $hitnum               = 0;
my $string;
my @search_hits;

foreach my $file (@xqresults)
{
	my $outputpath = $file->attr('outputpath');
	print "Outputpath is : " . $outputpath . "\n";
### Go through every hit of this resultfile
	@search_hits = map { $_ } $file->find('search_hit');
	foreach my $hit (@search_hits)
	{
### annotate native or with annotation array
		if ( $PARAMS->{'native'} )
		{
			my $string = $outputpath;
			unless ($string)
			{
				warn "No outputpath found! Hit not annotated!\n";
				$notannotationcounter++;
				next;
			}
			$hit->attr( 'annotation', $string );
			$annotationcounter++;
		}
### Go through all conditions
		if ( $PARAMS->{'annotationfile'} )
		{
			if ( $conditionshash->{$outputpath} )
			{

				#		print "outputpath is $outputpath\n";
				$annotationcounter++;
				$string = $conditionshash->{$outputpath};

				#		print "Annotate hit with $string\n";
				$hit->attr( 'annotation', $string );
			} else
			{
				$notannotationcounter++;
				print "Cannot annotate hit with outputpath $outputpath\n";
			}
		}
	}
}
print "Annotated $annotationcounter search hits, $notannotationcounter search hits could not be annotated\n";
#### WRITE THE XML FILE FILE
save_tree_to_file( $outfile, $tree );
$tree->delete;

#===  FUNCTION  ================================================================
#  NAME:        save_tree_to_file
#  PURPOSE:     save xml tree to a new file
#  DESCRIPTION: save xml tree to a new file
#  PARAMETERS:  filename, tree
#  RETURNS:     return
#===============================================================================
sub save_tree_to_file
{
	my $filename = shift;
	my $tree     = shift;
	print "Saving file $filename\n";
	open MYFILE, ">", $filename or die $!;
	print MYFILE $tree->as_XML;
	close(MYFILE);
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
	open FILE, $filename or die " Cannot open $filename $!";
	while ( my $line = <FILE> )
	{
		if ($verbose) { print "Reading line: $line" }
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

sub print_params
{
	my $hashref = shift;
	foreach my $key ( sort keys %$hashref )
	{
		my $value = $hashref->{$key};
		if ( ref($value) )
		{
			$value = $$value;
		}
		unless ($value) { $value = "not defined" }
		my $length = length($key);
		if ( $length > 6 )
		{
			print "$key \t  =>  $value\n";
		} else
		{
			print "$key \t\t  =>  $value\n";
		}
	}
}

sub usage()
{
	print "
	SOFTWARE: ", basename($0), " version $version
	
	AUTHOR: Thomas Walzthoeni

	INFORMATION: Program to annotate xquest.xml files 
 
 	USAGE: ", basename($0), " -Option [Parameter]

	REQUIRED OPTIONS: -annotationfile [filename] or -native
	
	-annotationfile [filename]: Annotates Search results with descriptions. 
	Use with a text file that contains: Filename and Description (one per line, tab-separated),
	Note that the Filenames have to be in the format: XXXX_matched (they must match the name of the search resultdirectory)
	
	-native: Annotate with the filenames used for the search (filenames are retrieved from the xquest xml file)
	
	OTHER OPTIONS [defaults]:
	-xmlfile [xquest.xml] input filename
	-out [annotated_xquest.xml] output filename
	-v verbose, prints program information
	
	";
	exit;
}
