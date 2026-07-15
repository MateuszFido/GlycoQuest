#!/usr/bin/env perl
use strict;
#---------------------------------------------------------------------------
# indexspecfiles.pl
# A software/script to index xQuest spec.xml files
# Execute indexspecfiles.pl -help to display information and usage options.
# AUTHOR: Thomas Walzthoeni
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
use Getopt::Long;
use File::Basename;
use XML::TreeBuilder;
use Storable;

#---------------------------------------------------------------------------
#  Variables
#---------------------------------------------------------------------------
my ( $input, $verbose, $version, $force, $help);

#---------------------------------------------------------------------------
#  Default values
#---------------------------------------------------------------------------
$version     = "1.0";
$input       = "xquest.xml";
#---------------------------------------------------------------------------
#  Options that are passed by arguments
#---------------------------------------------------------------------------
GetOptions(
			'in=s'  => \$input,          ## xQuest input filename"
			'force' =>\$force,
			'v'     => \$verbose,
			'help'  => \$help,
);

&usage() if $help;

unless (-e $input){
print "Cannot find input file $input\n";
exit;
}
parseXML($input);



# Parse xQuest xml file
sub parseXML
{
	my $xmlfilename   = shift;
### Parse the XML file
	unless ( -e $xmlfilename )
	{
		die "Cannot find xquest result file $xmlfilename\n";
	}
	my $tree = XML::TreeBuilder->new();
	$tree->parse_file($xmlfilename);
	print "Parsing xquest XML file\n";
### Parsing all headers // one or more if a it is a merged result
	my @resultsheader = map { $_ } $tree->find('xquest_results');
## the running id for spectra start with 1 (0 is undef)
## Parsing all specrum search results
	foreach my $header (@resultsheader)
	{
		my $headerid = $header->attr('outputpath');

		#print "Indexing Header: $headerid";
		## Define the filename of the spectrumxml file
		my $specfn               = $headerid . ".spec.xml";
		my $specxmlfile = File::Spec->catfile( $specfn );
		if ( -e $specxmlfile )
		{
			## parse the spectra into a separate hash and write a hashfile
			parse_spec_xml_file($specxmlfile);
		}else{
			warn "Cannot find the spec.xml file $specxmlfile\n";
		}
	}
	$tree->delete();
}

#---------------------------------------------------------------------------
#  SUB parseXML for PARSING THE XML FILE
#---------------------------------------------------------------------------
sub parse_spec_xml_file
{
	my $xmlfilename = shift;
	my $outfilename = $xmlfilename . ".hash";
	if ( -e $outfilename && !$force)
	{
		print "Specxml file $xmlfilename already exists, use -force to reindex.\n";
		return;
	}
	my $tree = XML::TreeBuilder->new();
	$tree->parse_file($xmlfilename);
	my @resultsheader = map { $_ } $tree->find('xquest_spectra');
	### index the spectra by filename
	my $spechash = {};
	foreach my $header (@resultsheader)
	{

		#	print "Parsing Spectrum xml file\n";
		my @spectra = map { $_ } $header->find('spectrum');
		foreach my $spectrum (@spectra)
		{
			my $filename = $spectrum->attr('filename');

			#my $type = $spectrum->attr('type');
			my $content = $spectrum->content();

			#$spechash->{$filename}=$spectrum->as_XML;
			$spechash->{$filename} = $content->[0];
		}
	}
	store $spechash, $outfilename;
	print "Indexed $outfilename\n";
	$tree->delete();
}

sub usage()
{
	print "
	SOFTWARE: ", basename($0), " version $version
	
	AUTHOR: Thomas Walzthoeni

	INFORMATION: A software/script to index specxml files.
 	
 	USAGE: ", basename($0), " -Option [Parameter]

	REQUIRED OPTIONS [defaults]: 
	-in [xquest.xml] xquest.xml input filename
	-force force reindexing if spec.xml files are already indexed

	OTHER OPTIONS:
 	-help print this information
 	

	EXAMPLE: ", basename($0), "
 		
	";
	exit;
}
