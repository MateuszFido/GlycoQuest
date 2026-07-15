#!/usr/bin/env perl
use strict;

#---------------------------------------------------------------------------
# pseudoSH.pl
# A software/script to ...
# Execute pseudoSH.pl -help to display information and usage options.
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
use Getopt::Long;
use File::Spec;
use File::Basename;
my $directoryBase = dirname(__FILE__);
use FindBin;
#use lib "$FindBin::Bin/../../perl5";
use XML::TreeBuilder;

##########################################################
# Include modules dir as lib that is relative to the Script path
##########################################################
#use FindBin;
use lib "$FindBin::Bin/../modules";
##########################################################
use xQUEST_mzXMLscan;
my ( $verbose, $help, $input );
my $version = "1.0";
my $outfile = 'MASTER_RUN_PSEUDO.txt';
GetOptions(
			'mz=s'    => \$input,
			'out=s'   => \$outfile,
			'verbose' => \$verbose,
			'help'    => \$help,
);
&usage() unless ($input);
&usage() if $help;
open INPUT, "<$input" or die $!;

#---------------------------------------------------------------------------
#  Parse the mzXML File
#---------------------------------------------------------------------------
my $MS2scans = parse_mzXML($input);

#---------------------------------------------------------------------------
#  Print the pseudo MasterMap
#---------------------------------------------------------------------------
print_pseudo_MM( $outfile, $MS2scans );
close(INPUT);

sub print_pseudo_MM
{
	my $outfile     = shift;
	my $MS2scans    = shift;
	my $numscans    = scalar(@$MS2scans);
	my @sortedscans = sort { $a->precursor_mz <=> $b->precursor_mz } @$MS2scans;
	open OUT, ">$outfile" or die "cannot open $outfile for writing $!\n";
	print OUT 'M/z	Tr	CHRG	AC	SQ	Pep.Prob	# Rep. Match	AREA_0	scannumber', "\n";
	foreach my $scan (@sortedscans)
	{

		unless ( $scan->precursorIntensity )
		{

			#warn "Scan ",$scan->id," contains no intensity value, setting scan intensity to 1\n";
			#next;
		}
		print OUT join "\t", $scan->precursor_mz, $scan->Tr / 60, ( '+' . $scan->FT_charge ), undef, undef, undef, 1, $scan->precursorIntensity || 1, $scan->scannumber;
		print OUT "\n";
	}

	#M/z	Tr	CHRG	AC	SQ	Pep.Prob	# Rep. Match	AREA_0	scannumber
	#350.4928894	29.3296666666667	+3				1	45632	2547
}

sub parse_mzXML
{
	my $mzXMLfile = shift;
	## get the basename
	( my $mzXMLfilename = $mzXMLfile ) =~ s/.mzXML//;

	#print $mzXMLfilename;
	my @MS2scans = ();
	print "parsing $mzXMLfile...\n";

	#( my $basename = $mzXMLfilename ) =~ s/.mzXML//;
	my $mzXMLtree = XML::TreeBuilder->new();
	$mzXMLtree->parse_file($mzXMLfile);
	foreach my $scan ( $mzXMLtree->find('scan') )
	{
		my $scannr = $scan->attr('num');

		#my $scanid = join ".", $mzXMLfilename, $scannr;
		my $scanid = $scannr;

		#print $scannr."\n";
		if ( $scan->attr_get_i('msLevel') == 2 )
		{
			## check if there are peaks in this scan
			my $peaknum = $scan->attr_get_i('peaksCount');
			unless ( $peaknum == 0 )
			{
				push @MS2scans, XQUEST_mzXMLscan->new( $scan, $scanid, $mzXMLfilename, 1 );
			} else
			{
				print "No peaks found in scan number $scannr, skip scan\n";
				next;
			}
		}
	}
	$mzXMLtree->delete;
	print "read ", scalar(@MS2scans), " MS/MS scans\n";
	return \@MS2scans;
}

sub usage()
{
	print "
	SOFTWARE: ", basename($0), " version $version
	
	AUTHOR: Thomas Walzthoeni, based on a previous version by Oliver Rinner.

	INFORMATION: Program to extract features based on the information in the MS/MS scanheaders.

 	USAGE: ", basename($0), " -Option [Parameter]

	REQUIRED OPTIONS:
	-mz [] mzXML input filename
	-out [$outfile] output filename
	-v print control information

	OTHER OPTIONS: 
	-help print this help
	
	EXAMPLE:: $0 -mz YYY-XXXXXX.mzXML -out MASTER_RUN_test.txt
";
	exit;
	print "\nwrites a pseudo MasterMap based on MS/MS scanheaders only
		usage: ", basename($0), " [-file input] [-options]
        options:
                        -mz mzXML input file 
                        -out [$outfile] output file
                        -v print control information


\n";
	exit;
}
