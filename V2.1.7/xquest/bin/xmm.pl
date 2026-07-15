#!/usr/bin/env perl
use strict;

#---------------------------------------------------------------------------
# xmm.pl
# A software/script to search for isotopic scan pairs.
# Execute xmm.pl -help to display information and usage options.
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
##########################################################
# Include modules dir as lib that is relative to the Script path
##########################################################
use FindBin;
use lib "$FindBin::Bin/../modules/xmmlibs";
#use lib "$FindBin::Bin/../../perl5";
##########################################################

use Getopt::Long;
use ParseMap;
use XMM_Statistics;
use PARAMS;
my ( $annotated_output, %transhash, $verbose, $normby, $translation, $help, $mzXMLs, $stdout, $output, @mzXMLfiles, @reference_mzxmls, $delta, $stats, $profileoutputname, $proteingrouplist, $reference_mzxml, $nhighest_quant_out_handle, $referenceprotein, $quantoutput_peptides, $quantoutput,
	 $nhighest_quant_out_handle_peps, $pseudoSH );
my $version = "1.0";
## the standard values
my $deffile       = "xmm.def";
my $inclusionlist = "inclusionlist.xls";
my $mastermap     = "MASTER_RUN.txt";
GetOptions(
			'trans=s'      => \$translation,
			'master=s'     => \$mastermap,
			'out=s'        => \$output,
			'profileout=s' => \$profileoutputname,
			'stats=s'      => \$stats,
			'ref=s'        => \$reference_mzxml,
			'isolist=s'    => \$delta,
			'mzXMLs=s'     => \$mzXMLs,
			'def=s'        => \$deffile,
			'inc=s'        => \$inclusionlist,
			'group=s'      => \$proteingrouplist,
			'verbose'      => \$verbose,
			'help'         => \$help,
			'norm=s'       => \$referenceprotein,
			'stdout'       => \$stdout,
			'quantout=s'   => \$quantoutput,
			'pseudoSH'     => \$pseudoSH,
);
&usage() unless ( -e $mastermap );
&usage() if $help;
my $PARAMS = PARAMS->new( $deffile, $verbose );

if ($mzXMLs)
{
	@mzXMLfiles = glob($mzXMLs);
	print "mzXMLfiles: @mzXMLfiles\n";
}
if ($reference_mzxml)
{
	@reference_mzxmls = glob($reference_mzxml);
	print "reference mzXML files: @reference_mzxmls\n";
}
unless ($output)
{
	( $output = $mastermap ) =~ s/\.\w+/_annotated.xls/;
}
unless ($stats)
{
	( $stats = $mastermap ) =~ s/\.\w+/_stats.xls/;
}
open INPUT,  "<$mastermap" or die $!;
open OUTPUT, ">$output"    or die $!;
open STATS,  ">$stats"     or die $!;
if ( $PARAMS->single_quant )
{
	unless ($quantoutput)
	{
		( $quantoutput          = $mastermap ) =~ s/\.\w+/_quant.xls/;
		( $quantoutput_peptides = $mastermap ) =~ s/\.\w+/_quant_peptides.xls/;
	}
	open QUANTOUTPUT, ">$quantoutput" or die $!;
	$nhighest_quant_out_handle = *QUANTOUTPUT;
	open QUANTOUTPUT_PEPS, ">$quantoutput_peptides" or die $!;
	$nhighest_quant_out_handle_peps = *QUANTOUTPUT_PEPS;
}
my $statsout     = *STATS;
my $resultoutput = *OUTPUT;

#my $delta_ETD_shiftout = *STDOUT;
if ($stdout)
{
	$resultoutput = *STDOUT;
}

#---------------------------------------------------------------------------
#  Parsing a MasterMap File
#---------------------------------------------------------------------------
print "\n\nparsing MasterMap $mastermap\n";
### Parses all lines of the MasterMap and creates Feature objects
my $MASTERMAP = ParseMap->new( $mastermap, $PARAMS );
if ( $PARAMS->normalize_MasterMap )
{
	$MASTERMAP->normalize_MasterMap;
}
if ( $PARAMS->profiler )
{
	my $profileoutfilehandle = my $profilestatsoutfilehandle = *STDOUT;
	if ($profileoutputname)
	{
		open PROFILEOUTPUT, ">$profileoutputname" . '_profilestats.xls' or die $!;
		$profilestatsoutfilehandle = *PROFILEOUTPUT;
		open PROFILES, ">$profileoutputname" . '_profiles.xls' or die $!;
		$profileoutfilehandle = *PROFILES;
	}
	if ($proteingrouplist)
	{
		my $grouphash = get_grouphash($proteingrouplist);
		$MASTERMAP->separate_protgroups_by_group_hash($grouphash);

		#  $MASTERMAP->quantify_protgroups_by_nhighest;
	} else
	{
		$MASTERMAP->separate_protgroups;
	}
	$MASTERMAP->calcprofiles;
	if ( $PARAMS->target )
	{
		$MASTERMAP->ScoreProtgroups;
	}
	if ($referenceprotein)
	{
		$MASTERMAP->normalize_profiles_to_target($referenceprotein);
	} elsif ( $PARAMS->levelto )
	{
		$MASTERMAP->normalize_profiles_to_target( $PARAMS->levelto );
	}
	$MASTERMAP->print_profile($profileoutfilehandle);
	$MASTERMAP->printprofilefile($profilestatsoutfilehandle);
}
if ( $PARAMS->single_quant )
{
	my $grouphash     = undef;
	my $proteingroups = undef;
	if ($proteingrouplist)
	{
		$grouphash     = get_grouphash($proteingrouplist);
		$proteingroups = $MASTERMAP->separate_protgroups_by_group_hash($grouphash);
		$MASTERMAP->quantify_protgroups_by_nhighest;
	} else
	{
		warn "a protein table is needed\n";
	}
	$MASTERMAP->print_nhighest_proteinquantifications( $nhighest_quant_out_handle, $nhighest_quant_out_handle_peps );
}

#---------------------------------------------------------------------------
#  Parse the mzXMLs
#---------------------------------------------------------------------------
if ( $mzXMLs && $pseudoSH )
{
	### Assign scannumbers from pseudosh
	$MASTERMAP->parse_pseudo_mzXML( \@mzXMLfiles, 0 );
} elsif ($mzXMLs)
{
	$MASTERMAP->parse_mzXML( \@mzXMLfiles, 1 );
	$MASTERMAP->assignscans;
}

#---------------------------------------------------------------------------
#  SEARCH FOR ISOTOPIC PAIRS
#---------------------------------------------------------------------------
if ( $PARAMS->deltashift )
{
	print "File to store deltashift: " . $delta . "\n";
	my $deltashiftout     = *STDOUT;
	my $deltafeatures_xls = *STDOUT;
	if ( $delta && !$stdout )
	{
		open DELTAOUT, ">$delta" or die $!;
		$deltashiftout = *DELTAOUT;
		open DELTAXLS, ">$delta" . "_isopairs.xls" or die $!;
		$deltafeatures_xls = *DELTAXLS;
	}
	if ($stdout)
	{
		$statsout = *STDOUT;
	}
	if ( $PARAMS->tripleshift )
	{
		print "Tripleshift is activated\n";
		$MASTERMAP->searchtripleshifts;
	}
	### searches for deltashifts and makes new isopair objects
	### isopair object consists 2 feature objects
	$MASTERMAP->searchdeltashifts;

	# exit;
	$MASTERMAP->delta_feature_stats( $PARAMS->featurestatIDs, $statsout, $deltafeatures_xls );
	if ($reference_mzxml)
	{
		open MAPPING_STATS, ">mapping.xls";
		$MASTERMAP->parse_ETD_mzXML( \@reference_mzxmls, 1 );
		$MASTERMAP->transferscans_from_reference( $PARAMS, *MAPPING_STATS );
		close(MAPPING_STATS);
	}
### here the scanpairs are printed, scanpairs are isopairs objects
	if ( $PARAMS->highresms2 )
	{
		print "Highres MS2 option selected\n";
	}
	if ( $PARAMS->printisotopicscanpairs )
	{
		my $twoscans = $MASTERMAP->get_isopairs_with_two_scans;
		foreach my $scanpair (@$twoscans)
		{
			if ( $PARAMS->highresms2 )
			{
				$scanpair->print_isopair_into_spectrumlist_highres($deltashiftout);
			} else
			{
				$scanpair->print_isopair_into_spectrumlist($deltashiftout);
			}
		}
		if ($reference_mzxml)
		{
			if ( $delta && !$stdout )
			{
				open DELTA_ETD_OUT_MGF, ">$delta" . "_etd_cid.mgf" or die $!;
				open DELTA_ETD_OUT_ISOPAIRS, ">$delta" . "_isopairs_etd_cid.txt"
				  or die $!;

				#		$delta_ETD_shiftout = *DELTA_ETD_OUT;
			}
			my $etdscanpairs = $MASTERMAP->get_isopairs_with_etd_cid_pair;
			foreach my $scanpair (@$etdscanpairs)
			{
				$scanpair->print_etd_isopair_into_spectrumlist( *DELTA_ETD_OUT_MGF, *DELTA_ETD_OUT_ISOPAIRS );
			}
			close(DELTA_ETD_OUT_MGF);
			close(DELTA_ETD_OUT_ISOPAIRS);
		}
	}
	if ( $PARAMS->printlightonlypairs )
	{
		my $lightscanonly = $MASTERMAP->get_isopairs_with_light_scan;
		foreach my $scanpair (@$lightscanonly)
		{
			if ( $PARAMS->highresms2 )
			{
				$scanpair->print_lightonly_isopair_into_spectrumlist_highres($deltashiftout);
			} else
			{
				$scanpair->print_lightonly_isopair_into_spectrumlist($deltashiftout);
			}

			#$scanpair->print_lightonly_isopair_into_spectrumlist($deltashiftout);
		}
		if ($reference_mzxml)
		{
			my $etdscanpairs = $MASTERMAP->get_isopairs_with_light_etd_cid_pair;
			foreach my $scanpair (@$etdscanpairs)
			{
				$scanpair->print_lightonly_etd_isopair_into_spectrumlist($deltashiftout);
			}
		}
	}
	if ( $PARAMS->printheavyonlypairs )
	{
		my $heavyscanonly = $MASTERMAP->get_isopairs_with_heavy_scan;
		foreach my $scanpair (@$heavyscanonly)
		{
			$scanpair->print_heavyonly_isopair_into_spectrumlist($deltashiftout);
		}
		if ($reference_mzxml)
		{
			my $etdscanpairs = $MASTERMAP->get_isopairs_with_heavy_etd_cid_pair;
			foreach my $scanpair (@$etdscanpairs)
			{
				$scanpair->print_heavyonly_etd_isopair_into_spectrumlist($deltashiftout);
			}
		}
	}
}
if ( $PARAMS->inclusion_list )
{
	$MASTERMAP->makeinclusionlist($inclusionlist);
}
$MASTERMAP->MapStats;
$MASTERMAP->printmastermap($resultoutput);
close(OUTPUT);
close(INPUT);

sub get_grouphash
{
	my %proteinhash = ();
	open LIST, "<$proteingrouplist" or die $!;
	while (<LIST>)
	{
		chomp;
		my ( $protein, $peps ) = split;
		my @peps = split /\+/, $peps;
		foreach my $pep (@peps)
		{
			$proteinhash{$pep} = $protein;
		}
	}
	close LIST;
	return \%proteinhash;
}

sub usage()
{
	print "
	SOFTWARE: ", basename($0), " version $version
	
	AUTHOR: Thomas Walzthoeni based on orginal work by Oliver Rinner.

	INFORMATION: A software/script to search for isotopic scan pairs.
 	
 	USAGE: ", basename($0), " -Option [Parameter]

	REQUIRED OPTIONS [defaults]:

  	-def [xmm.def] xmm.def file
	-mz [] mzXML filename, provided as \"mzXMLfilename.mzXML\"
	-iso [] output filename

	OTHER OPTIONS:	
	-pseudoSH, pseudoSH Master map is used, features are not assigned to scans (much faster), 
	scannumber is obtained from the master map, must be in combination with -mz to indicate mzxml

	EXAMPLE
	xmm.pl -master ./MASTER_RUN/MASTER_RUN.txt -iso mzXMLfilename_matched.txt -def xmm.def  -pseudoSH -mz \"mzXMLfilename.mzXML\"
	
";
	exit;
}
