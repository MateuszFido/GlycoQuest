#!/usr/bin/env perl
use strict;

#---------------------------------------------------------------------------
# xprophet.pl
# A software/script to estimate FDRs for cross-linked peptides identified by MS.
# Execute xprophet.pl -help to display information and usage options.
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
##########################################################
# Include modules dir as lib that is relative to the Script path
##########################################################
use FindBin;
use lib "$FindBin::Bin/../modules";
##########################################################
use Environment;
use Read_Params;
use Data::Dumper;
use Storable;
use Cwd;
use File::Path;
#use lib "$FindBin::Bin/../../perl5";
use XML::TreeBuilder;
use Getopt::Long;
use File::Basename;

#---------------------------------------------------------------------------
#  Variables
#---------------------------------------------------------------------------
my ( $input, $outfilename, $deffile, $help, $verbose, $version, $optifdr, $getdef, $optifile, $forceindex, $minionsmatched );

#---------------------------------------------------------------------------
#  Default values
#---------------------------------------------------------------------------
$version     = "2.5.5";
$input       = "xquest.xml";
$outfilename = "xproph_out.xml";
$deffile     = "xproph.def";
my $author      = "Thomas Walzthoeni";
my $affiliation = "ETH Zurich, Institute of Molecular Systems Biology, Wolfgang Pauli-Str. 16, CH-8093 Zurich";
my $mailto      = "walzthoeni\@imsb.biol.ethz.ch";

#---------------------------------------------------------------------------
#  Options that are passed by arguments
#---------------------------------------------------------------------------
GetOptions(
			'in=s'             => \$input,            ## input filename"
			'out=s'            => \$outfilename,      ## output filename
			'def=s'            => \$deffile,          ## deffilename
			'help'             => \$help,
			'v'                => \$verbose,
			'opti=s'           => \$optifdr,
			'optifile=s'       => \$optifile,
			'getdef'           => \$getdef,
			'forceindex'       => \$forceindex,
			'minionsmatched=s' => \$minionsmatched,
);

#---------------------------------------------------------------------------
#  Generate the paths
#---------------------------------------------------------------------------
my $dir         = getcwd;
my $xmlfilename = File::Spec->catfile( $dir, $input );
my $xmlbasename = basename( $input, ".xml" );
$deffile = File::Spec->catfile( $dir, $deffile );
my $resultshashfilename = File::Spec->catfile( $dir, "$xmlbasename.resulthash.db" );
&usage() if $help;
print "xProhpet version $version, by $author.\n";

#---------------------------------------------------------------------------
#  Check if a xproph.def file exists, otherwise generate it and exit
#---------------------------------------------------------------------------
if ( !-e $deffile || $getdef )
{
### Generate the Deffile
	print "No definition file found, will generate a definition file.\n" unless ($getdef);
	open DEFFILE, ">", $deffile or die "cannot open file $deffile $!";
	my $fh = *DEFFILE;
	print $fh "### xProphet v.$version definition file ###\n";
	print $fh "### Parameters for availaible filters, which are applied before calculating FDRs. ###\n";
	print $fh "minborder -5 # filter for minimum precursor mass error (ppm)\n";
	print $fh "maxborder 5  # filter for maximum precursor mass error (ppm)\n";
	print $fh "mindeltas 0.95  # filter for delta score, 0 is no filter, minimum delta score required, hits are rejected if larger or equal\n";
	print $fh "minionsmatched 0 # Filter for minimum matched ions per peptide\n";
	print $fh "uniquexl 1 # [1|0] calculate statistics based on unique IDs [1]\n";
	print $fh "qtransform 1 # transform simple FDR to q-FDR values\n";
	print $fh "### Internal parameters\n";
	print $fh "nidx 0 # reparse xml file\n";
	print $fh "minscore 0 # minimum ld-score to be considered\n";
	close(DEFFILE);
	print "Wrote definition file $deffile, please edit and rerun xProphet.\n";
	exit;
}

#---------------------------------------------------------------------------
#  Read the definition file
#---------------------------------------------------------------------------
my $PARAMS;
### Default parameters
$PARAMS->{'nranks'}         = 1;      ### Leave a 1, more ranks are not supported anymore
$PARAMS->{'minborder'}      = -5;
$PARAMS->{'maxborder'}      = 5;
$PARAMS->{'mindeltas'}      = 0.95;
$PARAMS->{'minionsmatched'} = 0;
$PARAMS->{'uniquexl'}       = 1;
$PARAMS->{'qtransform'}     = 1;
$PARAMS->{'minscore'}       = 15;

### Read param file and overwrite the default params
read_params( $deffile, $PARAMS );

### Overwrite by definition from input params
if ($forceindex)
{
	$PARAMS->{'nidx'} = 1;
}

if ($minionsmatched)
{
	$PARAMS->{'minionsmatched'} = $minionsmatched;
}

#---------------------------------------------------------------------------
#  Print params
#---------------------------------------------------------------------------
print "Filtering of precursor mass error from ", $PARAMS->{'minborder'}, " to ", $PARAMS->{'maxborder'}, " ppm is used.\n";
print "Filtering of hits by a deltascore of ", $PARAMS->{'mindeltas'}, " is used.\n";
print "Filtering of hits by minimum ions matched: ", $PARAMS->{'minionsmatched'}, " is used.\n" if ($PARAMS->{'minionsmatched'}>0);
print "Filtering of hits by minimum score of: ", $PARAMS->{'minscore'}, " is used.\n" if ($PARAMS->{'minscore'}>0);

if ( $PARAMS->{'uniquexl'} )
{
	print "Error model is generated based on unique cross-links.\n";
} else
{
	print "Error model is generated based on redundant cross-links.\n";
}

#---------------------------------------------------------------------------
# Parse xQuest XML file
#---------------------------------------------------------------------------
my $resultshash = {};
my $summaryhash = {};
my $nidx        = $PARAMS->{'nidx'};
if ( -e $resultshashfilename && !$nidx )
{
	print "NOTE: Retrieving hits from file $resultshashfilename\n";
	print "NOTE: If you changed your xquest.xml reparse the xquest.xml file using the nidx option in the param.def file\n";
	$resultshash = retrieve($resultshashfilename);
} else
{
	print "Parsing xQuest XML file $xmlfilename\n";
	parseXML( $xmlfilename, $resultshash, $summaryhash );
	store $resultshash, $resultshashfilename;
}

# Sorting of hits into classes and filtering of the hits (unique Ids, Mass error, deltaS, minionsmatched, minscore)
my ( $intralinks, $intradecoys, $interlinks, $interdecoys, $decoys, $targethits, $monolinks, $monodecoys, $looplinks, $looplinkdecoys, $hybriddecoysinterlinks, $fulldecoysinterlinks, $fulldecoysintralinks, $hybriddecoysintralinks, $intralinksIds, $spectrumidshash ) =
  _generate_result_hashs( $resultshash, $PARAMS );

#---------------------------------------------------------------------------
#  Generate Grouped Hashes // @groupsarray is used for the grouping
#---------------------------------------------------------------------------
my @groupsarray;
### The array for the score bins
### 0.1 steps from 0 to 100
### 0.1 steps are necassary since precision is 10
### Groupshash is 0, 0.1, 0.1 .. 100
foreach my $i ( 0 .. 1000 )
{
	push @groupsarray, $i * 0.1;
}

#generate grouped hashes with score and the number of hits above the corresponding score
my $intralinksgrouped             = _gen_grouped_hash( \@groupsarray, $intralinks,             $PARAMS );
my $intralinkdecoysgrouped        = _gen_grouped_hash( \@groupsarray, $intradecoys,            $PARAMS );
my $hybriddecoysintralinksgrouped = _gen_grouped_hash( \@groupsarray, $hybriddecoysintralinks, $PARAMS );
my $fulldecoysintralinksgrouped   = _gen_grouped_hash( \@groupsarray, $fulldecoysintralinks,   $PARAMS );
my $interlinksgrouped             = _gen_grouped_hash( \@groupsarray, $interlinks,             $PARAMS );
my $interlinkdecoysgrouped        = _gen_grouped_hash( \@groupsarray, $interdecoys,            $PARAMS );
my $hybriddecoysinterlinksgrouped = _gen_grouped_hash( \@groupsarray, $hybriddecoysinterlinks, $PARAMS );
my $fulldecoysinterlinksgrouped   = _gen_grouped_hash( \@groupsarray, $fulldecoysinterlinks,   $PARAMS );
my $mono_intra_grouped            = _gen_grouped_hash( \@groupsarray, $monolinks,              $PARAMS );
my $mono_intra_decoys_grouped     = _gen_grouped_hash( \@groupsarray, $monodecoys,             $PARAMS );
my $targethitsgrouped             = _gen_grouped_hash( \@groupsarray, $targethits,             $PARAMS );
my $decoyhitsgrouped              = _gen_grouped_hash( \@groupsarray, $decoys,                 $PARAMS );
## Generate a hash of all FP hits
## Is then directly used in the fdr estimation
my $fphitsgrouped = _genfphash( $intralinkdecoysgrouped, $fulldecoysintralinksgrouped, $interlinkdecoysgrouped, $fulldecoysinterlinksgrouped, $mono_intra_decoys_grouped );
### Generate a grouparray for the histograms with binsize 1
my @groupsarray2;

foreach my $i ( 0 .. 100 )
{
	push @groupsarray2, $i;
}
## Generate histograms hashes, with score and number of hits within the score bins
my $interhisto                   = _gen_histo_hash( \@groupsarray2, $interlinks,             $PARAMS );
my $interdchisto                 = _gen_histo_hash( \@groupsarray2, $interdecoys,            $PARAMS );
my $intrahisto                   = _gen_histo_hash( \@groupsarray2, $intralinks,             $PARAMS );
my $intradchisto                 = _gen_histo_hash( \@groupsarray2, $intradecoys,            $PARAMS );
my $mono_intra_histo             = _gen_histo_hash( \@groupsarray2, $monolinks,              $PARAMS );
my $mono_intra_decoys_histo      = _gen_histo_hash( \@groupsarray2, $monodecoys,             $PARAMS );
my $hybriddecoysinterlinks_histo = _gen_histo_hash( \@groupsarray2, $hybriddecoysinterlinks, $PARAMS );
my $fulldecoysinterlinks_histo   = _gen_histo_hash( \@groupsarray2, $fulldecoysinterlinks,   $PARAMS );
my $hybriddecoysintralinks_histo = _gen_histo_hash( \@groupsarray2, $hybriddecoysintralinks, $PARAMS );
my $fulldecoysintralinks_histo   = _gen_histo_hash( \@groupsarray2, $fulldecoysintralinks,   $PARAMS );

# generate cumulative grouped hashes with score and the number of hits above the corresponding score
my $intralinksgrouped2             = _gen_grouped_hash( \@groupsarray2, $intralinks,             $PARAMS );
my $intralinkdecoysgrouped2        = _gen_grouped_hash( \@groupsarray2, $intradecoys,            $PARAMS );
my $hybriddecoysintralinksgrouped2 = _gen_grouped_hash( \@groupsarray2, $hybriddecoysintralinks, $PARAMS );
my $fulldecoysintralinksgrouped2   = _gen_grouped_hash( \@groupsarray2, $fulldecoysintralinks,   $PARAMS );
my $interlinksgrouped2             = _gen_grouped_hash( \@groupsarray2, $interlinks,             $PARAMS );
my $interlinkdecoysgrouped2        = _gen_grouped_hash( \@groupsarray2, $interdecoys,            $PARAMS );
my $hybriddecoysinterlinksgrouped2 = _gen_grouped_hash( \@groupsarray2, $hybriddecoysinterlinks, $PARAMS );
my $fulldecoysinterlinksgrouped2   = _gen_grouped_hash( \@groupsarray2, $fulldecoysinterlinks,   $PARAMS );
my $mono_intra_grouped2            = _gen_grouped_hash( \@groupsarray2, $monolinks,              $PARAMS );
my $mono_intra_decoys_grouped2     = _gen_grouped_hash( \@groupsarray2, $monodecoys,             $PARAMS );
my $targethitsgrouped2             = _gen_grouped_hash( \@groupsarray2, $targethits,             $PARAMS );
my $decoyhitsgrouped2              = _gen_grouped_hash( \@groupsarray2, $decoys,                 $PARAMS );
##
my $fphitsgrouped2 = _genfphash( $intralinkdecoysgrouped2, $fulldecoysintralinksgrouped2, $interlinkdecoysgrouped2, $fulldecoysinterlinksgrouped2, $mono_intra_decoys_grouped2 );

#---------------------------------------------------------------------------
#  Create a statistics file
#---------------------------------------------------------------------------
my $resultfile = File::Spec->catfile( $dir, "xprophstat.xls" );
print "Writing statistics into $resultfile\n";
open RESULTS, ">", $resultfile or die "cannot open file $resultfile $!";
my $fh = *RESULTS;
### Intralink cumulative statistics
store_histogram( $fh, $intralinksgrouped2,           $PARAMS, "CUM STAT target Intra-protein cross links" );
store_histogram( $fh, $intralinkdecoysgrouped2,      $PARAMS, "CUM STAT decoy Intra-protein cross links" );
store_histogram( $fh, $fulldecoysintralinksgrouped2, $PARAMS, "CUM STAT full decoy Intra-protein cross links" );
### Interlink cumulative statistics
store_histogram( $fh, $interlinksgrouped2,           $PARAMS, "CUM STAT target Inter-protein cross links" );
store_histogram( $fh, $interlinkdecoysgrouped2,      $PARAMS, "CUM STAT decoy Inter-protein cross links" );
store_histogram( $fh, $fulldecoysinterlinksgrouped2, $PARAMS, "CUM STAT full decoy Inter-protein cross links" );
### Mono and Intralinks
store_histogram( $fh, $mono_intra_grouped2,        $PARAMS, "CUM STAT target mono- and intra-links" );
store_histogram( $fh, $mono_intra_decoys_grouped2, $PARAMS, "CUM STAT decoy mono- and intra-links" );
### Total EST FPS
store_histogram( $fh, $fphitsgrouped2,               $PARAMS, "TOTAL EST FPS" );
store_histogram( $fh, $interhisto,                   $PARAMS, "HISTOGRAM Inter-protein cross links" );
store_histogram( $fh, $interdchisto,                 $PARAMS, "HISTOGRAM Inter-protein decoy cross links" );
store_histogram( $fh, $hybriddecoysinterlinks_histo, $PARAMS, "HISTOGRAM Hybrid Inter-protein decoy cross links" );
store_histogram( $fh, $fulldecoysinterlinks_histo,   $PARAMS, "HISTOGRAM Full Inter-protein decoy cross links" );
store_histogram( $fh, $intrahisto,                   $PARAMS, "HISTOGRAM Intra-protein cross-links" );
store_histogram( $fh, $intradchisto,                 $PARAMS, "HISTOGRAM Intra-protein decoy cross-links" );
store_histogram( $fh, $hybriddecoysintralinks_histo, $PARAMS, "HISTOGRAM Hybrid Intra-protein decoy cross-links" );
store_histogram( $fh, $fulldecoysintralinks_histo,   $PARAMS, "HISTOGRAM Full Intra-protein decoy cross-links" );
store_histogram( $fh, $mono_intra_histo,             $PARAMS, "HISTOGRAM Mono and Looplinks" );
store_histogram( $fh, $mono_intra_decoys_histo,      $PARAMS, "HISTOGRAM decoy Mono and Looplinks" );
close(RESULTS);

#---------------------------------------------------------------------------
#  Generate the hash with FDRs
#---------------------------------------------------------------------------
print_hash();
my $intrafdr = calc_local_fdr_xlinks( $intralinksgrouped, $intralinkdecoysgrouped, $fulldecoysintralinksgrouped, $targethitsgrouped, $fphitsgrouped );
my $interfdr = calc_local_fdr_xlinks( $interlinksgrouped, $interlinkdecoysgrouped, $fulldecoysinterlinksgrouped, $targethitsgrouped, $fphitsgrouped );
## Calc FDR for mono and looplinks                 = 1;
my $monolooplinkfdr = calc_local_fdr_monoandlooplinks( $mono_intra_grouped, $mono_intra_decoys_grouped, $targethitsgrouped, $fphitsgrouped );

#---------------------------------------------------------------------------
#  qvalue transform
#---------------------------------------------------------------------------
if ( $PARAMS->{'qtransform'} )
{
	my $intrafdr2 = calcqvalue($intrafdr);
	print "Transformed FDR to q-FDR for intra-protein xlinks\n";
	my $interfdr2 = calcqvalue($interfdr);
	print "Transformed FDR to q-FDR for inter-protein xlinks\n";
	my $monolooplinkfdr2 = calcqvalue($monolooplinkfdr);
	print "Transformed FDR to q-FDR for mono- and loop-protein xlinks\n";
## replace hashes
	$intrafdr        = $intrafdr2;
	$interfdr        = $interfdr2;
	$monolooplinkfdr = $monolooplinkfdr2;
}
## Optimization run
## Count the number of the different links at FDR 0 to 0.1
if ($optifdr)
{
	print "Optimization run with FDR cutoff: $optifdr\n";
	open RES, ">>", $optifile or die "cannot open file $optifile $!";
	my $intralinks_at_fdr = count_hits_for_fdr( $intralinks, $intrafdr, $optifdr );
	print "Intra-protein cross-links @ FDR<= $optifdr : $intralinks_at_fdr\n";
	my $interlinks_at_fdr = count_hits_for_fdr( $interlinks, $interfdr, $optifdr );
	print "Inter-protein cross-links @ FDR<= $optifdr : $interlinks_at_fdr\n";
	my $mono_looplink_at_fdr = count_hits_for_fdr( $monolinks, $monolooplinkfdr, $optifdr );
	print "Mono and Looplinks cross-links @ FDR<= $optifdr : $mono_looplink_at_fdr\n";
## Prepare line for results
	my $mions  = $PARAMS->{'minionsmatched'};
	my $result = "$mions \t $optifdr \t $intralinks_at_fdr \t $interlinks_at_fdr \t $mono_looplink_at_fdr\n";
	print RES $result;
	close(RES);
## Exit here
	exit 0;
}

#---------------------------------------------------------------------------
#  Lookup the FDRs and write to an xml
#---------------------------------------------------------------------------
my $tree = parseXMLandAddFdr( $xmlfilename, $intrafdr, $interfdr, $monolooplinkfdr, $spectrumidshash );

#---------------------------------------------------------------------------
#  Save tree to file
#---------------------------------------------------------------------------
save_tree_to_file( $outfilename, $tree );

#
sub count_hits_for_fdr
{
	my $resulthash = shift;    # filtered resultshash, use the initial hashes
	my $fdrhash    = shift;    # corresponding fdr hash
	my $fdrcutoff  = shift;
## Check input
	unless ( $resulthash || $fdrhash || $fdrcutoff )
	{
		die "ERROR count_hits: Param missing\n";
	}
	my $count = 0;
## 	resultshashes are only a hash with the rank as key (only rank 1 avail) and score as arrayref
	my $resultsarray = $resulthash->{1};

	#print Dumper ($resultsarray);
## If there are no hits then the array is empty
	unless ( defined($resultsarray) )
	{
		return 0;
	}

	# Sort descending
	my @resultsarraysorted = sort { $b <=> $a } @$resultsarray;
	foreach my $score (@resultsarraysorted)
	{
		## transform to int with precision 10
		my $intscore = int( $score * 10 );
		unless ( defined( $fdrhash->{$intscore} ) )
		{
			warn " WARN: FDR not defined for $score $intscore!\n";
			next;
		}
		my $fdr = $fdrhash->{$intscore};
		## print "Score: $score FDR: $fdr\n";
		## count if lower than cutoff
		if ( $fdr <= $fdrcutoff )
		{
			$count++;
		}
	}

	# print "$count hits <= FDR cutoff $fdrcutoff \n";
	return $count;
}

sub print_hash
{
	my $hashref = shift;
	my @keys = sort { $a <=> $b } keys %{$hashref};
	foreach my $key (@keys)
	{
		print "Key:" . $key . "  Value: " . $hashref->{$key} . "\n";
	}
}

sub read_params
{
	my $deffile = shift;
	my $PARAMS  = shift;
	open DEFFILE, "<$deffile" or die "cannot open table $deffile $!";
	while ( my $line = <DEFFILE> )
	{
		chomp($line);
		if ($line)
		{
			my @results = split( " ", $line );
			$PARAMS->{ $results[0] } = $results[1];
		}
	}
	return $PARAMS;
}

sub calcqvalue
{
	my $fdrhash    = shift;
	my $qvaluehash = shift;
## sort by keys desc
	my @keys = sort { $b <=> $a } keys %$fdrhash;
	for ( my $i = 0 ; $i < @keys ; $i++ )
	{
		if ( $fdrhash->{ $keys[$i] } eq "n/a" )
		{
			next;
		}
		my $currentfdr = $fdrhash->{ $keys[$i] };

		#print "KEY:".$keys[$i]." FDR: ".$fdrhash->{$keys[$i]}."\n";
		my $smallestfdr = $currentfdr;
		## check if there is a smaller fdr below, start with the current value
		for ( my $y = $i ; $y < @keys ; $y++ )
		{
			my $fdrtocheck = $fdrhash->{ $keys[$y] };
			if ( $fdrtocheck < $smallestfdr )
			{
				$smallestfdr = $fdrtocheck;
			}
		}
		## if smallestfdr is smaller than currentfdr thant correct
		if ( $smallestfdr < $currentfdr )
		{

			#print "Correct FDR: $currentfdr to $smallestfdr\n";
			## replace value with min fdr
			$qvaluehash->{ $keys[$i] } = $smallestfdr;
		} else
		{
			## keep currentfdr because it is the smallest
			$qvaluehash->{ $keys[$i] } = $currentfdr;
		}
	}
	return $qvaluehash;
}

sub _genfphash
{
	my $intradecoys          = shift->{1};
	my $fulldecoysintralinks = shift->{1};
	my $interdecoys          = shift->{1};
	my $fulldecoysinterlinks = shift->{1};
	my $mono_intradecoys     = shift->{1};
	my @keys                 = sort { $a <=> $b } keys %{$intradecoys};
	my $fpgroupedhash        = {};
	foreach my $key (@keys)
	{
		my $intrafp = $intradecoys->{$key} - 2 * $fulldecoysintralinks->{$key};
		my $interfp = $interdecoys->{$key} - 2 * $fulldecoysinterlinks->{$key};
		my $monofp  = $mono_intradecoys->{$key};
		$fpgroupedhash->{1}->{$key} = $intrafp + $interfp + $monofp;
	}
	return $fpgroupedhash;
}

sub calc_local_fdr_monoandlooplinks
{
### get the hashes for rank 1
	my $targethitshash  = shift->{1};
	my $decoyhitshash   = shift->{1};
	my $targethitstotal = shift->{1};
	my $fphitstotal     = shift->{1};
	my $fdr;
	my $resultshash = {};
## get the binkeys (Bins of scores!)
	my @keys = sort { $a <=> $b } keys %$targethitshash;
	for ( my $i = 0 ; $i < @keys ; $i++ )
	{

		#print "Key is: ".$keys[$i];
		## check if there are any target hits in this bin
		unless ( ( $targethitshash->{ $keys[$i] } ) )
		{
			$fdr = "n/a";
		} else
		{
			### Calculation of the FDR
			my $t1;
			unless ( $fphitstotal->{ $keys[$i] } )
			{
				$t1 = 0;
			} else
			{
				### Est FP hits
				my $estfphits = $decoyhitshash->{ $keys[$i] };
				$t1 = $estfphits / $fphitstotal->{ $keys[$i] };
			}
			my $t2;
			unless ( $targethitstotal->{ $keys[$i] } )
			{
				$t2 = 0;
			} else
			{
				$t2 = $fphitstotal->{ $keys[$i] } / $targethitstotal->{ $keys[$i] };
			}
			my $t3;
			unless ( $targethitshash->{ $keys[$i] } )
			{
				$fdr = "n/a";
			} else
			{
				$t3 = $targethitshash->{ $keys[$i] } / $targethitstotal->{ $keys[$i] };
			}
			### Calculate the FDR
			$fdr = ( $t1 * $t2 ) / $t3;
		}
		### FDR * 10 is the precision
		$resultshash->{ $keys[$i] * 10 } = $fdr;
	}
	return $resultshash;
}

sub calc_local_fdr_xlinks
{
### get the hashes for rank 1
	my $targethitshash    = shift->{1};
	my $decoyhitshash     = shift->{1};
	my $fulldecoyshithash = shift->{1};
	my $targethitstotal   = shift->{1};
	my $fphitstotal       = shift->{1};
	my $fdr;
	my $resultshash = {};
	my $verbose     = 0;
## get the binkeys (Bins of scores!)
	my @keys = sort { $a <=> $b } keys %$targethitshash;

	for ( my $i = 0 ; $i < @keys ; $i++ )
	{
		## check if there are any target hits in this bin
		unless ( ( $targethitshash->{ $keys[$i] } ) )
		{
			$fdr = "n/a";
		} else
		{
			### Calculation of the FDR
			my $t1;
			### Est FP hits
			my $estfphits = $decoyhitshash->{ $keys[$i] } - 2 * $fulldecoyshithash->{ $keys[$i] };
			unless ( $fphitstotal->{ $keys[$i] } )
			{
				$t1 = 0;
			} else
			{
				$t1 = $estfphits / $fphitstotal->{ $keys[$i] };
			}
			my $t2;
			unless ( $targethitstotal->{ $keys[$i] } )
			{
				$t2 = 0;
			} else
			{
				$t2 = $fphitstotal->{ $keys[$i] } / $targethitstotal->{ $keys[$i] };
			}
			my $t3;
			unless ( $targethitshash->{ $keys[$i] } )
			{
				$fdr = "n/a";
			} else
			{
				$t3 = $targethitshash->{ $keys[$i] } / $targethitstotal->{ $keys[$i] };
			}
			### Calculate the FDR
			$fdr = ( $t1 * $t2 ) / $t3;
		}

		#$verbose && print "FDR: $fdr\n";
		### FDR: * 10 is the precision
		## check if fdr is not smaller than 0
		unless ( $fdr eq "n/a" )
		{
			if ( $fdr < 0 )
			{
				$fdr = 0;
			}
		}
		$resultshash->{ $keys[$i] * 10 } = $fdr;
	}
	return $resultshash;
}

sub get_stats
{
	my $score_to_test = shift;
	my $hashref       = shift;
	my $counter       = 0;
	foreach my $rank ( 1 .. $PARAMS->{'nranks'} )
	{
## extract the arrayref
		my $arrayref = $hashref->{$rank};
## make a cumulative statistik
## count all the elements greater or equal to a given score
		foreach my $score (@$arrayref)
		{
			if ( ( $score >= $score_to_test ) )
			{
				$counter++;
			}
		}
	}
	return $counter;
}

# Parse an XML file into a hash
sub parseXMLandAddFdr
{
	my $xmlfilename     = shift;
	my $intrafdr        = shift;
	my $interfdr        = shift;
	my $monolooplinkfdr = shift;
	my $spectrumidshash = shift;
	my $tree            = XML::TreeBuilder->new();
	$tree->parse_file($xmlfilename);
### Parsing all headers // one or more if a it is a merged result
	my @resultsheader = map { $_ } $tree->find('xquest_results');
## the running id for spectra start with 1 (0 is undef)
	my $spectrumid = 1;
	my $nseenbyxp  = 0;
## Parsing all specrum search results
	foreach my $header (@resultsheader)
	{
		my $headerid = $header->attr('outputpath');
		print "Annotating Header: $headerid\n";
		my $filespecindex = 0;
		## Define the filename of the spectrumxml file
		## parsing of all spectrum search results
		my @spectrumsearchelements = map { $_ } $header->find('spectrum_search');
		my $num = @spectrumsearchelements;

		#print "Found header element: $headerid with $num spectrum_search elements<br>";
		foreach my $spectrumsearchelements (@spectrumsearchelements)
		{
			my $spectrum_name  = $spectrumsearchelements->attr('spectrum');
			my $seenbyxprophet = 0;
			if ( $spectrumidshash->{$spectrum_name} )
			{
				$nseenbyxp++;
				$seenbyxprophet = 1;
			}
			### Get all search hits for this spectrum
			my @search_hits = map { $_ } $spectrumsearchelements->find('search_hit');
			my $i = 1;
			foreach my $hit (@search_hits)
			{
				### Gather information about the search hit
				my $rank  = $hit->attr('search_hit_rank');
				my $score = $hit->attr('score');
				my $error = $hit->attr('error_rel');
				my $type  = $hit->attr('type');
				my $id    = $hit->attr('id');

				#print "Search hit Score: $score\n";
				## Add the proteins
				my $spidp1 = $hit->attr('prot1');
				my $spidp2 = $hit->attr('prot2');
				### Determine the type of cross-link
				my $xlinktype = "-";
				if ( $type eq "xlink" )
				{
					$xlinktype = get_type_of_xlink( $spidp1, $spidp2 );
				} else
				{
					$xlinktype = get_type_of_monolink_looplink( $spidp1, $type );
				}
				#### Lookup the FDR
				my $fdr = "-";
				## 0.1 is the precision
				my $int = int( $score * 10 );

				#print "Score is $score, int to lookup is $int\n";
				if ( $xlinktype eq "intra/inter xl" )
				{

					#print "INTRALINK: ";
					#my $int= int ($score);
					$fdr = $intrafdr->{$int};

					#print "$fdr";
				}
				if ( $xlinktype eq "decoy intra/inter xl" )
				{

					#print "INTRALINK: ";
					#my $int= int ($score);
					$fdr = $intrafdr->{$int};

					#print "$fdr";
				}
				if ( $xlinktype eq "intra-protein xl" )
				{

					#print "INTRALINK: ";
					#my $int= int ($score);
					$fdr = $intrafdr->{$int};

					#print "$fdr";
				}
				if ( $xlinktype eq "decoy intra-protein xl" )
				{

					#print "INTRALINK: ";
					#my $int= int ($score);
					$fdr = $intrafdr->{$int};

					#print "$fdr";
				}
				if ( $xlinktype eq "inter-protein xl" )
				{

					#print "INTERLINK: ";
					#my $int= int ($score);
					$fdr = $interfdr->{$int};

					#print "$fdr";
				}
				if ( $xlinktype eq "decoy inter-protein xl" )
				{

					#print "INTERLINK: ";
					#my $int= int ($score);
					$fdr = $interfdr->{$int};

					#print "$fdr";
				}
				if ( $xlinktype eq "monolink" )
				{

					#print "INTERLINK: ";
					#my $int= int ($score);
					$fdr = $monolooplinkfdr->{$int};

					#print "$fdr";
				}
				if ( $xlinktype eq "decoy monolink" )
				{

					#print "INTERLINK: ";
					#my $int= int ($score);
					$fdr = $monolooplinkfdr->{$int};

					#print "$fdr";
				}
				if ( $xlinktype eq "intralink" )
				{

					#print "INTERLINK: ";
					#my $int= int ($score);
					$fdr = $monolooplinkfdr->{$int};

					#print "$fdr";
				}
				if ( $xlinktype eq "decoy intralink" )
				{

					#print "INTERLINK: ";
					#my $int= int ($score);
					$fdr = $monolooplinkfdr->{$int};

					#print "$fdr";
				}
				$hit->attr( 'fdr',        $fdr );
				$hit->attr( 'xprophet_f', $seenbyxprophet );
				$i++;
			}
			$spectrumid++;
			$filespecindex++;
		}
	}
	print "  Annotated xml with FDR; $nseenbyxp hits were seen by xProphet and flagged\n";
	return $tree;
}

sub store_fdr_hash
{
	my $fh      = shift;
	my $fdrhash = shift;
	my $desc    = shift;
	my @keys    = sort { $a <=> $b } keys %$fdrhash;
	print $fh $desc . "\t";
	for ( my $i = 0 ; $i < @keys ; $i++ )
	{
		print $fh $fdrhash->{ $keys[$i] } . "\t";
	}
	print $fh "\n\n";
	print "FDR results for $desc written\n";
}

sub calc_fdr
{
	my $tphash      = shift;
	my $fphash      = shift;
	my $PARAMS      = shift;
	my $f           = shift;
	my $resultshash = {};

	#	foreach my $rank ( 1 .. $PARAMS->{'nranks'} )
	#	{
	## This is only for rank 1
	my $tpsubhash = $tphash->{1};
	my $fpsubhash = $fphash->{1};
## get the binkeys (Bins of scores!)
	my @keys = sort { $a <=> $b } keys %$tpsubhash;
	for ( my $i = 0 ; $i < @keys ; $i++ )
	{
		my $fdr;
		unless ( defined( $fpsubhash->{ $keys[$i] } ) )
		{

			#die "Error: Key $i does not exist in the FP hash but in the TP hash!\n";
			$fdr = "n/a";
		} else
		{
			my $fp = $fpsubhash->{ $keys[$i] } * $f;
			if ( ( $fp + $tpsubhash->{ $keys[$i] } ) == 0 )
			{
				$fdr = "n/a";
			} else
			{

				#print "KEY:$i\n";
				#print "FP: ". $fpsubhash->{$i}."\n";
				#print "TP: ". $tpsubhash->{$i}."\n";
				$fdr = $fp / ( $fp + $tpsubhash->{ $keys[$i] } );
			}
		}
		$resultshash->{ $keys[$i] } = $fdr;
	}
	return $resultshash;
}

#---------------------------------------------------------------------------
#  SUB parseXML for PARSING THE XML FILE
#---------------------------------------------------------------------------
sub parse_spec_xml_file
{
	my $xmlfilename = shift;
	my $outfilename = $xmlfilename . ".hash";
	if ( -e $outfilename )
	{
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
	$tree->delete();
}

sub get_clean_id
{
	my $protid = shift;

	#print "Protein Id:$protid\n";
### Split the protein id by _ "eg: decoy_reverse_gi|147905534|ref|NP_001079812.1|"
### splits into decoy, reverse and the id: the id is always the last
	my @splittedp = split( /\_/, $protid );
## reconstitute if the part doesnt mach "reverse or decoy"
	my @reconstitute;
	foreach my $part (@splittedp)
	{

		#print "part: $part\n";
		unless ( $part =~ /decoy/ || $part =~ /reverse/ )
		{
			push @reconstitute, $part;
		}
	}
### reconstitute
	my $rec = join( "_", @reconstitute );

	#print "Reconstituted and clean id: $rec\n";
	return $rec;
}

sub get_type_of_monolink_looplink
{
	my $spidp1 = shift;
	my $type   = shift;
	my $typetoreport;
	my @prots1 = split( ",", $spidp1 );
	my $decoy = 0;
	foreach my $prot1 (@prots1)
	{
### Check if protein is a decoy protein
		if ( $prot1 =~ m/decoy/ )
		{
			$decoy = 1;
		}
	}
	if ($decoy)
	{
		$typetoreport = "decoy " . $type;
	} else
	{
		$typetoreport = $type;
	}
	return $typetoreport;
}

sub get_type_of_xlink
{
	my $spidp1 = shift;
	my $spidp2 = shift;
	my @prots1 = split( ",", $spidp1 );
	my @prots2 = split( ",", $spidp2 );
	my $type;
	my $intralink = 0;
	my $interlink = 0;
	my $decoy     = 0;
	foreach my $prot1 (@prots1)
	{
		my $p1 = get_clean_id($prot1);

		#exit;
		foreach my $prot2 (@prots2)
		{
			my $p2 = get_clean_id($prot2);

			#print "Prot1: $prot1 , Prot2: $prot2\n";
			if ( $p1 eq $p2 )
			{
				$intralink = 1;
			} else
			{
				$interlink = 1;
			}
			### Check if one protein is a decoy protein
			if ( $prot1 =~ m/decoy/ || $prot2 =~ m/decoy/ )
			{
				$decoy = 1;
			}
		}
	}
	if ( $intralink && $interlink )
	{
		$type = "intra/inter xl";
		if ($decoy)
		{
			$type = "decoy intra/inter xl";
		}
	}
	if ( $intralink && !$interlink )
	{
		$type = "intra-protein xl";
		if ($decoy)
		{
			$type = "decoy intra-protein xl";
		}
	}
	if ( !$intralink && $interlink )
	{
		$type = "inter-protein xl";
		if ($decoy)
		{
			$type = "decoy inter-protein xl";
		}
	}
	return $type;
}

sub store_histogram
{
	my $fh               = shift;
	my $hashref          = shift;
	my $PARAMS           = shift;
	my $desc             = shift;
	my $resheaderprinted = 0;
	foreach my $rank ( 1 .. $PARAMS->{'nranks'} )
	{
		my @resarray;
		my @resheader;
		my $a1 = "Rank (vert.),Bins(horiz.) $rank\t";
		push @resheader, $a1;
		push @resarray,  $rank . "\t";
		my $subhash = $hashref->{$rank};
## get the binkeys
		my @keys = sort { $a <=> $b } keys %$subhash;

		#print Dumper (@keys);
		#exit;
		my $sum = 0;

		#my @resvalues;
		#push @resarray, "\t";
		for ( my $i = 0 ; $i < @keys ; $i++ )
		{
			push @resheader, $keys[$i] . "\t";
			my $value = $subhash->{ $keys[$i] };
			$sum = $sum + $value;
			push @resarray, $value . "\t";
		}
		push @resheader, "SUM\n";
		push @resarray,  $sum . "\n";
		unless ($resheaderprinted)
		{
			print $fh $desc . "\n";
			print $fh @resheader;
			$resheaderprinted = 1;
		}
		print $fh @resarray;
	}
	print $fh "\n";
}
### generates hashes with the scores in it for all the diffenent types
sub _generate_result_hashs
{
	my $resulthash = shift;
	my $PARAMS     = shift;
	my @ukeys;
	foreach my $id ( keys %$resulthash )
	{
		push @ukeys, $id;
	}
	### Sort the keys descending by the score of the hits
	my @keys                   = sort { $resulthash->{$b}->{1}->{'score'} <=> $resulthash->{$a}->{1}->{'score'} } @ukeys;
	my $interlinks             = {};
	my $interdecoys            = {};
	my $intralinks             = {};
	my $intradecoys            = {};
	my $monolinks              = {};
	my $monodecoys             = {};
	my $looplinks              = {};
	my $looplinkdecoys         = {};
	my $decoys                 = {};                                                                                        ## all decoys
	my $targethits             = {};
	my $intralinksIds          = {};
	my $hybriddecoysinterlinks = {};                                                                                        ## the hybrid decoys
	my $fulldecoysinterlinks   = {};
	my $fulldecoysintralinks   = {};
	my $hybriddecoysintralinks = {};
	my $spectrumidshash        = {};
######## HASH Structure
	my $msg = "Number of elements:" . scalar(@keys) . "\n";
	print $msg;
	my $counter = 0;
	my $seen    = {};

	foreach my $key (@keys)
	{
		foreach my $rank ( 1 .. $PARAMS->{'nranks'} )
		{

			#print "KEY:".$key." Rank: ".$rank."\n";
			my $spectrumname = $resulthash->{$key}->{$rank}->{'spectrumname'};
			my $xltype       = $resulthash->{$key}->{$rank}->{'xlinktype'};
			my $score = $resulthash->{$key}->{$rank}->{'score'};
			
			### Filter by mass deviation
			my $delta = $resulthash->{$key}->{$rank}->{'error'};
			if ( ( $delta >= $PARAMS->{'minborder'} ) && ( $delta <= $PARAMS->{'maxborder'} ) )
			{

				#do nothing
			} else
			{
				next;
			}
			### Filter by deltaS
			my $deltaS = $resulthash->{$key}->{$rank}->{'deltaS'};
			if ( $PARAMS->{'mindeltas'} )
			{
				if ( $deltaS >= $PARAMS->{'mindeltas'} )
				{

					#print "DeltaS is: $deltaS, cutoff is ",$PARAMS->{'mindeltas'},"\n";
					next;
				}
			}
			### Filter by minmatched ions
			my $mions = $resulthash->{$key}->{$rank}->{'nminmatchedions'};
			if ( $PARAMS->{'minionsmatched'} )
			{
				if ( $mions < $PARAMS->{'minionsmatched'} )
				{
					next;
				}
			}
			
			## Filter by minscore
			if ( $PARAMS->{'minscore'}	)
			{
				if ( $score < $PARAMS->{'minscore'} )
				{
					next;
				}	
			}
			
			## Filter unique xls
			### Filter by unique Ids (starts now with the largest score)
			my $myid = $resulthash->{$key}->{$rank}->{'id'};
			if ( $PARAMS->{'uniquexl'} )
			{
				unless ( ( $seen->{$myid} ) )
				{
					$seen->{$myid} = 1;
				} else
				{
					next;
				}
			}
			
			
			### Get the type of cross-link
			my $type = $resulthash->{$key}->{$rank}->{'type'};
			unless ($type)
			{
				next;
			}

			#print "Type: " . $type . "\n";
			
			my $prot1 = $resulthash->{$key}->{$rank}->{'protein1'};
			my $prot2 = $resulthash->{$key}->{$rank}->{'protein2'};
			if ( ( $xltype =~ m/intra-protein xl/ ) || ( $xltype =~ m/intra\/inter xl/ ) )
			{
				if ( $xltype =~ m/decoy/ )
				{
					push @{ $intradecoys->{$rank} }, $score;
					push @{ $decoys->{$rank} },      $score;
					if ( $prot1 =~ m/decoy/ && $prot2 =~ m/decoy/ )
					{
						push @{ $fulldecoysintralinks->{$rank} }, $score;
					} else
					{
						push @{ $hybriddecoysintralinks->{$rank} }, $score;
					}
				} else
				{
					push @{ $intralinks->{$rank} },    $score;
					push @{ $intralinksIds->{$rank} }, "S: " . $score . " DeltsS:" . $deltaS . " Type: " . $xltype . " Id: " . $myid;
					push @{ $targethits->{$rank} },    $score;
				}
			}
			if ( ( $xltype =~ m/inter-protein xl/ ) )
			{
				## then check if its not a decoy
				if ( $xltype =~ m/decoy/ )
				{

					#	print "Decoy found: $prot1\n";
					push @{ $decoys->{$rank} },      $score;
					push @{ $interdecoys->{$rank} }, $score;
					if ( $prot1 =~ m/decoy/ && $prot2 =~ m/decoy/ )
					{
						push @{ $fulldecoysinterlinks->{$rank} }, $score;
					} else
					{
						push @{ $hybriddecoysinterlinks->{$rank} }, $score;
					}
				} else
				{

					#print "interlink found: $prot1\n";
					push @{ $interlinks->{$rank} }, $score;
					push @{ $targethits->{$rank} }, $score;

					#push @{ $inter_and_decoys->{$rank} }, $score;
				}
			}
			if ( $xltype =~ m/monolink/ || $xltype =~ m/intralink/ )
			{
				if ( $xltype =~ m/decoy/ )
				{
					push @{ $decoys->{$rank} },     $score;
					push @{ $monodecoys->{$rank} }, $score;
				} else
				{
					push @{ $monolinks->{$rank} },  $score;
					push @{ $targethits->{$rank} }, $score;
				}
			}
			$spectrumidshash->{$spectrumname} = 1;
		}
	}
	return $intralinks, $intradecoys, $interlinks, $interdecoys, $decoys, $targethits, $monolinks, $monodecoys, $looplinks, $looplinkdecoys, $hybriddecoysinterlinks, $fulldecoysinterlinks, $fulldecoysintralinks, $hybriddecoysintralinks, $intralinksIds, $spectrumidshash;
}
## generates a cummulative hash
sub _gen_grouped_hash
{
	my $groupsarray = shift;
	my $xlinkshash  = shift;
	my $PARAMS      = shift;
	my $resultarray = {};
	foreach my $rank ( 1 .. $PARAMS->{'nranks'} )
	{
## extract the arrayref
		my $arrayref = $xlinkshash->{$rank};
## make a cumulative statistik
## count all the elements greater or equal to a given score
		for ( my $i = 0 ; $i < scalar(@$groupsarray) ; $i++ )
		{
			my $counter  = 0;
			my $counter2 = 0;
			foreach my $score (@$arrayref)
			{
				$counter2++;
				if ( ( $score >= $groupsarray->[$i] ) )
				{
					$counter++;
				}
			}
			$resultarray->{$rank}->{ $groupsarray->[$i] } = $counter;
		}
	}
	return $resultarray;
}
## generates a histogram hash!
sub _gen_histo_hash
{
	my $groupsarray = shift;
	my $xlinkshash  = shift;
	my $PARAMS      = shift;
	my $resultarray = {};
	foreach my $rank ( 1 .. $PARAMS->{'nranks'} )
	{
## extract the arrayref
		my $arrayref = $xlinkshash->{$rank};
## make a histogram statistic
## count all the elements greater or equal to a given score
		for ( my $i = 0 ; $i < scalar(@$groupsarray) ; $i++ )
		{
			my $counter = 0;
			foreach my $score (@$arrayref)
			{
				unless ($score)
				{
					warn "It seems that $score is not set or is 0 \n";
					next;
				}
				if ( ( $score >= $groupsarray->[$i] ) && ( $score < $groupsarray->[ $i + 1 ] ) )
				{
					$counter++;
				}
			}
			$resultarray->{$rank}->{ $groupsarray->[$i] } = $counter;
		}
	}
	return $resultarray;
}

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

#---------------------------------------------------------------------------
#  SUB parseXML for PARSING THE XML FILE
#---------------------------------------------------------------------------
sub read_spec_xml_file
{
	my $xmlfilename = shift;
	my $spechash    = {};
	my $tree        = XML::TreeBuilder->new();
	print "Parsing mzXML  $xmlfilename\n";
	$tree->parse_file($xmlfilename);
	my @resultsheader = map { $_ } $tree->find('xquest_spectra');
	### index the spectra by filename
	foreach my $header (@resultsheader)
	{
		my @spectra = map { $_ } $header->find('spectrum');
		foreach my $spectrum (@spectra)
		{
			my $filename = $spectrum->attr('filename');
			my $content  = $spectrum->content();
			$spechash->{$filename} = $content->[0];
		}
	}
	$tree->delete();
	return $spechash;
}

# Parse an XML file into a hash
sub parseXML
{
	my $xmlfilename = shift;
	my $resultshash = shift;
	my $summaryhash = shift;
### Parse the XML file
	unless ( -e $xmlfilename )
	{
		die " Error: Cannot find xquest result file $xmlfilename\n";
	}
	$| = 1;    ## omit output buffering
	print "Parsing xquest XML file\n";
	my $tree = XML::TreeBuilder->new();
	$tree->parse_file($xmlfilename);
### Parsing all headers // one or more if a it is a merged result
	my @resultsheader = map { $_ } $tree->find('xquest_results');
## the running id for spectra start with 1 (0 is undef)
	my $spectrumid = 1;
## ref for files
	my $files = $summaryhash->{'xmlheader'};
## Parsing all specrum search results
	foreach my $header (@resultsheader)
	{
		my $headerid = $header->attr('outputpath');
		print "Indexing Header: $headerid\n";
		my $filespecindex = 0;
		## Define the filename of the spectrumxml file
		## parsing of all spectrum search results
		my @spectrumsearchelements = map { $_ } $header->find('spectrum_search');
		my $num = @spectrumsearchelements;

		#print "Found header element: $headerid with $num spectrum_search elements<br>";
		foreach my $spectrumsearchelements (@spectrumsearchelements)
		{
			## Get the spectrum name
			my $spectrum_name = $spectrumsearchelements->attr('spectrum');
			unless ($spectrum_name)
			{
				print "  Error ($0): Can't find the spectrum name\n";
				exit;
			}
			### Get all search hits for this spectrum
			my @search_hits = map { $_ } $spectrumsearchelements->find('search_hit');
			my $i = 1;
			### calculate the delta scores and set it as attribute for each hit
			calc_delta_score( \@search_hits );
			foreach my $hit (@search_hits)
			{
				### Gather information about the search hit
				my $rank        = $hit->attr('search_hit_rank');
				my $score       = $hit->attr('score');
				my $error       = $hit->attr('error_rel');
				my $type        = $hit->attr('type');
				my $id          = $hit->attr('id');
				my $deltas      = $hit->attr('deltaS');
				my $nmionsalpha = $hit->attr('num_of_matched_ions_alpha');
				my $nmionsbeta  = $hit->attr('num_of_matched_ions_beta');
				### Only store rank 1
				unless ( $rank == 1 )
				{
					next;
				}
				### Store in spec hash
				$resultshash->{$spectrumid}->{$rank}->{'spectrumname'}    = $spectrum_name;
				$resultshash->{$spectrumid}->{$rank}->{'score'}           = $score;
				$resultshash->{$spectrumid}->{$rank}->{'error'}           = $error;
				$resultshash->{$spectrumid}->{$rank}->{'type'}            = $type;
				$resultshash->{$spectrumid}->{$rank}->{'id'}              = $id;
				$resultshash->{$spectrumid}->{$rank}->{'deltaS'}          = $deltas;
				$resultshash->{$spectrumid}->{$rank}->{'nminmatchedions'} = get_n_min_ions_matched( $type, $nmionsalpha, $nmionsbeta );
				## Add the proteins
				my $spidp1 = $hit->attr('prot1');
				$resultshash->{$spectrumid}->{$rank}->{'protein1'} = $spidp1;
				my $spidp2 = $hit->attr('prot2');
				$resultshash->{$spectrumid}->{$rank}->{'protein2'} = $spidp2;
				### Determine the type of cross-link
				my $xlinktype = "-";
				if ( $type eq "xlink" )
				{
					$xlinktype = get_type_of_xlink( $spidp1, $spidp2 );
				} else
				{
					$xlinktype = get_type_of_monolink_looplink( $spidp1, $type );
				}
				$resultshash->{$spectrumid}->{$rank}->{'xlinktype'} = $xlinktype;
				$i++;
			}
			### STORE THE REF INTO THE TIED HASH
			$spectrumid++;
			$filespecindex++;
			$spectrumsearchelements->delete();
		}
		##
		push @{$files}, $headerid . "\t" . $filespecindex . "\n";
	}
	$summaryhash->{'xmlheader'} = $files;
	print "Parsed ", $spectrumid - 1, " spectra\n";
	$summaryhash->{'totalparsed'} .= "Parsed " . scalar( @{$files} ) . "\n" . ( $spectrumid - 1 ) . " spectra\n";

	#return $tree;
	$tree->delete();

	#print Dumper ($resultshash);
}

sub get_n_min_ions_matched
{
	my $type           = shift;
	my $nmionsalpha    = shift;
	my $nmionsbeta     = shift;
	my $minionsmatched = 0;
	if ( $type eq "xlink" )
	{
		if ( $nmionsalpha == $nmionsbeta )
		{
			$minionsmatched = $nmionsalpha;
		} elsif ( $nmionsalpha < $nmionsbeta )
		{
			$minionsmatched = $nmionsalpha;
		} elsif ( $nmionsbeta < $nmionsalpha )
		{
			$minionsmatched = $nmionsbeta;
		}
	} else
	{
		$minionsmatched = $nmionsalpha;
	}
	return $minionsmatched;
}

sub calc_delta_score
{
	my $search_hits = shift;
## Sort the search hits by the rank
	my @allhitssorted = sort { $a->attr('search_hit_rank') <=> $b->attr('search_hit_rank') } (@$search_hits);

	#---------------------------------------------------------------------------
	#  calculate the delta score here and set them as attribute
	#---------------------------------------------------------------------------
	#my @allhitssorted=sort { $b->attr('score') <=> $a->attr('score') } @search_hits;
	for ( my $i = 0 ; $i < ( scalar(@allhitssorted) ) ; $i++ )
	{

		#my $rank=$i+1;
		my $hit = $allhitssorted[$i];
		my $resultshashref;
		my $key;
		my $rank      = $hit->attr('search_hit_rank');
		my $score     = $hit->attr('score');
		my $error_rel = $hit->attr('error');
		my $id        = $hit->attr('structure');

		#print "Current Hit Rank $rank, Structure $id, Score $score, Error is $error_rel ppm\n";
		my $nextdifferentid;
		my $nextdifferentidscore;
### Get the next id that is different from the current start at the current id +1
		for ( my $z = $i ; $z < scalar(@allhitssorted) ; $z++ )
		{
			my $hit2 = $allhitssorted[ $z + 1 ];
			unless ($hit2)
			{

				#	print "Hit2 not set\n";
				last;
			}
			my $rank      = $hit2->attr('search_hit_rank');
			my $score     = $hit2->attr('score');
			my $error_rel = $hit2->attr('error');
			my $hit2id    = $hit2->attr('structure');

			#print "Current Hit Rank $rank, Structure $hit2id, Score $score, Error is $error_rel ppm\n";
### Check if the structure is the same
			if ( $id eq $hit2id )
			{

				#print "Hits are the same!\n";
				next;
			} else
			{
				$nextdifferentid      = $z;
				$nextdifferentidscore = $hit2->attr('score');
				last;
			}
		}
### Calculate the delta btw the hits
		my $delta;
		unless ($score)
		{
			$delta = "n/a";
		} else
		{
			if ( $score == 0 )
			{
				$delta = "n/a";
			} else
			{
				$delta = $nextdifferentidscore / $score;
			}
		}

		#print "DeltaScore for hit $rank: $delta";
		## Set the deltascore
		$hit->attr( 'deltaS', $delta );
	}
}

sub print_stats
{
	if ( @{ $interlinks->{1} } )
	{
		print "INTERLINKS (top rank):\t", @{ $interlinks->{1} } . "\n";
		$PARAMS->{'statusmsg'} .= "INTERLINKS (top rank):\t" . @{ $interlinks->{1} } . "\n";
	} else
	{
		print "INTERLINKS (top rank):\t0\n";
		$PARAMS->{'statusmsg'} .= "INTERLINKS (top rank):\t0\n";
	}
	if ( @{ $interlinks->{1} } )
	{
		print "Decoy INTERLINKS (top rank):\t", @{ $interdecoys->{1} } . "\n";
		$PARAMS->{'statusmsg'} .= "Decoy INTERLINKS (top rank):\t" . @{ $interdecoys->{1} } . "\n";
	} else
	{
		print "Decoy INTERLINKS (top rank):\t0\n";
		$PARAMS->{'statusmsg'} .= "Decoy INTERLINKS (top rank):\t0\n";
	}
	if ( @{ $intralinks->{1} } )
	{
		print "INTRALINKS (top rank):\t" . @{ $intralinks->{1} } . "\n";
		$PARAMS->{'statusmsg'} .= "INTRALINKS (top rank):\t" . @{ $intralinks->{1} } . "\n";
	} else
	{
		print "INTRALINKS (top rank):\t0\n";
		$PARAMS->{'statusmsg'} .= "INTRALINKS (top rank):\t0\n";
	}
	if ( @{ $intradecoys->{1} } )
	{
		print "Decoy INTRALINKS (top rank):\t" . @{ $intradecoys->{1} } . "\n";
		$PARAMS->{'statusmsg'} .= "Decoy INTRALINKS (top rank):\t" . @{ $intradecoys->{1} } . "\n";
	} else
	{
		print "Decoy INTRALINKS (top rank):\t0\n";
		$PARAMS->{'statusmsg'} .= "Decoy INTRALINKS (top rank):\t0\n";
	}
	if ( @{ $monolinks->{1} } )
	{
		print "Monolinks (top rank):\t" . @{ $monolinks->{1} } . "\n";
		$PARAMS->{'statusmsg'} .= "Monolinks (top rank):\t" . @{ $monolinks->{1} } . "\n";
	} else
	{
		print "Monolinks (top rank):\t0\n";
		$PARAMS->{'statusmsg'} .= "Monolinks (top rank):\t0\n";
	}
	if ( @{ $monodecoys->{1} } )
	{
		print "Decoy Monolinks (top rank):\t" . @{ $monodecoys->{1} } . "\n";
		$PARAMS->{'statusmsg'} .= "Decoy Monolinks (top rank):\t" . @{ $monodecoys->{1} } . "\n";
	} else
	{
		print "Decoy Monolinks (top rank):\t0\n";
		$PARAMS->{'statusmsg'} .= "Decoy Monolinks (top rank):\t0\n";
	}
	if ( @{ $looplinks->{1} } )
	{
		print "Looplinks (top rank):\t" . @{ $looplinks->{1} } . "\n";
		$PARAMS->{'statusmsg'} .= "Looplinks (top rank):\t" . @{ $looplinks->{1} } . "\n";
	} else
	{
		print "Looplinks (top rank):\t0\n";
		$PARAMS->{'statusmsg'} .= "Looplinks (top rank):\t0\n";
	}
	if (  @{ $looplinkdecoys->{1} } )
	{
		print "Decoy Looplinks (top rank):\t" . @{ $looplinkdecoys->{1} } . "\n";
		$PARAMS->{'statusmsg'} .= "Decoy Looplinks (top rank):\t" . @{ $looplinkdecoys->{1} } . "\n";
	} else
	{
		print "Decoy Looplinks (top rank):\t0\n";
		$PARAMS->{'statusmsg'} .= "Decoy Looplinks (top rank):\t0\n";
	}
	if ( @{ $decoys->{1} } )
	{
		print "DECOYS (top rank):\t" . @{ $decoys->{1} } . "\n";
		$PARAMS->{'statusmsg'} .= "DECOYS (top rank):\t" . @{ $decoys->{1} } . "\n";
	} else
	{
		print "DECOYS (top rank):\t0\n";
		$PARAMS->{'statusmsg'} .= "DECOYS (top rank):\t0\n";
	}
	if ( @{ $targethits->{1} } )
	{
		print "TARGET HITS (top rank):\t" . @{ $targethits->{1} } . "\n";
		$PARAMS->{'statusmsg'} .= "TARGET HITS (top rank):\t" . @{ $targethits->{1} } . "\n";
	} else
	{
		print "TARGET (top rank):\t0\n";
		$PARAMS->{'statusmsg'} .= "TARGET (top rank):\t0\n";
	}
}

sub usage()
{
	print "
	SOFTWARE: ", basename($0), " version $version
	
	AUTHOR: Thomas Walzthoeni

	INFORMATION: A software/script to estimate FDRs for cross-linked peptides identified by MS.
	Please note: Protein Identifiers must contain the string \"decoy\" to be recognized.
 	
 	USAGE: ", basename($0), " -Option [Parameter]

	REQUIRED OPTIONS [defaults]: 
	-in [xquest.xml] xquest.xml input filename
	-out [xproph_out.xml] output filename, can be used by the xQuest/xProphet viewer.
 	-def [xproph.def] xProphet definition file, see DEFINTION FILE for further information

	DEFINTION FILE
	To generate a definition file (xproph.def) template execute the program without any options (and no xproph.def in the working directory), or use -getdef.
	
	OTHER OPTIONS:
 	-help print this information
 	

	EXAMPLE: ", basename($0), " -in merged_xquest.xml -out xquest.xml
 		
	";
	exit;
}
