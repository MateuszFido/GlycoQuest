package ParseMap;
use strict;
#---------------------------------------------------------------------------
# Module: ParseMap.pm
# Author(s): Thomas Walzthoeni, xQuest specific modules are based on 
# original work by Oliver Rinner.
# Description: Module for parsing master map (feature map).
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
use Feature;
use XMM_Statistics;
use proteingroups;
use XML::TreeBuilder;
use mzXMLscan;
use XMM_ETD_CID;
use inclusion;
use Isopairs;
use File::Basename;
use Data::Dumper;

sub new {
	my $class     = shift();
	my $MasterMap = shift;
	my $PARAMS    = shift;

	my $self = {};
	bless $self, $class;

	my $translation = $PARAMS->translationtable;
	my $verbose     = $PARAMS->verbose;

	if ($translation) {
		print "using $translation as translationtable\n";
		$self->{'translation'} = $self->parsetranslationtable($translation);
	}
	else {
		$self->{'translation'} = undef;
	}

	$self->{'PARAMS'}     = $PARAMS;
	$self->{'verbose'}    = $verbose;
	$self->{'normfactor'} = 1;
	$self->{'runIDs'}     = $PARAMS->runIDs;

	open MASTER, "<$MasterMap" or die "cannot open Master Map $!";
	my $startmaster        = undef;
	my @runIDs             = ();
	my %runIDS_translation = ();
	my @features           = ();
	while (<MASTER>) {
		chomp;
		if ( $_ =~ /B0\d-\d+_p\[0\]/ ) {
			s/[\s+\(\)]//g;
			print $_, "\n";
			%runIDS_translation =
			  map { $2, $1 if /^(.*)\[(\d+)\]/g }
			  split /\|/;
			#exit;
			#print join "\t",map {$_} %runIDS_translation,"\n";
			next;
		}
		if (/M\/z\tTr.*(AREA)/) {
			while ( $_ =~ /AREA_(\d+)/gi ) {
				push @runIDs, $1;
			}
			print "ids in MasterMap: @runIDs\n";
			$startmaster = defined;

			next;

		}
		$self->{'definedIDs'}             = \@runIDs;
		$self->{'definedIDs_translation'} = \%runIDS_translation;

		if ($startmaster) {
		
			my @feature = split /\t/, $_;
			if ($translation) {
				if ( $feature[3] ) {
					$self->translationtable->{ $feature[3] }
					  && ( $feature[3] =
						$self->translationtable->{ $feature[3] } );
				}
			}
			#print join "\t", @runIDs, "\n";
			#print join "\t", @feature, "\n";
			$verbose && print join "\t", @feature, "\n";
			
			my $feature = Feature->new( $self, \@feature, \@runIDs );

			if ( $feature->matches_requirements ) {
				push @features, $feature;
			}

		}
		else {
			next;
		}
	}
	print "Read ", scalar(@features), " Features that match the requirements from MasterMap.\n";
	$self->{'features'} = \@features;
	return $self;
}

sub getrunIDs {
	my $self = shift;
	return $self->{'runIDs'};
}


sub transferscans_from_reference {
	my $self             = shift;
	my $PARAMS           = shift;
	my $outputfilehandle = shift;

	unless ($outputfilehandle) {
		$outputfilehandle = *STDOUT;
	}

	my $precursor_tolerance   = $PARAMS->etd_precursor_tolerance;
	my $etd_Tr_diff           = $PARAMS->etd_Tr_diff;
	my $dot_product_threshold = $PARAMS->etd_min_dotproduct;

	#my $isopairs_with_scans=$self->get_isopairs_with_any_scans;

	#	my $features_with_scan  =
	#	  $self->get_features_with_scan;    #features that got a primary scan

	my $features_with_scan =    #only features from isotopic pairs
	  $self->get_isopairs_features_with_scans;

	my $reference_scanhash = $self
	  ->get_reference_scan_mz_hash;    #scans from reference ordered by int m/z

	print $outputfilehandle "Feature_mz\tFeature_Tr\tassociated_scan_mz\tassociated_scan_Tr\tscan_ID\tdotproduct\n";

	print "mapping MS2 features from reference LC-MS runs to MS2 features assigned to isotopic pairs...";
	my $npairsassigned = 0;
	foreach my $feature (@$features_with_scan) {
		my @candidatescans = ();
		my $feature_mz     = $feature->mz;
		my $feature_Tr     = $feature->Tr;

		my $referencescan = $feature->associated_scan;

		if ( defined( $reference_scanhash->{ int($feature_mz) } ) ) {
			push @candidatescans,
			  @{ $reference_scanhash->{ int($feature_mz) } };
		}
		if ( defined( $reference_scanhash->{ int( $feature_mz + 1 ) } ) ) {
			push @candidatescans,
			  @{ $reference_scanhash->{ int( $feature_mz + 1 ) } };
		}
		if ( defined( $reference_scanhash->{ int( $feature_mz - 1 ) } ) ) {
			push @candidatescans,
			  @{ $reference_scanhash->{ int( $feature_mz - 1 ) } };
		}

		#print "feature mz: ", $feature->mz, " feature Tr ", $feature->Tr, "\n";

		foreach my $assigned_etd_scanobj (@candidatescans) {
			my $cid_scan = $assigned_etd_scanobj->cid_scan;
			if (
				(
					abs( $cid_scan->precursor_mz - $feature_mz ) <
					$precursor_tolerance
				)

				&& (
					( abs( $feature_Tr - $cid_scan->Tr ) / 60 ) < $etd_Tr_diff )
			  )

			{
				my $dotproduct =
				  Consensus::dotproduct( $referencescan->get_peaklist,
					$cid_scan->get_peaklist, $PARAMS );
				print $outputfilehandle $feature_mz, "\t", $feature_Tr, "\t",
				  $cid_scan->precursor_mz, "\t", $cid_scan->Tr, "\t",
				  $cid_scan->id, "\t$dotproduct\n";

				if ( $dotproduct > $dot_product_threshold ) {
					$npairsassigned++;
					$feature->assign_etd_cid_pair($assigned_etd_scanobj);
				}
			}
		}
	}
	print
	  " ...done, assigned $npairsassigned ETD MS2 scans to isotopic pairs\n";
}

sub makeinclusionlist {
	my $self              = shift;
	my $inclusionlistfile = shift;
	my $inclusionlist     = inclusion->new( $self, $self->params );
	$inclusionlist->make_list($inclusionlistfile);
}

sub get_isotopicpairs {
	my $self = shift;
	return $self->{'isotopicpairs'};
}

sub getdefinedIDs {
	my $self = shift;
	return $self->{'definedIDs'};
}

sub get_definedIDs_translation {
	my $self = shift;
	return $self->{'definedIDs_translation'};
}

sub getrunIDarray {
	my $self = shift;
	return split /,/, $self->{'runIDs'};
}

sub params {
	my $self = shift;
	return $self->{'PARAMS'};
}

sub min_Tr {
	my $self = shift;
	return $self->{'min_Tr'};
}

sub max_Tr {
	my $self = shift;
	$self->{'max_Tr'};
}

sub calc_min_max_Tr {
	my $self           = shift;
	my $features       = $self->getallfeatures;
	my @sortedfeatures = sort { $a->Tr <=> $b->Tr } @$features;
	$self->{'min_Tr'} = $sortedfeatures[0]->Tr;
	$self->{'max_Tr'} = $sortedfeatures[-1]->Tr;
}

sub getallfeatures {
	my $self = shift;
	return $self->{'features'};
}

sub get_features_with_scan {
	my $self             = shift;
	my $allfeatures      = $self->getallfeatures;
	my @featureswithscan = ();
	foreach my $feature (@$allfeatures) {
		if ( $feature->scan ) {
			push @featureswithscan, $feature;
		}
	}
	return \@featureswithscan;
}

sub get_features_wo_scan {
	my $self           = shift;
	my $allfeatures    = $self->getallfeatures;
	my @featureswoscan = ();
	foreach my $feature (@$allfeatures) {
		unless ( $feature->scan ) {
			push @featureswoscan, $feature;
		}
	}
	return \@featureswoscan;
}

sub normfactor {
	my $self = shift;
	return $self->{'normfactor'};
}

sub parse_ETD_mzXML {
	my $self       = shift;
	my $mzXMLfiles = shift;
	my $getpeaks   = shift;

	my $PARAMS          = $self->params;
	my $levels          = $PARAMS->scanlevels;
	my $msInstrumentIDs = $PARAMS->msInstrumentIDs;
	my %scanmzhash      = ();
	my %scan_pair_hash  = ();
	my $nscans          = 0;
	my $i               = 0;
	my $mzXMLfilename;

	if ( scalar @$mzXMLfiles == 0 ) {
		die "no mzXML file was indicated \n";
	}

	foreach my $mzXMLfile (@$mzXMLfiles) {

		my $mzXMLtree = XML::TreeBuilder->new();
		$mzXMLtree->parse_file($mzXMLfile);
		( $mzXMLfilename = basename($mzXMLfile) ) =~ s/\.mzXML//;
		print "parsing $mzXMLfilename...\n";
		my $dependent_ID = 0;
		foreach my $precursorscan (
			$mzXMLtree->look_down(
				_tag => 'scan',
				sub { $_[0]->attr('msLevel') < 2 }
			)
		  )
		{
			foreach my $scan ( $precursorscan->find('scan') ) {

				if ( number_is_in_array( $levels, $scan->attr_get_i('msLevel') )
				  )
				{
					$nscans++;
					my $mzxmlscan =
					  mzXMLscan->new( $scan, $i, $mzXMLfilename, $getpeaks );

					#$mzxmlscan->print_dta;
					my $precursor_mz        = $mzxmlscan->precursor_mz;
					my $precursor_intensity = $mzxmlscan->precursorIntensity;

					my $scanid = $precursor_mz . '::' . $precursor_intensity;
					my $retentiontime = $mzxmlscan->Tr;
					my $mz            = $mzxmlscan->FT_mz;

					#	push @{ $scanmzhash{ int($mz) }->{$scanid} }, $mzxmlscan;
					$scan_pair_hash{$scanid}
					  ->{ $mzxmlscan->fractionationtype } = $mzxmlscan;
				}
			}
		}

		$mzXMLtree->delete;
		$i++;
	}

	foreach my $scanpairid ( keys %scan_pair_hash ) {
		unless ( defined( $scan_pair_hash{$scanpairid}->{'cid'} )
			&& defined( $scan_pair_hash{$scanpairid}->{'etd'} ) )
		{
			warn "no pair of scans is given for $scanpairid $!";
			next;
		}

		my $ETD_CID_pair_obj =
		  XMM_ETD_CID->new( $scan_pair_hash{$scanpairid}->{'cid'},
			$scan_pair_hash{$scanpairid}->{'etd'}, $scanpairid );

#		my $cidscan=$ETD_CID_pair_obj->cid_scan;
#		my $etdscan=$ETD_CID_pair_obj->etd_scan;
#
#		print "$scanpairid: cid: ",$cidscan->scannumber,"\tetd: ",$etdscan->scannumber,"\n";
		my $precursor_mz = $ETD_CID_pair_obj->precursor_mz;
		push @{ $scanmzhash{ int($precursor_mz) } }, $ETD_CID_pair_obj;
	}

	$self->{'n_reference_scans'} = $nscans;
	print
"parsed $nscans MS scans from ETD experiment for levels @$levels; instrument ids @$msInstrumentIDs\n";

	$self->{'reference_scanmzhash'} = \%scanmzhash;
}


#### Emulates the parsing of a mzXML
#### in the pSH workflow the scans numbers are in the Mastermap stored
sub parse_pseudo_mzXML{
## used to emulate a XML for the pseudoSH workflow
my $self            = shift;
my $mzXMLfiles      = shift;
my $readpeaks       = 0;
my $PARAMS          = $self->params;	

my $levels         	= $PARAMS->scanlevels;
my $msInstrumentIDs = $PARAMS->msInstrumentIDs;
my $offset          = $PARAMS->Tr_offset;
my $Trtolerance     = $PARAMS->deltaTr;

my %scanTrhash;

my $nscans = 0;
my $i      = 0;
my $mzXMLfilename;
	
if ( scalar @$mzXMLfiles == 0 ) {
	die "no mzXML file name was indicated \n";
}

### Get all Features
my $features = $self->getallfeatures;

foreach my $mzXMLfile (@$mzXMLfiles) {

( $mzXMLfilename = basename($mzXMLfile) ) =~ s/\.mzXML//;
print "parsing $mzXMLfilename...\n";

## EVERY FEATURE is annotated with an pseudoscan
## Append to every feature the mzXMLscan object
	foreach my $scan (@$features) {
	$nscans++;
	my $mzxmlscan = mzXMLscan->new_pseudoscan( $scan, $i, $mzXMLfilename, $readpeaks );
	## annotate the scanobj
	$scan->annotate_with($mzxmlscan);
	## retentiontimehash
	my $retentiontime = $mzxmlscan->Tr;
	
		if ( $retentiontime > $offset ) {
		push @{ $scanTrhash{ int($retentiontime) } }, $mzxmlscan;
		}
		else {
		push @{ $scanTrhash{ int( $offset - $Trtolerance + 1 ) } },$mzxmlscan;
		}
	
	}
print "READ $nscans Pseudofeatures as Pseudoscans.\n";
}

}


sub parse_mzXML {
	my $self            = shift;
	my $mzXMLfiles      = shift;
	my $readpeaks       = shift;
	my $PARAMS          = $self->params;
	my $levels          = $PARAMS->scanlevels;
	my $msInstrumentIDs = $PARAMS->msInstrumentIDs;
	my $offset          = $PARAMS->Tr_offset;
	my $Trtolerance     = $PARAMS->deltaTr;

	my %scanTrhash;

	#my %fragmentationtypehash=();
	my $nscans = 0;
	my $i      = 0;
	my $mzXMLfilename;
	if ( scalar @$mzXMLfiles == 0 ) {
		die "no mzXML file was indicated \n";
	}
	foreach my $mzXMLfile (@$mzXMLfiles) {

		my $mzXMLtree = XML::TreeBuilder->new();
		$mzXMLtree->parse_file($mzXMLfile);
		( $mzXMLfilename = basename($mzXMLfile) ) =~ s/\.mzXML//;
		print "parsing $mzXMLfilename...\n";
		foreach my $scan ( $mzXMLtree->find_by_tag_name('scan') ) {
			if (
				number_is_in_array( $levels, $scan->attr_get_i('msLevel') )

			  )
			{
				$nscans++;
				my $mzxmlscan = mzXMLscan->new( $scan, $i, $mzXMLfilename, $readpeaks );

			 #				print $mzxmlscan->id, "\n";
			 #				my $normalizedpeaks = $mzxmlscan->get_binned_normalized_peaks;
			 #				foreach my $peakpairs (@$normalizedpeaks) {
			 #					if ( $peakpairs->[1] ) {
			 #						print $peakpairs->[0], "\t", $peakpairs->[1], "\n";
			 #					}
			 #				}

#			push @{$fragmentationtypehash{$mzxmlscan->fractionationtype}},$mzxmlscan->scannumber;
				my $retentiontime = $mzxmlscan->Tr;
				if ( $retentiontime > $offset ) {
					push @{ $scanTrhash{ int($retentiontime) } }, $mzxmlscan;
				}
				else {
					push @{ $scanTrhash{ int( $offset - $Trtolerance + 1 ) } },
					  $mzxmlscan;
				}
			}
		}
		$mzXMLtree->delete;
		$i++;
	}
	$self->{'nscans'} = $nscans;
	print
"parsed $nscans MS scans for levels @$levels; instrument ids @$msInstrumentIDs; fragmentationtype\n";

#foreach my $scantype(keys %fragmentationtypehash){
#print "number of $scantype scans: ",scalar(@{$fragmentationtypehash{$scantype}}),"\n";
#}
#print "\n";
	$self->{'scanTrhash'} = \%scanTrhash;
}

sub get_nscans {
	my $self = shift;
	return $self->{'nscans'};
}

sub get_isopairs_with_two_scans {
	my $self         = shift;
	my $all_isopairs = $self->get_isotopicpairs;
	my @pairswithscan;
	foreach my $isopair (@$all_isopairs) {
		if ( $isopair->light_feature->scan && $isopair->heavy_feature->scan ) {
			push @pairswithscan, $isopair;
		}
	}
	print "found ", scalar(@pairswithscan),
	  " isotopic pairs with scans for light and heavy partner\n";
	return \@pairswithscan;
}

sub get_isopairs_with_any_scans {
	my $self         = shift;
	my $all_isopairs = $self->get_isotopicpairs;
	my @pairswithscan;
	foreach my $isopair (@$all_isopairs) {
		if ( $isopair->light_feature->scan || $isopair->heavy_feature->scan ) {
			push @pairswithscan, $isopair;
		}
	}
	print "found ", scalar(@pairswithscan),
	  " isotopic pairs with scans for light or heavy partner\n";
	return \@pairswithscan;
}

sub get_isopairs_features_with_scans {
	my $self             = shift;
	my $all_isopairs     = $self->get_isotopicpairs;
	my @featureswithscan = ();
	foreach my $isopair (@$all_isopairs) {
		if ( $isopair->light_feature->scan ) {
			push @featureswithscan, $isopair->light_feature;
		}
		if ( $isopair->heavy_feature->scan ) {
			push @featureswithscan, $isopair->heavy_feature;
		}
	}
	print "found ", scalar(@featureswithscan),
	  " features from isotopic pairs with scans for light or heavy partner\n";
	return \@featureswithscan;
}

sub get_isopairs_with_etd_cid_pair {
	my $self         = shift;
	my $all_isopairs = $self->get_isotopicpairs;
	my @pairswithscan;
	foreach my $isopair (@$all_isopairs) {
		if (   $isopair->light_feature->etd_cid_pair
			&& $isopair->heavy_feature->etd_cid_pair )
		{
			push @pairswithscan, $isopair;
		}
	}
	print "found ", scalar(@pairswithscan),
	  " isotopic pairs with scans associated etd and cid scans\n";
	return \@pairswithscan;
}

sub get_isopairs_with_light_etd_cid_pair {
	my $self         = shift;
	my $all_isopairs = $self->get_isotopicpairs;
	my @pairswithscan;
	foreach my $isopair (@$all_isopairs) {
		if ( $isopair->light_feature->etd_cid_pair ) {
			push @pairswithscan, $isopair;
		}
	}
	print "found ", scalar(@pairswithscan),
" isotopic pairs with scans associated from etd scans to light precursor\n";
	return \@pairswithscan;
}

sub get_isopairs_with_heavy_etd_cid_pair {
	my $self         = shift;
	my $all_isopairs = $self->get_isotopicpairs;
	my @pairswithscan;
	foreach my $isopair (@$all_isopairs) {
		if ( $isopair->heavy_feature->etd_cid_pair ) {
			push @pairswithscan, $isopair;
		}
	}
	print "found ", scalar(@pairswithscan),
" isotopic pairs with scans associated from etd scans to heavy precursor\n";
	return \@pairswithscan;
}

sub getisopair_features {
	my $self             = shift;
	my $all_isopairs     = $self->get_isotopicpairs;
	my @selectedfeatures = ();

	if ( $self->params->inclusionlisttype eq "inclusion" ) {
		foreach my $isopair (@$all_isopairs) {
			if ( !$isopair->light_feature->scan ) {
				push @selectedfeatures, $isopair->light_feature;
			}
			if ( !$isopair->heavy_feature->scan ) {
				push @selectedfeatures, $isopair->heavy_feature;
			}

		}
	}

	if ( $self->params->inclusionlisttype eq "exclusion" ) {
		foreach my $isopair (@$all_isopairs) {
			if ( $isopair->light_feature->scan ) {
				push @selectedfeatures, $isopair->light_feature;
			}
			if ( $isopair->heavy_feature->scan ) {
				push @selectedfeatures, $isopair->heavy_feature;
			}

		}
	}

	elsif ( $self->params->inclusionlisttype eq "all" ) {
		foreach my $isopair (@$all_isopairs) {
			push @selectedfeatures, $isopair->light_feature;
			push @selectedfeatures, $isopair->heavy_feature;
		}
	}

	print "put ", scalar(@selectedfeatures),
	  " isotopic pair features into inclusionlist\n";
	return \@selectedfeatures;
}

sub get_isopairs_with_light_scan {
	my $self         = shift;
	my $all_isopairs = $self->get_isotopicpairs;
	my @pairswithscan;
	foreach my $isopair (@$all_isopairs) {
		if ( $isopair->light_feature->scan ) {
			push @pairswithscan, $isopair;
		}
	}
	print "found ", scalar(@pairswithscan),
	  " isotopic pairs with scans for light partner \n";
	return \@pairswithscan;
}

sub get_isopairs_with_heavy_scan {
	my $self         = shift;
	my $all_isopairs = $self->get_isotopicpairs;
	my @pairswithscan;
	foreach my $isopair (@$all_isopairs) {
		if ( $isopair->heavy_feature->scan ) {
			push @pairswithscan, $isopair;
		}
	}
	print "found ", scalar(@pairswithscan),
	  " isotopic pairs with scans for heavy partner \n";
	return \@pairswithscan;
}

sub getTrhash {
	my $self = shift;
	return $self->{'scanTrhash'};
}

sub get_reference_scan_mz_hash {
	my $self = shift;
	return $self->{'reference_scanmzhash'};
}

sub number_is_in_array {
	my $array  = shift;
	my $number = shift;
	foreach my $arrayentry (@$array) {
		if ( $number == $arrayentry ) {
			return 1;
		}
	}
	return 0;
}

sub printmastermap {
	my $self   = shift;
	my $output = shift;

	my $features       = $self->getallfeatures;
	my $inclusionlist  = $self->params->inclusion_list;
	my @sortedfeatures = sort { $a->mz <=> $b->mz } @$features;

	print $output "m/z\tTr\tcharge\tAC\tSQ\tPep.Prob\t#Rep.Match\tAREA_",
	  ( join "\tAREA_", @{ $self->getdefinedIDs } ),

"\tisotopic-shift\tisotopic-pair index\tmatched scan\tmatched mz\tmatched charge\tmatched Tr\tmatched bp intensity\tmatched runID\n";

	foreach my $feature (@sortedfeatures) {
		print $output $feature->mz, "\t", $feature->original_Tr, "\t",
		  $feature->charge, "\t", $feature->accession, "\t", $feature->seq,
		  "\t", $feature->p, "\t", $feature->nreplicates, "\t",
		  ( join "\t", @{ $feature->areaarray } ), "\t", $feature->is_isopair,
		  "\t", $feature->get_isopair_index;
		if ($inclusionlist) {
			my $scan = $feature->associated_scan;
			if ($scan) {
				print $output "\t", $scan->scannumber, "\t", $scan->FT_mz, "\t",
				  $scan->FT_charge, "\t", $scan->Tr / 60, "\t",
				  $scan->basePeakIntensity, "\t", $scan->runID, "\n";
			}
			else {
				print $output "\n";
			}
		}
		else {
			print $output "\n";
		}
	}

}


## assigning of scans
## by searching
sub assignscans {
	my $self     = shift;
	my $features = $self->getallfeatures;
	my $scanhash = $self->getTrhash;

	my $Troffset    = $self->params->Tr_offset;
	my $Trtolerance = $self->params->deltaTr;
	my $Trshift     = $self->params->Trshift;

	my $next2apex = undef;
	if ( $self->params->matchtype eq 'Tr_distance_2apex' ) {
		$next2apex = 1;
	}
	my $nassignedscans = 0;
	my %assignedscans  = ();
	my @deltastats     = ();
	print "assigning peaks to scans. offset: $Troffset Tr tolerance: $Trtolerance Trshift: $Trshift\n";
  
  Feature: foreach my $feature (@$features) {
		my @matchingscans = ();
		my $featuremz     = $feature->mz;
		my $featureTr     = $feature->Tr + $Trshift;
		
		if ( $featureTr <= $Troffset ) {
			$featureTr = $Troffset;
		}
		my $featurecharge = $feature->charge;

		#print "feature: $featureTr\n";
		for my $testTr ( 0 .. int($Trtolerance) ) {

			my @candidatescans = ();

			if ( $scanhash->{ int( $featureTr + $testTr ) } ) {
				push @candidatescans, @{ $scanhash->{ $featureTr + $testTr } };
			}
			if ( $scanhash->{ int( $featureTr - $testTr ) } ) {
				push @candidatescans, @{ $scanhash->{ $featureTr - $testTr } };
			}
			
			foreach my $scan (@candidatescans) {
				if ( $feature->matches($scan) ) {
					unless ( $assignedscans{ $scan->id }++ ) {
						$nassignedscans++;
						unless ( $Troffset == $featureTr ) {

#	print "featureTr: $featureTr ","scanTr: ",$scan->Tr," delta: ",$featureTr - $scan->Tr,"\n";
							push @deltastats, $featureTr - $Trshift - $scan->Tr;
						}
					}
					push @matchingscans, $scan;
					#print $scan;
				}
			}
		}
		
		if ( scalar(@matchingscans) == 1 ) {
			$feature->annotate_with( $matchingscans[0] );
		}
		elsif ( scalar(@matchingscans) > 1 ) {
			$feature->annotate_scan_by_defined_method(
				$self->params->annotationmethod, \@matchingscans );
		}
	}
	print "$nassignedscans out of ", $self->get_nscans,
	  " scans can be assigned to peaks (",
	sprintf( "%.2f", 100 * $nassignedscans / $self->get_nscans ), "%)\n";
	print "mean Tr difference: ",
	  sprintf( "%.2f", XMM_Statistics::mean( \@deltastats ) ),  " stdev: ",
	  sprintf( "%.2f", XMM_Statistics::stdev( \@deltastats ) ), "\n";
}

sub getproteingroups {
	my $self = shift;
	return $self->{'proteingroups'};
}

sub getproteingrouphash {
	my $self = shift;
	return $self->{'proteingrouphash'};
}

sub parsetranslationtable {
	my $self             = shift;
	my $translationtable = shift;
	my %transhash        = ();
	open TRANSLATION, "<$translationtable" or die $!;

	print "reading translationtable table $translationtable ...";
	while (<TRANSLATION>) {
		chomp;
		my @tmp = split /::/;

		# $verbose && print "@tmp\n";
		$transhash{ $tmp[0] } = $tmp[1];
	}
	print "... done \n";
	close(TRANSLATION);

	return \%transhash;
}

sub translationtable {
	my $self = shift;
	return $self->{'translation'};
}

sub print_profile{
	my $self          = shift;
	my $profileout    = shift;
	my $proteingroups = $self->getproteingroups;
		foreach my $protein (@$proteingroups) {
			$protein->print_profile($profileout);
		}
}

sub calcprofiles {
	my $self          = shift;
#	my $profileout    = shift;
	my $runIDs        = $self->getdefinedIDs;
	my $verbose       = $self->verbose;
	my $proteingroups = $self->getproteingroups;

	#my @runids = split /,/, $runIDs;

	if ( scalar(@$runIDs) == 2 ) {
		foreach my $protein (@$proteingroups) {
			$verbose && print "protein accession ", $protein->accession, "\n";
			$protein->calcratiostats($runIDs);
		#	$protein->print_profile($profileout);
			$verbose && $protein->printstats;
		}
	}
	elsif ( scalar(@$runIDs) > 2 ) {
		foreach my $protein (@$proteingroups) {
			$verbose && print "protein accession ", $protein->accession, "\n";
			$protein->calcprofilestats($runIDs);
		#	$protein->print_profile($profileout);
			$verbose && $protein->printstats;
		}
	}

	else {
		die "number of runIDs is < 2 ";
	}
}

sub quantify_protgroups_by_nhighest {
	my $self          = shift;
	my $definedids    = $self->getdefinedIDs;
	my $proteingroups = $self->getproteingroups;
	foreach my $proteingroup (@$proteingroups) {
		$proteingroup->quantify_by_nhighest($definedids);
	}
}

sub print_nhighest_proteinquantifications {
	my $self                    = shift;
	my $outfilehandle           = shift;
	my $pepinfofilehandle       = shift;
	my $defined_ids             = $self->getdefinedIDs;
	my $defined_ids_translation = $self->get_definedIDs_translation;

	( $outfilehandle     = *STDOUT ) unless $outfilehandle;
	( $pepinfofilehandle = *STDOUT ) unless $pepinfofilehandle;

	print $outfilehandle "protein\t",
	  ( join "\t", map { $defined_ids_translation->{$_} } @$defined_ids ), "\n";

	my $proteingroups = $self->getproteingroups;
	foreach my $proteingroup (@$proteingroups) {
		my $runID_areas = $proteingroup->get_nighest_quant_hash;

		print $outfilehandle $proteingroup->accession,     "\t";
		print $pepinfofilehandle $proteingroup->accession, "\t";

		foreach my $id (@$defined_ids) {
			print $outfilehandle $runID_areas->{$id}->{'area'}, "\t";
			my @features = @{ $runID_areas->{$id}->{'features'} };
			print $pepinfofilehandle ( join "::",
				map { $_->clean_seq } @features ), "\t";

		}
		print $outfilehandle "\n";
		print $pepinfofilehandle "\n";
	}

}

sub verbose {
	my $self = shift;
	return $self->{'verbose'};
}

sub printratiofile {
	my $self    = shift;
	my $outfile = shift;
	open OUTFILE, ">$outfile" or die $!;

	print OUTFILE
"Accession\t# features\t# unique features\t# valid ratios\tratio mean\tvariation factor\tratio stdev\tratio sem\tratio median\n";
	my $proteingroups = $self->getproteingroups;

	foreach my $protein (@$proteingroups) {
		print OUTFILE $protein->accession, "\t", $protein->nfeatures, "\t",
		  $protein->nuniquefeatures, "\t", $protein->nratios,      "\t",
		  $protein->ratiomean,       "\t", $protein->ratiopercent, "\t",
		  $protein->ratiostdev,      "\t", $protein->ratiosem,     "\t",
		  $protein->ratiomedian,     "\n";
	}
}

sub printprofilefile {
	my $self    = shift;
	my $outfile = shift;

	print $outfile
	  "Accession\t# features\t# unique features\t# valid profiles\t",
	  ( join $self->params->averaging."\t",  $self->getrunIDarray ), $self->params->averaging,"\t",
	  ( join "stdev\t", $self->getrunIDarray ), "stdev\t",
	  ( join "sem\t",   $self->getrunIDarray ), "sem\tprofilescore\tpeptides\n";
	my $proteingroups = $self->getproteingroups;

	my $averagingmethod = undef;
	if ( $self->params->averaging eq "mean" ) {
		$averagingmethod = \&proteingroups::meanprofile;
	}
	elsif ( $self->params->averaging eq "median" ) {
		$averagingmethod = \&proteingroups::medianprofile;
	}
	else {
		die "averaging method ", $self->params->averaging,
		  " is not implemented, use \"mean\" or \"median\" $!";
	}
	foreach my $protein (@$proteingroups) {
			my $featurestring = join ",", map {$_->seq.'_'.$_->charge} @{$protein->getfeatures_of_protgroup($self->params->average_nmostintense)};
		print $outfile $protein->accession, "\t", $protein->nfeatures, "\t",

	 		  $protein->nuniquefeatures, "\t", $protein->nprofiles, "\t", (join "\t",
	 		  @{ $protein->$averagingmethod },
	 		  ( join "\t", @{ $protein->profilestdev } ),
	 		  ( join "\t", @{ $protein->profilesem } )),
	 		  "\t",$protein->profilescore,"\t",$featurestring,
		  "\n";
	}
}

sub searchtripleshifts {
	my $self           = shift;
	my $PARAMS         = $self->params;
	my $deltaarray     = $self->params->triple; ## an array with the shift

	my $intensityratio = $self->params->pairratio;
	my $deltaTr        = $self->params->Triplepair_Tr_tolerance;

	my $MS1_Mr_tolerance_abs      = $self->params->Triplepair_Mr_tolerance;
	my $tolerancemeasure          = $self->params->Isopair_Mr_tolerance_measure;
	my $lightfeatureids           = $self->params->lightrunids;
	my $heavyfeatureids           = $self->params->heavyrunids;
	my $Isopair_require_same_lcid = $self->params->Isopair_require_same_lcid;

	my $features = $self->getallfeatures;
	my @sortedfeatures = sort { $a->Mr <=> $b->Mr } @$features;
	my @isopairs       = ();
	my $markedpairs	   = {};
	my $isopairindex   = 0;

#---------------------------------------------------------------------------
#  Search for isotopic pairs
#---------------------------------------------------------------------------	
	
	foreach my $delta (@$deltaarray) {
		print "Searching for features with $delta amu mass-shift.\n";
		print "Maximum allowed retention time shift for isotopic pairs is $deltaTr seconds.\n";
		my $j         = 0;
		my $i         = 0;
		my $lastindex = 0;

		for $i ( 0 .. $#sortedfeatures ) {
			my $MS1_Mr_tolerance = calcMS1tolerance( $MS1_Mr_tolerance_abs, $tolerancemeasure,$sortedfeatures[$i]->Mr );
			$j = $lastindex;
			#print "Checking Feature: $i\n";
			## EVAL aslong as the condition is false, aslong as it is a potential candidate for a pair
			until (	( ( $sortedfeatures[$j]->Mr - $sortedfeatures[$i]->Mr ) > ( $delta + $MS1_Mr_tolerance )) || $j >= $#sortedfeatures )
			{
				
				if ( abs( $sortedfeatures[$i]->Mr - $sortedfeatures[$j]->Mr + $delta ) <= $MS1_Mr_tolerance	
				&& abs( $sortedfeatures[$i]->Tr - $sortedfeatures[$j]->Tr )< $deltaTr 
				&& ( $sortedfeatures[$i]->charge == $sortedfeatures[$j]->charge )
					&& (
						$intensityratio < intensityratio(
							$sortedfeatures[$i]->totalarea,
							$sortedfeatures[$j]->totalarea
						)
					)
					&& contains( $sortedfeatures[$i]->lcid, $lightfeatureids )
					&& contains( $sortedfeatures[$j]->lcid, $heavyfeatureids )
					&& (
						!$Isopair_require_same_lcid
						|| overlaps(
							$sortedfeatures[$i]->lcid,
							$sortedfeatures[$j]->lcid
						)
					)
				  )
				{

					
					#push @isopairs, Isopairs->new( $sortedfeatures[$i], $sortedfeatures[$j],$PARAMS, $delta );
					$isopairindex++;
					# Mark as a heavy feature of a triple pair
					#print "Found Pair\n";
					# Mark the FEATURES THAT GIVE THE TRIPLE PAIR 
					$markedpairs->{$j}=1;
					$markedpairs->{$i}=1;
					#$sortedfeatures[$i] ->mark_as_isopair( -$delta, $isopairindex );
					#$sortedfeatures[$j] ->mark_as_isopair( $delta, $isopairindex );

				}
				if ( ( $sortedfeatures[$j]->Mr - $sortedfeatures[$i]->Mr ) <
					( $delta - $MS1_Mr_tolerance ) )
				{
					$lastindex++;
				}
				$j++;
			}
		}
	}
	#print Dumper ($markedpairs);
	my @keys = keys %$markedpairs;
	#exit;
	print "Found and marked ".scalar (@keys)." isotopic triple pairs\n";
	#exit;
	$self->{'triplepairs'} = $markedpairs;
	#print scalar(@isopairs);
}



### search for isotopic pairs
sub searchdeltashifts {
	my $self           = shift;
	my $PARAMS         = $self->params;
	my $deltaarray     = $self->params->delta; ## an array with the shift

	my $intensityratio = $self->params->pairratio;
	my $deltaTr        = $self->params->Isopair_Tr_tolerance;

	my $MS1_Mr_tolerance_abs      = $self->params->Isopair_Mr_tolerance;
	my $tolerancemeasure          = $self->params->Isopair_Mr_tolerance_measure;
	my $lightfeatureids           = $self->params->lightrunids;
	my $heavyfeatureids           = $self->params->heavyrunids;
	my $Isopair_require_same_lcid = $self->params->Isopair_require_same_lcid;

	my $features = $self->getallfeatures;

	my @sortedfeatures = sort { $a->Mr <=> $b->Mr } @$features;
	my @isopairs       = ();
	my $isopairindex   = 0;

	my $markedpairs = $self->{'triplepairs'};
	my $skipped = 0;
	#print Dumper ($markedpairs);
	#exit;
#---------------------------------------------------------------------------
#  Search for isotopic pairs
#---------------------------------------------------------------------------	
	
	foreach my $delta (@$deltaarray) {
		print "Searching for features with $delta amu mass-shift.\n";
		print "Maximum allowed retention time shift for isotopic pairs is $deltaTr seconds.\n";
		my $j         = 0;
		my $i         = 0;
		my $lastindex = 0;

		for $i ( 0 .. $#sortedfeatures ) {
			my $MS1_Mr_tolerance = calcMS1tolerance( $MS1_Mr_tolerance_abs, $tolerancemeasure,$sortedfeatures[$i]->Mr );
			$j = $lastindex;
			#print "Checking Feature: $i\n";
			## EVAL aslong as the condition is false, aslong as it is a potential candidate for a pair
			until (	( ( $sortedfeatures[$j]->Mr - $sortedfeatures[$i]->Mr ) > ( $delta + $MS1_Mr_tolerance )) || $j >= $#sortedfeatures )
			{
				
				if ( abs( $sortedfeatures[$i]->Mr - $sortedfeatures[$j]->Mr + $delta ) <= $MS1_Mr_tolerance	
				&& abs( $sortedfeatures[$i]->Tr - $sortedfeatures[$j]->Tr )< $deltaTr 
				&& ( $sortedfeatures[$i]->charge == $sortedfeatures[$j]->charge )
					&& (
						$intensityratio < intensityratio(
							$sortedfeatures[$i]->totalarea,
							$sortedfeatures[$j]->totalarea
						)
					)
					&& contains( $sortedfeatures[$i]->lcid, $lightfeatureids )
					&& contains( $sortedfeatures[$j]->lcid, $heavyfeatureids )
					&& (
						!$Isopair_require_same_lcid
						|| overlaps(
							$sortedfeatures[$i]->lcid,
							$sortedfeatures[$j]->lcid
						)
					)
				  )
				{

					$isopairindex++;
					
					if ($markedpairs->{$j} || $markedpairs->{$i} ){
					#print "Marked pair found, will skip this triple pair\n";
					$skipped++;
					}else{					
					push @isopairs, Isopairs->new( $sortedfeatures[$i], $sortedfeatures[$j],$PARAMS, $delta );
					$sortedfeatures[$i] ->mark_as_isopair( -$delta, $isopairindex );
					$sortedfeatures[$j] ->mark_as_isopair( $delta, $isopairindex );
					}
				}
				if ( ( $sortedfeatures[$j]->Mr - $sortedfeatures[$i]->Mr ) < ( $delta - $MS1_Mr_tolerance ) )
				{
					$lastindex++;
				}
				$j++;
			}
		}
	}
	$self->{'isotopicpairs'} = \@isopairs;
	print "Skipped: ". $skipped. " isotopic triplepairs\n";
	#print scalar(@isopairs);
}


sub overlaps {
	my $array1 = shift;
	my $array2 = shift;
	my $i;
	my $sum = 0;
	for $i ( 0 .. $#$array1 ) {
		$sum += $array1->[$i] * $array2->[$i];
	}
	return $sum;
}

sub normalize_MasterMap {
	my $self = shift;

	#my $stdev_vector   = $self->getMasterMap_stdev_vector;
	#my $average_vector = $self->getMasterMap_average_vector;

	#my $average_vector = $self->getMasterMap_median_vector;
	my $average_vector = $self->getMasterMap_topn_vector;
	

	my $features = $self->getallfeatures;
	my $i;
	foreach my $feature (@$features) {
		my $lcidvector = $feature->areaarray;
		for $i ( 0 .. $#$lcidvector ) {

			#			$lcidvector->[$i] =
			#			  ( $lcidvector->[$i] - $average_vector->[$i] ) /
			#			  $stdev_vector->[$i];
			$lcidvector->[$i] /= $average_vector->[$i];
		}
	}
}

sub getMasterMap_stdev_vector {
	my $self     = shift;
	my $features = $self->getallfeatures;
	my $i;
	my %lcids = ();
	foreach my $feature (@$features) {
		my $lcidvector = $feature->areaarray;
		for $i ( 0 .. $#$lcidvector ) {
			if ( $lcidvector->[$i] ) {
				push @{ $lcids{$i} }, $lcidvector->[$i];
			}
		}
	}
	my @stdevarray = ();
	foreach my $lcid ( sort { $a <=> $b } keys %lcids ) {
		my $stdev = XMM_Statistics::stdev( $lcids{$lcid} );
		push @stdevarray, $stdev;
		print "stev lcid $lcid = $stdev\n";
	}

	return \@stdevarray;
}

sub calcMS1tolerance {
	my $MS1_Mr_tolerance_abs = shift;
	my $tolerancemeasure     = shift;
	my $Mr                   = shift;

	if ( $tolerancemeasure eq "amu" ) {
		return $MS1_Mr_tolerance_abs;
	}
	elsif ( $tolerancemeasure eq "ppm" ) {
		return $MS1_Mr_tolerance_abs * 1e-6 * $Mr;
	}
	else {
		die "tolerance measure $tolerancemeasure is not defined $!";
	}

}


sub getMasterMap_topn_vector{

	my $self     = shift;
	my $features = $self->getallfeatures;
	my $i;
	my %lcids = ();
	foreach my $feature (@$features) {
		my $lcidvector = $feature->areaarray;
		for $i ( 0 .. $#$lcidvector ) {
			if ( $lcidvector->[$i] ) {
				push @{ $lcids{$i} }, $lcidvector->[$i];
			}
		}
	}
	my @meanarray = ();
	foreach my $lcid ( sort { $a <=> $b } keys %lcids ) {
		my $topten = XMM_Statistics::topn( $lcids{$lcid} );
		
	#	my $mean = XMM_Statistics::mean( $lcids{$lcid} );
		push @meanarray, $topten;
		print "average intensity lcid $lcid = $topten\n";
	}

	return \@meanarray;
}
	
sub getMasterMap_average_vector {
	my $self     = shift;
	my $features = $self->getallfeatures;
	my $i;
	my %lcids = ();
	foreach my $feature (@$features) {
		my $lcidvector = $feature->areaarray;
		for $i ( 0 .. $#$lcidvector ) {
			if ( $lcidvector->[$i] ) {
				push @{ $lcids{$i} }, $lcidvector->[$i];
			}
		}
	}
	my @meanarray = ();
	foreach my $lcid ( sort { $a <=> $b } keys %lcids ) {
		my $mean = XMM_Statistics::mean( $lcids{$lcid} );
		push @meanarray, $mean;
		print "average intensity lcid $lcid = $mean\n";
	}

	return \@meanarray;
}



sub getMasterMap_median_vector {
	my $self     = shift;
	my $features = $self->getallfeatures;
	my $i;
	my %lcids = ();
	foreach my $feature (@$features) {
		my $lcidvector = $feature->areaarray;
		for $i ( 0 .. $#$lcidvector ) {
			if ( $lcidvector->[$i] ) {
				push @{ $lcids{$i} }, $lcidvector->[$i];
			}
		}
	}
	my @meanarray = ();
	foreach my $lcid ( sort { $a <=> $b } keys %lcids ) {
		my $mean = XMM_Statistics::median( $lcids{$lcid} );
		push @meanarray, $mean;
		print "median intensity lcid $lcid = $mean\n";
	}

	return \@meanarray;
}

#sub delta_feature_stats {
#	my $self          = shift;
#	my $ids           = shift;
#	my $outfilehandle = shift;
#	unless ($outfilehandle) {
#		$outfilehandle = *STDOUT;
#	}
#
#	my $isotopicpairs = $self->get_isotopicpairs;
#	my %deltahash     = ();
#	foreach my $pair (@$isotopicpairs) {
#		$deltahash{ $pair->delta }++;
#	}
#
#	foreach my $delta ( sort { $a <=> $b } keys %deltahash ) {
#		print "found ", $deltahash{$delta}, " features with mass shift of ",
#		  $delta, " amu\n";
#	}
#	my $lightfeatureids = $self->params->lightrunids;
#	my $heavyfeatureids = $self->params->heavyrunids;
#
#	my $features = $self->getallfeatures;
#
#	my @featureintensities = ();
#	foreach my $feature (@$features) {
#		if ( $feature->ispresent($ids) ) {
#			push @featureintensities, $feature->areasum($ids);
#		}
#	}
#	my $topten = XMM_Statistics::topn( \@featureintensities, 10 );
#
#	my %chargedpairs;
#	foreach my $isopair (@$isotopicpairs) {
#		push @{ $chargedpairs{ $isopair->charge } }, $isopair;
#	}
#
#	print $outfilehandle
#"charge\tnpairs\tdelta Tr heavy - light\tsem\tmedian intensity light\tmedian intensity heavy\t light top10 ratio\theavy top10 ratio\n";
#	foreach my $charge ( sort { $a <=> $b } keys %chargedpairs ) {
#		my $npairs         = 0;
#		my @matched        = ();
#		my @Intensitylight = ();
#		my @Intensityheavy = ();
#		my @Trlight        = ();
#		my @Trheavy        = ();
#		my @Trdeltas       = ();
#
#		foreach my $isopair ( @{ $chargedpairs{$charge} } ) {
#			push @Trlight, $isopair->light_feature->Tr;
#			push @Trheavy, $isopair->heavy_feature->Tr;
#
#			push @Intensitylight,
#			  $isopair->light_feature->areasum($lightfeatureids);
#			push @Intensityheavy,
#			  $isopair->heavy_feature->areasum($heavyfeatureids);
#			push @Trdeltas,
#			  $isopair->light_feature->Tr - $isopair->heavy_feature->Tr;
#			$npairs++;
#		}
#
#		my $lightmean = XMM_Statistics::median( \@Intensitylight );
#		my $heavymean = XMM_Statistics::median( \@Intensityheavy );
#
#		print $outfilehandle "$charge\t$npairs\t";
#
#		print $outfilehandle XMM_Statistics::mean( \@Trdeltas ), "\t",
#		  XMM_Statistics::sem( \@Trdeltas ), "\t";
#		print $outfilehandle "$lightmean\t$heavymean\t";
#
#		print $outfilehandle XMM_Statistics::topn( \@Intensitylight, 10 ) /
#		  $topten, "\t";
#		print $outfilehandle XMM_Statistics::topn( \@Intensityheavy, 10 ) /
#		  $topten, "\n";
#	}
#}

sub delta_feature_stats {
	my $self              = shift;
	my $ids               = shift;
	my $outfilehandle     = shift;
	my $deltafeatures_xls = shift;
	
	unless ($outfilehandle) {
		$outfilehandle = *STDOUT;
	}
	unless ($deltafeatures_xls) {
		$deltafeatures_xls = *STDOUT;
	}

	print $deltafeatures_xls
"m/z\tTr_light\tTr_heavy\tdelta_Tr\tcharge\tMr_light\tMr_heavy\tdelta\tscan_light\tscan_heavy\trepeatcount\tAREA_",
	  ( join "\tAREA_", @{ $self->getdefinedIDs } ), "\n";
	
	my $isotopicpairs = $self->get_isotopicpairs;
	my %deltahash     = ();
	
	foreach my $pair (@$isotopicpairs) {
		
		$deltahash{ $pair->delta }++;
		$pair->profile;
		
		$pair->print_deltafeatures_xls($deltafeatures_xls);
		
	}

	foreach my $delta ( sort { $a <=> $b } keys %deltahash ) {
		print "found ", $deltahash{$delta}, " features with mass shift of ",
		  $delta, " amu\n";
	}
	
	my $lightfeatureids = $self->params->lightrunids;
	my $heavyfeatureids = $self->params->heavyrunids;

	my $features = $self->getallfeatures;

	my @featureintensities = ();
	foreach my $feature (@$features) {
		if ( $feature->ispresent($ids) ) {
			push @featureintensities, $feature->areasum($ids);
		}
	}
	my $topten = XMM_Statistics::topn( \@featureintensities, 50 )||1;

	my %chargedpairs;
	foreach my $isopair (@$isotopicpairs) {
		push @{ $chargedpairs{ $isopair->charge } }, $isopair;
	}

	print $outfilehandle
"charge\tnpairs\tnlightscans\tnheavyscans\tnbothscans\tdelta Tr heavy - light\tsem\tmedian intensity light\tmedian intensity heavy\t light top10 ratio\theavy top10 ratio\n";
	foreach my $charge ( sort { $a <=> $b } keys %chargedpairs ) {
		my $npairs              = 0;
		my $npairswithlightscan = 0;
		my $npairswithheavyscan = 0;
		my $npairswithbothscans = 0;

		my @matched        = ();
		my @Intensitylight = ();
		my @Intensityheavy = ();
		my @Trlight        = ();
		my @Trheavy        = ();
		my @Trdeltas       = ();

		foreach my $isopair ( @{ $chargedpairs{$charge} } ) {
			push @Trlight, $isopair->light_feature->Tr;
			push @Trheavy, $isopair->heavy_feature->Tr;

			push @Intensitylight,
			  $isopair->light_feature->areasum($lightfeatureids);
			push @Intensityheavy,
			  $isopair->heavy_feature->areasum($heavyfeatureids);
			push @Trdeltas,
			  $isopair->light_feature->Tr - $isopair->heavy_feature->Tr;
			$npairs++;

			if ( $isopair->light_feature->scan ) {
				$npairswithlightscan++;
			}
			if ( $isopair->heavy_feature->scan ) {
				$npairswithheavyscan++;
			}
			if (   $isopair->light_feature->scan
				&& $isopair->heavy_feature->scan )
			{
				$npairswithbothscans++;
			}

		}

		my $lightmean = XMM_Statistics::median( \@Intensitylight );
		my $heavymean = XMM_Statistics::median( \@Intensityheavy );

		print $outfilehandle
"$charge\t$npairs\t$npairswithlightscan\t$npairswithheavyscan\t$npairswithbothscans\t";

		print $outfilehandle XMM_Statistics::mean( \@Trdeltas ), "\t",
		  XMM_Statistics::sem( \@Trdeltas ), "\t";
		print $outfilehandle "$lightmean\t$heavymean\t";

		#print $outfilehandle XMM_Statistics::topn( \@Intensitylight, 50 ) /
		#  $topten||1, "\t";
		#print $outfilehandle XMM_Statistics::topn( \@Intensityheavy, 50 ) /
		#  $topten||1, "\n";
	}
}

sub contains {
	my $lcid     = shift;
	my $validids = shift;
	my %seen     = ();
	my $i;
	foreach my $validid (@$validids) {
		if ( $lcid->[$validid] ) {
			return 1;
		}
	}
	return 0;
}

sub intensityratio {
	my $i1          = shift;
	my $i2          = shift;
	my @intensities = sort { $a <=> $b } ( $i1, $i2 );
	return $intensities[0] / $intensities[1];
}

sub separate_protgroups {
	my $self              = shift;
	my %prothash          = ();
	my %proteingrouphash  = ();
	my $annotatedfeatures = $self->getfeatureswithseq;
	my $minfeaturecount   = $self->params->mincount;
	my @runIDarray        = $self->getrunIDarray;
	my @requiredids       = $self->params->requiredrunIDs;
	foreach my $feature (@$annotatedfeatures) {

		#print $feature->nrepeats(\@runIDarray),"\n";

		if ( $feature->ispresent( \@requiredids ) ) {
			push @{ $prothash{ $feature->accession } }, $feature;
			print $feature->accession, " ", $feature->seq, " ",
			  $feature->clean_seq, "\n";
		}
	}

	my @proteingroups = ();
	foreach my $protein ( keys %prothash ) {
		my $proteingroup =
		  proteingroups->new( $protein, $prothash{$protein}, $self );
		push @proteingroups, $proteingroup;
		$proteingrouphash{$protein} = $proteingroup;
	}
	$self->{'proteingroups'}    = \@proteingroups;
	$self->{'proteingrouphash'} = \%proteingrouphash;

	return \@proteingroups;
}

sub separate_protgroups_by_group_hash {
	my $self      = shift;
	my $grouphash = shift;

	my %prothash          = ();
	my %proteingrouphash  = ();
	my $annotatedfeatures = $self->getfeatureswithseq;
	my $minfeaturecount   = $self->params->mincount;
	my @runIDarray        = $self->getrunIDarray;
	my @requiredids       = $self->params->requiredrunIDs;
	my $strip_peptides    = $self->params->strip_peptides;
	my $verbose=$self->params->verbose;
	foreach my $feature (@$annotatedfeatures) {

		#if ( $feature->ispresent( \@requiredids ) ) {
		my $seq = '';
		if ($strip_peptides) {
			$seq = $feature->clean_seq;
		}
		else {
			$seq = $feature->seq;
		}
		if($grouphash->{$seq}) {
		push @{ $prothash{ $grouphash->{$seq} } },$feature;
		}else{
		$verbose && warn "sequence $seq is not assigned to any protein\n";
		}
#	print $grouphash->{$feature->clean_seq}," ",$feature->accession," ",$feature->seq," ", $feature->clean_seq,"\n";
#}
	}

	my @proteingroups = ();
	foreach my $protein ( keys %prothash ) {
		my $proteingroup =
		  proteingroups->new( $protein, $prothash{$protein}, $self );
		push @proteingroups, $proteingroup;
		$proteingrouphash{$protein} = $proteingroup;
	}
	$self->{'proteingroups'}    = \@proteingroups;
	$self->{'proteingrouphash'} = \%proteingrouphash;

	return \@proteingroups;
}

sub getfeatureswithseq {
	my $self        = shift;
	my $minprob     = $self->params->minprob;
	my $allfeatures = $self->getallfeatures;
	my @featureswithseq;
	foreach my $feature (@$allfeatures) {
		if ( $feature->p > $minprob ) {
			push @featureswithseq, $feature;
		}
	}
	return \@featureswithseq;
}

sub MapStats {
	my $self = shift;
	my $ids  = shift;
	unless ($ids) {
		$ids = $self->getdefinedIDs;
	}
	my $proteingroups = $self->getproteingroups;
	my $features      = $self->getallfeatures;

	my $nseqfeatures    = 0;
	my $nfeatures       = 0;
	my $nuniquefeatures = 0;

	foreach my $protein (@$proteingroups) {
		$nseqfeatures    += $protein->nfeatures;
		$nuniquefeatures += $protein->nuniquefeatures;
	}

	my @featureintensities = ();
	my %charges            = ();
	foreach my $feature (@$features) {
		if ( $feature->ispresent($ids) ) {
			push @featureintensities, $feature->areasum($ids);
		}
		$nfeatures++;
		$charges{ $feature->charge }++;
	}
	my $topten = XMM_Statistics::topn( \@featureintensities );

	foreach my $charge ( sort keys %charges ) {
		print "charge $charge n = ", $charges{$charge}, "\n";
	}
	print
"MasterMap contains $nfeatures features, $nuniquefeatures unique sequences \n";
	print "top ten intensity= $topten\n";
}


sub normalize_profiles_to_target{
	my $self               = shift;
	my $target             = shift;
	
	my $proteingroups      = $self->getproteingroups;
	my $proteingrouphash   = $self->getproteingrouphash;

	my $targetproteingroup = $proteingrouphash->{$target};

	my $averagingmethod = undef;
	if ( $self->params->averaging eq "mean" ) {
		$averagingmethod = \&proteingroups::meanprofile;
	}
	elsif ( $self->params->averaging eq "median" ) {
		$averagingmethod = \&proteingroups::medianprofile;
	}


	my $reference_profile=[];
	unless ( defined($targetproteingroup) ) {
		die
"target protein profile $target was not found. Check name in config file $!";
	}
	else {
		print "scoring against target profile of $target ", join "\t",
		  @{ $targetproteingroup->$averagingmethod }, "\n";
		  
		$reference_profile=$targetproteingroup->$averagingmethod;
	}

	my $targetprofile = $targetproteingroup->$averagingmethod;
	
	foreach my $protein (@$proteingroups) {

		my $profile =$protein->$averagingmethod;
			for my $i (0..$#$profile){
			$profile->[$i]/=$reference_profile->[$i];
			}
				if ( $self->params->renormalize ) {
			Norm::normalize($profile,
				"by_max");
				}
	}


}
sub ScoreProtgroups {
	my $self               = shift;
	my $proteingroups      = $self->getproteingroups;
	my $proteingrouphash   = $self->getproteingrouphash;
	my $target             = $self->params->target;
	my $targetproteingroup = $proteingrouphash->{$target};

	my $averagingmethod = undef;
	if ( $self->params->averaging eq "mean" ) {
		$averagingmethod = \&proteingroups::meanprofile;
	}
	elsif ( $self->params->averaging eq "median" ) {
		$averagingmethod = \&proteingroups::medianprofile;
	}

	unless ( defined($targetproteingroup) ) {
		die
"target protein profile $target was not found. Check name in config file $!";
	}
	else {
		print "scoring against target profile of $target ", join "\t",
		  @{ $targetproteingroup->$averagingmethod }, "\n";
	}

	my $targetprofile = $targetproteingroup->$averagingmethod;
	foreach my $protein (@$proteingroups) {

		my $score =
		  Score::score( $protein->$averagingmethod, $targetprofile,
			$self->params->score );

#		print "scoring ",
#		$protein->accession," ","@{$protein->$averagingmethod}", "target @$targetprofile score: $score\n";

		$protein->setscore($score);
	}
}

1;
