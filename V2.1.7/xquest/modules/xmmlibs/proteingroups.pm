package proteingroups;
use strict;
#---------------------------------------------------------------------------
# Module: proteingroups.pm
# Author(s): Thomas Walzthoeni, xQuest specific modules are based on 
# original work by Oliver Rinner.
# Description: xmm.pl specific module.
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
use XMM_Statistics;
use Score;
use Norm;

sub new {
	my $class            = shift();
	my $proteinaccession = shift;
	my $features         = shift;
	my $MasterMap        = shift;

	my $self = {};
	bless $self, $class;
	$self->{'accession'} = $proteinaccession;
	$self->{'features'}  = $features;
	$self->{'MasterMap'} = $MasterMap;

	return $self;
}

sub accession {
	my $self = shift;
	return $self->{'accession'};
}

sub normfactor {
	my $self = shift;
	return $self->MasterMap->normfactor;
}

sub MasterMap {
	my $self = shift;
	return $self->{'MasterMap'};
}

sub params {
	my $self = shift;
	return $self->MasterMap->params;
}

sub printfeatures {
	my $self     = shift;
	my $features = $self->getfeatures_of_protgroup;
	foreach my $feature (@$features) {
		print $feature->id, "\t", $feature->clean_seq, "\t",
		  ( join "\t", @{ $feature->areaarray } ), "\n";
	}
}

sub getfeatures_of_protgroup {
	my $self         = shift;
	my $nmostintense = shift;
	my $allfeatures  = $self->{'features'};

	unless ($nmostintense) {
		return $allfeatures;
	}
	else {
		my @returnfeatures       = ();
		my @nmostintensefeatures =
		  sort { $b->meanintensity <=> $a->meanintensity } @$allfeatures;
		my $i = 0;
		for (
			$i = 0 ;
			( $i < $nmostintense && $i <= $#nmostintensefeatures ) ;
			$i++
		  )
		{
			push @returnfeatures, $nmostintensefeatures[$i];
		}
		return \@returnfeatures;
	}
}

sub print_profile {
	my $self          = shift;
	my $outfilehandle = shift || *STDOUT;
	my $profile       = $self->meanprofile;
	print $outfilehandle $self->accession, "\t", ( join "\t", @$profile ), "\n";

}

sub nfeatures {
	my $self = shift;
	return scalar( @{ $self->getfeatures_of_protgroup } );
}

sub calcratiostats {
	my $self   = shift;
	my $runIDs = shift;

	my $normfactor = $self->normfactor;

	my $features = $self->getfeatures_of_protgroup;
	my @ratios   = ();
	foreach my $feature (@$features) {
		my $ratio = $feature->ratio($runIDs);

		#print "ratio: ",$runIDs,"\n";
		if ($ratio) {
			push @ratios, $ratio * $normfactor;
		}
	}

	$self->{'ratios'}    = \@ratios;
	$self->{'ratiomean'} = XMM_Statistics::mean( \@ratios );

	if ( $self->{'ratiomean'} > 1 ) {
		$self->{'ratiopercent'} = $self->{'ratiomean'} - 1;
	}
	else {
		if ( $self->{'ratiomean'} > 0 ) {
			$self->{'ratiopercent'} = -( 1 / $self->{'ratiomean'} - 1 );
		}
		else {
			$self->{'ratiopercent'} = 0;
		}
	}
	$self->{'ratiomean'} = XMM_Statistics::mean( \@ratios );

	$self->{'ratiomedian'} = XMM_Statistics::median( \@ratios );
	$self->{'ratiostdev'}  = XMM_Statistics::stdev( \@ratios );
	$self->{'ratiosem'}    = XMM_Statistics::sem( \@ratios );
	$self->{'nratios'}     = XMM_Statistics::n( \@ratios );
}

sub calcprofilestats {
	my $self   = shift;
	my $runIDs = shift;

	#my $normfactor                    = $self->normfactor;
	my $average_nmostintense          = $self->params->average_nmostintense;
	my $normalize_individual_profiles = $self->params->normalize_profiles;

	#my $features = $self->getfeatures_of_protgroup;
	my $features = $self->getfeatures_of_protgroup($average_nmostintense);

	my @profiles = ();
	foreach my $feature (@$features) {
		my $profile = $feature->profile($runIDs);

		#print "ratio: ",$runIDs,"\n";
		if ($profile) {
			if ($normalize_individual_profiles) {
				Norm::normalize( $profile, $normalize_individual_profiles );
			}
			push @profiles, $profile;
		}
	}

	$self->{'profiles'} = \@profiles;

	$self->{'meanprofile'} = XMM_Statistics::mean_array( \@profiles );
	if ( $self->params->renormalize ) {
	#	if ( $self->params->normalize_profiles ) {
			Norm::normalize( $self->{'meanprofile'},
				"by_max");
	#	}
	}

	$self->{'medianprofile'} = XMM_Statistics::median_array( \@profiles );
	if ( $self->params->renormalize ) {
	#	if ( $self->params->normalize_profiles ) {
			Norm::normalize( $self->{'medianprofile'},
				"by_max");
	#	}
	}
	$self->{'profilestdev'} = XMM_Statistics::stdev_array( \@profiles );
	$self->{'profilesem'}   = XMM_Statistics::sem_array( \@profiles );
	$self->{'nprofiles'}    = XMM_Statistics::n( \@profiles );
}

sub quantify_by_nhighest {
	my $self   = shift;
	my $runIDs = shift;

	my $PARAMS = $self->params;

	my $nhighest = $PARAMS->quantify_by_nhighest;
	my $features = $self->getfeatures_of_protgroup;

	my @profiles    = ();
	my %runID_areas = ();
	foreach my $id (@$runIDs) {
		my @areas = sort { $b <=> $a } map { $_->area($id) } @$features;
		my @peptides = sort { $b->area($id) <=> $a->area($id) } @$features;
		my @nhighest_areas = ();
		my @nhighest_peps  = ();

		my $i = 0;
		for ( $i = 0 ; ( $i < $nhighest && $i <= $#areas ) ; $i++ ) {
			push @nhighest_areas, $areas[$i];
			push @nhighest_peps,  $peptides[$i];

		}
		$runID_areas{$id}->{'area'} = XMM_Statistics::mean( \@nhighest_areas );
		$runID_areas{$id}->{'features'} = \@nhighest_peps;

		#		print "@areas\n";
		#		print "@nhighest_areas\n";
	}
	$self->{'nighest_quant'} = \%runID_areas;
}

sub get_nighest_quant_hash {
	my $self = shift;
	return $self->{'nighest_quant'};
}

sub meanprofile {
	my $self = shift;
	return $self->{'meanprofile'};
}

sub medianprofile {
	my $self = shift;
	return $self->{'medianprofile'};
}

sub profilestdev {
	my $self = shift;
	return $self->{'profilestdev'};
}

sub profilesem {
	my $self = shift;
	return $self->{'profilesem'};
}

sub setscore {
	my $self  = shift;
	my $score = shift;
	$self->{'profilescore'} = $score;
}

sub profilescore {
	my $self = shift;
	return $self->{'profilescore'};
}

sub nprofiles {
	my $self = shift;
	return $self->{'nprofiles'};
}

sub getratios {
	my $self = shift;
	return $self->{'ratios'};
}

sub ratiopercent {
	my $self = shift;
	return $self->{'ratiopercent'};
}

sub ratiomean {
	my $self = shift;
	return $self->{'ratiomean'};
}

sub ratiomedian {
	my $self = shift;
	return $self->{'ratiomedian'};
}

sub ratiostdev {
	my $self = shift;
	return $self->{'ratiostdev'};
}

sub ratiosem {
	my $self = shift;
	return $self->{'ratiosem'};
}

sub nratios {
	my $self = shift;
	return $self->{'nratios'};
}

sub nuniquefeatures {
	my $self         = shift;
	my %seen         = ();
	my $nuniquefeats = 0;
	my $features     = $self->getfeatures_of_protgroup;
	foreach my $feature (@$features) {
		unless ( $seen{ $feature->seq }++ ) {
			$nuniquefeats++;
		}
	}
	return $nuniquefeats;
}

sub printstats {
	my $self = shift;
	print "\#features = ",         $self->nfeatures,       "\n";
	print "\#unique sequences = ", $self->nuniquefeatures, "\n";
	print "\# of valid ratios = ", $self->nratios,         "\n";
	print "ratio mean ",   sprintf( "%.3f", $self->ratiomean ),   "\n";
	print "ratio median ", sprintf( "%.3f", $self->ratiomedian ), "\n";
	print "ratio stdev ",  sprintf( "%.3f", $self->ratiostdev ),  "\n";
	print "ratio sem ",    sprintf( "%.3f", $self->ratiosem ),    "\n";
}
1;
