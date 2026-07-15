package inclusion;
use strict;
#---------------------------------------------------------------------------
# Module: inclusion.pm
# Author(s): Thomas Walzthoeni, xQuest specific modules are based on 
# original work by Oliver Rinner.
# Description: Module for handling inclusionlists.
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
sub new {
	my $class     = shift();
	my $MasterMap = shift;
	my $PARAMS    = shift;

	my $self = {};
	bless $self, $class;

	$self->{'MasterMap'} = $MasterMap;
	$self->{'PARAMS'}    = $PARAMS;

	$self->calcbins;

	return $self;
}

sub make_list {
	my $self     = shift;
	my $listfile = shift;

	my $outfilehandle;
	unless ($listfile) {
		$outfilehandle = *STDOUT;
	}
	else {
		open OUT, ">$listfile" or die $!;
		$outfilehandle = *OUT;
	}

	my $allowedcharges=$self->params->IL_charges;
	my $bins_min = $self->bins_min;
	my $bins_max = $self->bins_max;
	my $features;
	my ( $sorted_Tr, $i );
	my $tr_tolerance = $self->tr_tolerance;
	if ( $self->params->delta) {
		$features = $self->MasterMap->getisopair_features;
	}
	elsif ( $self->params->inclusionlisttype eq "inclusion" ) {
		$features = $self->MasterMap->get_features_wo_scan;
	}
	elsif ( $self->params->inclusionlisttype eq "exclusion" ) {
		$features = $self->MasterMap->get_features_with_scan;
	}
	elsif ( $self->params->inclusionlisttype eq "all" ) {
		$features = $self->MasterMap->getallfeatures;
	}
	
	else {
		print "inclusion list type: ", $self->inclusionlisttype,
		  " is not defined";
		die;
	}

	

	if ($allowedcharges){
	print "putting only charges @$allowedcharges in inclusionlist\n";
		my @features=();
		foreach my $feature (@$features){
			if(number_is_in_array($feature->charge,$allowedcharges)){
				push @features,$feature;
			}
		}
	$features=\@features;
	}

	if ( $self->params->nmostintensefeatures ) {
		unless (
			$self->params->nmostintensefeatures->[1] >= scalar(@$features) )
		{
			@$features =
			  ( sort { $a->totalarea <=> $b->totalarea } @$features )
			  [ $self->params->nmostintensefeatures->[0]
			  .. $self->params->nmostintensefeatures->[1] ];
		}
		else {
			warn "there are less than ",
			  $self->params->nmostintensefeatures->[1],
			  "features available, ignoring intensity selection $!";
		}
	}

	print "Tr-tolerance = $tr_tolerance\n";
	my @sorted_Tr = sort { $a->Tr <=> $b->Tr } @$features;

	my $minindex = 0;
  feature: foreach my $feature (@sorted_Tr) {

		#		print $feature->Tr/60," ";
		for $i ( $minindex .. $#$bins_min ) {
			if ( $feature->Tr / 60 <= ( $bins_max->[$i] - $tr_tolerance ) ) {

		 #				print "lower ",$bins_min->[$i]," ","upper: ",$bins_max->[$i],"\n";
				$feature->setlower_bin( $bins_min->[$i] );
				$feature->setupper_bin( $bins_max->[$i] );
				$minindex = $i;
				next feature;
			}
		}
	}
	my %seen_ion_mz;
	my %seen_ion;
	my %nsegments;
	my $nfeatures              = 0;
	my $noriginalfeatures      = 0;
	my $redundant_features     = 0;
	my $halfredundant_features = 0;
	my $modified_features      = 0;

	#printout features with their segments
	my %segmentpopulation;
	foreach my $feature ( sort { $a->Tr <=> $b->Tr } @sorted_Tr ) {
		$noriginalfeatures++;
		unless ( defined( $feature->bin_min ) && defined( $feature->bin_max ) )
		{
			next;
		}
		if ( !$seen_ion_mz{ concatstring($feature) }++ ) {
			if ( $seen_ion{ $feature->mz }++ ) {
				print $outfilehandle sprintf( "%.4f",
					$feature->mz + 0.0001 * ( $seen_ion{ $feature->mz } - 1 ) ),
				  "\t", sprintf( "%.2f", $feature->bin_min ), "\t",
				  sprintf( "%.2f", $feature->bin_max ), "\t", $feature->Tr / 60,
				  "\t";
				print $outfilehandle $feature->charge, "\n";
				$modified_features++;
				my $segment = join "::", $feature->bin_min, $feature->bin_max;

				$segmentpopulation{$segment}++;
				$nfeatures++;

			}
			else {
				print $outfilehandle sprintf( "%.4f", $feature->mz ), "\t",
				  sprintf( "%.2f", $feature->bin_min ), "\t",
				  sprintf( "%.2f", $feature->bin_max ), "\t", $feature->Tr / 60,
				  "\t";
				print $outfilehandle $feature->charge, "\n";
				my $segment = join "::", $feature->bin_min, $feature->bin_max;

				$segmentpopulation{$segment}++;
				$nfeatures++;
			}
		}
		else {
			$redundant_features++;
		}

	}
	close($outfilehandle);

	my $max = 0;
	foreach ( keys %segmentpopulation ) {
		if ( $segmentpopulation{$_} > $max ) {
			$max = $segmentpopulation{$_};
		}
	}

	foreach my $segment (
		sort { ( split /::/, $a )[0] <=> ( split /::/, $b )[0] }
		keys %segmentpopulation
	  )
	{
		unless ( $segmentpopulation{$segment} ) {
			$segmentpopulation{$segment} = 0;
		}
		print $segment, " #features: ", $segmentpopulation{$segment}, "\t";
		if ( $segmentpopulation{$segment} ) {
			for $i ( 0 .. int( 30 * $segmentpopulation{$segment} / $max ) ) {
				print "#";
			}
		}
		print "\n";
	}

	my $notplaced_features =
	  scalar(@sorted_Tr) - $nfeatures - $redundant_features -
	  $halfredundant_features;
	print "\n--------------------------------------------\n",
	  "$nfeatures features were put into segments\n",,
"$redundant_features features discarded with identical m/z and Tr-segment\n",
	  "$modified_features features were modified with an increment of 0.01\n";
	if ($notplaced_features) {
		warn
"\n$notplaced_features features could not be placed in any segment (check segmentsize, limits, and overlap) \n";
	}
}

sub MasterMap {
	my $self = shift;
	return $self->{'MasterMap'};
}

sub concatstring {
	my $self = shift;
	my $concatstring = join "_", $self->mz, $self->Tr;
	return $concatstring;
}

sub tr_tolerance {
	my $self = shift;
	return $self->{'$tr_tolerance'};
}

sub number_is_in_array {
	my $number = shift;
	my $array  = shift;
	foreach my $arrayentry (@$array) {
		if ( $number == $arrayentry ) {
			return 1;
		}
	}
	return 0;
}


sub calcbins {
	my $self   = shift;
	my $PARAMS = $self->params;

	#	print "@{$PARAMS->bins_min}\n";
	if ( defined( $PARAMS->IL_bins_min ) && defined( $PARAMS->IL_bins_max ) ) {
		print "reading bins from xmm.def file\n";
		if ( scalar @{ $PARAMS->IL_bins_min } == scalar @{ $PARAMS->IL_bins_max } ) {
			$self->{'bins_min'} = $PARAMS->IL_bins_min;
			$self->{'bins_max'} = $PARAMS->IL_bins_max;
		}
		else {
			die
"number of lower bin boundaries needs to match number of upper bin boundaries $!";
		}
	}
	else {
		print "auto binning Tr in MasterMap \n";

		my @bins_min = ();
		my @bins_max = ();

		my $overlap     = $self->params->IL_bin_overlap;
		my $segmentsize = $self->params->IL_bin_segmentsize;
		my $limits      = $self->params->IL_bin_limits;
		my $MasterMap   = $self->MasterMap;

		my ( $mins_bin, $bins_max );

		my ( $trmin, $trmax );

		my $segmentshift = ( 1 - $overlap ) * $segmentsize;
		my $tr_tolerance = ( $segmentsize * ( 1 - $overlap ) ) / 2;
		$self->{'$tr_tolerance'} = $tr_tolerance;
		my $tr;
		my $nsegments = 0;

		if ($limits) {
			( $trmin, $trmax ) = split /,/, $limits;
		}
		else {
			$MasterMap->calc_min_max_Tr;

			if ( $PARAMS->IL_bin_offset ) {

				$trmin = $PARAMS->IL_bin_offset - $segmentshift;
				push @bins_min, sprintf( "%.2f", 0 );
				push @bins_max, sprintf( "%.2f", $PARAMS->IL_bin_offset  );
			}
			else {
				$trmin = ( $MasterMap->min_Tr - $tr_tolerance * 60 ) / 60;
				$trmin = 0 if $trmin < 0;
			}

			$trmax = $MasterMap->max_Tr / 60;
		}

		print "Tr min: ", $trmin, " Tr max: ", $trmax,
" Tr tolerance: $tr_tolerance  segmentsize: $segmentsize segmentshift: $segmentshift\n";

		for (
			$tr = $trmin ;
			$tr <= $trmax - $tr_tolerance ;
			$tr += $segmentshift
		  )
		{
			push @bins_min, sprintf( "%.2f", $tr+ 0.1 );
			push @bins_max, sprintf( "%.2f", $tr + $segmentsize );
			$nsegments++;
		}

		$self->{'bins_min'} = \@bins_min;
		$self->{'bins_max'} = \@bins_max;
	}

	print "bins start  @{ $self->{'bins_min'} } \n";
	print "bins end   @{ $self->{'bins_max'} }\n";
}

sub bins_min {
	my $self = shift;
	return $self->{'bins_min'};
}

sub bins_max {
	my $self = shift;
	return $self->{'bins_max'};
}

sub params {
	my $self = shift;
	return $self->{'PARAMS'};
}
1;
