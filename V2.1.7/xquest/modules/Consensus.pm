package Consensus;
use strict;
#---------------------------------------------------------------------------
# Module: Consensus.pm
# Author(s): Thomas Walzthoeni, xQuest specific modules are based on 
# original work by Oliver Rinner.
# Description: Specific functions.
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

sub get_binned_peaks {
	my $peaklist = shift;
	my $PARAMS  = shift;

	my $binsize = $PARAMS->etd_binsize;
	my $dynamic_range=$PARAMS->etd_dynamic_range;
	my $min_mz=$PARAMS->etd_min_mz;
	my $max_mz=$PARAMS->etd_max_mz;
	
	
	my @peaks = ();
	foreach my $peakpair (@$peaklist) {
		push @peaks, $peakpair->[1];
	}
	my $max = XMM_Statistics::max( \@peaks );

	my @filtered_peakpairs  = ();
	my @filteredintensities = ();
	my $minintensity        = $max / $dynamic_range;
	foreach my $peakpair (@$peaklist) {
		if ( $peakpair->[1] >= $minintensity ) {
			push @filtered_peakpairs,
			  [ $peakpair->[0], sqrt( $peakpair->[1] ) ];
			push @filteredintensities, $peakpair->[1];
		}
	}

	my %bins = ();
	for my $i ( int( $min_mz * $binsize ) .. int( $max_mz * $binsize ) ) {
		$bins{$i} = 0;
	}

	foreach my $peakpair (@filtered_peakpairs) {
		my $mz = $peakpair->[0];
		if($mz>=$min_mz && $mz <=$max_mz){
			$bins{ int( $mz * $binsize ) } += $peakpair->[1];
		}
	}

	my $normalizationfactor =
	  sqrt( XMM_Statistics::sum( \@filteredintensities ) );

	my @bins = ();

	foreach my $binned_mz ( sort { $a <=> $b } keys %bins ) {
		push @bins, [ $binned_mz, $bins{$binned_mz} / $normalizationfactor ];
	}
	return \@bins;
}

sub dotproduct {
	my $peaklist1  = shift;
	my $peaklist2  = shift;
	my $PARAMS=shift;
	my $normpeaks1 = get_binned_peaks($peaklist1,$PARAMS);
	my $normpeaks2 = get_binned_peaks($peaklist2,$PARAMS);

	unless ( $#$normpeaks1 == $#$normpeaks2 ) {
		die "vector lengths ",$#$normpeaks1," and ",$#$normpeaks2," do not match $!";
	}
	my $dotproduct = 0;
	for my $i ( 0 .. $#$normpeaks1 ) {
		$dotproduct += $normpeaks1->[$i]->[1] * $normpeaks2->[$i]->[1];
	}
	return $dotproduct;
}




1;
