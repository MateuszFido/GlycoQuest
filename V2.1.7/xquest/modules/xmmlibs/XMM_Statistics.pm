package XMM_Statistics;
use strict;
#---------------------------------------------------------------------------
# Module: XMM_Statistics.pm
# Author(s): Thomas Walzthoeni, xQuest specific modules are based on 
# original work by Oliver Rinner.
# Description: xmm specific statistic module.
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
sub mean {
	my $array = shift;
	my $sum   = 0;
	my $n     = 0;
	foreach (@$array) {
		$sum += $_;
		$n++;
	}
	if ($n) {
		return $sum / $n;
	}
	else {
		return 0;
	}
}

sub sum{
	my $array = shift;
	my $sum   = 0;
	foreach (@$array) {
		$sum += $_;
	}
return $sum;
}


sub average {
	my $array = shift;
	return mean($array);
}

sub mean_array {
	my $array = shift;

	my $n = 0;

	my %profilehash = ();
	my ( $i, $j );
	for $i ( 0 .. $#$array ) {
		my @meanarray = ();
		foreach $j ( 0 .. $#{ $array->[$i] } ) {
			push @{ $profilehash{$j} }, $array->[$i]->[$j];
		}
	}

	my @profilearray = ();
	foreach my $runID ( sort keys %profilehash ) {
		push @profilearray, mean( $profilehash{$runID} );
	}

	return \@profilearray;
}

sub median_array {
	my $array = shift;

	my $n = 0;

	my %profilehash = ();
	my ( $i, $j );
	for $i ( 0 .. $#$array ) {
		my @meanarray = ();
		foreach $j ( 0 .. $#{ $array->[$i] } ) {
			push @{ $profilehash{$j} }, $array->[$i]->[$j];
		}
	}

	my @profilearray = ();
	foreach my $runID ( sort keys %profilehash ) {
		push @profilearray, median( $profilehash{$runID} );
	}
	return \@profilearray;
}

sub sem_array {
	my $array = shift;

	my $n = 0;

	my %profilehash = ();
	my ( $i, $j );
	for $i ( 0 .. $#$array ) {
		my @meanarray = ();
		foreach $j ( 0 .. $#{ $array->[$i] } ) {
			push @{ $profilehash{$j} }, $array->[$i]->[$j];
		}
	}

	my @profilearray = ();
	foreach my $runID ( sort keys %profilehash ) {
		push @profilearray, sem( $profilehash{$runID} );
	}
	return \@profilearray;
}

sub stdev_array {
	my $array = shift;

	my $n = 0;

	my %profilehash = ();
	my ( $i, $j );
	for $i ( 0 .. $#$array ) {
		my @meanarray = ();
		foreach $j ( 0 .. $#{ $array->[$i] } ) {
			push @{ $profilehash{$j} }, $array->[$i]->[$j];
		}
	}

	my @profilearray = ();
	foreach my $runID ( sort keys %profilehash ) {
		push @profilearray, stdev( $profilehash{$runID} );
	}
	return \@profilearray;
}

sub max {
	my $vector = shift;
	my @sorted = sort { $b <=> $a } @$vector;
	return $sorted[0];
}

sub n {
	my $array = shift;
	return scalar(@$array);
}

sub var {
	my $vector = shift;
	my $mean   = mean($vector);
	my $n      = n($vector);

	my $squaresum = 0;
	foreach (@$vector) {
		$squaresum += ( $_ - $mean )**2;
	}
	if ( $n > 0 ) {
		return ( $squaresum / $n );
	}
	else {
		return 0;
	}
}

sub median {
	my $vector = shift;
	my @sorted = sort { $a <=> $b } @$vector;
	return $sorted[ int( scalar(@sorted) / 2 ) ];
}

sub getpercentile {
	my $vector = shift;
	my $value  = shift;
	my @sorted = sort { $a <=> $b } @$vector;
	my $i;
	for ( $i = 0 ; $i <= $#sorted ; $i++ ) {
		if ( $value < $sorted[$i] ) {

			#	print "$i: $value\n";
			return $i / scalar(@sorted);
		}
	}
}

sub topn {
	my $vector = shift;
	my $rank   = shift;
	unless ($rank) {
		$rank = 0.5;
	}
	my @sorted      = sort { $b <=> $a } @$vector;
	my $toptenindex = scalar(@$vector) *(1-$rank);
	my @toptenslice = @sorted[ 0 .. int($toptenindex) ];

	return mean( \@toptenslice );
}

sub stdev {
	my $vector = shift;
	my $n      = n($vector);

	if ( $n > 0 ) {
		return ( sqrt( var($vector) ) );
	}
	else {
		return 0;
	}
}

sub sem {
	my $vector = shift;
	my $n      = n($vector);
	if ( $n > 0 ) {
		return ( stdev($vector) / sqrt($n) );
	}
	else {
		return 0;
	}
}

1;
