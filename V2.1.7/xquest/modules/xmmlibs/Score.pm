package Score;
use strict;
#---------------------------------------------------------------------------
# Module: Score.pm
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

sub score {
	my $profile = shift;
	my $target  = shift;
	my $method  = shift;
	my $score=undef;

	if ( $method eq "abs" ) {
		$score=_abs($profile,$target);
	}
	elsif ( $method eq "mse" ) {
		$score=_mse($profile,$target);
	}
	else {
		die "scoring method $method is not implemented $!";
	}
	
return $score;
}

sub _abs {
	my $profile = shift;
	my $target  = shift;
	my $i;
	my $deltasum = 0;
	for $i ( 0 .. $#$profile ) {
		$deltasum += abs( $profile->[$i] - $target->[$i] );
	}
	return $deltasum /scalar(@$profile);
}

sub _mse {
	my $profile = shift;
	my $target  = shift;

	my $deltasum = 0;
	my $i;
	
	for $i ( 0 .. $#$profile ) {
		$deltasum += ( $profile->[$i] - $target->[$i] )**2;
	}
	return sqrt($deltasum / scalar(@$profile));
}

1;
