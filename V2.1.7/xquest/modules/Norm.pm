package Norm;
#---------------------------------------------------------------------------
# Module: Norm.pm
# Author(s): Thomas Walzthoeni, xQuest specific modules are based on 
# original work by Oliver Rinner.
# Description: Module for normalization.
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

sub normalize {
	my $array  = shift;
	my $method = shift;

	if ( $method eq "by_sum" ) {
		_by_sum($array);
	}
	elsif ( $method eq "by_max" ) {
		_by_max($array);
	}
	elsif ( $method eq "sqrt" ) {
		_sqrt($array);
	}

	else {
		die "normalization method $method is not implemented $!";
	}
}

sub _sqrt{
	my $array = shift;
	foreach $entry (@$array) {

		$entry = sqrt($entry);
	}
}

sub _by_sum {
	my $array = shift;
	my $sum   = 0;
	foreach my $entry (@$array) {
		$sum += $entry;
	}
	foreach $entry (@$array) {
if($sum){
		$entry /= $sum;
}
	}
}

sub _by_max {
	my $array = shift;
	my $max   = XMM_Statistics::max($array);
	foreach $entry (@$array) {
		if($max){
		$entry /= $max;
		}
	}
}

return 1;
