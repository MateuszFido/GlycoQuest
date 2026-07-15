package XMM_ETD_CID;
use strict;
#---------------------------------------------------------------------------
# Module: XMM_ETD_CID.pm
# Author(s): Thomas Walzthoeni, xQuest specific modules are based on 
# original work by Oliver Rinner.
# Description: xmm specific module for scan handling.
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
	my $class = shift;
	my $self  = {};
	bless $self, $class;
	
	$self->{'cid_scan'}     = shift;
	$self->{'etd_scan'}     = shift;
	$self->{'scanpair_id'}     = shift;
	
	
	unless($self->{'etd_scan'}->precursor_mz == $self->{'cid_scan'}->precursor_mz){
		die "peakpairs ",$self->{'cid_scan'}->id, " and ",$self->{'etd_scan'}->id, " do not have matching precursor ions $!";
	}
	$self->{'precursor_mz'} = $self->etd_scan->precursor_mz;
	#$self->{'precursor_charge'} = $self->etd_scan->precursor_charge;
	
	$self->{'precursor_Tr'} = $self->etd_scan->Tr;

	return $self;
}

sub scanpair_id{
	my $self = shift;
	return $self->{'scanpair_id'};
}


sub cid_scan {
	my $self = shift;
	return $self->{'cid_scan'};
}

sub etd_scan {
	my $self = shift;
	return $self->{'etd_scan'};
}

sub precursor_mz {
	my $self = shift;
	return $self->{'precursor_mz'};
}

sub precursor_Tr {
	my $self = shift;
	return $self->{'precursor_Tr'};
}




1;
