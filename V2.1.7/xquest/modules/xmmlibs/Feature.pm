package Feature;
use strict;
#---------------------------------------------------------------------------
# Module: Feature.pm
# Author(s): Thomas Walzthoeni, xQuest specific modules are based on 
# original work by Oliver Rinner.
# Description: Module for handling features.
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
use Norm;

my $Hatom = 1.0078250;

sub new {
	my $class       = shift();
	my $Master      = shift;	## is the Master obj ref
	my $featureinfo = shift;
	my $runIDs      = shift;
	my $verbose     = shift;

	my $self = {};

	bless $self, $class;
	#print $featureinfo->[1]."\n";
	my $id = join "::", @$featureinfo[ 0 .. 2 ];
	$self->{'id'}         = $id;
	$self->{'Master'}     = $Master;
	$self->{'m_z'}        = $featureinfo->[0];
	$self->{'originalTr'} = $featureinfo->[1];
	$self->{'Tr'}         = int( $featureinfo->[1] * 60 );
	( $self->{'charge'} = $featureinfo->[2] ) =~ s/\+//;
	unless ( $self->{'charge'} ) {
		$self->{'charge'} = 1;
	}
	$self->{'accession'}       = $featureinfo->[3];
#	$self->{'seq'}             = $featureinfo->[4];
	($self->{'seq'}=$featureinfo->[4]) =~s/^-\.//; 
	
	($self->{'seq'}=$self->{'seq'}) =~s/^\w\.//g; 
	($self->{'clean_seq'}=$self->{'seq'}) =~s/\[\d+\.\d+\]//g; 	
	$self->{'probability'}     = $featureinfo->[5];
	$self->{'nreplicates'}     = $featureinfo->[6];
	$self->{'areas'}           = [ @$featureinfo[ 7 .. ( 7 + $#$runIDs ) ] ];
	$self->{'associated_scan'} = undef;
	$self->{'isopair'}         = undef;
	$self->{'isopairindex'}    = undef;
	$self->{'pseudoSHscannum'}    = $featureinfo->[8];
	return $self;
}


sub get_pseudoSHscannum{
	my $self = shift;
	return $self->{'pseudoSHscannum'};	
}

sub mark_as_isopair {
	my $self = shift;
	$self->{'isopair'}      = shift;
	$self->{'isopairindex'} = shift;
}

sub is_isopair {
	my $self = shift;
	return $self->{'isopair'};
}

sub get_isopair_index {
	my $self = shift;
	return $self->{'isopairindex'};
}

sub associated_scan {
	my $self = shift;
	return $self->{'associated_scan'};
}

sub assign_etd_cid_pair {
	my $self = shift;
	$self->{'etd_cid_pair'} = shift;
}

sub etd_cid_pair {
	my $self = shift;
	return $self->{'etd_cid_pair'};
}

sub scan {
	my $self = shift;
	return $self->{'associated_scan'};
}

sub scanpath {
	my $self       = shift;
	my $scannumber = $self->scan->scannumber;
	my $scanID     = $self->scan->basename;
	my $charge     = $self->scan->FT_charge;
	my $scanstring;
	
	if ( $scannumber < 100 ) {
		$scanstring = join "", ".\/", $scanID, "\/", $scanID, ".000", $scannumber,
		   ".000", $scannumber, ".", $charge, ".dta";
	}
	elsif ( $scannumber < 1000 ) {
		$scanstring = join "", ".\/", $scanID, "\/", $scanID, ".00", $scannumber,
		  ".00", $scannumber, ".", $charge, ".dta";
	}elsif ( $scannumber < 10000 ){
	$scanstring = join "", ".\/", $scanID, "\/", $scanID, ".0", $scannumber,
		  ".0", $scannumber, ".", $charge, ".dta";	
		
	}elsif( $scannumber < 100000 ){
	$scanstring = join "", ".\/", $scanID, "\/", $scanID, ".", $scannumber,
	".", $scannumber, ".", $charge, ".dta";		
	}

#./B07-03342_c/B07-03342_c.1563.1563.3.dta       ./B07-03342_c/B07-03342_c.1575.1575.3.dta
	return $scanstring;
}


sub scanpath_highres {
	my $self       = shift;
	my $scannumber = $self->scan->scannumber;
	my $scanID     = $self->scan->basename;
	my $charge     = $self->scan->FT_charge;
	my $scanstring;
## PG MSMS File name contains no additional 0 values
$scanstring = join "", ".\/", $scanID, "\/", $scanID, ".", $scannumber,".", $scannumber, ".", $charge, ".dta";
#	if ( $scannumber >= 1000 ) {
#		$scanstring = join "", ".\/", $scanID, "\/", $scanID, ".", $scannumber,
#		   ".", $scannumber, ".", $charge, ".dta";
#	}
#	else {
#		$scanstring = join "", ".\/", $scanID, "\/", $scanID, ".0", $scannumber,
#		  ".0", $scannumber, ".", $charge, ".dta";
#	}

#./B07-03342_c/B07-03342_c.1563.1563.3.dta       ./B07-03342_c/B07-03342_c.1575.1575.3.dta
	return $scanstring;
}


sub original_Tr {
	my $self = shift;
	return $self->{'originalTr'};
}

sub annotate_with {
	my $self = shift;
	my $scan = shift;
	$self->{'associated_scan'} = $scan;
}

sub add_reference_scans {
	my $self  = shift;
	my $scans = shift;
	$self->{'reference_scans'} = $scans;
}

sub annotate_scan_by_defined_method {
	my $self             = shift;
	my $assignmentmethod = shift;
	my $scanarray        = shift;
	$self->{'associated_scan'} = $assignmentmethod->( $self, $scanarray );
}

sub setlower_bin {
	my $self = shift;
	$self->{'lower_bin'} = shift;
}

sub bin_min {
	my $self = shift;
	return $self->{'lower_bin'};
}

sub setupper_bin {
	my $self = shift;
	$self->{'upper_bin'} = shift;
}

sub bin_max {
	my $self = shift;

	return $self->{'upper_bin'};
}

sub matches {
	my $self           = shift;
	my $scan           = shift;			## retrieves the pseudoscanobject
	my $PARAMS         = $self->params;
	my $addtofeaturemz = 0;     #$self->params->feature_plusmz;
	my $mztolerance    = $PARAMS->MS1tolerance;
	my $feature_mz     = $self->mz;
	my $scan_mz        = $scan->FT_mz;
	my $fractionationtype = $scan->fractionationtype;

	my $feature_charge             = $self->charge;
	my $scan_charge                = $scan->FT_charge;
	my $require_same_charge        = $PARAMS->require_same_charge;
	my $require_fractionation_type = $PARAMS->require_fractionationtype;
	my $monotoggle                 = $PARAMS->monotoggle;

	if (
		(
			( abs( $feature_mz - $scan_mz ) < $mztolerance )
			|| (
				(
					$monotoggle
					&& (
						(
							abs(
								$feature_mz - ( 1 / $feature_charge ) - $scan_mz
							) < $mztolerance
						)
						|| (
							abs(
								$feature_mz + ( 1 / $feature_charge ) - $scan_mz
							) < $mztolerance
						)
					)
				)
			)
		)
		&& ( !$require_same_charge || ( $feature_charge == $scan_charge ) )
		&& ( !$require_fractionation_type
			|| ( $fractionationtype =~ /$require_fractionation_type/ ) )
	  )
	{
		return 1;
	}
	else {
		return 0;
	}

}

sub params {
	my $self = shift;
	return $self->MasterMap->params;
}

sub mz {
	my $self = shift;
	return $self->{'m_z'};
}

sub Mr {
	my $self = shift;
	return $self->{'m_z'} * $self->charge - $self->charge * $Hatom;
}

sub Tr {
	my $self = shift;
	return $self->{'Tr'};
}

sub charge {
	my $self = shift;
	return $self->{'charge'};
}

sub accession {
	my $self = shift;
	unless ( $self->{'accession'} ) {
		return "";
	}
	else {
		return $self->{'accession'};
	}
}

sub id {
	my $self = shift;
	return $self->{'id'};
}

sub translation {
	my $self = shift;
	return $self->MasterMap->translationtable;
}

sub MasterMap {
	my $self = shift;
	return $self->{'Master'};
}

sub seq {
	my $self = shift;
	return $self->{'seq'};
}


sub clean_seq {
	my $self = shift;
	return $self->{'clean_seq'};
}

sub p {
	my $self = shift;
	return $self->{'probability'};
}

sub nreplicates {
	my $self = shift;
	return $self->{'nreplicates'};
}

sub areaarray {
	my $self = shift;
#	if ( $self->params->normalize_profiles ) {
#		Norm::normalize( $self->{'areas'}, $self->params->normalize_profiles );
#	}

	return $self->{'areas'};
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

sub matches_requirements {
	my $self   = shift;
	my $PARAMS = $self->params;
	my $Tr     = $self->Tr;
	my $featurecount=$self->nreplicates;
	my $charge = $self->charge;

	if (   ( $Tr >= $PARAMS->min_Tr )
		&& ( $Tr <= $PARAMS->max_Tr )
		&& ($featurecount>=$PARAMS->minfeaturecount )
		&& ( number_is_in_array( $PARAMS->charge_states, $charge ) ) )
	{
		return 1;
	}
	else {
		return 0;
	}
}

sub area {
	my $self  = shift;
	my $runid = shift;
	return $self->areaarray->[$runid];
}

sub lcid {
	my $self       = shift;
	my $areas      = $self->areaarray;
	my $definedids = $self->MasterMap->getdefinedIDs;
	my @ids;
	foreach my $runid (@$definedids) {
		if ( $areas->[$runid] ) {
			push @ids, 1;
		}
		else {
			push @ids, 0;
		}
	}
	return \@ids;
}

sub totalarea {
	my $self   = shift;
	my $runids = $self->MasterMap->getdefinedIDs;
	my $total  = 0;
	foreach my $id (@$runids) {
		$total += $self->areaarray->[$id];
	}
	return $total;
}

sub meanintensity {
	my $self   = shift;
	my $runids = $self->MasterMap->getdefinedIDs;
	my $total  = 0;
	my $n      = 0;
	foreach my $id (@$runids) {
#		if ( $self->ispresent($id) ) {
#			$n++;
#		}
		$total += $self->areaarray->[$id];
	}
	#return $total / $n;
	return $total;
}

sub areasum {
	my $self   = shift;
	my $runids = shift;
	my $total  = 0;
	foreach my $id (@$runids) {
		$total += $self->areaarray->[$id];
	}
	return $total;
}

sub ratio {
	my $self   = shift;
	my $runids = shift;

	if ( $self->area( $runids->[0] ) && $self->area( $runids->[1] ) ) {
		return $self->area( $runids->[0] ) / $self->area( $runids->[1] );
	}
	else {
		return undef;
	}

}

sub profile {
	my $self     = shift;
	my $runids   = shift;
	my @profiles = ();
	foreach my $id (@$runids) {
		if ( $self->area($id) ) {
			push @profiles, $self->area($id);
		}
		else {
			push @profiles, 0;
		}
	}
	

	return \@profiles;
}

sub ispresent {
	my $self        = shift;
	my $requiredids = shift;

	my @profiles  = ();
	my $ispresent = 1;

	foreach my $id (@$requiredids) {
		$ispresent *= $self->area($id);
	}
	return $ispresent;
}

sub nrepeats {
	my $self     = shift;
	my $runids   = shift;
	my @profiles = ();
	my $nrepeats = 0;
	foreach my $id (@$runids) {
		if ( $self->area($id) ) {
			$nrepeats++;
		}
	}
	return $nrepeats;
}

1;

