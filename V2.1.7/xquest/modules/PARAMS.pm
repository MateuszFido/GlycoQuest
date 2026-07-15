package PARAMS;
use strict;
#---------------------------------------------------------------------------
# Module: PARAMS.pm
# Author(s): Thomas Walzthoeni, xQuest specific modules are based on 
# original work by Oliver Rinner.
# Description: Module for parameter handling.
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
	my $class   = shift();
	my $deffile = shift;
	my $verbose = shift;

	my $self = {};
	#$verbose=1;
	$self->{'PARAMS'} = readtables( $deffile, $verbose );
	bless $self, $class;
}

sub require_same_charge {
	my $self = shift;
	return $self->params->{'require_same_charge'};
}


sub levelto{
my $self=shift;
return 0 unless defined ($self->params->{'require_same_charge'});
return $self->params->{'levelto'};
}

sub average_nmostintense{
	my $self = shift;
	die "parameter average_nmostintense is not defined" unless defined($self->params->{'average_nmostintense'});
	return $self->params->{'average_nmostintense'};
}

############################# ETD parameters ########################
sub etd_min_mz {
	my $self = shift;
	if ( $self->params->{'etd_min_mz'} ) {
		return $self->params->{'etd_min_mz'};
	}
	else {
		return 200;
	}
}

sub etd_max_mz {
	my $self = shift;
	if ( $self->params->{'etd_max_mz'} ) {
		return $self->params->{'etd_max_mz'};
	}
	else {
		return 1600;
	}
}

sub etd_binsize {
	my $self = shift;
	if ( $self->params->{'etd_binsize'} ) {
		return $self->params->{'etd_binsize'};
	}
	else {
		return 1;
	}
}

sub etd_min_dotproduct {
	my $self = shift;
	if ( $self->params->{'etd_min_dotproduct'} ) {
		return $self->params->{'etd_min_dotproduct'};
	}
	else {
		return 1;
	}
}



sub etd_Tr_diff{
	my $self = shift;
	if ( $self->params->{'etd_Tr_diff'} ) {
		return $self->params->{'etd_Tr_diff'};
	}
	else {
		return 15;
	}
}


sub etd_precursor_tolerance {
	my $self = shift;
	if ( $self->params->{'etd_precursor_tolerance'} ) {
		return $self->params->{'etd_precursor_tolerance'};
	}
	else {
		return 2.0;
	}
}

sub etd_dynamic_range{
	my $self = shift;
	if ( $self->params->{'etd_dynamic_range'} ) {
		return $self->params->{'etd_dynamic_range'};
	}
	else {
		return 100;
	}
}

#################################################################################

sub require_fractionationtype {
	my $self = shift;
	return $self->params->{'require_fractionationtype'};
}

sub Isopair_require_same_lcid {
	my $self = shift;
	return $self->params->{'Isopair_require_same_lcid'};
}

sub IL_charges {
	my $self = shift;
	if ( $self->params->{'IL_charges'} == 0 ) {
		return 0;
	}
	elsif ( $self->params->{'IL_charges'} =~ /\A[\d+,]+\Z/ ) {
		my @bins = split /,/, $self->params->{'IL_charges'};
		return \@bins;
	}
	else {
		warn "parameter IL_charges ", $self->params->{'IL_charges'},
		  " is not correctly set ";
		die $!;
	}
}

sub readtables {
	my $deffile = shift;
	my $verbose = shift;
	my %PARAMS  = ();
	open DEFFILE, "<$deffile"
	  or die "cannot open xmm definition file $deffile $!";
	while (<DEFFILE>) {
		chomp;
		my @tmp = split;
		$PARAMS{ $tmp[0] } = $tmp[1];
		$verbose && print $tmp[0], "\t", $tmp[1], "\n";
	}
	return \%PARAMS;
}

sub scanlevels {
	my $self = shift;
	my @scanlevels = split /,/, $self->params->{'scanlevels'};
	return \@scanlevels;
}

sub msInstrumentIDs {
	my $self = shift;
	if ( defined( $self->params->{'msInstrumentIDs'} ) ) {
		my @instrumentids = split /,/, $self->params->{'msInstrumentIDs'};
		return \@instrumentids;
	}
	else {
		return (0);
	}
}

sub inclusion_list {
	my $self = shift;
	return $self->params->{'inclusion_list'};
}

sub IL_bin_overlap {
	my $self = shift;

	my $overlap = $self->params->{'IL_bin_overlap'} / 100;
	if ( $overlap >= 0 && $overlap <= 1 ) {
		return $overlap;
	}
	else {
		die "segement overlap of $overlap cannot be applied $!";
	}
}

sub inclusionlisttype {
	my $self = shift;
	return $self->params->{'inclusionlisttype'};
}

sub IL_bin_offset {
	my $self = shift;
	return $self->params->{'IL_bin_offset'};

}

sub IL_bins_min {
	my $self = shift;
	if ( !( $self->params->{'bins_min'} =~ /\d+,\d+/ ) ) {
		return undef;
	}
	else {
		my @bins = split /,/, $self->params->{'bins_min'};

		return \@bins;
	}

}

sub IL_bins_max {
	my $self = shift;
	if ( !( $self->params->{'bins_max'} =~ /\d+,\d+/ ) ) {
		return undef;
	}
	else {
		my @bins = split /,/, $self->params->{'bins_max'};
		return \@bins;
	}
}

sub IL_bin_segmentsize {
	my $self = shift;
	return $self->params->{'IL_bin_segmentsize'};
}

sub printisotopicscanpairs {
	my $self = shift;
	return $self->params->{'printisotopicscanpairs'};
}

sub printlightonlypairs {
	my $self = shift;
	return $self->params->{'printlightonlypairs'};
}

sub printheavyonlypairs {
	my $self = shift;
	return $self->params->{'printheavyonlypairs'};
}

sub nmostintensefeatures {
	my $self = shift;
	if ( $self->params->{'IL_nmostintense'} == 0 ) {
		return 0;
	}
	else {
		my @tmp = split /,/, $self->params->{'IL_nmostintense'};
		@tmp = sort { $a <=> $b } @tmp;
		return \@tmp;
	}
}

sub IL_bin_limits {
	my $self = shift;
	return $self->params->{'IL_bin_limits'};
}


sub strip_peptides{
	my $self = shift;
	return $self->params->{'strip_peptides'};
}

sub normalize_MasterMap {
	my $self = shift;
	return $self->params->{'normalize_MasterMap'};
}

sub profiler {
	my $self = shift;
	return $self->params->{'profiler'};
}


sub single_quant{
	my $self = shift;
	return $self->params->{'single_quant'};
}

sub ratiolizer {
	my $self = shift;
	return $self->params->{'ratiolizer'};
}

sub matchtype {
	my $self = shift;
	return $self->params->{'matchtype'};
}

sub monotoggle {
	my $self = shift;
	return $self->params->{'monotoggle'};
}

sub annotationmethod {
	my $self      = shift;
	my $matchtype = $self->matchtype;
	if ( $matchtype eq "Tr_distance_2apex" ) {
		return \&mzXMLscan::assign_by_distance2apex;
	}
	elsif ( $matchtype eq "most_intense" ) {
		return \&mzXMLscan::assign_by_basepeakIntensity;
	}
	else {
		die "peak matching method \"$matchtype\" is not defined $0";
	}
}

sub min_Tr {
	my $self = shift;
	return $self->params->{'min_Tr'} * 60;
}

sub max_Tr {
	my $self = shift;
	return $self->params->{'max_Tr'} * 60;
}

sub charge_states {
	my $self = shift;
	my @chargestates = split /,/, $self->params->{'charge_states'};
	return \@chargestates;
}

sub delta {
	my $self = shift;
	unless ( defined( $self->params->{'isotopeshift'} ) ) {
		die "parameter isotopeshift is not defined $!";
	}
	my @shifts = split /,/, $self->params->{'isotopeshift'};
	return \@shifts;
}

sub deltashift {
	my $self = shift;
	unless ( defined( $self->params->{'isotopeshift'} ) ) {
		die "parameter deltashift is not defined $!";
	}
	return $self->params->{'deltashift'};
}

sub tripleshift {
	my $self = shift;
	return $self->params->{'tripleshift'};
}

#
#
#sub isotopic_mass_shift {
#	my $self = shift;
#	unless(defined($self->params->{'isotopic_mass_shift'})){
#	 die "parameter isotopic_mass_shift is not defined $!";
#	}
#	return $self->params->{'isotopic_mass_shift'};
#}

sub lightrunids {
	my $self = shift;
	my @runids = split /,/, $self->params->{'lightrunids'};
	return \@runids;
}

sub heavyrunids {
	my $self = shift;
	my @runids = split /,/, $self->params->{'heavyrunids'};
	return \@runids;
}

sub featurestatIDs {
	my $self = shift;
	my @runids = split /,/, $self->params->{'featurestatIDs'};
	return \@runids;
}

sub MS1tolerance {
	my $self = shift;
	return $self->params->{'mz_tolerance'};
}

sub Isopair_Mr_tolerance {
	my $self = shift;
	unless ( defined( $self->params->{'Isopair_Mr_tolerance'} ) ) {
		die "parameter Isopair_Mr_tolerance is not defined $!";
	}
	return $self->params->{'Isopair_Mr_tolerance'};
}

sub Isopair_Mr_tolerance_measure {
	my $self             = shift;
	my $tolerancemeasure = $self->params->{'Isopair_Mr_tolerance_measure'};
	unless ( $tolerancemeasure eq "amu" || $tolerancemeasure eq "ppm" ) {
		warn
"Mr tolerance measure $tolerancemeasure is not defined, defaulting to \"amu\" \n";
		return "amu";
	}
	else {
		return $tolerancemeasure;
	}

}

sub feature_plusmz {
	my $self = shift;
	return $self->params->{'feature_plusmz'};
}

sub pairratio {
	my $self = shift;
	return $self->params->{'pairratio'};
}

sub deltaTr {
	my $self = shift;
	return int( $self->params->{'deltaTr'} * 60 );
}

sub Trshift {
	my $self = shift;
	unless ( defined( $self->params->{'Trshift'} ) ) {
		warn "parameter Trshift is not defined, assuming 0s\n";
		return 0;
	}
	else {
		return int( $self->params->{'Trshift'} * 60 );
	}
}

sub Isopair_Tr_tolerance {
	my $self = shift;
	unless ( defined( $self->params->{'Isopair_Tr_tolerance'} ) ) {
		die "parameter Isopair_Tr_tolerance is not defined $!";
	}
	return int( $self->params->{'Isopair_Tr_tolerance'} * 60 );
}

sub normalize_profiles {
	my $self = shift;
	return $self->params->{'normalize_profiles'};
}

sub averaging {
	my $self = shift;
	return $self->params->{'averaging'};
}

sub target {
	my $self = shift;
	return $self->params->{'target'};
}

sub score {
	my $self = shift;
	return $self->params->{'score'};
}

sub mincount {
	my $self = shift;
	return $self->params->{'minfeaturecount'};
}

sub minfeaturecount {
	my $self = shift;
	return $self->params->{'minfeaturecount'};
}

sub renormalize {
	my $self = shift;
	return $self->params->{'renormalize'};
}

sub params {
	my $self = shift;
	return $self->{'PARAMS'};
}

sub minprob {
	my $self = shift;
	return $self->{'PARAMS'}->{'minprob'};
}

sub Tr_offset {
	my $self = shift;
	return int( $self->{'PARAMS'}->{'Tr_offset'} * 60 );
}

sub translationtable {
	my $self = shift;
	return $self->{'PARAMS'}->{'translationtable'};
}

sub verbose {
	my $self = shift;
	return $self->{'PARAMS'}->{'verbose'};
}
sub quantify_by_nhighest{
	my $self = shift;
	return $self->{'PARAMS'}->{'quantify_by_nhighest'};
}
sub runIDs {
	my $self = shift;
	return $self->{'PARAMS'}->{'runIDs'};
}

sub requiredrunIDs {
	my $self = shift;

	return split /,/, $self->{'PARAMS'}->{'requiredrunIDs'};
}

1;
