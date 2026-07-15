package Isopairs;
use strict;
#---------------------------------------------------------------------------
# Module: Isopairs.pm
# Author(s): Thomas Walzthoeni, xQuest specific modules are based on 
# original work by Oliver Rinner.
# Description: Module for isotopic scan pair processing.
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
use File::Basename;
use File::Spec;

sub new {
	my $class = shift;
	my $self  = {};
	bless $self, $class;
	$self->{'light_feature'} = shift;
	$self->{'heavy_feature'} = shift;
	$self->{'PARAMS'}        = shift;
	$self->{'delta'}         = shift;
	return $self;
}

sub print_deltafeatures_xls {
	my $self          = shift;
	my $outfilehandle = shift;

	my $Mrlight   = $self->light_feature->Mr;
	my $Mrheavy   = $self->heavy_feature->Mr;
	
	my $delta     = sprintf( "%.4f", $Mrheavy - $Mrlight );
	my $Tr_light  = sprintf( "%.2f", $self->light_feature->Tr / 60 );
	my $Tr_heavy  = sprintf( "%.2f", $self->heavy_feature->Tr / 60 );
	my $delta_Tr  = $Tr_light - $Tr_heavy;
	my $mz        = $self->mz;
	my $charge    = $self->charge;
	my $spectrum1 = undef;
	my $spectrum2 = undef;
	if ( $self->light_feature->scan ) {
		$spectrum1 = $self->light_feature->scan->scannumber;
	}
	if ( $self->heavy_feature->scan ) {
		$spectrum2 = $self->heavy_feature->scan->scannumber;
	}
	my $profile     = $self->get_profile;
	my $repeatcount = 0;
	foreach my $lcid (@$profile) {
		if ($lcid) {
			$repeatcount++;
		}
	}

	print $outfilehandle
"$mz\t$Tr_light\t$Tr_heavy\t$delta_Tr\t$charge\t$Mrlight\t$Mrheavy\t$delta\t$spectrum1\t$spectrum2\t$repeatcount\t";
	print $outfilehandle join "\t", @$profile;
	print $outfilehandle "\n";

#1974.924633::50.320129  659.316551      3       ./B07-03342_c/B07-03342_c.1563.1563.3.dta       ./B07-03342_c/B07-03342_c.1575.1575.3.dta
}

sub print_isopair_into_spectrumlist {
	my $self          = shift;
	my $outfilehandle = shift;

	#my $id        = join "::", $self->Mr, $self->Tr;

	my $mz        = $self->mz;
	my $charge    = $self->charge;
	my $spectrum1 = $self->light_feature->scanpath;
	my $spectrum2 = $self->heavy_feature->scanpath;
	my $id        = join ",", basename($spectrum1), basename($spectrum2);
	my $scan1=$self->light_feature->scan->scannumber;
	my $scan2=$self->heavy_feature->scan->scannumber;
	my $Trscan1= $self->light_feature->scan->Tr;
	my $Trscan2= $self->heavy_feature->scan->Tr;
	my $mzscan1 = $self->light_feature->scan->FT_mz;
	my $mzscan2 = $self->heavy_feature->scan->FT_mz;
	
	print $outfilehandle $id, "\t", $mz, "\t", $charge, "\t", $spectrum1, "\t", $spectrum2, "\tlight\theavy","\t",$scan1.":".$scan2,"\t",$Trscan1.":".$Trscan2,"\t",$mzscan1.":".$mzscan2,"\n";

#1974.924633::50.320129  659.316551      3       ./B07-03342_c/B07-03342_c.1563.1563.3.dta       ./B07-03342_c/B07-03342_c.1575.1575.3.dta
}

sub print_isopair_into_spectrumlist_highres {
	my $self          = shift;
	my $outfilehandle = shift;

	#my $id        = join "::", $self->Mr, $self->Tr;

	my $mz        = $self->mz;
	my $charge    = $self->charge;
	my $spectrum1 = $self->light_feature->scanpath_highres;
	my $spectrum2 = $self->heavy_feature->scanpath_highres;
	my $id        = join ",", basename($spectrum1), basename($spectrum2);
	my $scan1=$self->light_feature->scan->scannumber;
	my $scan2=$self->heavy_feature->scan->scannumber;
	my $Trscan1= $self->light_feature->scan->Tr;
	my $Trscan2= $self->heavy_feature->scan->Tr;
	my $mzscan1 = $self->light_feature->scan->FT_mz;
	my $mzscan2 = $self->heavy_feature->scan->FT_mz;
		
	print $outfilehandle $id, "\t", $mz, "\t", $charge, "\t", $spectrum1, "\t",$spectrum2, "\tlight\theavy","\t",$scan1.":".$scan2,"\t",$Trscan1.":".$Trscan2,$mzscan1.":".$mzscan2,"\n";

#1974.924633::50.320129  659.316551      3       ./B07-03342_c/B07-03342_c.1563.1563.3.dta       ./B07-03342_c/B07-03342_c.1575.1575.3.dta
}




sub print_etd_isopair_into_spectrumlist {
	my $self          = shift;
	my $mgffilehandle = shift;
	my $isopairlistfilehandle   = shift;

	my $mz                 = $self->mz;
	my $charge             = $self->charge;
	my $etd_cid_pair_light = $self->light_feature->etd_cid_pair;
	my $etd_cid_pair_heavy = $self->heavy_feature->etd_cid_pair;
	my $etd_scan_light = $etd_cid_pair_light->etd_scan;
	my $etd_scan_heavy = $etd_cid_pair_heavy->etd_scan;
	my $cid_scan_light = $etd_cid_pair_light->cid_scan;
	my $cid_scan_heavy = $etd_cid_pair_heavy->cid_scan;
	
	


	$etd_scan_light->print_into_mgf($mgffilehandle,$mz,$charge,$etd_scan_light->id."_etd");
	$etd_scan_heavy->print_into_mgf($mgffilehandle,$mz,$charge,$etd_scan_heavy->id.'_etd');
	$cid_scan_light->print_into_mgf($mgffilehandle,$mz,$charge,$cid_scan_light->id."_cid");
	$cid_scan_heavy->print_into_mgf($mgffilehandle,$mz,$charge,$cid_scan_heavy->id.'_cid');

	print $isopairlistfilehandle $etd_scan_light->id,"_etd","\t",$etd_scan_heavy->id,"_etd","\n";
	print $isopairlistfilehandle $cid_scan_light->id,"_cid","\t",$cid_scan_heavy->id,"_cid","\n";
}

sub print_lightonly_isopair_into_spectrumlist_highres {
	my $self          = shift;
	my $outfilehandle = shift;

	#	my $id        = join "::", $self->Mr, $self->Tr;
	my $mz        = $self->mz;
	my $charge    = $self->charge;
	my $spectrum1 = $self->light_feature->scanpath_highres;
	my $spectrum2 = $self->light_feature->scanpath_highres;
	my $id        = join ",", basename($spectrum1), basename($spectrum2);
	my $scan1=$self->light_feature->scan->scannumber;
	my $scan2=$self->light_feature->scan->scannumber;
	my $Trscan1= $self->light_feature->scan->Tr;
	my $Trscan2= $self->light_feature->scan->Tr;
	my $mzscan1 = $self->light_feature->scan->FT_mz;
	my $mzscan2 = $self->light_feature->scan->FT_mz;
	
	print $outfilehandle $id, "\t", $mz, "\t", $charge, "\t", $spectrum1, "\t",
	  $spectrum2, "\tlight\tlight","\t",$scan1.":".$scan2,"\t",$Trscan1.":".$Trscan2,$mzscan1.":".$mzscan2,"\n";

#1974.924633::50.320129  659.316551      3       ./B07-03342_c/B07-03342_c.1563.1563.3.dta       ./B07-03342_c/B07-03342_c.1575.1575.3.dta
}

sub print_lightonly_isopair_into_spectrumlist {
	my $self          = shift;
	my $outfilehandle = shift;

	#	my $id        = join "::", $self->Mr, $self->Tr;
	my $mz        = $self->mz;
	my $charge    = $self->charge;
	my $spectrum1 = $self->light_feature->scanpath;
	my $spectrum2 = $self->light_feature->scanpath;
	my $id        = join ",", basename($spectrum1), basename($spectrum2);
	my $scan1=$self->light_feature->scan->scannumber;
	my $scan2=$self->light_feature->scan->scannumber;
	my $Trscan1= $self->light_feature->scan->Tr;
	my $Trscan2= $self->light_feature->scan->Tr;
	my $mzscan1 = $self->light_feature->scan->FT_mz;
	my $mzscan2 = $self->light_feature->scan->FT_mz;
	
	print $outfilehandle $id, "\t", $mz, "\t", $charge, "\t", $spectrum1, "\t",
	  $spectrum2, "\tlight\tlight","\t",$scan1.":".$scan2,"\t",$Trscan1.":".$Trscan2,$mzscan1.":".$mzscan2,"\n";

#1974.924633::50.320129  659.316551      3       ./B07-03342_c/B07-03342_c.1563.1563.3.dta       ./B07-03342_c/B07-03342_c.1575.1575.3.dta
}

sub print_heavyonly_isopair_into_spectrumlist {
	my $self          = shift;
	my $outfilehandle = shift;

	#	my $id        = join "::", $self->Mr, $self->Tr;
	my $mz        = $self->mz;
	my $charge    = $self->charge;
	my $spectrum1 = $self->heavy_feature->scanpath;
	my $spectrum2 = $self->heavy_feature->scanpath;
	my $id        = join ",", basename($spectrum1), basename($spectrum2);

	print $outfilehandle $id, "\t", $mz, "\t", $charge, "\t", $spectrum1, "\t",
	  $spectrum2, "\theavy\theavy\n";

#1974.924633::50.320129  659.316551      3       ./B07-03342_c/B07-03342_c.1563.1563.3.dta       ./B07-03342_c/B07-03342_c.1575.1575.3.dta
}

sub params {
	my $self = shift;
	return $self->{'PARAMS'};
}

sub mz {
	my $self = shift;
	return $self->light_feature->mz;
}

sub Tr {
	my $self = shift;
	return $self->light_feature->Tr;
}

sub Mr {
	my $self = shift;
	return $self->light_feature->Mr;
}

sub delta {
	my $self = shift;
	return $self->{'delta'};
}

sub charge {
	my $self = shift;
	if ( $self->light_feature->charge == $self->heavy_feature->charge ) {
		return $self->light_feature->charge;
	}
	else {
		warn "charges do not match for ", $self->id, "$!\n";
		return 0;
	}
}

sub FT_charge_light {
	my $self = shift;
	if ( $self->light_feature->scan ) {
		return $self->light_feature->scan->FT_charge;
	}
	else {
		return 0;
	}
}

sub FT_charge_heavy {
	my $self = shift;
	if ( $self->heavy_feature->scan ) {
		return $self->heavy_feature->scan->FT_charge;
	}
	else {
		return 0;
	}
}

sub profile {
	my $self    = shift;
	my $verbose = $self->params->verbose;

	my $lightfeature        = $self->light_feature;
	my $heavyfeature        = $self->heavy_feature;
	my $lightfeatureprofile = $lightfeature->areaarray;
	my $heavyfeatureprofile = $heavyfeature->areaarray;
	my @profilevector       = ();

	$verbose && print "____________________________\n";
	$verbose && print "light @$lightfeatureprofile\n";
	$verbose && print "heavy @$heavyfeatureprofile\n";

	my $meanprofile =
	  XMM_Statistics::mean_array(
		[ $lightfeatureprofile, $heavyfeatureprofile ] );

	if ( $self->params->normalize_profiles ) {
		Norm::normalize( $meanprofile, $self->params->normalize_profiles );
	}
	$verbose && print "mean: @$meanprofile\n";
	$self->{'profile'} = $meanprofile;
}

sub get_profile {
	my $self = shift;
	return $self->{'profile'};
}

sub light_feature {
	my $self = shift;
	return $self->{'light_feature'};
}

sub heavy_feature {
	my $self = shift;
	return $self->{'heavy_feature'};
}

1;
