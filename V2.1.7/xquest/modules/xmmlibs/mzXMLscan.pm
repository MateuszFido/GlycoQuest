package mzXMLscan;
use strict;
#---------------------------------------------------------------------------
# Module: mzXMLscan.pm
# Author(s): Thomas Walzthoeni, xQuest specific modules are based on 
# original work by Oliver Rinner.
# Description: Module for mzXML scanobjects.
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
use XML::TreeBuilder;
use Base64;
use Consensus;

my $Hatom = 1.0078250;

sub new_pseudoscan{
	my $class    = shift();
	my $scan     = shift;  # is a feature obj in this case
	my $id       = shift;
	my $basename = shift;
	my $getpeaks = shift;
	
	my $self = {};

	bless $self, $class;	
	my $scannumber = $scan->get_pseudoSHscannum;
	my $retentiontime = $scan->Tr;
	my $mz = $scan->mz;

	if ($getpeaks) {
		#$self->{'peaks'} = $scan->find('peaks')->content->[0];
	}
	my $basePeakIntensity = 1;
	$self->{'filterline'} = '';
	
	my $charge = $scan->charge;
	my $precursorIntensity = $scan->mz;
	my $activationMethod = "not defined -pseudoSH";

	$self->{'scannumber'}         = $scannumber;
	$self->{'basename'}           = $basename;
	$self->{'basePeakIntensity'}  = $basePeakIntensity;
	$self->{'FT_measured_mz'}     = $mz;
	$self->{'FT_measured_charge'} = $charge;
	$self->{'Tr'}                 = $retentiontime;
	$self->{'runID'}              = $id;
	$self->{'msInstrumentID'}     = "not set -pseudoSH";
	$self->{'precursorIntensity'} = $scan->area;

	return $self;
	
	
}

sub new {
	my $class    = shift();
	my $scan     = shift;
	my $id       = shift;
	my $basename = shift;
	my $getpeaks = shift;

	my $self = {};

	bless $self, $class;

	my $scannumber = $scan->attr('num');
	( my $retentiontime = $scan->attr('retentionTime') ) =~ s/PT(.*)S/$1/g;
	#print $retentiontime."\n";
	
	my $mz = $scan->find('precursorMz')->content->[0];

	if ($getpeaks) {
		$self->{'peaks'} = $scan->find('peaks')->content->[0];
	}

	my $basePeakIntensity = $scan->attr('basePeakIntensity');
	if ( $scan->attr('filterLine') ) {
		$self->{'filterline'} = $scan->attr('filterLine');
	}
	else {
		$self->{'filterline'} = '';
	}
	my $charge = $scan->find('precursorMz')->attr('precursorCharge');
	my $precursorIntensity =
	$scan->find('precursorMz')->attr('precursorIntensity');
	my $activationMethod = $scan->find('precursorMz')->attr('activationMethod');

	$self->{'scannumber'}         = $scannumber;
	$self->{'basename'}           = $basename;
	$self->{'basePeakIntensity'}  = $basePeakIntensity;
	$self->{'FT_measured_mz'}     = $mz;
	$self->{'FT_measured_charge'} = $charge;
	$self->{'Tr'}                 = $retentiontime;
	$self->{'runID'}              = $id;
	$self->{'msInstrumentID'}     = $scan->attr('msInstrumentID');
	$self->{'precursorIntensity'} = $precursorIntensity;

	return $self;
}

sub basename {
	my $self = shift;
	return $self->{'basename'};
}

sub runID {
	my $self = shift;
	return $self->{'runID'};
}

sub scannumber {
	my $self = shift;
	return $self->{'scannumber'};
}

sub basePeakIntensity {
	my $self = shift;
	return $self->{'basePeakIntensity'};
}

sub precursorIntensity {
	my $self = shift;
	return $self->{'precursorIntensity'};
}

sub filterline {
	my $self = shift;
	return $self->{'filterline'};
}

sub fractionationtype {
	my $self       = shift;
	my $filterline = $self->filterline;
	if ( $filterline =~ /\@etd\d+/ ) {
		return 'etd';
	}
	elsif ( $filterline =~ /\@cid\d+/ ) {
		return 'cid';
	}
	else {
		return 'nd';
	}
}

sub get_peaks_inbase64 {
	my $self = shift;
	return $self->{'peaks'};
}

sub print_dta {
	my $self          = shift;
	my $outfilehandle = shift;
	my $mz            = shift;
	my $charge        = shift;

	unless ($outfilehandle) {
		$outfilehandle = *STDOUT;
	}
	unless ($mz) {
		$mz = $self->FT_mz;
	}
	unless ($charge) {
		$charge = $self->FT_charge;
	}

	my $Mr    = $self->Mr;
	my $peaks = $self->get_peaklist;

	print $outfilehandle "$Mr\t$charge\n";
	foreach my $peakpair (@$peaks) {
		print $outfilehandle $peakpair->[0], "\t", $peakpair->[1], "\n";
	}
}

sub print_into_mgf {
	my $self          = shift;
	my $outfilehandle = shift;
	my $mz            = shift;
	my $charge        = shift;
	my $id=shift;
	my $peaks = $self->get_peaklist;

	print $outfilehandle "BEGIN IONS\n";
			print $outfilehandle "TITLE=", $id, "\n";
			print $outfilehandle "CHARGE=$charge+\n";
			print $outfilehandle "PEPMASS=$mz\n";

	foreach my $peakpair (@$peaks) {
		print $outfilehandle $peakpair->[0], "\t", $peakpair->[1], "\n";
	}
	print $outfilehandle "END IONS\n";
	
}


sub Mr {
	my $self   = shift;
	my $charge = $self->FT_charge;

	return $self->FT_mz * $charge - $charge * $Hatom;
}

sub get_peaklist {
	my $self = shift;

	my $base64string  = $self->get_peaks_inbase64;
	my $base64decoded = Base64::b64decode($base64string);
	my @hostOrder32   = unpack( "N*", $base64decoded );

	# unpack the binary data as host ordered 32 bit ints
	#foreach my $i (@hostOrder32) {
	#	my $float = unpack( "f", pack( "I", $i ) );
	#	print $float, " ";
	#	# The hostOrder32 array contains a list of
	#	# host ordered 32 bits entities which we want to re-interpret
	#	# as floats. In Perl this means we have to
	#	# pack it back as an int and then unpack it as a float
	#	# This would all have been simpler if only Perl
	#	# had had a network/host order option on unpack float
	#	# But we don't so alas we do the ordering operation
	#	# in the first unpack (N*) and then do the conversion
	#	# to float in the second
	#}
### or to clearly distinguish m/z from intensity:
	#my @mzs         = ();
	#my @intensities = ();
	my @peakpairs  = ();
	my $done       = 0;
	my $npeakpairs = scalar(@hostOrder32);

	#print "npeak=$npeakpairs\n";
	my $i = 0;
	for ( $i = 0 ; $i < $npeakpairs ; $i += 2 ) {
		my $mz        = unpack( "f", pack( "I", shift(@hostOrder32) ) );
		my $intensity = unpack( "f", pack( "I", shift(@hostOrder32) ) );
		push @peakpairs, [ $mz, $intensity ];
	}
	return \@peakpairs;
}

sub get_binned_normalized_peaks {
	my $self     = shift;
	my $peaklist = $self->get_peaklist;
	my $binned_peaks=Consensus::get_binned_peaks($peaklist);	
	return $binned_peaks;	
}

sub id {
	my $self = shift;
	return join "_", $self->basename, $self->scannumber;
}

sub FT_mz {
	my $self = shift;
	return $self->{'FT_measured_mz'};
}

sub precursor_mz {
	my $self = shift;
	return $self->{'FT_measured_mz'};
}

sub FT_charge {
	my $self   = shift;
	my $charge = $self->{'FT_measured_charge'};
	if ($charge) {
		return $self->{'FT_measured_charge'};
	}
	else {
		return 1;
	}
}

sub msInstrumentID {
	my $self = shift;
	return $self->{'msInstrumentID'};
}

sub Tr {
	my $self = shift;
	return $self->{'Tr'};
}

sub assign_by_basepeakIntensity {
	my $feature           = shift;
	my $scanarray         = shift;
	my @sortedbyintensity =
	  sort { $b->basePeakIntensity <=> $a->basePeakIntensity } @$scanarray;
	return $sortedbyintensity[0];
}

sub assign_by_distance2apex {
	my $feature         = shift;
	my $scanarray       = shift;
	my @sortedbydeltaTr =
	  sort { abs( $a->Tr - $feature->Tr ) <=> abs( $b->Tr - $feature->Tr ) }
	  @$scanarray;
	return $sortedbydeltaTr[0];
}

1;
