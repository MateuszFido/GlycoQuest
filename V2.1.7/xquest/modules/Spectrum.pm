package Spectrum;
use strict;
#---------------------------------------------------------------------------
# Module: Spectrum.pm
# Author(s): Thomas Walzthoeni, xQuest specific modules are based on 
# original work by Oliver Rinner.
# Description: Module for handling MS/MS spectra.
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
use Statistics;
use File::Basename;
use File::Copy;
use Xcorrelation;
use MIME::Base64;
use Data::Dumper;

sub new {
	my $class = shift();
	my $self  = {};
	bless $self, $class;
	my $spectrum     = shift;
	my $PARAMS       = shift;
	my $scantype     = shift;
	my $isotopeshift = shift;
	my $spechash	 = shift;
	$self->{'PARAMS'}       = $PARAMS;
	$self->{'scantype'}     = $scantype;
	$self->{'addedmass'}    = $isotopeshift;
	$self->{'specbasename'} = basename($spectrum);

	## common an xlinkspectrum
	my $commonspectrum = $self->{'commonspecname'} = join "", $spectrum,"_common.txt";
	my $xlinkspectrum = $self->{'xlinkspecname'} = join "", $spectrum,"_xlinker.txt";
	
	### Parsing the Spectra
	my %mz_intensityhash;
	my ( $intensities_common, $intensities_xlink );
	
	if ($PARAMS->{'specxml'}){
	#	print "SPECXML";
	print "Reading Spectra common ion Spectrum $commonspectrum from spectrum Hash";
	( $self->{'peaks'}, $self->{'ionpairs'}, $intensities_common ) = $self->_readspec_from_hash( $commonspectrum, $PARAMS, \%mz_intensityhash, $spechash);
	print "...done\n";
	print "Reading cross-link ion Spectrum $xlinkspectrum from spectrum Hash";
	( $self->{'xlinkerpeaks'}, $self->{'xlinkerionpairs'}, $intensities_xlink )  = $self->_readspec_from_hash( $xlinkspectrum, $PARAMS, \%mz_intensityhash, $spechash );		
	print "...done\n";
	}else{
	( $self->{'peaks'}, $self->{'ionpairs'}, $intensities_common ) = $self->_readspec( $commonspectrum, $PARAMS, \%mz_intensityhash );
	( $self->{'xlinkerpeaks'}, $self->{'xlinkerionpairs'}, $intensities_xlink )  = $self->_readspec( $xlinkspectrum, $PARAMS, \%mz_intensityhash );
	}

	
	unless ($PARAMS->{'specxml'}){
	if ( $PARAMS->{'outputpath'} ) {
	
		my $directory = $PARAMS->{'outputpath'};
		#print "COPY SPECTRUM\n";
		copy( $commonspectrum, $directory )
		  or warn
		  "copyying of $commonspectrum to $directory was not successfull: $!";
		copy( $xlinkspectrum, $directory )
		  or warn
		  "copyying of $commonspectrum to $directory was not successfull: $!";
	}
	}

	my @intensities=map {$mz_intensityhash{$_}} keys %mz_intensityhash;
#	print Dumper ($self->getcommonpairs);
	my @mzcommon=map { mz($_)} @{$self->getcommonpairs};
	my @mzxlink=map { mz($_)} @{$self->getxlinkpairs};
	#print Dumper (\@mzxlink);
	my @allmz= keys %mz_intensityhash;
	
	$self->{'mean_mz'}  = Statistics::mean( \@allmz );
	$self->{'mean_mz_common'}  = Statistics::mean( \@mzcommon );
	$self->{'mean_mz_xlink'}  = Statistics::mean( \@mzxlink );
	#print "M common:",$self->{'mean_mz_common'};
	#print "M xl:",$self->{'mean_mz_xlink'};
	#print "M all:",$self->{'mean_mz'};
	#print Dumper (\@allmz);
	$self->{'chargesorted_xlinkpeaks'} = $self->sortbycharge( $self->{'xlinkerionpairs'} );
	$self->{'chargesorted_xlinkpairs'} = $self->sort_ionpairs_bycharge( $self->{'xlinkerionpairs'} );

	$self->{'ionintensityhash'} = \%mz_intensityhash;
	$self->{'ncommonions'}      = $#{ $self->{'peaks'} } + 1;
	$self->{'nxlinkions'}       = $#{ $self->{'xlinkerpeaks'} } + 1;
	$self->{'ionintensitylist'}   = \@intensities;
	$self->{'mean_ionintensity'}  = Statistics::mean( \@intensities );
	$self->{'total_ionintensity'} = Statistics::sum( \@intensities );
	$self->{'ionintensity_stdev'} = Statistics::stdev( \@intensities );
	$self->{'ionintensity_sem'}   = Statistics::sem( \@intensities );
	$self->{'ionintensity_max'}   = Statistics::max( \@intensities );

	if ( $PARAMS->{'Iontagmode'} ) {
		$self->{'testions'} = $self->_gettestions($PARAMS);
	}
	if ( $PARAMS->{'normxcorr'} ) {
		$self->autocorrelation;
	}
	else {
		$self->{'autocorrxlink'}    = 1;
		$self->{'autocorrbackbone'} = 1;
		$self->{'autocorrall'}      = 1;
	}

	$self->calc_apriori_probabilities;
	#print "<br>Spectrum Summary:<br>";
	#print "Number of common ions:",$self->getncommonions,"<br>";
	#print "Number of xlink ions:",$self->getnxlinkions,"<br>";
	#print "Apriory pxlink:",$self->get_apriori_pxlink,"<br>";
	#print "Apriory pcommon:",$self->get_apriori_pcommon,"<br>";
	#print "COMMONPEAKS<br>";
	#print Dumper ($self->{'peaks'});
	#print "XL-PEAKS<br>";
	#print Dumper ($self->{'xlinkerpeaks'});
	return $self;
}

sub calc_apriori_probabilities {
	my $self     = shift;
	my $PARAMS   = $self->getParams;
	my $scantype = $self->get_scantype;

	my $ncommonpeaks = $self->getncommonions;
	my $nxlinkpeaks  = $self->getnxlinkions;

	unless ($nxlinkpeaks) {
		$nxlinkpeaks = 1;
	}
	unless ($ncommonpeaks) {
		$ncommonpeaks = 1;
	}
	my $ms2precision       = $PARAMS->{'ms2tolerance'};
	my $xlink_ms2tolerance = $PARAMS->{'xlink_ms2tolerance'};
	
	my $range = $PARAMS->{'maxionsize'} - $PARAMS->{'minionsize'};
	if ( $PARAMS->{'tolerancemeasure_ms2'} =~ /^ppm/i )
	{
		#$tolerance = $tolerance * 1e-6 * 1000;    #ppm to amu measure
		## take the mean of the range
		###my $meanmzcommon = $self->{'mean_mz_common'};
		### Calc the mean m/z of all peaks
		
		
		$ms2precision       = $self->{'mean_mz_common'} * 1e-6 * $PARAMS->{'ms2tolerance'};
		$xlink_ms2tolerance = $self->{'mean_mz_xlink'} * 1e-6 * $PARAMS->{'xlink_ms2tolerance'};
		#print "Mean mz:$self->{'mean_mz'}<br>";
		#print "Precision ms2:".$ms2precision."<br>";
		#print "Precision xlms2:".$xlink_ms2tolerance."<br>";
	}
#print "MS2precision:$ms2precision\n";

	
	$self->{'apriori_pcommon'} = ( 1 - ( ( 1 - 2 * $ms2precision / ( 0.5 * $range ) )**($ncommonpeaks) ) );

	if ( $scantype eq 'light_heavy' ) {
		## maxioncharge is the precursorion charge -1
		my $nxlinkcharges =  $self->maxioncharge_xlinks - $self->minioncharge_xlinks;

		unless ($nxlinkcharges) {
			$nxlinkcharges = 1;
		}
		$self->{'apriori_pxlink'} = (
			1 - (( 1 - 2 * $xlink_ms2tolerance / ( 0.5 * $range ) ) **( $nxlinkpeaks / $nxlinkcharges )	));
	}
	else {
		$self->{'apriori_pxlink'} = (
			1 - (( 1 - 2 * $xlink_ms2tolerance / ( 0.5 * $range ) )**($nxlinkpeaks)));

	}
	my $type;
	#print Dumper ($scantype);
	#print "apriori probabilities for matching type  pcommon ",$self->{'apriori_pcommon'}," pxlink: " ,$self->{'apriori_pxlink'},"\n";

}

sub get_scantype {
	my $self = shift;
	return $self->{'scantype'};
}

sub get_addedmass {
	my $self = shift;
	return $self->{'addedmass'};
}

sub get_apriori_pcommon {
	my $self = shift;
	return $self->{'apriori_pcommon'};
}

sub get_apriori_pxlink {
	my $self = shift;
	return $self->{'apriori_pxlink'};
}

sub maxioncharge_xlinks {
	my $self                = shift;
	my $definedmaxioncharge = $self->getParams->{'ioncharge_xlink'}->[-1];
	my $spectrumcharge      = $self->getprecursorCharge;
	if ( $definedmaxioncharge > $self->getprecursorCharge ) {
		return $spectrumcharge;
	}
	else {
		return $definedmaxioncharge;
	}
}

sub minioncharge_xlinks {
	my $self = shift;
	return $self->getParams->{'ioncharge_xlink'}->[0];
}

sub basepeakintensity {
	my $self = shift;
	return $self->{'ionintensity_max'};
}

sub sortbycharge {
	my $self      = shift;
	my $peakarray = shift;
	my %sortedchargehash;

	foreach my $peak (@$peakarray) {
		push @{ $sortedchargehash{ charge($peak) } }, mz($peak);

		#	print "peak: ",mz($peak)," charge: ",charge($peak),"\n";

	}
	return \%sortedchargehash;
}

sub sort_ionpairs_bycharge {
	my $self      = shift;
	my $peakarray = shift;
	my %sortedchargehash;
	foreach my $peak (@$peakarray) {
		push @{ $sortedchargehash{ charge($peak) } }, $peak;
	}
	return \%sortedchargehash;
}

sub mz {
	my $pairs = shift;
	return $pairs->[0];
}

sub getintensities {
	my $self = shift;
	return $self->{'ionintensitylist'};
}

sub intensity {
	my $pairs = shift;
	return $pairs->[1];
}

sub charge {
	my $pairs = shift;
	return $pairs->[2];
}

sub newmasslist {
	my $class = shift();
	my $self  = {};
	bless $self, $class;
	my $PARAMS    = shift;
	my $massvalue = shift;
	my $charge    = shift;

	$self->{'PARAMS'} = $PARAMS;
	if ( $massvalue && $charge ) {
		$self->_readmassonly_mz_charge( $massvalue, $charge );
	}
	elsif ($massvalue) {
		$self->_readmassonly($massvalue);
	}
	else {
		die "Mr value or mz + charge is required $0";
	}

	return $self;
}

sub getParams {
	my $self = shift;
	return $self->{'PARAMS'};
}

sub autocorrelation {
	my $self        = shift;
	my @commonpeaks = @{ $self->getcommonpeaks };
	my @commonpairs = @{ $self->getcommonpairs };
	my @xlinkpeaks  = @{ $self->getxlinkpeaks };
	my @xlinkpairs  = @{ $self->getxlinkpairs };
	my $PARAMS      = $self->getParams;
	my $maxionsize  = $self->getParams->{'maxionsize'};

	my $printout     = $PARAMS->{'printtables'};
	my $poolisotopes = $PARAMS->{'poolisotopes'};
	my $userealint   = $PARAMS->{'realintensities4xcorr'};

	my ( $i, $j, $k, @commonpeaktable, @xlinkpeaktable );
	### check if ppm measure is used ##added by TW
	my $tolerance;
	if ( $PARAMS->{'xcorrprecision'} )
	{
		$tolerance = $PARAMS->{'xcorrprecision'};
	} else
	{
		$tolerance = $PARAMS->{'ms2tolerance'};
	}

	my $precision = int( 1 / $tolerance);
	#print "Precision for bins: $precision\n";
#	my $precision = int( 1 / $PARAMS->{'ms2tolerance'} );

	my $delay = int( $PARAMS->{'xcorrdelay'} );

	my $weigthxlinker = $PARAMS->{'xlinkxcorrweigth'};
	my $weigthcommon  = $PARAMS->{'commonxcorrweigth'};
	#print "weight xl: $weigthxlinker\n";

	#todo: make reasonable matching function taking intensities into account
	if ($userealint) {
		foreach (@xlinkpairs) {
			$xlinkpeaktable[ int( $precision * $_->[0] ) ] = $_->[1];
		}

		foreach (@commonpairs) {
			$commonpeaktable[ int( $precision * $_->[0] ) ] = $_->[1];
		}
	}
	else {
		### transform the peaks to an array where the keys of the peaks are set to fixed values
		### the other slots are undef
		foreach (@xlinkpeaks) {
			#print "INT VALUE:".( int( $precision * $_ ))."\n";
			$xlinkpeaktable[ int( $precision * $_ ) ] = $weigthxlinker;
		}
		foreach (@commonpeaks) {
			$commonpeaktable[ int( $precision * $_ ) ] = $weigthcommon;
		}
	}

	my $rcommon =  Xcorrelation::calcxcorr( \@commonpeaktable, \@commonpeaktable, int($delay) );
	my $rxlink = Xcorrelation::calcxcorr( \@xlinkpeaktable, \@xlinkpeaktable, int($delay) );

	my ( $bbxcorr, $xlinkxcorr );

	if ( $delay == 0 ) {
		$bbxcorr    = $rcommon->[0]->[1];
		$xlinkxcorr = $rxlink->[0]->[1];
	}
	else {
		### if a delay was used then check out the value at delay 0
		for $i ( 0 .. $#$rcommon ) {
			if ( $rcommon->[$i]->[0] == 0 ) {
				$bbxcorr +=
				  $rcommon->[$i]->[1] + $rcommon->[ $i - 1 ]->[1] +
				  $rcommon->[ $i + 1 ]->[1];
			}
		}
		for $i ( 0 .. $#$rxlink ) {
			if ( $rxlink->[$i]->[0] == 0 ) {
				$xlinkxcorr +=
				  $rxlink->[$i]->[1] + $rxlink->[ $i - 1 ]->[1] +
				  $rxlink->[ $i + 1 ]->[1];
			}
		}
	}

	#print" xlink: $xlinkxcorr backbone $bbxcorr\n";
	$self->{'autocorrxlink'}    = $xlinkxcorr;
	$self->{'autocorrbackbone'} = $bbxcorr;

}

#sub _calcxcorr {
#	my $series1  = shift;
#	my $series2  = shift;
#	my $maxdelay = shift;
#	my @x        = @{$series1};
#	my @y        = @{$series2};
#	my $n        = $#$series1;
#	my ( $i, $j, $delay, $mx, $my, $sx, $sy, $sxy, $denom, @r );
#
#	#   /* Calculate the mean of the two series x[], y[] */
#	$my = 0;
#	$my = 0;
#	for ( $i = 0 ; $i <= $n ; $i++ ) {
#		$mx += $x[$i];
#		$my += $y[$i];
#	}
#
#	#print "n $n mean x $mx mean y =$my\n";
#	unless ($n) {
#		$n = 1;
#	}
#	$mx /= ($n);
#	$my /= ($n);
#
#	#   /* Calculate the denominator */
#	$sx = 0;
#	$sy = 0;
#	for ( $i = 0 ; $i <= $n ; $i++ ) {
#		$sx += ( $x[$i] - $mx ) * ( $x[$i] - $mx );
#		$sy += ( $y[$i] - $my ) * ( $y[$i] - $my );
#	}
#
#	$denom = sqrt( $sx * $sy );
#
#	#   /* Calculate the correlation series */
#	for ( $delay = -$maxdelay ; $delay < $maxdelay ; $delay++ ) {
#		$sxy = 0;
#		for ( $i = 0 ; $i < $n ; $i++ ) {
#			$j = $i + $delay;
#
#			unless ( $j < 0 || $j >= $n ) {
#				$sxy += ( $x[$i] - $mx ) * ( $y[$j] - $my );
#			}
#			else {
#				next;
#			}
#
#		}
#		if ( $denom != 0 ) {
#			push @r,
#			  [ $delay, ( $sxy / $denom )
#			  ]    #      /* r is the correlation coefficient at "delay" */
#			       #print "$delay ",$sxy / $denom,"\n";
#		}
#		else {
#			push @r, [ $delay,
#				0 ]    #      /* r is the correlation coefficient at "delay" */
#			           #print "$delay ",$sxy / $denom,"\n";
#		}
#	}
#
#	return \@r;
#}

sub autoxlink {
	my $self = shift;
	return $self->{'autocorrxlink'};

}

sub autobackbonecorr {
	my $self = shift;
	return $self->{'autocorrbackbone'};

}

sub autoallcorr {
	my $self = shift;
	return $self->{'autocorrall'};
}

sub getspecbasename {
	my $self = shift;
	return $self->{'specbasename'};
}

sub getncommonions {
	my $self = shift;
	return $self->{'ncommonions'};
}

sub getnxlinkions {
	my $self = shift;
	return $self->{'nxlinkions'};
}

sub DESTROY {
	my $this = shift;
}

sub _readspec {
	my $self             = shift;
	my $spectrum         = shift;
	my $PARAMS           = shift;
	my $mz_intensityhash = shift;
#print "Reading: $spectrum<br>";
	#my $minintensity     = $PARAMS->{'minionintensity'};
	my ( @peaks, @pairs, @intensities );
	open SPEC, "<$spectrum" or die "cannot open spectrum file $spectrum $!";
	my $header = <SPEC>;
	chomp($header);
	$self->{'header'} = $header;
	my $precursorMz = <SPEC>;
	chomp($precursorMz);
	$self->{'precursorMz'} = $precursorMz;

	#print "Mz: $precursorMz\n";

	my $precursorCharge = <SPEC>;
	chomp($precursorCharge);
	$self->{'precursorCharge'} = $precursorCharge;

	#print "Charge: $precursorCharge\n";

	while ( <SPEC>) {
#		chomp ($line);
#		unless ($line){
#		next;	
#		}
		
		my @tmp       = split;
		#print @tmp,"<br>";
		unless (@tmp){
			next;
		}
		my $mz        = $tmp[0];
		my $intensity = $tmp[1];
		my $charge    = $tmp[2];

		push @peaks, $mz;
		push @pairs, [ $mz, $intensity, $charge ];
#		push @intensities, sqrt($intensity);
#		$mz_intensityhash->{$mz} = sqrt($intensity);
		push @intensities, $intensity;
		$mz_intensityhash->{$mz} = $intensity;

	}

	close(SPEC);
	return \@peaks, \@pairs, \@intensities;

}

sub _readspec_from_hash {
	my $self             = shift;
	my $spectrum         = shift; ### must be the key in the spechash
	my $PARAMS           = shift;
	my $mz_intensityhash = shift;
	my $spechash=shift;

	## use the basename of the spectrum as key
	my $key=basename($spectrum);
	my ( @peaks, @pairs, @intensities );
	my $speccontent=$spechash->{$key};
	my $pl       = decode_base64($speccontent);
	
	#print Dumper ($pl);
	#exit;
	unless($pl){
	die "Empty Peaklist, cant proceed $!";	
	}
	
	my @plarray=split (/\n/,$pl);

	my $header = $plarray[0];
	chomp($header);
	$self->{'header'} = $header;
	my $precursorMz = $plarray[1];
	chomp($precursorMz);
	$self->{'precursorMz'} = $precursorMz;
#exit;
	#print "Mz: $precursorMz\n";
	#die;
	my $precursorCharge = $plarray[2];
	chomp($precursorCharge);
	
	$self->{'precursorCharge'} = $precursorCharge;

	#print "Charge: $precursorCharge\n";
	#my $lastelement=@#plarray;

	foreach my $line( 3 .. $#plarray ) {
		my $tmpline=$plarray[$line];
		chomp ($tmpline);
		unless ($tmpline){
		print "Warning: Empty line in spectrum found\n";
		next;	
		}
		#print "Index: $line: ".$tmpline."\n";		
		my @tmp       = split (/\t/,$tmpline);
		my $mz        = $tmp[0];
		my $intensity = $tmp[1];
		my $charge    = $tmp[2];
		#print "read mz: $mz, int: $intensity, z: $charge\n";
		
		push @peaks, $mz;
		push @pairs, [ $mz, $intensity, $charge ];
#		push @intensities, sqrt($intensity);
#		$mz_intensityhash->{$mz} = sqrt($intensity);
		push @intensities, $intensity;
		$mz_intensityhash->{$mz} = $intensity;

	}

	return \@peaks, \@pairs, \@intensities;

}







sub _readmassonly {
	my $self = shift;
	my $Mr   = shift;
	unless ( $Mr =~ /\d+\.\d+/ ) {
		die
"mass must be indicated in floating point numer e.g. 2343.231 not $Mr $0";
	}
	$self->{'precursorMz'}     = $Mr + 1.007825032;
	$self->{'precursorCharge'} = 1;
}

sub _readmassonly_mz_charge {
	my $self   = shift;
	my $mz     = shift;
	my $charge = shift;
	unless ( $mz =~ /\d+\.\d+/ && $charge =~ /\d/ ) {
		die
"mz must be indicated in floating point numer e.g. 2343.231 not $mz $0";
	}
	$self->{'precursorMz'}     = $mz;
	$self->{'precursorCharge'} = $charge;
}

sub _gettestions {
	my $self   = shift;
	my $PARAMS = shift;
	my $method;
	my ( $i, $j, @parentlist );

	if ( defined( $PARAMS->{'testionspick'} ) ) {
		$method = $PARAMS->{'testionspick'};
	}

	else { die "de novo method not defined (mzdescend|mzascend|intensity) $!"; }
	if ( $PARAMS->{'testionspick'} eq 'mzdescend' ) {
		my $peaks     = $self->getcommonpeaks;
		my $ntestions = 0;

		my $nions = $PARAMS->{'ntestions'};
		for ( $i = $#$peaks ; $i >= 0 && ( $#$peaks - $i ) < $nions ; $i-- )
		{    #pick peaks in descending m/z order
			if (   ( $peaks->[$i] >= $PARAMS->{'minionsize'} )
				&& ( $peaks->[$i] <= $PARAMS->{'maxionsize'} ) )
			{
				$ntestions++;
				push @parentlist, $peaks->[$i];
			}
		}

	}
	elsif ( $PARAMS->{'testionspick'} eq 'mzascend' ) {
		my $peaks     = $self->getcommonpeaks;
		my $nions     = $PARAMS->{'ntestions'};
		my $ntestions = 0;
		for ( $i = 0 ; $i <= $#$peaks && ( $i < $nions ) ; $i++ )
		{    #pick peaks in descending m/z order
			if (   ( $peaks->[$i] >= $PARAMS->{'minionsize'} )
				&& ( $peaks->[$i] <= $PARAMS->{'maxionsize'} ) )
			{
				$ntestions++;
				push @parentlist, $peaks->[$i];
			}
		}

	}
	elsif ( $PARAMS->{'testionspick'} eq 'intensity' ) {
		my $ntestions = 0;
		my $getnions  = $PARAMS->{'ntestions'};

		my @sortedpeaks = sort { $b->[1] <=> $a->[1] } @{ $self->getcommonpairs };
	
		foreach my $peak (@sortedpeaks) {
			if (   ( $peak->[0] >= $PARAMS->{'minionsize'} ) && ( $peak->[0] <= $PARAMS->{'maxionsize'} ) )
			{
				$ntestions++;
				push @parentlist, $peak->[0];
				#print $peak->[0], "\t", $peak->[1], "\n";
			}
			if ( $ntestions >= $getnions ) { last; }
		}
	}
	else {
		die "test ion peakpicking ", $PARAMS->{'testionspick'},
		  " not defined $!";
	}

	return \@parentlist;
}

sub getcommonpairs {
	my $self = shift;
	return $self->{'ionpairs'};
}

sub getxlinkpairs {
	my $self   = shift;
	my $charge = shift;
	if ( defined($charge) ) {
		if ( $self->{'chargesorted_xlinkpairs'}->{$charge} ) {
			return $self->{'chargesorted_xlinkpairs'}->{$charge};
		}
		else { return []; }
	}
	else {
		return $self->{'xlinkerionpairs'};
	}
}

sub getxlinkpeaks {
	my $self   = shift;
	my $charge = shift;
	if ( defined($charge) ) {
		if ( $self->{'chargesorted_xlinkpeaks'}->{$charge} ) {
			return $self->{'chargesorted_xlinkpeaks'}->{$charge};
		}
		else {
			return [];
		}
	}
	else {
		return $self->{'xlinkerpeaks'};
	}
}

sub get_ionintensity {
	my $self = shift;
	my $mz   = shift;
	return $self->get_ionintensityhash->{$mz};
}

#sub get_realionintensity {
#	my $self = shift;
#	my $mz   = shift;
#	return int( $self->get_ionintensityhash->{$mz}**2 );
#}

sub get_ionintensityhash {
	my $self = shift;
	return $self->{'ionintensityhash'};
}

sub get_average_ionintensity {
	my $self = shift;
	return sprintf( "%.2f", $self->{'mean_ionintensity'} );
}

sub get_total_ionintensity {
	my $self = shift;
	return sprintf( "%.2f", $self->{'total_ionintensity'} );
}

sub get_ionintensity_stdev {
	my $self = shift;
	return sprintf( "%.2f", $self->{'ionintensity_stdev'} );
}

sub get_ionintensity_sem {
	my $self = shift;
	return sprintf( "%.2f", $self->{'ionintensity_sem'} );
}

sub getSpecname {
	my $self = shift;
	return $self->{'specbasename'};
}

sub getcommonpeaks {
	my $self = shift;
	return $self->{'peaks'};
}

sub gettestions {
	my $self = shift;
	return $self->{'testions'};
}

sub getprecursorMz {
	my $self = shift;
	if ( defined( $self->{'precursorMz'} ) ) {
		return $self->{'precursorMz'};
	}
	else { warn "precursorMz is not defined, check spectrumheader!\n"; }
}

sub getspectrumheader{
	my $self = shift;
	if ( defined( $self->{'header'} ) ) {
		return $self->{'header'};
	}
	else { warn "header is not defined, check spectrum header!\n"; }
}


sub getprecursorCharge {
	my $self = shift;
	if ( defined( $self->{'precursorCharge'} ) ) {
		return $self->{'precursorCharge'};
	}
	else { warn "precursorCharge is not defined, check spectrum header!\n"; }
}

sub getms1mass {
	my $self    = shift;
	my $ms1mass = $self->getprecursorMz * $self->getprecursorCharge -
	  $self->getprecursorCharge * 1.007825032;

#print "getms1mass: Mz ",$self->getprecursorMz,"charge ",$self->getprecursorCharge," mass $ms1mass\n";
	return $ms1mass;
}

1;
