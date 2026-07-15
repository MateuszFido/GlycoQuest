package Xcorrelation;
use strict;
#---------------------------------------------------------------------------
# Module: Xcorrelation.pm
# Author(s): Thomas Walzthoeni, xQuest specific modules are based on 
# original work by Oliver Rinner.
# Description: Module for cross-correlation.
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
use Data::Dumper;

sub xcorrelation_common
{
	my $self            = shift;
	my $basefilename    = shift;
	my $printout        = shift;
	my $delay           = shift;
	my $precision       = shift;
	my $lossions        = shift;
	my $spectrum        = $self->getSpecObj;
	my @commonpeaks     = @{ $spectrum->getcommonpeaks };
	my @commonpairs     = @{ $spectrum->getcommonpairs };
	my @commonionstheor = @{ $self->getCommonIons };
	my $PARAMS          = $self->getParams;
	my $maxionsize      = $self->getParams->{'maxionsize'};
	unless ($printout)
	{
		$printout = $self->getParams->{'printtables'};
	}
	my $userealint = $self->getParams->{'realintensities4xcorr'};
	my $id         = $self->getid;
	$id =~ s/::/_/g;
	my $commoncorrfile = join "", $basefilename, "_commonxcorr.txt";
	$printout
	  && ( open XCOMMON, ">$commoncorrfile" or die "xcorrelation:cannot open $commoncorrfile\n" );
	my ( $i, $j, $k, @commoniontable, @commonpeaktable, );
	my $tolerance;

	if ( $PARAMS->{'xcorrprecision'} )
	{
		$tolerance = $PARAMS->{'xcorrprecision'};
	} else
	{
		$tolerance = $PARAMS->{'ms2tolerance'};
	}
	if ( $PARAMS->{'tolerancemeasure_ms2'} =~ /^ppm/i )
	{
		#$tolerance = $tolerance * 1e-6 * 1000;    #ppm to amu measure
	}
	#print "Tolerance: $tolerance<br>";
	#exit;
	unless ($precision)
	{
		$precision = int( 1 / $tolerance );
	}

	unless ($delay)
	{
		$delay = int( $PARAMS->{'xcorrdelay'} );
	}
	my $weigthcommon = $PARAMS->{'commonxcorrweigth'};
	for $i ( 0 .. $maxionsize * $precision )
	{
		$commoniontable[$i] = $commonpeaktable[$i] = 0;
	}
	#print "Delay: $delay<br>";
	#todo: make reasonable matching function taking intensities into account
	if ($userealint)
	{
		foreach (@commonpairs)
		{
			$commonpeaktable[ int( $precision * $_->[0] ) ] = $_->[1];
		}
	} else
	{
		foreach (@commonpeaks)
		{
			$commonpeaktable[ int( $precision * $_ ) ] = $weigthcommon;
		}
	}
	foreach (@commonionstheor)
	{
		$commoniontable[ int( $precision * $_ ) ] = $weigthcommon;
	}
##old
	my $rcommon = calcxcorr( \@commonpeaktable, \@commoniontable, int($delay) );
	my $bbxcorr = 0;

	#print Dumper ($rcommon);
#### ATTENTION: THE OLD XCORRELATION TOOK ALWAYS THE SUM OF ALL XCORR VALUES 
#### AROUND THE SPECIFIED DELAY!
#### THE DELAY IS USUALLY DEFINED IN THE XQUEST DEFFILE (5)
#### THE WEBINTERFACE REDEFINES THE DELAYS.
	
	for $i ( 0 .. $#$rcommon )
	{
		$bbxcorr += $rcommon->[$i]->[1];
		$printout && print XCOMMON $rcommon->[$i]->[0] / $precision, "\t", $rcommon->[$i]->[1], "\n";
	}

	$printout && drawxcorrgraph( $rcommon, $precision, ( join "_", $basefilename, "backbonecorr.png" ), "xcorr backbone-ions" );
	$printout && close(XCOMMON);

	#old: normalize with autocorrelation
	$self->{'backbonexcorr'} = $bbxcorr;
	$self->getSpecObj->autobackbonecorr  && ( $self->{'backbonexcorr'} = $bbxcorr / $self->getSpecObj->autobackbonecorr );
	#print "AUTOCORR: ",$self->getSpecObj->autobackbonecorr,"<br>";
	#print "<br>Xcorr Common ions (static): ",$self->{'backbonexcorr'} ,"<br>";

	#new: use normalized crosscorrelation
}

sub xcorrelation_common_normalized
{
	my $self            = shift;
	my $basefilename    = shift;
	my $printout        = shift;
	my $delay           = shift;
	my $precision       = shift;
	my $lossions        = shift;
	my $spectrum        = $self->getSpecObj;
	my @commonpeaks     = @{ $spectrum->getcommonpeaks };
	my @commonpairs     = @{ $spectrum->getcommonpairs };
	my @commonionstheor = @{ $self->getCommonIons };
	my $PARAMS          = $self->getParams;
	my $maxionsize      = $self->getParams->{'maxionsize'};
	unless ($printout)
	{
		$printout = $self->getParams->{'printtables'};
	}
	my $userealint = $self->getParams->{'realintensities4xcorr'};
	my $id         = $self->getid;
	$id =~ s/::/_/g;
	my $commoncorrfile = join "", $basefilename, "_commonxcorr.txt";
	$printout
	  && ( open XCOMMON, ">$commoncorrfile" or die "xcorrelation:cannot open $commoncorrfile\n" );
	my ( $i, $j, $k, @commoniontable, @commonpeaktable, );
	my $tolerance = 0.2;

	if ( $PARAMS->{'xcorrprecision'} )
	{
		$tolerance = $PARAMS->{'xcorrprecision'};
	} else
	{
		$tolerance = $PARAMS->{'ms2tolerance'};
	}
	if ( $PARAMS->{'tolerancemeasure_ms2'} =~ /^ppm/i )
	{
		$tolerance = $tolerance * 1e-6 * 1000;    #ppm to amu measure
	}
	unless ($precision)
	{
		$precision = int( 1 / $tolerance );
	}
	unless ($delay)
	{
		$delay = int( $PARAMS->{'xcorrdelay'} );
	}
	my $weigthcommon = $PARAMS->{'commonxcorrweigth'};
	for $i ( 0 .. $maxionsize * $precision )
	{
		$commoniontable[$i] = $commonpeaktable[$i] = 0;
	}

	#todo: make reasonable matching function taking intensities into account
	if ($userealint)
	{
		foreach (@commonpairs)
		{
			$commonpeaktable[ int( $precision * $_->[0] ) ] = $_->[1];
		}
	} else
	{
		foreach (@commonpeaks)
		{
			$commonpeaktable[ int( $precision * $_ ) ] = $weigthcommon;
		}
	}
	foreach (@commonionstheor)
	{
		$commoniontable[ int( $precision * $_ ) ] = $weigthcommon;
	}
	my $rcommon                = calcxcorr( \@commonpeaktable, \@commoniontable, int($delay) );
	my $bbxcorr                = 0;
	my $xcorr_tolerance_window = $PARAMS->{'xcorr_tolerance_window'} || 2;
	#my $xcorr_tolerance_window = 5;
	my @x                      = map { $_->[0] / $precision } @$rcommon;
	my @y                      = map { $_->[1] } @$rcommon;
	my @peaks                  = sort { $b->[1] <=> $a->[1] } map {
	if ( abs( $_->[0] ) <= $xcorr_tolerance_window ) { [ $_->[0], $_->[1] ] }
		else                                             { }
	} findpeak( \@x, \@y );
	my $max_xcorr_x = $peaks[0]->[0];
	my $max_xcorr_y = $peaks[0]->[1];
	$printout && print "<br>Xcorr Common ions: MaxX: $max_xcorr_x MaxY: $max_xcorr_y<br>";	
	$printout && drawxcorrgraph( $rcommon, $precision, ( join "_", $basefilename, "backbonecorr.png" ), "xcorr backbone-ions", $max_xcorr_x,$max_xcorr_y );
	$printout && close(XCOMMON);
	$self->{'backbonexcorr'} = $max_xcorr_y;
}

sub xcorrelation_xlink
{
	my $self         = shift;
	my $basefilename = shift;
	my $printout     = shift;
	my $delay        = shift;
	my $precision    = shift;
	my $lossions     = shift;
	my $spectrum     = $self->getSpecObj;
	my $PARAMS       = $self->getParams;
	my $mincharge    = $self->minioncharge_xlinks;
	my $maxcharge    = $self->maxioncharge_xlinks;
	my $maxionsize   = $self->getParams->{'maxionsize'};
	unless ($printout)
	{
		$printout = $self->getParams->{'printtables'};
	}
	my $poolisotopes = $self->getParams->{'poolisotopes'};
	my $userealint   = $self->getParams->{'realintensities4xcorr'};
	my $id           = $self->getid;
	$id =~ s/::/_/g;
	my $xlinkcorrfile = join "", $basefilename, "_xlinkxcorr.txt";
	$printout
	  && ( open XLINK, ">$xlinkcorrfile" or die "xcorrelation:cannot open $xlinkcorrfile\n" );
	my ( $i, $j, $k, %rxlink, @xlinkiontable, @xlinkpeaktable, );
	my $tolerance = 0.2;

	if ( $PARAMS->{'xcorrprecision'} )
	{
		$tolerance = $PARAMS->{'xcorrprecision'};
	} else
	{
		$tolerance = $PARAMS->{'xlink_ms2tolerance'};
	}
	if ( $PARAMS->{'tolerancemeasure_ms2'} =~ /^ppm/i )
	{
		#$tolerance = $tolerance * 1e-6 * 1000;    #ppm to amu measure
	}
	unless ($precision)
	{
		$precision = int( 1 / $tolerance );
	}
	unless ($delay)
	{
		$delay = int( $PARAMS->{'xcorrdelay'} );
	}
	my $weigthxlinker     = $self->getParams->{'xlinkxcorrweigth'};
	my $weigthxlinkerloss = $self->getParams->{'commonlossxcorrweigth'};

	#todo: make reasonable matching function taking intensities into account
	my $charge;
	for $charge ( $mincharge .. $maxcharge )
	{
		for $i ( 0 .. $maxionsize * $precision )
		{
			$xlinkiontable[$i] = $xlinkpeaktable[$i] = 0;
		}
		if ($userealint)
		{
			my @xlinkpairs = ( @{ $spectrum->getxlinkpairs($charge) }, @{ $spectrum->getxlinkpairs("0") } );
			foreach (@xlinkpairs)
			{
				$xlinkpeaktable[ int( $precision * $_->[0] ) ] = $_->[1];
			}
		} else
		{
			my @xlinkpeaks = ( @{ $spectrum->getxlinkpeaks($charge) }, @{ $spectrum->getxlinkpeaks("0") } );
			foreach (@xlinkpeaks)
			{
				$xlinkpeaktable[ int( $precision * $_ ) ] = $weigthxlinker;
			}
		}
		my @xlinkionstheor = @{ $self->getXlinkIons->{$charge} };
		foreach (@xlinkionstheor)
		{
			$xlinkiontable[ int( $precision * $_ ) ] = $weigthxlinker;
		}
		if ($lossions)
		{
			my @xlinklossionstheor = @{ $self->getxlinkLossIons($charge) };
			foreach (@xlinklossionstheor)
			{
				$xlinkiontable[ int( $precision * $_ ) ] = $weigthxlinkerloss;
			}
		}
		$rxlink{$charge} = calcxcorr( \@xlinkpeaktable, \@xlinkiontable, int($delay) );
	}
	my ( $xlinkxcorr, @xcorrarray );
	$xlinkxcorr = 0;
	my $nchargestates = 0;
	foreach my $chargestate ( keys %rxlink )
	{
		$nchargestates++;
		if ( $delay == 0 )
		{
			$xlinkxcorr += $rxlink{$chargestate}->[0]->[1];
		} else
		{
			my $xcorrpairs = $rxlink{$chargestate};
			for $i ( 0 .. $#$xcorrpairs )
			{
				$xlinkxcorr += $xcorrpairs->[$i]->[1];
				$xcorrarray[$i]->[1] += $xcorrpairs->[$i]->[1];
				$xcorrarray[$i]->[0] = $xcorrpairs->[$i]->[0];

				#				$printout
				#				  && print XLINK $xcorrpairs->[$i]->[0] / $precision,
				#				  "\t",$xcorrpairs->[$i]->[1],"\n";
			}
		}
	}
	my $max_xcorr = 0;
	for $i ( 0 .. $#xcorrarray )
	{
		$printout
		  && print XLINK $xcorrarray[$i]->[0] / $precision,
		  "\t", $xcorrarray[$i]->[1], "\n";
	}
	
	$printout && drawxcorrgraph( \@xcorrarray, $precision, ( join "_", $basefilename, "xcorr.png" ), "xcorr cross-linker ions" );
	$printout && close(XLINK);
	$self->getSpecObj->autoxlink && ( $xlinkxcorr = $xlinkxcorr / $self->getSpecObj->autoxlink );
	$self->{'xlinkxcorr'} = $xlinkxcorr / $nchargestates;
	#$printout && print "<br>Xcorr Xlinker ions (static): ",$self->{'xlinkxcorr'} ,"<br>";
}

sub xcorrelation_xlink_normalized
{
	my $self         = shift;
	my $basefilename = shift;
	my $printout     = shift;
	my $delay        = shift;
	my $precision    = shift;
	my $lossions     = shift;
	my $spectrum     = $self->getSpecObj;
	my $PARAMS       = $self->getParams;
	my $mincharge    = $self->minioncharge_xlinks;
	my $maxcharge    = $self->maxioncharge_xlinks;
	my $maxionsize   = $self->getParams->{'maxionsize'};
	unless ($printout)
	{
		$printout = $self->getParams->{'printtables'};
	}
	my $poolisotopes = $self->getParams->{'poolisotopes'};
	my $userealint   = $self->getParams->{'realintensities4xcorr'};
	my $id           = $self->getid;
	$id =~ s/::/_/g;
	my $xlinkcorrfile = join "", $basefilename, "_xlinkxcorr.txt";
	$printout
	  && ( open XLINK, ">$xlinkcorrfile" or die "xcorrelation:cannot open $xlinkcorrfile\n" );
	my ( $i, $j, $k, %rxlink, @xlinkiontable, @xlinkpeaktable, );
	my $tolerance = 0.2;

	if ( $PARAMS->{'xcorrprecision'} )
	{
		$tolerance = $PARAMS->{'xcorrprecision'};
	} else
	{
		$tolerance = $PARAMS->{'xlink_ms2tolerance'};
	}
	if ( $PARAMS->{'tolerancemeasure_ms2'} =~ /^ppm/i )
	{
		#$tolerance = $tolerance * 1e-6 * 1000;    #ppm to amu measure
	}
	unless ($precision)
	{
		$precision = int( 1 / $tolerance );
	}
	unless ($delay)
	{
		$delay = int( $PARAMS->{'xcorrdelay'} );
	}
	my $weigthxlinker     = $self->getParams->{'xlinkxcorrweigth'}      || 1;
	my $weigthxlinkerloss = $self->getParams->{'commonlossxcorrweigth'} || 0.1;

	#todo: make reasonable matching function taking intensities into account
	my $charge;
  chargestate: for $charge ( $mincharge .. $maxcharge )
	{
		for $i ( 0 .. $maxionsize * $precision )
		{
			$xlinkiontable[$i] = $xlinkpeaktable[$i] = 0;
		}
		if ($userealint)
		{
			my @xlinkpairs = ( @{ $spectrum->getxlinkpairs($charge) }, @{ $spectrum->getxlinkpairs("0") } );

			#			if(scalar(@xlinkpairs)<10){
			#								next chargestate;
			#			}
			foreach (@xlinkpairs)
			{
				$xlinkpeaktable[ int( $precision * $_->[0] ) ] = $_->[1];
			}
		} else
		{
			my @xlinkpeaks = ( @{ $spectrum->getxlinkpeaks($charge) }, @{ $spectrum->getxlinkpeaks("0") } );
			my $nxlinkpeaks = scalar(@xlinkpeaks);

			#			print "npeaks for charge $charge: $nxlinkpeaks\n";
			#			if(scalar(@xlinkpeaks)<10){
			#
			#				next chargestate;
			#			}
			foreach (@xlinkpeaks)
			{
				$xlinkpeaktable[ int( $precision * $_ ) ] = $nxlinkpeaks;
			}
		}
		my @xlinkionstheor = @{ $self->getXlinkIons->{$charge} };
		foreach (@xlinkionstheor)
		{
			$xlinkiontable[ int( $precision * $_ ) ] = $weigthxlinker;
		}
		if ($lossions)
		{
			my @xlinklossionstheor = @{ $self->getxlinkLossIons($charge) };
			foreach (@xlinklossionstheor)
			{
				$xlinkiontable[ int( $precision * $_ ) ] = $weigthxlinkerloss;
			}
		}
		$rxlink{$charge} = calcxcorr( \@xlinkpeaktable, \@xlinkiontable, int($delay) );
	}
	my ( $xlinkxcorr, @xcorrarray );
	$xlinkxcorr = 0;
	my @chargestates  = keys %rxlink;
	my $nchargestates = scalar(@chargestates);
	foreach my $chargestate (@chargestates)
	{
		$nchargestates++;
		if ( $delay == 0 )
		{
			$xlinkxcorr += $rxlink{$chargestate}->[0]->[1];
		} else
		{
			my $xcorrpairs = $rxlink{$chargestate};
			for $i ( 0 .. $#$xcorrpairs )
			{
				$xlinkxcorr += $xcorrpairs->[$i]->[1];
				$xcorrarray[$i]->[0] = $xcorrpairs->[$i]->[0];
				$xcorrarray[$i]->[1] += $xcorrpairs->[$i]->[1] / $nchargestates;

				#				$printout
				#				  && print XLINK $xcorrpairs->[$i]->[0] / $precision,
				#				  "\t",$xcorrpairs->[$i]->[1],"\n";
			}
		}
	}
	for $i ( 0 .. $#xcorrarray )
	{
		$printout  && print XLINK $xcorrarray[$i]->[0] / $precision,"\t", $xcorrarray[$i]->[1], "\n";
	}
	##### new find max within x delay for $rcommon
	my $xcorr_tolerance_window = $PARAMS->{'xcorr_tolerance_window'};
	my @x                      = map { $_->[0] / $precision } @xcorrarray;
	my @y                      = map { $_->[1] } @xcorrarray;
	my @peaks                  = sort { $b->[1] <=> $a->[1] } map {
		if ( abs( $_->[0] ) <= $xcorr_tolerance_window ) { [ $_->[0], $_->[1] ] }
		else                                             { }
	} findpeak( \@x, \@y, );
	#print Dumper (@xcorrarray);
	my $max_xcorr_x = $peaks[0]->[0];
	my $max_xcorr_y = $peaks[0]->[1];
	$printout && print "Xcorr Xlinker ions: MaxX: $max_xcorr_x MaxY: $max_xcorr_y<br>";
#####
	$printout && drawxcorrgraph( \@xcorrarray, $precision, ( join "_", $basefilename, "xcorr.png" ), "xcorr cross-linker ions" , $max_xcorr_x,$max_xcorr_y);
	$printout && close(XLINK);
	$self->{'xlinkxcorr'} = $max_xcorr_y;
}

sub calcxcorr
{
	my $series1  = shift;
	my $series2  = shift;
	my $maxdelay = shift;
	my @x        = @{$series1};
	my @y        = @{$series2};
	my $n        = $#$series1;
	my ( $i, $j, $delay, $mx, $my, $sx, $sy, $sxy, $denom, @r );

	#   /* Calculate the mean of the two series x[], y[] */
	$my = 0;
	$my = 0;
	for ( $i = 0 ; $i <= $n ; $i++ )
	{		
		$mx += $x[$i];
		$my += $y[$i];
	}

	#print "n $n mean x $mx mean y =$my\n";
	#if ($n==0){
	#print "warning, n was zero in xcorr.";
	#$n=1;
	#}
	$mx /= ($n);
	$my /= ($n);

	#   /* Calculate the denominator */
	$sx = 0;
	$sy = 0;
	for ( $i = 0 ; $i <= $n ; $i++ )
	{
		$sx += ( $x[$i] - $mx ) * ( $x[$i] - $mx );
		$sy += ( $y[$i] - $my ) * ( $y[$i] - $my );
	}
	$denom = sqrt( $sx * $sy );

	#   /* Calculate the correlation series */
	for ( $delay = -$maxdelay ; $delay <= $maxdelay ; $delay++ )
	{
		$sxy = 0;
		for ( $i = 0 ; $i < $n ; $i++ )
		{
			$j = $i + $delay;
			unless ( $j < 0 || $j >= $n )
			{
				$sxy += ( $x[$i] - $mx ) * ( $y[$j] - $my );
			} else
			{
				next;
			}

			#         /* Or should it be (?)
			#			if ( $j < 0 || $j >= $n ) { next; }
			#			else {
			#				$sxy += ( $x[$i] - $mx ) * ( $y[$j] - $my );
			#
			#				#$sxy += ($x[$i] - $mx) * (-$my);
			#				#$sxy += ( $x[$i] - $mx ) * ( $y[$j] - $my );
			#			}
			#*/;
		}
		if ( $denom != 0 )
		{
			push @r, [ $delay, ( $sxy / $denom ) ]    #      /* r is the correlation coefficient at "delay" */
			                                          #print "$delay ",$sxy / $denom,"\n";
		} else
		{
			push @r, [ $delay, 0 ]                    #      /* r is the correlation coefficient at "delay" */
			                                          #print "$delay ",$sxy / $denom,"\n";
		}
	}
	return \@r;
}

sub findpeak
{
	my $x            = shift;
	my $y            = shift;
	my $cutoff       = shift || 0;
	my $window       = shift || 1;
	my $derivative_y = Statistics::Derivative1( $x, $y );
	my %peaks        = ();
	#print "TRUE";
	for my $i ( 1 .. $#$x )
	{

		if ( ( $derivative_y->[$i] * $derivative_y->[ $i - 1 ] ) < 0 )
		{
			if ( $y->[$i] > $y->[ $i - 1 ] )
			{
				$peaks{$i}->{'peak'} = $y->[$i];
			} else
			{
				$peaks{$i}->{'peak'} = $y->[ $i - 1 ];
			}
		}
	}
	foreach my $index ( sort { $a <=> $b } keys %peaks )
	{
		my $minintensity = $cutoff * $peaks{$index}->{'peak'};
	  right: for my $j ( $index + 1 .. $#$x )
		{
			if ( ( $derivative_y->[$j] > 0 ) || ( $y->[$j] < $minintensity ) )
			{
				$peaks{$index}->{'right'} = $x->[$j];
				last right;
			} else
			{
				$peaks{$index}->{'right'} = $x->[$#$x];
			}
		}
	  left: for ( my $j = $index - 1 ; $index > 0 ; $j-- )
		{
			if ( ( $derivative_y->[$j] < 0 ) || ( $y->[$j] < $minintensity ) )
			{
				$peaks{$index}->{'left'} = $x->[ $j + 1 ];
				last left;
			} else
			{
				$peaks{$index}->{'left'} = $x->[0];
			}
		}
	}
	my @peakpairs = map { [ $x->[$_], $peaks{$_}->{'peak'}, $peaks{$_}->{'left'}, $peaks{$_}->{'right'} ] }
	  sort { $a <=> $b } keys %peaks;
	return @peakpairs;
}

sub drawxcorrgraph
{
	my $dat       = shift;
	my $precision = shift;
	my $filelabel = shift;
	my $title     = shift;
	my $maxpeak_x = shift;
	my $maxpeak_y=shift;
	

	### define the series that markes the picked peak
	my @z;
	
	if ($maxpeak_y){
	foreach my $arrayref (@$dat){
	my $xvalue=$arrayref->[0]/$precision;
	if ($xvalue == $maxpeak_x){
	push @z,$maxpeak_y;	
	}else{
	push @z,0;
	}	
	}
	}
	#print Dumper (\@z);
	
	if ($maxpeak_x && $maxpeak_y){

	}
	
	
	if ( defined($dat) )
	{
		require GD::Graph::lines;
		my ( @x, @y );
		@x = map { $_->[0] / $precision } @$dat;
		@y = map { $_->[1] } @$dat;
		#@z = ();    #map{if($_->[0]==$maxpeak_x){1}else{0}} @$dat;

		#print "@data\n";
		my $graph = GD::Graph::lines->new( 400, 300 );
		my @col = qw(black red blue cyan);
		$graph->set(
			x_label      => 'delta mz',
			x_label_skip => 5,
			zero_axis    => 1,
			t_margin     => 20,
			b_margin     => 20,
			l_margin     => 20,
			r_margin     => 20,
			line_width   => 3,

			#dclrs        => [qw(black red blue cyan)],
			dclrs   => \@col,
			y_label => "cross-link cross-correlation",
			title   => "$title",
		) or die $graph->error;
		
		
		my $gd = $graph->plot( [ \@x, \@y, \@z ] ) or die $graph->error;
		
		open( IMG, ">$filelabel" ) or die $!;
		binmode IMG;
		print IMG $gd->png;
		close IMG;
	}
}

1;
