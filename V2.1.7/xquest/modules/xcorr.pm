package xcorr;
use strict;
#---------------------------------------------------------------------------
# Module: xcorr.pm
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

use File::Spec;
use File::Basename;

sub xcorrelation {
	my $dtafile1 = shift;
	my $dtafile2 = shift;
	my $range=shift;
	my $select=shift;
	my $precision=shift;
	my $dir=shift;
	
	
	open DTA1, "<$dtafile1" or die "cannot open $dtafile1 $!";
	open DTA2, "<$dtafile2" or die "cannot open $dtafile2 $!";
	print "calc xcorrelation for $dtafile1 and $dtafile2 \n<br>";
	
	my ( @dtafile1, @dtafile2 );
	my $i = 0;
	<DTA1>; #crop dta header
	while (<DTA1>) {
		chomp;
		my @tmp = split;

		push @dtafile1, \@tmp;

	}
	my $i = 0;
	<DTA2>; #crop dta header
	while (<DTA2>) {
		chomp;
		my @tmp = split;
		push @dtafile2, \@tmp;

	}
	$dtafile1 =~ s/.dta//;
	$dtafile2 =~ s/.dta//;
	my $graphname=File::Spec->catfile($dir,(join "", basename($dtafile1), "_", basename($dtafile2), "_xcorr.png"));
	my $xcorrfile = File::Spec->catfile($dir,(join "", basename($dtafile1), "_", basename($dtafile2), "_xcorr.txt"));
	open XCORR, ">$xcorrfile" or die "cannot open $xcorrfile $! ";
	my ( @dat1, @dat2 );

	my ($i);
	for $i ( 0 .. int(2000*$precision )) {
		$dat1[$i] = $dat2[$i] = 0;
	}

	
	foreach (@dtafile1) {
		$dat1[ int( $precision*$_->[0] ) ] = sqrt($_->[1]);
	}

	foreach (@dtafile2) {
		$dat2[ int($precision* $_->[0] ) ] = sqrt($_->[1]);
	}

	my $xcorr = _calcxcorr( \@dat1, \@dat2, $range );
	my $peaks=peakpick($xcorr,6);

	my ( $xcorr0, $xcorr80 );
	foreach $i (@$xcorr) {
		print XCORR "", $i->[0], "\t", $i->[1], " \n";
		if ( $i->[0]==0) {
			$xcorr0 = $i->[1];
		}
		if ( $i->[0]==$select->[0]*$precision) {
			$xcorr80 = $i->[1];
		}
	}
##	if($xcorr80>$select->[1]){
#	print "xcorr0: $xcorr0\t";
#	print "xcorr at ",$select->[0],": $xcorr80\n";
##	}
drawgraph( $xcorr, $precision, $graphname,$peaks );
close(XCORR),close(DTA1),close(DTA2);
}


sub drawgraph {
	my $dat       = shift;
	my $precision = shift;
	my $filelabel = shift;
	my $peaks=shift;

	if ( defined($dat) ) {

		my ( @x, @y );
		foreach (@$dat) {
			push @x, $_->[0] / $precision;
			push @y, $_->[1];
		}

		# print "@data\n";
		require GD::Graph::lines;
		my $graph = GD::Graph::lines->new( 800, 600 );
		$graph->set(
					 x_label      => 'delta mz',
					 x_label_skip => 20,
					 zero_axis    => 1,
					 t_margin     => 20,
					 b_margin     => 20,
					 l_margin     => 20,
					 r_margin     => 20,
					 line_width   => 3,
					 dclrs        => [qw(black red blue cyan)],
					 y_label      => 'xcorr',
					 title        => "$filelabel",
		  )
		  or die $graph->error;
		my $gd = $graph->plot( [ \@x, \@y ] ) or die $graph->error;
		open( IMG, ">$filelabel" ) or die $!;
		binmode IMG;
		print IMG $gd->png;
		close IMG;
	}
}



sub _calcxcorr {
	my $series1  = shift;
	my $series2  = shift;
	my $maxdelay = shift;
	my @x        = @{$series1};
	my @y        = @{$series2};
	my $n        = $#$series1;
	my ( $i, $j, $delay, $mx, $my, $sx, $sy, $sxy, $denom, @r );

	#   /* Calculate the mean of the two series x[], y[] */
	$mx = 0;
	$my = 0;
	for ( $i = 0 ; $i <= $n ; $i++ ) {
		$mx += $x[$i];
		$my += $y[$i];
	}

	#print "n $n mean x $mx mean y =$my\n";
	$mx /= ($n);
	$my /= ($n);

	#   /* Calculate the denominator */
	$sx = 0;
	$sy = 0;
	for ( $i = 0 ; $i <= $n ; $i++ ) {
		$sx += ( $x[$i] - $mx ) * ( $x[$i] - $mx );
		$sy += ( $y[$i] - $my ) * ( $y[$i] - $my );
	}

	$denom = sqrt( $sx * $sy );

	#   /* Calculate the correlation series */
	for ( $delay = -$maxdelay ; $delay < $maxdelay ; $delay++ ) {
		$sxy = 0;
		for ( $i = 0 ; $i < $n ; $i++ ) {
			$j = $i + $delay;

			unless ( $j < 0 || $j >= $n ) {
				$sxy += ( $x[$i] - $mx ) * ( $y[$j] - $my );
			}
			else {
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
			*/

		}
		if ( $denom != 0 ) {
			push @r,
			  [ $delay, ( $sxy / $denom )
			  ]    #      /* r is the correlation coefficient at "delay" */
			       #print "$delay ",$sxy / $denom,"\n";
		}
		else {
			push @r, [ $delay,
					0 ] #      /* r is the correlation coefficient at "delay" */
			            #print "$delay ",$sxy / $denom,"\n";
		}
	}

	return \@r;
}

sub dotproduct {
	my $series1  = shift;
	my $series2  = shift;
	my $maxdelay = shift;
	my @x        = @{$series1};
	my @y        = @{$series2};
	my $n        = $#$series1;
	my ( $i,$mx, $my, $sx, $sy, $sxy, $denom, @r );

	#   /* Calculate the mean of the two series x[], y[] */
	$my = 0;
	$my = 0;
	for ( $i = 0 ; $i <= $n ; $i++ ) {
		$mx += $x[$i];
		$my += $y[$i];
	}

	#print "n $n mean x $mx mean y =$my\n";
	$mx /= ($n);
	$my /= ($n);

	#   /* Calculate the denominator */
	$sx = 0;
	$sy = 0;
	for ( $i = 0 ; $i <= $n ; $i++ ) {
		$sx += ( $x[$i] - $mx ) * ( $x[$i] - $mx );
		$sy += ( $y[$i] - $my ) * ( $y[$i] - $my );
	}

	$denom = sqrt( $sx * $sy );

	#   /* Calculate the correlation series */
		$sxy = 0;
		for ( $i = 0 ; $i < $n ; $i++ ) {

				$sxy += ( $x[$i] - $mx ) * ( $y[$i] - $my );

		
		}
		my $dotproduct=0;
		if ( $denom != 0 ) {
			my $dotproduct=$sxy / $denom ;
		}
	return $dotproduct;
}



sub peakpick{
my $vector=shift;
my $n=shift;
my @sortedpairs = sort {$b->[1]<=>$a->[1]} @$vector;
my @nhighest = @sortedpairs[0..int($n)];
return \@nhighest;
}

sub nmostintense{
my $vector=shift;
my $n=shift;
my @sorted= sort {$b<=>$a} @$vector;
return @sorted[0..int($n)];
}

1;