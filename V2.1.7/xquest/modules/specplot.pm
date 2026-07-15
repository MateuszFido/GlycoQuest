package specplot;
use strict;
#---------------------------------------------------------------------------
# Module: specplot.pm
# Author(s): Thomas Walzthoeni, xQuest specific modules are based on 
# original work by Oliver Rinner.
# Description: Module for spectrum plotting.
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
use GD;
sub new {
	my $class      = shift();
	my $self       = {};
	my $xdimension = shift;
	my $ydimension = shift;

	bless $self, $class;

	unless ( $xdimension && $ydimension ) {
		$xdimension = 800;
		$ydimension = 600;
	}
	$xdimension = $self->{'xsize'} = $xdimension;
	$ydimension = $self->{'ysize'} = $ydimension;

	# create a new image
	my $im = $self->{'image'} = new GD::Image( $xdimension, $ydimension );

	# allocate some colors
	$self->{'colorhash'}  = $self->_define_colorhash;
	$self->{'colorarray'} = [
		"green", "red",   "lightblue", "darkblue",
		"black", "black", "black",     "black",
		"black", "black", "black",     "black"
	];

	$self->{'plottypes'} = $self->_defineplotfunctions;

	# make the background transparent and interlaced

	$im->transparent('true');    #( $self->color('white') );
	$im->interlaced('true');

	return $self;
}

sub _getplottype {
	my $self = shift;
	return $self->{'plottypes'};
}

sub _getxsize {
	my $self = shift;
	return $self->{'xsize'};
}

sub _getysize {
	my $self = shift;
	return $self->{'ysize'};
}

sub _defineplotfunctions {
	my %functionhash = (
		"diamond"    => \&_drawdiamond,
		"cross"      => \&_drawcross,
		"annotation" => \&_drawannotation,
		"line" => \&_drawline
	);
	return \%functionhash;

}

sub setcolor {
	my $self       = shift;
	my @colorarray = @_;
	$self->{'colorarray'} = \@colorarray;
}

sub getcolorvector {
	my $self = shift;
	return $self->{'colorarray'};

}

sub _drawruler {
	my $self = shift;
	my $min  = shift;
	my $max  = shift;

	my $img   = $self->img;
	my $xsize = $self->_getxsize;
	my $ysize = $self->_getysize;

	$img->line(
		20,
		$ysize - 30,
		$xsize - 20,
		$ysize - 30,
		$self->color("black")
	);
	my $increment = ( $max - $min ) / 4;
	my ( $i, $j );

	for $i ( 0 .. 4 ) {
		for $j ( 0 .. 9 ) {
			$img->line(
				20 + $j * 19 + $i * 190, $ysize-30, 20 + $j * 19 + $i * 190, $ysize-26,
				$self->color("black")
			);
		}
		$img->line(
			20 + $i * 190,
			$ysize-30, 20 + $i * 190,
			$ysize-20, $self->color("black")
		);
		$img->string( gdSmallFont, 15 + $i * 190,
			$ysize-15, $min + $i * $increment,
			$self->color("black")
		);
	}
	$img->string( gdSmallFont, 300, $ysize-19, "mz", $self->color("black") );
}

sub drawlegend {
	my $self           = shift;
	my $xposition      = shift;
	my $yposition      = shift;
	my $labelcolorhash = shift;
	my $sortlist       = shift;

	my $gd = $self->img;
	my $i  = 0;
	if ($sortlist) {
		foreach my $label (@$sortlist) {
			$gd->string( gdSmallFont, $xposition, $yposition + $i * 20,
				$label, $self->color( $labelcolorhash->{$label} ) );
			$i++;
		}
	}
	else {
		foreach my $label ( keys %$labelcolorhash ) {
			$gd->string( gdSmallFont, $xposition, $yposition + $i * 20,
				$label, $self->color( $labelcolorhash->{$label} ) );
			$i++;
		}
	}
}

sub plotdata {
	my $self = shift;
	my $min  = shift;
	my $max  = shift;
	my $colorvector=shift||$self->getcolorvector;

	my @data        = @_;
	my $img         = $self->img;
	
	if ( !defined($min) && !defined($max) ) {
		( $min, $max ) = _findminmax( \@data );
	}
	if ( !defined($min) ) {
		$min = _findminx( \@data );
	}
	if ( !defined($max) ) {
		$max = _findmaxx( \@data );
	}

	my $minintensity = _findminintensity( \@data, $min, $max );
	my $maxintensity = _findmaxintensity( \@data, $min, $max );

	$self->_drawruler( $min, $max );

	#	print "minx: $min maxx: $max\n";
	#	print "maxi: $maxintensity\n";
	my ( $xscale, $yscale );
	my $offset = 0;
	if ( $maxintensity == 0 ) {
		$yscale = 0;
	}
	else {
		if ( $minintensity < 0 ) {
			$yscale = 500 / ( $maxintensity - $minintensity );
			$offset = 500 - abs($minintensity) * $yscale;

#print "max: $maxintensity min: $minintensity y-scale $yscale offset: $offset\n";
			$img->line(
				20, 570 - $offset,
				780, 570 - $offset,
				$self->color("black")
			);

		}
		else {
			$yscale = 500 / $maxintensity;
		}
	}
	if ( $max - $min == 0 ) {
		$xscale = 0;
	}
	else {
		$xscale = 760 / ( $max - $min );
	}

	$self->{'minx'}   = $min;
	$self->{'maxx'}   = $max;
	$self->{'minI'}   = $minintensity;
	$self->{'maxI'}   = $maxintensity;
	$self->{'yscale'} = $yscale;
	$self->{'xscale'} = $xscale;

	my $j = 0;
	foreach my $dataset (@data) {
		foreach my $intensitypair (@$dataset) {
			if (   $intensitypair->[0] >= $min
				&& $intensitypair->[0] <= $max )
			{
				my $x =
				  int( ( $xscale * ( $intensitypair->[0] - $min ) ) + 20 );
				my $y = int( 570 - $offset - $intensitypair->[1] * $yscale );
				$img->line( $x, 570 - $offset,
					$x, $y, $self->color( $colorvector->[$j] ) );

				#	print	$intensitypair->[0]," ",$intensitypair->[1]," $x $y\n";
			}
		}
		$j++;
	}
}

sub labelpeaks {
	my $self        = shift;
	my $plottypes   = shift;
	my $colorvector = shift;
	my $data        = shift;
	my $img         = $self->img;
	unless ($colorvector) {
		$colorvector = $self->getcolorvector;
	}

	my $drawsymbol   = $self->_getplottype;
	my $min          = $self->minx;
	my $max          = $self->maxx;
	my $maxintensity = $self->maxI;
	my $yscale       = $self->yscale;
	my $xscale       = $self->xscale;

	my $j = 0;
	foreach my $dataset (@$data) {
		foreach my $intensitypair (@$dataset) {
			if (   $intensitypair->[0] >= $min
				&& $intensitypair->[0] <= $max )
			{
				my $x =
				  int( ( $xscale * ( $intensitypair->[0] - $min ) ) + 20 );
				my $y        = int( 560 - $intensitypair->[1] * $yscale );
				my $label    = $intensitypair->[2];
				my $function = $drawsymbol->{ $plottypes->[$j] };
				if ($function) {

					#print "$function <br>";
					&$function( $self, $x, $y,
						$self->color( $colorvector->[$j] ), $label );
				}
			}
		}
		$j++;
	}
}


sub drawlineplots {
	my $self = shift;
	my $min  = shift;
	my $max  = shift;

	my @data        = @_;
	my $img         = $self->img;
	my $colorvector = $self->getcolorvector;

	if ( !defined($min) && !defined($max) ) {
		( $min, $max ) = _findminmax( \@data );
	}
	if ( !defined($min) ) {
		$min = _findminx( \@data );
	}
	if ( !defined($max) ) {
		$max = _findmaxx( \@data );
	}

	my $minintensity = _findminintensity( \@data, $min, $max );
	my $maxintensity = _findmaxintensity( \@data, $min, $max );

	$self->_drawruler( $min, $max );

	#	print "minx: $min maxx: $max\n";
	#	print "maxi: $maxintensity\n";
	my ( $xscale, $yscale );
	my $offset = 0;
	if ( $maxintensity == 0 ) {
		$yscale = 0;
	}
	else {
		if ( $minintensity < 0 ) {
			$yscale = 500 / ( $maxintensity - $minintensity );
			$offset = 500 - abs($minintensity) * $yscale;

#print "max: $maxintensity min: $minintensity y-scale $yscale offset: $offset\n";
#			$img->line(
#				20, 570 - $offset,
#				780, 570 - $offset,
#				$self->color("red")
#			);

		}
		else {
			$yscale = 500 / $maxintensity;
		}
	}
	if ( $max - $min == 0 ) {
		$xscale = 0;
	}
	else {
		$xscale = 760 / ( $max - $min );
	}

	$self->{'minx'}   = $min;
	$self->{'maxx'}   = $max;
	$self->{'minI'}   = $minintensity;
	$self->{'maxI'}   = $maxintensity;
	$self->{'yscale'} = $yscale;
	$self->{'xscale'} = $xscale;


	my $j = 0;
	foreach my $dataset (@data) {
		my $x_before=0;
		my $y_before=0;
		foreach my $intensitypair (@$dataset) {
				#print "<li>",$intensitypair->[0], " ",$intensitypair->[1],"</li>";
			if (   $intensitypair->[0] >= $min
				&& $intensitypair->[0] <= $max )
			{
				my $x =
				  int( ( $xscale * ( $intensitypair->[0] - $min ) ) + 20 );
				my $y        = int( 560 - $intensitypair->[1] * $yscale );
			

						_drawlineplot( $self, $x, $y,$x_before||$x,$y_before||$y,
						$self->color( $colorvector->[$j] ));

			$x_before=$x;
			$y_before=$y;
			}
		}
		$j++;
	}
}



sub _drawlineplot {
	my $self    = shift;
	my $x = shift;
	my $y = shift;
	my $x_before = shift;
	my $y_before = shift;
	
	my $color   = shift;
	
	my $img = $self->img;
	$img->line( $x_before, $y_before, $x, $y, $color );
}



sub _drawdiamond {
	my $self    = shift;
	my $xcenter = shift;
	my $ycenter = shift;
	my $color   = shift;
	my $label   = shift;

	my $img = $self->img;
	$img->line( $xcenter - 3, $ycenter, $xcenter, $ycenter - 5, $color );
	$img->line( $xcenter + 3, $ycenter, $xcenter, $ycenter - 5, $color );

	$img->line( $xcenter - 3, $ycenter, $xcenter, $ycenter + 5, $color );
	$img->line( $xcenter + 3, $ycenter, $xcenter, $ycenter + 5, $color );

}

sub _drawannotation {
	my $self  = shift;
	my $x     = shift;
	my $y     = shift;
	my $color = shift;

#           my $courier = GD::Font->load('./cour.ttf');# or die "Can't load font";
#	 my $courier = GD::Font->load('./cour.ttf') or die "Can't load font cour.ttf $!";
	my $label = shift;
	my $img   = $self->img;

	$img->stringUp( gdSmallFont , $x - 4, $y, $label, $color );

	#$img->stringUp( $courier, $x - 4, $y, $label, $color );
}

sub _drawcross {
	my $self    = shift;
	my $xcenter = shift;
	my $ycenter = shift;
	my $color   = shift;
	my $label   = shift;

	my $img = $self->img;
	$img->line( $xcenter - 4, $ycenter, $xcenter + 4, $ycenter, $color );
	$img->line( $xcenter, $ycenter - 4, $xcenter, $ycenter + 4, $color );
}

sub minx {
	my $self = shift;
	return $self->{'minx'};
}

sub maxx {
	my $self = shift;
	return $self->{'maxx'};
}

sub maxI {
	my $self = shift;
	return $self->{'maxI'};
}

sub minI {
	my $self = shift;
	return $self->{'maxI'};
}

sub yscale {
	my $self = shift;
	return $self->{'yscale'};
}

sub xscale {
	my $self = shift;
	return $self->{'xscale'};
}

sub _findmaxintensity {
	my $data = shift;
	my $minx = shift;
	my $maxx = shift;
	my $max  = 0;
	foreach my $dataset (@$data) {
		foreach my $intensitypair (@$dataset) {
			if (   $intensitypair->[1] > $max
				&& $intensitypair->[0] > $minx
				&& $intensitypair->[0] < $maxx )
			{
				$max = $intensitypair->[1];
			}
		}

	}
	return $max;
}

sub _findminintensity {
	my $data = shift;
	my $minx = shift;
	my $maxx = shift;
	my $min  = 1e99;
	foreach my $dataset (@$data) {
		foreach my $intensitypair (@$dataset) {
			if (   $intensitypair->[1] < $min
				&& $intensitypair->[0] > $minx
				&& $intensitypair->[0] < $maxx )
			{
				$min = $intensitypair->[1];
			}
		}

	}
	return $min;
}

sub _findminx {
	my $data = shift;

	my $min = undef;
	foreach my $dataset (@$data) {
		foreach my $intensitypair (@$dataset) {
			if ( !defined($min) || $intensitypair->[0] < $min ) {
				$min = $intensitypair->[0];
			}
		}

	}
	return $min;
}

sub _findminmax {
	my $data = shift;
	my $max  = undef;
	my $min  = undef;
	foreach my $dataset (@$data) {
		foreach my $intensitypair (@$dataset) {
			if ( !defined($max) || $intensitypair->[0] > $max ) {
				$max = $intensitypair->[0];
			}
			if ( !defined($min) || $intensitypair->[0] < $min ) {
				$min = $intensitypair->[0];
			}
		}

	}
	return $min, $max;
}

sub _findmaxx {
	my $data = shift;
	my $max  = undef;
	foreach my $dataset (@$data) {
		foreach my $intensitypair (@$dataset) {
			if ( !defined($max) || $intensitypair->[0] > $max ) {
				$max = $intensitypair->[0];
			}

		}

	}
	return $max;
}

sub color {
	my $self  = shift;
	my $color = shift;
	if ( $self->{'colorhash'}->{$color} ) {
		return $self->{'colorhash'}->{$color};
	}
	else {
		return $self->{'colorhash'}->{"black"};
	}

}

sub _define_colorhash {
	my $self  = shift;
	my $im    = $self->img;
	my $white = $im->colorAllocate( 255, 255, 255 );
	my $grey  = $im->colorAllocate( 123, 123, 123 );

	my $black     = $im->colorAllocate( 0,   0,   0 );
	my $red       = $im->colorAllocate( 255, 0,   0 );
	my $green     = $im->colorAllocate( 0,   255, 0 );
	
	my $yellow     = $im->colorAllocate( 223,   223, 0 );
	my $blue      = $im->colorAllocate( 0,   0,   255 );
	my $darkblue  = $im->colorAllocate( 0,   0,   255 );
	my $lightblue = $im->colorAllocate( 20,  29,  200 );

	my %colorhash = (
		'red'  => $red,
		'grey' => $grey,

		'white'     => $white,
		'green'     => $green,
		'blue'      => $blue,
		'black'     => $black,
		'yellow'     => $yellow,
		
		'lightblue' => $lightblue,
		'darkblue'  => $darkblue,
	);
	return \%colorhash;
}

sub printimage {
	my $self     = shift;
	my $filename = shift;
	open( IMG, ">$filename" ) or die "cannot write $filename $!";
	binmode IMG;
	print IMG $self->img->png;
	close IMG;
}

sub img {
	my $self = shift;
	return $self->{'image'};
}
1;
