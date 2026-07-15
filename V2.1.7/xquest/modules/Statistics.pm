package Statistics;
use strict;
#---------------------------------------------------------------------------
# Module: Statistics.pm
# Author(s): 
# Description: Module that provides statistics functions.
#---------------------------------------------------------------------------

sub Choose {    # Probability of getting $k heads is $n tosses
    my  $n=shift;
    my $k=shift;
    return undef unless defined $n && defined $k;
    my ( $result, $j ) = ( 1, 1 );
    if ( $k > $n || $k < 0 ) {
        return 0;
    }
    while ( $j <= $k ) {
        $result *= $n--;
        $result /= $j++;
    }
    return $result;
}
sub Derivative1 {
    my $x=shift;
    my $y=shift;
    my @y2;
    my $n=$#{$x};
    $y2[0]=($y->[1]-$y->[0])/($x->[1]-$x->[0]);
    $y2[$n]=($y->[$n]-$y->[$n-1])/($x->[$n]-$x->[$n-1]);
    my $i;
    for($i=1; $i<$n; $i++) {
    	my $denominator=($x->[$i+1]-$x->[$i-1]);
    	if($denominator){
	$y2[$i]=($y->[$i+1]-$y->[$i-1])/($x->[$i+1]-$x->[$i-1]);
    		}
    		else{
    			warn "Derivative1: denominator is 0\n";
    			$y2[$i]=0;
    		}
    }
    return \@y2;
}

sub Derivative2 {
    my ($x,$y,$yp1,$ypn)=@_;
    my $n=$#{$x};
    my ($i,@y2,@u);
    if (!defined $yp1) {
	$y2[0]=0; $u[0]=0;
    }
    else {
	$y2[0]=-0.5;
	$u[0]=(3/($x->[1]-$x->[0]))*(($y->[1]-$y->[0])/($x->[1]-$x->[0])-$yp1);
    }
    for($i=1; $i<$n; $i++) {
	my $sig=($x->[$i]-$x->[$i-1])/($x->[$i+1]-$x->[$i-1]);
	my $p=$sig*$y2[$i-1]+2.0;
	$y2[$i]=($sig-1.0)/$p;
	$u[$i]=(6.0*( ($y->[$i+1]-$y->[$i])/($x->[$i+1]-$x->[$i])-
		      ($y->[$i]-$y->[$i-1])/($x->[$i]-$x->[$i-1])
		     )/
		($x->[$i+1]-$x->[$i-1])-$sig*$u[$i-1])/$p;
    }
    my ($qn,$un);
    if (!defined $ypn) {
	$qn=0;
	$un=0;
    }
    else {
	$qn=0.5;
	$un=(3.0/($x->[$n]-$x->[$n-1]))*
	    ($ypn-($y->[$n]-$y->[$n-1])/($x->[$n]-$x->[$n-1]));
    }
    $y2[$n]=($un-$qn*$u[$n-1])/($qn*$y2[$n-1]+1.0);
    for($i=$n-1; $i>=0; --$i) {
	$y2[$i]=$y2[$i]*$y2[$i+1]+$u[$i];
    }
    return \@y2;
}

sub squaresum {
	my $vector = shift;
	my $mean   = shift || mean($vector);
	#print "vector: @$vector mean $mean\n";
	
	my $squaresum = 0;
	foreach my $value (@$vector) {
		$squaresum += ( $value - $mean )**2;
	}
	#print "squaresum: $squaresum";
	return $squaresum;
	
}

sub Binomial {    # probability of $k successes in $n attempts, given probability of $p
    my $n=shift;
    my $k=shift;
    my $p =shift;
    return $k == 0 if $p == 0;
    return $k != $n if $p == 1;
    return Choose( $n, $k ) * $p ** $k * ( 1 - $p ) ** ( $n - $k );
}

sub cumulativeBinomial {    # probability of $k successes in $n attempts, given probability of $p
    my $n=shift;
    my $k=shift;
    my $p =shift;
    return $k == 0 if $p == 0;
    return $k != $n if $p == 1;
    return 1  if $k > $n;
    my $j=0;
    my $pcum=0;
    for $j (0 .. $k){
      $pcum+=Choose( $n, $j ) * $p ** $j * ( 1 - $p ) ** ( $n - $j );
    }
    
return $pcum;
}


sub getpercentile{	
my $vector=shift;
my $value=shift;
my @sorted=sort {$a<=>$b} @$vector;
my $i;
for($i=0;$i<=$#sorted;$i++){
	if ($value<$sorted[$i]){
	#	print "$i: $value\n";
		return $i/scalar(@sorted); 
	}
 }
}


sub sum {
	my $array = shift;
	my $sum   = 0;
	foreach (@$array) {
		$sum += $_;
	}
	return $sum;
}

sub mean {
	my $array = shift;
	my $sum   = 0;
	my $n     = 0;
	foreach (@$array) {
		$sum += $_;
		$n++;
	}
	if ($n) {
		return $sum / $n;
	}
	else {
		return 0;
	}
}

sub average{
	my $array=shift;
return mean($array);
}

sub mean_array {
	my $array = shift;

	my $n         = 0;
	
	my %profilehash=();
	my ($i,$j);
	for $i (0..$#$array) {
	my @meanarray = ();	
		foreach $j (0..$#{$array->[$i]}) {
		push @{$profilehash{$j}},$array->[$i]->[$j];	
	  }
	}
	
	my @profilearray=();
#	@profilearray=map{mean($profilehash{$_})} sort keys %profilehash;

	foreach my $runID (sort keys %profilehash){
 	push @profilearray,mean($profilehash{$runID});
	}
	
	
return \@profilearray;
}




sub median_array {
	my $array = shift;

	my $n         = 0;
	
	my %profilehash=();
	my ($i,$j);
	for $i (0..$#$array) {
	my @meanarray = ();	
		foreach $j (0..$#{$array->[$i]}) {
		push @{$profilehash{$j}},$array->[$i]->[$j];	
	  }
	}
	
	my @profilearray=();
	foreach my $runID (sort keys %profilehash){
	push @profilearray,median($profilehash{$runID});
	}
return \@profilearray;
}

sub sem_array {
	my $array = shift;

	my $n         = 0;
	
	my %profilehash=();
	my ($i,$j);
	for $i (0..$#$array) {
	my @meanarray = ();	
		foreach $j (0..$#{$array->[$i]}) {
		push @{$profilehash{$j}},$array->[$i]->[$j];	
	  }
	}
	
	my @profilearray=();
	foreach my $runID (sort keys %profilehash){
	push @profilearray,sem($profilehash{$runID});
	}
return \@profilearray;
}

sub stdev_array {
	my $array = shift;

	my $n         = 0;
	
	my %profilehash=();
	my ($i,$j);
	for $i (0..$#$array) {
	my @meanarray = ();	
		foreach $j (0..$#{$array->[$i]}) {
		push @{$profilehash{$j}},$array->[$i]->[$j];	
	  }
	}
	
	my @profilearray=();
	foreach my $runID (sort keys %profilehash){
	push @profilearray,stdev($profilehash{$runID});
	}
return \@profilearray;
}



sub max{
my $vector=shift;
my @sorted=sort {$b<=>$a} @$vector;
 return $sorted[0];
}





sub n {
	my $array = shift;
	return scalar(@$array);
}

sub var {
	my $vector = shift;
	my $mean   = mean($vector);
	my $n      = n($vector);

	my $squaresum = 0;
	foreach (@$vector) {
		$squaresum += ( $_ - $mean )**2;
	}
	if ( $n > 0 ) {
		return ( $squaresum / $n );
	}
	else {
		return 0;
	}
}

sub median {
	my $vector = shift;
	my @sorted = sort { $a <=> $b } @$vector;
	return $sorted[ int( scalar(@sorted) / 2 ) ];
}



sub topten{
	my $vector = shift;
	my @sorted = sort { $b <=> $a } @$vector;
	my $toptenindex=scalar(@$vector)/10;
	my @toptenslice= @sorted[0..int($toptenindex)];
	
	return median(\@toptenslice);
}


sub topn{
	my $vector = shift;
	my $percentile = shift;
	my $fraction=100/$percentile;
	my @sorted = sort { $b <=> $a } @$vector;
	my $toptenindex=scalar(@$vector)/$fraction;
	my @toptenslice= @sorted[0..int($toptenindex)];
	
	return median(\@toptenslice);
}
sub stdev {
	my $vector = shift;
	my $n      = n($vector);

	if ( $n > 0 ) {
		return ( sqrt( var($vector) ) );
	}
	else {
		return 0;
	}
}

sub sem {
	my $vector = shift;
	my $n      = n($vector);
	if ( $n > 0 ) {
		return ( stdev($vector) / sqrt($n) );
	}
	else {
		return 0;
	}
}

1;