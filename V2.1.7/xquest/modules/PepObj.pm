package PepObj;
use strict;
#---------------------------------------------------------------------------
# Module: PepObj.pm
# Author(s): Thomas Walzthoeni, xQuest specific modules are based on 
# original work by Oliver Rinner.
# Description: Module for peptides.
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
use Ionseries;
use Data::Dumper;

sub new
{
	my $class = shift();
	my $self  = {};
	bless $self, $class;
	my $sequence  = shift;
	my $id        = shift;
	my $desc      = shift;
	my $MSTAB     = shift;
	my $PARAMS    = shift;
	my $verbose   = shift;
	my $proteinID = shift;
	my $seenion   = shift;
	$self->{'PARAMS'} = $PARAMS;
	$self->{'MSTAB'}  = $MSTAB;
	$self->{'Hatom'}  = $PARAMS->{'Hatom'};
	$verbose = $self->{'verbose'} = $PARAMS->{'verbose'};
	$self->{'seq'}       = $sequence;
	$self->{'desc'}      = $desc;
	$self->{'length'}    = length($sequence);
	$self->{'peptideID'} = $id;
	$self->{'seenion'}   = $seenion;

	#print Dumper ($PARAMS);
	if ($proteinID)
	{
		$self->{'proteinID'} = $proteinID;
	} else
	{
		$id =~ /(^\S+)::\d+::\d+/g;
		$self->{'proteinID'} = $1;
	}
	unless ( $PARAMS->{'redundant_peps'} )
	{
		$self->{'nonred'} = 1;
	}
	$self->{'description'} = $desc;
	my ( $parentmass, $Mplus1, $Mplus2, $Mplus3, $iontable, $fragmenttypes, $losstypes ) = Ionseries::calcions( $sequence, $MSTAB, $PARAMS, $verbose );
	$self->{'fragmenttypes'} = $fragmenttypes;
	$self->{'losstypes'}     = $losstypes;
	$self->{'parentmass'}    = $parentmass;
	$self->{'iontable'}      = $iontable;

	if ( $PARAMS->{'averageMS2'} )
	{
		$self->{'averagemass'} = calcaveragemass( $sequence, $MSTAB );
	}

	$self->{'chargestates'} = [ $Mplus1, $Mplus2, $Mplus3 ];
	$verbose && $self->printtable;
	return $self;
}

sub getfragmenttypes
{
	my $self = shift;
	return $self->{'fragmenttypes'};
}

sub getlosstypes
{
	my $self = shift;
	return $self->{'losstypes'};
}

sub calcaveragemass
{
	my $sequence    = shift;
	my $MSTAB       = shift;
	my @seq         = split //, $sequence;
	my $averagemass = 0;
	foreach my $aminoacid (@seq)
	{
		$averagemass += $MSTAB->{$aminoacid}->{'average'};
	}
	return $averagemass;
}

sub getaveragemass
{
	my $self = shift;
	return $self->{'averagemass'};
}

sub Mplus1
{
	my $self = shift;
	return $self->{'chargestates'}->[0];
}

sub Mplus2
{
	my $self = shift;
	return $self->{'chargestates'}->[1];
}

sub Mplus3
{
	my $self = shift;
	return $self->{'chargestates'}->[2];
}

sub molweight
{
	my $self = shift;
	return $self->{'parentmass'};
}

sub getmonoisotopicmass
{
	my $self = shift;
	return $self->{'parentmass'};
}

sub verbose
{
	my $self = shift;
	return $self->{'verbose'};
}

sub getlength
{
	my $self = shift;
	return $self->{'length'};
}

sub printtable
{
	my $self = shift;
	my $i;
	my $fwdiontypes = $self->fwdiontypes;
	my $reviontypes = $self->reviontypes;
	print $self->id, " ", $self->desc, "\t", $self->seq, "\n";
	print "parent mass = ", $self->molweight, "\n";
	print "[M+H]1+ = ",     $self->Mplus1,    "\n";
	print "[M+2H]2+ = ",    $self->Mplus2,    "\n";
	print "[M+3H]3+ = ",    $self->Mplus3,    "\n";
	my @seq           = split //, $self->seq;
	my $ionseries     = $self->getiontable;
	my @fragmenttypes = sort keys %$ionseries;
	my @losstypes     = sort keys %{ $ionseries->{ $fragmenttypes[0] } };
	my @index =
	  sort { $a <=> $b }
	  keys %{ $ionseries->{ $fragmenttypes[0] }->{'standard'} };
	print "fwd\t";

	foreach my $fragmenttype (@fragmenttypes)
	{
		if ( $fragmenttype =~ /$fwdiontypes/ )
		{
			foreach my $losstype (@losstypes)
			{
				print "$fragmenttype-$losstype ";
				print "\t";
			}
		}
	}
	print "AA\trev\t";
	foreach my $fragmenttype (@fragmenttypes)
	{
		if ( $fragmenttype =~ /$reviontypes/ )
		{
			foreach my $losstype (@losstypes)
			{
				print "$fragmenttype-$losstype ";
				print "\t";
			}
		}
	}
	print "\n";
	for $i ( 0 .. $#index )
	{

		#	print $seq[$i],"\t",$index[$i],"\t";
		print $index[$i], "\t";
		foreach my $fragmenttype (@fragmenttypes)
		{
			if ( $fragmenttype =~ /$fwdiontypes/ )
			{
				foreach my $losstype (@losstypes)
				{
					if ( defined( $ionseries->{$fragmenttype}->{$losstype}->{ $index[$i] } ) )
					{
						print $ionseries->{$fragmenttype}->{$losstype}->{ $index[$i] }, "\t";
					} else
					{
						print "-\t";
					}
				}
			}
		}
		print $seq[$i], "\t", $index[ $#index - $i ], "\t";
		foreach my $fragmenttype (@fragmenttypes)
		{
			if ( $fragmenttype =~ /$reviontypes/ )
			{
				foreach my $losstype (@losstypes)
				{
					if ( defined( $ionseries->{$fragmenttype}->{$losstype}->{ $index[ $#index - $i ] } ) )
					{
						print $ionseries->{$fragmenttype}->{$losstype}->{ $index[ $#index - $i ] }, "\t";
					} else
					{
						print "-\t";
					}
				}
			}
		}
		print "\n";
	}
	print "\n";
}

sub printpeps
{
	my $self    = shift;
	my $verbose = shift;
	my ( $i, $j );
	print $self->id, " ", $self->desc, "\t", $self->seq, "\n";
	if ($verbose)
	{
		print "parent mass = ", $self->molweight, "\n";
		print "[M+H]1+ = ",     $self->Mplus1,    "\n";
		print "[M+2H]2+ = ",    $self->Mplus2,    "\n";
		print "[M+3H]3+ = ",    $self->Mplus3,    "\n";
		for $i ( 1 .. $self->length - 1 )
		{    #ignore first and last b_ion
			print "b", $i + 1, "\t",    $self->get_b_ions->[$i],      "\t";
			print "b", $i + 1, "_2+\t", $self->get_b_ions2plus->[$i], "\n";
		}
		print "\n";
		for $i ( 0 .. $self->length )
		{
			print "y", $i + 1, "\t",    $self->get_y_ions->[$i],      "\t";
			print "y", $i + 1, "_2+\t", $self->get_y_ions2plus->[$i], "\n";
		}
		print "\n";

		#remove last b_ion
		foreach ( @{ $self->get_b_ions_waterloss } )
		{
			print "b-H20\t ", $_, "\n";
		}
		foreach ( @{ $self->get_b_ions_waterloss2plus } )
		{
			print "b2+-H20\t ", $_, "\n";
		}
	}
}

sub seq
{
	my $self = shift;
	return $self->{'seq'};
}

sub length
{
	my $self = shift;
	return length( $self->{'seq'} );
}

sub id
{
	my $self = shift;
	if ( $self->{'nonred'} )
	{
		return $self->{'seq'};
	} else
	{
		return $self->{'peptideID'};
	}
}

sub protid
{
	my $self = shift;
	return $self->{'proteinID'};
}

sub protidstring
{
	my $self = shift;
	return join ",", @{ $self->{'proteinID'} };
}

sub desc
{
	my $self = shift;
	return $self->{'description'};
}

sub DESTROY
{
	my $this = shift;
}

sub getfwdrev_lastxlinkpossibilities
{
	my $self            = shift;
	my $sequence        = $self->seq;
	my $PARAMS          = $self->getParams;
	my @seq             = split //, $sequence;
	my @revseq          = reverse(@seq);
	my $reversesequence = reverse($sequence);
	my $xlinkedAA       = $PARAMS->{'AArequired'};
	my $nocutatxlink    = $PARAMS->{'nocutatxlink'};
	my ( @Kpositionfwd, $lastKfwd, @Kpositionrev, $lastKrev );

	while ( $sequence =~ /$xlinkedAA/gi )
	{
		push @Kpositionfwd, pos($sequence);    ### Push all xlinksites into an array
	}

	#print Dumper(@Kpositionfwd);
	#print "KPOS-1".$Kpositionfwd[-1]."\n";
	#	print "LENGTH:".$self->getlength."\n";
	#if ( ( $Kpositionfwd[-1] == length($sequence) ) && $nocutatxlink )
	if ( ( $Kpositionfwd[-1] == ( $self->getlength ) ) && $nocutatxlink )    ### $Kpositionfwd[-1] ist the last array element
	{
		if ( defined( $Kpositionfwd[-2] ) )
		{
			$lastKfwd = $Kpositionfwd[-2];
		} else
		{
			$lastKfwd = 1;
		}
	} else
	{
		$lastKfwd = $Kpositionfwd[-1];
	}
	while ( $reversesequence =~ /$xlinkedAA/gi )
	{
		push @Kpositionrev, pos($reversesequence);
	}
	$lastKrev = $Kpositionrev[-1];    #cut at 5' K would be allowed, check that
	return ( $lastKfwd, $lastKrev );
}

sub getfwdrev_firstxlinkpossibilities
{
	my $self            = shift;
	my $sequence        = $self->seq;
	my $PARAMS          = $self->getParams;
	my @seq             = split //, $sequence;
	my @revseq          = reverse(@seq);
	my $reversesequence = reverse($sequence);
	my ( @Kpositionfwd, $firstKfwd, @Kpositionrev, $firstKrev );
	while ( $sequence =~ /($PARAMS->{'AArequired'})/gi )
	{
		push @Kpositionfwd, pos($sequence);
	}
	$firstKfwd = $Kpositionfwd[0];
	while ( $reversesequence =~ /($PARAMS->{'AArequired'})/gi )
	{
		push @Kpositionrev, pos($reversesequence);
	}
	$firstKrev = $Kpositionrev[0];
	return ( $firstKfwd, $firstKrev );
}

sub iontag_getionsforindex
{
	my $self            = shift;
	my $sequence        = $self->seq;
	my $PARAMS          = $self->getParams;
	my @seq             = split //, $sequence;
	my @revseq          = reverse(@seq);
	my $reversesequence = reverse($sequence);
	my $verbose         = $self->verbose;
	my $chargestates    = $PARAMS->{'indexcharges_common'};
	my ( $lastKfwd, $lastKrev ) = $self->getfwdrev_lastxlinkpossibilities;
	$verbose && print "indexing candidate common ions for $sequence from 1 to $lastKfwd fwd and 1 to $lastKrev reverse, chargestates @$chargestates\n";
	####### CHANGED 29/03/2010
	my @fwdions = @{ $self->getfwdions( 1, $lastKfwd - 1, $chargestates, undef ) };
	my @reverseions = @{ $self->getrevions( 1, $lastKrev - 1, $chargestates, undef ) };
	my @ions = ( @fwdions, @reverseions );
	return \@ions;
}

sub iontag_getpossiblecommonions
{
	my $self          = shift;
	my $sequence      = $self->seq;
	my $PARAMS        = $self->getParams;
	my $verbose       = $self->verbose;
	my $chargestates  = $PARAMS->{'ioncharge_common'};
	my $matchlossions = $PARAMS->{'uselossionsformatching'};
	my ( $lastKfwd, $lastKrev ) = $self->getfwdrev_lastxlinkpossibilities;
	$verbose
	  && print "get candidate common ions for $sequence from 1 to $lastKfwd fwd and 1 to $lastKrev reverse\n";
	my @fwdions = @{ $self->getfwdions( 1, $lastKfwd - 1, $chargestates, $matchlossions ) };
	my @reverseions = @{ $self->getrevions( 1, $lastKrev - 1, $chargestates, $matchlossions ) };
	my @ions = ( @fwdions, @reverseions );
	$verbose && print "@ions\n";
	return \@ions;
}

sub iontag_getpossiblexlinkions
{
	my $self            = shift;
	my $deltaMr         = shift;
	my $minioncharge    = shift;
	my $maxioncharge    = shift;
	my $verbose         = $self->verbose;
	my $Hatom           = $self->getHatom;
	my $sequence        = $self->seq;
	my $PARAMS          = $self->getParams;
	my @seq             = split //, $sequence;
	my @revseq          = reverse(@seq);
	my $seqlength       = scalar(@seq);
	my $reversesequence = reverse($sequence);
	my ( $firstKfwd, $firstKrev, @Kpositionfwd, @Kpositionrev );
	my $seqindexlength = $#seq;
	my $matchlossions  = $PARAMS->{'uselossionsformatching'};
	( $firstKfwd, $firstKrev ) = $self->getfwdrev_firstxlinkpossibilities;
	#$verbose=1;
	$verbose && print "get candidate xlink ions for $sequence from $firstKfwd to $seqlength fwd and $firstKrev to $seqlength reverse deltaMr = $deltaMr\n";
	########to do : change getfwdions to accomodate start and end of ionseries
	 
	
	### CHANGED BY TW
	### GET FWD IONS FROM START TO LAST XLINKABLE SITE-1
	#my @fwdions = @{ $self->getfwdions( $firstKfwd, undef, [1], $matchlossions ) };
	#print "Sequence: $sequence\n";
	#print "First K FROM FWD: $firstKfwd\n";
	#print "First K FROM REV: $firstKrev\n";
	#print "IONTABLE:\n";
	#$self->printtable;
	
	my @fwdions = @{ $self->getfwdions( $firstKfwd, $seqlength, [1], $matchlossions ) };
	
	#print Dumper (\@fwdions);
	my @reverseions = @{ $self->getrevions($firstKrev, $seqlength, [1], $matchlossions ) };
	#	print Dumper (\@reverseions);
	#	exit;
	my @ions       = ( @fwdions, @reverseions );
	my $minionsize = $PARAMS->{'minionsize'};
	my $maxionsize = $PARAMS->{'maxionsize'};
	my %xlinkionhash;
	for my $charge ( $minioncharge .. $maxioncharge )
	{
		my @tmp = ();
		$verbose && print "charge $charge: ";
		foreach my $ion (@ions)
		{
			my $xlinkion = ( $ion + $deltaMr + $Hatom * ( $charge - 1 ) ) / $charge;
			if ( $xlinkion >= $minionsize && $xlinkion <= $maxionsize )
			{
				push @tmp, $xlinkion;
			}
		}
		$verbose && print "@tmp\n";
		$xlinkionhash{$charge} = \@tmp;
	}
	$verbose && print "\n";
	return \%xlinkionhash;
}
### creates a hash where the fwd ions are indexed by the number of the aa in the sequence
### the ions are stored as an arrayreference.
sub get_fwdions_indexed
{
	my $self         = shift;
	my $chargestates = shift;
	my $Hatom        = $self->getHatom;
	my $PARAMS       = $self->getParams;
	my $minionsize   = $PARAMS->{'minionsize'};
	my $maxionsize   = $PARAMS->{'maxionsize'};
	my $fwdiontypes  = $self->fwdiontypes;
	unless ($chargestates)
	{
		$chargestates = [1];
	}
	my $ionseries = $self->getiontable;
	my $ions      = {};
	### Associate the fwd ions with AA positions
	foreach my $iontype ( keys %$ionseries )
	{

		#print $iontype."\n";
		unless ( $iontype =~ /$fwdiontypes/ )
		{
			next;
		}
		foreach my $losstype ( keys %{ $ionseries->{$iontype} } )
		{

			#print $losstype. "\n";
			## This is the subhash with AA positions as index and masses as values
			my $subhash = $ionseries->{$iontype}->{$losstype};
			foreach my $aa ( keys %$subhash )
			{
				### put the ions into a hash as an array where the index is the aa
				### Create the array if not exists
				my $arrayref;
				if ( $ions->{$aa} )
				{
					$arrayref = $ions->{$aa};
				} else
				{
					## create new arrayref
					$arrayref = [];
					$ions->{$aa} = $arrayref;
				}
				## add the ions for all acharge states to the arrayref
				my $ion = $ionseries->{$iontype}->{$losstype}->{$aa};
				foreach my $charge (@$chargestates)
				{

					#print "Charge $charge\n";
					my $chargedions = ( ( $ion + ( $charge - 1 ) * $Hatom ) / $charge );
					if ( ( $chargedions >= $minionsize ) && ( $chargedions <= $maxionsize ) )
					{
						push( @$arrayref, $chargedions );
					}
				}
			}
		}
	}

	#print Dumper ($ions);
	return $ions;
}

sub get_revions_indexed
{
	my $self         = shift;
	my $chargestates = shift;
	my $Hatom        = $self->getHatom;
	my $PARAMS       = $self->getParams;
	my $minionsize   = $PARAMS->{'minionsize'};
	my $maxionsize   = $PARAMS->{'maxionsize'};
	my $reviontypes  = $self->reviontypes;
	unless ($chargestates)
	{
		$chargestates = [1];
	}
	my $ionseries = $self->getiontable;
	my $ions      = {};
	### Associate the fwd ions with AA positions
	foreach my $iontype ( keys %$ionseries )
	{

		#print $iontype."\n";
		unless ( $iontype =~ /$reviontypes/ )
		{
			next;
		}
		foreach my $losstype ( keys %{ $ionseries->{$iontype} } )
		{

			#print $losstype. "\n";
			## This is the subhash with AA positions as index and masses as values
			my $subhash = $ionseries->{$iontype}->{$losstype};
			foreach my $aa ( keys %$subhash )
			{
				### put the ions into a hash as an array where the index is the aa
				### Create the array if not exists
				my $arrayref;
				if ( $ions->{$aa} )
				{
					$arrayref = $ions->{$aa};
				} else
				{
					## create new arrayref
					$arrayref = [];
					$ions->{$aa} = $arrayref;
				}
				## add the ions for all acharge states to the arrayref
				my $ion = $ionseries->{$iontype}->{$losstype}->{$aa};
				foreach my $charge (@$chargestates)
				{

					#print "Charge $charge\n";
					my $chargedions = ( ( $ion + ( $charge - 1 ) * $Hatom ) / $charge );
					if ( ( $chargedions >= $minionsize ) && ( $chargedions <= $maxionsize ) )
					{
						push( @$arrayref, $chargedions );
					}
				}
			}
		}
	}

	#print Dumper ($ions);
	return $ions;
}

sub get_iontag_ions_indexed
{
	my $self          = shift;
	my $ionsfromindex = shift;              ## ref on an array with aa position as index and an arrayref with the ions
	my $PARAMS        = shift;
	my $startindex    = shift;
	my $stopindex     = shift;
	my $ms1mass       = $self->molweight;
	unless ($ionsfromindex)
	{
		print "Error PepObj.pm sub get_iontag_ions_indexed(): no ion or hashreference was passed\n";
		exit;
	}
	unless ($startindex) { die "No startindex defined $!\n" }
	## make an index with integers of the ions as keys and the positions as values
	my $lookuphash = {};
	#print "Startindex: $startindex ";
	#print "Stopindex: $stopindex\n";
	foreach my $key ( keys %$ionsfromindex )
	{
		## check if within start and stop
		if ( $key > $stopindex || $key < $startindex ) { next }
		### Extract the array
		my $arrayref = $ionsfromindex->{$key};
		## add this key
		$lookuphash->{$key} = $arrayref;
	}
	return $lookuphash;
}

sub generate_lookup_hash
{
	my $self          = shift;
	my $ionsfromindex = shift;              ## ref on an array with aa position as index and an arrayref with the ions
	my $PARAMS        = shift;
	my $ms1mass       = $self->molweight;
	unless ($ionsfromindex)
	{
		print "Error PepObj.pm sub match_ion_from_index(): no ion or hashreference was passed\n";
		exit;
	}
	## Intprecision defines the precision of the array
	my $intprecision = $PARAMS->{'ionindexintprecision'};
	unless ($intprecision) { $intprecision = 10; }
	## Picktolerance is the tolerance the ions were picked.
	my $picktolerance = $PARAMS->{'picktolerance'};
	if ( $PARAMS->{'tolerancemeasure_ms2'} =~ /^ppm/i )
	{
		#$picktolerance = $picktolerance * 1e-6 * $ms1mass;    #ppm to amu measure
	}
	
	#print "Picktolerance: $picktolerance\n";
	## make an index with integers of the ions as keys and the positions as values
	my $lookuphash   = {};
	my $inttolerance = int( $intprecision * $picktolerance ) + 1;
	foreach my $key ( keys %$ionsfromindex )
	{
		## check if within start and stop
		### Extract the array
		my $arrayref = $ionsfromindex->{$key};
		foreach my $ion (@$arrayref)
		{

			#print $ion. "";
			### Set the bins including the tolerance with the aa
			#my $inttolerance = int( $intprecision * $picktolerance ) + 1;
			### Adjust the iontolerance if ppmeasure is selected
			if ( $PARAMS->{'picktolerance_measure'} =~ /^ppm/i )
			{
			$picktolerance=$PARAMS->{'picktolerance'} * 1e-6 * $ion;
			$inttolerance = int( $intprecision * $picktolerance ) + 1;	
			#print "PPM measure intolerance adjusted to: $inttolerance\n";
			}
			
			my $intionmz = int( $ion * $intprecision );    #get intprecision resolution


			#print "Inttolerance is:$inttolerance Setting bin: ";
			for my $rangex ( $intionmz - $inttolerance .. $intionmz + $inttolerance )
			{
				## Set the bin in the hash with the aa number
				$lookuphash->{$rangex} = $key;
			}
		}
	}
	return $lookuphash;
}

sub get_aapos_that_matched_index
{
	my $self         = shift;
	my $iontomatch   = shift;  ## is now an arrayreference
	#my $PARAMS       = shift;
	my $verbose      = shift;
	$verbose=0;
	my $PARAMS  = $self->getParams;
	my $ms1mass      = $self->molweight;
	my $chargestates = $PARAMS->{'indexcharges_common'};
	### fwd and rev ions are generated from the iontable
	### impossible ions are already sorted out
	my $fwdions = $self->get_fwdions_indexed($chargestates);
	my $revions = $self->get_revions_indexed($chargestates);
	my ( $lastKfwd, $lastKrev ) = $self->getfwdrev_lastxlinkpossibilities;

	## Then generate the ion tag ions (discount 1 position)
	my $fwdiontagions = $self->get_iontag_ions_indexed( $fwdions, $PARAMS, 1, $lastKfwd - 1 );
	my $reviontagions = $self->get_iontag_ions_indexed( $revions, $PARAMS, 1, $lastKrev - 1 );
	## Then generate the lookup hash
	my $fwdlookuphash = $self->generate_lookup_hash( $fwdiontagions, $PARAMS );
	my $revlookuphash = $self->generate_lookup_hash( $reviontagions, $PARAMS );

	unless ($iontomatch)
	{
		print "Error PepObj.pm get_aapos_that_matched_index(): no ion was passed\n";
		exit;
	}
	
	#$verbose = 1;
	if ($verbose)
	{
		#print "IONTABLE:\n";
		#print $self->printtable;
		print "FWD ion tag ions until residue $lastKfwd: \n";
		print Dumper($fwdiontagions);
		print "REV ion tag ions until residue $lastKrev: \n";
		print Dumper($reviontagions);
	}
	
	#### Intprecision defines the precision of the array
	my $intprecision = $PARAMS->{'ionindexintprecision'};
	unless ($intprecision) { $intprecision = 10; }
	#### Picktolerance is the tolerance the ions were picked.
	
	#my $picktolerance = $PARAMS->{'picktolerance'};
	#if ( $PARAMS->{'tolerancemeasure_ms2'} =~ /^ppm/i )
	#{
	#	$picktolerance = $picktolerance * 1e-6 * $ms1mass;    #ppm to amu measure
	#}
	
	$verbose && print "Seenions were: " , join("," ,@{$self->{'seenion'}}),"\n";
	my $result = {};
	
	foreach my $ionseen(@{$self->{'seenion'}}){
	#### converted int value of the ion
	my $intionmz = int( $ionseen * $intprecision );        #get intprecision resolution
	
	
	#### Matching of the ions
	if ( $fwdlookuphash->{$intionmz} )
	{
		my $pos = $fwdlookuphash->{$intionmz};
		my $realmass=$fwdiontagions->{$pos}->[0];
		#print $realmass;
		#exit;
		my $delta= $realmass- $ionseen;
		push @{$result->{'fwd'}},$pos;
		#$result->{'fwd'}           = $pos;
		#$result->{'ionmatchedbin'} = $intionmz;
		push @{$result->{'ionmatchedbin'}},$intionmz;
		push @{$result->{'deltaerror'}},$delta;
		$verbose && print "Fwd ion num $pos matches\n";
	}
	if ( $revlookuphash->{$intionmz} )
	{
		my $pos = $revlookuphash->{$intionmz};
		my $realmass=$reviontagions->{$pos}->[0];
		my $delta= $realmass- $ionseen;
		push @{$result->{'rev'}},$pos;
		#$result->{'rev'}           = $pos;
		#$result->{'ionmatchedbin'} = $intionmz;
		push @{$result->{'ionmatchedbin'}},$intionmz;
		push @{$result->{'deltaerror'}},$delta;
		$verbose && print "Rev ion num $pos matches\n";
	}
	
	
	}
	return $result;
}

sub print_lookuphash
{
	my $self    = shift;
	my $hashref = shift;
	print "IONS:POS ";
	foreach my $key ( sort keys %$hashref )
	{
		print $key. ":" . $hashref->{$key} . "\t";
	}
	print "\n";
}

sub fisher_yates_shuffle {
my $string = shift; 
				
if (CORE::length($string)==1){
return $string;
}

my @deck=split(//,$string);
return unless @deck; # must not be empty!		
my $i = @deck;		
while (--$i) {			
my $j = int rand ($i+1);			
@deck[$i,$j] = @deck[$j,$i];			
}
my $suffeledstring=join("",@deck);
return $suffeledstring;
}


sub generate_decoy_sequence
{
	my $self      = shift;
	my $foundions = shift;    ## hashref with 'fwd' and 'rev' as key, matched positions are arrayrefs with the positions
	my $topo 	  = shift;
	my $PARAMS  = $self->getParams;
	my $report    = {};
	my $sumflippedaa=0;
	my $frac;
	my $nionseen=0;
	my $seq = $self->seq;
	my $peplen = $self->length;
	
	$report->{'fwd__seq'}=$seq;
	$report->{'rev__seq'}=reverse $seq;
	#### Determine the position of the aa that was seen in the index (fwd direction)
	my @foundfwdpositions;
	my @foundrevpositions;
	
	if ( $foundions->{'fwd'} )
	{	
		foreach my $ipos(@{$foundions->{'fwd'}}){
		push @foundfwdpositions,$ipos;
		push @{$report->{'ionpos'}},"fwd".$ipos;
		$nionseen++;
		}
		#$ionseen = $foundions->{'fwd'};
		#$report->{'ionsseen'}++;	
		#$report->{'ionpos'} .= "fwd" . join(",",@{$foundions->{'fwd'}});
	}
	
	if ( $foundions->{'rev'} )
	{
		foreach my $ipos(@{$foundions->{'rev'}}){
		my $ionpos = $peplen - $ipos + 1;   ## if a revion was seen one pos may be added as the fixed position
		push @foundrevpositions,$ionpos;
		push @{$report->{'ionpos'}},"rev".$ipos;
		$nionseen++;
		}
		#$report->{'ionpos'} .= "rev" . join(",",@{$foundions->{'rev'}});
		#$ionseen = $peplen - $foundions->{'rev'}+1; ## if a revion was seen one pos may be added as the fixed position
		#$report->{'ionsseen'}++;
		#$report->{'ionpos'} .= "rev" . join(",",@{$foundions->{'rev'}});
	}
	
	unless ($nionseen)
	{
		print "Error PepObj: No ion found to generate decoy sequence.\n";
		exit;
	}
	$report->{'num_ionsseen'}=$nionseen;
	#print "Total Number of ions that were seen: $nionseen\n";
	
	
	## calc the real error of the ion that matched
	$report->{'errors'}=$foundions->{'deltaerror'};
	
	#### Get All Crosslinkable sites
	my $xlinkedAA    = $PARAMS->{'AArequired'};
	my $nocutatxlink = $PARAMS->{'nocutatxlink'};
	my ( $lastKfwd, $lastKrev ) = $self->getfwdrev_lastxlinkpossibilities;
	
	my $Kpositionfwd = {};

	#####################################
	# Define the fixed positions
	# The positions correspond to the AA positions
	# that should be preserved: eg. seq AQS, fixpos 2: AQ
	# Positions are always preserved c terminal to the selected residue
	#####################################
	
	## Keeping all Seenions for this peptide
	if ($PARAMS->{'Decoyfixseenionpos'}){
	foreach my $pos(@foundfwdpositions){
	$Kpositionfwd->{$pos} = defined;	
	#print "Defined Kpos: $pos\n";	
	}
	foreach my $pos(@foundrevpositions){
	### Discount one position from revions, the n term AA should be preserved
	$Kpositionfwd->{$pos-1} = defined;	
	#print "Defined Kpos: $pos\n";	
	}
	
	
	}

	### FIX the C terminus otherwise new xlinksite may emerge!
	$Kpositionfwd->{$peplen-1}  = defined;
	$Kpositionfwd->{$peplen}  = defined;
	### Report the fixed positions
	$report->{'fixedpositions'}=$Kpositionfwd;
	
	## make an array from the sequence
	my @seqarray = split( //, $seq );
	## Subseq all substrings from the sequence
	my $start = 0;
	my @subseqs;
	## Sort the positions hash into an array
	my @keys = sort { $a <=> $b } keys %$Kpositionfwd;
	#print "FWD sequencs: $seq\n";
	#print "Fixpositions:",join (",",@keys),"\n";
	### go through all the fixedpositions array
	for ( my $i = 0 ; $i <= $#keys ; $i++ )
	{
		## get the length of the subseq
		my $len=$keys[ $i ]-$start;
		## the length to the next pos
		#print "START IS $start, length is $len, i is $i, fixpos is $keys[ $i ] ";
		## Subseq always takes a startpos and samples then the length without taking the startpos itself
		my $subseq = substr( $seq, $start, $len );
		#print "Subseq is: $subseq\n";
		push @subseqs, $subseq;
		
		## then get the fixed res itself
		## Set the start to the position 
		#$start=$keys[ $i ]-1;
		$start=$keys[ $i ];
		#print "Start is now: $start\n";
		#$subseq = substr( $seq, $start, 1 );
		#push @subseqs, $subseq;
		#print "START IS $start, length is 1";
		#print "Subseq is: $subseq\n";
		### then set the start to the next pos	
#$start++;
}
	#print Dumper(@subseqs);
	my $decoy;
	### generate the decoysequence
	foreach my $substring (@subseqs)
	{
		unless ($substring)
		{
			next;
		}
		## Count the number of elements in the substring
		my $length = CORE::length($substring);
		unless ($length == 1){
		$sumflippedaa=$sumflippedaa+$length;
		}
		
		if ($PARAMS->{'Randomdecoys'}){
		my $dcstring = fisher_yates_shuffle($substring);
		$decoy .= $dcstring;
		}else{
		$decoy .= reverse($substring);
		}
	}
	#print "Target: $seq\n";
	#print "Decoy: $decoy\n";
	$report->{'num_altered_aas'}=$sumflippedaa;
	$report->{'peplength'}=$peplen;
	$report->{'amount_alteredAA'}=$sumflippedaa/$peplen;
	$report->{'decoyseq'} = $decoy;
	
	return $report;

}

sub make_iontag_decoy
{
	my $self    = shift;
	my $topo    = shift;  ## the AA that i cross-linked
	my $verbose = shift;
	my $PARAMS  = $self->getParams;
	my $MSTAB   = $self->{'MSTAB'};
	my $seenion = $self->{'seenion'};  ## is an arrayref
	my $ionpos  = $self->get_aapos_that_matched_index( $seenion, $verbose );
	my $report  = $self->generate_decoy_sequence( $ionpos,$topo );
	return $report;
}

sub getfwdions
{
	my $self         = shift;
	my $startindex   = shift;
	my $stopindex    = shift;
	my $chargestates = shift;
	my $getlossions  = shift;
	my $Hatom        = $self->getHatom;
	my $PARAMS       = $self->getParams;
	my $minionsize   = $PARAMS->{'minionsize'};
	my $maxionsize   = $PARAMS->{'maxionsize'};
	my $fwdiontypes  = $self->fwdiontypes;
	unless ($startindex)
	{
		$startindex = 0;
	}

	#	unless ($stopindex)
	#	{
	#		$stopindex = 10000;
	#	}
	#print "STARTINDEX: $startindex, STOPINDEX: $stopindex\n";
	unless ($chargestates)
	{
		$chargestates = [1];
	}
	my $ionseries = $self->getiontable;

	#print Dumper ($ionseries);
	#exit;
	my @ions = ();
	my $i;
	foreach my $iontype ( keys %$ionseries )
	{

		#if ( $iontype =~ /[abc]/ ) {
		if ( $iontype =~ /$fwdiontypes/ )
		{
			foreach my $losstype ( keys %{ $ionseries->{$iontype} } )
			{
				my @index = sort { $a <=> $b } keys %{ $ionseries->{$iontype}->{$losstype} };
				if ($getlossions)
				{
					for $i ( 0 .. $#index )
					{
						if ( ( $index[$i] >= $startindex ) && ( $index[$i] <= $stopindex ) )
						{
							push @ions, $ionseries->{$iontype}->{$losstype}->{ $index[$i] };
						}
					}
				} elsif ( $losstype eq "standard" )
				{
					for $i ( 0 .. $#index )
					{
						if ( ( $index[$i] >= $startindex ) && ( $index[$i] <= $stopindex ) )
						{
							push @ions, $ionseries->{$iontype}->{$losstype}->{ $index[$i] };
						}
					}
				}
			}
		}
	}
	my @ioncharges = ();
	foreach my $charge (@$chargestates)
	{
		foreach my $ion (@ions)
		{
			my $chargedions = ( ( $ion + ( $charge - 1 ) * $Hatom ) / $charge );
			if (    ( $chargedions >= $minionsize )
				 && ( $chargedions <= $maxionsize ) )
			{
				push @ioncharges, ( ( $ion + ( $charge - 1 ) * $Hatom ) / $charge );
			}
		}
	}
	return \@ioncharges;
}

sub getrevions
{
	my $self         = shift;
	my $startindex   = shift;
	my $stopindex    = shift;
	my $chargestates = shift;
	my $getlossions  = shift;
	my $Hatom        = $self->getHatom;
	my $PARAMS       = $self->getParams;
	my $minionsize   = $PARAMS->{'minionsize'};
	my $maxionsize   = $PARAMS->{'maxionsize'};
	my $reviontypes  = $self->reviontypes;
	unless ($startindex)
	{
		$startindex = 0;
	}

	#	unless ($stopindex)
	#	{
	#		$stopindex = 10000;
	#	}
	unless ($chargestates)
	{
		$chargestates = [1];
	}
	my $ionseries = $self->getiontable();
	my @ions      = ();
	my $i;
	foreach my $iontype ( keys %$ionseries )
	{

		#		if ( $iontype =~ /[xyz]/ ) {
		if ( $iontype =~ /$reviontypes/ )
		{
			foreach my $losstype ( keys %{ $ionseries->{$iontype} } )
			{
				my @index =
				  sort { $a <=> $b }
				  keys %{ $ionseries->{$iontype}->{$losstype} };
				if ($getlossions)
				{
					for $i ( 0 .. $#index )
					{
						if (    ( $index[$i] >= $startindex )
							 && ( $index[$i] <= $stopindex ) )
						{
							push @ions, $ionseries->{$iontype}->{$losstype}->{ $index[$i] };
						}
					}
				} elsif ( $losstype eq "standard" )
				{
					for $i ( 0 .. $#index )
					{
						if (    ( $index[$i] >= $startindex )
							 && ( $index[$i] <= $stopindex ) )
						{
							push @ions, $ionseries->{$iontype}->{$losstype}->{ $index[$i] };
						}
					}
				}
			}
		}
	}
	my @ioncharges = ();
	foreach my $charge (@$chargestates)
	{
		foreach my $ion (@ions)
		{
			my $chargedions = ( ( $ion + ( $charge - 1 ) * $Hatom ) / $charge );
			if (    ( $chargedions >= $minionsize )
				 && ( $chargedions <= $maxionsize ) )
			{
				push @ioncharges, ( ( $ion + ( $charge - 1 ) * $Hatom ) / $charge );
			}
		}
	}
	return \@ioncharges;
}

sub getParams
{
	my $self = shift;
	return $self->{'PARAMS'};
}

sub getHatom
{
	my $self = shift;
	return $self->{'Hatom'};
}

sub getplus1ions
{
	my $self     = shift;
	my $iontable = $self->getiontable();
}

sub getiontable
{
	my $self = shift;
	return $self->{'iontable'};
}

sub fwdiontypes
{
	my $self = shift;
	return $self->getParams->{'fwd_ions'};
}

sub reviontypes
{
	my $self = shift;
	return $self->getParams->{'rev_ions'};
}
1;
