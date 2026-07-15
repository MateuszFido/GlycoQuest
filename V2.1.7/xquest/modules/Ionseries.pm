package Ionseries;
use strict;
#---------------------------------------------------------------------------
# Module: Ionseries.pm
# Author(s): Thomas Walzthoeni, xQuest specific modules are based on 
# original work by Oliver Rinner.
# Description: Calculation of ionseries.
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
use Data::Dumper;
sub calcions
{
	my $sequence          = shift;
	my $MSTAB             = shift;
	my $PARAMS            = shift;
	my $ionseries         = $PARAMS->{'ionseries'};
	my $fragmentresiduals = $PARAMS->{'fragmentresiduals'};  ## Are the masses that have to be subtracted or added for the spcific fragments
	my @seq               = split //, $sequence;
	my @revseq            = reverse(@seq);
	my $reversesequence   = reverse($sequence);
	my ( $fragmenttype, @fragmenttypes );

	foreach $fragmenttype ( keys %$ionseries )
	{
		if ( $ionseries->{$fragmenttype} )
		{
			push @fragmenttypes, $fragmenttype;    ## fragmenttypes are y b ...
		}
	}

	#print "calc ions for @fragmenttypes\n";
	my %iontable = ();
	my ( $i, $j );
	my $monomass    = 0;
	my $ms2masstype = 'native';
	if ( $PARAMS->{'averageMS2'} )
	{
		$ms2masstype = 'average';
	}
	for $i ( 0 .. $#seq )
	{
		$monomass += $MSTAB->{ $seq[$i] }->{'native'};
	}

	my @losstypes = ('standard');
	my %seen      = ();
	my $fwdions   = $PARAMS->{'fwd_ions'};    ## defined in xquest def:  a|b|c
	my $revions   = $PARAMS->{'rev_ions'};
	
	foreach $fragmenttype (@fragmenttypes)
	{
		if ( $fragmenttype =~ /$fwdions/ )
		{
			my $mass             = 0;
			my $fragmentresidual = $fragmentresiduals->{$fragmenttype};

			#print "fragmentresidual $fragmenttype = $fragmentresidual\n";
			my $waterlosstoggle = 0;
			my $NH3losstoggle   = 0;
			my $waterlossmass   = 0;
			my $NH3lossmass     = 0;
			for $i ( 0 .. $#seq )
			{
			
			
				$mass += $MSTAB->{ $seq[$i] }->{$ms2masstype};
				$iontable{$fragmenttype}->{'standard'}->{ $i + 1 } = ( $mass + $fragmentresidual );
				if ( $PARAMS->{'waterloss'} )
				{
					if ( $MSTAB->{ $seq[$i] }->{'-water'} )
					{
						$waterlosstoggle = 1;
						$waterlossmass   = $MSTAB->{ $seq[$i] }->{'-water'};
						unless ( $seen{'-H2O'}++ )
						{
							push @losstypes, '-H20';
						}
					}
					if ($waterlosstoggle)
					{
						$iontable{$fragmenttype}->{'-H2O'}->{ $i + 1 } = ( $mass + $fragmentresidual ) + $waterlossmass;
					}
				}
				if ( $PARAMS->{'nh3loss'} )
				{
					if ( $MSTAB->{ $seq[$i] }->{'-NH3'} )
					{
						$NH3losstoggle = 1;
						$NH3lossmass   = $MSTAB->{ $seq[$i] }->{'-NH3'};
						unless ( $seen{'-NH3'}++ )
						{
							push @losstypes, '-NH3';
						}
					}
					if ($NH3losstoggle)
					{
						$iontable{$fragmenttype}->{'-NH3'}->{ $i + 1 } = ( $mass + $fragmentresidual ) + $NH3lossmass;
					}
				}
			}
		} elsif ( $fragmenttype =~ /$revions/ )
		{
			my $mass             = 0;
			my $fragmentresidual = $fragmentresiduals->{$fragmenttype};
			my $waterlosstoggle  = 0;
			my $NH3losstoggle    = 0;
			my $waterlossmass    = 0;
			my $NH3lossmass      = 0;

			#print "fragmentresidual $fragmenttype = $fragmentresidual\n";
			for $i ( 0 .. $#revseq )
			{			
				
				$mass += $MSTAB->{ $revseq[$i] }->{$ms2masstype};
				$iontable{$fragmenttype}->{'standard'}->{ $i + 1 } = ( $mass + $fragmentresidual );
				if ( $PARAMS->{'waterloss'} )
				{
					if ( $MSTAB->{ $revseq[$i] }->{'-water'} )
					{
						$waterlosstoggle = 1;
						$waterlossmass   = $MSTAB->{ $revseq[$i] }->{'-water'};
						unless ( $seen{'-H2O'}++ )
						{
							push @losstypes, '-H2O';
						}
					}
					if ($waterlosstoggle)
					{
						$iontable{$fragmenttype}->{'-H2O'}->{ $i + 1 } = ( $mass + $fragmentresidual ) + $waterlossmass;
					}
				}
				if ( $PARAMS->{'nh3loss'} )
				{
					if ( $MSTAB->{ $revseq[$i] }->{'-NH3'} )
					{
						$NH3losstoggle = 1;
						$NH3lossmass   = $MSTAB->{ $revseq[$i] }->{'-NH3'};
						unless ( $seen{'-NH3'}++ )
						{
							push @losstypes, '-NH3';
						}
					}
					if ($NH3losstoggle)
					{
						$iontable{$fragmenttype}->{'-NH3'}->{ $i + 1 } = ( $mass + $fragmentresidual ) + $NH3lossmass;
					}
				}
			}
		}
	}
	my $parentmass = $monomass + 2 * $MSTAB->{'Hatom'}->{'native'} + $MSTAB->{'Oatom'}->{'native'};
	my $Mplus1     = $parentmass + $MSTAB->{'Hatom'}->{'native'};
	my $Mplus2     = ( $parentmass + 2 * $MSTAB->{'Hatom'}->{'native'} ) / 2;
	my $Mplus3     = ( $parentmass + 3 * $MSTAB->{'Hatom'}->{'native'} ) / 3;
	return ( $parentmass, $Mplus1, $Mplus2, $Mplus3, \%iontable, \@fragmenttypes, \@losstypes );
}

sub printiontable
{
	my $iontable = shift;
	my $i;
	my @iontypes = sort keys %$iontable;
	foreach my $iontype (@iontypes)
	{
		my @lossions = keys %{ $iontable->{$iontype} };
		foreach my $lossiontype ( sort @lossions )
		{
			my @ions =
			  sort { $a <=> $b }
			  keys %{ $iontable->{$iontype}->{$lossiontype} };
			for $i ( 0 .. $#ions )
			{
				print "$iontype\t$lossiontype\t", $ions[$i], "\t", $iontable->{$iontype}->{$lossiontype}->{ $ions[$i] }, "\n";
			}
		}
	}
}
1;
