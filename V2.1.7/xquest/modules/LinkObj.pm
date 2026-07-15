package LinkObj;
use strict;
#---------------------------------------------------------------------------
# Module: LinkObj.pm
# Author(s): Thomas Walzthoeni, xQuest specific modules are based on 
# original work by Oliver Rinner.
# Description: Module for all peptide types.
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

#use lib '/cluster/apps/perl_modules/1209/lib64/perl5/';
use FindBin;
#use lib "$FindBin::Bin/../../perl5";
use GD;
use GD::Graph::linespoints;
use GD::Graph::bars;
use GD::Graph::mixed;
use GD::Graph::colour;
use File::Spec;
use File::Basename;
use Statistics;
use specplot;
use Xcorrelation;
use Data::Dumper;
use Statistics::Descriptive;

sub new
{
	my $class = shift();
	my $self  = {};
	bless $self, $class;
	$self->{'xlinktype'}  = shift;
	$self->{'hitpepObjs'} = shift;
	$self->{'xlinkermw'}  = shift;
	$self->{'topology'}   = shift;
	my $matchObj = shift;
	my $PARAMS   = shift;
	$self->{'MSTAB'} = shift;
	my $specobj   = shift;
	my $annotate  = shift;
	my $addedmass = shift;
	$self->{'addedmass'} = $addedmass;
	$self->{'matchobj'}  = $matchObj;
	$self->{'params'}    = $PARAMS;
	my $masslistonly = $PARAMS->{'massmatchonly'};

	if ($specobj)
	{
		$self->{'specObj'} = $specobj;
	} elsif ($matchObj)
	{
		$self->{'specObj'} = $matchObj->getspectrum;
	}
	$self->calc_minioncharge_common;
	$self->calc_maxioncharge_common;
	$self->calc_minioncharge_xlinks;
	$self->calc_maxioncharge_xlinks;
	
	#print "<br>Z states<br>Mincommon:",$self->calc_minioncharge_common, ", Maxcommon:", $self->calc_maxioncharge_common,"<br>";
	 
	$self->{'id'} = $self->makeid;
	## HERE THE MATCHING IS DONE
	unless ($masslistonly)
	{
		$self->makeiontable;
		$self->generate_full_table;

		#exit;
		#$self->printhtmliontable;
		#exit;
		$self->matchions;

		#exit;
		if ($matchObj)
		{
			## check if there are matched ions
			## for intralinks there are sometime no matched ions due to the topology
			if ( ( $self->get_number_of_XlinkIons + $self->get_number_of_CommonIons ) == 0 )
			{
				$self->{'prescore'}       = 0;
				$self->{'prescore_alpha'} = 0;
				$self->{'prescore_beta'}  = 0;
			} else
			{
				$self->{'prescore'}       = $self->calcprescore;
				$self->{'prescore_alpha'} = $self->calcprescore_alpha;
				$self->{'prescore_beta'}  = $self->calcprescore_beta;
			}
			$self->{'error'}     = $self->calcerrorppm;
			$self->{'error_rel'} = $self->calcerrorppmrelativ;
		}
	}
	return $self;
}

sub generate_full_table
{
	my $self      = shift;
	my $iontable  = $self->getiontable;
	my $fulltable = {};
	my @keys      = sort { $a cmp $b } keys %$iontable;    ## for strings
	foreach my $key (@keys)
	{

		#print "$key <br>";
		my $index;
		if ( $key =~ /alpha/ && $key =~ /_b_/ )
		{
			$index = "alpha_b";
		}
		if ( $key =~ /alpha/ && $key =~ /_y_/ )
		{
			$index = "alpha_y";
		}
		if ( $key =~ /beta/ && $key =~ /_b_/ )
		{
			$index = "beta_b";
		}
		if ( $key =~ /beta/ && $key =~ /_y_/ )
		{
			$index = "beta_y";
		}
		if ( $key =~ /alpha/ && $key =~ /_c_/ )
		{
			$index = "alpha_c";
		}
		if ( $key =~ /alpha/ && $key =~ /_z_/ )
		{
			$index = "alpha_z";
		}
		if ( $key =~ /beta/ && $key =~ /_c_/ )
		{
			$index = "beta_c";
		}
		if ( $key =~ /beta/ && $key =~ /_z_/ )
		{
			$index = "beta_z";
		}
		if ( $key =~ /alpha/ && $key =~ /_a_/ )
		{
			$index = "alpha_a";
		}
		if ( $key =~ /alpha/ && $key =~ /_x_/ )
		{
			$index = "alpha_x";
		}
		if ( $key =~ /beta/ && $key =~ /_a_/ )
		{
			$index = "beta_a";
		}
		if ( $key =~ /beta/ && $key =~ /_x_/ )
		{
			$index = "beta_x";
		}
		unless ($index)
		{
			die "No index set!";
		}

		#print "Key is: ", $key."<br>";
		#print "Ions are: ";
		#print Dumper ($iontable->{$key});
		#print "<br>";
### get the keys
		my @ionkeys = sort { $a <=> $b } keys %{ $iontable->{$key} };

		#print @ionkeys;
		foreach my $pos (@ionkeys)
		{
### get the ion
			my $ion = $iontable->{$key}->{$pos};

			#print "Ion is: $ion<br>";
			push @{ $fulltable->{$index}->{$pos} }, $ion;
		}
	}
	$self->{'full_iontable'} = $fulltable;

	#print Dumper ($fulltable);
}

sub get_full_iontable
{
	my $self = shift;
	return $self->{'full_iontable'};
}

sub set_target_dc_label
{
	my $self  = shift;
	my $label = shift;
	$self->{'target_decoy_label'} = $label;
}

sub get_target_decoy_label
{
	my $self = shift;
	return $self->{'target_decoy_label'};
}

sub getMatchingIonpositions
{
	my $self = shift;
	return $self->{'ionposition_matchhash'};
}

sub matchions
{
	my $self = shift;
	my %ion2peakmatchhash;
	my $errorhash = {};
	( $self->{'xlinkmatches'}, $self->{'xlinkmatchingpeaks'} ) = $self->matchxlinks( \%ion2peakmatchhash, $errorhash );
	( $self->{'commonmatches'}, $self->{'commonmatchingpeaks'} ) = $self->matchcommon( \%ion2peakmatchhash, $errorhash );    ## matched ions, and peaks

	#print Dumper ($errorhash);
	#exit;
	#	print "<br>Common Matches Ions:";
	#	print Dumper($self->{'commonmatches'});
	#print Dumper($errorhash);
	#print "<br>Common Matches Peaks:";
	#print Dumper($self->{'commonmatchingpeaks'});
	#print "<br>CommonMatches Alpha:<br>";
	#print Dumper ( $self->{'commonmatchingpeaks_alpha'} );
	#	print "<br>CommonMatches Beta:<br>";
	#	print Dumper ( $self->{'commonmatchingpeaks_beta'} );
	#	print "<br>XlinkMatches Alpha:<br>";
	#	print Dumper ( $self->{'xlinkmatchingpeaks_alpha'} );
	#	print "<br>XlinkMatches Beta:<br>";
	#	print Dumper ( $self->{'xlinkmatchingpeaks_beta'} );
	#	print "<br>Xlink matches Ions:";
	#print Dumper($self->{'xlinkmatches'});
	#print "<br>Xlink matches Peaks:";
	#print Dumper($self->{'xlinkmatchingpeaks'});
	#	print "<br>";
## Is now changed to a hash where the key is the ion (theoretical peak), but the value is an arrayref
## --> multiple peaks may match to one ion
	$self->{'ion2peakmatchhash'} = \%ion2peakmatchhash;
	## flat matchhash is a hash where all ions are in and the value is one
	$self->{'flat_ionmatchhash'}     = $self->make_flationmatchhash;
	$self->{'ionposition_matchhash'} = $self->calcMatchingIonpositions;

	#print "<br>ION to PEAK hash:<br>";
	#print Dumper ($self->{'ion2peakmatchhash'});
	### calculate the error of the matched peaks to ions in ppm
	$self->calc_error_matched_ions($errorhash);

	#print Dumper ($self->getiontable);
	#exit;
	#print Dumper ($self->get_ion_annotationtable);
	#print "<br>Xlinkions:<br>";
	#print Dumper ($self->getXlinkIons);
	#print Dumper ($self->{'flationtable'});
	#print "<br>";
	#$self->printiontable;
	## make peak to ion hash to lookup the the peak ion pairs
	my %tmp;
	foreach my $ion ( keys %ion2peakmatchhash )
	{
		my $arrayref = $ion2peakmatchhash{$ion};
		foreach my $peak (@$arrayref)
		{
			$tmp{$peak} = $ion;
		}
	}
	$self->{'peak2ionmatchhash'} = \%tmp;
}

sub calc_error_matched_ions
{
	my $self      = shift;
	my $errorhash = shift;
	my @errorarray;

	#print Dumper ($errorhash);
	foreach my $ion ( keys %$errorhash )
	{
		my $errors = $errorhash->{$ion};
		foreach my $error (@$errors)
		{
			push @errorarray, $error;
		}
	}
## calculate the mean
	my $stat = Statistics::Descriptive::Full->new();
	$stat->add_data(@errorarray);
	$self->{'mean_error_matchedpeaks'}       = $stat->mean();
	$self->{'mean_error_matchedpeaks_stdev'} = $stat->standard_deviation();
}

sub get_matcherror_mean
{
	my $self = shift;
	return $self->{'mean_error_matchedpeaks'};
}

sub get_matcherror_stdev
{
	my $self = shift;
	return $self->{'mean_error_matchedpeaks_stdev'};
}

sub makeiontable
{
	my $self      = shift;
	my $xlinktype = $self->getxlinktype;
	if ( $xlinktype eq "xlink" )
	{
		$self->_makeiontable_xlink;
	} elsif ( $xlinktype eq "monolink" )
	{
		$self->_makeiontable_monolink;
	} elsif ( $xlinktype eq "intralink" )
	{
		$self->_makeiontable_intralink;
	} else
	{
		print "xlinktype $xlinktype is not defined";
		die $!;
	}
}

sub calccommonions
{
	my $self      = shift;
	my $xlinktype = $self->getxlinktype;
	if ( $xlinktype eq "xlink" )
	{
		$self->_calccommonions_xlink;
	} elsif ( $xlinktype eq "monolink" )
	{
		$self->_calccommonions_monolink;
	} elsif ( $xlinktype eq "intralink" )
	{
		$self->_calccommonions_intralink;
	}
}

sub printxlinks
{
	my $self       = shift;
	my $filehandle = shift;
	my $xlinktype  = $self->getxlinktype;
	if ( $xlinktype eq "xlink" )
	{
		$self->_printxlinks_xlink($filehandle);
	} elsif ( $xlinktype eq "monolink" )
	{
		$self->_printxlinks_monolink($filehandle);
	} elsif ( $xlinktype eq "intralink" )
	{
		$self->_printxlinks_intralink($filehandle);
	} else
	{
		die "xlinktype $xlinktype is not defined $!";
	}
}

sub drawpepstructure
{
	my $self      = shift;
	my $dat1      = shift;
	my $dat2      = shift;
	my $filelabel = shift;
	my $logscale  = shift;
	my $lossions  = shift;
	my $gd        = shift;
	my $xlinktype = $self->getxlinktype;
	if ( $xlinktype eq "xlink" )
	{
		$self->_drawpepstructure_xlink( $dat1, $dat2, $filelabel, $logscale, $lossions, $gd, $xlinktype );
	} elsif ( $xlinktype eq "monolink" )
	{
		$self->_drawpepstructure_monolink( $dat1, $dat2, $filelabel, $logscale, $lossions, $gd, $xlinktype );
	} elsif ( $xlinktype eq "intralink" )
	{
		$self->_drawpepstructure_intralink( $dat1, $dat2, $filelabel, $logscale, $lossions, $gd, $xlinktype );
	} else
	{
		die "xlinktype $xlinktype is not defined $!";
	}
}

sub _makeiontable_xlink
{
	my $self             = shift;
	my $pepobjs          = $self->gethitObjs;
	my $pepObj1          = $pepobjs->[0];
	my $pepObj2          = $pepobjs->[1];
	my $verbose          = $self->verbose;
	my $xlinkpositions   = $self->gettopology;
	my $xlinkpos1        = $xlinkpositions->[0];
	my $xlinkpos2        = $xlinkpositions->[1];
	my $reversexlinkpos1 = $self->getreversetopology( $pepObj1, $xlinkpos1 );
	my $reversexlinkpos2 = $self->getreversetopology( $pepObj2, $xlinkpos2 );
	my $sequence1        = $pepObj1->seq;
	my $sequence2        = $pepObj2->seq;

	#	my $xlinkermass     = $self->getxlinkermass;
	my $xlinkermass         = $self->gettotalxlinkermass;
	my $minioncharge_common = $self->minioncharge_common;
	my $maxioncharge_common = $self->maxioncharge_common;
	my $minioncharge_xlinks = $self->minioncharge_xlinks;
	my $maxioncharge_xlinks = $self->maxioncharge_xlinks;
	my $Hatom               = $pepObj1->getHatom;
	my $PARAMS              = $self->getParams;
	my @commonions          = ();
	my @commonions_alpha    = ();
	my @commonions_beta     = ();
	my %xlinkions           = ();
	my %xlinkions_alpha     = ();
	my %xlinkions_beta      = ();
	my ( $alphax, $betax );

	if ( $PARAMS->{'averageMS2'} )
	{
		$alphax = $pepObj1->getaveragemass + $xlinkermass;
		$betax  = $pepObj2->getaveragemass + $xlinkermass;
	} else
	{
		$alphax = $pepObj1->getmonoisotopicmass + $xlinkermass;
		$betax  = $pepObj2->getmonoisotopicmass + $xlinkermass;
	}
	$verbose && print "make iontable for ", $self->getid, "\n";
	my $pepobjiontable1 = $pepObj1->getiontable;
	my $pepobjiontable2 = $pepObj2->getiontable;
	my $fragmenttypes1  = $pepObj1->getfragmenttypes;
	my $fragmenttypes2  = $pepObj2->getfragmenttypes;
	my $losstypes1      = $pepObj1->getlosstypes;       #sort ( keys %{ $pepobjiontable1->{ $fragmenttypes[0] } } );
	my $losstypes2      = $pepObj2->getlosstypes;       #sort ( keys %{ $pepobjiontable1->{ $fragmenttypes[0] } } );
	$verbose && print "losstypes1 @$losstypes1\n";
	$verbose && print "losstypes2 @$losstypes2\n";
	my @index1 =
	  sort { $a <=> $b }
	  keys %{ $pepobjiontable1->{ $fragmenttypes1->[0] }->{'standard'} };
	my @index2 =
	  sort { $a <=> $b }
	  keys %{ $pepobjiontable2->{ $fragmenttypes2->[0] }->{'standard'} };
	my %iontable;
	my %flationtable;
	my $i                       = 0;
	my $xlinkpos1_indexstart    = 0;
	my $xlinkpos2_indexstart    = 0;
	my $revxlinkpos1_indexstart = 0;
	my $revxlinkpos2_indexstart = 0;
	my $fwdiontypes             = $PARAMS->{'fwd_ions'};
	my $reviontypes             = $PARAMS->{'rev_ions'};
	my $minionsize              = $PARAMS->{'minionsize'};
	my $maxionsize              = $PARAMS->{'maxionsize'};

	for my $charge ( $minioncharge_xlinks .. $maxioncharge_xlinks )
	{
		$xlinkions{$charge}       = [];
		$xlinkions_alpha{$charge} = [];
		$xlinkions_beta{$charge}  = [];
	}
	my $ncommonionsalpha = 0;
	my $ncommonionsbeta  = 0;
	my $nxlinkionsalpha  = 0;
	my $nxlinkionsbeta   = 0;
########## common ions alpha ########################################
	for my $charge ( $minioncharge_common .. $maxioncharge_common )
	{
		foreach my $iontype (@$fragmenttypes1)
		{
			foreach my $lossiontype (@$losstypes1)
			{
				if ( $iontype =~ /$fwdiontypes/ )
				{
					for ( $i = 0 ; $index1[$i] < $xlinkpos1 ; $i++ )
					{
						if ( defined( $pepobjiontable1->{$iontype}->{$lossiontype}->{ $index1[$i] } ) )
						{
							my $ionmz = ( $pepobjiontable1->{$iontype}->{$lossiontype}->{ $index1[$i] } + ( $charge - 1 ) * $Hatom ) / $charge;
							my $iontypestring = join "_", "alpha", "common", $iontype, $lossiontype, "plus$charge";
							$iontable{$iontypestring}->{ $index1[$i] } = $ionmz;
							my $flationtypestring = join "_", "alpha", "common", $iontype, $index1[$i], $lossiontype, "plus$charge";
							$flationtable{$flationtypestring} = $ionmz;
							if ( $ionmz > $minionsize && $ionmz < $maxionsize )
							{
								push @commonions,       $ionmz;
								push @commonions_alpha, $ionmz;
								$ncommonionsalpha++;
							}
						}
					}
					$xlinkpos1_indexstart = $i;
				} elsif ( $iontype =~ /$reviontypes/ )
				{
					for ( $i = 0 ; $index1[$i] < $reversexlinkpos1 ; $i++ )
					{
						if ( defined( $pepobjiontable1->{$iontype}->{$lossiontype}->{ $index1[$i] } ) )
						{
							my $ionmz = ( $pepobjiontable1->{$iontype}->{$lossiontype}->{ $index1[$i] } + ( $charge - 1 ) * $Hatom ) / $charge;
							my $iontypestring = join "_", "alpha", "common", $iontype, $lossiontype, "plus$charge";
							$iontable{$iontypestring}->{ $index1[$i] } = $ionmz;
							my $flationtypestring = join "_", "alpha", "common", $iontype, $index1[$i], $lossiontype, "plus$charge";
							$flationtable{$flationtypestring} = $ionmz;
							if ( $ionmz > $minionsize && $ionmz < $maxionsize )
							{
								push @commonions,       $ionmz;
								push @commonions_alpha, $ionmz;
								$ncommonionsalpha++;
							}
						}
					}
					$revxlinkpos1_indexstart = $i;
				}
			}
		}
		foreach my $iontype (@$fragmenttypes2)
		{
			foreach my $lossiontype (@$losstypes2)
			{
				if ( $iontype =~ /$fwdiontypes/ )
				{
					for ( $i = 0 ; $index2[$i] < $xlinkpos2 ; $i++ )
					{
						if ( defined( $pepobjiontable2->{$iontype}->{$lossiontype}->{ $index2[$i] } ) )
						{
							my $ionmz = ( $pepobjiontable2->{$iontype}->{$lossiontype}->{ $index2[$i] } + ( $charge - 1 ) * $Hatom ) / $charge;
							my $iontypestring = join "_", "beta", "common", $iontype, $lossiontype, "plus$charge";
							$iontable{$iontypestring}->{ $index2[$i] } = $ionmz;
							my $flationtypestring = join "_", "beta", "common", $iontype, $index2[$i], $lossiontype, "plus$charge";
							$flationtable{$flationtypestring} = $ionmz;
							if ( $ionmz > $minionsize && $ionmz < $maxionsize )
							{
								push @commonions,      $ionmz;
								push @commonions_beta, $ionmz;
								$ncommonionsbeta++;
							}
						}
					}
					$xlinkpos2_indexstart = $i;
				} elsif ( $iontype =~ /$reviontypes/ )
				{
					for ( $i = 0 ; $index2[$i] < $reversexlinkpos2 ; $i++ )
					{
						if ( defined( $pepobjiontable2->{$iontype}->{$lossiontype}->{ $index2[$i] } ) )
						{
							my $ionmz = ( $pepobjiontable2->{$iontype}->{$lossiontype}->{ $index2[$i] } + ( $charge - 1 ) * $Hatom ) / $charge;
							my $iontypestring = join "_", "beta", "common", $iontype, $lossiontype, "plus$charge";
							$iontable{$iontypestring}->{ $index1[$i] } = $ionmz;
							my $flationtypestring = join "_", "beta", "common", $iontype, $index1[$i], $lossiontype, "plus$charge";
							$flationtable{$flationtypestring} = $ionmz;
							if ( $ionmz > $minionsize && $ionmz < $maxionsize )
							{
								push @commonions,      $ionmz;
								push @commonions_beta, $ionmz;
								$ncommonionsbeta++;
							}
						}
					}
					$revxlinkpos2_indexstart = $i;
				}
			}
		}
	}
	################### xlink ions
	for my $charge ( $minioncharge_xlinks .. $maxioncharge_xlinks )
	{
		foreach my $iontype (@$fragmenttypes1)
		{
			foreach my $lossiontype (@$losstypes1)
			{

				#print
				#"charge $charge iontype $iontype losstype $lossiontype \n";
				#print "index1 @index1 xlinkpos1 = $xlinkpos1\n";
				if ( $iontype =~ /$fwdiontypes/ )
				{
					for ( $i = $xlinkpos1_indexstart ; $i <= $#index1 ; $i++ )
					{

						#print "index1 @index1 xlinkpos1 = $xlinkpos1 $i index:", $index1[$i],"\n";
						if ( defined( $pepobjiontable1->{$iontype}->{$lossiontype}->{ $index1[$i] } ) )
						{
							my $ionmz = ( $pepobjiontable1->{$iontype}->{$lossiontype}->{ $index1[$i] } + $betax + ( $charge - 1 ) * $Hatom ) / $charge;
							my $iontypestring = join "_", "alpha", "xlink", $iontype, $lossiontype, "plus$charge";
							$iontable{$iontypestring}->{ $index1[$i] } = $ionmz;
							my $flationtypestring = join "_", "alpha", "xlink", $iontype, $index1[$i], $lossiontype, "plus$charge";
							$flationtable{$flationtypestring} = $ionmz;
							if ( $ionmz > $minionsize && $ionmz < $maxionsize )
							{
								push @{ $xlinkions_alpha{$charge} }, $ionmz;
								push @{ $xlinkions{$charge} },       $ionmz;
								$nxlinkionsalpha++;
							}
						}
					}
				} elsif ( $iontype =~ /$reviontypes/ )
				{
					for ( $i = $revxlinkpos1_indexstart ; $i <= $#index1 ; $i++ )
					{
						if ( defined( $pepobjiontable1->{$iontype}->{$lossiontype}->{ $index1[$i] } ) )
						{
							my $ionmz = ( $pepobjiontable1->{$iontype}->{$lossiontype}->{ $index1[$i] } + $betax + ( $charge - 1 ) * $Hatom ) / $charge;
							my $iontypestring = join "_", "alpha", "xlink", $iontype, $lossiontype, "plus$charge";
							$iontable{$iontypestring}->{ $index1[$i] } = $ionmz;
							my $flationtypestring = join "_", "alpha", "xlink", $iontype, $index1[$i], $lossiontype, "plus$charge";
							$flationtable{$flationtypestring} = $ionmz;
							if ( $ionmz > $minionsize && $ionmz < $maxionsize )
							{
								push @{ $xlinkions_alpha{$charge} }, $ionmz;
								push @{ $xlinkions{$charge} },       $ionmz;
								$nxlinkionsalpha++;
							}
						}
					}
				}
			}
		}
		foreach my $iontype (@$fragmenttypes2)
		{
			foreach my $lossiontype (@$losstypes2)
			{

				#print
				#"charge $charge iontype $iontype losstype $lossiontype \n";
				#print "index1 @index1 xlinkpos1 = $xlinkpos1\n";
				if ( $iontype =~ /$fwdiontypes/ )
				{
					for ( $i = $xlinkpos2_indexstart ; $i <= $#index2 ; $i++ )
					{

						#print "index1 @index1 xlinkpos1 = $xlinkpos1 $i index:", $index1[$i],"\n";
						if ( defined( $pepobjiontable2->{$iontype}->{$lossiontype}->{ $index2[$i] } ) )
						{
							my $ionmz = ( $pepobjiontable2->{$iontype}->{$lossiontype}->{ $index2[$i] } + $alphax + ( $charge - 1 ) * $Hatom ) / $charge;
							my $iontypestring = join "_", "beta", "xlink", $iontype, $lossiontype, "plus$charge";
							$iontable{$iontypestring}->{ $index2[$i] } = $ionmz;
							my $flationtypestring = join "_", "beta", "xlink", $iontype, $index2[$i], $lossiontype, "plus$charge";
							$flationtable{$flationtypestring} = $ionmz;
							if ( $ionmz > $minionsize && $ionmz < $maxionsize )
							{
								push @{ $xlinkions_beta{$charge} }, $ionmz;
								push @{ $xlinkions{$charge} },      $ionmz;
								$nxlinkionsbeta++;
							}
						}
					}
				} elsif ( $iontype =~ /$reviontypes/ )
				{
					for ( $i = $revxlinkpos2_indexstart ; $i <= $#index2 ; $i++ )
					{
						if ( defined( $pepobjiontable2->{$iontype}->{$lossiontype}->{ $index2[$i] } ) )
						{
							my $ionmz = ( $pepobjiontable2->{$iontype}->{$lossiontype}->{ $index2[$i] } + $alphax + ( $charge - 1 ) * $Hatom ) / $charge;
							my $iontypestring = join "_", "beta", "xlink", $iontype, $lossiontype, "plus$charge";
							$iontable{$iontypestring}->{ $index2[$i] } = $ionmz;
							my $flationtypestring = join "_", "beta", "xlink", $iontype, $index2[$i], $lossiontype, "plus$charge";
							$flationtable{$flationtypestring} = $ionmz;
							if ( $ionmz > $minionsize && $ionmz < $maxionsize )
							{
								push @{ $xlinkions_beta{$charge} }, $ionmz;
								push @{ $xlinkions{$charge} },      $ionmz;
								$nxlinkionsbeta++;
							}
						}
					}
				}
			}
		}
	}
	if ($verbose)
	{
		foreach my $entry ( sort keys %flationtable )
		{
			print "$entry\t", $flationtable{$entry}, "\n";
		}
	}

	#$self->{'commonions'}      = \@commonions;
	$self->{'commonions'}       = \@commonions;
	$self->{'commonions_alpha'} = \@commonions_alpha;
	$self->{'commonions_beta'}  = \@commonions_beta;
	$self->{'xlinkions'}        = \%xlinkions;
	$self->{'xlinkions_alpha'}  = \%xlinkions_alpha;
	$self->{'xlinkions_beta'}   = \%xlinkions_beta;
	$self->{'ncommonionsalpha'} = $ncommonionsalpha;
	$self->{'ncommonionsbeta'}  = $ncommonionsbeta;
	$self->{'nxlinkionsalpha'}  = $nxlinkionsalpha;
	$self->{'nxlinkionsbeta'}   = $nxlinkionsbeta;

	#$self->{'flationtable_norev'}        = \%flationtable;
	#print Dumper (\%flationtable);
	#exit;
	my %ionannotationtable = reverse %flationtable;
	$self->{'flationtable'}        = \%flationtable;
	$self->{'ion_annotationtable'} = \%ionannotationtable;
	$self->{'iontable'}            = \%iontable;
}

sub _makeiontable_monolink
{
	my $self             = shift;
	my $pepobjs          = $self->gethitObjs;
	my $pepObj1          = $pepobjs->[0];
	my $verbose          = $self->verbose;
	my $xlinkpositions   = $self->gettopology;
	my $xlinkpos1        = $xlinkpositions->[0];
	my $reversexlinkpos1 = $self->getreversetopology( $pepObj1, $xlinkpos1 );
	my $sequence1        = $pepObj1->seq;

	#my $xlinkermass             = $self->getxlinkermass;
	my $xlinkermass         = $self->gettotalxlinkermass;
	my $minioncharge_common = $self->minioncharge_common;
	my $maxioncharge_common = $self->maxioncharge_common;
	my $Hatom               = $pepObj1->getHatom;
	my $PARAMS              = $self->getParams;
	my $pepobjiontable1     = $pepObj1->getiontable;
	my $fragmenttypes1      = $pepObj1->getfragmenttypes;
	my $losstypes1          = $pepObj1->getlosstypes;       #sort ( keys %{ $pepobjiontable1->{ $fragmenttypes[0] } } );
	$verbose && print "losstypes1 @$losstypes1\n";
	my @index =
	  sort { $a <=> $b }
	  keys %{ $pepobjiontable1->{ $fragmenttypes1->[0] }->{'standard'} };
	my %iontable;
	my %flationtable;
	my @commonions              = ();
	my %xlinkions               = ();
	my $xlinkpos1_indexstart    = 0;
	my $xlinkpos2_indexstart    = 0;
	my $revxlinkpos1_indexstart = 0;
	my $rexlinkpos2_indexstart  = 0;
	my $fwdiontypes             = $PARAMS->{'fwd_ions'};
	my $reviontypes             = $PARAMS->{'rev_ions'};
	my $minionsize              = $PARAMS->{'minionsize'};
	my $maxionsize              = $PARAMS->{'maxionsize'};
	my $ncommonionsalpha        = 0;
	my $nxlinkionsalpha         = 0;

	for my $charge ( $minioncharge_common .. $maxioncharge_common )
	{
		$xlinkions{$charge} = [];
	}
	foreach my $iontype (@$fragmenttypes1)
	{
		if ( $iontype =~ /$fwdiontypes/ )
		{
			for my $charge ( $minioncharge_common .. $maxioncharge_common )
			{
				foreach my $lossiontype (@$losstypes1)
				{
					my $i = 0;
					while ( $index[$i] < $xlinkpos1 )
					{
						if ( defined( $pepobjiontable1->{$iontype}->{$lossiontype}->{ $index[$i] } ) )
						{
							my $ionmz = ( $pepobjiontable1->{$iontype}->{$lossiontype}->{ $index[$i] } + ( $charge - 1 ) * $Hatom ) / $charge;
							my $iontypestring = join "_", "alpha", "common", $iontype, $lossiontype, "plus$charge";
							$iontable{$iontypestring}->{ $index[$i] } = $ionmz;
							my $flationtypestring = join "_", "alpha", "common", $iontype, $index[$i], $lossiontype, "plus$charge";
							$flationtable{$flationtypestring} = $ionmz;
							if ( $ionmz > $minionsize && $ionmz < $maxionsize )
							{
								push @commonions, $ionmz;
								$ncommonionsalpha++;
							}
						}
						$i++;
					}
					for $i ( $i .. $#index )
					{
						if ( defined( $pepobjiontable1->{$iontype}->{$lossiontype}->{ $index[$i] } ) )
						{
							my $ionmz = ( $pepobjiontable1->{$iontype}->{$lossiontype}->{ $index[$i] } + $xlinkermass + ( $charge - 1 ) * $Hatom ) / $charge;
							my $iontypestring = join "_", "alpha", "xlink", $iontype, $lossiontype, "plus$charge";
							$iontable{$iontypestring}->{ $index[$i] } = $ionmz;
							my $flationtypestring = join "_", "alpha", "xlink", $iontype, $index[$i], $lossiontype, "plus$charge";
							$flationtable{$flationtypestring} = $ionmz;
							if ( $ionmz > $minionsize && $ionmz < $maxionsize )
							{
								push @{ $xlinkions{$charge} }, $ionmz;
								$nxlinkionsalpha++;
							}
						}
					}
				}
			}
		}
		if ( $iontype =~ /$reviontypes/ )
		{
			for my $charge ( $minioncharge_common .. $maxioncharge_common )
			{
				foreach my $lossiontype (@$losstypes1)
				{
					my $i = 0;
					while ( $index[$i] < $reversexlinkpos1 )
					{
						if ( defined( $pepobjiontable1->{$iontype}->{$lossiontype}->{ $index[$i] } ) )
						{
							my $ionmz = ( $pepobjiontable1->{$iontype}->{$lossiontype}->{ $index[$i] } + ( $charge - 1 ) * $Hatom ) / $charge;
							my $iontypestring = join "_", "alpha", "common", $iontype, $lossiontype, "plus$charge";
							$iontable{$iontypestring}->{ $index[$i] } = $ionmz;
							my $flationtypestring = join "_", "alpha", "common", $iontype, $index[$i], $lossiontype, "plus$charge";
							$flationtable{$flationtypestring} = $ionmz;
							if ( $ionmz > $minionsize && $ionmz < $maxionsize )
							{
								push @commonions, $ionmz;
								$ncommonionsalpha++;
							}
						}
						$i++;
					}
					for $i ( $i .. $#index )
					{
						if ( defined( $pepobjiontable1->{$iontype}->{$lossiontype}->{ $index[$i] } ) )
						{
							my $ionmz = ( $pepobjiontable1->{$iontype}->{$lossiontype}->{ $index[$i] } + $xlinkermass + ( $charge - 1 ) * $Hatom ) / $charge;
							my $iontypestring = join "_", "alpha", "xlink", $iontype, $lossiontype, "plus$charge";
							$iontable{$iontypestring}->{ $index[$i] } = $ionmz;
							my $flationtypestring = join "_", "alpha", "xlink", $iontype, $index[$i], $lossiontype, "plus$charge";
							$flationtable{$flationtypestring} = $ionmz;
							if ( $ionmz > $minionsize && $ionmz < $maxionsize )
							{
								push @{ $xlinkions{$charge} }, $ionmz;
								$nxlinkionsalpha++;
							}
						}
					}
				}
			}
		}
	}
	if ($verbose)
	{
		foreach my $entry ( sort keys %flationtable )
		{
			print "$entry\t", $flationtable{$entry}, "\n";
		}
	}
	my %ionannotationtable = reverse %flationtable;
	$self->{'flationtable'}        = \%flationtable;
	$self->{'commonions'}          = \@commonions;
	$self->{'xlinkions'}           = \%xlinkions;
	$self->{'ncommonionsalpha'}    = $ncommonionsalpha;
	$self->{'ncommonionsbeta'}     = 0;
	$self->{'nxlinkionsalpha'}     = $nxlinkionsalpha;
	$self->{'nxlinkionsbeta'}      = 0;
	$self->{'ion_annotationtable'} = \%ionannotationtable;
	$self->{'iontable'}            = \%iontable;
}

sub drawxlinkspec_alpha_beta
{
	my $self          = shift;
	my $logscale      = shift;
	my $lossions      = shift;
	my $min           = shift;
	my $max           = shift;
	my $xlinkspecfile = shift;
	my $labelpeaks    = shift;
	my $showstructure = shift;
	my $printtextfile = shift;
	return 0 unless $self->xlinktype eq "xlink";
	my $spectrum = $self->getSpecObj;
	my $PARAMS   = $self->getParams;

	if ( !defined($min) )
	{
		$min = $PARAMS->{'minionsize'};
	}
	if ( !defined($max) )
	{
		$max = $PARAMS->{'maxionsize'};
	}
	my $specplot  = specplot->new();
	my @xlinkions = ();
	my $charge;
	my $commonpairs = $spectrum->getcommonpairs;
	if ($logscale)
	{
		my @tmp = ();
		foreach my $pairs (@$commonpairs)
		{
			push @tmp, [ $pairs->[0], sqrt( $pairs->[1] ) ];
		}
		$commonpairs = \@tmp;
	}
	my $mincharge = $self->minioncharge_xlinks;
	my $maxcharge = $self->maxioncharge_xlinks;
	for $charge ( $mincharge .. $maxcharge )
	{
		push @xlinkions, @{ $spectrum->getxlinkpairs($charge) };
	}
	push @xlinkions, @{ $spectrum->getxlinkpairs(0) };
	my $xlinkpeakpairs;
	if ($logscale)
	{
		my @tmp = ();
		foreach my $pairs (@xlinkions)
		{
			push @tmp, [ $pairs->[0], sqrt( $pairs->[1] ) ];
		}
		$xlinkpeakpairs = \@tmp;
	} else
	{
		$xlinkpeakpairs = \@xlinkions;
	}

	#$specplot->plotdata( $min, $max,["black","black"], $commonpairs, $xlinkpeakpairs );
	my ( @matchpairs_alpha, @matchpairs_beta );
	my @matches_alpha = ( @{ $self->getcommonmatches_alpha }, @{ $self->getxlinkmatches_alpha } );
	foreach my $matchedionmz (@matches_alpha)
	{
		my $matchedpeak = $self->get_ion2peakmatch($matchedionmz);    ## is now an arrayref
		foreach my $matchedpeak (@$matchedpeak)
		{

			#my $label       = $self->annotate_matchedion($matchedionmz);
			my $label = $self->annotate_matchedpeak($matchedpeak);
			$label =~ s/plus/\+/;
			$label =~ s/alpha/a/;
			$label =~ s/beta/b/;
			$label =~ s/standard_//;
			my $peakintensity = 0;
			if ($logscale)
			{
				$peakintensity = sqrt( $spectrum->get_ionintensity($matchedpeak) );
			} else
			{
				$peakintensity = $spectrum->get_ionintensity($matchedpeak);
				#$peakintensity = $spectrum->get_realionintensity($matchedpeak);
			}
			#push @matchpairs_alpha, [ $matchedionmz, $peakintensity, int($matchedionmz) . "_" . $label ];
			push @matchpairs_alpha, [ $matchedpeak, $peakintensity, int($matchedpeak) . "_" . $label ];
		}
	}
	my @matches_beta = ( @{ $self->getcommonmatches_beta }, @{ $self->getxlinkmatches_beta } );
	foreach my $matchedionmz (@matches_beta)
	{
		my $matchedpeak = $self->get_ion2peakmatch($matchedionmz);
		foreach my $matchedpeak (@$matchedpeak)
		{

			#my $label       = $self->annotate_matchedion($matchedionmz);
			my $label = $self->annotate_matchedpeak($matchedpeak);
			$label =~ s/plus/\+/;
			$label =~ s/alpha/a/;
			$label =~ s/beta/b/;
			$label =~ s/standard_//;
			my $peakintensity = 0;
			if ($logscale)
			{
				$peakintensity = sqrt( $spectrum->get_ionintensity($matchedpeak) );
			} else
			{
				$peakintensity = $spectrum->get_ionintensity($matchedpeak);
				#$peakintensity = $spectrum->get_realionintensity($matchedpeak);
			}
			#push @matchpairs_beta, [ $matchedionmz, $peakintensity, int($matchedionmz) . "_" . $label ];
			push @matchpairs_beta, [ $matchedpeak, $peakintensity, int($matchedpeak) . "_" . $label ];
		}
	}

	#$specplot->setcolor( "black", "black", "blue", "yellow" );
	$specplot->plotdata( $min, $max, [ "black", "black", "blue", "yellow" ], $commonpairs, $xlinkpeakpairs, \@matchpairs_alpha, \@matchpairs_beta );
	if ($labelpeaks)
	{
		$specplot->labelpeaks( [ "annotation", "annotation", "annotation", "annotation" ], [ "blue", "yellow", "blue", "lightblue" ], [ \@matchpairs_alpha, \@matchpairs_beta, ] );
	} else
	{
		$specplot->labelpeaks( [ "diamond", "diamond", "annotation", "annotation" ], [ "blue", "yellow", "blue", "lightblue" ], [ \@matchpairs_alpha, \@matchpairs_beta, ] );
	}
	my %colorhash = (
					  "alpha" => "blue",
					  "beta"  => "yellow",
	);
	my @sortlist = ( "alpha", "beta" );
	$specplot->drawlegend( 600, 20, \%colorhash, \@sortlist );

	#
	#	if ($showstructure) {
	#		$self->drawpepstructure(
	#			\@xlinkions, $commonpairs, $xlinkspecfile,
	#			$logscale,   $lossions,    $specplot->img
	#		);
	#	}
	unless ( defined($xlinkspecfile) )
	{
		my $id = $self->getid;
		$id =~ s/::/_/g;
		$xlinkspecfile = File::Spec->catfile( $self->getParams->{'outputpath'}, ( join "", $spectrum->getspecbasename, $id, ".png" ) );
	}

	#	print "open $xlinkspecfile <br>";
	$specplot->printimage($xlinkspecfile);
	### chmod the image
	return $xlinkspecfile;
}

sub _makeiontable_intralink
{
	my $self                = shift;
	my $pepobjs             = $self->gethitObjs;
	my $pepObj1             = $pepobjs->[0];
	my $verbose             = $self->verbose;
	my $xlinkpositions      = $self->gettopology;
	my $xlinkpos1           = $xlinkpositions->[0];
	my $xlinkpos2           = $xlinkpositions->[1];
	my $reversexlinkpos1    = $self->getreversetopology( $pepObj1, $xlinkpos2 );
	my $reversexlinkpos2    = $self->getreversetopology( $pepObj1, $xlinkpos1 );
	my $sequence1           = $pepObj1->seq;
	my $xlinkermass         = $self->gettotalxlinkermass;
	my $minioncharge_common = $self->minioncharge_common;
	my $maxioncharge_common = $self->maxioncharge_common;
	my $Hatom               = $pepObj1->getHatom;
	my $PARAMS              = $self->getParams;
	my $pepobjiontable1     = $pepObj1->getiontable;
	my $fragmenttypes1      = $pepObj1->getfragmenttypes;
	my $losstypes1          = $pepObj1->getlosstypes;                              #sort ( keys %{ $pepobjiontable1->{ $fragmenttypes[0] } } );
	$verbose && print "losstypes1 @$losstypes1\n";
	my @index =
	  sort { $a <=> $b }
	  keys %{ $pepobjiontable1->{ $fragmenttypes1->[0] }->{'standard'} };
	my %iontable;
	my %flationtable;
	my @commonions              = ();
	my %xlinkions               = ();
	my $xlinkpos1_indexstart    = 0;
	my $xlinkpos2_indexstart    = 0;
	my $revxlinkpos1_indexstart = 0;
	my $rexlinkpos2_indexstart  = 0;
	my $fwdiontypes             = $PARAMS->{'fwd_ions'};
	my $reviontypes             = $PARAMS->{'rev_ions'};
	my $minionsize              = $PARAMS->{'minionsize'};
	my $maxionsize              = $PARAMS->{'maxionsize'};
	my $ncommonionsalpha        = 0;
	my $nxlinkionsalpha         = 0;

	for my $charge ( $minioncharge_common .. $maxioncharge_common )
	{
		$xlinkions{$charge} = [];
	}

	#print $self->xlinktype," mincharge $minioncharge_common maxcharge $maxioncharge_common\n";
	foreach my $iontype (@$fragmenttypes1)
	{
		if ( $iontype =~ /$fwdiontypes/ )
		{
			for my $charge ( $minioncharge_common .. $maxioncharge_common )
			{
				foreach my $lossiontype (@$losstypes1)
				{
					my $i = 0;
					while ( $index[$i] < $xlinkpos1 )
					{
						if ( defined( $pepobjiontable1->{$iontype}->{$lossiontype}->{ $index[$i] } ) )
						{
							my $ionmz = ( $pepobjiontable1->{$iontype}->{$lossiontype}->{ $index[$i] } + ( $charge - 1 ) * $Hatom ) / $charge;
							my $iontypestring = join "_", "alpha", "common", $iontype, $lossiontype, "plus$charge";
							$iontable{$iontypestring}->{ $index[$i] } = $ionmz;
							my $flationtypestring = join "_", "alpha", "common", $iontype, $index[$i], $lossiontype, "plus$charge";
							$flationtable{$flationtypestring} = $ionmz;
							if ( $ionmz > $minionsize && $ionmz < $maxionsize )
							{
								push @commonions, $ionmz;
								$ncommonionsalpha++;
							}
						}
						$i++;
					}
					while ( $index[$i] < $xlinkpos2 )
					{
						$i++;
					}
					for $i ( $i .. $#index )
					{
						if ( defined( $pepobjiontable1->{$iontype}->{$lossiontype}->{ $index[$i] } ) )
						{
							my $ionmz = ( $pepobjiontable1->{$iontype}->{$lossiontype}->{ $index[$i] } + $xlinkermass + ( $charge - 1 ) * $Hatom ) / $charge;
							my $iontypestring = join "_", "alpha", "xlink", $iontype, $lossiontype, "plus$charge";
							$iontable{$iontypestring}->{ $index[$i] } = $ionmz;
							my $flationtypestring = join "_", "alpha", "xlink", $iontype, $index[$i], $lossiontype, "plus$charge";
							$flationtable{$flationtypestring} = $ionmz;
							if ( $ionmz > $minionsize && $ionmz < $maxionsize )
							{
								push @{ $xlinkions{$charge} }, $ionmz;
								$nxlinkionsalpha++;
							}
						}
					}
				}
			}
		}
		if ( $iontype =~ /$reviontypes/ )
		{
			for my $charge ( $minioncharge_common .. $maxioncharge_common )
			{
				foreach my $lossiontype (@$losstypes1)
				{
					my $i = 0;
					while ( $index[$i] < $reversexlinkpos1 )
					{
						if ( defined( $pepobjiontable1->{$iontype}->{$lossiontype}->{ $index[$i] } ) )
						{
							my $ionmz = ( $pepobjiontable1->{$iontype}->{$lossiontype}->{ $index[$i] } + ( $charge - 1 ) * $Hatom ) / $charge;
							my $iontypestring = join "_", "alpha", "common", $iontype, $lossiontype, "plus$charge";
							$iontable{$iontypestring}->{ $index[$i] } = $ionmz;
							my $flationtypestring = join "_", "alpha", "common", $iontype, $index[$i], $lossiontype, "plus$charge";
							$flationtable{$flationtypestring} = $ionmz;
							if ( $ionmz > $minionsize && $ionmz < $maxionsize )
							{
								push @commonions, $ionmz;
								$ncommonionsalpha++;
							}
						}
						$i++;
					}
					while ( $index[$i] < $reversexlinkpos2 )
					{
						$i++;
					}
					for $i ( $i .. $#index )
					{
						if ( defined( $pepobjiontable1->{$iontype}->{$lossiontype}->{ $index[$i] } ) )
						{
							my $ionmz = ( $pepobjiontable1->{$iontype}->{$lossiontype}->{ $index[$i] } + $xlinkermass + ( $charge - 1 ) * $Hatom ) / $charge;
							my $iontypestring = join "_", "alpha", "xlink", $iontype, $lossiontype, "plus$charge";
							$iontable{$iontypestring}->{ $index[$i] } = $ionmz;
							my $flationtypestring = join "_", "alpha", "xlink", $iontype, $index[$i], $lossiontype, "plus$charge";
							$flationtable{$flationtypestring} = $ionmz;
							if ( $ionmz > $minionsize && $ionmz < $maxionsize )
							{
								push @{ $xlinkions{$charge} }, $ionmz;
								$nxlinkionsalpha++;
							}
						}
					}
				}
			}
		}
	}
	if ($verbose)
	{
		foreach my $entry ( sort keys %flationtable )
		{
			print "$entry\t", $flationtable{$entry}, "\n";
		}
	}

	#foreach my $charge (sort {$a<=>$b} keys %xlinkions){
	# my @tmp=@{$xlinkions{$charge}};
	#
	# print "charge $charge @tmp\n";
	#}
	my %ionannotationtable = reverse %flationtable;
	$self->{'flationtable'}        = \%flationtable;
	$self->{'commonions'}          = \@commonions;
	$self->{'xlinkions'}           = \%xlinkions;
	$self->{'ncommonionsalpha'}    = $ncommonionsalpha;
	$self->{'ncommonionsbeta'}     = 0;
	$self->{'nxlinkionsalpha'}     = $nxlinkionsalpha;
	$self->{'nxlinkionsbeta'}      = 0;
	$self->{'ion_annotationtable'} = \%ionannotationtable;
	$self->{'iontable'}            = \%iontable;
}

sub _makeiontable_intralink_old
{
	my $self    = shift;
	my $pepobjs = $self->gethitObjs;
	my $pepObj1 = $pepobjs->[0];
	my $verbose = $self->verbose;
	print "_makeiontable_intralink: \n";
	my $xlinkpositions      = $self->gettopology;
	my $xlinkpos1           = $xlinkpositions->[0];
	my $xlinkpos2           = $xlinkpositions->[1];
	my $reversexlinkpos1    = $self->getreversetopology( $pepObj1, $xlinkpos1 );
	my $reversexlinkpos2    = $self->getreversetopology( $pepObj1, $xlinkpos2 );
	my $sequence1           = $pepObj1->seq;
	my $deltaMr             = $self->getxlinkermass;
	my $minioncharge_common = $self->minioncharge_common;
	my $maxioncharge_common = $self->maxioncharge_common;
	my $Hatom               = $pepObj1->getHatom;
	my $PARAMS              = $self->getParams;
	my $pepobjiontable1     = $pepObj1->getiontable;
	my $fragmenttypes1      = $pepObj1->getfragmenttypes;
	my $losstypes1          = $pepObj1->getlosstypes;                              #sort ( keys %{ $pepobjiontable1->{ $fragmenttypes[0] } } );
	$verbose && print "losstypes1 @$losstypes1\n";
	my @index =
	  sort { $a <=> $b }
	  keys %{ $pepobjiontable1->{ $fragmenttypes1->[0] }->{'standard'} };
	my %iontable;
	my %flationtable;
	my @commonions              = ();
	my %xlinkions               = ();
	my $xlinkpos1_indexstart    = 0;
	my $xlinkpos2_indexstart    = 0;
	my $revxlinkpos1_indexstart = 0;
	my $rexlinkpos2_indexstart  = 0;
	my $fwdiontypes             = $PARAMS->{'fwd_ions'};
	my $reviontypes             = $PARAMS->{'rev_ions'};
	my $minionsize              = $PARAMS->{'minionsize'};
	my $maxionsize              = $PARAMS->{'maxionsize'};

	foreach my $iontype (@$fragmenttypes1)
	{
		if ( $iontype =~ /[abc]/ )
		{
			for my $charge ( $minioncharge_common .. $maxioncharge_common )
			{
				foreach my $lossiontype (@$losstypes1)
				{
					my $i = 0;
					while ( $index[$i] < $xlinkpos1 )
					{
						if ( defined( $pepobjiontable1->{$iontype}->{$lossiontype}->{ $index[$i] } ) )
						{
							my $ionmz = ( $pepobjiontable1->{$iontype}->{$lossiontype}->{ $index[$i] } + ( $charge - 1 ) * $Hatom ) / $charge;
							my $iontypestring = join "_", "alpha", "common", $iontype, $lossiontype, "plus$charge";
							$iontable{$iontypestring}->{ $index[$i] } = $ionmz;
							my $flationtypestring = join "_", "alpha", "common", $iontype, $index[$i], $lossiontype, "plus$charge";
							$flationtable{$flationtypestring} = $ionmz;
							if ( $ionmz > $minionsize && $ionmz < $maxionsize )
							{
								push @commonions, $ionmz;
							}
						}
						$i++;
					}
					while ( $index[$i] < $xlinkpos2 )
					{
						$i++;
					}
					for $i ( $i .. $#index )
					{
						if ( defined( $pepobjiontable1->{$iontype}->{$lossiontype}->{ $index[$i] } ) )
						{
							my $ionmz = ( $pepobjiontable1->{$iontype}->{$lossiontype}->{ $index[$i] } + $deltaMr + ( $charge - 1 ) * $Hatom ) / $charge;
							my $iontypestring = join "_", "alpha", "xlink", $iontype, $lossiontype, "plus$charge";
							$iontable{$iontypestring}->{ $index[$i] } = $ionmz;
							my $flationtypestring = join "_", "alpha", "xlink", $iontype, $index[$i], $lossiontype, "plus$charge";
							$flationtable{$flationtypestring} = $ionmz;
							if (    $ionmz > $minionsize
								 && $ionmz < $maxionsize )
							{
								push @{ $xlinkions{$charge} }, $ionmz;
							}
						}
					}
				}
			}
		}
		if ( $iontype =~ /$reviontypes/ )
		{
			for my $charge ( $minioncharge_common .. $maxioncharge_common )
			{
				foreach my $lossiontype (@$losstypes1)
				{
					my $i = 0;
					while ( $index[$i] < $reversexlinkpos1 )
					{
						if ( defined( $pepobjiontable1->{$iontype}->{$lossiontype}->{ $index[$i] } ) )
						{
							my $ionmz = ( $pepobjiontable1->{$iontype}->{$lossiontype}->{ $index[$i] } + ( $charge - 1 ) * $Hatom ) / $charge;
							my $iontypestring = join "_", "alpha", "common", $iontype, $lossiontype, "plus$charge";
							$iontable{$iontypestring}->{ $index[$i] } = $ionmz;
							my $flationtypestring = join "_", "alpha", "common", $iontype, $index[$i], $lossiontype, "plus$charge";
							$flationtable{$flationtypestring} = $ionmz;
							if ( $ionmz > $minionsize && $ionmz < $maxionsize )
							{
								push @commonions, $ionmz;
							}
						}
						$i++;
					}
					while ( $index[$i] < $reversexlinkpos2 )
					{
						$i++;
					}
					for $i ( $i .. $#index )
					{
						if ( defined( $pepobjiontable1->{$iontype}->{$lossiontype}->{ $index[$i] } ) )
						{
							my $ionmz = ( $pepobjiontable1->{$iontype}->{$lossiontype}->{ $index[$i] } + ( $charge - 1 ) * $Hatom ) / $charge;
							my $iontypestring = join "_", "alpha", "xlink", $iontype, $lossiontype, "plus$charge";
							$iontable{$iontypestring}->{ $index[$i] } = $ionmz;
							my $flationtypestring = join "_", "alpha", "xlink", $iontype, $index[$i], $lossiontype, "plus$charge";
							$flationtable{$flationtypestring} = $ionmz;
							if (    $ionmz > $minionsize
								 && $ionmz < $maxionsize )
							{
								push @{ $xlinkions{$charge} }, $ionmz;
							}
						}
					}
				}
			}
		}
	}
	if ($verbose)
	{
		foreach my $entry ( sort keys %flationtable )
		{
			print "$entry\t", $flationtable{$entry}, "\n";
		}
	}
	my %ionannotationtable = reverse %flationtable;
	$self->{'commonions'}          = \@commonions;
	$self->{'xlinkions'}           = \%xlinkions;
	$self->{'ion_annotationtable'} = \%ionannotationtable;
	$self->{'iontable'}            = \%iontable;
}

sub makeid
{
	my $self      = shift;
	my $xlinktype = $self->getxlinktype;
	if ( $xlinktype eq "xlink" )
	{
		$self->_makeid_xlink;
	} elsif ( $xlinktype eq "monolink" )
	{
		$self->_makeid_monolink;
	} elsif ( $xlinktype eq "intralink" )
	{
		$self->_makeid_intralink;
	}
}

sub _makeid_xlink
{
	my $self     = shift;
	my $topology = $self->gettopology;
	my $hitobjs  = $self->gethitObjs;
	my $hit1     = $hitobjs->[0];
	my $hit2     = $hitobjs->[1];
	my @tmp      = ( ( join ".", $hit1->id, $topology->[0] ), ( join ".", $hit2->id, $topology->[1] ) );
	@tmp = sort @tmp;
	return join "", $hit1->id, "-", $hit2->id, "-", $self->topologystring;
}

sub _makeid_xlink_andid
{
	my $self     = shift;
	my $topology = $self->gettopology;
	my $hitobjs  = $self->gethitObjs;
	my $hit1     = $hitobjs->[0];
	my $hit2     = $hitobjs->[1];
	my @tmp      = ( ( join ".", $hit1->id, $topology->[0] ), ( join ".", $hit2->id, $topology->[1] ) );
	@tmp = sort @tmp;
	return join "", $hit1->id, "-", $hit2->id, "-", $self->topologystring, " ID1: ", @{ $hit1->protid }, " ID2: ", @{ $hit2->protid },;
}

sub _makeid_monolink
{
	my $self             = shift;
	my $monolinkposition = $self->gettopology->[0];
	my $monolinkweight   = $self->getxlinkermass;
	my $hit              = $self->gethitObjs->[0];
	my @seq              = split //, $hit->seq;
	return join "", $hit->id, "-", $seq[ $monolinkposition - 1 ], $monolinkposition, "-", int($monolinkweight);
}

sub _makeid_intralink
{
	my $self     = shift;
	my $topology = $self->gettopology;
	my $hitobjs  = $self->gethitObjs;
	my $hit      = $hitobjs->[0];
	my @seq      = split //, $hit->seq;
	return join "", $hit->id, "-", $seq[ $topology->[0] - 1 ], $topology->[0], "-", $seq[ $topology->[1] - 1 ], $topology->[1];
}

sub getreversetopology
{
	my $self          = shift;
	my $pepobj        = shift;
	my $xlinkposition = shift;
	my $sequence      = $pepobj->seq;
	return length($sequence) - $xlinkposition + 1;
}

sub nonred
{
	my $array = shift;
	my %seen  = ();
	my $entry;
	my @unique = ();
	foreach $entry (@$array)
	{
		push( @unique, $entry ) unless $seen{$entry}++;
	}
	@$array = @unique;
}

sub get_ion2peakmatchhash
{
	my $self = shift;
	return $self->{'ion2peakmatchhash'};
}

sub get_ion2peakmatch
{
	my $self  = shift;
	my $ionmz = shift;
	return $self->{'ion2peakmatchhash'}->{$ionmz};    ## is now an arrayref
}
################scoring related functions ##########
## calculation of the wTIC score #added by TW v 2.0.1a
sub calc_wTIC
{
	my $self      = shift;
	my $xlinktype = $self->getxlinktype;
	my $pepobjs   = $self->gethitObjs;
	my $wTIC      = 0;
	my $wFracA;
	my $wFracB;
	my $maxdigestlength = $self->get_maxaa;
	my $mindigestlength = $self->get_minaa;
	my $pepAlength;
	my $pepBlength;
	my $aatotal;
	my ( $ticA, $ticB );
	my ( $FracA, $FracB, $InvFracA, $InvFracB, $FracMin, $FracMax, $InvMin, $InvMax );
## calculate the min (not needed) an max inverse values
	$FracMin = $mindigestlength / ( $mindigestlength + $maxdigestlength );
	$FracMax = $maxdigestlength / ( $mindigestlength + $maxdigestlength );
	$InvMin  = 1 / $FracMax;
	$InvMax  = 1 / $FracMin;
## get the length of the peptides
	if ( $xlinktype eq "xlink" )
	{
		my $pepObj1   = $pepobjs->[0];
		my $pepObj2   = $pepobjs->[1];
		my $sequence1 = $pepObj1->seq;
		my $sequence2 = $pepObj2->seq;

		#print "Sequence A: $sequence1\n";
		#print "Sequence B: $sequence2\n";
		# calculate the length of the peptides
		$pepAlength = length($sequence1);
		$pepBlength = length($sequence2);

		#print "pepA length: $pepAlength\n";
		#print "pepB length: $pepBlength\n";
## get the TIC % for peptide A and peptide B
		$ticA = $self->get_TICperc_alpha;
		$ticB = $self->get_TICperc_beta;
	} else
	{
## in case of mono and intralinks the lenght of peptide B is set equal as peptide A, but TIC_pepB is set to 0
		my $pepObj1   = $pepobjs->[0];
		my $sequence1 = $pepObj1->seq;
		$pepAlength = length($sequence1);
		$pepBlength = ( $maxdigestlength + $mindigestlength ) - $pepAlength;
## get the TIC % for peptide A and set the TIC for peptide B = 0
## the TIC A is the full TIC
		$ticA = $self->get_TICperc;
		$ticB = 0;
	}
## calculate the wTIC score
## get the number of total AAs
	$aatotal = $pepAlength + $pepBlength;
## calculate the fraction of pep A and b
	$FracA = $pepAlength / $aatotal;
	$FracB = $pepBlength / $aatotal;
## calculate the invers of the fractions
	$InvFracA = 1 / $FracA;
	$InvFracB = 1 / $FracB;
## calculate the weights ==> normalize Inv values to 1
	$wFracA = $InvFracA / $InvMax;
	$wFracB = $InvFracB / $InvMax;
## calculation of the final wTIC
	$wTIC = $wFracA * $ticA + $wFracB * $ticB;
	$self->{'wTIC'} = $wTIC;
## for debugging ##
	#print "################\n";
	#print "xlinktype: $xlinktype \n";
	#print "aa length total: $aatotal pepA: $pepAlength pepB: $pepBlength\n";
	#print "Mindigestlength: $mindigestlength Maxdigestlength: $maxdigestlength FracMax: $FracMax FracMax: $FracMin IvnMax: $InvMax IvnMin: $InvMin \n";
	#print "FracA: $FracA FracB: $FracB InvFraA: $InvFracA InvFracB: $InvFracB\n";
	#print "TIC A: $ticA weightA: $wFracA TIC B: $ticB weightB: $wFracB \n";
	#print "wTIC: $wTIC\n";
	#print "################\n";
	return $wTIC;
}
## may be used for normalization for values from zero to x to values 0 to y
sub normalize
{
	my $maxbefore = shift;
	my $maxafter  = shift;
	my $value     = shift;
	my $normalized;
	if ( $maxbefore == 0 )
	{
		$normalized = 0;
	} else
	{
		$normalized = ( $value * $maxafter ) / $maxbefore;
	}
	return $normalized or die "Cannot normalize value: sub normalize";
}

sub get_maxaa
{
	my $self = shift;
	return $self->getParams->{'maxdigestlength'};
}

sub get_minaa
{
	my $self = shift;
	return $self->getParams->{'mindigestlength'};
}
## EOF calc wTIC
sub get_wTIC
{
	my $self = shift;
	return $self->{'wTIC'};
}

sub get_TICperc
{
	my $self = shift;
	return $self->{'percTICmatched'};
}

sub get_TICperc_alpha
{
	my $self = shift;
	return $self->{'percTICmatched_alpha'};
}

sub get_TICperc_beta
{
	my $self = shift;
	return $self->{'percTICmatched_beta'};
}

sub calc_TICperc
{
	my $self          = shift;
	my $specObj       = $self->getSpecObj;
	my $intensityhash = $specObj->get_ionintensityhash;
	my @matchingpeaks = ( @{ $self->{'xlinkmatchingpeaks'} }, @{ $self->{'commonmatchingpeaks'} }, );
	nonred( \@matchingpeaks );
	my $nmatchinpeaks           = scalar(@matchingpeaks);
	my @matchingpeakintensities = map { $intensityhash->{$_} } @matchingpeaks;
	my $totalmatchintensity     = Statistics::sum( \@matchingpeakintensities );

	#print Dumper (\@matchingpeakintensities);
	my $totalintensity = $specObj->get_total_ionintensity;
	$self->{'totalmatchintensity'} = $totalmatchintensity;
	my $percTICmatched = $totalmatchintensity / $totalintensity;

	#print "Matching intensities total: $totalmatchintensity<br>";
	$self->{'percTICmatched'} = $percTICmatched;
	return $percTICmatched;
}
### TIC alpha and TIC beta have to be corrected if one several ions matches
### to the same peak. The intensity is then splited to A and B
sub calc_TICperc_alpha
{
	my $self          = shift;
	my $specObj       = $self->getSpecObj;
	my $intensityhash = $specObj->get_ionintensityhash;
	my @matchingpeaks = ( @{ $self->{'xlinkmatchingpeaks_alpha'} }, @{ $self->{'commonmatchingpeaks_alpha'} }, );
	nonred( \@matchingpeaks );
	my $nmatchinpeaks = scalar(@matchingpeaks);

	#print $nmatchinpeaks;
	my @matchingpeakintensities = map { $intensityhash->{$_} } @matchingpeaks;
	my $totalmatchintensity = Statistics::sum( \@matchingpeakintensities );

	#print "Matching intensities alpha: $totalmatchintensity<br>";
	my $totalintensity = $specObj->get_total_ionintensity;
	my $percTICmatched = $totalmatchintensity / $totalintensity;
	$self->{'percTICmatched_alpha'} = $percTICmatched;
	return $percTICmatched;
}

sub calc_TICperc_beta
{
	my $self          = shift;
	my $specObj       = $self->getSpecObj;
	my $intensityhash = $specObj->get_ionintensityhash;
	my @matchingpeaks = ( @{ $self->{'xlinkmatchingpeaks_beta'} }, @{ $self->{'commonmatchingpeaks_beta'} }, );
	nonred( \@matchingpeaks );
	my $nmatchinpeaks           = scalar(@matchingpeaks);
	my @matchingpeakintensities = map { $intensityhash->{$_} } @matchingpeaks;
	my $totalmatchintensity     = Statistics::sum( \@matchingpeakintensities );

	#print "Matching intensities beta: $totalmatchintensity<br>";
	my $totalintensity = $specObj->get_total_ionintensity;
	my $percTICmatched = $totalmatchintensity / $totalintensity;
	$self->{'percTICmatched_beta'} = $percTICmatched;
	return $percTICmatched;
}

sub normalize_TICab
{
	my $self   = shift;
	my $ticA   = $self->{'percTICmatched_alpha'};
	my $ticB   = $self->{'percTICmatched_beta'};
	my $TICtot = $self->{'percTICmatched'};
	my $TICsum = $ticA + $ticB;
## normalize function
## param 1 max before normalization
## param 2 max after normalization
## param 3 value
	my $ticAnorm = normalize( $TICsum, $TICtot, $ticA );
	my $ticBnorm = normalize( $TICsum, $TICtot, $ticB );

	#print "Normalize TICa from sum of TICA+B $TICsum to TIC tot value $TICtot value: $ticA to: $ticAnorm\n ";
	#print "Normalize TICb from sum of TICA+B $TICsum to TIC tot value $TICtot value: $ticB to: $ticBnorm\n ";
## overrides tic
	$self->{'percTICmatched_alpha'} = $ticAnorm;
	$self->{'percTICmatched_beta'}  = $ticBnorm;
}

sub get_total_matchIntensity
{
	my $self = shift;
	return $self->{'totalmatchintensity'};
}

sub calcNmatches
{
	my $self     = shift;
	my $nmatches = ( $#{ $self->getxlinkmatches } + $#{ $self->getcommonmatches } + 2 );
	return $nmatches;
}

sub get_meanmatchintensity
{
	my $self = shift;
	return $self->{'meanmatchintensity'};
}

#sub getPRMI {
#	my $self = shift;
#	return $self->{'percentrank_matchintensity'};
#}
sub transform_peaks_to_log
{
	my $self            = shift;
	my $peakarray       = shift;
	my $specobj         = $self->getSpecObj;
	my $transformedhash = {};
	foreach my $peak (@$peakarray)
	{
## get intensity
		#print "<br>";
		my $int = $specobj->get_ionintensity($peak);

		#print "Peak: $peak Int: $int<br>";
		my $intsqrt = sqrt($int);

		#my $intsqrt = $int;
		#print "Peak: $peak Intsqrt: $intsqrt<br>";
		$transformedhash->{$intsqrt} = $peak;
	}
	return $transformedhash;
}

sub sort_hash_into_bins
{
	my $self         = shift;
	my $binarray     = shift;
	my $intensithash = shift;
### sort the bins
	my @sortedbins = sort { $b <=> $a } @$binarray;
### get the hashkeys
	my @intensityarray = sort { $a <=> $b } ( keys %$intensithash );
	my $sortedhash = {};
	for ( my $i = 0 ; $i < scalar(@sortedbins) ; $i++ )
	{

		foreach my $int (@intensityarray)
		{
			if ( ( $int <= $sortedbins[$i] ) && ( $int > $sortedbins[ $i + 1 ] ) )
			{
## then count the peak into that bin
				$sortedhash->{ $sortedbins[$i] }++;
			}
		}
	}
	return $sortedhash;
}

sub sort_ions_into_bins_xlink
{
	my $self = shift;
	my $bins = shift;
### Transform all matched peaks, keys are intensities
	my $commonalphalog = $self->transform_peaks_to_log( $self->{'commonmatchingpeaks_alpha'} );
	my $xlalphalog     = $self->transform_peaks_to_log( $self->{'xlinkmatchingpeaks_alpha'} );
	my $commonbetalog  = $self->transform_peaks_to_log( $self->{'commonmatchingpeaks_beta'} );
	my $xlbetalog      = $self->transform_peaks_to_log( $self->{'xlinkmatchingpeaks_beta'} );
### Sort into bins
	my $commonalphabins = $self->sort_hash_into_bins( $bins, $commonalphalog );
	my $xlalphabins     = $self->sort_hash_into_bins( $bins, $xlalphalog );
	my $commonbetabins  = $self->sort_hash_into_bins( $bins, $commonbetalog );
	my $xlbetabins      = $self->sort_hash_into_bins( $bins, $xlbetalog );
	return ( $commonalphabins, $xlalphabins, $commonbetabins, $xlbetabins );
}

sub sort_ions_into_bins_mono_intra
{
	my $self = shift;
	my $bins = shift;
### Transform all matched peaks, keys are intensities
	my $commonlog = $self->transform_peaks_to_log( $self->{'commonmatchingpeaks'} );
	my $xllog     = $self->transform_peaks_to_log( $self->{'xlinkmatchingpeaks'} );
### Sort into bins
	my $commonbins  = $self->sort_hash_into_bins( $bins, $commonlog );
	my $xlalphabins = $self->sort_hash_into_bins( $bins, $xllog );
	return ( $commonbins, $xlalphabins );
}

sub weighted_matchodds
{
	my $self         = shift;
	my $matchbin     = shift;
	my $binarray     = shift;
	my $weightsarray = shift;
	my $possibleions = shift;
	my $aprioryprob  = shift;
	my $wmatchodds   = 0;
### Go through the bins and calculate the weighted matchodds
	for ( my $i = 0 ; $i < scalar(@$binarray) ; $i++ )
	{

		#print "checking bin: $binarray->[$i]<br>";
		if ( $matchbin->{ $binarray->[$i] } )
		{
			my $prob = -log( 1 - Statistics::cumulativeBinomial( $possibleions, $matchbin->{ $binarray->[$i] }, $aprioryprob ) + 1E-5 );
			$wmatchodds += ( $prob * $weightsarray->[$i] );

			#print "Weighted Modds for bin: $binarray->[$i] weight: $weightsarray->[$i]<br>";
		}
	}
	return $wmatchodds;
}
## Added by TW
sub get_intsum_score
{
	my $self = shift;
## is calculated during TIC perc
	return $self->{'totalmatchintensity'};
}
## Added by TW
sub calc_apriory_probs
{
	my $self             = shift;
	my $pepobjs          = $self->gethitObjs;
	my $PARAMS           = $self->getParams;
	my $matchobject      = $self->getmatchObj;
	my $INFOINDEX_target = $matchobject->get_infoindex('infoindex_target');
	my $INFOINDEX_decoy  = $matchobject->get_infoindex('infoindex_decoy');

	#$INFOINDEX_decoy->{'npeps'}=0;
	#print Dumper ($INFOINDEX_target);
	if ( $PARAMS->{'RuntimeDecoys'} )
	{
		my $INFOINDEX_decoy = $matchobject->get_infoindex('infoindex_decoy');
	}
	my $xlinktype = $self->getxlinktype;
##
	my $protein1;
	my $protein2;
	my $npeptidesA;
	my $npeptidesB;
	my $npeptidestotal = $INFOINDEX_target->{'npeps'} + $INFOINDEX_decoy->{'npeps'};

	#print "Total number of peptides: $npeptidestotal\n";
	if ( $xlinktype eq "xlink" )
	{
### get the protein ids
		my $pepObj1 = $pepobjs->[0];
		my $pepObj2 = $pepobjs->[1];
		my $pepid1  = $pepObj1->protid;    ## is an array-reference
		my $pepid2  = $pepObj2->protid;
		my ( $type, $combination ) = $self->get_type_of_xlink( $pepid1, $pepid2 );

		#print "Type of xlink: ". ($type)."\n";
		if ( $type eq "intra/inter xl" || $type eq "intra-protein xl" )
		{
## check which protein is equal
			$protein1 = $combination->{'intra'}->{'p1'};
			$protein2 = $combination->{'intra'}->{'p2'};
		} else
		{
## check which protein is equal
			$protein1 = $combination->{'inter'}->{'p1'};
			$protein2 = $combination->{'inter'}->{'p2'};
		}
### Get the numbers of peptides
		if ( $protein1 =~ m/decoy/ )
		{
			$npeptidesA = $INFOINDEX_decoy->{'npeptideshash'}->{$protein1};
		} else
		{
			$npeptidesA = $INFOINDEX_target->{'npeptideshash'}->{$protein1};
		}
		if ( $protein2 =~ m/decoy/ )
		{
			$npeptidesB = $INFOINDEX_decoy->{'npeptideshash'}->{$protein2};
		} else
		{
			$npeptidesB = $INFOINDEX_target->{'npeptideshash'}->{$protein2};
		}
### Calc the prob for the first draw
		my $apriori_prob_intra = ( $npeptidesA / $npeptidestotal ) * ( $npeptidesB / $npeptidestotal );
		my $apriori_prob_inter = 1 - $apriori_prob_intra;
		if ( $type eq "intra/inter xl" || $type eq "intra-protein xl" )
		{
			return $apriori_prob_intra;

			#print "Apriory probability: $apriori_prob_intra\n";
		} else
		{
			return $apriori_prob_inter;

			#print "Apriory probability: $apriori_prob_inter\n";
		}

		#print "NPeptides Protein A: ". $npeptidesA."\n";
		#print "NPeptides Protein B: ". $npeptidesB."\n";
	} else
	{

		#print "Apriory probability: 1\n";
## For mono and intralinks
		return 1;
	}
}

sub get_type_of_xlink
{
	my $self   = shift;
	my $spidp1 = shift;
	my $spidp2 = shift;

	#my @prots1 = split( ",", $spidp1 );
	#my @prots2 = split( ",", $spidp2 );
	my $type;
	my $intralink   = 0;
	my $interlink   = 0;
	my $combination = {};

	#print Dumper (\@prots1);
	foreach my $prot1 (@$spidp1)
	{

		#print "HERE!\n";
		#print $prot1."\n";
		foreach my $prot2 (@$spidp2)
		{

			#print $prot1."\n";
			if ( $prot1 eq $prot2 )
			{
				$intralink                      = 1;
				$combination->{'intra'}->{'p1'} = $prot1;
				$combination->{'intra'}->{'p2'} = $prot2;
			} else
			{
				$interlink                      = 1;
				$combination->{'inter'}->{'p1'} = $prot1;
				$combination->{'inter'}->{'p2'} = $prot2;
			}
		}
	}
	if ( $intralink && $interlink )
	{
		$type = "intra/inter xl";
	}
	if ( $intralink && !$interlink )
	{
		$type = "intra-protein xl";
	}
	if ( !$intralink && $interlink )
	{
		$type = "inter-protein xl";
	}
	return $type, $combination;
}
## Added by TW
sub calc_weighted_matchodds
{
	my $self        = shift;
	my $xlinktype   = $self->getxlinktype;
	my $wmatchodds  = 0;
	my $pcommon     = $self->get_apriori_pcommon;
	my $pxlink      = $self->get_apriori_pxlink;
	my $verbose     = 0;
	my $wm_odds_sum = 0;
	$verbose && print "<br>Matchratio Score params: <br>";
	$verbose && print "Apriory common: $pcommon <br>";
	$verbose && print "Apriory xlink: $pxlink <br>";
### define the bins in which the spectrum should be separated
	my $nbins = 5;
## bins have to be sorted from largest to lowest
	my @bins = ( 10, 5, 2.5, 1.25, 0.625, 0.3125 );
## weights have to have the same order than bins
	my @weights = ( 1 / 2, 1 / 4, 1 / 8, 1 / 16, 1 / 32, 1 / 64 );

	if ( $xlinktype eq "xlink" )
	{
### Sort the ions into the bins, transform intensities to log scale
		my ( $commonalphabins, $xlalphalogbins, $commonbetalogbins, $xlbetalogbins ) = $self->sort_ions_into_bins_xlink( \@bins );

		#print "<br>Ca Bins:<br>";
		#print Dumper($commonalphabins);
		#print "<br>Xa Bins:<br>";
		#print Dumper($xlalphalogbins);
		#print "<br>Cb Bins:<br>";
		#print Dumper($commonbetalogbins);
		#print "<br>Xb Bins:<br>";
		#print Dumper($xlbetalogbins);
### calculation of weighted matchodds for the alpha chain
		my $ncionsalpha      = $self->get_number_of_CommonIons_alpha;
		my $nxionsalpha      = $self->get_number_of_XlinkIons_alpha;
		my $commonalphamodds = $self->weighted_matchodds( $commonalphabins, \@bins, \@weights, $ncionsalpha, $pcommon );
		my $xlalphamodds     = $self->weighted_matchodds( $xlalphalogbins, \@bins, \@weights, $nxionsalpha, $pxlink );
		$verbose && print "wMatchodds aC: $commonalphamodds<br>";
		$verbose && print "wMatchodds aX: $xlalphamodds<br>";
### calculation of weighted matchodds for the beta chain
		my $ncionsbeta      = $self->get_number_of_CommonIons_beta;
		my $nxionsbeta      = $self->get_number_of_XlinkIons_beta;
		my $commonbetamodds = $self->weighted_matchodds( $commonbetalogbins, \@bins, \@weights, $ncionsbeta, $pcommon );
		my $xlbetamodds     = $self->weighted_matchodds( $xlbetalogbins, \@bins, \@weights, $nxionsbeta, $pxlink );
		$verbose && print "wMatchodds bC: $commonbetamodds<br>";
		$verbose && print "wMatchodds bX: $xlbetamodds<br>";
		$wmatchodds = ( ( $commonalphamodds + $xlalphamodds ) / 2 + ( $commonbetamodds + $xlbetamodds ) / 2 );
		$wm_odds_sum = ( $commonalphamodds / 4 + $xlalphamodds / 4 + $commonbetamodds / 4 + $xlbetamodds / 4 );
		$verbose && print "wMatchodds Sum: $wmatchodds<br>";
	} else
	{
		### for mono and intralinks
		my ( $commonbins, $xlalphabins ) = $self->sort_ions_into_bins_mono_intra( \@bins );
		my $ncommonions = $self->get_number_of_CommonIons;
		my $commonmodds = $self->weighted_matchodds( $commonbins, \@bins, \@weights, $ncommonions, $pcommon );
		my $nxlinkions  = $self->get_number_of_XlinkIons;
		my $xlinkmodds  = $self->weighted_matchodds( $xlalphabins, \@bins, \@weights, $nxlinkions, $pxlink );
		$wmatchodds  = ( $commonmodds + $xlinkmodds ) / 2;
		$wm_odds_sum = ( $commonmodds / 2 + $xlinkmodds / 2 );
		$verbose && print "wMatchodds Mean: $wmatchodds<br>";
		$verbose && print "wMatchodds Sum: $wmatchodds<br>";
	}
	$self->{'weighted_matchodds_mean'} = $wmatchodds;
	$self->{'weighted_matchodds_sum'}  = $wm_odds_sum;
	return $wmatchodds;
}
### matchration is now matchodds score
sub calcmatchratio
{
	my $self       = shift;
	my $xlinktype  = $self->getxlinktype;
	my $matchratio = 0;
	my $pcommon    = $self->get_apriori_pcommon;
	my $pxlink     = $self->get_apriori_pxlink;
	if ( $xlinktype eq "xlink" )
	{
		my $nxmatchesalpha         = $self->get_number_of_Xlinkmatches_alpha;
		my $ncmatchesalpha         = $self->get_number_of_Commonmatches_alpha;
		my $ncionsalpha            = $self->get_number_of_CommonIons_alpha;
		my $nxionsalpha            = $self->get_number_of_XlinkIons_alpha;
		my $ncionsbeta             = $self->get_number_of_CommonIons_beta;
		my $nxionsbeta             = $self->get_number_of_XlinkIons_beta;
		my $nxmatchesbeta          = $self->get_number_of_Xlinkmatches_beta;
		my $ncmatchesbeta          = $self->get_number_of_Commonmatches_beta;
		my $probabilityalphacommon = -log( 1 - Statistics::cumulativeBinomial( $ncionsalpha, $ncmatchesalpha, $pcommon ) + 1E-5 );
		$self->{'matchratio_alphacommon'} = $probabilityalphacommon;
		my $probabilityalphaxlink = -log( 1 - Statistics::cumulativeBinomial( $nxionsalpha, $nxmatchesalpha, $pxlink ) + 1E-5 );
		$self->{'matchratio_alphaxlink'} = $probabilityalphaxlink;
		my $probabilitybetacommon = -log( 1 - Statistics::cumulativeBinomial( $ncionsbeta, $ncmatchesbeta, $pcommon ) + 1E-5 );
		$self->{'matchratio_betacommon'} = $probabilitybetacommon;
		my $probabilitybetaxlink = -log( 1 - Statistics::cumulativeBinomial( $nxionsbeta, $nxmatchesbeta, $pxlink ) + 1E-5 );
		$self->{'matchratio_betaxlink'} = $probabilitybetaxlink;

		#print "nxmatchesalpha:$nxmatchesalpha of $nxionsalpha<br>";
		#print "ncmachesalpha: $ncmatchesalpha of $ncionsalpha <br>";
		#print "nxmatchesbeta: $nxmatchesbeta of $nxionsbeta <br>";
		#print "ncmatchesbeta: $ncmatchesbeta of $ncionsbeta <br>";
		$matchratio = ( $probabilityalphacommon + $probabilityalphaxlink + $probabilitybetaxlink + $probabilitybetacommon ) / 4;

		#print "Matchratio: $matchratio<br>";
	} else
	{
		my $nxlinkmatches          = $self->get_number_of_Xlinkmatches;
		my $ncommonmatches         = $self->get_number_of_Commonmatches;
		my $nxlinkions             = $self->get_number_of_XlinkIons;
		my $ncommonions            = $self->get_number_of_CommonIons;
		my $probabilityalphacommon = -log( 1 - Statistics::cumulativeBinomial( $ncommonions, $ncommonmatches, $pcommon ) + 1E-5 );
		$self->{'matchratio_alphacommon'} = $probabilityalphacommon;

		#print "nxlinkions $nxlinkions\n";
		#print "nxlinkmatches $nxlinkmatches\n";
		#print "Pxlink $pxlink\n";
		#print Statistics::cumulativeBinomial( $nxlinkions, $nxlinkmatches, $pxlink );
		my $probabilityalphaxlink = -log( 1 - Statistics::cumulativeBinomial( $nxlinkions, $nxlinkmatches, $pxlink ) + 1E-5 );

		#print "<br>Prob xl: $probabilityalphaxlink<br>";
		$self->{'matchratio_alphaxlink'} = $probabilityalphaxlink;
		$self->{'matchratio_betacommon'} = 0;
		$self->{'matchratio_betaxlink'}  = 0;
		$matchratio                      = ( $probabilityalphacommon + $probabilityalphaxlink ) / 2;
	}
	$self->{'matchratio'} = $matchratio;
	return $matchratio;
}

sub get_matchratio_alphacommon
{
	my $self = shift;
	return $self->{'matchratio_alphacommon'};
}

sub get_matchratio_alphaxlink
{
	my $self = shift;
	return $self->{'matchratio_alphaxlink'};
}

sub get_matchratio_betaxlink
{
	my $self = shift;
	return $self->{'matchratio_betaxlink'};
}

sub get_matchratio_betacommon
{
	my $self = shift;
	return $self->{'matchratio_betacommon'};
}

sub calcprescore
{
	my $self      = shift;
	my $xlinktype = $self->getxlinktype;
	my $prescore  = 0;
	if ( $xlinktype eq "xlink" )
	{
		my $alpha = ( $self->get_number_of_Xlinkmatches_alpha + $self->get_number_of_Commonmatches_alpha ) / ( $self->get_number_of_CommonIons_alpha + $self->get_number_of_XlinkIons_alpha );
		my $beta  = ( $self->get_number_of_Xlinkmatches_beta + $self->get_number_of_Commonmatches_beta ) /   ( $self->get_number_of_CommonIons_beta + $self->get_number_of_XlinkIons_beta );
		$prescore = sqrt( $alpha * $beta );
	} else
	{
		$prescore = ( $self->get_number_of_Xlinkmatches + $self->get_number_of_Commonmatches ) / ( $self->get_number_of_XlinkIons + $self->get_number_of_CommonIons );
	}
	return $prescore;
}

sub get_num_of_matched_ions_alpha
{
	my $self      = shift;
	my $xlinktype = $self->getxlinktype;
	my $numalpha  = 0;
	if ( $xlinktype eq "xlink" )
	{
		$numalpha = $self->get_number_of_Xlinkmatches_alpha + $self->get_number_of_Commonmatches_alpha;
	} else
	{
		$numalpha = $self->get_number_of_Xlinkmatches + $self->get_number_of_Commonmatches;
	}
	return $numalpha;
}

sub get_num_of_matched_ions_beta
{
	my $self      = shift;
	my $xlinktype = $self->getxlinktype;
	my $numalpha  = 0;
	if ( $xlinktype eq "xlink" )
	{
		$numalpha = $self->get_number_of_Xlinkmatches_beta + $self->get_number_of_Commonmatches_beta;
	} else
	{

		#my $numalpha=( $self->get_number_of_Xlinkmatches + $self->get_number_of_Commonmatches);
		$numalpha = "-";
	}
	return $numalpha;
}

sub calcprescore_alpha
{
	my $self      = shift;
	my $xlinktype = $self->getxlinktype;
	my $prescore  = 0;
	if ( $xlinktype eq "xlink" )
	{
		my $alpha = ( $self->get_number_of_Xlinkmatches_alpha + $self->get_number_of_Commonmatches_alpha ) / ( $self->get_number_of_CommonIons_alpha + $self->get_number_of_XlinkIons_alpha );
		return $alpha;
	} else
	{
		$prescore = ( $self->get_number_of_Xlinkmatches + $self->get_number_of_Commonmatches ) / ( $self->get_number_of_XlinkIons + $self->get_number_of_CommonIons );
		return $prescore;
	}
}

sub calcprescore_beta
{
	my $self      = shift;
	my $xlinktype = $self->getxlinktype;
	my $prescore  = 0;
	if ( $xlinktype eq "xlink" )
	{
		my $beta = ( $self->get_number_of_Xlinkmatches_beta + $self->get_number_of_Commonmatches_beta ) / ( $self->get_number_of_CommonIons_beta + $self->get_number_of_XlinkIons_beta );
		return $beta;
	} else
	{
		$prescore = ( $self->get_number_of_Xlinkmatches + $self->get_number_of_Commonmatches ) / ( $self->get_number_of_XlinkIons + $self->get_number_of_CommonIons );
	}
	return 0;
}

sub check_for_series
{
	my $self            = shift;
	my $series          = shift;
	my $keysions        = shift;
	my $fulltable       = $self->get_full_iontable;
	my $matchedionshash = $self->get_flat_ionmatchhash;

	#print Dumper ($fulltable);
	# Iontable: hash with 'alpha_y' => { '6' => [ '712.3742108', '356.691017916', '238.129953621333' ],
	# 6 is the position
### get all keys, keys are alpha_y, alpha_b, beta_y, and beta_b
	my @keysa_b = sort { $a <=> $b } keys %{ $fulltable->{$keysions} };
## keys are the ionpositions
	#print Dumper (@keysa_b);
### get all hits // ions
### Iterate from one key to the next
### count if a hit is found
	my $lastfound   = 0;
	my $seriesfound = 0;
### go through the whole ionseries
	foreach my $key (@keysa_b)
	{
### check if lastfound is one before the one right now
		if ( $lastfound == ( $key - 1 ) )
		{

			#print "Series may continue!<br>";
## the series may continue
		} else
		{
			$seriesfound = 0;
		}

		#print "KEY IS $key<br>";
### Get the th ions
		my $ions = $fulltable->{$keysions}->{$key};

		#print Dumper ($ions);
		my $matched = 0;
		foreach my $ion (@$ions)
		{

			#print "Check ion : $ion<br>";
## check if this is matched
			if ( $matchedionshash->{$ion} )
			{

				#print "Ion $ion is matched on position $key!<br>";
				$matched   = 1;
				$lastfound = $key;
				$seriesfound++;
				### Go to the next position if one ion of a certain position is found to be matched
				### Do not count multiple matched z states as additional ions in the series.
				### Counts series as seen in the tickmarks on the peptide plot
				last;
			}
		}

		#print "last key is:", $#keysa_b, "key is $key<br>";
		if ( !$matched || $key == $#keysa_b + 1 )
		{
## then it is a gap here or it is the last key of the peptides ions then also check if there is a series
			#print "GAP on position: $key<br>";
			#print "The series before was $seriesfound but is now interrupted<br>";
### Here the series can be reported if a gap is found or if it is the last element
### only record if it is a non 0 series
			unless ( $seriesfound == 0 )
			{
				$series->{$seriesfound}++;
			}
### reset the series
			$seriesfound = 0;
		} else
		{

			#print "Series so far: $seriesfound<br>";
		}
	}

	#print Dumper ($series);
}

sub calc_series_score_sub
{
	my $self     = shift;
	my $numalpha = shift;
	my $series   = shift;
	my $max      = ( $numalpha**2 ) * 2;    ## 2 ionseries, maxscore for one series is: num^2

	#print "MAX is: $numalpha^2 * 2: $max<br>";
	if ( $max == 0 )
	{
		return 0;
	}
	my $seriesscorealpha = 0;
	foreach my $key ( keys %$series )
	{

		#print "KEY IS: $key";
		my $score = $key * $key;
## check how often seen
		my $scoretotal = $score * $series->{$key};

		#print "Subscore: $scoretotal, = $score * numberseen<br>";
		$seriesscorealpha = $seriesscorealpha + $scoretotal;
	}

	#print "Score total is: $seriesscorealpha<br>";
## Normalize to 100
	my $normalizedscore = $seriesscorealpha * 100 / $max;

	#print "Score normalized is: ",$normalizedscore;
	return $normalizedscore;
}

sub calc_series_score
{
	my $self      = shift;
	my $xlinktype = $self->getxlinktype;
	my $fulltable = $self->get_full_iontable;

	#print Dumper ($fulltable);
	if ( $xlinktype eq "xlink" )
	{

		#print Dumper ($hash);
		my $series = {};
		$self->check_for_series( $series, "alpha_b" );
		$self->check_for_series( $series, "alpha_y" );
		### calculate the series score for the first peptide
		### get the number of AAs
		my @keysa_b = sort { $a <=> $b } keys %{ $fulltable->{alpha_b} };

		#print Dumper ($fulltable);
		my $numalpha = scalar(@keysa_b);
		my $scorealpha = $self->calc_series_score_sub( $numalpha, $series );

		#print "<br>Score alpha: ",$scorealpha,"<br>";
		#print "<br>Series alpha:<br>";
		#print Dumper ($series);
		#print "Maxnumalpha: $numalpha<br>";
		$series = {};
		$self->check_for_series( $series, "beta_b" );
		$self->check_for_series( $series, "beta_y" );

		#print "Series beta:<br>";
		#print Dumper ($series);
		#exit;
		my @keysb_b   = sort { $a <=> $b } keys %{ $fulltable->{beta_b} };
		my $numbeta   = scalar(@keysb_b);
		my $scorebeta = $self->calc_series_score_sub( $numbeta, $series );

		#print "Score beta: ",$scorebeta,"<br>";
		#print "Series beta:<br>";
		#print Dumper ($series);
		#print "<br>";
		my $fullscore = ( $scorealpha + $scorebeta ) / 2;
		$self->{'series_score_mean'} = $fullscore;

		#print "<br>Full Score: $fullscore<br>";
	} else
	{
		### Calc the score for monolinks
		my $series = {};
		$self->check_for_series( $series, "alpha_b" );
		$self->check_for_series( $series, "alpha_y" );
		### calculate the series score for the first peptide
		### get the number of AAs
		my @keysa_b    = sort { $a <=> $b } keys %{ $fulltable->{alpha_b} };
		my $numalpha   = scalar(@keysa_b);
		my $scorealpha = $self->calc_series_score_sub( $numalpha, $series );

		#print "<br>Score alpha: ",$scorealpha,"<br>";
		$self->{'series_score_mean'} = $scorealpha;
	}
}

sub get_series_score
{
	my $self = shift;
	return $self->{'series_score_mean'};
}

sub calcfullscore
{
	my $self                   = shift;
	my $xcorr_window_tolerance = shift;
	$self->calcmatchratio;    ## matchodds score
	$self->calc_weighted_matchodds;
	### Calc The Series Score
	#print "Here:<br>";
	$self->calc_series_score;

	# xcorr score
	if ( $xcorr_window_tolerance > 0 )
	{
		$self->xcorrelation_common_normalized;
		$self->xcorrelation_xlink_normalized;
	} else
	{
		$self->xcorrelation_common;
		$self->xcorrelation_xlink;
	}

	# TIC score
	$self->calc_TICperc;
	my $xlinktype = $self->getxlinktype;
	if ( $xlinktype eq "xlink" )
	{
		$self->calc_TICperc_alpha;
		$self->calc_TICperc_beta;
		## new function to normalice the TICa and TICb to TICtotal.
		## is necassary because some ions match to same peaks.
		#exit;
		$self->normalize_TICab;
	} else
	{

		# for mono links and loop links
		# TODO:
		$self->{'percTICmatched_alpha'} = $self->{'percTICmatched'};
	}
	## calculate the wTIC
	$self->calc_wTIC;
	## get the int_sum
	#print "Sum of matched peaks: ", $self->get_intsum_score, "\n";
	## calc the apriory score
	unless ( $self->{'apriory_match_probs_log'} )
	{    ### is already set when used in offline mode
		    #print "TRUE";
		    #exit;
		$self->{'apriory_match_probs'} = $self->calc_apriory_probs;
		unless ( $self->{'apriory_match_probs'} )
		{
			$self->{'apriory_match_probs_log'} = 0;
		} else
		{
			$self->{'apriory_match_probs_log'} = log10( $self->{'apriory_match_probs'} );
		}
	} else
	{

		#$self->{'apriory_match_probs_log'}=log10($self->{'apriory_match_probs'});
	}

	#print "Apriory match probabilities: ", $self->calc_apriory_probs,"\n";
	#print "Apriory match probabilities log: ", $self->{'apriory_match_probs_log'},"\n";
	## calculates the weighted score
	$self->calcscore;
}

sub set_lapS
{
	my $self = shift;
	my $aps  = shift;
	$self->{'apriory_match_probs_log'} = $aps;

	#print "Score Set";
	return;
}

sub log10
{
	my $n = shift;
	return log($n) / log(10);
}

sub get_apriory_match_probs_log
{
	my $self = shift;
	return $self->{'apriory_match_probs_log'};
}

sub get_apriory_match_probs
{
	my $self = shift;
	return $self->{'apriory_match_probs'};
}

sub getprescore
{
	my $self = shift;
	return $self->{'prescore'};
}

sub getprescore_alpha
{
	my $self = shift;
	return $self->{'prescore_alpha'};
}

sub getprescore_beta
{
	my $self = shift;
	return $self->{'prescore_beta'};
}

sub get_weightedmatchodds_mean
{
	my $self = shift;
	return $self->{'weighted_matchodds_mean'};
}

sub get_weightedmatchodds_sum
{
	my $self = shift;
	return $self->{'weighted_matchodds_sum'};
}

sub getmatchratio
{
	my $self = shift;
	return $self->{'matchratio'};
}
##############################################
sub getiontable
{
	my $self = shift;
	return $self->{'iontable'};
}

sub match_offline
{
	my $self    = shift;
	my $specobj = shift;
	my %matchhash;
	$self->{'specObj'} = $specobj;
	my $dummy;
	( $self->{'xlinkmatches'},  $dummy ) = $self->matchxlinks( \%matchhash );
	( $self->{'commonmatches'}, $dummy ) = $self->matchcommon( \%matchhash );
	$self->{'ion2peakmatchhash'} = \%matchhash;
	## my %tmp = reverse %matchhash;
	## make peak to ion hash to lookup the the peak ion pairs
	my %tmp;
	foreach my $ion ( keys %matchhash )
	{
		my $arrayref = $matchhash{$ion};
		foreach my $peak (@$arrayref)
		{
			$tmp{$peak} = $ion;
		}
	}
	$self->{'peak2ionmatchhash'} = \%tmp;

	#$self->{'peak2ionmatchhash'} = \%tmp;
}

sub getxlinktype
{
	my $self = shift;
	return $self->{'xlinktype'};
}

sub getxlinkmatches
{
	my $self = shift;
	return $self->{'xlinkmatches'};
}

sub getxlinkmatches_alpha
{
	my $self = shift;
	return $self->{'xlinkmatches_alpha'};
}

sub getxlinkmatches_beta
{
	my $self = shift;
	return $self->{'xlinkmatches_beta'};
}

sub getcommonmatches
{
	my $self = shift;
	return $self->{'commonmatches'};
}

sub getcommonmatches_alpha
{
	my $self = shift;
	return $self->{'commonmatches_alpha'};
}

sub getcommonmatches_beta
{
	my $self = shift;
	return $self->{'commonmatches_beta'};
}

sub gethithash
{
	my $self = shift;
	return $self->{'hithash'};
}

sub getSpecObj
{
	my $self = shift;
	return $self->{'specObj'};
}

sub printiontable
{
	my $self     = shift;
	my $iontable = $self->getiontable;
	foreach my $iontype ( sort keys %$iontable )
	{
		print $iontype, "\t";
		foreach my $position ( sort { $a <=> $b } keys %{ $iontable->{$iontype} } )
		{
			if ( $iontable->{$iontype}->{$position} )
			{
				print "$position:", sprintf( "%.2f", $iontable->{$iontype}->{$position} ), " ";
			} else
			{
				print "$position: - ";
			}
		}
		print "\n";
	}
	print "\n";
}

sub topologystring
{
	my $self         = shift;
	my $combinations = $self->gettopology;
	my $topology     = join "", "a", $combinations->[0];
	$combinations->[1] && ( $topology .= join "", "-b", $combinations->[1] );
	return $topology;
}

sub getid
{
	my $self = shift;
	return $self->{'id'};
}

sub getuniqueid2
{
	my $self = shift;
	return $self->{'uniqueid'};
}

sub getmatchObj
{
	my $self = shift;
	return $self->{'matchobj'};
}

sub gethitObjs
{
	my $self = shift;
	return $self->{'hitpepObjs'};
}

sub verbose
{
	my $self = shift;
	return $self->getParams->{'verbose'};
}

sub xlinktype
{
	my $self = shift;
	return $self->{'xlinktype'};
}

sub getMSTAB
{
	my $self = shift;
	if ( $self->{'MSTAB'} )
	{
		return $self->{'MSTAB'};
	} else
	{
		return $self->getmatchObj->getMSTAB;
	}
}

sub matchcommon
{    #match single peaks and keep theoretical peaks
	my $self      = shift;
	my $matchhash = shift;
	my $errorhash = shift;
	my $PARAMS    = $self->getParams;
	my $spectrum  = $self->getSpecObj;
	my ( $ion, $peak );
	my $tolerance      = $PARAMS->{'ms2tolerance'};
	my $peaks          = $spectrum->getcommonpeaks;
	my $candidate_ions = $self->getCommonIons;
	my ( $matchingions, $matchingpeaks );
	( $matchingions, $matchingpeaks ) = Match::matchpeaks( $candidate_ions, $peaks, $PARAMS, $tolerance, $matchhash, undef, $errorhash );
	my $xlinktype = $self->getxlinktype;

	if ( $xlinktype eq "xlink" )
	{
		my $dummy                = {};
		my $candidate_ions_alpha = $self->getCommonIons_alpha;
		( $self->{'commonmatches_alpha'}, $self->{'commonmatchingpeaks_alpha'} ) = Match::matchpeaks( $candidate_ions_alpha, $peaks, $PARAMS, $tolerance, $dummy );
		my $candidate_ions_beta = $self->getCommonIons_beta;
		$dummy = {};
		( $self->{'commonmatches_beta'}, $self->{'commonmatchingpeaks_beta'} ) = Match::matchpeaks( $candidate_ions_beta, $peaks, $PARAMS, $tolerance, $dummy );
	}
	### Added by TW
	## Filter redundancy out of the peaklist (could happen if alpha and beta share ions)
	nonred($matchingpeaks);
	return $matchingions, $matchingpeaks;
}

sub matchxlinks
{
	my $self      = shift;
	my $matchhash = shift;
	my $errorhash = shift;
	my $PARAMS    = $self->getParams;
	my $specobj   = $self->getSpecObj;
	my ( $tolerance, $ion, $charge, $peak, @matchions, @matchpeaks, $match2ndisotope );
	if ( $PARAMS->{'xlink_ms2tolerance'} )
	{
		$tolerance = $PARAMS->{'xlink_ms2tolerance'};
	} else
	{
		$tolerance = $PARAMS->{'ms2tolerance'};
	}
	if ( $PARAMS->{'CID_match2ndisotope'} )
	{
		$match2ndisotope = 1;
	} else
	{
		$match2ndisotope = 0;
	}
	if ( $PARAMS->{'tolerancemeasure_ms2'} =~ /^ppm/i )
	{

		# CHANGED	$tolerance = $tolerance * 1e-6 * 1000;    #ppm to amu measure
	}
	my $maxcharge = $self->maxioncharge_xlinks;
	my $mincharge = $self->minioncharge_xlinks;
	for $charge ( $mincharge .. $maxcharge )
	{
		my @peaks = ( @{ $self->getxlinkpeaks($charge) }, @{ $self->getxlinkpeaks("0") } );

		#get defined chargestate and undefined peaks
		my @candidate_ions = ( @{ $self->getXlinkIons->{$charge} } );
		my ( $matchingions, $matchingpeaks );
		( $matchingions, $matchingpeaks ) = Match::matchpeaks( \@candidate_ions, \@peaks, $PARAMS, $tolerance, $matchhash, $match2ndisotope / $charge, $errorhash );
		push @matchions,  @$matchingions;
		push @matchpeaks, @$matchingpeaks;
	}
	my $xlinktype = $self->getxlinktype;
	### This is the separate matching for alpha and beta
	if ( $xlinktype eq "xlink" )
	{
		my ( @matchions_alpha, @matchions_beta, @matchpeaks_alpha, @matchpeaks_beta );
		for $charge ( $mincharge .. $maxcharge )
		{
			my @peaks = ( @{ $self->getxlinkpeaks($charge) }, @{ $self->getxlinkpeaks("0") } );

			#get defined chargestate and undefined peaks
			my @candidate_ions_alpha = ( @{ $self->getXlinkIons_alpha->{$charge} } );
			my ( $matchingions_alpha, $matchingpeaks_alpha );
			my $dummy = {};
			( $matchingions_alpha, $matchingpeaks_alpha ) = Match::matchpeaks( \@candidate_ions_alpha, \@peaks, $PARAMS, $tolerance, $dummy, $match2ndisotope / $charge );
			push @matchions_alpha,  @$matchingions_alpha;
			push @matchpeaks_alpha, @$matchingpeaks_alpha;
			my @candidate_ions_beta = ( @{ $self->getXlinkIons_beta->{$charge} } );
			my ( $matchingions_beta, $matchingpeaks_beta );
			$dummy = {};
			( $matchingions_beta, $matchingpeaks_beta ) = Match::matchpeaks( \@candidate_ions_beta, \@peaks, $PARAMS, $tolerance, $dummy, $match2ndisotope / $charge );
			push @matchions_beta,  @$matchingions_beta;
			push @matchpeaks_beta, @$matchingpeaks_beta;
		}
		$self->{'xlinkmatches_alpha'}       = \@matchions_alpha;
		$self->{'xlinkmatchingpeaks_alpha'} = \@matchpeaks_alpha;
		$self->{'xlinkmatches_beta'}        = \@matchions_beta;
		$self->{'xlinkmatchingpeaks_beta'}  = \@matchpeaks_beta;
	}
	### Added by TW
	## Filter redundancy out of the peaklist (could happen if alpha and beta share ions)
	nonred( \@matchpeaks );
	return \@matchions, \@matchpeaks;
}

sub getxlinkpeaks
{
	my $self     = shift;
	my $charge   = shift;
	my $PARAMS   = $self->getParams;
	my $spectrum = $self->getSpecObj;
	if ( defined($charge) )
	{
		return $spectrum->getxlinkpeaks($charge);
	} else
	{
		my $mincharge = $self->minioncharge;
		my $maxcharge = $self->maxioncharge;
		my @tmp       = ();
		my $charge;
		for $charge ( $mincharge .. $maxcharge )
		{
			push @tmp, @{ $spectrum->getxlinkpeaks($charge) };
		}
		push @tmp, @{ $spectrum->getxlinkpeaks(0) };    #undefined charge peaks
		return \@tmp;
	}
}

sub getxlinkpairs
{
	my $self     = shift;
	my $charge   = shift;
	my $PARAMS   = $self->getParams;
	my $spectrum = $self->getSpecObj;
	if ( defined($charge) )
	{
		return $spectrum->getxlinkpairs($charge);
	} else
	{
		my $mincharge = $self->minioncharge_xlinks;
		my $maxcharge = $self->maxioncharge_xlinks;
		my @tmp       = ();
		my $charge;
		for $charge ( $mincharge .. $maxcharge )
		{
			push @tmp, @{ $spectrum->getxlinkpairs($charge) };
		}
		push @tmp, @{ $spectrum->getxlinkpairs(0) };    #undefined charge peaks
		return \@tmp;
	}
}

sub getXlinkIons
{
	my $self = shift;
	return $self->{'xlinkions'};
}

sub getXlinkIons_alpha
{
	my $self = shift;
	return $self->{'xlinkions_alpha'};
}

sub getXlinkIons_beta
{
	my $self = shift;
	return $self->{'xlinkions_beta'};
}

sub getCommonIons
{
	my $self = shift;
	return $self->{'commonions'};
}

sub getCommonIons_alpha
{
	my $self = shift;
	return $self->{'commonions_alpha'};
}

sub getCommonIons_beta
{
	my $self = shift;
	return $self->{'commonions_beta'};
}

sub getParams
{
	my $self = shift;
	return $self->{'params'};
}

sub getms1mass
{
	my $self = shift;
	return $self->{'matchobj'}->getms1mass;
}

sub get_flat_ionmatchhash
{
	my $self = shift;
	return $self->{'flat_ionmatchhash'};
}

sub get_flat_ionmatchhash_alpha
{
	my $self = shift;
	return $self->{'flat_ionmatchhash_alpha'};
}

sub get_flat_ionmatchhash_beta
{
	my $self = shift;
	return $self->{'flat_ionmatchhash_beta'};
}

sub get_flat_iontable
{
	my $self = shift;
	return $self->{'flationtable'};
}

sub getprecursorMz
{
	my $self = shift;
	return $self->getSpecObj->getprecursorMz;
}

sub getprecursorCharge
{
	my $self = shift;
	return $self->getSpecObj->getprecursorCharge;
}

sub get_apriori_pcommon
{
	my $self = shift;
	return $self->getSpecObj->get_apriori_pcommon;
}

sub get_apriori_pxlink
{
	my $self = shift;
	return $self->getSpecObj->get_apriori_pxlink;
}

sub getspectrumname
{
	my $self = shift;
	return $self->getSpecObj->getSpecname;
}

sub getxlinkermass
{
	my $self = shift;
	return $self->{'xlinkermw'};
}

sub gettotalxlinkermass
{
	my $self = shift;
	if ( $self->{'addedmass'} )
	{
		return $self->{'xlinkermw'} + $self->{'addedmass'};
	} elsif ( $self->getParams->{'addedmass'} )
	{
		return $self->{'xlinkermw'} + $self->getParams->{'addedmass'};
	} else
	{
		$self->{'xlinkermw'};
	}
}

sub get_ion_annotationtable
{
	my $self = shift;
	return $self->{'ion_annotationtable'};
}

sub getMatchingIonsIDs_for_pepstructure
{
	my $self             = shift;
	my $commonmatches    = $self->getcommonmatches;
	my $xlinkmatches     = $self->getxlinkmatches;
	my $iontable         = $self->getiontable;
	my $verbose          = $self->verbose;
	my $flationmatchhash = $self->get_flat_ionmatchhash;
	my %hittable;
	my %transiontable = (
						  'a'     => 'fwd',
						  'b'     => 'fwd',
						  'c'     => 'fwd',
						  'x'     => 'rev',
						  'y'     => 'rev',
						  'z'     => 'rev',
						  'beta'  => 'beta',
						  'alpha' => 'alpha',
	);

	foreach my $iontype ( sort ( keys %{$iontable} ) )
	{
		my @iontags = split /_/, $iontype;
		$verbose && print "$iontype<br>";
		foreach my $ionposition ( keys %{ $iontable->{$iontype} } )
		{
			my $mz = $iontable->{$iontype}->{$ionposition};
			if ( $flationmatchhash->{$mz} )
			{
				$hittable{ $transiontable{ $iontags[0] } }->{ $transiontable{ $iontags[2] } }->{$ionposition} = defined;
			}
		}
	}
	return \%hittable;
}

sub printhtmliontable_old
{
	my $self    = shift;
	my $hithash = shift;
	unless ($hithash)
	{
		$hithash = $self->getMatchingIonpositions;
	}
	my $iontable = $self->getiontable;
	my $ionnames = $self->getionnamehash;
	my $PARAMS   = $self->getParams;
	my %seqhash;
	if ( $PARAMS->{'variable_mod'} )
	{
	#	my ( $AA, $delta ) = split /,|:/, $PARAMS->{'variable_mod'};
	#	my $roundeddelta = sprintf( "%.0f", $delta );
	#	$seqhash{'X'} = join "", $AA, "[", int($roundeddelta), "]";
		my @AAdelta = split /,|:/, $PARAMS->{'variable_mod'};
		my @AAmods = ('X', 'U', 'B', 'J');
		for my $i (0..($#AAdelta / 2 - 1) ){
			my $roundeddelta = sprintf( "%.0f", $AAdelta[2*$i + 1]);
			$seqhash{$AAmods[$i]} = join "", $AAdelta[2*$i], "[", int($roundeddelta), "]";
		}
	}
	foreach my $iontype ( keys %$hithash )
	{
		my @positions = keys %{ $hithash->{$iontype} };

		#print "$iontype @positions <br>\n";
	}
	my ( @bionsalpha, @yionsalpha, @bionsbeta, @yionsbeta );
	my $fwdiontypes = $PARAMS->{'fwd_ions'};
	my $reviontypes = $PARAMS->{'rev_ions'};
	foreach ( sort keys %$iontable )
	{
		/alpha_(common|xlink)_($fwdiontypes)/ && ( push @bionsalpha, $_ );
		/alpha_(xlink|common)_($reviontypes)/ && ( push @yionsalpha, $_ );
		/beta_(xlink|common)_($fwdiontypes)/  && ( push @bionsbeta,  $_ );
		/beta_(xlink|common)_($reviontypes)/  && ( push @yionsbeta,  $_ );
	}
	my $seq1             = $self->gethitObjs->[0]->seq;
	my $seq2             = $self->gethitObjs->[1]->seq;
	my @seq1             = split //, $seq1;
	my @seq2             = split //, $seq2;
	my $xlinkerposition1 = $self->gettopology->[0];
	my $xlinkerposition2 = $self->gettopology->[1];
	my $i                = 0;
	print '<h2>alpha-chain</h2>';
	print ' <table border="1">';

	foreach my $iontype (@bionsalpha)
	{
		print '<tr>';
		if ( defined( $ionnames->{$iontype} ) )
		{
			print '<th class="smallfont">', $ionnames->{$iontype}, '</th>';
		} else
		{
			print '<th class="smallfont">', $iontype, '</th>';
		}
		for my $position ( 1 .. $#seq1 + 1 )
		{
			if ( $iontable->{$iontype}->{$position} )
			{
				if ( $hithash->{$iontype}->{$position} )
				{
					if ( $iontype =~ /H2O|NH3/ )
					{
						print '<th class=blue>', sprintf( "%.2f", $iontable->{$iontype}->{$position} ), '</th>';
					} elsif ( $iontype =~ /xlink/ )
					{
						print '<th class=red>', sprintf( "%.2f", $iontable->{$iontype}->{$position} ), '</th>';
					} else
					{
						print '<th class=green>', sprintf( "%.2f", $iontable->{$iontype}->{$position} ), '</th>';
					}
				} else
				{
					print '<th>', sprintf( "%.2f", $iontable->{$iontype}->{$position} ), '</th>';
				}
			} else
			{
				print '<th class="smallfont">-</th>';
			}
		}
		print '</tr>';
	}
	print '<tr>';
	print '<div class="smallfont">';
	print '<th >AA</th>';
	foreach my $AA (@seq1)
	{
		if ( $seqhash{$AA} )
		{
			$AA = $seqhash{$AA};
		}
		$i++;
		if ( $i == $xlinkerposition1 )
		{
			print '<th class="orange">', $AA, '</th>';
		} else
		{
			print '<th >', $AA, '</th>';
		}
	}
	print '</div>';
	print '</tr>';
	foreach my $iontype (@yionsalpha)
	{
		print '<tr>';
		if ( defined( $ionnames->{$iontype} ) )
		{
			print '<th class="smallfont">', $ionnames->{$iontype}, '</th>';
		} else
		{
			print '<th class="smallfont">', $iontype, '</th>';
		}
		for my $position ( 1 .. $#seq1 + 1 )
		{
			if ( $iontable->{$iontype}->{ $#seq1 - $position + 2 } )
			{
				if ( $hithash->{$iontype}->{ $#seq1 - $position + 2 } )
				{
					if ( $iontype =~ /H2O|NH3/ )
					{
						print '<th class=blue>', sprintf( "%.2f", $iontable->{$iontype}->{ $#seq1 - $position + 2 } ), '</th>';
					} elsif ( $iontype =~ /xlink/ )
					{
						print '<th class=red>', sprintf( "%.2f", $iontable->{$iontype}->{ $#seq1 - $position + 2 } ), '</th>';
					} else
					{
						print '<th class=green>', sprintf( "%.2f", $iontable->{$iontype}->{ $#seq1 - $position + 2 } ), '</th>';
					}
				} else
				{
					print '<th>', sprintf( "%.2f", $iontable->{$iontype}->{ $#seq1 - $position + 2 } ), '</th>';
				}
			} else
			{
				print '<th class="smallfont">-</th>';
			}
		}
		print '</tr>';
	}
	print '</table>';
	if ( $self->xlinktype eq "xlink" )
	{
		my $i = 0;
		print '<hr> <h2>beta-chain</h2>';
		print '<table border="1">';
		foreach my $iontype (@bionsbeta)
		{
			print '<tr>';
			if ( defined( $ionnames->{$iontype} ) )
			{
				print '<th class="smallfont">', $ionnames->{$iontype}, '</th>';
			} else
			{
				print '<th class="smallfont">', $iontype, '</th>';
			}
			for my $position ( 1 .. $#seq2 + 1 )
			{
				if ( $iontable->{$iontype}->{$position} )
				{
					if ( $hithash->{$iontype}->{$position} )
					{
						if ( $iontype =~ /H2O|NH3/ )
						{
							print '<th class=blue>', sprintf( "%.2f", $iontable->{$iontype}->{$position} ), '</th>';
						} elsif ( $iontype =~ /xlink/ )
						{
							print '<th class=red>', sprintf( "%.2f", $iontable->{$iontype}->{$position} ), '</th>';
						} else
						{
							print '<th class=green>', sprintf( "%.2f", $iontable->{$iontype}->{$position} ), '</th>';
						}
					} else
					{
						print '<th>', sprintf( "%.2f", $iontable->{$iontype}->{$position} ), '</th>';
					}
				} else
				{
					print '<th class="smallfont">-</th>';
				}
			}
			print '</tr>';
		}
		print '<tr>';
		print '<div class="smallfont">';
		print '<th >AA</th>';
		foreach my $AA (@seq2)
		{
			if ( $seqhash{$AA} )
			{
				$AA = $seqhash{$AA};
			}
			$i++;
			if ( $i == $xlinkerposition2 )
			{
				print '<th class="orange">', $AA, '</th>';
			} else
			{
				print '<th >', $AA, '</th>';
			}
		}
		print '</div>';
		print '</tr>';
		foreach my $iontype (@yionsbeta)
		{
			print '<tr>';
			if ( defined( $ionnames->{$iontype} ) )
			{
				print '<th class="smallfont">', $ionnames->{$iontype}, '</th>';
			} else
			{
				print '<th class="smallfont">', $iontype, '</th>';
			}
			for my $position ( 1 .. $#seq2 + 1 )
			{
				if ( $iontable->{$iontype}->{ $#seq2 - $position + 2 } )
				{
					if ( $hithash->{$iontype}->{ $#seq2 - $position + 2 } )
					{
						if ( $iontype =~ /NH3|H2O/ )
						{
							print '<th class=blue>', sprintf( "%.2f", $iontable->{$iontype}->{ $#seq2 - $position + 2 } ), '</th>';
						} elsif ( $iontype =~ /xlink/ )
						{
							print '<th class=red>', sprintf( "%.2f", $iontable->{$iontype}->{ $#seq2 - $position + 2 } ), '</th>';
						} else
						{
							print '<th class=green>', sprintf( "%.2f", $iontable->{$iontype}->{ $#seq2 - $position + 2 } ), '</th>';
						}
					} else
					{
						print '<th>', sprintf( "%.2f", $iontable->{$iontype}->{ $#seq2 - $position + 2 } ), '</th>';
					}
				} else
				{
					print '<th class="smallfont">-</th>';
				}
			}
			print '</tr>';
		}
		print '</table>';
	}
}

sub getminionsize
{
	my $self = shift;
	return $self->getParams->{'minionsize'};
}

sub getmaxionsize
{
	my $self = shift;
	return $self->getParams->{'maxionsize'};
}

sub printhtmliontable
{
	my $self       = shift;
	my $minionsize = shift;
	my $maxionsize = shift;
	my $iontable   = $self->getiontable;
	my $PARAMS     = $self->getParams;
	unless ($minionsize)
	{
		$minionsize = $PARAMS->{'minionsize'};
	}
	unless ($maxionsize)
	{
		$maxionsize = $PARAMS->{'maxionsize'};
	}
	my %seqhash;
	if ( $PARAMS->{'variable_mod'} )
	{
	#	my ( $AA, $delta ) = split /,|:/, $PARAMS->{'variable_mod'};
	#	my $roundeddelta = sprintf( "%.0f", $delta );
	#	$seqhash{'X'} = join "", $AA, "[", int($roundeddelta), "]";
		my @AAdelta = split /,|:/, $PARAMS->{'variable_mod'};
		my @AAmods = ('X', 'U', 'B', 'J');
		for my $i (0..($#AAdelta / 2 - 1) ){
			my $roundeddelta = sprintf( "%.0f", $AAdelta[2*$i + 1]);
			$seqhash{$AAmods[$i]} = join "", $AAdelta[2*$i], "[", int($roundeddelta), "]";
		}
	}
	my ( @bionsalpha, @yionsalpha, @bionsbeta, @yionsbeta );
	my $fwdiontypes = $PARAMS->{'fwd_ions'};
	my $reviontypes = $PARAMS->{'rev_ions'};
	foreach ( sort keys %$iontable )
	{
		/alpha_(common|xlink)_($fwdiontypes)/ && ( push @bionsalpha, $_ );
		/alpha_(xlink|common)_($reviontypes)/ && ( push @yionsalpha, $_ );
		/beta_(xlink|common)_($fwdiontypes)/  && ( push @bionsbeta,  $_ );
		/beta_(xlink|common)_($reviontypes)/  && ( push @yionsbeta,  $_ );
	}
	my $seq1             = $self->gethitObjs->[0]->seq;
	
	my $seq2;
	if ($self->xlinktype eq "xlink"){
	$seq2             = $self->gethitObjs->[1]->seq;
	}
	
	my @seq1             = split //, $seq1;
	my @seq2             = split //, $seq2;
	my $xlinkerposition1 = $self->gettopology->[0];
	my $xlinkerposition2 = $self->gettopology->[1];
	my $i                = 0;
	my $commonmatches    = $self->getcommonmatches;
	my $xlinkmatches     = $self->getxlinkmatches;
	my %matchhash        = ();

	foreach my $mz ( @$commonmatches, @$xlinkmatches )
	{
		$matchhash{$mz} = 1;
	}

	#print "min: $minionsize max: $maxionsize <br>";
	print '<h2>alpha-chain</h2>';
	print ' <table border="1">';
	foreach my $iontype (@bionsalpha)
	{
		print '<tr>';
		print '<th class="smallfont">', $iontype, '</th>';
		for my $position ( 1 .. $#seq1 + 1 )
		{
			if ( $iontable->{$iontype}->{$position} )
			{
				my $mz = $iontable->{$iontype}->{$position};
				if ( $matchhash{$mz} )
				{
					if ( $iontype =~ /H2O|NH3/ )
					{
						print '<th class=blue>', sprintf( "%.2f", $mz ), '</th>';
					} elsif ( $iontype =~ /xlink/ )
					{
						print '<th class=red>', sprintf( "%.2f", $mz ), '</th>';
					} else
					{
						print '<th class=green>', sprintf( "%.2f", $mz ), '</th>';
					}
				} elsif ( ( $mz < $minionsize ) || ( $mz > $maxionsize ) )
				{
					print '<th class=low>', sprintf( "%.2f", $mz ), '</th>';
				} else
				{
					print '<th class=normal>', sprintf( "%.2f", $mz ), '</th>';
				}
			} else
			{
				print '<th class=smallfont>-</th>';
			}
		}
		print '</tr>';
	}
	print '<tr>';
	print '<div class="smallfont">';
	print '<th >AA</th>';
	foreach my $AA (@seq1)
	{
		if ( $seqhash{$AA} )
		{
			$AA = $seqhash{$AA};
		}
		$i++;
		if (
			 $i == $xlinkerposition1
			 || (    ( $self->xlinktype eq "intralink" )
				  && ( $i == $xlinkerposition2 ) )
		  )
		{
			print '<th class="orange">', $AA, '</th>';
		} else
		{
			print '<th >', $AA, '</th>';
		}
	}
	print '</div>';
	print '</tr>';
	foreach my $iontype (@yionsalpha)
	{
		print '<tr>';
		print '<th class="smallfont">', $iontype, '</th>';
		for my $position ( 1 .. $#seq1 + 1 )
		{
			if ( $iontable->{$iontype}->{ $#seq1 - $position + 2 } )
			{
				my $mz = $iontable->{$iontype}->{ $#seq1 - $position + 2 };
				if ( $matchhash{$mz} )
				{
					if ( $iontype =~ /H2O|NH3/ )
					{
						print '<th class=blue>', sprintf( "%.2f", $mz ), '</th>';
					} elsif ( $iontype =~ /xlink/ )
					{
						print '<th class=red>', sprintf( "%.2f", $mz ), '</th>';
					} else
					{
						print '<th class=green>', sprintf( "%.2f", $mz ), '</th>';
					}
				} elsif ( ( $mz < $minionsize ) || ( $mz > $maxionsize ) )
				{
					print '<th class=low>', sprintf( "%.2f", $mz ), '</th>';
				} else
				{
					print '<th>', sprintf( "%.2f", $mz ), '</th>';
				}
			} else
			{
				print '<th class="smallfont">-</th>';
			}
		}
		print '</tr>';
	}
	print '</table>';
	if ( $self->xlinktype eq "xlink" )
	{
		my $i = 0;
		print '<hr> <h2>beta-chain</h2>';
		print '<table border="1">';
		foreach my $iontype (@bionsbeta)
		{
			print '<tr>';
			print '<th class="smallfont">', $iontype, '</th>';
			for my $position ( 1 .. $#seq2 + 1 )
			{
				if ( $iontable->{$iontype}->{$position} )
				{
					my $mz = $iontable->{$iontype}->{$position};
					if ( $matchhash{$mz} )
					{
						if ( $iontype =~ /H2O|NH3/ )
						{
							print '<th class=blue>', sprintf( "%.2f", $mz ), '</th>';
						} elsif ( $iontype =~ /xlink/ )
						{
							print '<th class=red>', sprintf( "%.2f", $mz ), '</th>';
						} else
						{
							print '<th class=green>', sprintf( "%.2f", $mz ), '</th>';
						}
					} elsif ( ( $mz < $minionsize ) || ( $mz > $maxionsize ) )
					{
						print '<th class=low>', sprintf( "%.2f", $mz ), '</th>';
					} else
					{
						print '<th>', sprintf( "%.2f", $mz ), '</th>';
					}
				} else
				{
					print '<th class="smallfont">-</th>';
				}
			}
			print '</tr>';
		}
		print '<tr>';
		print '<div class="smallfont">';
		print '<th >AA</th>';
		foreach my $AA (@seq2)
		{
			if ( $seqhash{$AA} )
			{
				$AA = $seqhash{$AA};
			}
			$i++;
			if ( $i == $xlinkerposition2 )
			{
				print '<th class="orange">', $AA, '</th>';
			} else
			{
				print '<th >', $AA, '</th>';
			}
		}
		print '</div>';
		print '</tr>';
		foreach my $iontype (@yionsbeta)
		{
			print '<tr>';
			print '<th class="smallfont">', $iontype, '</th>';
			for my $position ( 1 .. $#seq2 + 1 )
			{
				if ( $iontable->{$iontype}->{ $#seq2 - $position + 2 } )
				{
					my $mz = $iontable->{$iontype}->{ $#seq2 - $position + 2 };
					if ( $matchhash{$mz} )
					{
						if ( $iontype =~ /NH3|H2O/ )
						{
							print '<th class=blue>', sprintf( "%.2f", $mz ), '</th>';
						} elsif ( $iontype =~ /xlink/ )
						{
							print '<th class=red>', sprintf( "%.2f", $mz ), '</th>';
						} else
						{
							print '<th class=green>', sprintf( "%.2f", $mz ), '</th>';
						}
					} elsif ( ( $mz < $minionsize ) || ( $mz > $maxionsize ) )
					{
						print '<th class=low>', sprintf( "%.2f", $mz ), '</th>';
					} else
					{
						print '<th>', sprintf( "%.2f", $mz ), '</th>';
					}
				} else
				{
					print '<th class="smallfont">-</th>';
				}
			}
			print '</tr>';
		}
		print '</table>';
	}
}

sub calcMatchingIonpositions
{    #get matchin ions
	my $self      = shift;
	my $matchhash = shift;

	#	my $commonmatches = $self->getcommonmatches;
	#	my $xlinkmatches  = $self->getxlinkmatches;
	unless ($matchhash)
	{
		$matchhash = $self->get_flat_ionmatchhash;
	}
	my $iontable = $self->getiontable;
	my %reversetable;
	foreach my $iontype ( keys %{$iontable} )
	{
		foreach my $ionposition ( keys %{ $iontable->{$iontype} } )
		{
			my $mz = $iontable->{$iontype}->{$ionposition};
			if ( $matchhash->{$mz} )
			{
				$reversetable{$iontype}->{$ionposition} = 1;
			}
		}
	}
	return \%reversetable;
}

sub make_flationmatchhash
{
	my $self          = shift;
	my %matchhash     = ();
	my $commonmatches = $self->getcommonmatches;
	my $xlinkmatches  = $self->getxlinkmatches;
	foreach my $mz ( @$commonmatches, @$xlinkmatches )
	{
		$matchhash{$mz} = 1;
	}

	#	my $xlinktype = $self->getxlinktype;
	#
	#	if ( $xlinktype eq "xlink" && $self->params->{'matchsinglepeps'} ) {
	#	$self->{'flat_ionmatchhash_alpha'}     = $self->make_flationmatchhash_alpha;
	#	$self->{'flat_ionmatchhash_beta'}     = $self->make_flationmatchhash_beta;
	#	}
	return \%matchhash;
}

sub printmatchingions
{
	my $self              = shift;
	my $matchingionsalpha = $self->getcommonmatches;
	print "commonmatches all\n";
	foreach my $ion (@$matchingionsalpha)
	{
		print $ion, "\n";
	}
}

sub printmatchingionsalpha
{
	my $self              = shift;
	my $matchhash         = $self->get_ion2peakmatchhash;
	my $matchingionsalpha = $self->getcommonmatches_alpha;
	print "commonmatches alpha<br>\n";
	foreach my $ion (@$matchingionsalpha)
	{
		print $ion, "\t", $matchhash->{$ion}, "\n<br>";
	}
}

sub printmatchingionsbeta
{
	my $self             = shift;
	my $matchingionsbeta = $self->getcommonmatches_beta;
	my $matchhash        = $self->get_ion2peakmatchhash;
	print "commonmatches beta<br>\n";
	foreach my $ion (@$matchingionsbeta)
	{
		print $ion, "\t", $matchhash->{$ion}, "\n<br>";
	}
}

sub make_flationmatchhash_alpha
{
	my $self          = shift;
	my %matchhash     = ();
	my $commonmatches = $self->getcommonmatches_alpha;
	my $xlinkmatches  = $self->getxlinkmatches_alpha;
	foreach my $mz ( @$commonmatches, @$xlinkmatches )
	{
		$matchhash{$mz} = 1;
	}
	return \%matchhash;
}

sub make_flationmatchhash_beta
{
	my $self          = shift;
	my %matchhash     = ();
	my $commonmatches = $self->getcommonmatches_beta;
	my $xlinkmatches  = $self->getxlinkmatches_beta;
	foreach my $mz ( @$commonmatches, @$xlinkmatches )
	{
		$matchhash{$mz} = 1;
	}
	return \%matchhash;
}

sub makehithash_for_xml
{
	my $self   = shift;
	my $ionmzs = $self->getiontable;
	my $table  = $self->getMatchingIonpositions;
	my %hithash;
	foreach my $iontype ( sort keys %$ionmzs )
	{
		$hithash{$iontype} = [];
		foreach my $position ( sort keys %{ $table->{$iontype} } )
		{
			if ( defined( $table->{$iontype}->{$position} ) )
			{
				push @{ $hithash{$iontype} }, $position;
			}
		}
	}
	return \%hithash;
}

sub printhit_attributes
{
	my $self          = shift;
	my $outfilehandle = shift;
	my $hithash       = $self->makehithash_for_xml;
	foreach my $iontype ( sort keys %$hithash )
	{
		my $hitstring = join ",", @{ $hithash->{$iontype} };
		print $outfilehandle " ", $iontype, '="', $hitstring, '"';
	}

	#	print $outfilehandle " >";
}

sub printhittable
{
	my $self                 = shift;
	my $matchingionpositions = shift;
	unless ($matchingionpositions)
	{
		$matchingionpositions = $self->getMatchingIonpositions;
	}
	my $iontable          = $self->getiontable;
	my $ion2peakmatchhash = $self->get_ion2peakmatchhash;    ## ion is key, peaks are arrayref
	my $specObj           = $self->getSpecObj;

	#print Dumper ($ion2peakmatchhash);
	my %colorcode = (
					  "common"     => "green",
					  "xlink"      => "red",
					  "loss"       => "blue",
					  "2ndisotope" => "orange",
	);
	print "<br>";
	print '<table border="2">';
	print '<tr>
       <th>type</th>
       <th>position</th>
       <th>ion th</th>
       <th>peak</th>
       <th>delta mz</th>  
       <th>delta ppm</th>  
       <th>intensity</th>
     </tr>';
	foreach my $iontype ( sort keys %$iontable )
	{

		foreach my $position ( sort keys %{ $iontable->{$iontype} } )
		{
			if ( defined( $matchingionpositions->{$iontype}->{$position} ) )
			{
				print '<tr>';
				if ( $iontype =~ /H2O|NH3/ )
				{
					print '<th class=blue>', $iontype, '</th>';
				} elsif ( $iontype =~ /xlink/ )
				{
					print '<th class=red>', $iontype, '</th>';
				} else
				{
					print '<th class=green>', $iontype, '</th>';
				}
				my $ion_mz = $iontable->{$iontype}->{$position};
				my @matchedpeakarray;
				my @deltaarray;
				my @peakintensityarray;
				my $peakanno;
				my @ppmarray;
				foreach my $peak ( @{ $ion2peakmatchhash->{$ion_mz} } )
				{
					### prepare the arrays
					push @matchedpeakarray, sprintf( "%.3f", $peak );
					push @deltaarray,       sprintf( "%.3f", $ion_mz - $peak );
					my $deltamz = $ion_mz - $peak;
					my $deltappm = 1e6 * ( $deltamz / $ion_mz );
					push @ppmarray,           sprintf( "%.0f", $deltappm );
					push @peakintensityarray, sprintf( "%.3f", $specObj->get_ionintensity($peak) );

					#$peakanno      = $self->annotate_matchedpeak($peak);
				}

				#my $matchedpeak   = sprintf( "%.3f", @{$ion2peakmatchhash->{$ion_mz}} );
				#my $delta         = $ion_mz - $matchedpeak;
				#my $peakintensity = $specObj->get_ionintensity($matchedpeak);
				#my $peakanno      = $self->annotate_matchedpeak($matchedpeak);
				print '<th>', $position, '</th>', '<th>', sprintf( "%.3f", $ion_mz ), '</th>', '<th>', join( ",", @matchedpeakarray ), '</th>', '<th>', join( ",", @deltaarray ), '</th>', '<th>', join( ",", @ppmarray ), '</th>', '<th>', join( ",", @peakintensityarray ), '</th>', '<th>',

				  #sprintf( "%.3f", $matchedpeak ), '</th>', '<th>',
				  #($matchedpeak ), '</th>', '<th>',
				  #sprintf( "%.3f", $delta ),       '</th>',
				  #'<th>', sprintf( "%.3f", $peakintensity ),
				  #'</th>';    #, " $peakanno ";
			}
			print '</tr>';
		}
	}
	print '</table>';
	print "Please note: if delta values are larger than the selected tolerance then this due to matching of the second isotope peak.<br>";
}

sub generate_subspec
{
	my $self                 = shift;
	my $subtractpeptide      = shift;
	my $matchingionpositions = $self->getMatchingIonpositions;
	my $iontable             = $self->getiontable;
	my $ion2peakmatchhash    = $self->get_ion2peakmatchhash;     ## ion is key, peaks are arrayref
	my $specObj              = $self->getSpecObj;
	my $matchedpeaksinspec;
    my $Hatom=1.007825032;
	#print Dumper ($ion2peakmatchhash);
	my %colorcode = (
					  "common"     => "green",
					  "xlink"      => "red",
					  "loss"       => "blue",
					  "2ndisotope" => "orange",
	);
	print "<br>";
	print '<table border="2">';
	print '<tr>
       <th>type</th>
       <th>position</th>
       <th>ion th</th>
       <th>peak</th>
       <th>delta mz</th>  
       <th>delta ppm</th>  
       <th>intensity</th>
     </tr>';
	foreach my $iontype ( sort keys %$iontable )
	{

		foreach my $position ( sort keys %{ $iontable->{$iontype} } )
		{
			if ( defined( $matchingionpositions->{$iontype}->{$position} ) )
			{
				print '<tr>';
				if ( $iontype =~ /$subtractpeptide/ )
				{

					#print "True\n";
				} else
				{
					next;
				}
				if ( $iontype =~ /H2O|NH3/ )
				{
					print '<th class=blue>', $iontype, '</th>';
				} elsif ( $iontype =~ /xlink/ )
				{
					print '<th class=red>', $iontype, '</th>';
				} else
				{
					print '<th class=green>', $iontype, '</th>';
				}
				my $ion_mz = $iontable->{$iontype}->{$position};
				my @matchedpeakarray;
				my @deltaarray;
				my @peakintensityarray;
				my $peakanno;
				my @ppmarray;
				foreach my $peak ( @{ $ion2peakmatchhash->{$ion_mz} } )
				{
					### prepare the arrays
					push @matchedpeakarray, sprintf( "%.3f", $peak );
					$matchedpeaksinspec->{$peak} = 1;
					push @deltaarray, sprintf( "%.3f", $ion_mz - $peak );
					my $deltamz = $ion_mz - $peak;
					my $deltappm = 1e6 * ( $deltamz / $ion_mz );
					push @ppmarray,           sprintf( "%.0f", $deltappm );
					push @peakintensityarray, sprintf( "%.3f", $specObj->get_ionintensity($peak) );

					#$peakanno      = $self->annotate_matchedpeak($peak);
				}

				#my $matchedpeak   = sprintf( "%.3f", @{$ion2peakmatchhash->{$ion_mz}} );
				#my $delta         = $ion_mz - $matchedpeak;
				#my $peakintensity = $specObj->get_ionintensity($matchedpeak);
				#my $peakanno      = $self->annotate_matchedpeak($matchedpeak);
				print '<th>', $position, '</th>', '<th>', sprintf( "%.3f", $ion_mz ), '</th>', '<th>', join( ",", @matchedpeakarray ), '</th>', '<th>', join( ",", @deltaarray ), '</th>', '<th>', join( ",", @ppmarray ), '</th>', '<th>', join( ",", @peakintensityarray ), '</th>', '<th>',

				  #sprintf( "%.3f", $matchedpeak ), '</th>', '<th>',
				  #($matchedpeak ), '</th>', '<th>',
				  #sprintf( "%.3f", $delta ),       '</th>',
				  #'<th>', sprintf( "%.3f", $peakintensity ),
				  #'</th>';    #, " $peakanno ";
			}
			print '</tr>';
		}
	}
	print '</table>';
	print "Please note: if delta values are larger than the selected tolerance then this due to matching of the second isotope peak.<br>";
	#print Dumper ($matchedpeaksinspec);


## Get peptides
my $hitobjects = $self->gethitObjs;
	
my $precursormass=$self->getprecursorMz;
my $precursorcharge= $self->getprecursorCharge;
my $xlinkermass= $self->getxlinkermass;
my $ms1mass=$self->getms1mass;
my $theoreticalmass= $self->calc_mr;
my $ppmerror = $self->calcerrorppm;

print "Mr: $precursormass<br>";
print "Charge: $precursorcharge<br>";
print "XL-Mass: $xlinkermass<br>";
print "MS1 mass: $ms1mass<br>";
print "MS1 mass theoretical: $theoreticalmass<br>";
print "Error ppm: $ppmerror<br>";
#$hit->molweight;
	
my $hitobjecttosubstract;
my $targethitobj;
	if ( $subtractpeptide =~ /alpha/ )
	{
		$hitobjecttosubstract = $hitobjects->[0];
		$targethitobj = $hitobjects->[1];
	} else
	{
		$hitobjecttosubstract = $hitobjects->[1];
		$targethitobj = $hitobjects->[0];
	}

my $pepobject1mass= $hitobjects->[0]->molweight;
my $pepobject2mass= $hitobjects->[1]->molweight;
print "Pepobject1mass: $pepobject1mass<br>";
print "Pepobject2mass: $pepobject2mass<br>";
my $masstosubtract= $hitobjecttosubstract->molweight + $xlinkermass;
print "Subtraction mass: $masstosubtract<br>";
my $z1target=$targethitobj->molweight+$Hatom;
$targethitobj->printtable;
#exit;


print "Z1mass of target: $z1target<br>";


#print Dumper($hitobjecttosubstract);
## Get peaks from specobject
my $commonions= $self->getcommonpairs;
my $xlinkions = $self->getxlinkpairs;
print Dumper($commonions);
#exit;
my $subtractedpeaks={};


foreach my $commonpeak (@$commonions){
print $commonpeak->[0];
my $mz=$commonpeak->[0];
my $int = $commonpeak->[1];
if ($matchedpeaksinspec->{$commonpeak->[0]}){
print "-->seen ";	
}else{
### Check mass and filter out impossible masses
if (checkmass($mz, 200, $targethitobj->molweight+$Hatom)){
$subtractedpeaks->{$mz}=$int;
}
}
print "<br>";	
}

foreach my $xlinkpeak (@$xlinkions){
print $xlinkpeak->[0];
my $mz=$xlinkpeak->[0];
my $int=$xlinkpeak->[1];
my $z= $xlinkpeak->[2];
my $nominalmass = $mz * $z - $z * $Hatom;
my $subtractedmass= $nominalmass-$masstosubtract;
my $z1mass=$subtractedmass+$Hatom;
print " z: $z nomianlmass: $nominalmass  converted mass: $subtractedmass ";

if ($matchedpeaksinspec->{$xlinkpeak->[0]}){
print "-->seen ";	
}else{
### Check mass and filter out impossible masses
if (checkmass($z1mass, 200, $targethitobj->molweight+$Hatom)){
$subtractedpeaks->{$z1mass}=$int;}	
}
print "<br>";	
}

## Print spectrum
print Dumper ($subtractedpeaks);
print "<br>";
print "BEGIN IONS<br>";
print "TITLE=testspectrum.dta<br>";
print "PEPMASS=$z1target<br>";
print "CHARGE=1+<br>";
foreach my $value( sort {$a <=> $b} keys %$subtractedpeaks){
print "$value $subtractedpeaks->{$value}<br>";
}
print "END IONS<br>";











}

sub checkmass{
my $mass=shift;
my $lowerlimit=shift;
my $uperlimit=shift;

if ($mass > $lowerlimit && $mass < $uperlimit){

#print "checking mass: $mass .. ok<br>";
return 1;
}else{
return 0;
}

}



sub get_peak2ionmatchhash
{
	my $self = shift;
	return $self->{'peak2ionmatchhash'};
}

sub get_peak2ionmatch
{
	my $self = shift;
	my $peak = shift;
	return $self->{'peak2ionmatchhash'}->{$peak};
}

sub annotate_matchedpeak
{
	my $self              = shift;
	my $peakmz            = shift;
	my $ionannotation     = $self->get_ion_annotationtable;
	my $peak2ionmatchhash = $self->get_peak2ionmatchhash;

	#	print Dumper ($peak2ionmatchhash);
	#	print "PEAK: $peakmz, ION: $peak2ionmatchhash->{$peakmz}<br>";
	return $ionannotation->{ $peak2ionmatchhash->{$peakmz} };
}

sub annotate_matchedion
{
	my $self  = shift;
	my $ionmz = shift;
	return $self->get_ion_annotationtable->{$ionmz};
}

sub getionintensity
{
	my $self = shift;
	my $mz   = shift;
	return $self->getSpecObj->get_ionintensity->{$mz};
}

sub gettopology
{
	my $self = shift;
	return $self->{'topology'};
}

sub check_topology{
my $self=shift;
my $pepobjs          = $self->gethitObjs;
my $PARAMS           = $self->getParams;
my $xlinktype = $self->getxlinktype;
my $toplogy=$self->gettopology;
my $possibletopologies=$PARAMS->{'possibleTopology'};

my ($AA1, $AA2);
## if xlinktype is a monolink return 1, nothing to check
if ( $xlinktype eq "monolink" ){
return 1;
} elsif ( $xlinktype eq "xlink" )
{
$AA1=$self->get_AA($pepobjs->[0],$toplogy->[0]);
$AA2=$self->get_AA($pepobjs->[1],$toplogy->[1]);
#print "Xlink: Linked AAs are: $AA1 : $AA2\n";	
} elsif ( $xlinktype eq "intralink" )
{
$AA1=$self->get_AA($pepobjs->[0],$toplogy->[0]);
$AA2=$self->get_AA($pepobjs->[0],$toplogy->[1]);
#print "Looplink: Linked AAs are: $AA1 : $AA2\n";	
}

unless ($AA1 && $AA2 && $possibletopologies){
die "Variable not det in check_topology, in LinkObj.pm\n";
}

my $topostring=$AA1.":".$AA2;
#print "Checking Topology: $topostring\n";
## check if Topology is valid
if ($possibletopologies->{$topostring}){
#print "Topology possibel\n";	
return 1;
}else{
return 0;
#print "Topology not possibel\n";	
}


	
	
}

sub get_AA{
my $self=shift;
my $pepobject=shift;
my $topo=shift;
my $AA = substr($pepobject->seq, $topo-1, 1);
return $AA;
	
}

sub calcerrorppm
{
	my $self         = shift;
	my $hits         = $self->gethitObjs;
	my $ms1parent    = $self->getms1mass;
	my $xlinkermass  = $self->getxlinkermass;
	my $xlinkpepmass = 0;
	foreach my $hit (@$hits)
	{
		$xlinkpepmass += $hit->molweight;
	}
	$xlinkpepmass += $xlinkermass;
	return 1e6 * abs( $ms1parent - $xlinkpepmass ) / $ms1parent;
}

sub calcerrorppmrelativ
{
	my $self         = shift;
	my $hits         = $self->gethitObjs;
	my $ms1parent    = $self->getms1mass;
	my $xlinkermass  = $self->getxlinkermass;
	my $xlinkpepmass = 0;
	foreach my $hit (@$hits)
	{
		$xlinkpepmass += $hit->molweight;
	}
	$xlinkpepmass += $xlinkermass;
	return 1e6 * ( $ms1parent - $xlinkpepmass ) / $ms1parent;
}

sub calc_specific_ions{
my $self=shift;
my $hitobjects         = $self->gethitObjs;
my $pepobject1mass= $hitobjects->[0]->molweight;
my $pepobject2mass= $hitobjects->[1]->molweight;
my $xlinkermass= $self->getxlinkermass;
my $theoreticalmass= $self->calc_mr;

my $Hatom=1.007825032;
my $nh2mass=17.02655;
my $specificion_b_pep1=$theoreticalmass-$pepobject1mass-$nh2mass;
my $specificion_b_pep2=$pepobject2mass-$nh2mass;
my $specificion_b_pep1_z1=$specificion_b_pep1+$Hatom;
my $specificion_b_pep2_z1=$specificion_b_pep2+$Hatom;
my $specificion_b_pep1_z2=($specificion_b_pep1+2*$Hatom)/2;
my $specificion_b_pep2_z2=($specificion_b_pep2+2*$Hatom)/2;
my $specificion_b_pep1_z3=($specificion_b_pep1+3*$Hatom)/3;
my $specificion_b_pep2_z3=($specificion_b_pep2+3*$Hatom)/3;
my $specificion_b_pep1_z4=($specificion_b_pep1+4*$Hatom)/4;
my $specificion_b_pep2_z4=($specificion_b_pep2+4*$Hatom)/4;

my @ions;
push @ions,  $specificion_b_pep1_z1;
push @ions,  $specificion_b_pep2_z1;
push @ions,  $specificion_b_pep1_z2;
push @ions,  $specificion_b_pep2_z2;
push @ions,  $specificion_b_pep1_z3;
push @ions,  $specificion_b_pep2_z3;
push @ions,  $specificion_b_pep1_z4;
push @ions,  $specificion_b_pep2_z4;
my $k=128.094963;
print "Linkermass:$xlinkermass<br>";
my $specificion_y_pep1=$pepobject1mass+$xlinkermass+$k+$Hatom;
my $specificion_y_pep2=$pepobject2mass+$xlinkermass+$k+$Hatom;
my $specificion_y_pep1_z1=$specificion_y_pep1+$Hatom;
my $specificion_y_pep2_z1=$specificion_y_pep2+$Hatom;
my $specificion_y_pep1_z2=($specificion_y_pep1+2*$Hatom)/2;
my $specificion_y_pep2_z2=($specificion_y_pep2+2*$Hatom)/2;
my $specificion_y_pep1_z3=($specificion_y_pep1+3*$Hatom)/3;
my $specificion_y_pep2_z3=($specificion_y_pep2+3*$Hatom)/3;
my $specificion_y_pep1_z4=($specificion_y_pep1+4*$Hatom)/4;
my $specificion_y_pep2_z4=($specificion_y_pep2+4*$Hatom)/4;
my @xlions;
push @xlions,  $specificion_y_pep1_z1;
push @xlions,  $specificion_y_pep1_z2;
push @xlions,  $specificion_y_pep1_z3;
push @xlions,  $specificion_y_pep1_z4;
push @xlions,  $specificion_y_pep2_z1;
push @xlions,  $specificion_y_pep2_z2;
push @xlions,  $specificion_y_pep2_z3;
push @xlions,  $specificion_y_pep2_z4;

print "Specific Mass b pep1: $specificion_b_pep1<br>";
print "Specific Mass b pep2: $specificion_b_pep2<br>";
print "Specific Mass y pep1: $specificion_y_pep1<br>";
print "Specific Mass y pep2: $specificion_y_pep2<br>";
# Match the ions
my $PARAMS    = $self->getParams;
my $spectrum  = $self->getSpecObj;
my ( $ion, $peak );
my $tolerance      = $PARAMS->{'ms2tolerance'};
my $peaks          = $spectrum->getcommonpeaks;
my ( $matchingions, $matchingpeaks );
my $errorhash;
my $matchhash;

my $candidate_ions         = \@ions;

( $matchingions, $matchingpeaks ) = Match::matchpeaks( $candidate_ions, $peaks, $PARAMS, $tolerance, $matchhash, undef, $errorhash );
#print Dumper($candidate_ions);
#print Dumper($peaks);
#print Dumper($matchingpeaks);
my @xlinkpeaks;
	my $maxcharge = 5;
	my $mincharge = 1;
	for my $charge ( $mincharge .. $maxcharge )
	{
	my  @tmppeaks=( @{ $self->getxlinkpeaks($charge) });		
	push @xlinkpeaks, @tmppeaks;
	}
#print Dumper(@xlinkpeaks);
( $matchingions, $matchingpeaks ) = Match::matchpeaks( \@xlions, \@xlinkpeaks, $PARAMS, $tolerance, $matchhash, undef, $errorhash );
#print Dumper($peaks);
print Dumper($matchingpeaks);

}

sub calc_mr
{
	my $self         = shift;
	my $hits         = $self->gethitObjs;
	my $xlinkermass  = $self->getxlinkermass;
	my $xlinkpepmass = 0;
	foreach my $hit (@$hits)
	{
		$xlinkpepmass += $hit->molweight;
	}
	$xlinkpepmass += $xlinkermass;
	return $xlinkpepmass;
}

sub _printxlinks_intralink
{
	my $self    = shift;
	my $outfile = shift;
	unless ( defined($outfile) )
	{
		$outfile = *STDOUT;
	}
	my $hit              = $self->gethitObjs->[0];
	my $intralinkmass    = $self->getxlinkermass;
	my $monolinkpepmass  = $hit->molweight + $intralinkmass;
	my $ms1parent        = $self->getms1mass;
	my $xlinkmatchratio  = join "\/", $self->get_number_of_Xlinkmatches, $self->get_number_of_XlinkIons;
	my $commonmatchratio = join "\/", $self->get_number_of_Commonmatches, $self->get_number_of_CommonIons;
	print $outfile $self->getspectrumname, ": ";
	print $outfile $self->getid,           " ";
	print $outfile $self->topologystring,  " ";
	print $outfile " ", $hit->protidstring, " Mr: ", sprintf( "%.4f", $monolinkpepmass );
	print $outfile "| Mr measured: ", sprintf( "%.4f error: ", $ms1parent ), sprintf( "%.1f ppm ", $self->calcerrorppm );
	print $outfile " xlinkions matched: ",    $xlinkmatchratio,  "\t";
	print $outfile "backbone ions matched: ", $commonmatchratio, "\t";
	print $outfile "xcorrx: ",   $self->{'xlinkxcorr'},    "\t";
	print $outfile "xcorrbb: ",  $self->{'backbonexcorr'}, "\t";
	print $outfile "xcorrall: ", $self->{'allxcorr'},      "\t";
	print $outfile " | prescore: ",    sprintf( "%.2f ", $self->getprescore );
	print $outfile " |  matchratio: ", sprintf( "%.6f ", $self->getmatchratio );

	#	print $outfile " |  PRMI: ",  sprintf( "%.2f ", $self->getPRMI );
	#	print $outfile " |  wPRMI: ", sprintf( "%.2f ", $self->getwPRMI );
	print $outfile " |  TICperc: ", sprintf( "%.2f ", $self->get_TICperc );
	print $outfile "score: ", $self->getscore, "\n";
}

sub _printxlinks_monolink
{
	my $self    = shift;
	my $outfile = shift;
	unless ( defined($outfile) )
	{
		$outfile = *STDOUT;
	}
	my $hit              = $self->gethitObjs->[0];
	my $combinations     = $self->gettopology->[0];
	my $monolinkmass     = $self->getxlinkermass;
	my $monolinkpepmass  = $hit->molweight + $monolinkmass;
	my $ms1parent        = $self->getms1mass;
	my $xlinkmatchratio  = join "\/", $self->get_number_of_Xlinkmatches, $self->get_number_of_XlinkIons;
	my $commonmatchratio = join "\/", $self->get_number_of_Commonmatches, $self->get_number_of_CommonIons;

	#
	print $outfile $self->getspectrumname, " | ";
	print $outfile $self->getid,           " | ";
	print $outfile $self->topologystring,  " | ";
	print $outfile " ", $hit->protidstring,, " | Mr: ", sprintf( "%.4f", $monolinkpepmass );
	print $outfile " | Mr measured: ", sprintf( " | %.4f error: ", $ms1parent ), sprintf( "%.1f ppm ", $self->calcerrorppm );
	print $outfile " | xlinkions matched: ",     $xlinkmatchratio,;
	print $outfile " | backbone ions matched: ", $commonmatchratio, " ";
	print $outfile " | xcorrx: ",                sprintf( "%.3f ", $self->{'xlinkxcorr'} );
	print $outfile " | xcorrbb: ",               sprintf( "%.3f ", $self->{'backbonexcorr'} );
	print $outfile " | xcorrall: ",              sprintf( "%.3f ", $self->{'allxcorr'} );
	print $outfile " | prescore: ",              sprintf( "%.2f ", $self->getprescore );
	print $outfile " |  matchratio: ",           sprintf( "%.6f ", $self->getmatchratio );

	#
	#	print $outfile " |  PRMI: ",  sprintf( "%.2f ", $self->getPRMI );
	#	print $outfile " |  wPRMI: ", sprintf( "%.2f ", $self->getwPRMI );
	print $outfile " |  TICperc: ", sprintf( "%.2f ", $self->get_TICperc );
	print $outfile " | score: ",    sprintf( "%.2f ", $self->getscore );
	print $outfile "\n";
}

sub _printxlinks_xlink
{
	my $self    = shift;
	my $outfile = shift;
	unless ( defined($outfile) )
	{
		$outfile = *STDOUT;
	}
	my $xlinkmatchratio        = join "\/", $self->get_number_of_Xlinkmatches,        $self->get_number_of_XlinkIons;
	my $xlinkmatchratio_alpha  = join "\/", $self->get_number_of_Xlinkmatches_alpha,  $self->get_number_of_XlinkIons_alpha;
	my $xlinkmatchratio_beta   = join "\/", $self->get_number_of_Xlinkmatches_beta,   $self->get_number_of_XlinkIons_beta;
	my $commonmatchratio       = join "\/", $self->get_number_of_Commonmatches,       $self->get_number_of_CommonIons;
	my $commonmatchratio_alpha = join "\/", $self->get_number_of_Commonmatches_alpha, $self->get_number_of_CommonIons_alpha;
	my $commonmatchratio_beta  = join "\/", $self->get_number_of_Commonmatches_beta,  $self->get_number_of_CommonIons_beta;
	my $hits                   = $self->gethitObjs;
	my $hit1                   = $hits->[0];
	my $hit2                   = $hits->[1];
	my $combinations           = $self->gettopology;
	my $xlinkermass            = $self->getxlinkermass;
	my $xlinkpepmass           = $hit1->molweight + $hit2->molweight + $xlinkermass;
	my $ms1parent              = $self->getms1mass;
	print $outfile $self->getspectrumname, " | ";
	print $outfile $self->getid,           " | ";
	print $outfile $self->topologystring,  " | ";
	print $outfile " ", $hit1->protidstring, "-", $hit2->protidstring, " | Mr: ", sprintf( "%.4f", $xlinkpepmass );
	print $outfile "| Mr measured: ", sprintf( " | %.4f error: ", $ms1parent ), sprintf( " | %.1f ppm ", $self->calcerrorppm );
	print $outfile " | xlinkions matched: ",           $xlinkmatchratio;
	print $outfile " | xlinkions_alpha matched: ",     $xlinkmatchratio_alpha;
	print $outfile " | xlinkions_beta matched: ",      $xlinkmatchratio_beta;
	print $outfile " | backbone ions matched: ",       $commonmatchratio;
	print $outfile " | backbone ions alpha matched: ", $commonmatchratio_alpha;
	print $outfile " | backbone ions beta matched: ",  $commonmatchratio_beta;
	print $outfile " | xcorrx: ",                      sprintf( "%.3f ", $self->{'xlinkxcorr'} );
	print $outfile " | xcorrbb: ",                     sprintf( "%.3f ", $self->{'backbonexcorr'} );
	print $outfile " |xcorrall: ",                     sprintf( "%.3f ", $self->{'allxcorr'} );
	print $outfile " |  prescore: ",                   sprintf( "%.2f ", $self->getprescore );
	print $outfile " |  matchratio: ",                 sprintf( "%.6f ", $self->getmatchratio );
	print $outfile " |  meanmatchintensity: ",         sprintf( "%.2f ", $self->get_meanmatchintensity );

	#	print $outfile " |  PRMI: ",  sprintf( "%.2f ", $self->getPRMI );
	#	print $outfile " |  wPRMI: ", sprintf( "%.2f ", $self->getwPRMI );
	print $outfile " |  TICperc: ", sprintf( "%.2f ", $self->get_TICperc );
	print $outfile " | score: ",    sprintf( "%.2f ", $self->getscore );
	print $outfile "\n";
}

sub formatseq
{
	my $self      = shift;
	my $seqstring = shift;
	my @seq       = split //, $seqstring;
	my $formatseq = "";
	my %seqhash;
	if ( $self->getParams->{'variable_mod'} )
	{
	#	my ( $AA, $delta ) = split /,|:/, $self->getParams->{'variable_mod'};
	#	my $roundeddelta = sprintf( "%.0f", $delta );
	#	$seqhash{'X'} = join "", $AA, "[", int($roundeddelta), "]";
		my @AAdelta = split /,|:/, $self->getParams->{'variable_mod'};
		my @AAmods = ('X', 'U', 'B', 'J');
		for my $i (0..($#AAdelta / 2 - 1) ){
			my $roundeddelta = sprintf( "%.0f", $AAdelta[2*$i + 1]);
			$seqhash{$AAmods[$i]} = join "", $AAdelta[2*$i], "[", int($roundeddelta), "]";
		}
		foreach my $aminoacid (@seq)
		{
			if ( $seqhash{$aminoacid} )
			{
				$formatseq .= $seqhash{$aminoacid};
			} else
			{
				$formatseq .= $aminoacid;
			}
		}
	} else
	{
		$formatseq = $seqstring;
	}
	return $formatseq;
}

sub get_number_of_Xlinkmatches
{
	my $self = shift;
	return $#{ $self->{'xlinkmatches'} } + 1;
}

sub get_number_of_Xlinkmatches_alpha
{
	my $self = shift;
	return $#{ $self->{'xlinkmatches_alpha'} } + 1;
}

sub get_number_of_Xlinkmatches_beta
{
	my $self = shift;
	return $#{ $self->{'xlinkmatches_beta'} } + 1;
}

sub get_number_of_XlinkIons
{
	my $self = shift;

	#	my $xlinkionhash = $self->{'xlinkions'};
	#	my $nions        = 0;
	#	foreach my $charge ( keys %$xlinkionhash ) {
	#		$nions += scalar( @{ $xlinkionhash->{$charge} } );
	#	}
	#	return $nions;
	return $self->{'nxlinkionsalpha'} + $self->{'nxlinkionsbeta'};
}

sub get_number_of_XlinkIons_alpha
{
	my $self = shift;

	#	my $xlinkionhash = $self->getXlinkIons_alpha;
	#	my $nions        = 0;
	#	foreach my $charge ( keys %$xlinkionhash ) {
	#		$nions += scalar( @{ $xlinkionhash->{$charge} } );
	#	}
	#	return $nions;
	return $self->{'nxlinkionsalpha'};
}

sub get_number_of_XlinkIons_beta
{
	my $self         = shift;
	my $xlinkionhash = $self->getXlinkIons_beta;
	my $nions        = 0;
	foreach my $charge ( keys %$xlinkionhash )
	{
		$nions += scalar( @{ $xlinkionhash->{$charge} } );
	}
	return $nions;
	return $self->{'nxlinkionsbeta'};
}

sub get_number_of_CommonIons
{
	my $self = shift;
	return $self->{'ncommonionsalpha'} + $self->{'ncommonionsbeta'};

	#return $#{ $self->{'commonions'} } + 1;
}

sub get_number_of_CommonIons_alpha
{
	my $self = shift;
	return $self->{'ncommonionsalpha'};

	#	my $commonionsalpha=$self->getCommonIons_alpha;
	#
	#	return $#$commonionsalpha + 1;
}

sub get_number_of_CommonIons_beta
{
	my $self = shift;

	#	my $commonionsbeta=$self->getCommonIons_beta;
	#	return $#$commonionsbeta + 1;
	return $self->{'ncommonionsbeta'};
}

sub get_number_of_CommonPeaks
{
	my $self = shift;
	return $self->getSpecObj->getncommonions;
}

sub get_number_of_XlinkPeaks
{
	my $self = shift;
	return $self->getSpecObj->getnxlinkions;
}

sub get_number_of_Commonmatches
{
	my $self = shift;
	return $#{ $self->{'commonmatches'} } + 1;
}

sub get_number_of_Commonmatches_alpha
{
	my $self = shift;
	return $#{ $self->{'commonmatches_alpha'} } + 1;
}

sub get_number_of_Commonmatches_beta
{
	my $self = shift;
	return $#{ $self->{'commonmatches_beta'} } + 1;
}

sub printxlinksXML
{
	my $self            = shift;
	my $outfile         = shift;
	my $rank            = shift;
	my $specfilename    = shift;
	my $webparams       = shift;
	my $printionmatches = shift;
	my $hits            = $self->gethitObjs;
	my $hit1            = $hits->[0];
	my $hit2            = $hits->[1];
	my $combinations    = $self->gettopology;
	my $xlinkermass     = $self->getxlinkermass; ## is defined individually for each linkobject
	my $xlinkpepmass    = 0;
	foreach my $hit (@$hits)
	{
		$xlinkpepmass += $hit->molweight;
	}
	$xlinkpepmass += $xlinkermass;

	#	my $xlinkpepmass = $hit1->molweight + $hit2->molweight + $xlinkermass;
	my $ms1parent        = $self->getms1mass;
	my $ms1mz            = $self->getprecursorMz;
	my $ms1charge        = $self->getprecursorCharge;
	my $xlinkmatchratio  = join "\/", $self->get_number_of_Xlinkmatches, $self->get_number_of_XlinkIons;
	my $commonmatchratio = join "\/", $self->get_number_of_Commonmatches, $self->get_number_of_CommonIons;
	print $outfile "<search_hit ";
	print $outfile "search_hit_rank=\"", $rank, "\" ";
	print $outfile "id=\"",   $self->getid,        "\" ";
	print $outfile "type=\"", $self->getxlinktype, "\" ";
	print $outfile "structure=\"", $self->formatseq( $hit1->seq ), $hit2 && ( "-", $self->formatseq( $hit2->seq ) ), "\" ";
	print $outfile "seq1=\"", $hit1->seq, "\" ";
	print $outfile "seq2=\"", $hit2 && $hit2->seq, "\" ";
	print $outfile "prot1=\"", $hit1->protidstring, "\" ";
	print $outfile "prot2=\"", $hit2 && $hit2->protidstring, "\" ";
	print $outfile "topology=\"", $self->topologystring, "\" ";
	print $outfile "xlinkposition=\"", ( join ",", @{ $self->gettopology } ), "\" ";
	print $outfile "Mr=\"", sprintf( "%.5f", $xlinkpepmass ), "\" ";
	print $outfile "mz=\"", sprintf( "%.5f", $ms1mz ),        "\" ";
	print $outfile "charge=\"",      $ms1charge,   "\" ";
	print $outfile "xlinkermass=\"", $xlinkermass, "\" ";
	print $outfile "measured_mass=\"", sprintf( "%.4f", $ms1parent ), "\" ";
	print $outfile "error=\"",     sprintf( "%.1f", $self->calcerrorppm ),        "\" ";
	print $outfile "error_rel=\"", sprintf( "%.1f", $self->calcerrorppmrelativ ), "\" ";
	print $outfile "xlinkions_matched=\"",    $self->{'nxlinkmatchratio'}, "\" ";
	print $outfile "backboneions_matched=\"", $self->{'commonmatchratio'}, "\" ";

	#	print $outfile "xlinkions_match_error_mean=\"",
	#	  $self->get_xlink_matcherror_mean, "\" ";
	#	print $outfile "xlinkions_match_error_stdev=\"",
	#	  $self->get_xlink_matcherror_stdev, "\" ";
	#
	#	print $outfile "backbone_match_error_mean=\"",
	#	  $self->get_common_matcherror_mean, "\" ";
	#	print $outfile "backbone_match_error_stdev=\"",
	#	  $self->get_common_matcherror_stdev, "\" ";
## Added by TW
	print $outfile "weighted_matchodds_mean=\"", $self->get_weightedmatchodds_mean, "\" ";
	print $outfile "weighted_matchodds_sum=\"",  $self->get_weightedmatchodds_sum,  "\" ";
## Added by TW
	print $outfile "match_error_mean=\"",  $self->get_matcherror_mean,  "\" ";
	print $outfile "match_error_stdev=\"", $self->get_matcherror_stdev, "\" ";
	print $outfile "xcorrx=\"", sprintf( "%.5f", $self->{'xlinkxcorr'} ),    "\" ";
	print $outfile "xcorrb=\"", sprintf( "%.5f", $self->{'backbonexcorr'} ), "\" ";
	print $outfile "match_odds=\"",             sprintf( "%.5f", $self->getmatchratio ),              "\" ";
	print $outfile "prescore=\"",               sprintf( "%.5f", $self->getprescore ),                "\" ";
	print $outfile "prescore_alpha=\"",         sprintf( "%.5f", $self->getprescore_alpha ),          "\" ";
	print $outfile "prescore_beta=\"",          sprintf( "%.5f", $self->getprescore_beta ),           "\" ";
	print $outfile "match_odds_alphacommon=\"", sprintf( "%.5f", $self->get_matchratio_alphacommon ), "\" ";
	print $outfile "match_odds_betacommon=\"",  sprintf( "%.5f", $self->get_matchratio_betacommon ),  "\" ";
	print $outfile "match_odds_alphaxlink=\"",  sprintf( "%.5f", $self->get_matchratio_alphaxlink ),  "\" ";
	print $outfile "match_odds_betaxlink=\"",   sprintf( "%.5f", $self->get_matchratio_betaxlink ),   "\" ";
	print $outfile "num_of_matched_ions_alpha=\"",        $self->get_num_of_matched_ions_alpha,     "\" ";
	print $outfile "num_of_matched_ions_beta=\"",         $self->get_num_of_matched_ions_beta,      "\" ";
	print $outfile "num_of_matched_common_ions_alpha=\"", $self->get_number_of_Commonmatches_alpha, "\" ";
	print $outfile "num_of_matched_common_ions_beta=\"",  $self->get_number_of_Commonmatches_beta,  "\" ";
	print $outfile "num_of_matched_xlink_ions_alpha=\"",  $self->get_number_of_Xlinkmatches_alpha,  "\" ";
	print $outfile "num_of_matched_xlink_ions_beta=\"",   $self->get_number_of_Xlinkmatches_beta,   "\" ";
	print $outfile "xcorrall=\"", sprintf( "%.5f", $self->{'allxcorr'} ), "\" ";

	#	print $outfile "PRMI=\"",  sprintf( "%.3f", $self->getPRMI ),  "\" ";
	#	print $outfile "wPRMI=\"", sprintf( "%.3f", $self->getwPRMI ), "\" ";
	print $outfile "TIC=\"", sprintf( "%.5f", $self->get_TICperc ), "\" ";
	print $outfile "TIC_alpha=\"", sprintf( "%.5f", $self->get_TICperc_alpha || 0 ), "\" ";
	print $outfile "TIC_beta=\"",  sprintf( "%.5f", $self->get_TICperc_beta  || 0 ), "\" ";
	print $outfile "wTIC=\"",      sprintf( "%.5f", $self->get_wTIC          || 0 ), "\" ";
## Added by TW
	print $outfile "intsum=\"",                  sprintf( "%.5f", $self->get_intsum_score            || 0 ), "\" ";
	print $outfile "apriori_match_probs=\"",     sprintf( "%.5f", $self->get_apriory_match_probs     || 0 ), "\" ";
	print $outfile "apriori_match_probs_log=\"", sprintf( "%.5f", $self->get_apriory_match_probs_log || 0 ), "\" ";
	print $outfile "series_score_mean=\"",       sprintf( "%.5f", $self->get_series_score            || 0 ), "\" ";
	if ($specfilename)
	{
		print $outfile "annotated_spec=\"", basename($specfilename), "\" ";
	} else
	{
		print $outfile "annotated_spec=\"\" ";
	}
	print $outfile "score=\"", sprintf( "%.2f", $self->getscore ), "\" >\n";
	if ($printionmatches)
	{
		$self->printmatchedionsXML($outfile);
	}
	print $outfile "<\/search_hit>\n";
}

sub printmatchedionsXML
{
	my $self                 = shift;
	my $outfile              = shift;
	my $matchingionpositions = $self->getMatchingIonpositions;
	my $iontable             = $self->getiontable;
	my $ion2peakmatchhash    = $self->get_ion2peakmatchhash;
	my $specObj              = $self->getSpecObj;

	print $outfile "<matched_ions source=\"xquest_v2.1.7\">\n";
	foreach my $iontype ( sort keys %$iontable )
	{
		foreach my $position ( sort keys %{ $iontable->{$iontype} } )
		{
			next unless defined( $matchingionpositions->{$iontype}->{$position} );
			my $ion_mz = $iontable->{$iontype}->{$position};
			next unless defined( $ion2peakmatchhash->{$ion_mz} );
			foreach my $peak ( @{ $ion2peakmatchhash->{$ion_mz} } )
			{
				my $deltamz  = $ion_mz - $peak;
				my $deltappm = $ion_mz ? 1e6 * ( $deltamz / $ion_mz ) : 0;
				my $intensity = $specObj->get_ionintensity($peak);
				print $outfile "<matched_ion ";
				print $outfile "ion_type=\"", xml_escape_gq($iontype), "\" ";
				print $outfile "position=\"", xml_escape_gq($position), "\" ";
				print $outfile "label=\"", xml_escape_gq( $iontype . $position ), "\" ";
				print $outfile "theoretical_mz=\"", sprintf( "%.6f", $ion_mz ), "\" ";
				print $outfile "observed_mz=\"", sprintf( "%.6f", $peak ), "\" ";
				print $outfile "delta_mz=\"", sprintf( "%.6f", $deltamz ), "\" ";
				print $outfile "delta_ppm=\"", sprintf( "%.3f", $deltappm ), "\" ";
				print $outfile "intensity=\"", sprintf( "%.6f", $intensity ), "\"";
				print $outfile "/>\n";
			}
		}
	}
	print $outfile "<\/matched_ions>\n";
}

sub xml_escape_gq
{
	my $value = shift;
	return "" unless defined($value);
	$value =~ s/&/&amp;/g;
	$value =~ s/"/&quot;/g;
	$value =~ s/</&lt;/g;
	$value =~ s/>/&gt;/g;
	return $value;
}

sub drawxlinkspec
{
	my $self          = shift;
	my $logscale      = shift;
	my $lossions      = shift;
	my $min           = shift;
	my $max           = shift;
	my $xlinkspecfile = shift;
	my $labelpeaks    = shift;
	my $showstructure = shift;
	my $printtextfile = shift;
	my $spectrum      = $self->getSpecObj;
	my $PARAMS        = $self->getParams;
	if ( !defined($min) )
	{
		$min = $PARAMS->{'minionsize'};
	}
	if ( !defined($max) )
	{
		$max = $PARAMS->{'maxionsize'};
	}
	my $specplot  = specplot->new();
	my @xlinkions = ();
	my $charge;
	my $commonpairs = $spectrum->getcommonpairs;
	if ($logscale)
	{
		my @tmp = ();
		foreach my $pairs (@$commonpairs)
		{
			push @tmp, [ $pairs->[0], sqrt( $pairs->[1] ) ];
		}
		$commonpairs = \@tmp;
	}
	my $mincharge = $self->minioncharge_xlinks;
	my $maxcharge = $self->maxioncharge_xlinks;
	for $charge ( $mincharge .. $maxcharge )
	{
		push @xlinkions, @{ $spectrum->getxlinkpairs($charge) };
	}
	push @xlinkions, @{ $spectrum->getxlinkpairs(0) };
	my $xlinkpeakpairs;
	if ($logscale)
	{
		my @tmp = ();
		foreach my $pairs (@xlinkions)
		{
			push @tmp, [ $pairs->[0], sqrt( $pairs->[1] ) ];
		}
		$xlinkpeakpairs = \@tmp;
	} else
	{
		$xlinkpeakpairs = \@xlinkions;
	}

	#$specplot->setcolor( "green", "red" );
	$specplot->plotdata( $min, $max, undef, $commonpairs, $xlinkpeakpairs );
	my ( @commonmatchpairs, @commonlosspairs, @xlinkmatchpairs, @xlinklosspairs );
## get the commonmatches
	my $commonmatches = $self->getcommonmatches;
	foreach my $matchedionmz (@$commonmatches)
	{
		my $matchedpeak = $self->get_ion2peakmatch($matchedionmz);    ## returns now an arrayref
		foreach my $matchedpeak (@$matchedpeak)
		{

			#print "<br>".$matchedionmz."\n";
			#print $matchedpeak."\n";
			#my $matchedpeak=$matchedionmz;
			#my $label       = $self->annotate_matchedion($matchedionmz);
			my $label = $self->annotate_matchedpeak($matchedpeak);

			#print Dumper ($label);
			$label =~ s/plus/\+/;
			$label =~ s/alpha/a/;
			$label =~ s/beta/b/;
			$label =~ s/standard_//;
			my $peakintensity = 0;
			if ($logscale)
			{
				$peakintensity = sqrt( $spectrum->get_ionintensity($matchedpeak) );
			} else
			{
				$peakintensity = $spectrum->get_ionintensity($matchedpeak);

				#$peakintensity = $spectrum->get_realionintensity($matchedpeak);
			}
			if ( $label =~ /(H2O|NH3)/ )
			{
				push @commonlosspairs, [ $matchedpeak, $peakintensity, int($matchedpeak) . "_" . $label ];
			} else
			{
				push @commonmatchpairs, [ $matchedpeak, $peakintensity, int($matchedpeak) . "_" . $label ];
			}
		}
	}
	my $xlinkmatches = $self->getxlinkmatches;

	#print Dumper ($xlinkmatches);
	foreach my $matchedionmz (@$xlinkmatches)
	{
		my $matchedpeak = $self->get_ion2peakmatch($matchedionmz);    ## is now an arrayreference -> more peaks can match to 1 ion
		foreach my $matchedpeak (@$matchedpeak)
		{

			#print "Looking for matchedpeak:: $matchedpeak, ion: $matchedionmz<br>";
			#print Dumper($self->get_ion_annotationtable->{$matchedionmz});
			my $label = $self->annotate_matchedpeak($matchedpeak);

			#print Dumper ($label)."<br>";
			#exit;
			$label =~ s/plus/\+/;
			$label =~ s/alpha/a/;
			$label =~ s/beta/b/;
			$label =~ s/standard_//;
			my $peakintensity = 0;
			if ($logscale)
			{
				$peakintensity = sqrt( $spectrum->get_ionintensity($matchedpeak) );
			} else
			{
				$peakintensity = $spectrum->get_ionintensity($matchedpeak);

				#$peakintensity = $spectrum->get_realionintensity($matchedpeak);
			}
			if ( $label =~ /(H2O|NH3)/ )
			{
				push @xlinklosspairs, [ $matchedpeak, $peakintensity, int($matchedpeak) . "_" . $label ];
			} else
			{
				push @xlinkmatchpairs, [ $matchedpeak, $peakintensity, int($matchedpeak) . "_" . $label ];
			}
		}
	}
	if ($labelpeaks)
	{
		$specplot->labelpeaks( [ "annotation", "annotation", "annotation", "annotation" ], [ "green", "red", "blue", "lightblue" ], [ \@commonmatchpairs, \@xlinkmatchpairs, \@commonlosspairs, \@xlinklosspairs ] );
		my %colorhash = (
						  "common-peaks"                  => "green",
						  "xlink-peaks"                   => "red",
						  "labels indicate matched peaks" => "grey",
		);
		my @sortlist = ( "common-peaks", "xlink-peaks", "labels indicate matched peaks" );
		$specplot->drawlegend( 600, 20, \%colorhash, \@sortlist );
	} else
	{
		$specplot->labelpeaks( [ "diamond", "diamond", "cross", "cross" ], [ "green", "red", "blue", "lightblue" ], [ \@commonmatchpairs, \@xlinkmatchpairs, \@commonlosspairs, \@xlinklosspairs ] );
		my %colorhash = (
						  "common-peaks"                    => "green",
						  "xlink-peaks"                     => "red",
						  "diamonds indicate matched peaks" => "grey",
		);
		my @sortlist = ( "common-peaks", "xlink-peaks", "diamonds indicate matched peaks" );
		$specplot->drawlegend( 600, 20, \%colorhash, \@sortlist );
	}
	if ($printtextfile)
	{
		my $xlinktextfile_common_matches = $printtextfile . '_commonmatches.txt';
		my $xlinktextfile_xlink_matches  = $printtextfile . '_xlinkmatches.txt';
		my $xlinktextfile_common_peaks   = $printtextfile . '_commonpeaks.txt';
		my $xlinktextfile_xlink_peaks    = $printtextfile . '_xlinkpeaks.txt';
		open COMMONPEAKS, ">$xlinktextfile_common_peaks" or die $!;

		#	print "open $xlinktextfile_common_peaks <br>";
		open XLINKPEAKS,    ">$xlinktextfile_xlink_peaks"    or die $!;
		open COMMONMATCHES, ">$xlinktextfile_common_matches" or die $!;
		open XLINKMATCHES,  ">$xlinktextfile_xlink_matches"  or die $!;
		foreach my $pair (@commonmatchpairs)
		{
			print COMMONMATCHES join "\t", @$pair;
			print COMMONMATCHES "\n";
		}
		foreach my $pair (@xlinkmatchpairs)
		{
			print XLINKMATCHES join "\t", @$pair;
			print XLINKMATCHES "\n";
		}
		foreach my $pair (@$commonpairs)
		{
			print COMMONPEAKS join "\t", @$pair;
			print COMMONPEAKS "\n";
		}
		foreach my $pair (@$xlinkpeakpairs)
		{
			print XLINKPEAKS join "\t", @$pair;
			print XLINKPEAKS "\n";
		}
		close(COMMONMATCHES);
		close(XLINKMATCHES);
		close(COMMONPEAKS);
		close(XLINKPEAKS);
	}
	if ($showstructure)
	{
		$self->drawpepstructure( \@xlinkions, $commonpairs, $xlinkspecfile, $logscale, $lossions, $specplot->img );
	}
	unless ( defined($xlinkspecfile) )
	{
		my $id = $self->getid;
		$id =~ s/::/_/g;
		$xlinkspecfile = File::Spec->catfile( $self->getParams->{'outputpath'}, ( join "", $spectrum->getspecbasename, $id, ".png" ) );
	}

	#	print "open $xlinkspecfile <br>";
	$specplot->printimage($xlinkspecfile);
	return $xlinkspecfile;
}

sub xcorrelation_common
{
	my $self         = shift;
	my $basefilename = shift;
	my $printout     = shift;
	my $delay        = shift;
	my $precision    = shift;
	my $lossions     = shift;
	Xcorrelation::xcorrelation_common( $self, $basefilename, $printout, $delay, $precision, $lossions );
}

sub xcorrelation_common_normalized
{
	my $self         = shift;
	my $basefilename = shift;
	my $printout     = shift;
	my $delay        = shift;
	my $precision    = shift;
	my $lossions     = shift;
	Xcorrelation::xcorrelation_common_normalized( $self, $basefilename, $printout, $delay, $precision, $lossions );
}

sub xcorrelation_xlink
{
	my $self         = shift;
	my $basefilename = shift;
	my $printout     = shift;
	my $delay        = shift;
	my $precision    = shift;
	my $lossions     = shift;
	Xcorrelation::xcorrelation_xlink( $self, $basefilename, $printout, $delay, $precision, $lossions );
}

sub xcorrelation_xlink_normalized
{
	my $self         = shift;
	my $basefilename = shift;
	my $printout     = shift;
	my $delay        = shift;
	my $precision    = shift;
	my $lossions     = shift;
	Xcorrelation::xcorrelation_xlink_normalized( $self, $basefilename, $printout, $delay, $precision, $lossions );
}

sub calcscore
{
	my $self             = shift;
	my $PARAMS           = $self->getParams;
	my $xcorrweigth      = $PARAMS->{'xcorrxweight'};
	my $backboneweigth   = $PARAMS->{'xcorrbweight'};
	my $ppmweigth        = $PARAMS->{'ppmweight'};
	my $matchratioweigth = $PARAMS->{'matchratioweight'};
	my $wPRMIweight      = $PARAMS->{'wPRMIweight'};
	my $matchoddsweight  = $PARAMS->{'matchoddsweight'};
	my $wmatchoddsweight = $PARAMS->{'wmatchoddsweight'};
	my $wTICweight       = $PARAMS->{'wTICweight'};
	my $TICweight        = $PARAMS->{'TICweight'};
	my $seriesweight     = $PARAMS->{'seriesweight'};
	my $intsumweight     = $PARAMS->{'intsumweight'};
	## old ppmweight (never used)
	## $ppmweigth * log( 1 + ( 1 / $self->calcerrorppm ) )
	#print "Error rel ppm: ", $self->calcerrorppmrelativ, "<br>";
	$self->{'score'} =
	  ( $ppmweigth * $self->calcerrorppmrelativ ) +
	  ( $xcorrweigth * $self->{'xlinkxcorr'} ) +
	  ( $backboneweigth * $self->{'backbonexcorr'} ) +
	  ( $matchoddsweight * $self->getmatchratio ) +
	  ( $wmatchoddsweight * $self->get_weightedmatchodds_mean ) +
	  ( $TICweight * $self->get_TICperc ) +
	  ( $wTICweight * $self->get_wTIC ) +
	  ( $seriesweight * $self->get_series_score ) +
	  ( $intsumweight * $self->get_intsum_score );

	#print "Score: ", $self->{'score'}, "<br>";
	return $self->{'score'};
}

sub print_subscores_cmdline
{
	my $self             = shift;
	my $PARAMS           = $self->getParams;
	my $xcorrweigth      = $PARAMS->{'xcorrxweight'};
	my $backboneweigth   = $PARAMS->{'xcorrbweight'};
	my $ppmweigth        = $PARAMS->{'ppmweight'};
	my $matchratioweigth = $PARAMS->{'matchratioweight'};
	my $wPRMIweight      = $PARAMS->{'wPRMIweight'};
	my $matchoddsweight  = $PARAMS->{'matchoddsweight'};
	my $wTICweight       = $PARAMS->{'wTICweight'};
	my $TICweight        = $PARAMS->{'TICweight'};
	my $mr               = $self->calc_mr;
	my $ms1parent        = $self->getms1mass;
	my $errorrel         = $self->calcerrorppmrelativ;
	print "\nFullscore:" . $self->calcscore . "\n";
	print "Subscores:\n";
	print "Xcorrx:$self->{'xlinkxcorr'}\n";
	print "Xcorrb:$self->{'backbonexcorr'}\n";
	print "MatchOdds:" . $self->getmatchratio . "\n";
	print "Weighted MatchOdds Mean: ", $self->get_weightedmatchodds_mean . "\n";

	#print "Weighted MatchOdds Sum: ",$self->get_weightedmatchodds_sum . "\n";
	print "TIC:" . $self->get_TICperc . "\n";
	print "wTIC:" . $self->get_wTIC . "\n";
	print "Precursor mass:\t $ms1parent\n";
	print "MR calculated:\t $mr\n";
	print "Error [ppm]:\t $errorrel\n";
	return;
}

sub print_subscores
{
	my $self             = shift;
	my $PARAMS           = $self->getParams;
	my $xcorrweigth      = $PARAMS->{'xcorrxweight'};
	my $backboneweigth   = $PARAMS->{'xcorrbweight'};
	my $ppmweigth        = $PARAMS->{'ppmweight'};
	my $matchratioweigth = $PARAMS->{'matchratioweight'};
	my $wPRMIweight      = $PARAMS->{'wPRMIweight'};
	my $matchoddsweight  = $PARAMS->{'matchoddsweight'};
	my $wTICweight       = $PARAMS->{'wTICweight'};
	my $TICweight        = $PARAMS->{'TICweight'};
	my $mr               = $self->calc_mr;
	my $ms1parent        = $self->getms1mass;
	my $errorrel         = $self->calcerrorppmrelativ;
	
	print "<br>Fullscore:" . $self->calcscore . "<br>";
	print "Subscores:<br>";
	print "Xcorrx:$self->{'xlinkxcorr'}<br>";
	print "Xcorrb:$self->{'backbonexcorr'}<br>";
	print "MatchOdds:" . $self->getmatchratio . "<br>";
	#print "Weighted MatchOdds Mean: ", $self->get_weightedmatchodds_mean . "<br>";
	#print "Weighted MatchOdds Sum: ",$self->get_weightedmatchodds_sum . "<br>";
	print "TIC:" . $self->get_TICperc . "<br>";
	print "wTIC:" . $self->get_wTIC . "<br>";

	
	## Get peptides
	my $hitobjects = $self->gethitObjs;		
	print "Precursor mass:\t $ms1parent<br>";
	my $xlinkermass= $self->getxlinkermass;
	print "X-LINKER mass:\t $xlinkermass<br>";

	if ($self->xlinktype eq "monolink" or $self->xlinktype eq "intralink"){
	my $pepobject1mass= $hitobjects->[0]->molweight;
	print "Pepobject1mass: $pepobject1mass<br>";	
	}else{
	my $pepobject1mass= $hitobjects->[0]->molweight;
	print "Pepobject1mass: $pepobject1mass<br>";
	my $pepobject2mass= $hitobjects->[1]->molweight;
	print "Pepobject2mass: $pepobject2mass<br>";		
	}	

	print "MR calculated:\t $mr<br>";
	print "Error [ppm]:\t $errorrel<br>";	
	
	return;
}

sub getscore
{
	my $self = shift;
	return $self->{'score'};
}

sub drawxcorrgraph
{
	my $dat       = shift;
	my $precision = shift;
	my $filelabel = shift;
	my $title     = shift;
	if ( defined($dat) )
	{
		my ( @x, @y );
		foreach (@$dat)
		{
			push @x, $_->[0] / $precision;
			push @y, $_->[1];
		}

		# print "@data\n";
		my $graph = GD::Graph::lines->new( 400, 300 );
		$graph->set(
					 x_label      => 'delta mz',
					 x_label_skip => 5,
					 zero_axis    => 1,
					 t_margin     => 20,
					 b_margin     => 20,
					 l_margin     => 20,
					 r_margin     => 20,
					 line_width   => 3,
					 dclrs        => [qw(black red blue cyan)],
					 y_label      => 'cross-link  cross-correlation',
					 title        => "$title",
		) or die $graph->error;
		my $gd = $graph->plot( [ \@x, \@y ] ) or die $graph->error;
		open( IMG, ">$filelabel" ) or die $!;
		binmode IMG;
		print IMG $gd->png;
		close IMG;
	}
}

sub _drawpepstructure_xlink
{
	my $self      = shift;
	my $dat1      = shift;
	my $dat2      = shift;
	my $filelabel = shift;
	my $logscale  = shift;
	my $lossions  = shift;
	my $gd        = shift;
	my $specObj   = $self->getSpecObj;
	my $score     = sprintf( "%.2f", $self->getscore );
	my ( $i, $j );
	my $iontable  = $self->getMatchingIonsIDs_for_pepstructure;
	my $pepObj1   = $self->gethitObjs->[0];
	my $pepObj2   = $self->gethitObjs->[1];
	my $seq1      = $pepObj1->seq;
	my $seq2      = $pepObj2->seq;
	my $topology  = $self->gettopology;
	my @seqarray1 = split //, $seq1;
	my @seqarray2 = split //, $seq2;
	my ( $x1, $x2, $y1, $y2 );
	my $xposition        = 20;
	my $yposition1       = 20;
	my $yposition2       = 60;
	my $xlinkpositions   = $self->gettopology;
	my $xlinkpos1        = $xlinkpositions->[0];
	my $xlinkpos2        = $xlinkpositions->[1];
	my $reversexlinkpos1 = $self->getreversetopology( $pepObj1, $xlinkpos1 );
	my $reversexlinkpos2 = $self->getreversetopology( $pepObj2, $xlinkpos2 );
	my $blue             = $gd->colorAllocate( 0, 0, 255 );
	my $green            = $gd->colorAllocate( 0, 255, 0 );
	my $red              = $gd->colorAllocate( 255, 0, 0 );
	my $white            = $gd->colorAllocate( 255, 255, 255 );
	my %colorcode = (
					  'common' => $green,
					  'xlink'  => $red,
	);
	my $type = "common";

	for $i ( 0 .. $#seqarray1 )
	{
		$gd->string( gdLargeFont, $xposition + $i * 12, $yposition1, $seqarray1[$i], $blue );
		if ( defined( $iontable->{'alpha'}->{'fwd'}->{ $i + 1 } ) )
		{
			if ( $i + 1 < $xlinkpos1 )
			{
				$type = 'common';
			} else
			{
				$type = 'xlink';
			}
			$gd->line( $xposition + 2 + $i * 12,     $yposition1 - 6, $xposition + 2 + 4 + $i * 12, $yposition1 - 6, $colorcode{$type} );
			$gd->line( $xposition + 2 + 4 + $i * 12, $yposition1 - 6, $xposition + 2 + 4 + $i * 12, $yposition1 - 2, $colorcode{$type} );

			#		print "balpha ",$i+1,"\n";
		}
		if ( defined( $iontable->{'alpha'}->{'rev'}->{ $i + 1 } ) )
		{
			if ( $i + 1 < $reversexlinkpos1 )
			{
				$type = 'common';
			} else
			{
				$type = 'xlink';
			}
			$gd->line( $xposition + ( $#seqarray1 - $i ) * 12, $yposition1 + 24, $xposition + 4 + ( $#seqarray1 - $i ) * 12, $yposition1 + 24, $colorcode{$type} );
			$gd->line( $xposition + ( $#seqarray1 - $i ) * 12, $yposition1 + 24, $xposition + ( $#seqarray1 - $i ) * 12, $yposition1 + 20, $colorcode{$type} );

			#		print "yalpha ",$i,"\n";
		}
		if ( $i == ( $topology->[0] - 1 ) )
		{
			$x1 = $xposition + $i * 12 + 2;
			$y1 = $yposition1 + 15;
		}
	}
	for $i ( 0 .. $#seqarray2 )
	{
		$gd->string( gdLargeFont, $xposition + $i * 12, $yposition2, $seqarray2[$i], $blue );
		if ( defined( $iontable->{'beta'}->{'fwd'}->{ $i + 1 } ) )
		{
			if ( $i + 1 < $xlinkpos2 )
			{
				$type = 'common';
			} else
			{
				$type = 'xlink';
			}
			$gd->line( $xposition + 2 + $i * 12,     $yposition2 - 6, $xposition + 2 + 4 + $i * 12, $yposition2 - 6, $colorcode{$type} );
			$gd->line( $xposition + 2 + 4 + $i * 12, $yposition2 - 6, $xposition + 2 + 4 + $i * 12, $yposition2 - 2, $colorcode{$type} );

			#		print "balpha ",$i+1,"\n";
		}
		if ( defined( $iontable->{'beta'}->{'rev'}->{ $i + 1 } ) )
		{
			if ( $i + 1 < $reversexlinkpos2 )
			{
				$type = 'common';
			} else
			{
				$type = 'xlink';
			}
			$gd->line( $xposition + ( $#seqarray2 - $i ) * 12, $yposition2 + 24, $xposition + 4 + ( $#seqarray2 - $i ) * 12, $yposition2 + 24, $colorcode{$type} );
			$gd->line( $xposition + ( $#seqarray2 - $i ) * 12, $yposition2 + 24, $xposition + ( $#seqarray2 - $i ) * 12, $yposition2 + 20, $colorcode{$type} );
		}
		if ( $i == ( $topology->[1] - 1 ) )
		{
			$x2 = $xposition + $i * 12 + 2;
			$y2 = $yposition2;
		}
	}

	#$gd->dashedLine( $x1, $y1, $x2, $y2, $blue );
	$gd->line( $x1, $y1, $x2, $y2, $blue );

	#	$gd->string( gdSmallFont, $xposition + 400, 20, "x-linker",      $red );
	#	$gd->string( gdSmallFont, $xposition + 400, 30, "backbone ions", $green );
	#	$gd->string( gdSmallFont, $xposition + 400,
	#		40, ( join "", "intensity = ", $self->getSpecObj->basepeakintensity ),
	#		$red );
	#	open( IMG, ">$filelabel" ) or die $!;
	#	binmode IMG;
	#	print IMG $gd->png;
	#	close IMG;
}

sub _drawpepstructure_monolink
{
	my $self      = shift;
	my $dat1      = shift;
	my $dat2      = shift;
	my $filelabel = shift;
	my $logscale  = shift;
	my $lossions  = shift;
	my $gd        = shift;
	my $specObj   = $self->getSpecObj;
	my $score     = sprintf( "%.2f", $self->getscore );
	my ( $i, $j );
	my $iontable  = $self->getMatchingIonsIDs_for_pepstructure;
	my $pepObj1   = $self->gethitObjs->[0];
	my $seq1      = $pepObj1->seq;
	my @seqarray1 = split //, $seq1;
	my ( $x1, $y1, $x2, $y2 );
	my $xposition        = 20;
	my $yposition1       = 20;
	my $yposition2       = 60;
	my $xlinkpositions   = $self->gettopology;
	my $xlinkpos1        = $xlinkpositions->[0];
	my $reversexlinkpos1 = $self->getreversetopology( $pepObj1, $xlinkpos1 );
	my $blue             = $gd->colorAllocate( 0, 0, 255 );
	my $green            = $gd->colorAllocate( 0, 255, 0 );
	my $red              = $gd->colorAllocate( 255, 0, 0 );
	my $white            = $gd->colorAllocate( 255, 255, 255 );
	my %colorcode = (
					  'common' => $green,
					  'xlink'  => $red,
	);
	my $type = "common";

	for $i ( 0 .. $#seqarray1 )
	{
		$gd->string( gdLargeFont, $xposition + $i * 12, $yposition1, $seqarray1[$i], $blue );
		if ( defined( $iontable->{'alpha'}->{'fwd'}->{ $i + 1 } ) )
		{
			if ( $i + 1 < $xlinkpos1 )
			{
				$type = 'common';
			} else
			{
				$type = 'xlink';
			}
			$gd->line( $xposition + 2 + $i * 12,     $yposition1 - 6, $xposition + 2 + 4 + $i * 12, $yposition1 - 6, $colorcode{$type} );
			$gd->line( $xposition + 2 + 4 + $i * 12, $yposition1 - 6, $xposition + 2 + 4 + $i * 12, $yposition1 - 2, $colorcode{$type} );

			#		print "balpha ",$i+1,"\n";
		}
		if ( defined( $iontable->{'alpha'}->{'rev'}->{ $i + 1 } ) )
		{
			if ( $i + 1 < $reversexlinkpos1 )
			{
				$type = 'common';
			} else
			{
				$type = 'xlink';
			}
			$gd->line( $xposition + ( $#seqarray1 - $i ) * 12, $yposition1 + 24, $xposition + 4 + ( $#seqarray1 - $i ) * 12, $yposition1 + 24, $colorcode{$type} );
			$gd->line( $xposition + ( $#seqarray1 - $i ) * 12, $yposition1 + 24, $xposition + ( $#seqarray1 - $i ) * 12, $yposition1 + 20, $colorcode{$type} );

			#		print "yalpha ",$i,"\n";
		}
		if ( $i == ( $xlinkpos1 - 1 ) )
		{
			$x1 = $xposition + $i * 12 + 2;
			$x2 = $x1;
			$y1 = $yposition1 + 15;
			$y2 = $y1 + 30;
		}
	}
	$gd->dashedLine( $x1, $y1, $x2, $y2, $blue );
	$gd->rectangle( $x1 - 2, $y2, $x2 + 2, $y2 + 5, $blue );

	#	$gd->string( gdSmallFont, $xposition + 400, 20, "x-linker",      $red );
	#	$gd->string( gdSmallFont, $xposition + 400, 30, "backbone ions", $green );
	#	$gd->string( gdSmallFont, $xposition + 400,
	#		40, ( join "", "intensity = ", $self->getSpecObj->basepeakintensity ),
	#		$red );
	#	open( IMG, ">$filelabel" ) or die $!;
	#	binmode IMG;
	#	print IMG $gd->png;
	#	close IMG;
}

sub _drawpepstructure_intralink
{
	my $self      = shift;
	my $dat1      = shift;
	my $dat2      = shift;
	my $filelabel = shift;
	my $logscale  = shift;
	my $lossions  = shift;
	my $gd        = shift;
	my $specObj   = $self->getSpecObj;
	my $score     = sprintf( "%.2f", $self->getscore );
	my ( $i, $j );
	my $iontable = $self->getMatchingIonsIDs_for_pepstructure;
	my $pepObj1  = $self->gethitObjs->[0];

	#	my $pepObj2  = $self->gethitObjs->[1];
	my $seq1 = $pepObj1->seq;

	#my $seq2 = $pepObj2->seq;
	my $topology = $self->gettopology;
	my @seqarray1 = split //, $seq1;

	#my @seqarray2 = split //, $seq2;
	my ( $x1a, $x2a, $y1a, $y2a, $x1b, $x2b, $y1b, $y2b );
	my $xposition        = 20;
	my $yposition1       = 20;
	my $yposition2       = 60;
	my $xlinkpositions   = $self->gettopology;
	my $xlinkpos1        = $xlinkpositions->[0];
	my $xlinkpos2        = $xlinkpositions->[1];
	my $reversexlinkpos1 = $self->getreversetopology( $pepObj1, $xlinkpos1 );

	#	my $reversexlinkpos2 = $self->getreversetopology( $pepObj2, $xlinkpos2 );
	my $blue  = $gd->colorAllocate( 0,   0,   255 );
	my $green = $gd->colorAllocate( 0,   255, 0 );
	my $red   = $gd->colorAllocate( 255, 0,   0 );
	my $white = $gd->colorAllocate( 255, 255, 255 );
	my %colorcode = (
					  'common' => $green,
					  'xlink'  => $red,
	);
	my $type = "common";
	for $i ( 0 .. $#seqarray1 )
	{
		$gd->string( gdLargeFont, $xposition + $i * 12, $yposition1, $seqarray1[$i], $blue );
		if ( defined( $iontable->{'alpha'}->{'fwd'}->{ $i + 1 } ) )
		{
			if ( $i + 1 < $xlinkpos1 )
			{
				$type = 'common';
			} else
			{
				$type = 'xlink';
			}
			$gd->line( $xposition + 2 + $i * 12,     $yposition1 - 6, $xposition + 2 + 4 + $i * 12, $yposition1 - 6, $colorcode{$type} );
			$gd->line( $xposition + 2 + 4 + $i * 12, $yposition1 - 6, $xposition + 2 + 4 + $i * 12, $yposition1 - 2, $colorcode{$type} );

			#		print "balpha ",$i+1,"\n";
		}
		if ( defined( $iontable->{'alpha'}->{'rev'}->{ $i + 1 } ) )
		{
			if ( $i + 1 < $reversexlinkpos1 )
			{
				$type = 'common';
			} else
			{
				$type = 'xlink';
			}
			$gd->line( $xposition + ( $#seqarray1 - $i ) * 12, $yposition1 + 24, $xposition + 4 + ( $#seqarray1 - $i ) * 12, $yposition1 + 24, $colorcode{$type} );
			$gd->line( $xposition + ( $#seqarray1 - $i ) * 12, $yposition1 + 24, $xposition + ( $#seqarray1 - $i ) * 12, $yposition1 + 20, $colorcode{$type} );

			#		print "yalpha ",$i,"\n";
		}
		if ( $i == ( $topology->[0] - 1 ) )
		{
			$x1a = $xposition + $i * 12 + 2;
			$x2a = $x1a;
			$y1a = $yposition1 + 15;
			$y2a = $y1a + 25;
			$gd->dashedLine( $x1a, $y1a, $x2a, $y2a, $blue );
		}
		if ( $i == ( $topology->[1] - 1 ) )
		{
			$x1b = $xposition + $i * 12 + 2;
			$x2b = $x1b;
			$y1b = $yposition1 + 15;
			$y2b = $y1b + 25;
			$gd->dashedLine( $x1b, $y1b, $x2b, $y2b, $blue );
		}
	}
	$gd->dashedLine( $x1a, $y2a, $x1b, $y2b, $blue );

	#
	#	$gd->string( gdSmallFont, $xposition + 400, 20, "x-linker",      $red );
	#	$gd->string( gdSmallFont, $xposition + 400, 30, "backbone ions", $green );
	#	$gd->string( gdSmallFont, $xposition + 400,
	#		40, ( join "", "intensity = ", $self->getSpecObj->basepeakintensity ),
	#		$red );
	#	open( IMG, ">$filelabel" ) or die $!;
	#	binmode IMG;
	#	print IMG $gd->png;
	#	close IMG;
}

sub get_common_matcherror_mean
{
	my $self = shift;
	return $self->{'common_matcherror_mean'};
}

sub get_common_matcherror_stdev
{
	my $self = shift;
	return $self->{'common_matcherror_stdev'};
}
#######get charges for ionseries xlinks charges are like for common ions in monolinks
sub calc_maxioncharge_xlinks
{
	my $self = shift;
	if ( $self->getxlinktype eq "xlink" )
	{
		my $definedmaxioncharge = $self->getParams->{'ioncharge_xlink'}->[-1];
		my $spectrumcharge      = $self->getSpecObj->getprecursorCharge;
		if ( $definedmaxioncharge > $self->getSpecObj->getprecursorCharge )
		{
			$self->{'maxioncharge_xlinks'} = $spectrumcharge;
		} else
		{
			$self->{'maxioncharge_xlinks'} = $definedmaxioncharge;
		}
	} else
	{
		$self->{'maxioncharge_xlinks'} = $self->maxioncharge_common;
	}
}

sub maxioncharge_xlinks
{
	my $self = shift;
	return $self->{'maxioncharge_xlinks'};
}

sub calc_minioncharge_xlinks
{
	my $self = shift;
	if ( $self->getxlinktype eq "xlink" )
	{
		$self->{'minioncharge_xlinks'} = $self->getParams->{'ioncharge_xlink'}->[0];
	} else
	{
		$self->{'minioncharge_xlinks'} = $self->minioncharge_common;
	}
}

sub minioncharge_xlinks
{
	my $self = shift;
	return $self->{'minioncharge_xlinks'};
}

sub calc_maxioncharge_common
{
	my $self                = shift;
	my $definedmaxioncharge = $self->getParams->{'ioncharge_common'}->[-1];
	my $spectrumcharge      = $self->getSpecObj->getprecursorCharge;
	if ( $definedmaxioncharge > $self->getSpecObj->getprecursorCharge )
	{
		$self->{'maxioncharge_common'} = $spectrumcharge;
	} else
	{
		$self->{'maxioncharge_common'} = $definedmaxioncharge;
	}
}

sub maxioncharge_common
{
	my $self = shift;
	return $self->{'maxioncharge_common'};
}

sub calc_minioncharge_common
{
	my $self = shift;
	$self->{'minioncharge_common'} = $self->getParams->{'ioncharge_common'}->[0];
}

sub minioncharge_common
{
	my $self = shift;
	return $self->{'minioncharge_common'};
}
#####################################################33
sub get_xlink_matcherror_mean
{
	my $self = shift;
	return $self->{'xlink_matcherror_mean'};
}

sub get_xlink_matcherror_stdev
{
	my $self = shift;
	return $self->{'xlink_matcherror_stdev'};
}

sub getcommonpairs
{
	my $self = shift;
	return $self->getSpecObj->{'ionpairs'};
}

sub getcommonpeaks
{
	my $self = shift;
	return $self->getSpecObj->{'peaks'};
}
1;
