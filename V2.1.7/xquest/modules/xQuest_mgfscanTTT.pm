package xQuest_mgfscanTTT;
use strict;
#---------------------------------------------------------------------------
# Module: xQuest_mgfscanTTT.pm
# Author(s): Thomas Walzthoeni, xQuest specific modules are based on 
# original work by Oliver Rinner.
# Description: Module for handling mgf scans of TTT instrument.
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
use Carp;

my $Hatom = 1.0078250;

sub new
{
	my $class     = shift;
	my $scanarray = shift;
	my $self      = {};
	bless $self, $class;
	$self->{'scanarray'} = $scanarray;
	unless ( $self->_init == 1 )
	{
		confess "Error: Scan not correctly initialized\n";
		exit;
	}
	return $self;
}

sub _init
{
	my $self      = shift;
	my $scanarray = $self->{'scanarray'};

	#print Dumper ($scanarray);
	my $s = @$scanarray;
	my $peaks;
	for ( my $i = 0 ; $i < $s ; $i++ )
	{
		chomp($scanarray->[$i]);
		
		if ( $i == 0 )
		{

			# do nothing
		} elsif ( $i == 1 )
		{			
			my @split = split( /=/, $scanarray->[$i] );
			my @split2=split( /:/, $split[1] );
			$self->{'title'} = $split2[1];
		} elsif ( $i == 2 )
		{
			my @split  = split( /=/, $scanarray->[$i] );
			my @split2 = split( //,  $split[1] );
			$self->{'z'} = $split2[0];
		} elsif ( $i == 3 )
		{
			my @split  = split( /=/, $scanarray->[$i] );
			my @split2 = split( / /, $split[1] );
			$self->{'mz'} = $split2[0];
		} elsif ( $i == 4 )
		{
			my @split  = split( /=/, $scanarray->[$i] );
			$self->{'rtsec'} = $split[1];
		}else
		{
			unless ( $scanarray->[$i] =~ m/END IONS/ )
			{
				my @split = split( / /, $scanarray->[$i] );
				unless($split[0] && $split[1]){ # empty line
				next;
				}
				unless($split[0] || $split[1]){
				warn "Parsing mgf: No mz or intensity value found, line skipped\n";
				next;
				}
				push @{ $self->{'peaks'} }, [ $split[0], $split[1] ];
			}
		}
		#print $scanarray->[$i] . "\n";
	}
	return 1;
}

sub print_matchlist_light{
my $self = shift;
my $id = $self->title.",".$self->title;
my $matchlist = $id."\t".$self->precursor_mz."\t".$self->charge."\t". $self->title."\t". $self->title."\tlight\tlight\n";
return $matchlist;
}

sub scannumber {
	my $self = shift;
	return $self->{'title'};
}


sub Tr {
	my $self = shift;
	return $self->{'rtsec'};
}

sub precursor_mz {
	my $self = shift;
	return $self->{'mz'};
}

sub charge {
	my $self = shift;
	return $self->{'z'};
}

sub FT_charge {
	my $self = shift;
	return $self->{'z'};
}


sub title {
	my $self = shift;
	return $self->{'title'};
} 

sub get_peaklist{
my $self=shift;
return $self->{'peaks'};
}

sub Mr {
	my $self   = shift;
	my $charge = $self->charge;
	return $self->precursor_mz * $charge - $charge * $Hatom;
}




1;
