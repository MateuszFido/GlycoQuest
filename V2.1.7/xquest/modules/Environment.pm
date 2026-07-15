package Environment;
use strict;
use warnings;

#---------------------------------------------------------------------------
# Module: Environment.pm
# Author(s): Thomas Walzthoeni
# Description: Module to set xQuest/xProphet specific environment variables.
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
use Cwd;

sub new
{
	my $class = shift;
	my $self  = {};
	my $cwd   = getcwd;
	my %machines;

	# Define paths to xQuest on the server
	# Machines/Servers Hash: for Windows machines use the name that is stored in $ENV{'LOGONSERVER'}
	# for Linux machines use the hostname that you get with the command >hostname
	#---------------------------------------------------------------------------
	# To define a server please edit the hostname as the key e.g. ('xquest-desktop')
	# and as value as a useful name.
	#---------------------------------------------------------------------------
	$machines{'xqxp-desktop'} = "xquestvm";
        $machines{'brutus'}      = "brutus";
		$machines{'euler'}      = "euler";

	#---------------------------------------------------------------------------
	$self->{'machines'} = \%machines;
	my %serverpaths;

	#---------------------------------------------------------------------------
	#  Edit the paths
	#---------------------------------------------------------------------------
	# Edit the paths for the specific server, as key use the name of the server as defined above
	# Paths for the xQuest VM
	$serverpaths{'xquestvm'}{'xquest_stable'} = "\/home\/xquest\/xquest";
	$serverpaths{'xquestvm'}{'web.config'}    = "\/home\/xquest\/xquest\/conf\/web.config";
	$serverpaths{'xquestvm'}{'mass.def'} = "\/home\/xquest\/xquest\/deffiles\/mass_table.def";
	
#    $serverpaths{'brutus'}{'xquest_stable'} = "\/cluster\/apps\/xquest\/V2_1_4\/xquest";
    $serverpaths{'brutus'}{'web.config'} = "\/cluster\/apps\/xquest\/V2_1_4\/xquest\/conf\/web.config";
    $serverpaths{'brutus'}{'mass.def'} = "\/cluster\/apps\/xquest\/V2_1_4\/xquest\/deffiles\/mass_table.def";
    $serverpaths{'brutus'}{'xquest_stable'} = "\/IMSB\/ra\/mewing\/xquestWorkingVersion\/V2_1_4\/xquest";

	$serverpaths{'euler'}{'web.config'} = "\/cluster\/apps\/imsb\/xquest/\/V2.1.5\/xquest\/conf\/web.config"; 
    $serverpaths{'euler'}{'mass.def'} = "\/cluster\/apps\/imsb\/xquest\/V2.1.5\/xquest\/deffiles\/mass_table.def";
    #$serverpaths{'euler'}{'xquest_stable'} = "\/nfs\/nas21.ethz.ch\/nas\/fs2102\/biol_ibt_usr_s1\/mewing\/xquestWorkingVersion\/V2_1_4/xquest";
	$serverpaths{'euler'}{'xquest_stable'} = "\/cluster\/apps\/imsb\/xquest/\/V2.1.5\/xquest";


	#---------------------------------------------------------------------------
	#  Note: If you work on a cluster where the hostnames are different for the
	#        individual nodes, then you can use a regex to define the environment,
	# 		 an example is listed at line 107
	#---------------------------------------------------------------------------
	$self->{'paths'} = \%serverpaths;
### Debug vars
	$self->{'error'} = undef;
### Server var: the name of the machine
	$self->{'server'} = undef;

	# Bless
	bless( $self, $class );
### Call these functions to initialize $self->{'server'}
	unless ( $self->_set_environment() == 1 )
	{
		my $test = Dumper( \%ENV );
		confess "There was an error:", "Errormsg: $test;", $self->{'error'};
	}

	# Return the reference
	return $self;
}

sub _set_environment
{
	my $self            = shift;

	if ( my $root = $ENV{'XQUEST_ROOT'} ) {
		$root =~ s/\/$//;
		if ( -d $root ) {
			$self->{'paths'}{'glycoquest_local'} = {
				'xquest_stable' => $root,
				'web.config'    => "$root/conf/web.config",
				'mass.def'      => "$root/deffiles/mass_table.def",
			};
			$self->{'machines'}{'glycoquest_local'} = 'glycoquest_local';
			$self->{'server'} = 'glycoquest_local';
			return 1;
		}
	}

	my $operatingsystem = $^O;
	my $server;
	my $name;
	if ( $operatingsystem eq "MSWin32" )
	{
## get the Computer name (this var is set for pl scripts)
		$name = $ENV{'LOGONSERVER'};
		unless ($name)
		{
## use the server_admin var for cgi scripts
			$name = $ENV{'SERVER_ADMIN'};
		}
		chomp($name);
		if ( ( $self->{'machines'}->{$name} ) )
		{
			$server = $self->{'machines'}->{$name};
		}
	} elsif ( $operatingsystem eq "linux" || $operatingsystem eq "darwin" )
	{
## get the network node name
		$name = `hostname`;
		unless ($name)
		{
## use the server_admin var for cgi scripts
			$name = $ENV{'SERVER_ADMIN'};
		}
		chomp($name);
		if ( $name =~ m/hpc-net.ethz.ch/ || $name =~ m/brutus/ )
		{
			$name = "brutus";
			$server = $self->{'machines'}->{$name};
		}
		elsif ( $name =~ m/e\d{4}/ || $name =~ m/euler/ || $name =~ m/eu-ms/ || $name =~ m/eu-c7/ || $name =~ m/eu-login/ )
		{
			$name = "euler";
			$server = $self->{'machines'}->{$name};
		}
		elsif ( ( $self->{'machines'}->{$name} ) )
		{
			$server = $self->{'machines'}->{$name};
		}
		else
		{
			# This (and the above corresponding elsif is where you will make changes if the cluster nocdes change name)
			$name = "euler";
			$server = $self->{'machines'}->{$name};
		}
	} else
	{
		$self->{'error'} = "Unknown Operating System: " . $name;
		return 0;
	}
## check if a server is set
	unless ($server)
	{
		$self->{'error'} = "No server set: " . $name;
		return 0;
	}
	$self->{'server'} = $server;
	return 1;
}

sub get_env
{
	my $self = shift;
	return $self->{'server'};
}
## pass array content
sub get_path
{
	my $self   = shift;
	my $key    = shift;
	my $server = $self->{'server'};
	my $path   = $self->{'paths'}->{$server}->{$key};
	unless ($path)
	{
		$self->{'error'} = "Error: No path found in paths hash: server: $server, Key:$key";
		die confess( $self->{'error'} );
	}
	return $path;
}
1;
