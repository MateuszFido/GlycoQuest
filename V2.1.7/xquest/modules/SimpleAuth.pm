package SimpleAuth;
use strict;
use warnings;
#---------------------------------------------------------------------------
# Module: SimpleAuth.pm
# Author(s): Thomas Walzthoeni
# Description: Simple module for user authentication.
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
use Common;
use File::Spec;

sub new
{
	my $class = shift;
	my $id = shift;   ### ID is the path to the folder and can be wathomas/test or test only
	
	my $self  = {};
	# Bless
	bless( $self, $class );

#---------------------------------------------------------------------------
#  PARAMETERS
#---------------------------------------------------------------------------
## standard log file, if you use another then redefine with the set_logfile method
my $logfile = File::Spec->catfile( getcwd, "logs","error.log" );
$self -> {'logfile'} = $logfile;

$self->{'status'}=1; ## 0 not authenticated, 1 auth ok, -1 some error
#---------------------------------------------------------------------------
#  Check if the Rootfolder is set
#---------------------------------------------------------------------------
### Get the rootfoldername
my @splitet = split ("/", $id);
my $rootfoldername = $splitet[0];
#print $rootfoldername;

unless ($rootfoldername){
$self->{'status'}=-1;
$self->{'error'}="No root folder set for authentification.<br>";
return $self;
}

#---------------------------------------------------------------------------
#  Get the Username
#---------------------------------------------------------------------------
my $currentusername=$ENV{'REMOTE_USER'};

unless ($currentusername){
$self->{'status'}=-1;
$self->{'error'}="Error: User is not authenticated.<br>";
$self->{'user'}="Unknown";
}else{
$self->{'user'}=$currentusername;
}

#---------------------------------------------------------------------------
#  Define the Protected folders here
#---------------------------------------------------------------------------
my $protectedfolders = {};
## Add a user or more users to a certain folder e.g.
#push @{$protectedfolders->{'foldername'}}, "username";

#---------------------------------------------------------------------------
#  Check if the folder is protected
#---------------------------------------------------------------------------
### Check if the user has access to the rootfolder
if ($protectedfolders->{$rootfoldername}){
$self->{'status'}=0;
$self->{'error'}="Error: Folder is not availaible<br>";

## check if the user is in the list
foreach my $user (@{$protectedfolders->{$rootfoldername} }){
if ($user eq $currentusername){
$self->{'status'}=1;
$self->{'error'}=0;
}
}
}
### initialize the authentication
### if not auth then return a bad status and an error msg.
return $self;
}

sub set_logfile{
my $self=shift;
my $logfile=shift;
$self -> {'logfile'} = $logfile;
}


sub get_status{
my $self=shift;
return $self->{'status'};
}

sub get_error{
my $self=shift;
return $self->{'error'};
}

sub writelog{
my $self =shift;
my $msg =shift;
my $logfile = $self -> {'logfile'};
$msg= Common::prettytime().": ". $msg;
Common::save_to_file($logfile,$msg,1);
return;
}; 


1;