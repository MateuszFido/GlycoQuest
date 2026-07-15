#!C:/Perl/bin/perl.exe
use strict;

#---------------------------------------------------------------------------
# dispdef.cgi
# A software/script to display definition files.
# Author: Thomas Walzthoeni
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
use strict;
use CGI qw/:standard :html3/;
use CGI::Carp qw(fatalsToBrowser);
use CGI::FastTemplate;
use File::Basename;
use File::Spec;
use File::Path;
use Data::Dumper;
use Cwd;
use File::stat;
use Time::localtime;
use Time::Piece;
use Storable;

# include the directory for the xquest modules ###########
use File::Spec::Functions qw(rel2abs);
use File::Basename;
use FindBin;
use lib "$FindBin::Bin/../modules";
##########################################################
use Environment;
use SimpleAuth;

#---------------------------------------------------------------------------
#  Version info
#---------------------------------------------------------------------------
my $version = "1.0";

#---------------------------------------------------------------------------
#  Create the cgi object and initialize parameters
#---------------------------------------------------------------------------
my $form = new CGI;

#---------------------------------------------------------------------------
#  Read the webparams: defined in the web.config file, path must be defined in Environment.pm
#---------------------------------------------------------------------------
my $env       = Environment->new();
my $webconfig = $env->get_path('web.config');
my $WEBPARAMS = readtables($webconfig);
my $id        = $form->param('id');             ## You can use this if you want to pass the param by URL, just uncomment
## Get the path
# The resutdirectory base e.g. /home/xquest/
my $resultdirbase = $WEBPARAMS->{'resultdirbase'};    # defined in the webparams
## Concatenate the paths, the path is then e.g. /home/xquest/username/xquestresults
my $file = File::Spec->catfile( $resultdirbase, $id );

#---------------------------------------------------------------------------
#  Print the HTML header
#---------------------------------------------------------------------------
print_header();
print "<h1>xQuest/xProphet deffile viewer</h1><hr>";
## Read the file
#print $file."<br>";
if ( -e $file )
{
	print_def($file);
} else
{
	print '<div class="error">' . "Error: cannot find $file<br>" . "</div>";
}
print $form->end_form;
print_footer();
print $form->end_html;

sub print_def
{
	my $deffile = shift;
	my %PARAMS;
	print "<div>";
	open DEF, "<$deffile" or die "cannot open table $deffile $!";
	while ( my $line = <DEF> )
	{
		print $line. "<br>";
	}
	print "</div>";

	#print "Path to file: $file<br>";
}

sub print_header
{
	print $form->header('text/html');
	print $form->start_html(
							 -title   => 'xQuest/xProphet deffile viewer',
							 -author  => 'Thomas Walzthoeni',
							 -expires => '+30s',
							 -meta    => { 'keywords' => 'MS cross-link xQuest xProhet', },
							 -style   => { 'src' => $WEBPARAMS->{'css_stylesheet'}, },
	);
}

sub readtables
{
	my $webparam = shift;
	my %WEBPARAMS;
	open WEBCONFIG, "<$webparam" or warn "could not open web config file $! ignoring";
	while (<WEBCONFIG>)
	{
		chomp;
		my @keyvalue = split /::/;
		$WEBPARAMS{ $keyvalue[0] } = $keyvalue[1];
	}
	return \%WEBPARAMS;
}

sub print_footer
{
	my $htmltable;
	$htmltable .= "<div style=\"position:relative;border:0px solid #00ff00;bottom:-25px; text-align: right\" >";
	$htmltable .= "xQuest/xProphet deffile viewer, version $version, Thomas Walzthoeni / ETH Zurich";
	$htmltable .= '</div>';
	print $htmltable;
}
