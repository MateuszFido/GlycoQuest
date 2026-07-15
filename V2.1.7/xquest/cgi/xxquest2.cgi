#!C:/Perl/bin/perl.exe
use strict;
#---------------------------------------------------------------------------
# xxquest2.cgi
# A software/script to display xQuest/xProphet results.
# Author: Thomas Walzthoeni based on orginal version of Oliver Rinner
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
use CGI qw/:standard :html3/;
use CGI::Carp qw(fatalsToBrowser);
use CGI::FastTemplate;
use Socket;
use File::Basename;
use File::Spec;
use File::Path;

#use Mail::Sender;
use XML::Parser;
use Data::Dumper;
use XML::TreeBuilder;
use Cwd;
use Archive::Tar;
use Bio::Index::Fasta;
### Packages for database use
use MLDBM qw(DB_File Storable);
use DB_File;
use Fcntl;
use Data::Dumper;
use Storable;
use CGI::Session;
use HTML::PageIndex;
# Include the directory for the xquest modules ###########
use File::Spec::Functions qw(rel2abs);
use File::Basename;
use FindBin;
use lib "$FindBin::Bin/../modules";
##########################################################
use Environment;
use SimpleAuth;
my $env       = Environment->new();
my $webconfig = $env->get_path('web.config');

#---------------------------------------------------------------------------
#  Debugging of Params
#---------------------------------------------------------------------------
my $debug = 0;    # to display param fields for debugging
my $debugmsg;     # write sth into this var, is printed at the beginning of the table

#---------------------------------------------------------------------------
#  Create the cgi object and initialize parmeters
#---------------------------------------------------------------------------
my $form              = new CGI;
my $OS                = $^O;
my $myself            = self_url();
my $time              = localtime();
my $id                = param('id');
my $firstrun          = $form->param('firstrun');
my $xmlfilename       = $form->param('resultxml');
my $peptidelistname   = $form->param('peptidelistname');
my $urlbase           = $form->param('urlbase');
my $maxrank           = $form->param('maxrank');
my $WEBPARAMS         = readwebparams($webconfig);
my $resultdirbase     = $WEBPARAMS->{'resultdirbase'};
my $resulturlbase     = $WEBPARAMS->{'resulturlbase'};
my $xlinkionparser    = $WEBPARAMS->{'xlinkionparser'};
my $monolinkionparser = $WEBPARAMS->{'monolinkionparser'};
my $xxquest           = basename( $WEBPARAMS->{'xmlparser'} );
my $version           = $form->param('v');
my $viewerversion     = "2.2.3";
### Can be used as an Urlparam to create a newindex
### Is not implemented in the form
my $newindex   = $form->param('newindex');
my $newdbindex = $form->param('newdbidx');

unless ($version)
{
	$version = "current";
}
my $hashfilename;
my $resultdirectory;

#---------------------------------------------------------------------------
#  Prints the html header
#---------------------------------------------------------------------------
# print_header();
# setparams();
#---------------------------------------------------------------------------
# Check if the id was posted (id is the url to the directory)
# and set some paths
#---------------------------------------------------------------------------
## Define the file names for the result files
if ($id)
{
	if ( ($xmlfilename) )
	{
		$hashfilename = File::Spec->catfile( $resultdirbase, $id, "$xmlfilename.hash" );
		$xmlfilename  = File::Spec->catfile( $resultdirbase, $id, $xmlfilename );
		$resultdirectory = File::Spec->catfile( $resultdirbase, $id );
	} else
	{
		$hashfilename = File::Spec->catfile( $resultdirbase, $id, 'xquest.xml.hash' );
		$xmlfilename  = File::Spec->catfile( $resultdirbase, $id, 'xquest.xml' );
		$resultdirectory = File::Spec->catfile( $resultdirbase, $id );
	}
	$peptidelistname = File::Spec->catfile( $resultdirbase, $id, 'database_peptides.txt' );
	$urlbase = join "\/", $resulturlbase, $id;
} else
{
	print "Search id was not posted <br>";
	exit;
}
my $xmlurl             = join "\/", $resulturlbase, $id, 'xquest.xml';
my $statusurl          = join "\/", $resulturlbase, $id, 'xquest.stat.done';
my $peptidelisturl     = join "\/", $resulturlbase, $id, 'database_peptides.txt';
my $local_xquestdefurl = join "\/", $resulturlbase, $id, 'xquest.def';
my $statusfile = File::Spec->catfile( $resultdirbase, $id, 'xquest.stat' );
my $xquestdef  = File::Spec->catfile( $resultdirbase, $id, 'xquest.def' );
my $xqprophdef = File::Spec->catfile( $resultdirbase, $id, 'xproph.def' );
my $folder     = File::Spec->catfile( $resultdirbase, $id );

#---------------------------------------------------------------------------
#  Error msg variable
#---------------------------------------------------------------------------
my $errormsg;
my $htmlheader;
my $warnings;
$htmlheader .= '<h1>xQuest/xProphet results viewer ' . basename($id) . '</h1><hr>' . "\n";
my $htmlerror = 0;

#---------------------------------------------------------------------------
#  ### FOLDER PROTECTION ###
#---------------------------------------------------------------------------
my $noauth  = 1;
my $authobj = SimpleAuth->new($id);
my $auth    = $authobj->get_status;
### CHECK FOR AUTH STATUS
unless ($noauth)
{
	unless ( $auth == 1 )
	{
		$errormsg .= $authobj->get_error;
		print_header();
		print $htmlheader;
		print $errormsg;
		my $msg = "$0: Access to $id denied.\n";
		$authobj->writelog($msg);
		exit 0;
	}
}

#---------------------------------------------------------------------------
#  Print html header
#---------------------------------------------------------------------------
### check if the resultdirectory exists
unless ( -e $resultdirectory )
{
### check if a tgz exists
	my $tgz = File::Spec->catfile( $resultdirbase, $id . ".tgz" );
	if ( -e $tgz )
	{
### untar the tgz!
		my $tar = Archive::Tar->new;
### Switch the tars working directory!
		$tar->setcwd($resultdirbase);
		$tar->read($tgz);
### Extract
		my $list = $tar->extract();
	}
}
unless ( -e $xquestdef )
{
	$errormsg .= "<div class\"error\">Software error: xQuest $xquestdef definition file was not found!</div><br>";
}
unless ( -e $xmlfilename )
{
	$errormsg .= "<div class=\"error\">Software error: xQuest XML file $xmlfilename was not found.</div><br>";
}
if ($errormsg)
{
	print_header();
	print $htmlheader;
	print $errormsg;
	print $form->end_form;
	print $form->end_html;
	exit 0;
}
my $PARAMS = read_def($xquestdef);

#---------------------------------------------------------------------------
#  Initialize a Session
#---------------------------------------------------------------------------
my $session;
my $cookie;
my $sid;
my $urlparam;
## Set the session name
CGI::Session->name("sid");
## Look if the client sent a cookie
my $cookieanswer = $form->cookie('sid');
## try to load the session (loads if one is found but creates a new),
## pass the cgi_obj, this will load session also from the cookie id
#  $session = CGI::Session->load( "driver:File", $form, { Directory => "E:/Apache2/tmp" } ) or die CGI::Session->errstr;
$session = CGI::Session->load("driver:File") or die CGI::Session->errstr;
## Check if the session is empty, if yes create a new one
if ( $session->is_empty )
{
## Create a new session
	#$session = new CGI::Session( "driver:File", $form, { Directory => "E:/Apache2/tmp" } );
	$session = new CGI::Session( "driver:File", $form );
## Create a cookie
	$cookie = $form->cookie( -name => $session->name, -value => $session->id );
}
## Set set sid variable
$sid = $session->id();
if ("$cookieanswer")
{

	#print "Cookie found: ".$cookieanswer."<br>";
} else
{
	$urlparam = "sid=$sid";

	#print "No cookie found! Generate all links with GET param $urlparam.<br>";
}

#---------------------------------------------------------------------------
#  Save parameters of cgi object to the session and define defaults
#  if they were not set
#---------------------------------------------------------------------------
$session->save_param($form);
my $maxerrorfrom;
my $maxerrorto;
my $unique_ids;
my $deltas;
my $mions;
my $minscore;
## Read the xproph.def file and set the parameters accordingly
if ( -e $xqprophdef )
{
	read_params( $xqprophdef, $PARAMS );
	$maxerrorfrom = $PARAMS->{'minborder'};
	$maxerrorto   = $PARAMS->{'maxborder'};
	$unique_ids   = $PARAMS->{'uniquexl'};
	$deltas       = $PARAMS->{'mindeltas'};
	$mions        = $PARAMS->{'minionsmatched'};
	$minscore     = $PARAMS->{'minscore'};
	unless ( $session->param('submit_form') || $session->param('saveselected') || $session->param('page') )
	{
		$debugmsg = '<div id="highlight">Filters were loaded from xproph.def</div>';
	}
}
## Check if the user has set other parameters, override xproph parameters
if ( $session->param("maxerrorfrom") )
{
	$maxerrorfrom = $session->param("maxerrorfrom");
}
if ( $session->param("maxerrorto") )
{
	$maxerrorto = $session->param("maxerrorto");
}
unless ( $session->param('unique_ids') )
{
## Check if form was submitted
	if ( $session->param('submit_form') )
	{
		$unique_ids = undef;
	}
} else
{
	$unique_ids = 1;
}
unless ( $session->param('deltas') )
{
## check if form was submitted
	if ( $session->param('submit_form') )
	{
		$deltas = undef;
	}
} else
{
	$deltas = $session->param('deltas');
}
unless ( $session->param('mions') )
{
## Check if form was submitted
	if ( $session->param('submit_form') )
	{
		$mions = undef;
	}
} else
{
	$mions = $session->param('mions');
}
unless ( $session->param('minscore') )
{
## Check if form was submitted
	if ( $session->param('submit_form') )
	{
		$minscore = undef;
	}
} else
{
	$minscore = $session->param('minscore');
}
my $xprophet_flagged;
unless ( $session->param('xprophet_flagged') )
{
## check if form was submitted
	if ( $session->param('submit_form') )
	{
		$xprophet_flagged = undef;
	}
} else
{
	$xprophet_flagged = 1;
}
my $uniquerestraints;
unless ( $session->param('uniquerestraints') )
{
## check if form was submitted
	if ( $session->param('submit_form') )
	{
		$uniquerestraints = undef;
	}
} else
{
	$uniquerestraints = 1;
}

#unless ($mions) { $mions = 0 }
## If parameters are not set use params from the xquest def
unless ($maxerrorfrom) { $maxerrorfrom = -$PARAMS->{'ms1tolerance'} }
unless ($maxerrorto)   { $maxerrorto   = $PARAMS->{'ms1tolerance'} }
### Initialize the parameters if they were not set
my $reporttype = $session->param('reporttype');
unless ($reporttype) { $reporttype = 1 }
my $disptype = $session->param('disptype');
unless ($disptype) { $disptype = 1 }
my $dispmax = $session->param("dispmax");
unless ($dispmax) { $dispmax = 100 }
my $dispnhits = $session->param("dispnhits");
unless ($dispnhits) { $dispnhits = 1 }
my $filter_error_on_all_hits = $session->param("filter_error_on_all_hits");
unless ($filter_error_on_all_hits) { $filter_error_on_all_hits = undef }
my $maxscore = $session->param('maxscore');
unless ($maxscore) { $maxscore = 0 }
### Create new index (initialize the hash again)
### Check if the param was set
my $make_new_index;

if ( $session->param('make_new_index') )
{
	$make_new_index = $session->param('make_new_index');
}
unless ($make_new_index) { $make_new_index = undef }
my $deltaAAmin = $session->param('deltaAAmin');
unless ($deltaAAmin) { $deltaAAmin = 0 }
my $deltaAAmax = $session->param('deltaAAmax');
unless ($deltaAAmax) { $deltaAAmax = 0 }
my $nseenmin = $session->param('nseenmin');
unless ($nseenmin) { $nseenmin = 0 }
my $nseenmax = $session->param('nseenmax');
unless ($nseenmax) { $nseenmax = 0 }
my $sequencefilter = $session->param('sequencefilter');
unless ($sequencefilter) { $sequencefilter = 0 }
my $annotationfilter = $session->param('annotationfilter');
unless ($annotationfilter) { $annotationfilter = 0 }
my $spectrum = $session->param('spectrum');
unless ($spectrum) { $spectrum = 0 }
my $fdrcutoff = $session->param('fdrcutoff');
unless ($fdrcutoff) { $fdrcutoff = 0 }
my $pagedisp = $session->param('page');
unless ($pagedisp) { $pagedisp = 1 }
my $selectallhits = $session->param('selectallhits');
unless ($selectallhits) { $selectallhits = 0 }
my $unselectallhits = $session->param('unselectallhits');
unless ($unselectallhits) { $unselectallhits = 0 }

#---------------------------------------------------------------------------
#  Verbose
#---------------------------------------------------------------------------
my $verbose = 0;

#---------------------------------------------------------------------------
#  Index the fasta database
#---------------------------------------------------------------------------
my $fasta = basename( $PARAMS->{'database'} );

#$debugmsg .= "Basename db:$fasta<br>";
my $database = File::Spec->catfile( $resultdirbase, $id, $fasta );
my $databasedecoy;
if ( $PARAMS->{'database_dc'} )
{
	my $fastadc = basename( $PARAMS->{'database_dc'} );
	$databasedecoy = File::Spec->catfile( $resultdirbase, $id, $fastadc );
	unless ( -e $databasedecoy )
	{

		#$debugmsg .= '<div class="error">Warning: database ' . $databasedecoy . ' does not exist.</div><p>';
	}
}
my $idx;

unless ( -e $database )
{
	print_header();
	$errormsg .= '<div class="error">Software error: database ' . $database . ' does not exist.</div><p>';
}

if ($errormsg)
{
	print_header();
	print $htmlheader;
	print $errormsg;
	print $form->end_form;
	print $form->end_html;
	exit 0;
}
## Unlink the index file if newindex is selected
## Create a version file if not defined
my $index_version_file = File::Spec->catfile( $resultdirbase, $id, "fasta_index.version" );
my $indexversion = Bio::Index::Fasta->_version();
unless ($indexversion)
{
	$indexversion = "undefined";
}
my $index_version;
unless ( -e $index_version_file )
{
## write file
	write_version_to_file( $index_version_file, $indexversion );
} else
{
	$index_version = get_version_from_file($index_version_file);
}
## Compare
unless ( $index_version eq $indexversion )
{
	$newdbindex = 1;
} else
{
}
if ( $make_new_index || $newdbindex )
{
	unlink( $database . ".idx" );
}
unless ( -e $database . ".idx" )
{
	$newdbindex = 1;
}
if ( $make_new_index || $newdbindex )
{
	($verbose) && print "Creating index <br>";
### generate an index file
	$idx = Bio::Index::Fasta->new( -filename => $database . ".idx", -write_flag => 1 );
### parse the databases
	$idx->make_index($database);
### parse the decoy db if availaible
	if ( -e $databasedecoy )
	{
		$idx->make_index($databasedecoy);
	}
} else
{
	## read from indexfile
	($verbose) && print "Database index exists! <br>";
	$idx = Bio::Index::Fasta->new( $database . ".idx" );
}

#---------------------------------------------------------------------------
#  Parsing and storing or retrieving the resulthash from the xml file
#---------------------------------------------------------------------------
my $dbfilesemaphore;

#my %resultshash;
my $reshashref = {};
($verbose) && print "Generated hashref $reshashref<br>";
if ( ( -e $hashfilename && -r $hashfilename ) && ( !$make_new_index ) )
{
	$dbfilesemaphore = 1;
## reads the Hash from the DB
	$reshashref = retrieve($hashfilename);
	($verbose) && print "Result hash found $reshashref<br>";
} else
{
	$dbfilesemaphore = 0;
### Parse the file and store the hash to the db file
	parseXML( $xmlfilename, $reshashref, $idx, $resultdirbase, $id );
	store $reshashref, $hashfilename;
	($verbose) && print "Created new index, parsed xml and stored to hash<br>";
	## undef the makenewindex afterwards
	($verbose) && print $make_new_index;
	### reset the form parameter
	$form->param( 'make_new_index', 0 );
}

sub print_form
{

	#---------------------------------------------------------------------------
	#  GENERATE THE FORM FOR FILTERING AND SORTING
	#---------------------------------------------------------------------------
	print start_multipart_form;
	print '<TABLE BORDER=0 CELLSPACING=1 CELLPADDING=3>';
	print '<TR BGCOLOR=#EEEEFF NOWRAP>';
	print '<th style="color: #0000FF; font-size: 12pt" align="left" colspan="2">General settings</td>';
	print '<th style="color: #0000FF; font-size: 12pt" align="left" colspan="6">Filter settings</td>';
	print '</TR>';

	#---------------------------------------------------------------------------
	#  ROW 1
	#---------------------------------------------------------------------------
	print '<TR BGCOLOR=#EEEEFF NOWRAP>';
	print '<TD>';
	my @values = ( 1, 2, 4, 5 );
	my %labels = ( 1 => 'Html Table', 2 => 'TSV file', 3 => "idXML", 4 => "xTract csv", 5 => "IL csv" );
	print "Select type of report:</TD><TD> ",
	  $form->popup_menu(
						 -name    => 'reporttype',
						 -values  => \@values,
						 -labels  => \%labels,
						 -default => $reporttype,
	  );
	print $form->br;
	print '</TD>';
	print '<TD>';
	@values = ( 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22 );

	#@values = ( 1, 2, 3, 4, 9,11,12);
	my %typelabels = (
					   1  => 'all',
					   2  => 'cross-links(all)',
					   3  => 'mono-links (all)',
					   4  => 'loop-links (all)',
					   5  => 'decoy mono-links',
					   6  => 'target mono-links',
					   7  => 'decoy loop-links',
					   8  => 'target loop-links',
					   9  => 'inter-protein cross-links (all)',
					   10 => 'inter-protein cross-links (target)',
					   11 => 'intra-protein cross-links (target)',
					   12 => 'intra / inter xl (target)',
					   13 => 'decoy inter-protein xl',
					   14 => 'decoy intra-protein xl',
					   15 => 'decoy intra/inter xl',
					   16 => 'selected hits',
					   17 => 'target hits',
					   18 => 'decoy hits',
					   19 => 'target mono- and loop-links',
					   20 => 'decoy mono- and loop-links',
						21 => 'target cross-links',
					   22 => 'decoy cross-links',					   
					   
					   
	);
	print "Filter by type (top hit):</TD><TD> ",
	  $form->popup_menu(
						 -name    => 'disptype',
						 -values  => \@values,
						 -labels  => \%typelabels,
						 -default => $disptype,
	  );
	print '</TD>';

	#print $form->br;
	print '<TD>';
	print "Filter by min # seen >:</TD><TD> ",
	  $form->textfield(
						-name      => 'nseenmin',
						-value     => $nseenmin,
						-size      => 2,
						-maxlength => 5
	  );
	print " and <  ",
	  $form->textfield(
						-name      => 'nseenmax',
						-value     => $nseenmax,
						-size      => 2,
						-maxlength => 5
	  );
	print $form->br;
	print '</TD>';
	print '<TD>';
	print "Filter by min ions matched>:</TD><TD> ",
	  $form->textfield(
						-name      => 'mions',
						-value     => $mions,
						-size      => 2,
						-maxlength => 3
	  );
	print $form->br;
	print '</TD>';
	print '</TR>';

	#---------------------------------------------------------------------------
	#  ROW 2
	#---------------------------------------------------------------------------
	print '<TR BGCOLOR=#EEEEFF NOWRAP>';
	print '<TD>';
	@values = ( 'all', 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 );
	print "Show n ranks per spectrum:</TD><TD> ",
	  $form->popup_menu(
						 -name    => 'dispnhits',
						 -values  => \@values,
						 -default => $dispnhits,
	  );
	print '</TD>';
	print '<TD>';
	print "Filter hits by max $PARAMS->{'tolerancemeasure'} (Range):</TD>";
	print "<TD> From ",
	  $form->textfield(
						-name      => 'maxerrorfrom',
						-value     => $maxerrorfrom,
						-size      => 2,
						-maxlength => 5
	  );
	print " to ",
	  $form->textfield(
						-name      => 'maxerrorto',
						-value     => $maxerrorto,
						-size      => 2,
						-maxlength => 5
	  );
	print "&nbsp $PARAMS->{'tolerancemeasure'} ";
	print "All hits:",
	  $form->checkbox(
					   -name    => 'filter_error_on_all_hits',
					   -checked => $filter_error_on_all_hits,
					   -value   => 1,
					   -label   => "",
	  );
	print '</TD>';
	print '<TD>';
	print "Filter by sequence:</TD>";
	print "<TD>",
	  $form->textfield(
						-name      => 'sequencefilter',
						-value     => $sequencefilter,
						-size      => 25,
						-maxlength => 100,
	  );
	print "</TD>";
	print '<TD>';
	print "Select all hits:</TD>";
	print "<TD>",
	  $form->checkbox(
					   -name    => 'selectallhits',
					   -checked => undef,
					   -value   => 1,
					   -label   => "",
	  );
	print '</TD>';
	print '</TR>';

	#---------------------------------------------------------------------------
	#  ROW 3
	#---------------------------------------------------------------------------
	print '<TR BGCOLOR=#EEEEFF NOWRAP>';
	print '<TD>';
	@values = ( 25, 50, 75, 100, 500 );
	print "Number of hits per page:</TD><TD> ",
	  $form->popup_menu(
						 -name    => 'dispmax',
						 -values  => \@values,
						 -default => $dispmax,
	  );
	print '</TD>';
	print '<TD>';
	print "Filter by unique ids (top hit):</TD><TD> ",
	  $form->checkbox(
					   -name    => 'unique_ids',
					   -checked => $unique_ids,
					   -value   => 1,
					   -label   => "",
	  );
	print '</TD>';
	print '<TD>';
	print "Filter by annotation:</TD>";
	print "<TD>",
	  $form->textfield(
						-name      => 'annotationfilter',
						-value     => $annotationfilter,
						-size      => 25,
						-maxlength => 100,
	  );
	print "</TD>";
	print '<TD>';
	print "Unselect all hits:</TD>";
	print "<TD>",
	  $form->checkbox(
					   -name    => 'unselectallhits',
					   -checked => undef,
					   -value   => 1,
					   -label   => "",
	  );
	print '</TD>';
	print '</TR>';

	#---------------------------------------------------------------------------
	#  Row 4
	#---------------------------------------------------------------------------
	print '<TR BGCOLOR=#EEEEFF NOWRAP>';
	print '<TD>';
	print "Create new index:</TD><TD> ",
	  $form->checkbox(
					   -name    => 'make_new_index',
					   -checked => undef,
					   -value   => 1,
					   -label   => "",
	  );
	print '</TD>';
	print '<TD>';
	print "Show scores (top hit) >:</TD><TD> ",
	  $form->textfield(
						-name      => 'minscore',
						-value     => $minscore,
						-size      => 2,
						-maxlength => 5
	  );
	print " and < ",
	  $form->textfield(
						-name      => 'maxscore',
						-value     => $maxscore,
						-size      => 2,
						-maxlength => 5
	  );
	print '</TD>';
	print '<TD>';
	print "Filter by deltaS (top hit) <:</TD><TD> ",
	  $form->textfield(
						-name      => 'deltas',
						-value     => $deltas,
						-size      => 2,
						-maxlength => 5
	  );
	print '</TD>';
	print '<TD>';
	print "xProphet flagged:</TD>";
	print "<TD>",
	  $form->checkbox(
					   -name    => 'xprophet_flagged',
					   -checked => $xprophet_flagged,
					   -value   => 1,
					   -label   => "",
	  );
	print '</TD>';
	print '</TR>';

	#---------------------------------------------------------------------------
	#  Row 5
	#---------------------------------------------------------------------------
	print '<TR BGCOLOR=#EEEEFF NOWRAP>';
	print '<TD>';
	print "Refresh:</TD><TD> ",
	  $form->submit( -name  => 'submit_form',
					 -value => 'Update', );
	print $form->submit( -name  => 'saveselected',
						 -value => 'Save Selected', );
	print '</TD>';
	print '<TD>';
	print "Filter by &Delta;AA (top hit, Range):</TD><TD>From  ",
	  $form->textfield(
						-name      => 'deltaAAmin',
						-value     => $deltaAAmin,
						-size      => 2,
						-maxlength => 5
	  );
	print " to ",
	  $form->textfield(
						-name      => 'deltaAAmax',
						-value     => $deltaAAmax,
						-size      => 2,
						-maxlength => 5
	  );
	print '</TD>';
	print '<TD>';
	print "Filter by FDR <:</TD><TD>",
	  $form->textfield(
						-name      => 'fdrcutoff',
						-value     => $fdrcutoff,
						-size      => 2,
						-maxlength => 5
	  );
	print '</TD>';
	print "<TD/>Filter by unique restraints:</TD>";
	print "<TD>",
	  $form->checkbox(
					   -name    => 'uniquerestraints',
					   -checked => $uniquerestraints,
					   -value   => 1,
					   -label   => "",
	  );
	print '</TD>';
	print '</TR>';

	#---------------------------------------------------------------------------
	#  END OF TABLE
	#---------------------------------------------------------------------------
	print '</TABLE>';

	#---------------------------------------------------------------------------
	#  Hidden fields ->search id
	#---------------------------------------------------------------------------
	print $form->hidden( 'id', $id );

	#---------------------------------------------------------------------------
	#  Hidden fields ->resultxmlfilename
	#---------------------------------------------------------------------------
	print $form->hidden( 'resultxml', basename($xmlfilename) );

	#---------------------------------------------------------------------------
	#  Hidden fields ->pagedisp
	#---------------------------------------------------------------------------
	print $form->hidden( 'page', $pagedisp );

	#---------------------------------------------------------------------------
	#  END THE FORM
	#---------------------------------------------------------------------------
	#print endform; ## is now ended at the end of the html table (checkboxes validate)
}

#---------------------------------------------------------------------------
#  Sorting functions applied if html output is selected
#---------------------------------------------------------------------------
#print "Resulthash: $reshashref<br>";
my $ids               = get_all_ids($reshashref);
my $numberoftotalhits = scalar(@$ids);
my $sorted;
my $debugarray;
if ( ( $reporttype == 1 ) || ( $reporttype == 2 ) || ( $reporttype == 3 ) || ( $reporttype == 4 ) || ( $reporttype == 5 ) )
{
	($verbose) && print "Total number of spectra: $numberoftotalhits<br>";
### sort decending by score
	$sorted = sort_hash_desc( $reshashref, $ids );
## Filter by type 1 is show all
	unless ( ( $disptype == 1 ) )
	{
		my %types = (
					  1  => 'all',
					  2  => 'xlink',
					  3  => 'monolink',
					  4  => 'intralink',
					  5  => 'decoy monolink',
					  6  => 'monolink',
					  7  => 'decoy intralink',
					  8  => 'intralink',
					  9  => 'inter-protein xl',
					  10 => 'inter-protein xl',
					  11 => 'intra-protein xl',
					  12 => 'intra/inter xl',
					  13 => 'decoy inter-protein xl',
					  14 => 'decoy intra-protein xl',
					  15 => 'decoy intra/inter xl',
					  17 => 'decoy',
					  18 => 'decoy',
		);
		if ( $disptype > 4 )
		{

			#print $disptype;
			if ( $disptype == 9 )
			{

				# type 5 is all inter protein cross-links
				$sorted = sort_xlinks_by_type_match_string( $reshashref, $sorted, "inter-protein xl" );
			} elsif ( $disptype == 16 )
			{
				$sorted = sort_xlinks_by_selected( $reshashref, $sorted, $unique_ids );
			} elsif ( $disptype == 17 or $disptype == 18 )
			{
				if ( $disptype == 17 )
				{
					$sorted = sort_xlinks_by_type_exclude_match_string( $reshashref, $sorted, $types{$disptype} );
				}
				if ( $disptype == 18 )
				{
					$sorted = sort_xlinks_by_type_match_string( $reshashref, $sorted, $types{$disptype} );
				}
			} elsif ( $disptype == 19 or $disptype == 20 )
			{
				$sorted = sort_xlinks_by_2types_and_decoy( $reshashref, $sorted, "monolink", "intralink" ) if ( $disptype == 19 );
				$sorted = sort_xlinks_by_2types_and_decoy( $reshashref, $sorted, "monolink", "intralink", 1 ) if ( $disptype == 20 );
			} elsif( $disptype == 21 or $disptype == 22 ){
			$sorted = sort_xlinks_by_type_and_decoy( $reshashref, $sorted, "xlink" ) if ( $disptype == 21 );
			$sorted = sort_xlinks_by_type_and_decoy( $reshashref, $sorted, "xlink", 1 ) if ( $disptype == 22 );	
			}else
			{
				$sorted = sort_xlinks_by_type( $reshashref, $sorted, $types{$disptype} );
			}
		} else
		{
			$sorted = sort_by_type( $reshashref, $sorted, $types{$disptype} );
		}
	}
	
	## Filter by ppm error if filtering is selected on top hit (all hits is not selected)
	#  print "MS1 Tolerance: ", $PARAMS->{'ms1tolerance'};
	my $ppmstd = $PARAMS->{'ms1tolerance'};
	if ( ( ( -$ppmstd != $maxerrorfrom ) || ( $ppmstd != $maxerrorto ) ) && ( !$filter_error_on_all_hits ) )
	{
		$sorted = filter_by_ppm( $reshashref, $sorted, $maxerrorfrom, $maxerrorto );
	}
## Filter by deltaAA
	if ( $deltaAAmin || $deltaAAmax )
	{
		($verbose) && print "Sorting by deltaAA<br>";
		$sorted = filter_by_deltaAA( $reshashref, $sorted, $deltaAAmin, $deltaAAmax );
	}
## Filter by nseen
	if ( $nseenmin || $nseenmax )
	{
		($verbose) && print "Sorting by nseen<br>";
		$sorted = filter_by_nseen( $reshashref, $sorted, $nseenmin, $nseenmax );
	}
## Filter by sequence
	if ($sequencefilter)
	{
		($verbose) && print "Filter by sequence<br>";
		### check for allowed chars
		#unless ($useraddress =~ /^([-\@\w.]+)$/)
		if ( $sequencefilter =~ /^[\w .!?-]+$/ )
		{
			$sorted = filter_by_sequence( $reshashref, $sorted, $sequencefilter );
		}
	}
	($verbose) && print "Number of spectra after sorting: " . scalar(@$sorted);

	#$numberoftotalhits = scalar(@$sorted);
## Filter by annotation
	if ($annotationfilter)
	{
		#print "TRUE";
		($verbose) && print "Filter by annotation<br>";
		### check for allowed chars
		$sorted = filter_by_annotation( $reshashref, $sorted, $annotationfilter );
	}
## Filter by spectrum id
	if ($spectrum)
	{
		($verbose) && print "Filter by spectrum id<br>";
		### check for allowed chars
		$sorted = filter_by_spectrum( $reshashref, $sorted, $spectrum );
	}
### calculate the deltascore // updates the reshhashref, also when filter on all is used
	$reshashref = calc_delta_score( $reshashref, $sorted );
### filter by delta score
	if ($deltas)
	{
		$sorted = filter_by_deltas( $reshashref, $sorted, $deltas );
	}
	if ($fdrcutoff)
	{
		$sorted = filter_by_fdr( $reshashref, $sorted, $fdrcutoff );
	}

	# filter by mions (min matched ions per peptide)
	if ($mions)
	{
		$sorted = filter_by_mions( $reshashref, $sorted, $mions );
	}
## Filter by unique ids
	if ($unique_ids)
	{
		my $groups;
		( $sorted, $groups, $debugarray ) = sort_by_unique_ids( $reshashref, $sorted );
	}
## Filter by score
	if ( $minscore > 0 || $maxscore > 0 )
	{
		($verbose) && print "Filter by score: $minscore!";
		$sorted = filter_by_score( $reshashref, $sorted, $minscore, $maxscore );
	}
## Filter by xprophet flag
	if ($xprophet_flagged)
	{
		$sorted = filter_by_xprophet_flag( $reshashref, $sorted );
	}
## Filter by unique restraints
	if ($uniquerestraints)
	{
		my $nduplicates;
		my $nchecked;
		( $sorted, $nduplicates, $nchecked, $debugmsg ) = filter_by_unique_restraints( $reshashref, $sorted );
		$debugmsg .= "Filtered $nduplicates redundant restraints, checked $nchecked IDs.";
	}
}
($verbose) && print "Number of spectra after sorting: " . scalar(@$sorted);
$numberoftotalhits = scalar(@$sorted);
## Update the selected hits
$debugmsg .= update_validated_hits( $reshashref, $sorted, $session, $hashfilename, $dispnhits );

## remove hits where no hit is left when filtering on all hits is used
if ($filter_error_on_all_hits)
{
	$sorted = filter_out_empty( $reshashref, $sorted );
}

#---------------------------------------------------------------------------
#  Displaying the page indexes
#---------------------------------------------------------------------------
## How many pages are there in total
my $pagestotal = int( $numberoftotalhits / $dispmax );
my $rest       = $numberoftotalhits % $dispmax;
if ($rest) { $pagestotal++ }
## Which page should be displayed?
## my $pagedisp = $session->param('page');
## unless ($pagedisp) { $pagedisp = 1 }
### Check if the user has changed the result per page and if the page is possible
if ( $pagedisp > $pagestotal )
{
	$pagedisp = 1;
}
## Generate the html for the navigation
my $pagenavi = new HTML::PageIndex;
## Make the base url
my $baseurl = File::Spec->catfile( $resulturlbase, $id );
my $pagenavig = $pagenavi->makeindex( $pagestotal, $pagedisp, "$xxquest?id=$id&$urlparam", "page", 0 );
my $pagenavihtml;
if ($pagenavig)
{
	$pagenavihtml .= $form->hr;
	$pagenavihtml .= '<TABLE><tr>';
	$pagenavihtml .= '<TD>Display Page</TD><TD>';
	$pagenavihtml .= $pagenavig;
	$pagenavihtml .= '</TD></TR></TABLE>';
	$pagenavihtml .= $form->hr;
}
if ( $reporttype == 1 )
{
	$pagenavihtml .= "<div>Number of spectra after sorting: " . scalar(@$sorted) . "<div>";

	#$pagenavihtml .= "$disptype";
}

#---------------------------------------------------------------------------
#  DEBUG HASH THAT MAY BE PRINTED
#---------------------------------------------------------------------------
if ($debugarray)
{
	foreach my $line (@$debugarray)
	{
		$pagenavihtml .= "$line<br>";
	}
}
### Set the offset parameters
my $offset = ( $pagedisp - 1 ) * $dispmax;
my $limit  = ( $offset + $dispmax ) - 1;

#---------------------------------------------------------------------------
#  Generate and Print the HTML table
#---------------------------------------------------------------------------
#slice the array according to the offset and limit if HTML output is selected
# DONT SLICE IF EXCEL OUTPUT IS SELECTED
if ( $reporttype == 1 )
{
	$sorted = slice_array( $offset, $limit, $sorted );
}
my $htmltable;
my $ruid = $offset + 1;
if ( $reporttype == 1 )
{
	$htmltable .= '<div class="tabletext" width="100%">' . "\n";
	$htmltable .= '<table border="1"  width="100%">';
	
}
my $excel = 0;
my $exelfh;
my $header;
my $outfile;
my $outfilebasename;

if ( $reporttype == 2 )
{
	$outfile = File::Spec->catfile( $resultdirbase, $id, "results.xls" );
	open EXCEL, ">$outfile" or die $!;
	
	$outfilebasename = basename($id) . ".xls";
	# $debugmsg .= "FN: $outfilebasename<br>";
	$exelfh          = *EXCEL;
}
if ( $reporttype == 3 )
{
	$outfile = File::Spec->catfile( $resultdirbase, $id, "results.idXML" );
	open EXCEL, ">$outfile" or die $!;
	$outfilebasename = basename($id) . ".idXML";
	$exelfh          = *EXCEL;
}
if ( $reporttype == 4 )
{
	$outfile = File::Spec->catfile( $resultdirbase, $id, "xtract.csv" );
	$outfilebasename = basename($id) . "_xtract.csv";
	open EXCEL, ">$outfile" or die $!;
	$exelfh = *EXCEL;
}
if ( $reporttype == 5 )
{
	$outfile = File::Spec->catfile( $resultdirbase, $id, "xquest_inclusionlist.xls" );
	open EXCEL, ">$outfile" or die $!;
	$outfilebasename = basename($id) . "_IL.xls";
	$exelfh          = *EXCEL;
}
my $storehash = {};
## a hash for storing the filtered results
my @hitheader = qw(Rank Sequence Protein1 Protein2 Matchodds TIC wTIC xcorrx xcorrb intsum deltaS mions ld-Score rel_error[ppm] type xl-type nAA1 nAA2 &Delta;AA nseen view Sel comment annotation a FDR);
## ITERATE THROUGH THE FILTERED RESULTS AND GENERATE EXCEL OR HTML OUTPUT
foreach my $key (@$sorted)
{
### Get the xml element
	my $xmlelm = ( $reshashref->{$key}->{'xml'} );
	unless ($xmlelm)
	{

		#print "ERROR: NO XML Element found<br>";
		next;
	}
	my $specsearch = XML::TreeBuilder->new;
	$specsearch->parse($xmlelm);
### Get all search hits for this spectrum
	my @search_hits  = map { $_ } $specsearch->find('search_hit');
	my $specname     = $specsearch->attr('spectrum');
	my $scantype     = $specsearch->attr('scantype');
	my @msid         = split( "\\.", $specname );
	my $precursormz  = sprintf( "%.4f", $specsearch->attr('mz_precursor') );
	my $charge       = $specsearch->attr('charge_precursor');
	my $measuredmass = sprintf( "%.4f", $specsearch->attr('Mr_precursor') );
	my $rts          = $specsearch->attr('rtsecscans');                           ## is in seconds
	my @rtimes       = split( /:/, $rts );
	my $rt1          = $rtimes[0];
	my $rt2          = $rtimes[1];
	my $rtstring     = sprintf( "%.2f", $rt1 ) . ":" . sprintf( "%.2f", $rt2 );
	my $mzs          = $specsearch->attr('mzscans');
	my @mzvalues     = split( /:/, $mzs );
	my $mz1          = $mzvalues[0];
	my $mz2          = $mzvalues[1];
	my $mzstring     = sprintf( "%.3f", $mz1 ) . ":" . sprintf( "%.3f", $mz2 );
	if ( $reporttype == 1 )
	{
		$htmltable .= generate_spectrum_header_tr( \@hitheader, $specname, $scantype, $msid[0], $precursormz, $charge, $measuredmass, $rtstring, $mzstring );
		### Print the Header for the search hits
		$htmltable .= '<tr>' . "\n";
		foreach my $string (@hitheader)
		{
			$htmltable .= "<td>$string</td>";
		}
		$htmltable .= '</tr>' . "\n";
	}
	### print the search hits
	#---------------------------------------------------------------------------
	#  ITERATE through all search hits
	#---------------------------------------------------------------------------
	foreach my $hit ( sort { $b->attr('score') <=> $a->attr('score') } @search_hits )
	{
		my $rank      = $hit->attr('search_hit_rank');
		my $type      = $hit->attr('type');
		my $error     = $hit->attr('error');
		my $error_rel = $hit->attr('error_rel');
		## Check if validation is set for this hit
		my $valid   = $reshashref->{$key}->{'hits'}->{$rank}->{'validated'};
		my $comment = $reshashref->{$key}->{'hits'}->{$rank}->{'comment'};
		## Check if disptype 16 (Selected hits) is on and check if the hit has been validated (must not be always rank 1)
		## If unique ids filter is not set then bypass this since the redundant IDs are not validated
		if ( $disptype == 16 && ( $unique_ids == 1 ) )
		{
			unless ($valid) { next }
		}
		if ( $selectallhits == 1 )
		{
			$valid = 1;
		}
		if ( $unselectallhits == 1 )
		{
			$valid = 0;
		}
		############# FILTERS ############
		### CHECK IF NRANKS FILTER IS ON
		if ( $dispnhits ne "all" )
		{
			## Only filter if seleted hits is not on and filter for unique restraints is on
			unless ( $disptype == 16 && ( $unique_ids == 1 ) )
			{
				if ( $rank > $dispnhits ) { next; }
			} else
			{
				## Check if it is not rank 1
				## and is validated
				if ( $rank > 1 )
				{
					next unless ($valid);
				}
			}
		}
		### CHECK IF FILTERING OF ERROR IS ENABLED ON ALL HITS, IF YES FILTER
		if ( ($filter_error_on_all_hits) )
		{
			if ( $error_rel < $maxerrorfrom || $error_rel > $maxerrorto )
			{
				next;
			}
		}
		if ( $reporttype == 1 )
		{
			$htmltable .= "<tr>\n";
		}
		### The id, sequences, and proteins
		my $id             = $hit->attr('id');
		my $seq1           = $hit->attr('seq1');
		my $seq2           = $hit->attr('seq2');
		my $concatsequence = $seq1 . $seq2;
		##
		#---------------------------------------------------------------------------
		#  Define the Modifications
		#---------------------------------------------------------------------------
		my $modseq1 = get_mods( $seq1, $PARAMS, "-a" );
		my $modseq2 = get_mods( $seq2, $PARAMS, "-b" );
		my @modarray;
		foreach my $modstring (@$modseq1)
		{
			push @modarray, $modstring;
		}
		foreach my $modstring (@$modseq2)
		{
			push @modarray, $modstring;
		}
		my $modstring = join( ":", @modarray );
		unless ($seq2)
		{
			$seq2 = "-";
		}
		my $spidp1 = $hit->attr('prot1');
		my $spidp2 = $hit->attr('prot2');
		unless ($spidp2) { $spidp2 = "-" }
		#### Mr, Mz, z ####
		my $mr = sprintf( "%.3f", $hit->attr('measured_mass') );
		my $mz = sprintf( "%.3f", $hit->attr('mz') );
		my $z  = sprintf( "%.3f", $hit->attr('charge') );
		my $z2 = $hit->attr('charge');
		#### Scoring ####
		my $score             = sprintf( "%.2f", $hit->attr('score') );
		my $matchodds         = sprintf( "%.2f", $hit->attr('match_odds') );
		my $wmatchoddssum     = sprintf( "%.2f", $hit->attr('weighted_matchodds_sum') );
		my $wmatchoddsmean    = sprintf( "%.2f", $hit->attr('weighted_matchodds_mean') );
		my $tic               = sprintf( "%.2f", $hit->attr('TIC') );
		my $tica              = sprintf( "%.2f", $hit->attr('TIC_alpha') );
		my $ticb              = sprintf( "%.2f", $hit->attr('TIC_beta') );
		my $wtic              = sprintf( "%.2f", $hit->attr('wTIC') );
		my $xcorrx            = sprintf( "%.2f", $hit->attr('xcorrx') );
		my $xcorrb            = sprintf( "%.2f", $hit->attr('xcorrb') );
		my $intsumscore       = sprintf( "%.0f", $hit->attr('intsum') );
		my $apriorymatchscore = sprintf( "%.1f", $hit->attr('apriori_match_probs_log') );
		my $seriesscore       = sprintf( "%.1f", $hit->attr('series_score_mean') );
		my $fdr               = "-";

		if ( defined( $hit->attr('fdr') ) )
		{
			$fdr = $hit->attr('fdr');
			$fdr = sprintf( "%.3f", $fdr );
		}
		my $ePvalue = "-";
		if ( defined( $hit->attr('ePvalue') ) )
		{
			$ePvalue = $hit->attr('ePvalue');
			$ePvalue = sprintf( "%.3f", $ePvalue );
		}
		my $apriorymatchscorefull = $hit->attr('apriori_match_probs_log');
		my $deltascore            = $reshashref->{$key}->{'hits'}->{$rank}->{'deltascore'};
		my $nminmatchedions       = $reshashref->{$key}->{'hits'}->{$rank}->{'nminmatchedions'};
		$deltascore = sprintf( "%.2f", $deltascore );
		### Determine the type of cross-link
		my $xlinktype = "-";
		if ( $type eq "xlink" )
		{
			$xlinktype = get_type_of_xlink( $spidp1, $spidp2 );
		} else
		{
			$xlinktype = get_type_of_mono_intra_link( $spidp1, $type );
		}
		### Determine the absolute positions of the cross-linked residues
		my $xlinkpos = $hit->attr('xlinkposition');    ## E.g. 3,5
		my ( $res1, $res2 ) = get_absposition( $idx, $xlinkpos, $hit );
		unless ($res2) { $res2 = "n/a" }
		unless ($res1) { $res1 = "n/a" }
		my $distinsequence = $reshashref->{$key}->{'hits'}->{$rank}->{'distinsequence'};
		my @pdbarray;
		my $pdbstr;
		my @structureannotation = map { $_ } $hit->find('structure_annotation');
		my $smallestdist;
		### Print the structure annotation
		if (@structureannotation)
		{
			if ( $reporttype == 1 )
			{
				$pdbstr = generate_dist_annotation( \@structureannotation );
				push( @pdbarray, $pdbstr );
			}
			if ( $reporttype == 2 )
			{
				( $pdbstr, $smallestdist ) = generate_dist_annotation_excel( \@structureannotation );
				push( @pdbarray, $pdbstr );
			}
		} else
		{
			## Check if the hit is annotated
			my $annotation = $hit->attr('annotation');
			if ($annotation)
			{
				$pdbstr = $annotation;
			} else
			{
				$pdbstr = "-";
			}
			push( @pdbarray, $pdbstr );
		}
		my $structureannotaition = join( ",", @pdbarray );
		my $annotationurl = join "", 'xannotate.cgi?', 'protein1=', $hit->attr('prot1'), ";protein2=", $hit->attr('prot2'), ';plottype=profile', ";id=$id", ";hitid=", $hit->attr('id'), ";database=$database", ";xquestdef=", $xquestdef, ";topology=", $hit->attr('topology'), ';seq1=',
		  $hit->attr('seq1'), ';seq2=', $hit->attr('seq2'), ";xlinkposition=", $hit->attr('xlinkposition');
		my $annotation = '<a target="_blank" href="' . $annotationurl . '">seq</a>';
		my $nseen      = $reshashref->{$key}->{'hits'}->{$rank}->{'nseen'};
		my $specfn     = $reshashref->{$key}->{'hits'}->{$rank}->{'specfile'};
		my $specurl    = get_xions_url( $hit, $key, $rank, "spectrum", $specname, $specfn, $scantype, "view" );
		my $xcorrburl  = get_xions_url( $hit, $key, $rank, "xcorr", $specname, $specfn, $scantype, $xcorrb );
		my $xcorrxurl  = get_xions_url( $hit, $key, $rank, "xcorr", $specname, $specfn, $scantype, $xcorrx );

		#my $ePvalueurl = get_epvalue_url($ePvalue);
		my $uxID = "-";
		if ( defined( $reshashref->{$key}->{'hits'}->{$rank}->{'uxID'} ) )
		{
			$uxID = $reshashref->{$key}->{'hits'}->{$rank}->{'uxID'};
		}
		my $xp_f = $reshashref->{$key}->{'hits'}->{1}->{'xprophet_f'};
		## write to excel
		if ( $reporttype == 2 )
		{
			my @excelheader = qw(Rank Id Protein1 Protein2 Type XLType Spectrum AbsPos1 AbsPos2 deltaAA Annotation Mr Mz z Error_rel[ppm] nseen mions MatchOdds Xcorrx Xcorrb TIC TicA TicB WTIC intsum deltaS ld-Score FDR uxID comment);
			unless ($header)
			{
				print_to_excel( $exelfh, \@excelheader );
				$header = 1;
			}
			my @excelarray =
			  ( $rank, $id, $spidp1, $spidp2, $type, $xlinktype, $specname, $res1, $res2, $distinsequence, $structureannotaition, $mr, $mz, $z, $error_rel, $nseen, $nminmatchedions, $matchodds, $xcorrx, $xcorrb, $tic, $tica, $ticb, $wtic, $intsumscore, $deltascore, $score, $fdr, $uxID, $comment );
			print_to_excel( $exelfh, \@excelarray );
		}
		## write table
		if ( $reporttype == 1 )
		{
			## format the protein ids so that there can be line breaks
			my @prots1 = split( ",", $spidp1 );
			my @prots2 = split( ",", $spidp2 );
			my $spidp1str       = get_prot_link($spidp1);
			my $spidp2str       = get_prot_link($spidp2);
			my $htmltablefields = {};
			my $formathash      = {};
### INFO: IF COLUMNS SHOULD APPEAR IN THE TABLE ADD THE HEADER COLUMN TO THE @hitheader
###		  the hash index must match with the header
###		  THE formathash may be used to add additional tags
			$formathash->{'ld-Score'}            = 'BGCOLOR="#FFBE9F"';
			$htmltablefields->{'Rank'}           = $rank;
			$htmltablefields->{'Sequence'}       = $id;
			$htmltablefields->{'Protein1'}       = $spidp1str;
			$htmltablefields->{'Protein2'}       = $spidp2str;
			$htmltablefields->{'Matchodds'}      = $matchodds;
			$htmltablefields->{'TIC'}            = $tic;
			$htmltablefields->{'wTIC'}           = $wtic;
			$htmltablefields->{'xcorrx'}         = $xcorrxurl;
			$htmltablefields->{'xcorrb'}         = $xcorrburl;
			$htmltablefields->{'intsum'}         = $intsumscore;
			$htmltablefields->{'deltaS'}         = $deltascore;
			$htmltablefields->{'mions'}          = $nminmatchedions;
			$htmltablefields->{'ld-Score'}       = $score;
			$htmltablefields->{'score'}          = $score;
			$htmltablefields->{'rel_error[ppm]'} = $error_rel;
			$htmltablefields->{'type'}           = $type;
			$htmltablefields->{'xl-type'}        = $xlinktype;
			$htmltablefields->{'nAA1'}           = $res1;
			$htmltablefields->{'nAA2'}           = $res2;
			$htmltablefields->{'&Delta;AA'}      = $distinsequence;
			$htmltablefields->{'nseen'}          = $nseen;
			$htmltablefields->{'view'}           = $specurl;
			$htmltablefields->{'annotation'}     = $structureannotaition;
			$htmltablefields->{'a'}              = $annotation;
			$htmltablefields->{'FDR'}            = $fdr;

			#$htmltablefields->{'eP'}             = $ePvalueurl;
			#$htmltablefields->{'xp_f'}           = $xp_f;
			#$htmltablefields->{'uxID'}           = $uxID;
			my $validationname = "validation_" . $key . "_" . $rank;
			my $commentname    = "comment_" . $key . "_" . $rank;
			if ( $valid == 1 )
			{
				$form->param( $validationname => 1 );
			} else
			{
				$form->param( $validationname => -1 );
			}
			$htmltablefields->{'Sel'} = $form->checkbox( -name => $validationname, -value => 1, -checked => $valid, -label => '' );
			$htmltablefields->{'comment'} = $form->textfield( -name => $commentname, -value => $comment, -size => 2, -maxlength => 25 );
			$htmltable .= print_td_hash( $htmltablefields, $formathash, \@hitheader );
			$htmltable .= "</tr>\n";
		}
		if ( $reporttype == 5 )
		{
			### Generate Inclusion List
			# 1. collect hits here, 2.generate IL, 3. store to excelfile
			if ( $scantype eq "light_heavy" )
			{
				## then make 2 entries, light and heavy, make a hash with rt and an array for the mz values, allowas also several mzs at the same rt
				#push @{ $storehash->{$rt1} }, $mz1;
				#push @{ $storehash->{$rt2} }, $mz2;
				push @{ $storehash->{$id}->{"light"} }, [ $mz1, $rt1 / 60, $score, $z2, $id . "-light", $type ];
				push @{ $storehash->{$id}->{"heavy"} }, [ $mz2, $rt2 / 60, $score, $z2, $id . "-heavy", $type ];
			} else
			{

				#push @{ $storehash->{$rt1} }, $mz1;
				push @{ $storehash->{$id}->{"light"} }, [ $mz1, $rt1 / 60, $score, $z2, $id . "-light", $type ];
			}
		}
		if ( $reporttype == 3 )
		{
			### Generate Id.xml
			#print $scantype."<br>";
			# 1. collect hits here, 2.generate Id.xml
			## Collect light and heavy as a separate ID
			if ( $scantype eq "light_heavy" )
			{
				## then make 2 entries, light and heavy, make a hash with rt and an array for the mz values, allowas also several mzs at the same rt
				## cat an L and an H for the light / heavy at the end
				my $concatsequencel = $id . ":L";
				my $concatsequenceh = $id . ":H";
				push @{ $storehash->{$id}->{"light"} }, [ $mz1, $rt1, $score, $z2, $id, $type ];
				push @{ $storehash->{$id}->{"heavy"} }, [ $mz2, $rt2, $score, $z2, $id, $type ];
			} else
			{
				$concatsequence = $id . ":L";
				push @{ $storehash->{$id}->{"light"} }, [ $mz1, $rt1, $score, $z2, $id, $type ];
			}
		}
		if ( $reporttype == 4 )
		{
			### Generate xTract csv
			# 1. collect hits here, 2.generate csv file
			## Collect light and heavy as a separate IDs
			## Allow to use redundant IDs
			## cut the spectrumname
			my ( $lightscanname, $heavyscanname ) = cut_spectrum_filename($specname);
			
			# Generate the protein identifier, sep proteins by ":"
			my $protidcsv="-";
			$protidcsv =  $spidp1 . ":" . $spidp2;
			$protidcsv =  $spidp1 if ( $spidp2 eq "-");
			
			# Mask if a comma is in the protein string
			if ($protidcsv =~ m/,/){
			$protidcsv = '"'.$protidcsv.'"';
			}
			
			
			if ( $scantype eq "light_heavy" )
			{
				## then make 2 entries, light and heavy, make a hash with rt and an array for the mz values, allow also several mzs at the same rt
				my $peptidehash = {};
				## Light Feature
				$peptidehash->{'seq'}   = $id;
				$peptidehash->{'score'} = $score;
				$peptidehash->{'fdr'}   = $fdr;
				$peptidehash->{'prot'}  = $protidcsv;
				$peptidehash->{'type'}  = $type . ":light";
				$peptidehash->{'mz'}    = $mz1;
				$peptidehash->{'tr'}    = $rt1;
				$peptidehash->{'z'}     = $z2;
				$peptidehash->{'scan'}  = $lightscanname;
				$peptidehash->{'mod'}   = $modstring;
				$peptidehash->{'uxID'}  = $uxID;
				$peptidehash->{'preINT'}  = "";
				push @{ $storehash->{$id}->{"light"} }, $peptidehash;
				## Heavy Feature
				$peptidehash            = {};
				$peptidehash->{'seq'}   = $id;
				$peptidehash->{'score'} = $score;
				$peptidehash->{'fdr'}   = $fdr;
				$peptidehash->{'prot'}  = $protidcsv;
				$peptidehash->{'type'}  = $type . ":heavy";
				$peptidehash->{'mz'}    = $mz2;
				$peptidehash->{'tr'}    = $rt2;
				$peptidehash->{'z'}     = $z2;
				$peptidehash->{'scan'}  = $heavyscanname;
				$peptidehash->{'mod'}   = $modstring;
				$peptidehash->{'uxID'}  = $uxID;
				$peptidehash->{'preINT'}  = "";
				push @{ $storehash->{$id}->{"heavy"} }, $peptidehash;
			} else
			{
				my $peptidehash = {};

				# ToDo: SCAN Numbers are currently not reported in xquest.xml
				## Light Feature
				$peptidehash->{'seq'}   = $id;
				$peptidehash->{'score'} = $score;
				$peptidehash->{'fdr'}   = $fdr;
				$peptidehash->{'prot'}  = $protidcsv;
				$peptidehash->{'type'}  = $type . ":light";
				$peptidehash->{'mz'}    = $mz1;
				$peptidehash->{'tr'}    = $rt1;
				$peptidehash->{'z'}     = $z2;
				$peptidehash->{'scan'}  = $lightscanname;
				$peptidehash->{'mod'}   = $modstring;
				$peptidehash->{'uxID'}  = $uxID;
				$peptidehash->{'preINT'}  = "";
				push @{ $storehash->{$id}->{"light"} }, $peptidehash;
			}
		}
	}
	$specsearch->delete();
	$ruid++;
}


if ( $reporttype == 1 )
{
	$htmltable .= '</table>';
	$htmltable .= '</div>';
	$htmltable .= "<div style=\"position:relative;border:0px solid #00ff00;bottom:-25px; text-align: right\" >";
	$htmltable .= "xQuest/xProphet viewer version $viewerversion, Thomas Walzthoeni / ETH Zurich";
	$htmltable .= '</div>';
}

#---------------------------------------------------------------------------
# PRINT THE TABLE, THE SUMMARY OR SEND THE EXCEL FILE
#---------------------------------------------------------------------------
if ( $reporttype == 1 )
{
	print_header();
	print $htmlheader;
	print_form();    # Header
	if ($debug)
	{
		debug_param($session);
	}
	if ($debugmsg)
	{
		print $debugmsg;
	}
	print $pagenavihtml;
	print $htmltable;
	print $form->end_form;
	print $form->end_html;
}
if ( $reporttype == 2 )
{
	close($exelfh);
	send_excel_file();
}

# idXML file
if ( $reporttype == 3 )
{

	#print_header();
	#print '<h1>xQuest results viewer - Summary</h1><hr>' . "\n";
	#print_form();
	#print Dumper($storehash);
	my @proteinelements;
	my @peptideelements;
	my $idxmlheader;
	$idxmlheader .= '<?xml version="1.0" encoding="UTF-8"?>' . "\n";
	$idxmlheader .= '<?xml-stylesheet type="text/xsl" href="file:////OpenMS/share/OpenMS/XSL/IdXML.xsl"?>' . "\n";
	$idxmlheader .= '<IdXML version="1.2" xsi:noNamespaceSchemaLocation="http://open-ms.sourceforge.net/schemas/IdXML_1_2.xsd" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">' . "\n";
	$idxmlheader .= '<SearchParameters id="SP_0" db="./current.fasta" db_version="" taxonomy="" mass_type="monoisotopic" charges="" enzyme="trypsin" missed_cleavages="0" precursor_peak_tolerance="0" peak_mass_tolerance="0" >' . "\n";
	$idxmlheader .= '<FixedModification name="C+57.0215" />' . "\n";
	$idxmlheader .= '</SearchParameters>' . "\n";
	$idxmlheader .= '<IdentificationRun date="" search_engine="xQuest" search_engine_version="" search_parameters_ref="SP_0" >' . "\n";
	print $exelfh $idxmlheader . "\n";
	## get all ids
	my @ids         = keys %$storehash;
	my $numelements = @ids;
	for ( my $i = 0 ; $i < $numelements ; $i++ )
	{
		my $id            = $ids[$i];
		my $lightfeatures = $storehash->{$id}->{"light"};
		my $heavyfeatures = $storehash->{$id}->{"heavy"};

		#print "<br>$id<br>";
		my $proteinidxml = '<ProteinHit id="PH_' . $i . '" accession="XL-Id-' . $i . '" score="0" sequence="' . $id . '" ></ProteinHit>';
		push @proteinelements, $proteinidxml;
		my $pepidexml;
		foreach my $lightfeaturearray (@$lightfeatures)
		{
			my $mz        = @$lightfeaturearray[0];
			my $rt        = @$lightfeaturearray[1];
			my $score     = @$lightfeaturearray[2];
			my $z         = @$lightfeaturearray[3];
			my $concatseq = @$lightfeaturearray[4];
			my $type      = @$lightfeaturearray[5];

			#print "$rt, $mz, $score<br>";
			$pepidexml = '<PeptideIdentification score_type="ld-score" higher_score_better="true" significance_threshold="0" MZ="' . $mz . '" RT="' . $rt . '" >"' . "\n";
			$pepidexml .= '<PeptideHit score="' . $score . '" sequence="' . $concatseq . '" label="light" charge="' . $z . '" type="' . $type . '" protein_refs="PH_' . $i . '" ></PeptideHit>' . "\n";
			$pepidexml .= '</PeptideIdentification>' . "\n";
			push @peptideelements, $pepidexml;
		}
		foreach my $heavyfeaturearray (@$heavyfeatures)
		{
			my $mz        = @$heavyfeaturearray[0];
			my $rt        = @$heavyfeaturearray[1];
			my $score     = @$heavyfeaturearray[2];
			my $z         = @$heavyfeaturearray[3];
			my $concatseq = @$heavyfeaturearray[4];
			my $type      = @$heavyfeaturearray[5];

			#print "$rt, $mz, $score<br>";
			$pepidexml = '<PeptideIdentification score_type="ld-score" higher_score_better="true" significance_threshold="0" MZ="' . $mz . '" RT="' . $rt . '" >"' . "\n";
			$pepidexml .= '<PeptideHit score="' . $score . '" sequence="' . $concatseq . '" label="heavy" charge="' . $z . '" type="' . $type . '" protein_refs="PH_' . $i . '" ></PeptideHit>' . "\n";
			$pepidexml .= '</PeptideIdentification>' . "\n";
			push @peptideelements, $pepidexml;
		}
	}
	print $exelfh '<ProteinIdentification score_type="" higher_score_better="true" significance_threshold="0" >' . "\n";
	foreach my $protidxml (@proteinelements)
	{
		print $exelfh $protidxml . "\n";
	}
	print $exelfh "</ProteinIdentification>\n";
	foreach my $peptideidxml (@peptideelements)
	{
		print $exelfh $peptideidxml . "\n";
	}
	my $idxmlfooter = '</IdentificationRun>
</IdXML>';
	print $exelfh $idxmlfooter;
	close($exelfh);
	## send the excel file
	send_idxml_file();
}

# xtract csv format
if ( $reporttype == 4 )
{
	my @proteinelements;
	my @peptideelements;
	## Print the header to the output file
	## Write the headerline
	#my $headerline = "scan,seq,prot,mod,score,type,mz,z,tr,fdr, uxID";
	my @header = qw(scan seq prot mod score type mz z tr fdr uxID preINT);
	my $headerline = join( ",", @header );
	print $exelfh $headerline . "\n";
	## get all ids, id is the seq, filter out redundant spectra
	my @ids         = keys %$storehash;
	my $numelements = @ids;
	my $seen_scan_seq = {};
	
	for ( my $i = 0 ; $i < $numelements ; $i++ )
	{
		my $id                = $ids[$i];
		my $lightfeaturearray = $storehash->{$id}->{"light"};
		my $heavyfeaturearray = $storehash->{$id}->{"heavy"};
		
		if ($lightfeaturearray)
		{
			foreach my $lightfeaturehash (@$lightfeaturearray)
			{
				my $scan = $lightfeaturehash->{'scan'};
				my $seq = $lightfeaturehash->{'seq'};
				my $scan_id = $scan."_".$seq;
				
				if ( $seen_scan_seq->{$scan_id} )
				{
					next;
				} else
				{
					$seen_scan_seq->{$scan_id} = 1;
				}
				my @resultsarray;
				foreach my $key (@header)
				{
					push @resultsarray, $lightfeaturehash->{$key};
				}
				my $string = join( ",", @resultsarray );
				print $exelfh $string . "\n";
			}
		}
		if ($heavyfeaturearray)
		{
			foreach my $heavyfeaturehash (@$heavyfeaturearray)
			{
				my $scan = $heavyfeaturehash->{'scan'};
				my $seq = $heavyfeaturehash->{'seq'};
				my $scan_id = $scan."_".$seq;		
						
				if ( $seen_scan_seq->{$scan_id} )
				{
					next;
				} else
				{
					$seen_scan_seq->{$scan_id} = 1;
				}
				
				my @resultsarray;
				foreach my $key (@header)
				{
					push @resultsarray, $heavyfeaturehash->{$key};
				}
				my $string = join( ",", @resultsarray );
				print $exelfh $string . "\n";
			}
		}
	}
	close($exelfh);
	## send the excel file
	send_xtractcsv_file();
}
## IL
if ( $reporttype == 5 )
{

	#print_header();
	#print '<h1>xQuest results viewer - Inclusion List Generator</h1><hr>' . "\n";
	#print_form();
	## Generation of IL
	# 1. check the number of masses / if > 2000 make more ILs
	# 2. make one entry for each peptide +/-5minutes rt
	# FORMAT
	# MZ \t START(min) \t END (MIN) \t NAME
	my @rtkeys      = keys %$storehash;
	my @ids         = keys %$storehash;
	my $nids        = @ids;
	my $numelements = @ids;
	my $maxsizeofIL = 2000;
	my $rthash;
	my $deltart        = 5;                                 # in minutes
	my $seen           = {};
	my $nbins          = 10;
	my $gradientlength = 90;
	my $binlength      = int( $gradientlength / $nbins );
	my $results        = {};

	for ( my $i = 0 ; $i < $nids ; $i++ )
	{
		## Get light and heavy feature array
		my $id = $ids[$i];

		#print "ID: $id<br>";
		my $lightarray = $storehash->{$id}->{"light"};

		#print Dumper ($lightarray);
		my $heavyarray = $storehash->{$id}->{"heavy"};

		# array structure [ $mz1, $rt1, $score, $z2, $id, $type ]
		if ($lightarray)
		{
			foreach my $array (@$lightarray)
			{

				#print Dumper ($array);
				## get the result hash (mz starttr endtr id)
				## key trstart is the start time
				## key resarray is the result array
				my $reshash  = generate_result_hash( $nbins, $binlength, $array );
				my $trstart  = $reshash->{'trstart'};
				my $resarray = $reshash->{'resarray'};
				push @{ $results->{$trstart} }, $resarray;
			}
		}
		if ($heavyarray)
		{
			foreach my $array (@$heavyarray)
			{
				## get the result hash (mz starttr endtr id)
				## key trstart is the start time
				## key resarray is the result array
				my $reshash = generate_result_hash( $nbins, $binlength, $array );

				#print Dumper ($reshash);
				#exit;
				my $trstart  = $reshash->{'trstart'};
				my $resarray = $reshash->{'resarray'};
				push @{ $results->{$trstart} }, $resarray;
			}
		}
	}
	## Sort the keys
	my @sortedkeys = sort { $a <=> $b } keys %$results;
	foreach my $key (@sortedkeys)
	{
		my $idarray = $results->{$key};
		foreach my $array (@$idarray)
		{

			#print Dumper ($array);
			my $reline = join( "\t", @$array );

			#print "$reline<br>";
			print $exelfh $reline . "\n";
		}
	}
	close($exelfh);
	send_il_file();
}

sub get_version_from_file
{
	my $filename = shift;
	open FILE, "<$filename" or die " No such file: $filename $!";
	my $firstLine = <FILE>;
	close(FILE);
	return $firstLine;
}

sub write_version_to_file
{
	my $filename = shift;
	my $version  = shift;
	open FILE, ">$filename" or die $!;
	print FILE $version . "\n";
	close(FILE);
}

sub generate_result_hash
{
	my $nbins        = shift;
	my $binlength    = shift;
	my $featurearray = shift;

	#print Dumper ($featurearray);
	#exit;
	my $mz = $featurearray->[0];

	#print $mz."<br>";
	#exit;
	my $tr = $featurearray->[1];
	my $id = $featurearray->[4];
	my @resultsar;
	my $resul    = {};
	my $foundbin = 0;
	foreach my $b ( 0 .. ( $nbins - 1 ) )
	{
		my $startrt = $b * $binlength;
		my $endrt   = ( $b + 1 ) * $binlength;

		#print "Bin $b: $startrt to $endrt<br>";
		if ( $tr >= $startrt && $tr < $endrt )
		{
			$foundbin = 1;

			#print "in bin!<br>";
			## Generate line;
			$resul->{'trstart'} = $startrt;
			push @resultsar, $mz;
			push @resultsar, $startrt;
			push @resultsar, $endrt;
			push @resultsar, $id;
			$resul->{'resarray'} = \@resultsar;
		}
	}
	unless ($foundbin)
	{
		print "NO BIN FOUND FOR\n";
		print Dumper ($featurearray);
		exit;
	}

	#print Dumper ($resul);
	return $resul;
}

sub get_type_of_monolink_looplink
{
	my $spidp1 = shift;
	my $type   = shift;
	my $typetoreport;
	my @prots1 = split( ",", $spidp1 );
	my $decoy = 0;
	foreach my $prot1 (@prots1)
	{
### Check if protein is a decoy protein
		if ( $prot1 =~ m/decoy/ )
		{
			$decoy = 1;
		}
	}
	if ($decoy)
	{
		$typetoreport = "decoy " . $type;
	} else
	{
		$typetoreport = $type;
	}
	return $typetoreport;
}

#---------------------------------------------------------------------------
#  Functions
#---------------------------------------------------------------------------
sub cut_spectrum_filename
{
	my $spectrumfilename = shift;
	my @specfilenames;

	#print "Specfilename: $spectrumfilename<br>";
## split the filename, get splitsites: e.g.: wathomas_M1111_179.c.02974.02974.4_wathomas_M1111_179.c.02776.02776.4
	my @cutsites;
	my $bnfn = basename($spectrumfilename);
	while ( $bnfn =~ /_/gi )
	{
		push @cutsites, pos($bnfn);
	}
## substring CUT STRING IN THE MIDDLE _
	my $cutat;
	my $numcutsites = @cutsites;
	for ( my $i = 1 ; $i < $numcutsites + 1 ; $i++ )
	{

		#print "Checking $i of $numcutsites cutsites<br>";
		if ( ( $i * 2 > $numcutsites ) && !( $i > ( $numcutsites + 1 ) / 2 ) )
		{
			$cutat = $i;

			#print "Cutsite to split:$i<br>";
		}
	}
	my $cut1    = $cutsites[ $cutat - 1 ];
	my $string1 = substr( $bnfn, 0, $cut1 - 1 );
	my $string2 = substr( $bnfn, $cut1 );
	return ( $string1, $string2 );
}

sub send_excel_file
{
	print "Content-type: application/octet-stream\n";
	print "Content-Disposition: attachment;filename=$outfilebasename\n\n";
	print_file($outfile);
}

sub send_idxml_file
{
	print "Content-type: application/octet-stream\n";
	print "Content-Disposition: attachment;filename=$outfilebasename\n\n";
	print_file($outfile);
}

sub send_xtractcsv_file
{
	print "Content-type: application/octet-stream\n";
	print "Content-Disposition: attachment;filename=$outfilebasename\n\n";
	print_file($outfile);
}

sub send_il_file
{
	print "Content-type: application/octet-stream\n";
	print "Content-Disposition: attachment;filename=$outfilebasename\n\n";
	print_file($outfile);
}

sub calc_delta_score
{
	my $resultshashref = shift;
	my $sorted         = shift;

	#my $hit;
### ITERATE THROUGH THE FILTERED RESULTS
	foreach my $key (@$sorted)
	{
### Get the xml element
		my $xmlelm = ( $reshashref->{$key}->{'xml'} );
		unless ($xmlelm)
		{
			print "ERROR: NO XML Element found<br>";
			next;
		}

		#my $specsearch = XML::TreeBuilder->new;
		#$specsearch->parse($xmlelm);
### Get all search hits for this spectrum
		my $subhash = $resultshashref->{$key}->{'hits'};
		my @allhitssorted = sort { $a <=> $b } keys %$subhash;

		#print Dumper (@allhitssorted);
		#exit;
		#my @search_hits  = map { $_ } $specsearch->find('search_hit');
		#---------------------------------------------------------------------------
		#  calculate the delta score here and set them as attribute
		#---------------------------------------------------------------------------
		#my @allhitssorted=sort { $b->attr('score') <=> $a->attr('score') } @search_hits;
		for ( my $i = 0 ; $i < ( scalar(@allhitssorted) ) ; $i++ )
		{

			#my $rank=$i+1;
			my $hit       = $allhitssorted[$i];
			my $score     = $resultshashref->{$key}->{'hits'}->{$hit}->{'score'};
			my $error_rel = $resultshashref->{$key}->{'hits'}->{$hit}->{'error'};
			my $id        = $resultshashref->{$key}->{'hits'}->{$hit}->{'structure'};

			#print "Score of current hit $id nr. $hit: $score error is $error_rel ppm\n";
			#exit;
############# FILTERS ############
### CHECK IF FILTERING OF ERROR IS ENABLED ON ALL HITS, IF YES FILTER
			if ( ($filter_error_on_all_hits) )
			{
				if ( $error_rel < $maxerrorfrom || $error_rel > $maxerrorto )
				{
					next;
				}
			}
			my $nextdifferentid;
			my $nextdifferentidscore;
### Get the next id that is different from the current start at the current id +1
			for ( my $z = $i ; $z < scalar(@allhitssorted) ; $z++ )
			{
				my $hit2 = $allhitssorted[ $z + 1 ];
				unless ($hit2)
				{

					#	print "Hit2 not set\n";
					last;
				}
				my $score     = $resultshashref->{$key}->{'hits'}->{$hit2}->{'score'};
				my $error_rel = $resultshashref->{$key}->{'hits'}->{$hit2}->{'error'};
############# FILTERS ############
### CHECK IF FILTERING OF ERROR IS ENABLED ON ALL HITS, IF YES FILTER
				if ( ($filter_error_on_all_hits) )
				{
					if ( $error_rel < $maxerrorfrom || $error_rel > $maxerrorto )
					{
						next;
					}
				}
				my $hit2id = $resultshashref->{$key}->{'hits'}->{$hit2}->{'structure'};

				#print "Score of next hit $hit2id nr. $hit2: $score error is $error_rel ppm\n";
### Check if the structure is the same
				if ( $id eq $hit2id )
				{

					#print "Hits are the same!\n";
					next;
				} else
				{
					$nextdifferentid      = $z;
					$nextdifferentidscore = $resultshashref->{$key}->{'hits'}->{$hit2}->{'score'};
					last;
				}
			}
### Calculate the delta btw the hits
			my $delta;
			unless ($score)
			{
				$delta = "n/a";
			} else
			{
				if ( $score == 0 )
				{
					$delta = "n/a";
				} else
				{
					$delta = $nextdifferentidscore / $score;
				}
			}
			$reshashref->{$key}->{'hits'}->{$hit}->{'deltascore'} = $delta;
		}
	}
	return $resultshashref;
}

#---------------------------------------------------------------------------
#  Update the reshashref for the validated hits
#---------------------------------------------------------------------------
sub update_validated_hits
{
	my $hashref        = shift;
	my $sorted         = shift;
	my $session        = shift;
	my $hashreffile    = shift;
	my $shownranks     = shift;
	my @parameters     = $session->param();
	my $validationhash = {};
	my $commenthash    = {};
	my $debug;
	
    #$debug = Dumper($sorted);
	
## Don't update unless the Form was submitted via "Save Selected"
## But: if selected hits are 
	unless ( $session->param('saveselected') )
	{
		return $debug;
	}

	# Unset the saveselected param, otherwise it is stored and if the
	# page is changed, the hits are updated
	$session->param( 'saveselected', 0 );

	# unset the selectallhits and the unselectallhits
	$session->param( 'selectallhits',   0 );
	$session->param( 'unselectallhits', 0 );
	$form->param( 'selectallhits',   0 );
	$form->param( 'unselectallhits', 0 );

	# Update only on the hits that were displayed
	my $dispmax  = $session->param("dispmax");
	my $pagedisp = $session->param('page');
### Set the offset parameters
	my $offset = ( $pagedisp - 1 ) * $dispmax;
	my $limit  = ( $offset + $dispmax ) - 1;
### Create a hash with entries that can be modified
	$sorted = slice_array( $offset, $limit, $sorted );
	my $addcount;
	my $rmcount;
## Parse all validation parameters
	for ( my $i = 0 ; $i < @parameters ; $i++ )
	{
		my $name = $parameters[$i];
		if ( $name =~ m/validation/ )
		{
			my @validationarray = split( /_/, $name );

			# set key and rank in the validationhash, only one per spectrum is possible
			$validationhash->{ $validationarray[1] } = $validationarray[2];
		}
		if ( $name =~ m/comment/ )
		{
			my @validationarray = split( /_/, $name );

			# set key and rank in the validationhash, only one per spectrum is possible
			$commenthash->{ $validationarray[1] } = $validationarray[2];
		}
	}
	my @validationhashkeys = keys %$validationhash;
## Add the validation if the box is checked
	foreach my $id (@$sorted)
	{
		if ( $validationhash->{$id} )
		{
			## add permanently to validated hits
			$hashref->{$id}->{'hits'}->{ $validationhash->{$id} }->{'validated'} = 1;
			$addcount++;
			## Update the comment for this hit
			#my $commentname= "comment_" . $id . "_" . $validationhash->{$id};
			#$hashref->{$id}->{'hits'}->{ $validationhash->{$id} }->{'comment'} = $session->param($commentname);
		}
	}
## Add the comments
	my @commentshashkeys = keys %$commenthash;
	foreach my $id (@$sorted)
	{
		if ( $commenthash->{$id} )
		{
			## add permanently to validated hits
			#$hashref->{$id}->{'hits'}->{ $validationhash->{$id} }->{'validated'} = 1;
			## Update the comment for this hit
			my $commentname = "comment_" . $id . "_" . $commenthash->{$id};
			$hashref->{$id}->{'hits'}->{ $commenthash->{$id} }->{'comment'} = $session->param($commentname);
		}
	}
## Remove the validation for all hits where the box is not checked
	foreach my $id (@$sorted)
	{
		## get all search hits
		my @ranks = keys %{ $hashref->{$id}->{'hits'} };
		foreach my $rank (@ranks)
		{
			## check if validation is on
			if ( $hashref->{$id}->{'hits'}->{$rank}->{'validated'} == 1 )
			{
				## check if checkbox is not on
				my $rankon = $validationhash->{$id};
				## uncheck if not set
				unless ( $rankon == $rank )
				{
					## only unset if nranks is at least as large as the searchhit!
					if ( $shownranks >= $rank )
					{
						$hashref->{$id}->{'hits'}->{$rank}->{'validated'} = 0;
						$rmcount++;
						## unset the comment
						#$hashref->{$id}->{'hits'}->{$rank}->{'comment'} = 0;
					}
				}
			}
		}
	}
	## Store after adding of validation tags
	$debug = "<div id=\"highlight\"> Added: $addcount Removed: $rmcount selected Ids</div>";

	store $hashref, $hashreffile;
	return $debug;
}

sub generate_spectrum_header_tr
{
	my $hitheaderref = shift;
	my $specname     = shift;
	my $scantype     = shift;
	my $basename     = shift;
	my $mz           = shift;
	my $charge       = shift;
	my $measuredmass = shift;
	my $rts          = shift;
	my $mzs          = shift;
	my $num          = scalar(@$hitheaderref);
	my $html;
	$html .= '<tr class="tableheading">' . "\n";
	$html .= "<td>id $ruid </td>";
	$html .= '<td colspan="' . ( $num - 1 ) . '" align="right">';
	$html .= '<span align="right">Spectrum info: ';

	#$html .= '<span>';
	#$html .= "SpectrumId:$specname,n";
	#$html .= "Scantype: $scantype, ";
	#$html .= "MS basename: $basename<br>";
	$html .= "MS m/z: $mz, ";
	$html .= "Charge: $charge, ";
	$html .= "Precursor mass: $measuredmass, ";
	$html .= "RTs: $rts, ";
	$html .= "MZs: $mzs";

	#$html .= '</span>'.$specname.'</td>';
	$html .= '</tr>' . "\n";
	return $html;
}

sub generate_dist_annotation
{
	my $structureannotation = shift;    ## arrayref
	my $html;
	$html .= '<a align="right" class="infobox2" href="#INFO">Show annotation';
	$html .= '<span>';
	foreach my $structureannotation (@$structureannotation)
	{
		my $pdb       = $structureannotation->attr('pdb');
		my $distinseq = $structureannotation->attr('distanceinsequence');
		my $euklid    = $structureannotation->attr('eukliddist');
		my $griddist  = $structureannotation->attr('griddist');
		my $atopo     = $structureannotation->attr('alpha_topology');
		my $btopo     = $structureannotation->attr('beta_topology');
		my $pdbstr    = $pdb . "|" . $distinseq . "| euklid:" . $euklid . "A | grid:" . $griddist . "A |" . $atopo . "|" . $btopo . "<br>";
		$html .= $pdbstr;

		#<structure_annotation alpha_topology="LYS-41-F" beta_topology="LYS-34-F" distanceinsequence="7" eukliddist="6.2" griddist="6.3" pdb="3DW8.pdb">
	}
	$html .= '</span></a></td>';
	$html .= '</tr>' . "\n";
	return $html;
}

sub generate_dist_annotation_excel
{
	my $structureannotation = shift;    ## arrayref
	my $html;
	my $smallestdistance = 0;
	$html .= '"';
	foreach my $structureannotation (@$structureannotation)
	{
		my $pdb       = $structureannotation->attr('pdb');
		my $distinseq = $structureannotation->attr('distanceinsequence');
		my $euklid    = $structureannotation->attr('eukliddist');
		if ($smallestdistance)
		{
			## check if this is the smallest so far
			if ( $euklid < $smallestdistance )
			{
				$smallestdistance = $euklid;
			}
		} else
		{
			$smallestdistance = $euklid;
		}
		my $griddist = $structureannotation->attr('griddist');
		my $atopo    = $structureannotation->attr('alpha_topology');
		my $btopo    = $structureannotation->attr('beta_topology');
		my $pdbstr   = $pdb . "|" . $distinseq . "| euklid:" . $euklid . "A | grid:" . $griddist . "A |" . $atopo . "|" . $btopo . "\n";
		$html .= $pdbstr;

		#<structure_annotation alpha_topology="LYS-41-F" beta_topology="LYS-34-F" distanceinsequence="7" eukliddist="6.2" griddist="6.3" pdb="3DW8.pdb">
	}

	#$html .= '</span></a></td>';
	$html .= '"';
	return ( $html, $smallestdistance );
}

sub print_file
{
	my $filename = shift;
	my $verbose  = shift;
	my @array;
	if ($verbose) { print "Reading from file $filename\n" }
	open FILE, $filename or die $!;
	while ( my $line = <FILE> )
	{
		if ($verbose) { print "Reading line $line" }

		#chomp($line);
		#push( @array, $line );
		print $line;
	}

	#print "\n";
	close FILE;

	#return @array;
}

sub print_hash
{
	my $hashref = shift;
	foreach my $key ( keys %$hashref )
	{
		print "Protein: " . $key . " Score:" . $hashref->{$key} . "<br>";
	}
}

sub get_epvalue_url
{
	my $epvalue    = shift;
	my $epvalueurl = "displaypvalueplots.cgi?id=$id";
	my $html       = "<a target=\"_blank\" href=\"$epvalueurl\">$epvalue</a>";
	return $html;
}

sub get_xions_url
{
	my $hit               = shift;
	my $specid            = shift;
	my $rank              = shift;
	my $plottype          = shift;
	my $spectrum          = shift;
	my $specfilename      = shift;
	my $scantype          = shift;
	my $annotation        = shift;
	my $type              = $hit->attr('type');
	my $xlid              = $hit->attr('id');
	my $seq1              = $hit->attr('seq1');
	my $seq2              = $hit->attr('seq2');
	my $xlmass            = $hit->attr('xlinkermass');
	my $xlpos             = $hit->attr('xlinkposition');
	my $apriorymatchscore = $hit->attr('apriori_match_probs_log');
	my $xcorrurl          = "xions2.cgi?id=$id;plottype=$plottype;spectrum=/tmp/$spectrum;specfilename=$specfilename;type=$type;xlid=$xlid;seq1=$seq1;seq2=$seq2;xlpos=$xlpos;xlmass=$xlmass;scantype=$scantype;lapS=$apriorymatchscore";
	my $html              = "<a target=\"_blank\" href=\"$xcorrurl\">$annotation</a>";
	return $html;
}

sub print_table_header
{
}

sub print_table_row
{
}

sub print_table_footer
{
}

sub get_mods
{
	my $seq    = shift;
	my $PARAMS = shift;
	my $modid  = shift;
	my @modificationresults;

	#$debugmsg .= "Seq: $seq<br>";
	## Check the sequence for fixed and static modifications
	# Static modifications
	my @modaas = keys %{ $PARAMS->{'fixedmod'} };
	foreach my $modaa (@modaas)
	{
		my $modmass = $PARAMS->{'fixedmod'}->{$modaa};
		if ( $seq =~ /$modaa/ )
		{
			my @tmpseq1 = split //, $seq;
			my $length = scalar(@tmpseq1);

			#$debugmsg .= "Length: $length<br>";
			for ( my $i = 0 ; $i <= $#tmpseq1 ; $i++ )
			{

				#$debugmsg .= "i: $i aa: $tmpseq1[$i]<br>";
				if ( $tmpseq1[$i] eq $modaa )
				{
					my $pos = $i + 1;

					#$debugmsg .= "--> mod pos: $pos<br>";
					my $modstr = $pos . $modid . "=" . $modmass;
					push @modificationresults, $modstr;

					#$modificationresults->{$pos.$modid} = $modmass;
				}
			}
		}
	}

	# Variable modifications
	if ( $PARAMS->{'variable_mod'} )
	{
		my ( $modAA, $modmass ) = split /,|:/, $PARAMS->{'variable_mod'};
		### check if there is a modified AA (X) in the sequence
		if ( $seq =~ /X/ )
		{
			my @tmpseq1 = split //, $seq;
			for ( my $i = 0 ; $i <= $#tmpseq1 ; $i++ )
			{
				if ( $tmpseq1[$i] eq "X" )
				{
					my $pos    = $i + 1;
					my $modstr = $pos . $modid . "=" . $modmass;
					push @modificationresults, $modstr;

					#$modificationresults->{$pos.$modid} = $shift;
				}
			}
		}
	}
	return \@modificationresults;
}

sub get_absposition
{
	my $idx       = shift;
	my $xlinksite = shift;
	my $hit       = shift;
	my $seq1      = $hit->attr('seq1');
	my $seq2      = $hit->attr('seq2');
	my $prot1     = $hit->attr('prot1');
	my $prot2     = $hit->attr('prot2');
## Get the type
	my $type   = $hit->attr('type');
	my @prots1 = split( ",", $prot1 );
	my @prots2 = split( ",", $prot2 );

	#print ref($idx);
	#print "xlpos: $xlinksite, Prot1 $prot1, Prot2 $prot2, Seq1 $seq1, Seq2 $seq2<br>";
	### check if there is a variable modification defined
	if ( $PARAMS->{'variable_mod'} )
	{
		my ( $modAA, $shift ) = split /,|:/, $PARAMS->{'variable_mod'};
		### check if there is a modified AA (X) in the sequence
		if ( $seq1 =~ /X/ )
		{
			my @tmpseq1 = split //, $seq1;
			for ( my $i = 0 ; $i <= $#tmpseq1 ; $i++ )
			{
				if ( $tmpseq1[$i] eq "X" )
				{
					## modify the AA
					$tmpseq1[$i] = $modAA;
				}
			}
			### mod seq1
			$seq1 = join "", @tmpseq1;
		}
		if ( $seq2 =~ /X/ )
		{
			my @tmpseq2 = split //, $seq2;
			for ( my $i = 0 ; $i <= $#tmpseq2 ; $i++ )
			{
				if ( $tmpseq2[$i] eq "X" )
				{
					## modify the AA
					$tmpseq2[$i] = $modAA;
				}
			}
			### mod seq1
			$seq2 = join "", @tmpseq2;
		}
	}
	my @resultspep1 = ();
	my @resultspep2 = ();
	my $restring1;
	my $restring2;

	# E.g. 3,5
	# my $type      = $hit->attr('type');
	my @xlinkpositions = split( ",", $xlinksite );
	my $xlpos1         = $xlinkpositions[0];
	my $xlpos2         = $xlinkpositions[1];
	if (@prots1)
	{
		foreach my $prot1 (@prots1)
		{
			my $seqobj = $idx->fetch($prot1);
			if ($seqobj)
			{
				my $sequence = $seqobj->seq;

				#	print $sequence;
				my @res = match_sequence( $seq1, $sequence, $xlpos1 );
				$restring1 = join( "+", @res );
				push @resultspep1, $restring1;
			}
		}
		## Get the second position for intralinks
		if ( $type eq "intralink" )
		{
			foreach my $prot1 (@prots1)
			{
				my $seqobj = $idx->fetch($prot1);
				if ($seqobj)
				{
					my $sequence = $seqobj->seq;
					my @res = match_sequence( $seq1, $sequence, $xlpos2 );
					$restring2 = join( "+", @res );
					push @resultspep2, $restring2;
				}
			}
		}
	}
	if ($prot2)
	{
		foreach my $prot2 (@prots2)
		{
			my $seqobj = $idx->fetch($prot2);
			if ($seqobj)
			{
				my $sequence = $seqobj->seq;
				my @res = match_sequence( $seq2, $sequence, $xlpos2 );
				$restring2 = join( "+", @res );
				push @resultspep2, $restring2;
			}
		}
	}
	return join( ",", @resultspep1 ), join( ",", @resultspep2 );
}

sub match_sequence
{
	my $peptidesequence = shift;
	my $proteinsequence = shift;
	my $xlinksite       = shift;    # the number of the xlinked  aa in the peptide
	my @xlinksites;
	my @resxlinksites;
	chomp($peptidesequence);
	while ( $proteinsequence =~ /$peptidesequence/gi )
	{
		push @xlinksites, pos($proteinsequence);

		#print "position: $pos\n";
	}
	foreach my $site (@xlinksites)
	{
		my $numberaa = $site - length($peptidesequence) + $xlinksite;
		push @resxlinksites, $numberaa;
	}
	return @resxlinksites;
}

sub get_desc_of_id
{
	my $spidp1 = shift;
	my $spidp2 = shift;
	my @prots1 = split( ",", $spidp1 );
	my @prots2 = split( ",", $spidp2 );
	my $desc;
}

sub get_clean_id
{
	my $protid = shift;

	#print "Protein Id:$protid\n";
### Split the protein id by _ "eg: decoy_reverse_gi|147905534|ref|NP_001079812.1|"
### splits into decoy, reverse and the id: the id is always the last
	my @splittedp = split( /\_/, $protid );
## reconstitute if the part doesnt mach "reverse or decoy"
	my @reconstitute;
	foreach my $part (@splittedp)
	{

		#print "part: $part\n";
		unless ( $part =~ /decoy/ || $part =~ /reverse/ )
		{
			push @reconstitute, $part;
		}
	}
### reconstitute
	my $rec = join( "_", @reconstitute );

	#print "Reconstituted and clean id: $rec\n";
	return $rec;
}

sub get_type_of_mono_intra_link
{
	my $spidp1 = shift;
	my $type   = shift;
	my @prots1 = split( ",", $spidp1 );
	my $decoy  = 0;
	my $annotatedtype;
	foreach my $prot1 (@prots1)
	{
### Check the protein is a decoy protein
		if ( $prot1 =~ m/decoy/ )
		{
			$decoy = 1;
		}
	}
	if ($decoy)
	{
		$annotatedtype = "decoy " . $type;
	} else
	{
		$annotatedtype = $type;
	}
	return $annotatedtype;
}

sub get_type_of_xlink
{
	my $spidp1 = shift;
	my $spidp2 = shift;
	my @prots1 = split( ",", $spidp1 );
	my @prots2 = split( ",", $spidp2 );
	my $type;
	my $intralink = 0;
	my $interlink = 0;
	my $decoy     = 0;
	foreach my $prot1 (@prots1)
	{
		my $p1 = get_clean_id($prot1);

		#exit;
		foreach my $prot2 (@prots2)
		{
			my $p2 = get_clean_id($prot2);

			#print "Prot1: $prot1 , Prot2: $prot2\n";
			if ( $p1 eq $p2 )
			{
				$intralink = 1;
			} else
			{
				$interlink = 1;
			}
			### Check if one protein is a decoy protein
			if ( $prot1 =~ m/decoy/ || $prot2 =~ m/decoy/ )
			{
				$decoy = 1;
			}
		}
	}
	if ( $intralink && $interlink )
	{
		$type = "intra/inter xl";
		if ($decoy)
		{
			$type = "decoy intra/inter xl";
		}
	}
	if ( $intralink && !$interlink )
	{
		$type = "intra-protein xl";
		if ($decoy)
		{
			$type = "decoy intra-protein xl";
		}
	}
	if ( !$intralink && $interlink )
	{
		$type = "inter-protein xl";
		if ($decoy)
		{
			$type = "decoy inter-protein xl";
		}
	}
	return $type;
}

sub slice_array
{
	my $from     = shift;
	my $to       = shift;
	my $arrayref = shift;
	my @resultarray;
	my $numelements = scalar(@$arrayref);

	#print "num elements: $numelements";
	## make sure that the to is not larger than the number of total elements
	if ( $to > $numelements )
	{
		$to = $numelements - 1;
	}
	if ( $from > $numelements )
	{
		print "Error: Selection not possible!<br>";
		exit;
	}
	@resultarray = @$arrayref[ $from .. $to ];
	return \@resultarray;
}

sub get_all_ids
{
	my $resultshashref = shift;
	my @resultarray;
	foreach my $id ( keys %$resultshashref )
	{
		push @resultarray, $id;
	}
	return \@resultarray;
}

sub sort_by_type
{
	my $resultshashref   = shift;
	my $hashkeysarrayref = shift;
	my $type             = shift;
	my @resultarray;
### select all ids from the list
	foreach my $id (@$hashkeysarrayref)
	{
		my $mytype = $resultshashref->{$id}->{'hits'}->{1}->{'type'};
		unless ($mytype)
		{
			next;
			print "Error: mytype was not set found in resultsarray.<br>";

			#exit;
		}

		#print "mytype is: $mytype sorting fot $type";
		if ( $mytype eq "$type" )
		{

			#print "true";
			push( @resultarray, $id );
		}
	}
	return \@resultarray;
}

sub print_th
{
	my $arrayref = shift;
	my $string;
	foreach my $text (@$arrayref)
	{
		$string .= "<th>$text</th>";
	}
	print $string;
}

sub print_td_hash
{
	my $fields     = shift;    #HR
	my $formathash = shift;    #HR
	my $hitheader  = shift;    #AR
	my $string;
	foreach my $field (@$hitheader)
	{
		### check is there is a format string in the formathash
		if ( $formathash->{$field} )
		{
			$string .= "<td " . $formathash->{$field} . ">" . $fields->{$field} . "</td>";
		} else
		{
			$string .= "<td>" . $fields->{$field} . "</td>";
		}

		#$string .= "<td>$text</td>";
	}
	return $string;
}

sub print_td
{
	my @array = @_;
	my $string;
	foreach my $text (@array)
	{
		$string .= "<td>$text</td>";
	}
	return $string;
}

sub print_to_excel
{
	my $filehandle = shift;
	my $arrayref   = shift;
	my $string     = join( "\t", @$arrayref );
	print $filehandle $string . "\n";
}

sub make_groups
{
	my $resultshashref = shift;
	my $listref        = shift;
	my $groups         = shift;
	my @resultarray;
	foreach my $id (@$listref)
	{
		my $myid = $resultshashref->{$id}->{'hits'}->{1}->{'id'};
		## select array with ids from group
		my $array = $groups->{$myid};
		push @resultarray, @$array;
	}
	return ( \@resultarray );
}

sub generate_unique_restraints_hash_bub
{
	my $type           = shift;
	my $proteinids     = shift;    # Array with Protein IDs [P1,P2]
	my $AAposinprotein = shift;    # Array with the cross-linked AA pos in the proteins [XLpos1,XLpos2]
	my @debugarray;
	my $resulthash = {};
	my @prots1     = split( ",", $proteinids->[0] );
	my @prots2     = split( ",", $proteinids->[1] );
	my @positions1 = split( ",", $AAposinprotein->[0] );
	my @positions2 = split( ",", $AAposinprotein->[1] );
	if (@prots1)
	{

		for ( my $i = 0 ; $i < @prots1 ; $i++ )
		{

			# get the positions (if multiple positions in the prot are matched they are separated by a +)
			my @positionsP1 = split( '\+', $positions1[$i] );
			foreach my $position (@positionsP1)
			{
				my $id = $prots1[$i] . "-" . $position;
				push @{ $resulthash->{'P1'} }, $id;
			}
		}
	}
	if (@prots2)
	{
		for ( my $i = 0 ; $i < @prots2 ; $i++ )
		{

			# get the positions (if multiple positions in the prot are matched they are separated by a +)
			my @positionsP2 = split( '\+', $positions2[$i] );
			foreach my $position (@positionsP2)
			{
				my $id = $prots2[$i] . "-" . $position;
				push @{ $resulthash->{'P2'} }, $id;
			}
		}
	}
## Combine the ids, forward and reverse
	my @reshashkeys = keys %$resulthash;
	my $results     = {};
	my $idsP1       = $resulthash->{'P1'};
	my $idsP2       = $resulthash->{'P2'};
	foreach my $id1 (@$idsP1)
	{
		if ($idsP2)
		{
			foreach my $id2 (@$idsP2)
			{
				my $fwdid = $id1 . "-" . $id2;
				my $revid = $id2 . "-" . $id1;
				$results->{'fwd'}->{$fwdid} = 1;
				$results->{'rev'}->{$revid} = 1;
			}
		} else
		{
			my $fwdid = $id1;
			$results->{'fwd'}->{$fwdid} = 1;
		}
	}
	return $results;
}

sub generate_unique_restraint_ID
{
	my $type           = shift;
	my $proteinids     = shift;    # Array with Protein IDs [P1,P2]
	my $AAposinprotein = shift;    # Array with the cross-linked AA pos in the proteins [XLpos1,XLpos2] mult hits within a protein are joined with a +, different prots are separated by ,
	my @debugarray;
	my $resulthash = {};
	my @prots1     = split( ",", $proteinids->[0] );
	my @prots2     = split( ",", $proteinids->[1] );
	my @positions1 = split( ",", $AAposinprotein->[0] );
	my @positions2 = split( ",", $AAposinprotein->[1] );
	## Description
	## Usually there is one protein ID per peptide and one position in the protein sequence ID.
	## This can be complicated if the peptide matches to multiple proteins or positions within a protein.
	## Note: In the xquest viewer multiple protein matches are separated by ","
	## Multiple positions within a protein are separated by "+"; If matched to multiple proteins the separator for the matched positions is ","
	## The cross-links are sorted by the length of the sequences (P1-P2)
	## protocol for defining a unique identifier (uxID) for a restraint independent of the length of the peptide:
	## 1. Collect the protein IDs and the positions of the cross-linked residues
	## (Protein-ID string) for the individual peptides of a cross-link separately
	## - Multiple matched positions within a protein are separated by a "+" sign
	## Example:
	## Peptide 1:
	## e.g. Simple (one position): PROTEIN_A:53
	## e.g. Complex (multiple positions within P1): PROTEIN_A:53+165
	## Peptide 2:
	## e.g. Simple: PROTEIN_B:78
	## e.g. Complex (multiple positions within P1): PROTEIN_B:78+189
	## 1.b. If a peptide maps to multiple proteins collect the protein ID strings analog to 1.
	## 1.c.
	## Sort the individual protein ID strings in an ascending lexicographical order (only if 1.b is given) (function perl sort { lc($a) cmp lc($b) } )
	## and join the sorted protein ID strings using ":|:" as separator
	## 2. Sort the individual protein ID strings of the individual peptides in an ascending lexicographical order
	## and join them to a single string using ":x:" as separator
	## Restraints from monolinks and loop-links are generated analog to cross-links
	## In case of monolinks there is only 1 peptide and 1 modified residue
	## In case of looplinks there is 1 peptide but 2 residues, the residues are considered independendly as for cross-links.
	## E.g.
	## Cross-link: sp|P48607|SPZ_DROME:123:x:sp|P48607|SPZ_DROME:456
	## Loop-link:  sp|P48607|SPZ_DROME:321:x:sp|P48607|SPZ_DROME:323
	## Monolink:   sp|P48607|SPZ_DROME:181
	## 1. Sort the protein IDs by name
	my @protstrings1;
	my @protstrings2;
	if (@prots1)
	{
		for ( my $i = 0 ; $i < @prots1 ; $i++ )
		{

			# Get the positions (if multiple positions in the prot are matched they are separated by a +)
			my $protstring = $prots1[$i] . ":" . $positions1[$i];
			push @protstrings1, $protstring;
		}
		## Looplinks only have 1 protein but 2 sites
		## Generate a second string for the second xl position analog to xls
		if ( $type eq "intralink" )
		{
			for ( my $i = 0 ; $i < @prots1 ; $i++ )
			{
				# Get the positions (if multiple positions in the prot are matched they are separated by a +)
				my $protstring = $prots1[$i] . ":" . $positions2[$i];
				push @protstrings2, $protstring;
			}
		}
	}
	if (@prots2)
	{
		for ( my $i = 0 ; $i < @prots2 ; $i++ )
		{

			# get the positions (if multiple positions in the prot are matched they are separated by a +)
			my $protstring = $prots2[$i] . ":" . $positions2[$i];
			push @protstrings2, $protstring;
		}
	}
## Sort
	my @sorted_protstrings_1 = sort { lc($a) cmp lc($b) } @protstrings1;    # alphabetical sort
	my @sorted_protstrings_2 = sort { lc($a) cmp lc($b) } @protstrings2;    # alphabetical sort
## Join these
	my $joined_prot1 = join( ":|:", @sorted_protstrings_1 );
	my $joined_prot2 = join( ":|:", @sorted_protstrings_2 );
## Collect
	my @allprot;
	if ( $joined_prot1 && $joined_prot2 )
	{
		push @allprot, $joined_prot1, $joined_prot2;
	} elsif ($joined_prot1)
	{
		push @allprot, $joined_prot1;
	} else
	{
		die "  ERROR: CAnnot defined uxID";
	}
## Sort again
	my @sorted_final = sort { lc($a) cmp lc($b) } @allprot;    # alphabetical sort
## uxID string
	my $uxiD = join( ":x:", @sorted_final );
	
	# Remove whitespace if there are any
	$uxiD =~ tr/ //ds;
	
	return $uxiD;
}

sub sort_by_unique_restraints
{
	my $resultshashref = shift;
	my $listref        = shift;
	my @resultarray;
	my @debugarray,
### select all ids from the list
	  my $seen = {};
	my $restraintgroups = {};
	foreach my $id (@$listref)
	{
		my $myid = $resultshashref->{$id}->{'hits'}->{1}->{'id'};
		unless ($myid)
		{
			next;
		}
		unless ( ( $seen->{$myid} ) )
		{
			$seen->{$myid} = 1;
			push @resultarray, $id;
			$restraintgroups->{$myid} = [$id];
		} else
		{
			push( @{ $restraintgroups->{$myid} }, $id );
		}
	}
	return \@resultarray, $restraintgroups, \@debugarray;
}

sub sort_by_unique_ids
{
	my $resultshashref = shift;
	my $listref        = shift;
	my @resultarray;
	my @debugarray,
### select all ids from the list
	  my $seen = {};
	my $groups = {};
	foreach my $id (@$listref)
	{
		my $myid = $resultshashref->{$id}->{'hits'}->{1}->{'id'};
		unless ($myid)
		{
			next;
		}
		unless ( ( $seen->{$myid} ) )
		{
			$seen->{$myid} = 1;
			push @resultarray, $id;
			$groups->{$myid} = [$id];
		} else
		{
			push( @{ $groups->{$myid} }, $id );
		}
	}
	return \@resultarray, $groups, \@debugarray;
}

sub sort_xlinks_by_type_and_decoy
{
	my $hashref = shift;
	my $listref = shift;
	my $type_s   = shift;
	my $decoy   = shift;
	my @resultarray;
### select all ids from the list
	foreach my $id (@$listref)
	{
		my $type = $hashref->{$id}->{'hits'}->{1}->{'type'};
		my $xltype = $hashref->{$id}->{'hits'}->{1}->{'xltype'};
		## Check if it matches type1 or type 2
		if ( $type eq $type_s )
		{
			## if decoy is defined only use decoys
			if ( defined($decoy) )
			{
				if ( $xltype =~ m/decoy/ )
				{
					push @resultarray, $id;
				}
			} else
			{
				## exclude decoy hits
				next if ( $xltype =~ m/decoy/ );
				push @resultarray, $id;
			}
		}
	}
	return \@resultarray;
}


sub sort_xlinks_by_2types_and_decoy
{
	my $hashref = shift;
	my $listref = shift;
	my $type1   = shift;
	my $type2   = shift;
	my $decoy   = shift;
	my @resultarray;
### select all ids from the list
	foreach my $id (@$listref)
	{
		my $xltype = $hashref->{$id}->{'hits'}->{1}->{'xltype'};
		## Check if it matches type1 or type 2
		if ( $xltype =~ m/$type1/ or $xltype =~ m/$type2/ )
		{
			## if decoy is defined only use decoys
			if ( defined($decoy) )
			{
				if ( $xltype =~ m/decoy/ )
				{
					push @resultarray, $id;
				}
			} else
			{
				## exclude decoy hits
				next if ( $xltype =~ m/decoy/ );
				push @resultarray, $id;
			}
		}
	}
	return \@resultarray;
}

sub sort_xlinks_by_type_match_string
{
	my $hashref = shift;
	my $listref = shift;
	my $typestr = shift;
	my @resultarray;
### select all ids from the list
	foreach my $id (@$listref)
	{

		#print $id;
		my $xltype = $hashref->{$id}->{'hits'}->{1}->{'xltype'};
		if ( $xltype =~ m/$typestr/ )
		{
			push @resultarray, $id;
		}
	}
	return \@resultarray;
}

sub sort_xlinks_by_type_exclude_match_string
{
	my $hashref = shift;
	my $listref = shift;
	my $typestr = shift;
	my @resultarray;
### select all ids from the list
	foreach my $id (@$listref)
	{

		#print $id;
		my $xltype = $hashref->{$id}->{'hits'}->{1}->{'xltype'};
		if ( $xltype =~ m/$typestr/ )
		{
			next;
		} else
		{
			push @resultarray, $id;
		}
	}
	return \@resultarray;
}

sub sort_xlinks_by_selected
{
	my $hashref   = shift;
	my $listref   = shift;
	my $uniqueids = shift;
	
	## If unique IDs is not selected then get all selected Ids and the corresponding redundant IDs.
	my @resultarray;
	unless ( $uniqueids == 1 )
	{
		# 1. Get the Ids of the selected hits and put into seenhash
		my $seen   = {};
		my $seenid = {};
		foreach my $id (@$listref)
		{
			## Get all ranks
			my @ranks = %{ $hashref->{$id}->{'hits'} };
			foreach my $rank (@ranks)
			{
				## check if there is a validated hit
				if ( $hashref->{$id}->{'hits'}->{$rank}->{'validated'} )
				{
					my $myid = $hashref->{$id}->{'hits'}->{$rank}->{'id'};
					$seen->{$myid}++;
					$seenid->{$id}++;
					push @resultarray, $id;
				}
			}
		}
		## 2. go through all again and collect all IDs that are stored in the seenhash; use only rank 1 for the redundant ones
#		foreach my $id (@$listref)
#		{
#			my $myid = $hashref->{$id}->{'hits'}->{1}->{'id'};
#			unless ($myid)
#			{
#				next;
#			}
#			if ( ( $seen->{$myid} ) && !( $seenid->{$id} ) )
#			{
#				# push @resultarray, $id;
#			}
#		}
	} else
	{
		# Select all ids from the list
		foreach my $id (@$listref)
		{
			# get all ranks
			my @ranks = %{ $hashref->{$id}->{'hits'} };
			foreach my $rank (@ranks)
			{

				# check if there is a validated hit
				if ( $hashref->{$id}->{'hits'}->{$rank}->{'validated'} )
				{
					push @resultarray, $id;
				}
			}
		}
	}
	return \@resultarray;
}

sub sort_xlinks_by_type
{
	my $hashref = shift;
	my $listref = shift;
	my $type    = shift;
	my @resultarray;
### select all ids from the list
	foreach my $id (@$listref)
	{

		#print $id;
		my $xltype = $hashref->{$id}->{'hits'}->{1}->{'xltype'};
		if ( $xltype eq $type )
		{
			push @resultarray, $id;
		}
	}
	return \@resultarray;
}

sub sort_decoy_xlinks_by_type
{
	my $hashref = shift;
	my $listref = shift;
	my $type    = shift;
	my @resultarray;
### select all ids from the list
	foreach my $id (@$listref)
	{

		#print $id;
		my $xltype = $hashref->{$id}->{'hits'}->{1}->{'xltype'};
		if ( $xltype eq $type )
		{
			push @resultarray, $id;
		}
	}
	return \@resultarray;
}

sub sort_hash_desc
{
	my $hashref = shift;
	my $listref = shift;
	my %unsorted;
	my @resultarray;
### select all ids from the list
	foreach my $id (@$listref)
	{

		#print $id;
		my $highscore = $hashref->{$id}->{'hits'}->{1}->{'score'};
		$unsorted{$id} = $highscore;
	}
	@resultarray = sort { $unsorted{$b} <=> $unsorted{$a} } ( keys(%unsorted) );
	return \@resultarray;
}

sub filter_by_ppm
{
	my $hashref   = shift;
	my $listref   = shift;
	my $limitfrom = shift;
	my $limitto   = shift;
	my @resultarray;

	#print "Filter hash by $limit ppm<br>";
	foreach my $id (@$listref)
	{
		## get the error
		my $error = $hashref->{$id}->{'hits'}->{1}->{'error'};

		#	print "Top hit of ID: $id has error of $error ppm<br>";
		if ( ( $error >= $limitfrom ) && ( $error <= $limitto ) )
		{
			push @resultarray, $id;
		}
	}
	return \@resultarray;
}

sub filter_by_nseen
{
	my $hashref  = shift;
	my $listref  = shift;
	my $nseenmin = shift;
	my $nseenmax = shift;
	my @resultarray;
	foreach my $id (@$listref)
	{
		my $numseen = $hashref->{$id}->{'hits'}->{1}->{'nseen'};
		unless ( ( $nseenmin && $nseenmin > $numseen ) || ( $nseenmax && $nseenmax < $numseen ) )
		{
			push @resultarray, $id;
		}
	}
	return \@resultarray;
}

sub filter_by_score
{
	my $hashref  = shift;
	my $listref  = shift;
	my $minscore = shift;
	my $maxscore = shift;
	my @resultarray;
	foreach my $id (@$listref)
	{
		## get the score
		my $score = $hashref->{$id}->{'hits'}->{1}->{'score'};

		#	print $score."<br>";
		if ( $score > $minscore )
		{
			unless ( $maxscore == 0 )
			{
				if ( $maxscore < $score )
				{
					next;
				}
			}
			push @resultarray, $id;
		}
	}
	return \@resultarray;
}

sub filter_by_deltas
{
	my $hashref = shift;
	my $listref = shift;
	my $delta   = shift;
	my @resultarray;
	foreach my $id (@$listref)
	{
		## get the score
		my $deltas = $hashref->{$id}->{'hits'}->{1}->{'deltascore'};

		#	print $deltas."<br>";
		if ( $deltas < $delta )
		{
			push @resultarray, $id;
		}
	}
	return \@resultarray;
}

sub filter_by_fdr
{
	my $hashref = shift;
	my $listref = shift;
	my $delta   = shift;
	my @resultarray;
	foreach my $id (@$listref)
	{
		## get the fdr
		my $fdr = $hashref->{$id}->{'hits'}->{1}->{'fdr'};
		if ( $fdr eq "-" )
		{
			next;
		}
		if ( $fdr < $delta )
		{
			push @resultarray, $id;
		}
	}
	return \@resultarray;
}

sub filter_by_mions
{
	my $hashref        = shift;
	my $listref        = shift;
	my $minionsmatched = shift;
	my @resultarray;
	foreach my $id (@$listref)
	{
		## get the mions
		my $mions = $hashref->{$id}->{'hits'}->{1}->{'nminmatchedions'};
		unless ($minionsmatched)
		{
			next;
		}
		if ( $mions >= $minionsmatched )
		{
			push @resultarray, $id;
		}
	}
	return \@resultarray;
}

sub filter_by_xprophet_flag
{
	my $hashref        = shift;
	my $listref        = shift;
	my $minionsmatched = shift;
	my @resultarray;
	foreach my $id (@$listref)
	{
		## get the mions
		my $flag = $hashref->{$id}->{'hits'}->{1}->{'xprophet_f'};
		unless ($flag)
		{
			next;
		} else
		{
			push @resultarray, $id;
		}
	}
	return \@resultarray;
}

sub filter_by_unique_restraints
{
	my $resultshash    = shift;
	my $listref        = shift;
	my $minionsmatched = shift;
	my @resultarray;
	my $seenhash    = {};
	my $nduplicates = 0;
	my $nchecked    = 0;
	my $debugmsg;

	foreach my $id (@$listref)
	{
		## get the uxID
		my $uxID = $resultshash->{$id}->{'hits'}->{1}->{'uxID'};

		#$debugmsg .= "$uxID<br>";
		#$debugmsg.="<pre>".Dumper ($restrainthash)."<pre>";
		#$debugmsg.="<br>";
		if ( $seenhash->{$uxID} == 1 )
		{
			$nchecked++;
			$nduplicates++;
		} else
		{
			push @resultarray, $id;
			$seenhash->{$uxID} = 1;
			$nchecked++;
		}
	}
	return \@resultarray, $nduplicates, $nchecked, $debugmsg;
}

sub filter_by_deltaAA
{
	my $hashref = shift;
	my $listref = shift;
	my $min     = shift;
	my $max     = shift;
	my @resultarray;
	foreach my $id (@$listref)
	{
		## get the score
		my $dist = $hashref->{$id}->{'hits'}->{1}->{'distinsequence'};

		#	print $dist."<br>";
		unless ( ( $max && $dist > $max ) || ( $min && $dist < $min ) )
		{
			push @resultarray, $id;
		}
	}
	return \@resultarray;
}

sub filter_by_sequence
{
	my $hashref  = shift;
	my $listref  = shift;
	my $sequence = shift;
	my @resultarray;
	foreach my $id (@$listref)
	{
		## get the seq
		my $seq = $hashref->{$id}->{'hits'}->{1}->{'id'};

		#	print $dist."<br>";
		if ( $seq =~ m/$sequence/ )
		{
			push @resultarray, $id;
		}
	}
	return \@resultarray;
}

sub filter_by_annotation
{
	my $hashref    = shift;
	my $listref    = shift;
	my $annotation = shift;
	my @resultarray;
	foreach my $id (@$listref)
	{
		## get the seq
		my $anno = $hashref->{$id}->{'hits'}->{1}->{'annotation'};
		if ( $anno =~ m/$annotation/ )
		{
			push @resultarray, $id;
		}
	}
	return \@resultarray;
}

sub filter_by_spectrum
{
	my $hashref    = shift;
	my $listref    = shift;
	my $annotation = shift;
	my @resultarray;
	foreach my $id (@$listref)
	{
		## get the seq
		my $anno = $hashref->{$id}->{'hits'}->{1}->{'spectrum'};
		if ( $anno =~ m/$annotation/ )
		{
			push @resultarray, $id;
		}
	}
	return \@resultarray;
}

sub sort_hash
{
	my $hashref = shift;
	my %hash    = %$hashref;
	my %sorted;
	foreach my $key ( sort { $hashref->{$b} <=> $hashref->{$a} } keys %$hashref )
	{
		print "Key: " . $key . " Value:" . $hashref->{$key} . "<br>";
		$sorted{$key} = $hashref->{$key};
	}
	return %sorted;
}

sub debug
{

	# for debuging of param()
	if ($debug)
	{
		my @parameters = $form->param();
		print "<table border=\"1\">";
		print "<tr><td>name</td><td>value</td></tr>";
		for ( my $i = 0 ; $i < @parameters ; $i++ )
		{
			print "<tr><td>";
			print "$parameters[$i] </td>";
			print "<td>" . $form->param( $parameters[$i] ) . "</td></tr>";
		}
		print "</table>";
	}
}

sub get_prot_link
{
	my $protsstring = shift;                  ## arrayref
	my @prots = split( ",", $protsstring );
	unless (@prots)
	{
		return "-";
	}
	my $spidp1str = join( ", ", @prots );
	my @htmlarray;
	my $seen = {};
	foreach my $protein (@prots)
	{
		my $seqobj = $idx->fetch($protein);

		#$html.= Dumper ($seqobj);
		my $html;

		#http://www.uniprot.org/uniprot/
		unless ( $seen->{$protein} )
		{
			if ($seqobj)
			{
				my $ids = $seqobj->id;
				my @ids = split( /\|/, $ids );

				#$html.= Dumper (\@ids);
				my $uniplink = "http://www.uniprot.org/uniprot/$ids[1]";
				$html .= '<a align="right" class="infobox2" href=' . $uniplink . '>' . $seqobj->id;
				$html .= '<span>';
				$html .= $seqobj->desc;
			} else
			{
				my @ids = split( /\|/, $protein );
				my $uniplink = "http://www.uniprot.org/uniprot/$ids[1]";
				$html .= '<a align="right" class="infobox2" href="$uniplink">' . $protein;
				$html .= '<span>';
				$html .= 'No description found. Probably the id is not in the database.';
			}
			$html .= '</span></a>';
			$seen->{$protein}++;
			push @htmlarray, $html;
		}
	}
	my $html = join( ", ", @htmlarray );
	return $html;
}

sub filter_out_empty
{
	my $resultshashref = shift;
	my $sorted         = shift;
	my @resultarray;
### ITERATE THROUGH THE FILTERED RESULTS
	foreach my $key (@$sorted)
	{
## setn nhits to 0
		my $nhits = 0;
### Get the xml element
		my $xmlelm = ( $reshashref->{$key}->{'xml'} );
		unless ($xmlelm)
		{
			print "ERROR: NO XML Element found<br>";
			next;
		}
### Get all search hits for this spectrum
		my $subhash = $resultshashref->{$key}->{'hits'};
		my @allhitssorted = sort { $a <=> $b } keys %$subhash;
		for ( my $i = 0 ; $i < ( scalar(@allhitssorted) ) ; $i++ )
		{
			my $hit       = $allhitssorted[$i];
			my $score     = $resultshashref->{$key}->{'hits'}->{$hit}->{'score'};
			my $error_rel = $resultshashref->{$key}->{'hits'}->{$hit}->{'error'};
			my $id        = $resultshashref->{$key}->{'hits'}->{$hit}->{'structure'};
############# FILTERS ############
### CHECK IF FILTERING OF ERROR IS ENABLED ON ALL HITS, IF YES FILTER
			if ( ($filter_error_on_all_hits) )
			{
				if ( $error_rel < $maxerrorfrom || $error_rel > $maxerrorto )
				{
					next;
				} else
				{
					$nhits++;
				}
			}
		}
		## Push into hash if there are still hits
		if ( $nhits > 0 )
		{
			push @resultarray, $key;
		}
	}
	return \@resultarray;
}

#---------------------------------------------------------------------------
#  SUB parseXML for PARSING THE XML FILE
#---------------------------------------------------------------------------
sub parse_spec_xml_file
{
	my $xmlfilename = shift;
	my $outfilename = $xmlfilename . ".hash";
	if ( -e $outfilename )
	{
		return;
	}
	my $tree = XML::TreeBuilder->new();
	$tree->parse_file($xmlfilename);
	my @resultsheader = map { $_ } $tree->find('xquest_spectra');
	### index the spectra by filename
	my $spechash = {};
	foreach my $header (@resultsheader)
	{

		#	print "Parsing Spectrum xml file\n";
		my @spectra = map { $_ } $header->find('spectrum');
		foreach my $spectrum (@spectra)
		{
			my $filename = $spectrum->attr('filename');

			#my $type = $spectrum->attr('type');
			my $content = $spectrum->content();

			#$spechash->{$filename}=$spectrum->as_XML;
			$spechash->{$filename} = $content->[0];
		}
	}
	store $spechash, $outfilename;
	$tree->delete();
}

sub parseXML
{
	my $xmlfilename   = shift;
	my $resultshash   = shift;
	my $idx           = shift;
	my $resultdirbase = shift;
	my $id            = shift;
### Parse the XML file
	unless ( -e $xmlfilename )
	{
		die "Cannot find xquest result file $xmlfilename\n";
	}
	my $tree = XML::TreeBuilder->new();
	$tree->parse_file($xmlfilename);

	#print "Parsing xquest XML file<br>";
### Parsing all headers // one or more if a it is a merged result
	my @resultsheader = map { $_ } $tree->find('xquest_results');
## the running id for spectra start with 1 (0 is undef)
	my $spectrumid = 1;
## Counting of redundant hits
	my $seen          = {};
	my $seenannotaion = {};
## Parsing all specrum search results
	foreach my $header (@resultsheader)
	{
		my $headerid = $header->attr('outputpath');
		## Define the filename of the spectrumxml file
		my $specfn               = $headerid . ".spec.xml";
		my $exportspectrafromxml = 0;
		if ($exportspectrafromxml)
		{
			### check if there are spectra in there (old verision)
			my @spectra = map { $_ } $header->find('spectrum');
			if (@spectra)
			{

				#print "generate spec.hash file\n";
				my $outfilename = $specfn . ".hash";
				my $specxmlfile = File::Spec->catfile( $resultdirbase, $id, $outfilename );
				my $spechash    = {};
				foreach my $spectrum (@spectra)
				{
					my $filename = basename( $spectrum->attr('filename') );
					my $content  = $spectrum->content();

					#$spechash->{$filename}=$spectrum->as_XML;
					$spechash->{$filename} = $content->[0];
					$spectrum->delete_content();
				}
				## Store it as hash
				#print "Store hash $specxmlfile\n";
				store $spechash, $specxmlfile;
			}
		}
		my $specxmlfile = File::Spec->catfile( $resultdirbase, $id, $specfn );
		if ( -e $specxmlfile )
		{
			## parse the spectra into a separate hash and write a hashfile
			parse_spec_xml_file($specxmlfile);
		}
		## parsing of all spectrum search results
		my @spectrumsearchelements = map { $_ } $header->find('spectrum_search');
		my $num = @spectrumsearchelements;

		#print "Found header element: $headerid with $num spectrum_search elements<br>";
		foreach my $spectrumsearchelements (@spectrumsearchelements)
		{
			$resultshash->{$spectrumid}->{'xml'} = $spectrumsearchelements->as_XML;
			## save also the spectum identifier assocoated with the spectrum id in the hash
			my $specidentifier = $spectrumsearchelements->attr('spectrum');
			## get rt and mz
			#my $rts = $spectrumsearchelements->attr('rtsecscans');
			#my $mzs = $spectrumsearchelements->attr('mzscans');
			### Get all search hits for this spectrum
			my @search_hits = map { $_ } $spectrumsearchelements->find('search_hit');
			my $i = 1;
			foreach my $hit (@search_hits)
			{
				### Get the rank
				my $rank = $hit->attr('search_hit_rank');
				### ONLY INDEX UNTIL rank 5
				#if ( $rank > 5 )
				#{
				#	next;
				#}
				my $score         = $hit->attr('score');
				my $error         = $hit->attr('error_rel');
				my $type          = $hit->attr('type');
				my $id            = $hit->attr('id');
				my $annotation    = $hit->attr('annotation');
				my $structure     = $hit->attr('structure');
				my $fdr           = $hit->attr('fdr');
				my $nmionsalpha   = $hit->attr('num_of_matched_ions_alpha');
				my $nmionsbeta    = $hit->attr('num_of_matched_ions_beta');
				my $xprophet_flag = $hit->attr('xprophet_f');
				$resultshash->{$spectrumid}->{'hits'}->{$rank}->{'score'}           = $score;
				$resultshash->{$spectrumid}->{'hits'}->{$rank}->{'error'}           = $error;
				$resultshash->{$spectrumid}->{'hits'}->{$rank}->{'type'}            = $type;
				$resultshash->{$spectrumid}->{'hits'}->{$rank}->{'id'}              = $id;
				$resultshash->{$spectrumid}->{'hits'}->{$rank}->{'specfile'}        = $specfn;
				$resultshash->{$spectrumid}->{'hits'}->{$rank}->{'structure'}       = $structure;
				$resultshash->{$spectrumid}->{'hits'}->{$rank}->{'nseen'}           = 1;
				$resultshash->{$spectrumid}->{'hits'}->{$rank}->{'fdr'}             = $fdr;
				$resultshash->{$spectrumid}->{'hits'}->{$rank}->{'spectrum'}        = $specidentifier;
				$resultshash->{$spectrumid}->{'hits'}->{$rank}->{'xprophet_f'}      = $xprophet_flag;
				$resultshash->{$spectrumid}->{'hits'}->{$rank}->{'nminmatchedions'} = get_n_min_ions_matched( $type, $nmionsalpha, $nmionsbeta );
				### check if the hit has an annotation
				if ($annotation)
				{
					$resultshash->{$spectrumid}->{'hits'}->{$rank}->{'annotation'} = $annotation;
				}
				###
				# check if seen on this rank
				if ( $seen->{$rank}->{$id} )
				{
					### it has been seen already
					### add one to the counter
					$seen->{$rank}->{$id}->{'count'}++;
					## add this hit to the seen ids
					push @{ $seen->{$rank}->{$id}->{'specids'} }, $spectrumid;
					## update the counts
					foreach my $specid ( @{ $seen->{$rank}->{$id}->{'specids'} } )
					{
						$resultshash->{$specid}->{'hits'}->{$rank}->{'nseen'} = $seen->{$rank}->{$id}->{'count'};
					}
					### now update all entries ids where it has been seen
				} else
				{
					### it has not been seen so far.
					$seen->{$rank}->{$id}->{'count'} = 1;
					push @{ $seen->{$rank}->{$id}->{'specids'} }, $spectrumid;
				}
				## Add the type of xlink
				my $spidp1 = $hit->attr('prot1');
				my $spidp2 = $hit->attr('prot2');
				### Determine the type of cross-link
				my $xlinktype = "-";
				if ( $type eq "xlink" )
				{
					$xlinktype = get_type_of_xlink( $spidp1, $spidp2 );
				} else
				{
					## get the type of mono or intra-links
					$xlinktype = get_type_of_mono_intra_link( $spidp1, $type );
				}
				$resultshash->{$spectrumid}->{'hits'}->{$rank}->{'xltype'} = $xlinktype;
				### DETERMINE THE POSITIONS
				### Determine the absolute positions of the cross-linked residues
				my $xlinkpos = $hit->attr('xlinkposition');
				my ( $res1, $res2 ) = get_absposition( $idx, $xlinkpos, $hit );
				my @aa1 = split( ",", $res1 );
				my @aa2 = split( ",", $res2 );
				my $distinsequence = "-";
				if ( $type eq "xlink" )
				{
					$distinsequence = abs( $aa1[0] - $aa2[0] );
				}
				if ( $type eq "intralink" )
				{
					$distinsequence = abs( $aa1[0] - $aa2[0] );
				}
				$resultshash->{$spectrumid}->{'hits'}->{$rank}->{'distinsequence'} = $distinsequence;
				## Generate restraint IDs --> arrayrefs fwd IDs and rev IDS
				## xlinks P1-P2-Pos1-Pos2
				## monolinks:
				## Define the restraints hash
				my $uxID = generate_unique_restraint_ID( $type, [ $spidp1, $spidp2 ], [ $res1, $res2 ] );
				$resultshash->{$spectrumid}->{'hits'}->{$rank}->{'uxID'} = $uxID;
				$i++;
			}
			$spectrumid++;
			$spectrumsearchelements->delete();
		}
	}
	$tree->delete();

	#exit;
}

#---------------------------------------------------------------------------
#  EOF xml parsing
#---------------------------------------------------------------------------
sub get_n_min_ions_matched
{
	my $type           = shift;
	my $nmionsalpha    = shift;
	my $nmionsbeta     = shift;
	my $minionsmatched = 0;
	if ( $type eq "xlink" )
	{
		if ( $nmionsalpha == $nmionsbeta )
		{
			$minionsmatched = $nmionsalpha;
		} elsif ( $nmionsalpha < $nmionsbeta )
		{
			$minionsmatched = $nmionsalpha;
		} elsif ( $nmionsbeta < $nmionsalpha )
		{
			$minionsmatched = $nmionsbeta;
		}
	} else
	{
		$minionsmatched = $nmionsalpha;
	}
	return $minionsmatched;
}

sub setparams
{
	print $form->hidden( -name  => 'id',
						 -value => $form->param('id'), );

	#	print $form->hidden(
	#		-name  => 'maxrank',
	#		-value => $form->param('maxrank'),
	#	);
	print $form->hidden( -name  => 'reload',
						 -value => 1, );
}

sub gethighscore
{
	my $specsearch = shift;
	my @results    = $specsearch->find('search_hit');
	my $tophit     = $results[0];
	if ( defined($tophit) )
	{
		return $tophit->attr('score');
	} else
	{
		return 0;
	}
}

sub getxlinktype
{
	my $specsearch = shift;
	my @results    = $specsearch->find('search_hit');
	my $tophit     = $results[0];
	if ( defined($tophit) )
	{
		return $tophit->attr('type');
	} else
	{
		return 0;
	}
}

sub getxlinktypename
{
	my $specsearch = shift;
	my $type       = shift;

	#print $type;
	my @results = $specsearch->find('search_hit');
	my $tophit  = $results[0];
	print $tophit->attr('type');
	if ( defined($tophit) )
	{

		#print "true";
		if ( $type eq $tophit->attr('type') )
		{
			return $tophit->attr('type');
		}
	} else
	{
		return 0;
	}
}

sub getfirstname
{
	my $specsearch = shift;
	my @results    = $specsearch->find('search_hit');
	my $tophit     = $results[0];
	if ( defined($tophit) )
	{
		return $tophit->attr('id');
	} else
	{
		return 0;
	}
}

sub printsortselection
{
	print '<table class="greybox" border="thin" bordercolor="#F0F0F0" width="780"   >';
	print '<tr>';
	print '<td>';
	my %sortlabels = (
					   'byscore' => 'score',
					   'byid'    => 'id',
					   'byspec'  => 'spectrum name',
					   'bytype'  => 'cross-link type',
	);
	print "sort by: ",
	  radio_group(
				   -name    => 'sortby',
				   -values  => [ 'byscore', 'byid', 'byspec', 'bytype' ],
				   -default => 'byscore',
				   -labels  => \%sortlabels,
	  );
	print "    show:   ",
	  popup_menu(
				  {
					name    => "show",
					values  => [ "all", "cross-links", "mono-links", "intra-links" ],
					default => ["all"]
				  }
	  );
	print '</td>';
	print '</tr>';
	print '</table>';
}

sub printhitselection
{
	print '<table class="greybox" border="thin" bordercolor="#F0F0F0" width="780"   >';
	print '<tr>';
	print '<td>';
	my %labels = (
				   'tophit'          => 'show top hit only',
				   'selectedhit'     => 'show all hits for selected spectra',
				   'allhits'         => 'show all hits for all spectra',
				   'selectedhits'    => 'show only selected hits',
				   'nonselectedhits' => 'show only non-selected hits',
	);
	print '</tr>';
	print '<tr>';
	print '<td>';
	print checkbox( -name => 'excel', -label => 'create Excel file' );
	print '</td>';
	print '<td>';
	print "show scores > ";
	print textfield(
					 -name      => 'minscore',
					 -value     => 0,
					 -size      => 2,
					 -maxlength => 5
	);
	print " show rank1 to ";
	print textfield(
					 -name      => 'maxrank',
					 -value     => 3,
					 -size      => 2,
					 -maxlength => 2
	);
	print '</td>';
	print '</tr>';
	print '<tr>';
	print '<td>';
	print submit( -name => 'sortby', -label => 'refresh' );
	print '</td>';
	print '</tr>';
	print '</table>';
}

sub print_header
{
	print $form->header('text/html');
	print $form->start_html(
							 -title   => 'xQuest/xProphet results viewer',
							 -author  => 'walzthoeni@imsb.biol.ethz.ch',
							 -base    => 'true',
							 -expires => '+30s',
							 -meta    => { 'keywords' => 'Ms cross-link Walzthoeni Rinner' },
							 -style   => { 'src' => $WEBPARAMS->{'css_stylesheet'}, },
	);
	print <<EOF;
<style type="text/css">
  a.infobox { border-bottom: 1px; text-decoration:none; }
  a.infobox:hover { color:#c30; background:white; }
  a.infobox span { visibility:hidden; position:absolute; left:-99em;
    margin-top:1.5em; padding:1em; text-decoration:none; }
  a.infobox:hover span, a.infobox:focus span, a.infobox:active span {
    visibility:visible; left:1em;
    border:1px solid #c30; color:blue; background:white; }
  
  a.infobox2 { color:black; border-bottom: 1px; text-decoration:none; }
  a.infobox2:hover { color:#c30; background:white; }
  a.infobox2 span { visibility:hidden; position:absolute; left:-99em;
    margin-top:1.5em; padding:1em; text-decoration:none; }
  a.infobox2:hover span, a.infobox:focus span, a.infobox:active span {
    visibility:visible; left:1em;
    border:1px solid #c30; color:blue; background:white; }
</style>
<!--[if IE 5]><style type="text/css">
  a.infobox span { display:none; }
  a.infobox:hover span { display:block; }
</style><![endif]-->
EOF
}

sub printrestart_button
{
	print '<input  type="submit" value="new search"></form>';
}

sub readwebparams
{
	my $webparam = shift;
	my %WEBPARAMS;
	open WEBCONFIG, "<$webparam"
	  or warn "could not open web config file $! ignoring";
	while (<WEBCONFIG>)
	{
		chomp;
		my @keyvalue = split /::/;
		$WEBPARAMS{ $keyvalue[0] } = $keyvalue[1];
	}
	return \%WEBPARAMS;
}

sub readtables
{
	my $webparam = shift;
	my %WEBPARAMS;
	open WEBCONFIG, "<$webparam"
	  or warn "could not open web config file $! ignoring";
	while (<WEBCONFIG>)
	{
		chomp;
		my @keyvalue = split /::/;
		$WEBPARAMS{ $keyvalue[0] } = $keyvalue[1];
	}
	return \%WEBPARAMS;
}

sub read_params
{
	my $deffile = shift;
	my $PARAMS  = shift;
	open DEFFILE, "<$deffile" or die "cannot open table $deffile $!";
	while ( my $line = <DEFFILE> )
	{
		chomp($line);
		if ($line)
		{
			my @results = split( " ", $line );
			$PARAMS->{ $results[0] } = $results[1];
		}
	}
	return $PARAMS;
}

sub read_def
{
	my $deffile = shift;
	my %PARAMS;
	open DEF, "<$deffile" or die "cannot open table $deffile $!";

	#read in definitions
	#read in definitions
	my $configdef     = undef;
	my $enzymedef     = undef;
	my $modifications = undef;
	while ( my $line = <DEF> )
	{
		chomp($line);

		#my $line = $_;
		if ( $line =~ /^digestdef/ )
		{
			$configdef     = defined;
			$enzymedef     = undef;
			$modifications = undef;
			next;
		} elsif ( $line =~ /^enzymedef/ )
		{
			$configdef     = undef;
			$modifications = undef;
			$enzymedef     = defined;
			next;
		} elsif ( $line =~ /^modifications/ )
		{
			$configdef     = undef;
			$enzymedef     = undef;
			$modifications = defined;
			next;
		}
		if ( defined($configdef) )
		{
			if ($line)
			{
				my @results = split( ' ', $line );
				$PARAMS{ $results[0] } = $results[1];
			}
		}
		## Only the fixed mods
		if ( defined($modifications) )
		{
			if ($line)
			{
				my @results = split( ' ', $line );
				if ( $results[1] )
				{
					$PARAMS{'fixedmod'}{ $results[0] } = $results[1];
				}
			}
		}
	}
	close DEF;
	return \%PARAMS;
}

sub debug_param
{

	#my $var=self;
	my $var = shift;

	# for debuging of param()
	my @parameters = $var->param();
	@parameters = sort (@parameters);
	print "<table border=\"1\">";
	print "<tr><td>name</td><td>value</td></tr>";
	for ( my $i = 0 ; $i < @parameters ; $i++ )
	{
		print "<tr><td>";
		print "$parameters[$i] </td>";
		print "<td>" . $var->param( $parameters[$i] ) . "</td></tr>";
	}
	print "</table>";
	return;
}
