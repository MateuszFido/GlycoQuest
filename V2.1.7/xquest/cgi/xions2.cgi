#!C:/Perl/bin/perl.exe
use strict;
#---------------------------------------------------------------------------
# xions2.cgi
# A software/script to display xQuest MS/MS spectra.
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
use CGI::Carp 'fatalsToBrowser';
use CGI::FastTemplate;
use Socket;
use File::Basename;
use File::Spec;
use File::Path;
#use Mail::Sender;
use XML::Parser;
use XML::TreeBuilder;

# include the directory for the xquest modules ###########
use File::Spec::Functions qw(rel2abs);
use File::Basename;
use FindBin;
use lib "$FindBin::Bin/../modules";
##########################################################
use Environment;
use PepObj;
use Match;
use Spectrum;
use Xcorrelation;
use LinkObj;
use Data::Dumper;
use Storable;
use MIME::Base64;
use Read_Params;

# in silico matching of peptide to spectrum
# sub getspecobj creates a spectrum object from the spectrum
# sub print_table does all the matching
my $env       = Environment->new();
my $webconfig = $env->get_path('web.config');
my $form      = new CGI;
my $myself    = self_url();
my $time      = localtime();
### PARSE THE PARAMS FROM THE URL
my $sessionid        = $form->param('id');
my $spectrumfilename = $form->param('spectrum');
my $specfilename     = $form->param('specfilename');
my $xltype           = $form->param('type');
my $xlid             = $form->param('xlid');
my $seq1             = $form->param('seq1');
my $seq2             = $form->param('seq2');
my $xlpos            = $form->param('xlpos');
my $xlmass           = $form->param('xlmass');
my $scantype         = $form->param('scantype');
my $apScore 		 = $form->param('lapS');
my $showScores 		 = $form->param('showscores');
### Build a hash of the searchhit
my $searchhit={};
$searchhit->{'seesionid'}=$sessionid;
$searchhit->{'spectrum'}=$spectrumfilename;
$searchhit->{'specfilename'}=$specfilename;
$searchhit->{'type'}=$xltype;
$searchhit->{'xlid'}=$xlid;
$searchhit->{'seq1'}=$seq1;
$searchhit->{'seq2'}=$seq2;
$searchhit->{'xlpos'}=$xlpos;
$searchhit->{'xlmass'}=$xlmass;
$searchhit->{'scantype'}=$scantype;
$searchhit->{'lapS'}=$apScore;

### Other not used params
my $hitid         = $form->param('hitid');
my $proteinids    = $form->param('proteins');
my $xmlfilename   = $form->param('resultxml');
my $normxcorr     = $form->param('norm_xcorr');
my $urlbase       = $form->param('urlbase');
my $plottype      = $form->param('plottype');
my $xquestdeffile = $form->param('xquestdef');
### WEBPARAMS
my $WEBPARAMS      = readwebparams($webconfig);
my $resultdirbase  = $WEBPARAMS->{'resultdirbase'};
my $resulturlbase  = $WEBPARAMS->{'resulturlbase'};
my $xlinkionparser = $WEBPARAMS->{'xlinkionparser'};
my ( $waterloss, $pngurl, $lightspectrum );
print_header();
print '<h1>xQuest results</h1><h3>', $xlid,      '</h3>';
#print "Specfilename: $specfilename<br>";
my $debug = 0;

if ($debug)
{
debug_param($form);
}
print $form->start_form;

my $xquestdir;
if ( defined( $ENV{'XQUEST_DIR'} ) )
{
	$xquestdir = $ENV{'XQUEST_DIR'};
} else
{
	## the cgi-bin is the xquest dir
	#warn "environment variable XQUEST_DIR is not defined $!, setting it to .";
	$xquestdir = ".";
}

#### LOAD THE CONFIG FILES ####
#open masstable and definition table
my $masstable = File::Spec->catfile( $xquestdir, "mass_table.def" );
my $masstable = $env->get_path('mass.def');
my $xquestdefcopy = File::Spec->catfile( $resultdirbase, $sessionid, "xquest.def" );
my $masslist;
my ( $MSTAB, $PARAMS, $ENZ, $basename ) = Read_Params::readtables( $masstable, $xquestdefcopy, $webconfig, $masslist,$xquestdir,0,$sessionid,$resultdirbase );

##### store posted parameters ######
setparams();

if ($sessionid)
{
	$xmlfilename = File::Spec->catfile( $resultdirbase, $sessionid, 'xquest.xml' );
	my $pngname = File::Spec->catfile( $resultdirbase, $sessionid, 'tmp.png' );
	$spectrumfilename = File::Spec->catfile( $resultdirbase, $sessionid, $form->param('spectrum') );
	$urlbase = join "\/", $resulturlbase, $sessionid;
	$lightspectrum = join "",   $spectrumfilename, "_common.txt";
} else
{
	die "no sessionid $!";
}

my $tmpdir = File::Spec->catfile( $resultdirbase, $sessionid, "tmp" );
### Create tmp folder for spectra and imgs
unless ( -e $tmpdir )
{
	mkdir($tmpdir) or die "cannot create tmp directory $tmpdir $!";
}
my @pngfiles = glob("$tmpdir/*.png");
my @dtafiles=glob("$tmpdir/*.dta");
my @txtfiles=glob("$tmpdir/*.txt");

## Deletes all png/dta/txt files so that they dont pile up in the tmp dir
foreach my $file (@pngfiles)
{
	#print "Will delete file: $file<br>";
	unlink($file);
}
foreach my $file (@dtafiles)
{
	#print "Will delete file: $file<br>";
	unlink($file);
}
foreach my $file (@txtfiles)
{
	#print "Will delete file: $file<br>";
	unlink($file);
}


#print "Plottype: $plottype\n";
print_ionselect( $plottype, $PARAMS );
print p, submit( -name => 'reload' );

#my $xmlurl = join "\/", $resulturlbase, $sessionid, 'xquest.xml';
unless ( -e $xmlfilename )
{
	print "<em>results not ready. Check back later (you can bookmark this page)<em>";
	print '<form method="get" action="javascript:bookmark()">';
	print '<input type="submit" value=" bookmark "></form>';
	exit 0;
}

#---------------------------------------------------------------------------
#  Read the spectrumfilehash
#---------------------------------------------------------------------------
### Initialize the XML
### Use the hash and pull out the searchhit as xml
##  Check if a db is already there
my $filename = File::Spec->catfile( $resultdirbase, $sessionid, 'resultdatabase.hash' );
### Filename of spectrumfile
my $specfilenamehash = File::Spec->catfile( $resultdirbase, $sessionid, "$specfilename.hash" );
#print $specfilenamehash. "<br>";
my $dbfilesemaphore;
my $hashref;
if ( -e $specfilenamehash && -r $specfilenamehash )
{
	$dbfilesemaphore = 1;
## reads the Hash from the DB
	$hashref = retrieve($specfilenamehash);

	#%resultshash = %$hashref;
	#print "Result hash found<br>";
	#print Dumper($hashref);
} else
{

	#Todo: parse the xml
	die "Cannot find resultshash $specfilenamehash $! <br>";
}

#---------------------------------------------------------------------------
#  Get & Store the Spectra in the resultfolder
#---------------------------------------------------------------------------
my @specfilenames;
#print "Specfilename: $spectrumfilename<br>";
## split the filename, get splitsites: e.g.: wathomas_M1111_179.c.02974.02974.4_wathomas_M1111_179.c.02776.02776.4
my @cutsites;
my $bnfn=basename ($spectrumfilename);

while ( $bnfn =~ /_/gi ) {
push @cutsites, pos($bnfn);
}

#print Dumper (@cutsites);

## substring CUT STRING IN THE MIDDLE _
my $cutat;
my $numcutsites=@cutsites;
for (my $i=1; $i<$numcutsites+1; $i++){
#print "Checking $i of $numcutsites cutsites<br>";
if (($i*2>$numcutsites) && !($i>($numcutsites+1)/2) ){
$cutat=$i;
#print "Cutsite to split:$i<br>";	
}
}


my $cut1=$cutsites[$cutat-1];
my $string1 = substr($bnfn, 0,$cut1-1);
my $string2 = substr($bnfn, $cut1);


my $light_specname   = $string1 . ".dta";
my $heavy_specname   = $string2 . ".dta";

#print $light_specname."<br>";
#print $heavy_specname."<br>";
my $commonspecfilename  = basename($spectrumfilename) . "_common.txt";
my $xlinkerspecfilename = basename($spectrumfilename) . "_xlinker.txt";
push @specfilenames, $light_specname, $heavy_specname, $commonspecfilename, $xlinkerspecfilename;
foreach my $spectrum (@specfilenames)
{
	my $filename = File::Spec->catfile( $resultdirbase, $sessionid, "tmp", $spectrum );
	my $content  = $hashref->{$spectrum};
	my $pl       = decode_base64($content);

	#print "Decoded: ".$pl."<br>";
	unless ( -e $filename )
	{
		#print "Write file $filename<br>";
		save_to_file( $filename, $pl );
	}
}
#exit;
# $hitid = "LCVLHEKTPVSEK-CASIQKFGER-a7-b6";
### HIT ID MUST BE DEFINED
my $addedmass;

#$addedmass = $specsearch->attr('addedMass');
$addedmass = 0;
#print Dumper($searchhit);
#exit;
## Function where spectra is printed
#print "nh3loss:".$PARAMS->{'nh3loss'}."<br>";
print_table( $searchhit, $addedmass );
printtail();
#$searchhit->delete;

#$tree->delete;
sub printtail
{
	print '</div>';
	print $form->end_html;
	print $form->end_form;
}

#===  FUNCTION  ================================================================
#  NAME:        save_to_file
#  PURPOSE:     save_to_file
#  DESCRIPTION: Save a string to a file
#  PARAMETERS:  $filename, $string, $append (1 if string should be appended)
#  RETURNS:     1 if done correctly
#===============================================================================
sub save_to_file
{
	my $filename = shift;
	my $text     = shift;
	my $append   = shift;
	if ($append)
	{
		open MYOUTFILE, ">>", "$filename"
		  or die "cannot open file $filename $!";
	} else
	{
		open MYOUTFILE, ">>", "$filename"
		  or die "cannot open file $filename $!";
	}
	print MYOUTFILE $text . "\n";
	close MYOUTFILE;
	return 1;
}

sub print_ionselect
{
	my $plottype = shift;
	my $PARAMS   = shift;

	### Default parameters
	### If the page is loaded the first time (no reload was clicked)
	### then get the default parameters form the xquest def
	
	my ($matchsecond,$waterloss,$nhloss);
	unless ($form->param('reload')){
	$matchsecond=$PARAMS->{'CID_match2ndisotope'};	
	$waterloss=$PARAMS->{'waterloss'};
	$nhloss=$PARAMS->{'nh3loss'};
		
	}



	printrangeselect( $plottype, $PARAMS->{'minionsize'}, $PARAMS->{'maxionsize'}, $PARAMS->{'ms2tolerance'}, $PARAMS->{'xlink_ms2tolerance'}, $PARAMS->{'xcorrprecision'}, $PARAMS->{'xcorrdelay'} );
	unless ( $plottype eq "xcorr" )
	{
		print_logscale_select();
		print "     ", checkbox( -name => 'labels', );
		print "     ", checkbox( -name => 'hide peptide structure' );
		print '<br>';
	}
	print hr, checkbox(
		-name  => 'waterloss',
		-label => 'water-loss',
		-value => 1,
		-checked => $waterloss,
	);
	print "     ", checkbox(
		-name  => 'nh3loss',
		-label => 'NH3-loss',
		-value => 1,
		-checked => $nhloss,
	);
	unless ( $plottype eq "xcorr" )
	{
		print "     ", checkbox(
			-name  => 'CID_match2ndisotope',
			-label => 'allow matching of 2nd isotopic peak',		
			-value =>1,
			-checked => $matchsecond,
		);
		print "     ",
		  checkbox( -name  => 'hittable',
					-label => 'show table with matching ions', );
		print "     ",
		  checkbox( -name  => 'spectrumcomparison',
					-label => 'show comparison of original spectra', );
	}
	print checkbox_group(
						  -name    => 'ionseries',
						  -values  => [ 'a', 'b', 'c', 'x', 'y', 'z' ],
						  -default => $PARAMS->{'ionseries_array'},
						  -rows    => 3,
						  -columns => 2
	);
}

sub printrangeselect
{
	my $plottype          = shift;
	my $min               = shift;
	my $max               = shift;
	my $ms2tolerance      = shift;
	my $xlinkms2tolerance = shift;
	my $xcorrbinsize      = shift;
	my $xcorrdelay        = shift;
	if ( $plottype eq "xcorr" )
	{
		print "  x-corr bin-size [Da]";
		print textfield(
						 -name      => 'xcorrprecision',
						 -value     => $xcorrbinsize,
						 -size      => 6,
						 -maxlength => 6
		);
		print " x-corr delay [Da]";
		print textfield(
						 -name      => 'plotxcorrdelay',
						 -value     => $xcorrdelay,
						 -size      => 6,
						 -maxlength => 6
		);
		print "<br>";
	} else
	{
		print "mz min. ";
		print textfield(
						 -name      => 'minionsize',
						 -value     => $min,
						 -size      => 6,
						 -maxlength => 6
		);
		print "   mz max. ";
		print textfield(
						 -name      => 'maxionsize',
						 -value     => $max,
						 -size      => 6,
						 -maxlength => 6
		);
		print "<br>MS2 tolerance for common-ions [Da]";
		print textfield(
						 -name      => 'ms2tolerance',
						 -value     => $ms2tolerance,
						 -size      => 6,
						 -maxlength => 6
		);
		print "     MS2 tolerance for xlink-ions [Da]";
		print textfield(
						 -name      => 'xlink_ms2tolerance',
						 -value     => $xlinkms2tolerance,
						 -size      => 6,
						 -maxlength => 6
		);
	}
}

sub print_logscale_select
{
	print "     ", checkbox( -name => 'logscale', );
}

sub print_header
{
	print $form->header('text/html');
	print $form->start_html(
							 -title   => 'xQuest ionseries',
							 -author  => 'rinner@xquest.org',
							 -base    => 'true',
							 -pragma  => "no-cache",
							 -expires => 'now',
							 -meta    => {
										'keywords'  => 'ms cross-link rinner',
										'copyright' => 'copyright 2006 Oliver Rinner'
							 },
							 -style => { 'src' => $WEBPARAMS->{'css_stylesheet'}, },
	);
}

sub print_table
{
	my $searchhit = shift;
	my $plusMr    = shift;
	my $xquestdir;
	if ( defined( $ENV{'XQUEST_DIR'} ) )
	{
		$xquestdir = $ENV{'XQUEST_DIR'};
	} else
	{
		$xquestdir = ".";
	}
	my $masstable = File::Spec->catfile( $xquestdir, "mass_table.def" );
	my $webconfig = File::Spec->catfile( $xquestdir, "web.config" );

	#open masstable and definition table
	#my ( $MSTAB, $PARAMS, $ENZ, $basename, $WEBPARAMS ) = readtables( $masstable, $xquestdefcopy, $webconfig );

	#my $matchions = $searchhit->find('match_ions');
	my $matchions;
	if ( $plottype eq "spectrum" )
	{
		printxlinkform( $matchions, $form->param('minionsize'), $form->param('maxionsize'), $form->param('labels'), $plusMr );
	} elsif ( $plottype eq "xcorr" )
	{
		printxlinkxcorr( $matchions, $PARAMS, $xquestdefcopy, $plusMr );
	}
}

sub printxlinkform
{
	my $matchions = shift;
	my $min       = shift;
	my $max       = shift;
	my $labels    = shift;
	my $plusMr    = shift;
	my $specname  = shift;
	my $hithash;

	my %buffer;
	
	my $specobj = getspecobj( $spectrumfilename, $xquestdefcopy, \%buffer );
	## The xlink object is a LinkObj
	## The function getxlinkions initializes and returns a xlinkobject!
	my $xlinkobj = getxlinkions( $searchhit->{'type'}, $searchhit->{'seq1'}, $searchhit->{'seq2'}, $searchhit->{'xlmass'}, $searchhit->{'xlpos'}, $xquestdefcopy, undef, $specobj, $plusMr );
#print "Lightspectrum: $lightspectrum<br>";
	if ( -e $lightspectrum )
	{
		my $commonpeaks = $xlinkobj->getcommonpairs;
		my $xlinkpeaks  = $xlinkobj->getxlinkpairs;
		## offline matching ##
		#$xlinkobj->match_offline($specobj);
		## only for testing ###
		## make a matchobject dummy
		my $matchobjdummy=Match->new_matchobj_dummy($specobj);
		$xlinkobj->{'matchobj'}=$matchobjdummy;
		#print "APS", $searchhit->{'apS'};
		#exit;
		$xlinkobj->set_lapS($searchhit->{'lapS'});	
		$xlinkobj->calcfullscore ($PARAMS->{'xcorr_tolerance_window'});
		
		if ($showScores){
		$xlinkobj->print_subscores;	
		}
		#$xlinkobj->printiontable;
		#$xlinkobj->calc_specific_ions;
		my $currenttime = time;
		my $pngfilename = join ".",basename($spectrumfilename), $currenttime, "png";
		$pngfilename =~ s/:/_/g;
		my $pngfilename_cp = join ".",basename($spectrumfilename), $currenttime, 'cp', "png";
		$pngfilename =~ s/:/_/g;
		my $specfilename = join ".",basename($spectrumfilename), $currenttime, "speccomparision.png";
		$specfilename =~ s/:/_/g;
		my $specfilename_cp = join ".",basename($spectrumfilename), $currenttime, 'cp', "speccomparision.png";
		$specfilename =~ s/:/_/g;
		my $pngname = File::Spec->catfile( $resultdirbase, $sessionid, "tmp", $pngfilename );
		$pngurl = join "\/", $resulturlbase, $sessionid, "tmp", $pngfilename;
		my $specname = File::Spec->catfile( $resultdirbase, $sessionid, "tmp", $specfilename );
		my $specurl = join "\/", $resulturlbase, $sessionid, "tmp", $specfilename;
		my ( $lossions, $logscale );

		if ( $form->param('logscale') )
		{
			$logscale = 1;
		}
		if ( $form->param('waterloss') || $form->param('nh3loss') )
		{
			$lossions = 1;
		}
		$xlinkobj->drawxlinkspec( $logscale, $lossions, $min, $max, $pngname, $labels, !$form->param('hide peptide structure'), $specname );
		print '<p><img width="800" height="600" src=', $pngurl, ' alt="Spectrum"></p>';
		if ( $xlinkobj->xlinktype eq "xlink" )
		{
			my $pngname_cp = File::Spec->catfile( $resultdirbase, $sessionid, "tmp", $pngfilename_cp );
			my $pngurl_cp = join "\/", $resulturlbase, $sessionid, "tmp", $pngfilename_cp;
			my $specname_cp = File::Spec->catfile( $resultdirbase, $sessionid, "tmp", $specfilename_cp );
			my $specurl = join "\/", $resulturlbase, $sessionid, "tmp", $specfilename_cp;
			$xlinkobj->drawxlinkspec_alpha_beta( $logscale, $lossions, $min, $max, $pngname_cp, $labels, !$form->param('hide peptide structure'), $specname_cp );
			print '<p><img width="800" height="600" src=', $pngurl_cp, ' alt="Spectrum"></p>';
		}
		my ( $lightspec, $heavyspec ) =  map { File::Spec->catfile( $resultdirbase, $sessionid, "tmp", $_ ) }
		  split /,/, $specobj->getspectrumheader;
		if ( $lightspec eq $heavyspec )
		{
			print '<p class="textblock">Only one spectrum was provided. Spectrum sorting cannot be applied. All peaks are matched against common- and xlink-ions</p>';
		}
		if ( $form->param('spectrumcomparison') )
		{
			
			if ( -e $lightspec && -e $heavyspec )
			{
				spectrum_comparison( $specobj, $specname, $lightspec, $heavyspec, $min, $max );
				print '<p><img width="800" height="600" src=', $specurl, ' alt="Spectrum not available"></p>';
			}
		}
	}
	if ( $form->param('original matches') )
	{
		$hithash = gethithash($matchions);
	} else
	{
		$hithash = $xlinkobj->getMatchingIonpositions;
	}
	$xlinkobj->printhtmliontable( $form->param('minionsize'), $form->param('maxionsize'), );
	
	print '<br>';
	
		print ' <table border="1">';
		print '<tr>';
		print "<td>Matched  common ions alpha:</td><td>", $xlinkobj->get_number_of_Commonmatches_alpha, "</td>";
		print "<td>Matched xlink ions alpha:</td><td>", $xlinkobj->get_number_of_Xlinkmatches_alpha, "</td>";
		print "<td>Matched ions alpha: </td><td>", $xlinkobj->get_num_of_matched_ions_alpha, "</td>";
		print '</tr><tr>';
		print "<td>Matched  common ions beta: </td><td>", $xlinkobj->get_number_of_Commonmatches_beta, "</td>";
		print "<td>Matched xlink ions beta: </td><td>", $xlinkobj->get_number_of_Xlinkmatches_beta, "</td>";
		print "<td>Matched ions beta: </td><td>", $xlinkobj->get_num_of_matched_ions_beta, "</td>";
		print '</tr>';
		print '</table>';
	
	if ( $form->param('hittable') )
	{
		$xlinkobj->printhittable($hithash);
	}
	
	
	#$xlinkobj->generate_subspec("beta");
	
	
	
}

sub printxlinkxcorr
{
	my $matchions = shift;
	my $PARAMS    = shift;
	my $xquestdef = shift;
	my $plusMr    = shift;
	my $hithash;
	my %buffer;
	my $specobj = getspecobj( $spectrumfilename, $xquestdefcopy,\%buffer );
	#my $xlinkobj = getxlinkions( $searchhit->{'type'}, $searchhit->{'seq1'}, $searchhit->{'seq2'}, $searchhit->{'xlmass'}, $searchhit->{'xlinkpos'}, $xquestdef, undef, $specobj, $plusMr );
	my $xlinkobj = getxlinkions( $searchhit->{'type'}, $searchhit->{'seq1'}, $searchhit->{'seq2'}, $searchhit->{'xlmass'}, $searchhit->{'xlpos'}, $xquestdefcopy, undef, $specobj, $plusMr );
		
	if ( $form->param('original matches') )
	{
		$hithash = gethithash($matchions);
	} else
	{
		$hithash = $xlinkobj->getMatchingIonpositions;
	}
	my $commonpeaks = $specobj->getcommonpairs;
	my $xlinkpeaks  = $specobj->getxlinkpairs;
	#$xlinkobj->match_offline($specobj);
	my $xcorrbase = join "_", $hitid, "xcorr", time;
	$xcorrbase =~ s/:/_/g;
	my $xcorrbasefilename = File::Spec->catfile( $resultdirbase, $sessionid, "tmp",$xcorrbase );
	my $xcorrburl = join "", $resulturlbase, "\/", $sessionid, "\/tmp\/", $xcorrbase, "_backbonecorr.png";
	my $xcorrxurl = join "", $resulturlbase, "\/", $sessionid, "\/tmp\/", $xcorrbase, "_xcorr.png";

	if ( $form->{'ms2tolerance'} )
	{
		$PARAMS->{'ms2tolerance'} = $form->param('ms2tolerance');
	}
	if ( $form->{'xlink_ms2tolerance'} )
	{
		$PARAMS->{'xlink_ms2tolerance'} = $form->param('xlink_ms2tolerance');
	}
	if ( $form->{'xcorrprecision'} )
	{
		$PARAMS->{'xcorrprecision'} = $form->param('xcorrprecision');
	}
	if ( $form->param('plotxcorrdelay') )
	{
		$PARAMS->{'plotxcorrdelay'} = $form->param('plotxcorrdelay');
	}else{
	#$PARAMS->{'plotxcorrdelay'} =10;	
	}
	my $xcorrprecision = int( 1 / $PARAMS->{'xcorrprecision'} );
	
	if ($PARAMS->{'xcorr_tolerance_window'})
	{
		print "<br> using normalized cross-correlation\n";
		$xlinkobj->xcorrelation_common_normalized( $xcorrbasefilename, 1, $PARAMS->{'plotxcorrdelay'} , $xcorrprecision, 0 );
		$xlinkobj->xcorrelation_xlink_normalized( $xcorrbasefilename, 1, $PARAMS->{'plotxcorrdelay'} , $xcorrprecision, 0 );
	} else
	{
		$xlinkobj->xcorrelation_common( $xcorrbasefilename, 1, $PARAMS->{'plotxcorrdelay'}, $xcorrprecision, 0 );
		$xlinkobj->xcorrelation_xlink( $xcorrbasefilename, 1, $PARAMS->{'plotxcorrdelay'} , $xcorrprecision, 0 );
	}

	print '<p><img src=', $xcorrburl, ' alt="xcorr_common"></p>';
	print '<p><img src=', $xcorrxurl, ' alt="xcorr_xlink"></p>';
}
### set the params that have to be remembered
sub setparams
{
	
	### STORE PARAMS
	### THESE ARE ALL THE NON-DYNAMIC PARAMS
	### THE ARE ALWAYS SUBMITTED AND STORED AS HIDDEN FIELDS
	print $form->hidden( -name => 'id',           -value => param('id') );
	print $form->hidden( -name => 'plottype',     -value => param('plottype') );
	print $form->hidden( -name => 'spectrum',     -value => param('spectrum') );
	print $form->hidden( -name => 'specfilename', -value => param('specfilename') );
	print $form->hidden( -name => 'type',         -value => param('type') );
	print $form->hidden( -name => 'xlid',         -value => param('xlid') );
	print $form->hidden( -name => 'seq1',         -value => param('seq1') );
	print $form->hidden( -name => 'seq2',         -value => param('seq2') );
	print $form->hidden( -name => 'xlpos',        -value => param('xlpos') );
	print $form->hidden( -name => 'xlmass',       -value => param('xlmass') );
	print $form->hidden( -name => 'scantype',       -value => param('scantype') );
	print $form->hidden( -name => 'lapS',       -value => param('lapS') );
	## SET THOSE PARAMS THAT ARE NOT PASSED BY THE URL AND THAT ARE DYNAMIC
	## ONLY STORE PARAMS AT INITIAL LOADING OF THE CGI
	unless ($form->param('reload')){
	print $form->hidden( -name => 'CID_match2ndisotope', -value => $PARAMS->{'CID_match2ndisotope'} );	
	print $form->hidden( -name => 'waterloss', -value => $PARAMS->{'waterloss'} );
	print $form->hidden( -name => 'nh3loss', -value => $PARAMS->{'nh3loss'} );
	}
	
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
## creates an insilico spectrum object from
## spectrumfile
sub getspecobj
{
	my $spectrumfile = shift;
	my $xquestdef    = shift;
	my $xquestdir;

	#print $spectrumfile;
	if ( defined( $ENV{'XQUEST_DIR'} ) )
	{
		$xquestdir = $ENV{'XQUEST_DIR'};
	} else
	{
		$xquestdir = ".";
	}
	my $masstable = File::Spec->catfile( $xquestdir, "mass_table.def" );
	#my $webconfig = File::Spec->catfile( $xquestdir, $webconfig );
	my $masstable = $env->get_path('mass.def');
	#open masstable and definition table
	my ( $MSTAB, $PARAMS, $ENZ, $basename, $WEBPARAMS ) = readtables( $masstable, $xquestdef, $webconfig );
	my %dummy;
	#print Dumper ($searchhit);
	my $specobj = Spectrum->new( $spectrumfile, $PARAMS,$searchhit->{'scantype'}, \%dummy );
	return $specobj;
}

sub gethithash
{
	my $matchions = shift;
	my $hithash;
	foreach my $ionseries ( $matchions->all_attr )
	{
		unless ( $ionseries = !/^_/ )
		{
			foreach my $position ( split /,/, $matchions->attr($ionseries) )
			{
				chomp;
				$hithash->{$ionseries}->{$position} = 1;
			}
		}
	}
	return $hithash;
}

sub getxlinkions
{
	my $xlinktype   = shift;
	my $pep1        = shift;
	my $pep2        = shift;
	my $xlinkermass = shift;
	my $topology    = shift;
	my $xquestdef   = shift;
	my $hithash     = shift;
	my $specobj     = shift;
	my $plusMr      = shift;
	my $xquestdir;
	my @topology = split /,/, $topology;

	if ( defined( $ENV{'XQUEST_DIR'} ) )
	{
		$xquestdir = $ENV{'XQUEST_DIR'};
	} else
	{
		$xquestdir = ".";
	}
	my $masstable = File::Spec->catfile( $xquestdir, "mass_table.def" );
	my $masstable = $env->get_path('mass.def');
	
	#my $webconfig = File::Spec->catfile( $xquestdir, $webconfig );

	#open masstable and definition table
	#my ( $MSTAB, $PARAMS, $ENZ, $basename, $WEBPARAMS ) = readtables( $masstable, $xquestdef, $webconfig );
	
	if ( $form->param('ionseries') )
	{
		my @ionseries = $form->param('ionseries');
		my %ionseries = ();
		foreach my $ion (@ionseries)
		{
			$ionseries{$ion} = 1;
		}
		$PARAMS->{'ionseries'} = \%ionseries;
	}
#	print "nh3loss:".$PARAMS->{'nh3loss'}."<br>";
	
	#### REVISE THE PARAMS IF RELOAD WAS USED ####
	#### OTHERWISE THE PARAMS FROM THE DEF ARE USED ####
	if ($form->param('reload')){
	if ( $form->param('CID_match2ndisotope') )
	{
	$PARAMS->{'CID_match2ndisotope'} = 1;
	}else{
	$PARAMS->{'CID_match2ndisotope'} = 0;
	}
	
	if ( $form->param('waterloss') )
	{
	$PARAMS->{'waterloss'} = 1;
	}else{
	$PARAMS->{'waterloss'} = 0;
	}
	
	if ( $form->param('nh3loss') )
	{
	$PARAMS->{'nh3loss'} = 1;
	}else{
	$PARAMS->{'nh3loss'} = 0;
	}	
	}
	
	my $pepobj1 = PepObj->new( $pep1, 'alpha', 'pep1', $MSTAB, $PARAMS, 0 );
	my $pepobj2 = PepObj->new( $pep2, 'beta',  'pep2', $MSTAB, $PARAMS, 0 );

	if ( $form->param('ms2tolerance') )
	{
		$PARAMS->{'ms2tolerance'} = $form->param('ms2tolerance');
	}
	if ( $form->param('xlink_ms2tolerance') )
	{
		$PARAMS->{'xlink_ms2tolerance'} = $form->param('xlink_ms2tolerance');
	}
	
my $xlinkobj;

if ($xlinktype eq "monolink" or $xlinktype eq "intralink"){
$xlinkobj = LinkObj->new( $xlinktype, [ $pepobj1 ], $xlinkermass, \@topology, undef, $PARAMS, $MSTAB, $specobj, 1, $plusMr );
}else{
$xlinkobj = LinkObj->new( $xlinktype, [ $pepobj1, $pepobj2 ], $xlinkermass, \@topology, undef, $PARAMS, $MSTAB, $specobj, 1, $plusMr );	
}


	return $xlinkobj;
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
	}
	close DEF;
	return \%PARAMS;
}



sub readtables
{
	my $table     = shift;
	my $xquestdef = shift;
	my $webconfig = shift;
	my $masslist  = shift;
	open TABLE, "<$table" or die "cannot open table $table $!";
	my ( %MSTAB, %PARAMS, %ENZ, %MOD, %WEBPARAMS );
	while ( my $line = <TABLE> )
	{
		chomp($line);
		if ($line)
		{
			my @results = split( " ", $line );
			$MSTAB{ $results[0] }->{'native'}  = $results[1];
			$MSTAB{ $results[0] }->{'average'} = $results[2];
			$MSTAB{ $results[0] }->{'-water'}  = $results[3];
			$MSTAB{ $results[0] }->{'-NH3'}    = $results[4];
		}
	}
	open DEF, "<$xquestdef" or die "cannot open xquest definition file $xquestdef $!";

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
			if ( $line ne "" )
			{

				#print $line."\n";
				my @results = split( ' ', $line );
				$PARAMS{ $results[0] } = $results[1];

				#print "key:".$results[0]."value ".$results[1]."\n";
			}

			#print $_[0], "\t", $PARAMS{ $_[0] }, "\n";
		} elsif ( defined($enzymedef) )
		{
			if ( $line ne "" )
			{
				my @results = split( ' ', $line );
				$ENZ{ $results[0] }->{'name'}     = $results[1];
				$ENZ{ $results[0] }->{'cutAA'}    = $results[3];
				$ENZ{ $results[0] }->{'notcutAA'} = $results[4];
			}
		} elsif ( defined($modifications) )
		{
			if ( $line ne "" )
			{
				my @results = split( ' ', $line );
				$MOD{ $results[0] } = $results[1];
			}
		}
	}

	#modify mass table for fixed modifications
	foreach my $key ( keys %MSTAB )
	{
		if ( $MOD{$key} )
		{

			#	$verbose && print "Fixed Modification defined: ", $key, " ", $MSTAB{$key}->{'native'}, "+ ", $MOD{$key}, "\n";
			$MSTAB{$key}->{'native'} += $MOD{$key};

			#print $MSTAB{$_}->{'native'},"\n";
		}
	}

	#define variable modification X
	if ( $PARAMS{'variable_mod'} )
	{
		my ( $AA, $delta ) = split /,|:/, $PARAMS{'variable_mod'};
		$MSTAB{'X'}->{'native'} = $MSTAB{$AA}->{'native'} + $delta;
	}

	#$verbose && print "modificaton X: ", $MSTAB{'X'}->{'native'}, "\n";
	if ( $PARAMS{'AArequired'} )
	{
		chomp( $PARAMS{'AArequired'} );
		my @xlinktargets = split /,|:|\|/, $PARAMS{'AArequired'};
		my $AAstring;
		foreach my $AA (@xlinktargets)
		{
			$AAstring = join "|", @xlinktargets;
		}
		$PARAMS{'AArequired'} = $AAstring;
	}
	if ( $PARAMS{'usenprescores'} )
	{
		my $use_nprescores = $PARAMS{'usenprescores'};
	}

	#$verbose && ( print "xlink targets: ", $PARAMS{'AArequired'}, "\n" );
	my $dbname;
	my $dbbasename=basename($PARAMS{'database'});

	if ( -e $PARAMS{'database'} )
	{
		$dbname = $PARAMS{'database'};
	} elsif ( -e File::Spec->catfile( $xquestdir, $PARAMS{'database'} ) )
	{
		$dbname = File::Spec->catfile( $xquestdir, $PARAMS{'database'} );
	} elsif( -e File::Spec->catfile( $resultdirbase, $sessionid, $dbbasename ) ){
		$dbname=File::Spec->catfile( $resultdirbase, $sessionid, $dbbasename );
	}else
	{
		die "cannot open database file $dbname $!";
	}
	## Get the Db path/basename is used for the db indices
	$dbname =~ s/\.\w+//;
	open WEBCONFIG, "<$webconfig"
	  or die "could not open web config file $webconfig $! ignoring";
	while (<WEBCONFIG>)
	{
		chomp;
		my @keyvalue = split /::/;    #:: alows for white space in file paths
		unless (/^#/)
		{
			$WEBPARAMS{ $keyvalue[0] } = $keyvalue[1];
		}
	}
	if ( $PARAMS{'xlinktypes'} )
	{                                 #compatibility stuff
		my @tmp = split //, $PARAMS{'xlinktypes'};
		$PARAMS{'search_monolinks'}       = $tmp[0];
		$PARAMS{'search_intralinks'}      = $tmp[1];
		$PARAMS{'search_intracrosslinks'} = $tmp[2];
		$PARAMS{'search_intercrosslinks'} = $tmp[3];
	}

	#MALDIspecific settings
	unless ( $PARAMS{'minpepmr'} )
	{
		$PARAMS{'minpepmr'} = 110 * $PARAMS{'mindigestlength'};
	}
	unless ( $PARAMS{'maxpepmr'} )
	{
		$PARAMS{'maxpepmr'} = 110 * $PARAMS{'maxdigestlength'};
	}
	unless ( $PARAMS->{'writetodiskaftern'} )
	{
		$PARAMS->{'writetodiskaftern'} = 100;
	}
	unless ( $PARAMS->{'tolerancemeasure_ms2'} )
	{
		$PARAMS->{'tolerancemeasure_ms2'} = "Da";    #$PARAMS->{'tolerancemeasure'};
	}
	if ( $PARAMS{'Iontagmode'} )
	{
		$PARAMS{'enumerate'} = 0;
		if ( $PARAMS{'Iontag_charges_for_index'} )
		{
			my @ioncharges = split /,/, $PARAMS{'Iontag_charges_for_index'};
			$PARAMS{'indexcharges_common'} = \@ioncharges;
		} else
		{
			$PARAMS{'indexcharges_common'} = [1];
		}
	}
	if ( $PARAMS{'ioncharge_common'} )
	{
		my @ioncharges = sort { $a <=> $b } split /,/, $PARAMS{'ioncharge_common'};
		$PARAMS{'ioncharge_common'} = \@ioncharges;
	} else
	{
		$PARAMS{'ioncharge_common'} = [ 1, 2 ];
	}
	if ( $PARAMS{'ioncharge_xlink'} )
	{
		my @ioncharges = sort { $a <=> $b } split /,/, $PARAMS{'ioncharge_xlink'};
		$PARAMS{'ioncharge_xlink'} = \@ioncharges;
	} else
	{
		$PARAMS{'ioncharge_xlink'} = [ 1, 2, 3, 4 ];
	}
	if ( $masslist && !$PARAMS{'massmatchonly'} )
	{
		warn "masslist $masslist indicated: changing searchmode to massmatchonly\b";
		$PARAMS{'massmatchonly'} = 1;
		$PARAMS{'enumerate'}     = 0;
		$PARAMS{'Iontagmode'}    = 0;
	}
	if ( $PARAMS{'ionseries'} )
	{
		my @ions = split //, $PARAMS{'ionseries'};
		my %ionseries = (
						  'a' => $ions[0],
						  'b' => $ions[1],
						  'c' => $ions[2],
						  'x' => $ions[3],
						  'y' => $ions[4],
						  'z' => $ions[5],
		);
		my @ionseries_array = ();
		foreach my $ion ( keys %ionseries )
		{
			if ( $ionseries{$ion} )
			{
				push @ionseries_array, $ion;
			}
		}
		$PARAMS{'ionseries'}       = \%ionseries;
		$PARAMS{'ionseries_array'} = \@ionseries_array;
	}
	my $ms2masstype = 'native';
	if ( $PARAMS{'averageMS2'} )
	{
		$ms2masstype = 'average';
	}
	if ( $PARAMS{'define_enzyme'} )
	{
		print $PARAMS{'define_enzyme'}, "\n";
		my ( $cutterm, $AAdefinitions ) = split /@/,  $PARAMS{'define_enzyme'};
		my ( $cutsite, $noncutsite )    = split /\^/, $AAdefinitions;
		my @cutsites    = split //, $cutsite;
		my @noncutsites = split //, $noncutsite;
		$PARAMS{'cutAA'}    = join "|", @cutsites;
		$PARAMS{'notcutAA'} = join "|", @noncutsites;
		$PARAMS{'cutterm'}  = $cutterm;
		print "AA to cut: ",     $PARAMS{'cutAA'},    "\n";
		print "AA not to cut: ", $PARAMS{'notcutAA'}, "\n";
		print "cut terminus: ",  $PARAMS{'cutterm'},  "\n";
	}
	my %fragmentresiduals = (
							  'a' => $PARAMS{'a_ion'},
							  'b' => $PARAMS{'b_ion'},
							  'c' => $PARAMS{'c_ion'},
							  'x' => $PARAMS{'x_ion'},
							  'y' => $PARAMS{'y_ion'},
							  'z' => $PARAMS{'z_ion'},
	);
	$PARAMS{'fragmentresiduals'} = \%fragmentresiduals;
	$PARAMS{'Hatom'}             = $MSTAB{'Hatom'}->{$ms2masstype};
	return ( \%MSTAB, \%PARAMS, \%ENZ, $dbname, \%WEBPARAMS );
}

sub spectrum_comparison
{
	my $specobj     = shift;
	my $pngfilename = shift;
	my $lightspec   = shift;
	my $heavyspec   = shift;
	my $min         = shift;
	my $max         = shift;

	#	my $base=shift;
	unless ( -e $lightspec && -e $heavyspec )
	{
		return 0;
	}
	my $commonpeaks = $specobj->getcommonpairs;

#		foreach my $peak (@$commonpeaks){
#			print "Common peak: @$peak<br>";
#		}
	
	my $xlinkpeaks = $specobj->getxlinkpairs;
	
#		foreach my $peak (@$xlinkpeaks){
#			print "Xlink peak: @$peak<br>";
#		}
	
	my @PL1        = ();
	my @PL2        = ();
	open PL1, "<$lightspec" or die $!;
	<PL1>;    #purge first line;
	while (<PL1>)
	{
		chomp;
		my @tmp = split;

		#	print "@tmp<br>";
		push @PL1, \@tmp;
	}
	close(PL1);
	open PL2, "<$heavyspec" or die $!;
	<PL2>;    #purge first line;
	while (<PL2>)
	{
		chomp;
		my @tmp = split;
		push @PL2, \@tmp;
	}
	close(PL2);
	foreach my $peakpair (@PL2)
	{
		$peakpair->[1] *= -1;
	}
	my $specobj1 = specplot->new();
	$specobj1->setcolor( "black", "black", "green", "red" );
	$specobj1->plotdata( $min, $max, [ "black", "black", "green", "red" ], \@PL1, \@PL2, $commonpeaks, $xlinkpeaks );

	#	$specobj1->drawlegend( 200, 10, { basename($pngfilename), "black" } );
	my %colorhash = (
					  "positive: light spectrum peaks" => "grey",
					  "negative: heavy spectrum peaks" => "grey",
					  "matched common-spectrum peaks"  => "green",
					  "matched xlink-spectrum peaks"   => "red",
	);
	my @sortlist = ( "positive: light spectrum peaks", "negative: heavy spectrum peaks", "matched common-spectrum peaks", "matched xlink-spectrum peaks" );
	$specobj1->drawlegend( 600, 20, \%colorhash, \@sortlist );
	$specobj1->printimage($pngfilename);
}

sub scaleintensity
{
	my $pairs  = shift;
	my $scale  = shift;
	my $factor = shift;
	unless ($factor)
	{
		$factor = 1;
	}
	$pairs->[1] = $factor * ( $pairs->[1] / $scale );
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
