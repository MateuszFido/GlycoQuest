#!/usr/bin/env perl
use strict;
use warnings;

#---------------------------------------------------------------------------
# compare_peaks3.pl
# Program to compare light and heavy MS/MS scan pairs of cross-linked peptides.
# Execute compare_peaks3.pl -help to display information and usage options.
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

# include the directory for the xquest modules ###########
use File::Spec::Functions qw(rel2abs);
use File::Basename;
use FindBin;
use lib "$FindBin::Bin/../modules";
#use lib "$FindBin::Bin/../../perl5";
##########################################################
use Getopt::Long;
use File::Basename;
use File::Spec;
use File::Copy;

use Statistics;
use XML::TreeBuilder;
use xQUEST_mzXMLscan;
use xQuest_mgfscanTTT;
use MIME::Base64;
use Data::Dumper;

# Default definitions
my $Hatom                = 1.007825032;
my $C13shift             = 1.0033548378;
my $threshold            = 1;
my $tolerance            = 0.2;
my $tolerancexl          = 0.3;
my $peakratio            = 0.1;
my $usenhighest          = 150;
my $average_common_peaks = 0;
my $outfile              = "test";
my $isotopediff          = 12.075321;
my $path                 = ".";
my $directory            = ".";
my $precision            = 10;
my $range                = 15;
my $minchargestate       = 1;
my $mincommonions        = 0;
my $scaleby              = 'max';
my $dynamic_range        = 0;
my $highresolutionms2    = 0;
my $minionsize           = 200;
my $maxionsize           = 1600;
my ( $xquestdef, $jslist, $id, $filelist, $mz, $charge, $superverbose, $input1, $input2, $verbose, $help, $i, $j, $matchlist, $xcorrelation, $dtafiles, $lightonlyfile, $addcommon2xlink, $weboutput, $scaleintensity, $printspectra, $resultdir, $genxml, $genmgf, $force, $filenamespec,
	 $tolerancemeasure, $deconvolutexlinkions );
GetOptions(
			'match=s'           => \$matchlist,              ## param is passed by runxquest
			'resultdir=s'       => \$resultdir,              ## param is passed by runxquest
			'mincharge=s'       => \$minchargestate,
			'filelist=s'        => \$filelist,
			'xrange=s'          => \$range,
			'av'                => \$average_common_peaks,
			'xprecision=s'      => \$precision,
			'mz=s'              => \$mz,
			'charge=s'          => \$charge,
			'xcorr=s'           => \$xcorrelation,
			'pr=s'              => \$peakratio,
			'path=s'            => \$path,
			'dtafiles=s'        => \$dtafiles,
			'lightonly=s'       => \$lightonlyfile,
			'jslist=s'          => \$jslist,
			'def=s'             => \$xquestdef,
			'thr=s'             => \$threshold,
			'tol=s'             => \$tolerance,
			'out=s'             => \$outfile,
			'verbose'           => \$verbose,
			'sv'                => \$superverbose,
			'web'               => \$weboutput,
			'dir=s'             => \$directory,              ## param is passed by runxquest as machedfilename+dir
			'help'              => \$help,
			'genxml=s'          => \$genxml,                 ## XML mode, write spec to an xml file, param is the filename of the original xml file
			'specfilenameout=s' => \$filenamespec,           ## if given as parameter this filename is used for the specxml file (paralell processing) otherwise the fn is created from the ml name
			'cpforce'           => \$force,
			'genmgf=s'          => \$genmgf,                 ## parameter for parsing an mgf file
);

#---------------------------------------------------------------------------
# Information
#---------------------------------------------------------------------------
my $version    = "3.4";
my $scriptinfo = {};
$scriptinfo->{'version'} = $version;
$scriptinfo->{'author'}  = "Thomas Walzthoeni modified from orginal script by Oliver Rinner";

## prints some information
my $report = 1;

# omit output buffering ## useful on cluster
$| = 1;

&usage unless ( $matchlist || $filelist || $dtafiles || $lightonlyfile );
&usage if $help;
if ($lightonlyfile)
{
	$addcommon2xlink = 1;
}
$report && print "Comparepeaks Parameters:\n -match $matchlist\n -dir $directory\n -resultdir $resultdir\n";
if ( $genxml && $report ) { print "-XMLmode: mzXML: $genxml\n" }
$range *= $precision;
my $defhash;

#---------------------------------------------------------------------------
#  Read the parameters from xquest.def if it is provided
#  THE standard defined VARS are OVERWRITTEN by the Params!
#---------------------------------------------------------------------------
if ($xquestdef)
{
	!$weboutput && ( print "Spectrum sorting: reading parameters from $xquestdef\n" );
	$defhash = readtables($xquestdef);
	$defhash->{'cp_threshold'} && ( $threshold = $defhash->{'cp_threshold'} );
	!$weboutput && ( print "cp_threshold ", $defhash->{'cp_threshold'}, "\n" );
	$defhash->{'cp_nhighest'} && ( $usenhighest = $defhash->{'cp_nhighest'} );
	!$weboutput && print "cp_nhighest ", $defhash->{'cp_nhighest'}, "\n";
	$defhash->{'cp_tolerancemeasure'} && ( $tolerancemeasure = $defhash->{'cp_tolerancemeasure'} );
	!$weboutput && print "cp_tolerancemeasure ", $defhash->{'cp_tolerancemeasure'}, "\n";
	$defhash->{'cp_tolerance'} && ( $tolerance = $defhash->{'cp_tolerance'} );
	!$weboutput && print "cp_tolerance ", $defhash->{'cp_tolerance'}, "\n";
	$defhash->{'cp_tolerancexl'} && ( $tolerancexl = $defhash->{'cp_tolerancexl'} );
	!$weboutput && print "cp_tolerancexl ", $defhash->{'cp_tolerancexl'}, "\n";
	$defhash->{'cp_peakratio'} && ( $peakratio = $defhash->{'cp_peakratio'} );
	!$weboutput && print "cp_peakratio ", $defhash->{'cp_peakratio'}, "\n";
	$defhash->{'cp_isotopediff'}  && ( $isotopediff  = $defhash->{'cp_isotopediff'} );
	$defhash->{'cp_spectralplot'} && ( $printspectra = $defhash->{'cp_spectralplot'} );
	!$weboutput && print "cp_isotopediff ", $defhash->{'cp_isotopediff'}, "\n";
	$defhash->{'cp_xcorrelation'} && ( $xcorrelation = $defhash->{'cp_xcorrelation'} );
	!$weboutput && print "cp_xcorrelation ", $defhash->{'cp_xcorrelation'}, "\n";
	$defhash->{'cp_dynamic_range'} && ( $dynamic_range = $defhash->{'cp_dynamic_range'} );
	!$weboutput && print "cp_dynamic_range ", $defhash->{'cp_dynamic_range'}, "\n";
	$defhash->{'cp_scaleby'} && ( $scaleby = $defhash->{'cp_scaleby'} );
	!$weboutput && print "cp_scaleby ", $defhash->{'cp_scaleby'}, "\n";
	$defhash->{'cp_scaleintensity'} && ( $scaleintensity = $defhash->{'cp_scaleintensity'} );
	!$weboutput && print "cp_scaleintensity ", $defhash->{'cp_scaleintensity'}, "\n";
	$defhash->{'cp_minpeaknumber'} && ( $mincommonions = $defhash->{'cp_minpeaknumber'} );
	!$weboutput && print "cp_minpeaknumber ", $defhash->{'cp_minpeaknumber'}, "\n";
	## ADDED PARAMETER FOR HIGH RESOLUTION MS2
	$defhash->{'cp_highresolutionms2'} && ( $highresolutionms2 = $defhash->{'cp_highresolutionms2'} );
	!$weboutput && print "cp_highresolutionms2 ", $defhash->{'cp_highresolutionms2'}, "\n";
	## ADDED PARAMETERs FOR THE SPECTRUM PLOT
	$defhash->{'minionsize'} && ( $minionsize = $defhash->{'minionsize'} );
	!$weboutput && print "minionsize ", $defhash->{'minionsize'}, "\n";
	$defhash->{'maxionsize'} && ( $maxionsize = $defhash->{'maxionsize'} );
	!$weboutput && print "maxionsize ", $defhash->{'maxionsize'}, "\n";
	## Added for deconvolution
	$defhash->{'cp_deconvolute_xlink_ions'} && ( $deconvolutexlinkions = $defhash->{'cp_deconvolute_xlink_ions'} );
	!$weboutput && print "deconvolute xlink ions: ", $defhash->{'cp_deconvolute_xlink_ions'}, "\n";
} else
{
	die "xquest definition file $xquestdef was not found $!";
}

#---------------------------------------------------------------------------
#  Create resultdirectories
#---------------------------------------------------------------------------
### if not defined, then default is "."
unless ( -e $directory )
{
	$report && print "cannot find spectrum directory $directory ... create new\n";
	mkdir $directory;
} else
{
	$report && print "Spectrum directory $directory exists.\n";
}
### create the xquest result directory if it is do not exists, otherwise the orginal dta will
## not be copied to the resultdirectory
unless ( -e $resultdir )
{
	$report && print "cannot find result directory $resultdir ... create new\n";
	mkdir $resultdir;
} else
{
	$report && print "Result directory $resultdir exists.\n";
}

#---------------------------------------------------------------------------
#  EOF Create resultdirectories
#---------------------------------------------------------------------------
my $lower_ratio = $peakratio;
my $upper_ratio = 1 / $peakratio;
my ( @MATCH, $matchpair );
if ($matchlist)
{
	open MATCH, "<$matchlist" or die "cannot open matchlist $matchlist $!";
	$matchlist =~ s/\.txt//;
	### check if listfile already exists
	my $listfilname = $matchlist . "_isotopepairs.txt";
	print $listfilname,"\n";
	if ($genxml)
	{
		if ( -e $listfilname && !$force )
		{
			if ($genxml)
			{
				unless ($filenamespec)
				{
					$filenamespec = "$resultdir" . ".spec." . "xml";
				}
				my $xmlspectfile = File::Spec->catfile( $resultdir, $filenamespec );
				if ( -e $xmlspectfile )
				{
					print "Filelist $listfilname and specxml $xmlspectfile file already exists, will skip compare_peaks3.pl, use -cpforce to reextract\n";
					exit;
				}
			} else
			{
				print "Filelist $listfilname already exists, will skip compare_peaks3.pl, use -cpforce to reextract\n";
				exit;
			}
		}
	}
	open FILELIST, ">$listfilname" or die "cannot open matchlist $listfilname $!";    ## create a FH for the isopairs

	#---------------------------------------------------------------------------
	#  PARSE THE MATCHLIST
	#---------------------------------------------------------------------------
	while (<MATCH>)
	{
		chomp;
		my @tmp = split;
		push @MATCH, \@tmp;
	}
	close(MATCH);
} elsif ( $filelist && $mz && $charge )
{
	chomp($filelist);
	my @files = map { basename($_) } split /,/, $filelist;
	push @MATCH, [ ( join ",", @files ), $mz, $charge, ( split /,/, $filelist ), 'light', 'heavy' ];
} elsif ($dtafiles)
{
	chomp($dtafiles);
	my ( $light, $heavy ) = split /,/, $dtafiles;
	open LIGHT, "<$light" or die $!;
	my $lightheader = <LIGHT>;
	chomp($lightheader);
	my ( $Mrplus1, $charge ) = split /\s+/, $lightheader;
	my $mz = ( $Mrplus1 - $Hatom + ( $charge * $Hatom ) ) / $charge;
	my @files = map { basename($_) } split /,/, $dtafiles;
	print "mz: $mz charge $charge\n";
	push @MATCH, [ ( join ",", @files ), $mz, $charge, ( split /,/, $dtafiles ), 'light', 'heavy' ];
} elsif ($lightonlyfile)
{
	chomp($lightonlyfile);
	open LIGHT, "<$lightonlyfile" or die $!;
	my $lightheader = <LIGHT>;
	chomp($lightheader);
	my ( $Mrplus1, $charge ) = split /\s+/, $lightheader;
	my $mz = ( $Mrplus1 - $Hatom + ( $charge * $Hatom ) ) / $charge;
	print "calculated from dta file m/z: $mz charge $charge\n";
	push @MATCH, [ ( join ",", basename($lightonlyfile), basename($lightonlyfile) ), $mz, $charge, $lightonlyfile, $lightonlyfile, 'light', 'light' ];
} else
{
	&usage();
}
my %seen        = ();
my $nonredpairs = 0;
my $npairs      = 0;
my $notfound    = 0;
my $rejected    = 0;
#### XML MODE: Parse XML if genxml parameter is selected
if ($genxml)
{
	print "Youre in xml mode using xml for saving spectra\n";
	#### define the outputfilename
	## outputpath is already the resultdirectory
	my $filenamespec = "$resultdir" . ".spec." . "xml";
	my $xmlspectfile = File::Spec->catfile( $resultdir, $filenamespec );
	## ADD a CHECK HERE IF THE FILE ALREADY EXISTS!!
	## IF YES USE A FORCE SEMAPHORE
	if ( -e $xmlspectfile && !$force )
	{
		print "Spec XML already exists, use -cpforce to rerun compare_peaks3\n";
		exit;
	}
	open INPUT, "<$genxml" or die "Cannot open xml file $!";
	my $MS2scans = parse_mzXML($genxml);
	close(INPUT);

	#print "XMLSPECOUTFILE is $xmlspectfile\n";
	open XMLSPECFILE, ">$xmlspectfile" or die "cannot open $xmlspectfile for writing! $!\n";
	openXMLspecHeader(*XMLSPECFILE);
### matchpair
	foreach $matchpair (@MATCH)
	{
		$npairs++;
		$nonredpairs += process_matched_scans_gen_xml( $matchpair, *FILELIST, *XMLSPECFILE, $MS2scans );
	}
	closeXMLspecHeader(*XMLSPECFILE);
	close(XMLSPECFILE);
} elsif ($genmgf)
{
	print "Mgf mode parsing an mgf file\n";
	#### define the outputfilename
	## outputpath is already the resultdirectory
	my $filenamespec = "$resultdir" . ".spec." . "xml";
	my $xmlspectfile = File::Spec->catfile( $resultdir, $filenamespec );
	## ADD a CHECK HERE IF THE FILE ALREADY EXISTS!
	## IF YES USE A FORCE SEMAPHORE
	if ( -e $xmlspectfile && !$force )
	{
		print "Spec XML already exists, use -cpforce to rerun compare_peaks3\n";
		exit;
	}
	my $MS2scans = parse_mgf($genmgf);

	#print "XMLSPECOUTFILE is $xmlspectfile\n";
	open XMLSPECFILE, ">$xmlspectfile" or die "cannot open $xmlspectfile for writing! $!\n";
	openXMLspecHeader(*XMLSPECFILE);
### Matchpair
	foreach $matchpair (@MATCH)
	{
		$npairs++;
		$nonredpairs += process_matched_scans_gen_xml_from_mgf( $matchpair, *FILELIST, *XMLSPECFILE, $MS2scans );
	}
	closeXMLspecHeader(*XMLSPECFILE);
	close(XMLSPECFILE);
} else
{
## old peakmatching
	foreach $matchpair (@MATCH)
	{
		$npairs++;
		$nonredpairs += process_matched_scans( $matchpair, *FILELIST );
	}
}
print "Processed $npairs spectral pairs, $nonredpairs are non-redundant\n";
print "$notfound Spectra were not found\n";
print "$rejected Spectra were rejected (not enough matched peaks)\n";
if ( $nonredpairs < 1 )
{
	print "<div class=\"error\">Number of non-redundant valid spectra is 0 - nothing to search </div>";
	exit 0;
}
close(FILELIST);
#### END OF CP
sub openXMLspecHeader
{
	my $specoutfile = shift;
	my $date        = localtime;
	print $specoutfile "<?xml version=\"1.0\" encoding=\"UTF-8\"?>";
	print $specoutfile "<xquest_spectra compare_peaks_version=\"$version\" date=\"$date\" author=\"Thomas Walzthoeni,Oliver Rinner\" homepage=\"http://proteomics.ethz.ch\" resultdir=\"", $resultdir, "\" deffile=\"", basename($xquestdef), "\" ";
	print $specoutfile ">\n";
}

sub closeXMLspecHeader
{
	my $specoutfile = shift;
	print $specoutfile "<\/xquest_spectra>";

}
### returns a hash with all ms2 scans as xQuest_mzXMLscam Objects
### the scannumber is is used as an index.
sub parse_mzXML
{
	my $mzXMLfile = shift;
	my $scanshash = {};
	## get the basename
	( my $mzXMLfilename = $mzXMLfile ) =~ s/.mzXML//;

	#print $mzXMLfilename;
	print "parsing $mzXMLfile...\n";
	my $mzXMLtree = XML::TreeBuilder->new();
	$mzXMLtree->parse_file($mzXMLfile);
	## Extract all scans and put them into a hash with the dta filename as index.
	foreach my $scan ( $mzXMLtree->find('scan') )
	{
		my $scannr = $scan->attr('num');

		#print $scannr."\n";
		unless ($scannr)
		{
			die "Cannot index scans: No scannumber found!\n";
		}
		if ( $scan->attr_get_i('msLevel') == 2 )
		{
			## check if there are peaks in this scan
			my $peaknum = $scan->attr_get_i('peaksCount');
			unless ( $peaknum == 0 )
			{
				$scanshash->{$scannr} = XQUEST_mzXMLscan->new( $scan, $scannr, $mzXMLfilename, 1 );
			} else
			{
				print "No peaks found in scan number $scannr, skip scan\n";
				next;
			}
		}
	}
	$mzXMLtree->delete;
	my $numberofscans = scalar( keys %$scanshash );
	print "read ", $numberofscans, " MS/MS scans\n";
	return $scanshash;
}

#---------------------------------------------------------------------------
#  Parse mgf file and generate a hash with all scanobjects and scannumber as index
#---------------------------------------------------------------------------
sub parse_mgf
{
	my $mgffile   = shift;
	my $mgf       = read_mgf_file( $mgffile, 0 );
	my $scanshash = {};
	my @keys      = keys %$mgf;
	my @MS2scans  = ();
	foreach my $key (@keys)
	{
		my $mgfarray = $mgf->{$key};
		my $scan     = xQuest_mgfscanTTT->new($mgfarray);
		my $title    = $scan->title;

		#print "$title\n";
		$scanshash->{$title} = $scan;
	}
	my $numberofscans = scalar( keys %$scanshash );
	print "read ", $numberofscans, " MS/MS scans from mgf file\n";
	return $scanshash;
}

#===  FUNCTION  ================================================================
#  NAME:        read_mgf_file
#  PURPOSE:     read_mgf_file
#  DESCRIPTION: Read mgf file by scans and put results into an hash with arrays
#  PARAMETERS:  $filename (including path if script is not executed in the folder)
#  RETURNS:     Hash
#===============================================================================
sub read_mgf_file
{
	my $filename = shift;
	my $verbose  = shift;
	my @array;
	if ($verbose) { print "Reading from file $filename\n" }
	open FILE, $filename or die $!;
	my $start = 0;
	my $num   = 0;
	my $scans = {};

	while ( my $line = <FILE> )
	{
		if ($verbose) { print "Reading line $line" }
		if ( $line =~ m/BEGIN IONS/ )
		{
			$num++;
		}
		chomp($line);
		unless ($line)
		{

			#next;
		}
		push @{ $scans->{$num} }, $line;
	}
	($verbose) && print "\n";
	close FILE;
	return $scans;
}

#---------------------------------------------------------------------------
#  FUNCTION FOR PROCESSING SCANS from MGF
#---------------------------------------------------------------------------
sub process_matched_scans_gen_xml_from_mgf
{
	my $matchpair      = shift;
	my $filelisthandle = shift;                  ### Filelist is the isotopelist outfilename
	my $XMloutfile     = shift;                  ### xml specfile
	my $MS2Scans       = shift;                  ### xQuest_mzXMLscanobjects
	my $id             = $matchpair->[0];        #M09-10341.c.03492.03492.3.dta,M09-10341.c.03368.03368.3.dta
	my $mz             = $matchpair->[1];
	my $charge         = $matchpair->[2];
	my $input1         = $matchpair->[3];        #./M09-10341.c/M09-10341.c.03492.03492.3.dta
	my $input2         = $matchpair->[4];        #./M09-10341.c/M09-10341.c.03368.03368.3.dta
	my $scantype1      = $matchpair->[5];        # light
	my $scantype2      = $matchpair->[6];        # heavy
	my $scans          = $matchpair->[7];        # 3492:3368
	my $rttimes        = $matchpair->[8];        # 1234:1238 ## rt in seconds
	my $mzofscans		= $matchpair->[9];		 # 531.53973389:535.55877686

	my @scanssplit     = split( /:/, $scans );
	my $scan1          = $scanssplit[0];
	my $scan2          = $scanssplit[1];

	#print Dumper ($MS2Scans);
	unless ( $scan1 or $scan2 )
	{
		die "Found no scannumber for id $id in matchlist\n";
	}
	$report && print "\nMatching Scans: Scan1: $scan1 and Scan2:$scan2\n";
	## check if the scannumbers exists
	##
	unless ( ( $MS2Scans->{$scan1} ) && ( $MS2Scans->{$scan2} ) )
	{
		warn "Scanobject for id $id does not exists!\n";
		$notfound++;
		next;
	}
	## get the scan objects
	my $scan1obj                = $MS2Scans->{$scan1};
	my $scan2obj                = $MS2Scans->{$scan2};
	my $isotopicshift4xlinkions = 0;
	$addcommon2xlink = ();
	if ( ( $scantype1 eq 'heavy' ) && ( $scantype2 eq 'heavy' ) )
	{
		$isotopicshift4xlinkions = $isotopediff;
	}
	### IS USED IF LIGHTONLY IS SELECTED
	if ( $scantype1 eq $scantype2 )
	{
		$addcommon2xlink = 1;
	}
	$verbose && print "id: $id mz: $mz charge: $charge dta-file1: $input1 dta-file2: $input2 scan1: $scantype1 scan2: $scantype2\n";
	if ( ( $input1 eq "void" ) || ( $input2 eq "void" ) )
	{
		warn "There seems to be an empty line in the matchlist, will skip it\n";
		next;
	}
	$xcorrelation && xcorr::xcorrelation( $input1, $input2, $range, undef, $precision, $directory );
	my ( @PL1, @PL2, @intensities1, @intensities2 );
#### READ THE PEAKLISTS FROM THE XML HASH
	#print Dumper ($scan1obj);
	## SPECTRUM 1
	my $scan1peaks = $scan1obj->get_peaklist;    ## is a array reference / only the peaaks
	                                             # print Dumper($scan1peaks);
	## get the basename of input 1
	$input1 =~ s/\.dta//;
	### generate the header (same as for a dta)
	my $mr1 = $scan1obj->Mr;
	my $z1  = $scan1obj->FT_charge;
	unless ( $mr1 or $z1 )
	{
		die "Charge or Mr not set!\n";
	}
	my $datheader1 = "$mr1\t$z1\n";
	my $specline   = 0;
	foreach my $peakpair (@$scan1peaks)
	{    ## peakpair is an array ref
		my $intensity = $peakpair->[1];
		my $mz        = $peakpair->[0];
		$specline++;
		unless ( defined( $intensity && $mz ) )
		{
			die "Error Specline $specline is not valid\n";
		}
		## check if the intensity is above the threshold
		if ( $intensity >= $threshold )
		{
			push @PL1,          $peakpair;
			push @intensities1, $peakpair->[1];
		}
	}
	my $numpeaksPL1 = scalar(@PL1);
	$report && print "#Ions Scan1 (after sorting intens. > $threshold): $numpeaksPL1\n";
	## SPECTRUM 2
	my $scan2peaks = $scan2obj->get_peaklist;    ## is a array reference / only the peaaks
	                                             #print Dumper($scan1peaks);
	## get the basename of input 1
	$input2 =~ s/\.dta//;
	### generate the header (same as for a dta)
	my $mr2 = $scan2obj->Mr;
	my $z2  = $scan2obj->FT_charge;
	unless ( $mr2 or $z2 )
	{
		die "Charge or Mr not set!\n";
	}
	my $datheader2 = "$mr2\t$z2\n";
	$specline = 0;
	foreach my $peakpair (@$scan2peaks)
	{    ## peakpair is an array ref
		my $intensity = $peakpair->[1];
		my $mz        = $peakpair->[0];
		$specline++;
		unless ( defined( $intensity && $mz ) )
		{
			die "Error Specline $specline is not valid\n";
		}
		## check if the intensity is above the threshold
		if ( $intensity >= $threshold )
		{
			push @PL2,          $peakpair;
			push @intensities2, $peakpair->[1];
		}
	}
	my $numpeaksPL2 = scalar(@PL2);

	#$report && print "Number of ions in spectrum2: $specline, after sorting (above intensity threshold:$threshold): $numpeaksPL2\n";
	$report && print "#Ions Scan2 (after sorting intens. > $threshold): $numpeaksPL2\n";

	# EOF LOADING THE PEAKS OF SPECTRUM2
	## set the path and basename for the result files
	$outfile = File::Spec->catfile( $directory, ( join "_", basename($input1), basename($input2) ) );
	if ( $seen{$outfile}++ )
	{
		return 0;
	}
	my $max1 = Statistics::max( \@intensities1 );
	my $max2 = Statistics::max( \@intensities2 );

	# Filtering of the peaklists by dynamic range
	if ($dynamic_range)
	{
		my $threshold1 = $max1 / $dynamic_range;
		my $threshold2 = $max2 / $dynamic_range;
		my @tmp1       = ();
		foreach my $peak (@PL1)
		{
			if ( intensity($peak) > $threshold1 )
			{
				push @tmp1, $peak;
			}
		}
		my @tmp2 = ();
		foreach my $peak (@PL2)
		{
			if ( intensity($peak) > $threshold2 )
			{
				push @tmp2, $peak;
			}
		}
		@PL1         = @tmp1;
		@PL2         = @tmp2;
		$numpeaksPL1 = scalar(@PL1);
		$numpeaksPL2 = scalar(@PL2);

		#$report && print "Number of ions in Spectrum 1 after filtereng by dynamic range (): $numpeaksPL1\n";
		#$report && print "Number of ions in Spectrum 2 after filtereng by dynamic range (): $numpeaksPL2\n";
	}
	## Rescaleing the Intensities ##############################################################################
	if ($scaleintensity)
	{
		my $topten1 = 0;
		my $topten2 = 0;
		if ( $scaleby eq 'top5' )
		{
			$topten1 = Statistics::topn( \@intensities1, 5 );
			!$weboutput
			  && print "mean intensity of  5\% most intense peaks for $input1: $topten1<br>\n";
			$topten2 = Statistics::topn( \@intensities2, 5 );
			!$weboutput && print "mean intensity of  5\% most intense peaks for $input2: $topten2<br>\n";
		} else
		{
			$topten1 = $max1;
			$topten2 = $max2;
			!$weboutput && print "<br>maximum intensity spectrum1: $topten1<br>\n";
			!$weboutput && print "maximum intensity spectrum2: $topten2<br>\n";
		}
		foreach my $peak (@PL1)
		{
			scaleintensity( $peak, $topten1, 100 );
			$superverbose && print mz($peak), "\t", intensity($peak), "\n";
		}
		$superverbose && print "\nspectrum2\n";
		foreach my $peak (@PL2)
		{
			scaleintensity( $peak, $topten2, 100 );
			$superverbose && print mz($peak), "\t", intensity($peak), "\n";
		}
	}
	## EOF Rescaleing the Intensities ##############################################################################
	## GENERATE ARRRAYS FOR THE PLOTS -> orginal peaks#####
	my @PL1_plot = @PL1;
	my @PL2_plot = @PL2;
	####
	# Modify the peaklist to include mispicked transformed peaks for peakmatching
	# by the deconvolution/deistoping algorithm, includes -C13 and -2*C13 peak
	if ( ($highresolutionms2) && !($addcommon2xlink) )
	{
		my @templist    = ();
		my @orginallist = ();
		print "Add -1C13 shift to every peak of the light peaklist\n";
		foreach my $peak (@PL1)
		{
			my $oldmz;
			my $newmz1;
			my $newmz2;
			my $newmz3;
			my $newmz4;
			$oldmz = mz($peak);
			my $intensity = intensity($peak);

			#print $C13shift;
			$newmz1 = $oldmz - $C13shift;
			$newmz2 = $oldmz - 2 * $C13shift;
			$newmz3 = $oldmz + $C13shift;
			$newmz4 = $oldmz + 2 * $C13shift;

			#print $oldmz." ".$newmz1."\n";
			my @orginalpeak = ( $oldmz,  $intensity );
			my @newpeak1    = ( $newmz1, $intensity );
			my @newpeak2    = ( $newmz2, $intensity );
			my @newpeak3    = ( $newmz3, $intensity );
			my @newpeak4    = ( $newmz4, $intensity );

			#print "oldmz: $oldmz orginalmz:$oldmz\n";
			push( @templist, \@orginalpeak, \@newpeak1, \@newpeak2, \@newpeak3, \@newpeak4 );
		}
		@PL1 = @templist;
	}
	## Saving the normalized files into the  XML as dta base 64 encoded#####
	my @peaklist_light = sort { mz($a) <=> mz($b) } @PL1;
	my @peaklist_heavy = sort { mz($a) <=> mz($b) } @PL2;
	my $normalized_dta1 = $input1 . ".dta";
	my $normalized_dta2 = $input2 . ".dta";

	#open PL1norm, ">$normalized_dta1" or die $!;
	#open PL2norm, ">$normalized_dta2" or die $!;
	#print PL1norm $datheader1, "\n";
	#print PL2norm $datheader2, "\n";
	my $plstring1 = $datheader1;
	foreach my $peak (@peaklist_light)
	{
		$plstring1 .= join "\t", @$peak, "\n";

		#print PL1norm "\n";
	}

	#close(PL1norm);
	my $plstring2 = $datheader2;
	foreach my $peak (@peaklist_heavy)
	{
		$plstring2 .= join "\t", @$peak, "\n";

		#print PL2norm join "\t", @$peak;
		#print PL2norm "\n";
	}

	#close(PL2norm);
	### Save the normalized spectra to the xml // done after peakmatching
	#save_spectra_to_XML_base64enc($normalized_dta1,$scantype1,$XMloutfile,$plstring1);
	#save_spectra_to_XML_base64enc($normalized_dta2,$scantype2,$XMloutfile,$plstring2);
	# For Debuging
	#$verbose =1;
	if ($verbose)
	{

		# print the peaklists
		print "\n#### PEAKLIST LIGHT: ####\n";
		foreach my $peak (@peaklist_light)
		{
			printpeak($peak);
		}
		print "#### END OF PEAKLIST LIGHT ###\n";
		print "\n#### PEAKLIST HEAVY: ####\n";
		foreach my $peak (@peaklist_heavy)
		{
			printpeak($peak);
		}
		print "#### END OF PEAKLIST HEAVY ###\n";
	}
	## PEAK MATCHING AND SEPARATION OF COMMON and Cross-linking peaks
	## search common peaks list1
	my ( @meancommon_mz, @meancommon_int, @meancommon_delta, @xlink, $pair, @commonpeaks_light, @commonpeaks_heavy );
	if ($verbose) { print "\n#### MATCHING of common ions of light to heavy spectrum: ####\n" }
	@commonpeaks_light = matchpeaks( \@peaklist_light, \@peaklist_heavy, $tolerance, $lower_ratio, $upper_ratio, 0 );
	if ( $average_common_peaks || $verbose )
	{
		if ($verbose) { print "\n#### MATCHING of common ions of heavy to light spectrum: ####\n" }
		@commonpeaks_heavy = matchpeaks( \@peaklist_heavy, \@peaklist_light, $tolerance, $lower_ratio, $upper_ratio, 0 );
	}
	## check if enough common peaks were found
	unless ( scalar(@commonpeaks_light) >= $mincommonions )
	{
		if ($weboutput)
		{
			print "<div class=\"error\"><br>",
			  "number of common peaks is smaller than the required minimum of $mincommonions: ",
			  scalar(@commonpeaks_light), " in $input1 and $input2\: #peaks1 > ", $threshold, ": ",
			  scalar(@peaklist_light), " | #peaks2 > ", $threshold, ": ",
			  scalar(@peaklist_heavy),, "->ignored</div>";
			$rejected++;
		} else
		{
			warn "number of common peaks is smaller than the required minimum of $mincommonions: ", scalar(@commonpeaks_light), " in $input1 and $input2\: #peaks1 > ", $threshold, ": ", scalar(@peaklist_light), " | #peaks2 > ", $threshold, ": ", scalar(@peaklist_heavy),, "->ignored\n";
			$rejected++;
		}
		return 0;
	}
#### SAVE THE NORMALIZED SPECTRA TO THE XML FILE
	### Save the normalized spectra to the xml
	#$report && print "Saving normalized spectrum $normalized_dta1\n";
	save_spectra_to_XML_base64enc( $normalized_dta1, $scantype1, $XMloutfile, $plstring1 );

	#$report && print "Saving normalized spectrum $normalized_dta2\n";
	save_spectra_to_XML_base64enc( $normalized_dta2, $scantype2, $XMloutfile, $plstring2 );
#### HERE THE LINE IS PRINTED TO THE ISOTOPEPAIR LIST
## add the rttimes and the mz values
	print $filelisthandle "$outfile\t$isotopicshift4xlinkions\t", $scantype1, "_", "$scantype2\t$rttimes\t$mzofscans\n";
	my $spline = "$outfile\t$isotopicshift4xlinkions\t" . $scantype1 . "_" . "$scantype2\t$rttimes\t$mzofscans\n";

	#print "Writing line to file: $filelisthandle: $spline";
#### PEAK MATCHING AND SEPARATION OF CROSSLINKER PEAKS $spline;
	#copy the orginal chargestate
	my $charge_orginal = $charge;
	if ($highresolutionms2)
	{

		# if highresolution ms2 are used set the charge state to 1
		# all peaks should be deisotoped and should have charge 1
		$charge = 1;
		$report && print "HIGH RESOLUTION MS2 charge is set to 1\n";
	}

	#  search xlink peaks list
	my ( @xlinkpeaks_light, @xlinkpeaks_heavy );
	### replace the ions of the light precursor to the xlinkpeaks if light or heavyonly is compared, no peakmatching of x-linker ions
	### all charge states are then also 0, which makes sense
	### the common and xlinker ion lists are then equal
	if ($addcommon2xlink)
	{
		@xlinkpeaks_light = @commonpeaks_light;
		@xlinkpeaks_heavy = @commonpeaks_heavy;

		# SET THE CHARGE STATE TO 1 in high res mode, peaks are already deconvoluted
		if ($highresolutionms2)
		{
			$report && print "High res ms2: charge of x-linker ions is set to 1";
			foreach my $peak (@xlinkpeaks_light)
			{
				setcharge( $peak, 1 );
			}
			foreach my $peak (@xlinkpeaks_heavy)
			{
				setcharge( $peak, 1 );
			}
		}
	} else
	{
		### Search for x-linker ions for all possible chage states, minchargestate=1, max = precursor chargestate
		#TODO: Should max charge not z-1 be?
		for my $chargestate ( $minchargestate .. $charge )
		{
			$verbose && print "searching for +$chargestate charged x-link-ions <br>\n";
			my @transformedpeaks = ();
			## makes a transformed peaklist where the heavy peaks are shifted to the left by the isodiff/chargestates
			## if these are compared with the light list they should match if a xlink is there
			foreach my $heavypeak (@peaklist_heavy)
			{
				push @transformedpeaks, [ $heavypeak->[0] - ( $isotopediff / $chargestate ), $heavypeak->[1] ];
			}
			push @xlinkpeaks_light, matchpeaks( \@peaklist_light, \@transformedpeaks, $tolerancexl, $lower_ratio, $upper_ratio, $chargestate );
			$verbose && ( push @xlinkpeaks_heavy, matchpeaks( \@transformedpeaks, \@peaklist_light, $tolerancexl, $lower_ratio, $upper_ratio, $chargestate ) );
		}
	}

	# END OF PEAKMATCHING
	if ($verbose)
	{

		# print the matched peaks
		print "\n#### Common PEAKS MATCHED LIGHT: ####\n";
		foreach my $peak (@commonpeaks_light)
		{
			printpeak($peak);
		}
		print "#### END OF MATCHED COMMON PEAKLIST LIGHT ####\n";
		print "\n#### Common PEAKS MATCHED HEAVY:####\n";
		foreach my $peak (@commonpeaks_heavy)
		{
			printpeak($peak);
		}
		print "#### END OF MATCHED COMMON PEAKLIST HEAVY ###\n";
	}
	my $commonmin = 0;
	my $xlinkmin  = 0;

	# sort out if nhighest param is used
	if ( scalar(@commonpeaks_light) > $usenhighest )
	{
		my @tmp = sort { intensity($b) <=> intensity($a) } @commonpeaks_light;
		$commonmin = intensity( $tmp[$usenhighest] );
	}
	if ( scalar(@xlinkpeaks_light) > $usenhighest )
	{
		my @tmp = sort { intensity($b) <=> intensity($a) } @xlinkpeaks_light;
		$xlinkmin = intensity( $tmp[$usenhighest] );
	}

	# AVERAGE the common peaks from the light and heavy spectrum
	if ($average_common_peaks)
	{
		for $i ( 0 .. $#commonpeaks_light )
		{
			$meancommon_mz[$i] = ( mz( $commonpeaks_light[$i] ) + mz( $commonpeaks_heavy[$i] ) ) / 2;
			$verbose && ( $meancommon_delta[$i] = sprintf( "%.4f", ( mz( $commonpeaks_light[$i] ) - mz( $commonpeaks_heavy[$i] ) ) ) );
			$meancommon_int[$i] = ( intensity( $commonpeaks_light[$i] ) + intensity( $commonpeaks_heavy[$i] ) ) / 2;
		}
	}
	$verbose && print "Common peaks tolerance=$tolerance, required peakratio = $peakratio\n";
	$verbose && print "mz1\tintensity1\tmz2\tintensity2\tmeanmz\tmeanintensity\tdelta\n";
	my $ncommonpeaks = 0;
	## GENERATE THE FILES FOR THE COMMON AND XLINKER Files
	my $commonheader = "$id\n$mz\n$charge_orginal\n";
	my $xlheader     = "$id\n$mz\n$charge_orginal\n";

	#print COMMONPEAKS "$id\n$mz\n$charge_orginal\n";
	#print XLINKERPEAKS "$id\n$mz\n$charge_orginal\n";
	my $commoncontent = $commonheader;
	my $xlcontent     = $xlheader;
	for ( $i = 0 ; $i <= $#commonpeaks_light ; $i += 1 )
	{
		if ( intensity( $commonpeaks_light[$i] ) >= $commonmin )
		{
			$ncommonpeaks++;
			$verbose && print mz( $commonpeaks_light[$i] ), "\t", intensity( $commonpeaks_light[$i] ), "\t\t", mz( $commonpeaks_heavy[$i] ), "\t", intensity( $commonpeaks_heavy[$i] ), "\t", $meancommon_mz[$i], "\t", $meancommon_int[$i], "\t", $meancommon_delta[$i], "\n";
			if ($average_common_peaks)
			{
				$commoncontent .= $meancommon_mz[$i] . "\t" . $meancommon_int[$i] . "\t0" . "\n";
			} else
			{
				$commoncontent .= mz( $commonpeaks_light[$i] ) . "\t" . intensity( $commonpeaks_light[$i] ) . "\t0" . "\n";
			}
		}
	}
	if ($verbose)
	{
		print "xlinker peaks\n";
		for ( $i = 0 ; $i <= $#xlinkpeaks_light ; $i += 1 )
		{
			print mz( $xlinkpeaks_light[$i] ), "\t", intensity( $xlinkpeaks_light[$i] ), "\tcharge: ", charge( $xlinkpeaks_light[$i] ), "\t";
			print mz( $xlinkpeaks_heavy[$i] ), "\t", intensity( $xlinkpeaks_heavy[$i] ), "\tcharge: ", charge( $xlinkpeaks_heavy[$i] ), "\t";
			print "\n";
		}
	}
	my $nxlinkpeaks = 0;
	foreach my $xlinkpeak (@xlinkpeaks_light)
	{
		if ( intensity($xlinkpeak) >= $xlinkmin )
		{
			$nxlinkpeaks++;
			### Deconvolute if the option is selected
			if ($deconvolutexlinkions)
			{
				my $deconvolutedpeak = deconvolute($xlinkpeak);
				$xlcontent .= mz($deconvolutedpeak) . "\t" . intensity($deconvolutedpeak) . "\t" . charge($deconvolutedpeak) . "\n";
			} else
			{
				$xlcontent .= mz($xlinkpeak) . "\t" . intensity($xlinkpeak) . "\t" . charge($xlinkpeak) . "\n";
			}
		}
	}
	#### SAVE THE COMMON AND XL SPECTRA TO THE XML FILE
	# Open the filehandles for common and xlinker peklists
	my $fncommon = $outfile . "_common.txt";
	my $fnxl     = $outfile . "_xlinker.txt";
	save_spectra_to_XML_base64enc( $fncommon, "common",  $XMloutfile, $commoncontent );
	save_spectra_to_XML_base64enc( $fnxl,     "xlinker", $XMloutfile, $xlcontent );
	### PRINT OVERLAY SPECTRUM ################################################
	if ($printspectra)
	{
		require specplot;
		my $specfile = join "", $outfile, "_specplot.png";
		foreach my $peakpair (@PL2_plot)
		{
			$peakpair->[1] *= -1;
		}
		foreach my $peakpair (@commonpeaks_heavy)
		{
			$peakpair->[1] *= -1;
		}
		foreach my $peakpair (@xlinkpeaks_heavy)
		{
			$peakpair->[1] *= -1;
		}
		my $specobj1 = specplot->new();

		#$specobj1->setcolor( "grey", "grey", "green", "red" );
		$specobj1->plotdata( $minionsize, $maxionsize, [ "grey", "grey", "green", "red" ], \@PL1_plot, \@PL2_plot, \@commonpeaks_light, \@xlinkpeaks_light );
		$specobj1->drawlegend( 200, 10, { basename($outfile), "black" } );
		my %colorhash = (
						  "positive: light spectrum peaks" => "grey",
						  "negative: heavy spectrum peaks" => "grey",
						  "matched common-spectrum peaks"  => "green",
						  "matched xlink-spectrum peaks"   => "red",
		);
		$specobj1->drawlegend( 600, 20, \%colorhash );
		$specobj1->printimage($specfile);
	}
	if ($weboutput)
	{
		print '<br><b>', basename($input1), ",", basename($input2), ': sorting peaks into common peaks and x-linker peaks:</b><br>';
		print 'number of peaks spectrum1  > ', $threshold, ': ',
		  scalar(@peaklist_light), '<br>number of peaks spectrum2 > ',
		  $threshold, ': ', scalar(@peaklist_heavy),
		  '<br>number of  common-peaks > ', sprintf( "%.2f", $commonmin ), ' ',
		  $ncommonpeaks, '<br>number of xlink-peaks > ',
		  sprintf( "%.2f", $xlinkmin ), ' ', $nxlinkpeaks, '<br>';
	} else
	{
		print "#peaks1 > ", $threshold, ": ", scalar(@peaklist_light), "|  #peaks2 > ", $threshold, ": ", scalar(@peaklist_heavy), "| #commonpeaks > ", sprintf( "%.2f", $commonmin ), "% ", $ncommonpeaks, "| #xlinkpeaks > ", sprintf( "%.2f", $xlinkmin ), "% ", $nxlinkpeaks, "<br>\n";
	}

	#exit;
	return 1;
}

#---------------------------------------------------------------------------
#  FUNCTION FOR PROCESSING WITH PSEUDOSH WORKFLOW
#---------------------------------------------------------------------------
sub process_matched_scans_gen_xml
{
	my $matchpair      = shift;
	my $filelisthandle = shift;                  ### Filelist is the isotopelist outfilename
	my $XMloutfile     = shift;                  ### xml specfile
	my $MS2Scans       = shift;                  ### xQuest_mzXMLscanobjects
	my $id             = $matchpair->[0];        #M09-10341.c.03492.03492.3.dta,M09-10341.c.03368.03368.3.dta
	my $mz             = $matchpair->[1];
	my $charge         = $matchpair->[2];
	my $input1         = $matchpair->[3];        #./M09-10341.c/M09-10341.c.03492.03492.3.dta
	my $input2         = $matchpair->[4];        #./M09-10341.c/M09-10341.c.03368.03368.3.dta
	my $scantype1      = $matchpair->[5];        # light
	my $scantype2      = $matchpair->[6];        # heavy
	my $scans          = $matchpair->[7];        # 3492:3368
	my $rttimes        = $matchpair->[8];        # 1234:1238 ## rt in seconds
	my $mzofscans		= $matchpair->[9];		 # 531.53973389:535.55877686
	my @scanssplit     = split( /:/, $scans );
	my $scan1          = $scanssplit[0];
	my $scan2          = $scanssplit[1];

	
	unless ( $scan1 or $scan2 )
	{
		die "Found no scannumber for id $id in matchlist\n";
	}
	$report && print "\nMatching Scans: Scan1: $scan1 and Scan2:$scan2\n";
	## check if the scannumbers exists
	##
	unless ( ( $MS2Scans->{$scan1} ) && ( $MS2Scans->{$scan2} ) )
	{
		warn "Scanobject for id $id does not exists!\n";
		$notfound++;
		next;
	}
	## get the scan objects
	my $scan1obj                = $MS2Scans->{$scan1};
	my $scan2obj                = $MS2Scans->{$scan2};
	my $isotopicshift4xlinkions = 0;
	$addcommon2xlink = ();
	if ( ( $scantype1 eq 'heavy' ) && ( $scantype2 eq 'heavy' ) )
	{
		$isotopicshift4xlinkions = $isotopediff;
	}
	### IS USED IF LIGHTONLY IS SELECTED
	if ( $scantype1 eq $scantype2 )
	{
		$addcommon2xlink = 1;
	}
	$verbose && print "id: $id mz: $mz charge: $charge dta-file1: $input1 dta-file2: $input2 scan1: $scantype1 scan2: $scantype2\n";
	if ( ( $input1 eq "void" ) || ( $input2 eq "void" ) )
	{
		warn "There seems to be an empty line in the matchlist, will skip it\n";
		next;
	}
	$xcorrelation && xcorr::xcorrelation( $input1, $input2, $range, undef, $precision, $directory );
	my ( @PL1, @PL2, @intensities1, @intensities2 );
#### READ THE PEAKLISTS FROM THE XML HASH
	## SPECTRUM 1
	my $scan1peaks = $scan1obj->get_peaklist;    ## is a array reference / only the peaaks
	## print Dumper($scan1peaks);
	## get the basename of input 1
	$input1 =~ s/\.dta//;
	### generate the header (same as for a dta)
	my $mr1 = $scan1obj->Mr;
	my $z1  = $scan1obj->FT_charge;
	unless ( $mr1 or $z1 )
	{
		die "Charge or Mr not set!\n";
	}
	my $datheader1 = "$mr1\t$z1\n";
	my $specline   = 0;
	foreach my $peakpair (@$scan1peaks)
	{    ## peakpair is an array ref
		my $intensity = $peakpair->[1];
		my $mz        = $peakpair->[0];
		$specline++;
		unless ( defined( $intensity && $mz ) )
		{
			die "Error Specline $specline is not valid\n";
		}
		## check if the intensity is above the threshold
		if ( $intensity >= $threshold )
		{
			push @PL1,          $peakpair;
			push @intensities1, $peakpair->[1];
		}
	}
	my $numpeaksPL1 = scalar(@PL1);
	$report && print "#Ions Scan1 (after sorting intens. > $threshold): $numpeaksPL1\n";
	## SPECTRUM 2
	my $scan2peaks = $scan2obj->get_peaklist;    ## is a array reference / only the peaaks
	                                             #print Dumper($scan1peaks);
	## get the basename of input 1
	$input2 =~ s/\.dta//;
	### generate the header (same as for a dta)
	my $mr2 = $scan2obj->Mr;
	my $z2  = $scan2obj->FT_charge;
	unless ( $mr2 or $z2 )
	{
		die "Charge or Mr not set!\n";
	}
	my $datheader2 = "$mr2\t$z2\n";
	$specline = 0;
	foreach my $peakpair (@$scan2peaks)
	{    ## peakpair is an array ref
		my $intensity = $peakpair->[1];
		my $mz        = $peakpair->[0];
		$specline++;
		unless ( defined( $intensity && $mz ) )
		{
			die "Error Specline $specline is not valid\n";
		}
		## check if the intensity is above the threshold
		if ( $intensity >= $threshold )
		{
			push @PL2,          $peakpair;
			push @intensities2, $peakpair->[1];
		}
	}
	my $numpeaksPL2 = scalar(@PL2);

	#$report && print "Number of ions in spectrum2: $specline, after sorting (above intensity threshold:$threshold): $numpeaksPL2\n";
	$report && print "#Ions Scan2 (after sorting intens. > $threshold): $numpeaksPL2\n";

	# EOF LOADING THE PEAKS OF SPECTRUM2
	## set the path and basename for the result files
	$outfile = File::Spec->catfile( $directory, ( join "_", basename($input1), basename($input2) ) );
	if ( $seen{$outfile}++ )
	{
		return 0;
	}
	my $max1 = Statistics::max( \@intensities1 );
	my $max2 = Statistics::max( \@intensities2 );

	# Filtering of the peaklists by dynamic range
	if ($dynamic_range)
	{
		my $threshold1 = $max1 / $dynamic_range;
		my $threshold2 = $max2 / $dynamic_range;
		my @tmp1       = ();
		foreach my $peak (@PL1)
		{
			if ( intensity($peak) > $threshold1 )
			{
				push @tmp1, $peak;
			}
		}
		my @tmp2 = ();
		foreach my $peak (@PL2)
		{
			if ( intensity($peak) > $threshold2 )
			{
				push @tmp2, $peak;
			}
		}
		@PL1         = @tmp1;
		@PL2         = @tmp2;
		$numpeaksPL1 = scalar(@PL1);
		$numpeaksPL2 = scalar(@PL2);

		#$report && print "Number of ions in Spectrum 1 after filtereng by dynamic range (): $numpeaksPL1\n";
		#$report && print "Number of ions in Spectrum 2 after filtereng by dynamic range (): $numpeaksPL2\n";
	}
	## Rescaleing the Intensities ##############################################################################
	if ($scaleintensity)
	{
		my $topten1 = 0;
		my $topten2 = 0;
		if ( $scaleby eq 'top5' )
		{
			$topten1 = Statistics::topn( \@intensities1, 5 );
			!$weboutput
			  && print "mean intensity of  5\% most intense peaks for $input1: $topten1<br>\n";
			$topten2 = Statistics::topn( \@intensities2, 5 );
			!$weboutput && print "mean intensity of  5\% most intense peaks for $input2: $topten2<br>\n";
		} else
		{
			$topten1 = $max1;
			$topten2 = $max2;
			!$weboutput && print "<br>maximum intensity spectrum1: $topten1<br>\n";
			!$weboutput && print "maximum intensity spectrum2: $topten2<br>\n";
		}
		foreach my $peak (@PL1)
		{
			scaleintensity( $peak, $topten1, 100 );
			$superverbose && print mz($peak), "\t", intensity($peak), "\n";
		}
		$superverbose && print "\nspectrum2\n";
		foreach my $peak (@PL2)
		{
			scaleintensity( $peak, $topten2, 100 );
			$superverbose && print mz($peak), "\t", intensity($peak), "\n";
		}
	}
	## EOF Rescaleing the Intensities ##############################################################################
	## GENERATE ARRRAYS FOR THE PLOTS -> orginal peaks#####
	my @PL1_plot = @PL1;
	my @PL2_plot = @PL2;
	####
	# Modify the peaklist to include mispicked transformed peaks for peakmatching
	# by the deconvolution/deistoping algorithm, includes -C13 and -2*C13 peak
	if ( ($highresolutionms2) && !($addcommon2xlink) )
	{
		my @templist    = ();
		my @orginallist = ();
		print "Add -1C13 shift to every peak of the light peaklist\n";
		foreach my $peak (@PL1)
		{
			my $oldmz;
			my $newmz1;
			my $newmz2;
			my $newmz3;
			my $newmz4;
			$oldmz = mz($peak);
			my $intensity = intensity($peak);

			#print $C13shift;
			$newmz1 = $oldmz - $C13shift;
			$newmz2 = $oldmz - 2 * $C13shift;
			$newmz3 = $oldmz + $C13shift;
			$newmz4 = $oldmz + 2 * $C13shift;

			#print $oldmz." ".$newmz1."\n";
			my @orginalpeak = ( $oldmz,  $intensity );
			my @newpeak1    = ( $newmz1, $intensity );
			my @newpeak2    = ( $newmz2, $intensity );
			my @newpeak3    = ( $newmz3, $intensity );
			my @newpeak4    = ( $newmz4, $intensity );

			#print "oldmz: $oldmz orginalmz:$oldmz\n";
			push( @templist, \@orginalpeak, \@newpeak1, \@newpeak2, \@newpeak3, \@newpeak4 );
		}
		@PL1 = @templist;
	}
	## Saving the normalized files into the  XML as dta base 64 encoded#####
	my @peaklist_light = sort { mz($a) <=> mz($b) } @PL1;
	my @peaklist_heavy = sort { mz($a) <=> mz($b) } @PL2;
	my $normalized_dta1 = $input1 . ".dta";
	my $normalized_dta2 = $input2 . ".dta";

	#open PL1norm, ">$normalized_dta1" or die $!;
	#open PL2norm, ">$normalized_dta2" or die $!;
	#print PL1norm $datheader1, "\n";
	#print PL2norm $datheader2, "\n";
	my $plstring1 = $datheader1;
	foreach my $peak (@peaklist_light)
	{
		$plstring1 .= join "\t", @$peak, "\n";

		#print PL1norm "\n";
	}

	#close(PL1norm);
	my $plstring2 = $datheader2;
	foreach my $peak (@peaklist_heavy)
	{
		$plstring2 .= join "\t", @$peak, "\n";

		#print PL2norm join "\t", @$peak;
		#print PL2norm "\n";
	}

	#close(PL2norm);
	### Save the normalized spectra to the xml // done after peakmatching
	#save_spectra_to_XML_base64enc($normalized_dta1,$scantype1,$XMloutfile,$plstring1);
	#save_spectra_to_XML_base64enc($normalized_dta2,$scantype2,$XMloutfile,$plstring2);
	# For Debuging
	#$verbose =1;
	if ($verbose)
	{

		# print the peaklists
		print "\n#### PEAKLIST LIGHT: ####\n";
		foreach my $peak (@peaklist_light)
		{
			printpeak($peak);
		}
		print "#### END OF PEAKLIST LIGHT ###\n";
		print "\n#### PEAKLIST HEAVY: ####\n";
		foreach my $peak (@peaklist_heavy)
		{
			printpeak($peak);
		}
		print "#### END OF PEAKLIST HEAVY ###\n";
	}
	## PEAK MATCHING AND SEPARATION OF COMMON and Cross-linking peaks
	## search common peaks list1
	my ( @meancommon_mz, @meancommon_int, @meancommon_delta, @xlink, $pair, @commonpeaks_light, @commonpeaks_heavy );
	if ($verbose) { print "\n#### MATCHING of common ions of light to heavy spectrum: ####\n" }
	@commonpeaks_light = matchpeaks( \@peaklist_light, \@peaklist_heavy, $tolerance, $lower_ratio, $upper_ratio, 0 );
	if ( $average_common_peaks || $verbose )
	{
		if ($verbose) { print "\n#### MATCHING of common ions of heavy to light spectrum: ####\n" }
		@commonpeaks_heavy = matchpeaks( \@peaklist_heavy, \@peaklist_light, $tolerance, $lower_ratio, $upper_ratio, 0 );
	}
	## check if enough common peaks were found
	unless ( scalar(@commonpeaks_light) >= $mincommonions )
	{
		if ($weboutput)
		{
			print "<div class=\"error\"><br>",
			  "number of common peaks is smaller than the required minimum of $mincommonions: ",
			  scalar(@commonpeaks_light), " in $input1 and $input2\: #peaks1 > ", $threshold, ": ",
			  scalar(@peaklist_light), " | #peaks2 > ", $threshold, ": ",
			  scalar(@peaklist_heavy),, "->ignored</div>";
			$rejected++;
		} else
		{
			warn "number of common peaks is smaller than the required minimum of $mincommonions: ", scalar(@commonpeaks_light), " in $input1 and $input2\: #peaks1 > ", $threshold, ": ", scalar(@peaklist_light), " | #peaks2 > ", $threshold, ": ", scalar(@peaklist_heavy),, "->ignored\n";
			$rejected++;
		}
		return 0;
	}
#### SAVE THE NORMALIZED SPECTRA TO THE XML FILE
	### Save the normalized spectra to the xml
	#$report && print "Saving normalized spectrum $normalized_dta1\n";
	save_spectra_to_XML_base64enc( $normalized_dta1, $scantype1, $XMloutfile, $plstring1 );

	#$report && print "Saving normalized spectrum $normalized_dta2\n";
	save_spectra_to_XML_base64enc( $normalized_dta2, $scantype2, $XMloutfile, $plstring2 );
#### HERE THE LINE IS PRINTED TO THE ISOTOPEPAIR LIST
## add the rttimes and the mz values

	print $filelisthandle "$outfile\t$isotopicshift4xlinkions\t", $scantype1, "_", "$scantype2\t$rttimes\t$mzofscans\n";
	my $spline = "$outfile\t$isotopicshift4xlinkions\t" . $scantype1 . "_" . "$scantype2\t$rttimes\t$mzofscans\n";

	#print "Writing line to file: $filelisthandle: $spline";
#### PEAK MATCHING AND SEPARATION OF CROSSLINKER PEAKS $spline;
	#copy the orginal chargestate
	my $charge_orginal = $charge;
	if ($highresolutionms2)
	{

		# if highresolution ms2 are used set the charge state to 1
		# all peaks should be deisotoped and should have charge 1
		$charge = 1;
		$report && print "HIGH RESOLUTION MS2 charge is set to 1\n";
	}

	#  search xlink peaks list
	my ( @xlinkpeaks_light, @xlinkpeaks_heavy );
	### Lightonly searches: replace the ions of the light precursor to the xlinkpeaks if light or heavyonly is compared, no peakmatching of x-linker ions
	### all charge states are then also 0, which makes sense
	### the common and xlinker ion lists are then equal
	if ($addcommon2xlink)
	{
		@xlinkpeaks_light = @commonpeaks_light;
		@xlinkpeaks_heavy = @commonpeaks_heavy;

		# SET THE CHARGE STATE TO 1 in high res mode, peaks are already deconvoluted
		if ($highresolutionms2)
		{
			$report && print "High res ms2: charge of x-linker ions is set to 1";
			foreach my $peak (@xlinkpeaks_light)
			{
				setcharge( $peak, 1 );
			}
			foreach my $peak (@xlinkpeaks_heavy)
			{
				setcharge( $peak, 1 );
			}
		}
	} else
	{
		### Search for x-linker ions for all possible chage states, minchargestate=1, max = precursor chargestate
		for my $chargestate ( $minchargestate .. $charge )
		{
			$verbose && print "searching for +$chargestate charged x-link-ions <br>\n";
			my @transformedpeaks = ();
			## makes a transformed peaklist where the heavy peaks are shifted to the left by the isodiff/chargestates
			## if these are compared with the light list they should match if a xlink is there
			foreach my $heavypeak (@peaklist_heavy)
			{
				push @transformedpeaks, [ $heavypeak->[0] - ( $isotopediff / $chargestate ), $heavypeak->[1] ];
			}
			push @xlinkpeaks_light, matchpeaks( \@peaklist_light, \@transformedpeaks, $tolerancexl, $lower_ratio, $upper_ratio, $chargestate );
			$verbose && ( push @xlinkpeaks_heavy, matchpeaks( \@transformedpeaks, \@peaklist_light, $tolerancexl, $lower_ratio, $upper_ratio, $chargestate ) );
		}
	}

	# END OF PEAKMATCHING
	if ($verbose)
	{

		# print the matched peaks
		print "\n#### Common PEAKS MATCHED LIGHT: ####\n";
		foreach my $peak (@commonpeaks_light)
		{
			printpeak($peak);
		}
		print "#### END OF MATCHED COMMON PEAKLIST LIGHT ####\n";
		print "\n#### Common PEAKS MATCHED HEAVY:####\n";
		foreach my $peak (@commonpeaks_heavy)
		{
			printpeak($peak);
		}
		print "#### END OF MATCHED COMMON PEAKLIST HEAVY ###\n";
	}
	my $commonmin = 0;
	my $xlinkmin  = 0;

	# sort out if nhighest param is used
	if ( scalar(@commonpeaks_light) > $usenhighest )
	{
		my @tmp = sort { intensity($b) <=> intensity($a) } @commonpeaks_light;
		$commonmin = intensity( $tmp[$usenhighest] );
	}
	if ( scalar(@xlinkpeaks_light) > $usenhighest )
	{
		my @tmp = sort { intensity($b) <=> intensity($a) } @xlinkpeaks_light;
		$xlinkmin = intensity( $tmp[$usenhighest] );
	}

	# AVERAGE the common peaks from the light and heavy spectrum
	if ($average_common_peaks)
	{
		for $i ( 0 .. $#commonpeaks_light )
		{
			$meancommon_mz[$i] = ( mz( $commonpeaks_light[$i] ) + mz( $commonpeaks_heavy[$i] ) ) / 2;
			$verbose && ( $meancommon_delta[$i] = sprintf( "%.4f", ( mz( $commonpeaks_light[$i] ) - mz( $commonpeaks_heavy[$i] ) ) ) );
			$meancommon_int[$i] = ( intensity( $commonpeaks_light[$i] ) + intensity( $commonpeaks_heavy[$i] ) ) / 2;
		}
	}
	$verbose && print "Common peaks tolerance=$tolerance, required peakratio = $peakratio\n";
	$verbose && print "mz1\tintensity1\tmz2\tintensity2\tmeanmz\tmeanintensity\tdelta\n";
	my $ncommonpeaks = 0;
	## GENERATE THE FILES FOR THE COMMON AND XLINKER Files
	my $commonheader = "$id\n$mz\n$charge_orginal\n";
	my $xlheader     = "$id\n$mz\n$charge_orginal\n";

	#print COMMONPEAKS "$id\n$mz\n$charge_orginal\n";
	#print XLINKERPEAKS "$id\n$mz\n$charge_orginal\n";
	my $commoncontent = $commonheader;
	my $xlcontent     = $xlheader;
	for ( $i = 0 ; $i <= $#commonpeaks_light ; $i += 1 )
	{
		if ( intensity( $commonpeaks_light[$i] ) >= $commonmin )
		{
			$ncommonpeaks++;
			$verbose && print mz( $commonpeaks_light[$i] ), "\t", intensity( $commonpeaks_light[$i] ), "\t\t", mz( $commonpeaks_heavy[$i] ), "\t", intensity( $commonpeaks_heavy[$i] ), "\t", $meancommon_mz[$i], "\t", $meancommon_int[$i], "\t", $meancommon_delta[$i], "\n";
			if ($average_common_peaks)
			{
				$commoncontent .= $meancommon_mz[$i] . "\t" . $meancommon_int[$i] . "\t0" . "\n";
			} else
			{
				$commoncontent .= mz( $commonpeaks_light[$i] ) . "\t" . intensity( $commonpeaks_light[$i] ) . "\t0" . "\n";
			}
		}
	}
	if ($verbose)
	{
		print "xlinker peaks\n";
		for ( $i = 0 ; $i <= $#xlinkpeaks_light ; $i += 1 )
		{
			print mz( $xlinkpeaks_light[$i] ), "\t", intensity( $xlinkpeaks_light[$i] ), "\tcharge: ", charge( $xlinkpeaks_light[$i] ), "\t";
			print mz( $xlinkpeaks_heavy[$i] ), "\t", intensity( $xlinkpeaks_heavy[$i] ), "\tcharge: ", charge( $xlinkpeaks_heavy[$i] ), "\t";
			print "\n";
		}
	}
	my $nxlinkpeaks = 0;
	foreach my $xlinkpeak (@xlinkpeaks_light)
	{
		if ( intensity($xlinkpeak) >= $xlinkmin )
		{
			$nxlinkpeaks++;
			### Deconvolute if the option is selected
			if ($deconvolutexlinkions)
			{
				my $deconvolutedpeak = deconvolute($xlinkpeak);
				$xlcontent .= mz($deconvolutedpeak) . "\t" . intensity($deconvolutedpeak) . "\t" . charge($deconvolutedpeak) . "\n";
			} else
			{
				$xlcontent .= mz($xlinkpeak) . "\t" . intensity($xlinkpeak) . "\t" . charge($xlinkpeak) . "\n";
			}
		}
	}
	#### SAVE THE COMMON AND XL SPECTRA TO THE XML FILE
	# Open the filehandles for common and xlinker peklists
	my $fncommon = $outfile . "_common.txt";
	my $fnxl     = $outfile . "_xlinker.txt";
	save_spectra_to_XML_base64enc( $fncommon, "common",  $XMloutfile, $commoncontent );
	save_spectra_to_XML_base64enc( $fnxl,     "xlinker", $XMloutfile, $xlcontent );
	### PRINT OVERLAY SPECTRUM ################################################
	if ($printspectra)
	{
		require specplot;
		my $specfile = join "", $outfile, "_specplot.png";
		foreach my $peakpair (@PL2_plot)
		{
			$peakpair->[1] *= -1;
		}
		foreach my $peakpair (@commonpeaks_heavy)
		{
			$peakpair->[1] *= -1;
		}
		foreach my $peakpair (@xlinkpeaks_heavy)
		{
			$peakpair->[1] *= -1;
		}
		my $specobj1 = specplot->new();

		#$specobj1->setcolor( "grey", "grey", "green", "red" );
		$specobj1->plotdata( $minionsize, $maxionsize, [ "grey", "grey", "green", "red" ], \@PL1_plot, \@PL2_plot, \@commonpeaks_light, \@xlinkpeaks_light );
		$specobj1->drawlegend( 200, 10, { basename($outfile), "black" } );
		my %colorhash = (
						  "positive: light spectrum peaks" => "grey",
						  "negative: heavy spectrum peaks" => "grey",
						  "matched common-spectrum peaks"  => "green",
						  "matched xlink-spectrum peaks"   => "red",
		);
		$specobj1->drawlegend( 600, 20, \%colorhash );
		$specobj1->printimage($specfile);
	}
	if ($weboutput)
	{
		print '<br><b>', basename($input1), ",", basename($input2), ': sorting peaks into common peaks and x-linker peaks:</b><br>';
		print 'number of peaks spectrum1  > ', $threshold, ': ',
		  scalar(@peaklist_light), '<br>number of peaks spectrum2 > ',
		  $threshold, ': ', scalar(@peaklist_heavy),
		  '<br>number of  common-peaks > ', sprintf( "%.2f", $commonmin ), ' ',
		  $ncommonpeaks, '<br>number of xlink-peaks > ',
		  sprintf( "%.2f", $xlinkmin ), ' ', $nxlinkpeaks, '<br>';
	} else
	{
		print "#peaks1 > ", $threshold, ": ", scalar(@peaklist_light), "|  #peaks2 > ", $threshold, ": ", scalar(@peaklist_heavy), "| #commonpeaks > ", sprintf( "%.2f", $commonmin ), "% ", $ncommonpeaks, "| #xlinkpeaks > ", sprintf( "%.2f", $xlinkmin ), "% ", $nxlinkpeaks, "<br>\n";
	}

	#exit;
	return 1;
}

#===  FUNCTION  ================================================================
#  NAME:        deconvolute
#  PURPOSE:     deconvolute a cross-linker peak
#  DESCRIPTION: deconvolute a cross-linker peak
#  PARAMETERS:  $peak
#  RETURNS:     deconvoluted peak
#===============================================================================
sub deconvolute
{
	my $peak        = shift;
	my $mz          = mz($peak);
	my $int         = intensity($peak);
	my $z           = charge($peak);
	my $nominalmass = $mz * $z - $z * $Hatom;
	my $z1mz        = $nominalmass + $Hatom;

	#print "PEAK: mz: $mz, int: $int, z:$z, nominal mass: $nominalmass z1mz: $z1mz\n";
	my $z1peak = [ $z1mz, $int, 1 ];
	return $z1peak;
}

#===  FUNCTION  ================================================================
#  NAME:        save_spectra_to_XML_base64enc
#  PURPOSE:     save_spectra_to_XML_base64enc
#  DESCRIPTION: describtion
#  PARAMETERS:  $filename (used as attr. in xml), $typeattrybute,
#  RETURNS:     return
#===============================================================================
sub save_spectra_to_XML_base64enc
{
	my $filename      = shift;    ##specfilename
	my $typeattribute = shift;
	my $xmlfilename   = shift;    ## ref to filehandle *SPECFILE
	my $buffertosave  = shift;
	my $text;
	## load the file
	## encode the content
	my $encoded_text = encode_base64($buffertosave);
	## create an xml element
	my $filenamebn = basename($filename);
	$text = "<spectrum filename=\"$filenamebn\" type=\"$typeattribute\">";
	$text .= $encoded_text;
	$text .= "</spectrum>\n";
	### Store the spectra in the xml file
	printToXML( $xmlfilename, $text );
	return;
}
## added to save the spectra in the xmlresult
## Save the spectra in a separate file
## print spectrum light/heavy/common/xlink into one tag base 64 encoded
sub printToXML
{
	my $xlinkfile = shift;
	my $content   = shift;
	print $xlinkfile $content;
}

sub process_matched_scans
{
	my $matchpair               = shift;
	my $filelisthandle          = shift;
	my $id                      = $matchpair->[0];
	my $mz                      = $matchpair->[1];
	my $charge                  = $matchpair->[2];
	my $input1                  = $matchpair->[3];
	my $input2                  = $matchpair->[4];
	my $scantype1               = $matchpair->[5];
	my $scantype2               = $matchpair->[6];
## use the dta spectra	
	my $rttimes        = $matchpair->[8];        # 1234:1238 ## rt in seconds
	my $mzofscans		= $matchpair->[9];		 # 531.53973389:535.55877686
	

	my $isotopicshift4xlinkions = 0;
	$addcommon2xlink = ();

	if ( ( $scantype1 eq 'heavy' ) && ( $scantype2 eq 'heavy' ) )
	{
		$isotopicshift4xlinkions = $isotopediff;
	}
	if ( $scantype1 eq $scantype2 )
	{
		$addcommon2xlink = 1;
	}
	$verbose && print "id: $id mz: $mz charge: $charge dta-file1: $input1 dta-file2: $input2 scan1: $scantype1 scan2: $scantype2\n";
	if ( ( $input1 eq "void" ) || ( $input2 eq "void" ) )
	{
		next;
	}
	unless ( ( -e $input1 ) && ( -e $input2 ) )
	{
		warn "cannot find $input1 and $input2\n";
		$notfound++;
		next;
	}
	$xcorrelation && xcorr::xcorrelation( $input1, $input2, $range, undef, $precision, $directory );
	my ( @PL1, @PL2, @intensities1, @intensities2 );
	open XML1, "<$input1" or die "cannot open dat file $input1 $!\n";
	$input1 =~ s/\.dta//;
	my $datheader1 = <XML1>;
	chomp($datheader1);    #get dta header containing m/z and charge

	#LOADING THE PEAKS OF SPECTRUM1 #########################################################
	#SORTS OUT THE PEAKS THAT ARE BELOW THE INTENSITY THRESHOLD cp_threshold ################
	#Each Peak is splitted into an mz/intensity array and a REFERENCE is pushed into @PL1 ###
	my $specline = 0;
	while (<XML1>)
	{
		chomp;
		my @tmp = split;
		$specline++;
		print "<div class=\"error\"><br>spectrum line $specline \"$_\" is not valid. I read m/z = $tmp[0] intensity= $tmp[1] </div>\n"
		  unless ( defined( $tmp[0] && $tmp[1] ) );
		if ( $tmp[1] >= $threshold )
		{
			push @PL1,          \@tmp;
			push @intensities1, $tmp[1];
		}
	}
	print "<div class=\"error\"><br>warning: read only $specline peakpairs from spectrum </div>\n"
	  unless ( $specline > 10 );
	close(XML1);
	my $numpeaksPL1 = scalar(@PL1);
	$report && print "Number of ions in spectrum1: $specline, after sorting (above threshold:$threshold): $numpeaksPL1\n";

	# EOF LOADING THE PEAKS OF SPECTRUM1
	#LOADING THE PEAKS OF SPECTRUM2 #########################################################
	#SORTS OUT THE PEAKS THAT ARE BELOW THE INTENSITY THRESHOLD cp_threshold ################
	#Each Peak is splitted into an mz/intensity array and a REFERENCE is pushed into @PL2 ###
	$specline = 0;
	open XML2, "<$input2" or die "cannot open dat file $input2 $!\n";
	$input2 =~ s/\.dta//;
	my $datheader2 = <XML2>;
	chomp($datheader2);    #get dta header containing m/z and charge
	while (<XML2>)
	{
		chomp;
		my @tmp = split;
		$specline++;
		print "<div class=\"error\"><br>spectrum line $specline \"$_\" is not valid. I read m/z = $tmp[0] intensity= $tmp[1] </div>\n"
		  unless ( defined( $tmp[0] && $tmp[1] ) );
		if ( $tmp[1] >= $threshold )
		{
			push @PL2,          \@tmp;
			push @intensities2, $tmp[1];
		}
	}
	print "<div class=\"error\"><br>warning: read only $specline peakpairs from spectrum </div>\n" unless ( $specline > 10 );
	close(XML2);
	my $numpeaksPL2 = scalar(@PL2);
	$report && print "Number of ions in spectrum2: $specline, after sorting (above threshold:$threshold): $numpeaksPL2\n";

	# EOF LOADING THE PEAKS OF SPECTRUM2
	## set the path and basename for the result files
	$outfile = File::Spec->catfile( $directory, ( join "_", basename($input1), basename($input2) ) );
	if ( $seen{$outfile}++ )
	{
		return 0;
	}
	my $max1 = Statistics::max( \@intensities1 );
	my $max2 = Statistics::max( \@intensities2 );

	# Filtering of the peaklists by dynamic range
	if ($dynamic_range)
	{
		my $threshold1 = $max1 / $dynamic_range;
		my $threshold2 = $max2 / $dynamic_range;
		my @tmp1       = ();
		foreach my $peak (@PL1)
		{
			if ( intensity($peak) > $threshold1 )
			{
				push @tmp1, $peak;
			}
		}
		my @tmp2 = ();
		foreach my $peak (@PL2)
		{
			if ( intensity($peak) > $threshold2 )
			{
				push @tmp2, $peak;
			}
		}
		@PL1         = @tmp1;
		@PL2         = @tmp2;
		$numpeaksPL1 = scalar(@PL1);
		$numpeaksPL2 = scalar(@PL2);
		$report && print "Number of ions in Spectrum 1 after filtereng by dynamic range (): $numpeaksPL1\n";
		$report && print "Number of ions in Spectrum 2 after filtereng by dynamic range (): $numpeaksPL2\n";
	}
	## Rescaleing the Intensities ##############################################################################
	if ($scaleintensity)
	{
		my $topten1 = 0;
		my $topten2 = 0;
		if ( $scaleby eq 'top5' )
		{
			$topten1 = Statistics::topn( \@intensities1, 5 );
			!$weboutput
			  && print "mean intensity of  5\% most intense peaks for $input1: $topten1<br>\n";
			$topten2 = Statistics::topn( \@intensities2, 5 );
			!$weboutput
			  && print "mean intensity of  5\% most intense peaks for $input2: $topten2<br>\n";
		} else
		{
			$topten1 = $max1;
			$topten2 = $max2;
			!$weboutput
			  && print "<br>maximum intensity spectrum1: $topten1<br>\n";
			!$weboutput && print "maximum intensity spectrum2: $topten2<br>\n";
		}
		foreach my $peak (@PL1)
		{
			scaleintensity( $peak, $topten1, 100 );
			$superverbose && print mz($peak), "\t", intensity($peak), "\n";
		}
		$superverbose && print "\nspectrum2\n";
		foreach my $peak (@PL2)
		{
			scaleintensity( $peak, $topten2, 100 );
			$superverbose && print mz($peak), "\t", intensity($peak), "\n";
		}
	}
	## EOF Rescaleing the Intensities ##############################################################################
	## GENERATE ARRRAYS FOR THE PLOTS -> orginal peaks#####
	my @PL1_plot = @PL1;
	my @PL2_plot = @PL2;
	####
	# Modify the peaklist to include mispicked transformed peaks for peakmatching
	# by the deconvolution/deistoping algorithm, includes -C13 and -2*C13 peak
	if ( ($highresolutionms2) && !($addcommon2xlink) )
	{
		my @templist    = ();
		my @orginallist = ();
		print "Add -1C13 shift to every peak of the light peaklist\n";
		foreach my $peak (@PL1)
		{
			my $oldmz;
			my $newmz1;
			my $newmz2;
			my $newmz3;
			my $newmz4;
			$oldmz = mz($peak);
			my $intensity = intensity($peak);

			#print $C13shift;
			$newmz1 = $oldmz - $C13shift;
			$newmz2 = $oldmz - 2 * $C13shift;
			$newmz3 = $oldmz + $C13shift;
			$newmz4 = $oldmz + 2 * $C13shift;

			#print $oldmz." ".$newmz1."\n";
			my @orginalpeak = ( $oldmz,  $intensity );
			my @newpeak1    = ( $newmz1, $intensity );
			my @newpeak2    = ( $newmz2, $intensity );
			my @newpeak3    = ( $newmz3, $intensity );
			my @newpeak4    = ( $newmz4, $intensity );

			#print "oldmz: $oldmz orginalmz:$oldmz\n";
			push( @templist, \@orginalpeak, \@newpeak1, \@newpeak2, \@newpeak3, \@newpeak4 );
		}
		@PL1 = @templist;
	}
	## Saving the normalized files into the spectra folder e.g O09-08xxx_p #####
	my @peaklist_light = sort { mz($a) <=> mz($b) } @PL1;
	my @peaklist_heavy = sort { mz($a) <=> mz($b) } @PL2;
	my $normalized_dta1 = $input1 . "_normalized.dta";
	my $normalized_dta2 = $input2 . "_normalized.dta";
	open PL1norm, ">$normalized_dta1" or die $!;
	open PL2norm, ">$normalized_dta2" or die $!;
	print PL1norm $datheader1, "\n";
	print PL2norm $datheader2, "\n";

	foreach my $peak (@peaklist_light)
	{
		print PL1norm join "\t", @$peak;
		print PL1norm "\n";
	}
	close(PL1norm);
	foreach my $peak (@peaklist_heavy)
	{
		print PL2norm join "\t", @$peak;
		print PL2norm "\n";
	}
	close(PL2norm);

	# Copy the normalized spectra also to the spectrumdirectory ###################
	# The normalized spectra are copied into the xQuest result directory, not in the
	# e.g O09-07106_shortlist, not in the O09-07106_shortlistdir where the processend
	# spectra are copied to.
	if ( -e $resultdir )
	{
		my $destination_dtafilename1 = "$resultdir\/" . basename($input1) . ".dta";
		my $destination_dtafilename2 = "$resultdir\/" . basename($input2) . ".dta";
		print "Copying $normalized_dta1 to $destination_dtafilename1\n";
		print "Copying $normalized_dta2 to $destination_dtafilename2\n";
		copy( $normalized_dta1, $destination_dtafilename1 ) or die $!;
		copy( $normalized_dta2, $destination_dtafilename2 ) or die $!;
	} else
	{
		print $resultdir, "does not exist!";
	}

	# For Debuging
	if ($verbose)
	{

		# print the peaklists
		print "\n#### PEAKLIST LIGHT: ####\n";
		foreach my $peak (@peaklist_light)
		{
			printpeak($peak);
		}
		print "#### END OF PEAKLIST LIGHT ###\n";
		print "\n#### PEAKLIST HEAVY: ####\n";
		foreach my $peak (@peaklist_heavy)
		{
			printpeak($peak);
		}
		print "#### END OF PEAKLIST HEAVY ###\n";
	}
	## PEAK MATCHING AND SEPARATION OF COMMON PEAKS
	## search common peaks list1
	my ( @meancommon_mz, @meancommon_int, @meancommon_delta, @xlink, $pair, @commonpeaks_light, @commonpeaks_heavy );
	if ($verbose) { print "\n#### MATCHING of common ions of light to heavy spectrum: ####\n" }
	@commonpeaks_light = matchpeaks( \@peaklist_light, \@peaklist_heavy, $tolerance, $lower_ratio, $upper_ratio, 0 );
	if ( $average_common_peaks || $verbose )
	{
		if ($verbose) { print "\n#### MATCHING of common ions of heavy to light spectrum: ####\n" }
		@commonpeaks_heavy = matchpeaks( \@peaklist_heavy, \@peaklist_light, $tolerance, $lower_ratio, $upper_ratio, 0 );
	}
	## check if enough common peaks were found
	unless ( scalar(@commonpeaks_light) >= $mincommonions )
	{
		if ($weboutput)
		{
			print "<div class=\"error\"><br>",
			  "number of common peaks is smaller than the required minimum of $mincommonions: ",
			  scalar(@commonpeaks_light), " in $input1 and $input2\: #peaks1 > ", $threshold, ": ",
			  scalar(@peaklist_light), " | #peaks2 > ", $threshold, ": ",
			  scalar(@peaklist_heavy),, "->ignored</div>";
			$rejected++;
		} else
		{
			warn "number of common peaks is smaller than the required minimum of $mincommonions: ", scalar(@commonpeaks_light), " in $input1 and $input2\: #peaks1 > ", $threshold, ": ", scalar(@peaklist_light), " | #peaks2 > ", $threshold, ": ", scalar(@peaklist_heavy),, "->ignored\n";
			$rejected++;
		}
		return 0;
	}

	# Open the filehandles for common and xlinker peklists
	open COMMONPEAKS,  ">$outfile" . "_common.txt";
	open XLINKERPEAKS, ">$outfile" . "_xlinker.txt";
	print $filelisthandle "$outfile\t$isotopicshift4xlinkions\t", $scantype1, "_", "$scantype2\t$rttimes\t$mzofscans\n";
	## PEAK MATCHING AND SEPARATION OF CROSSLINKER PEAKS
	# copy the orginal chargestate
	my $charge_orginal = $charge;
	if ($highresolutionms2)
	{

		# if highresolution ms2 are used set the charge state to 1
		# all peaks should be deisotoped and should have charge 1
		$charge = 1;
		print "HIGH RESOLUTION MS2 charge is set to 1\n";
	}

	#  search xlink peaks list
	my ( @xlinkpeaks_light, @xlinkpeaks_heavy );
	### replace the ions of the light precursor to the xlinkpeaks if light or heavyonly is compared, no peakmatching of x-linker ions
	### all charge states are then also 0, which makes sense
	### the common and xlinker ion lists are then equal
	if ($addcommon2xlink)
	{
		@xlinkpeaks_light = @commonpeaks_light;
		@xlinkpeaks_heavy = @commonpeaks_heavy;

		# SET THE CHARGE STATE TO 1 in high res mode, peaks are already deconvoluted
		if ($highresolutionms2)
		{
			print "High res ms2: charge of x-linker ions is set to 1";
			foreach my $peak (@xlinkpeaks_light)
			{
				setcharge( $peak, 1 );
			}
			foreach my $peak (@xlinkpeaks_heavy)
			{
				setcharge( $peak, 1 );
			}
		}
	} else
	{
		### Search for x-linker ions for all possible chage states, minchargestate=1, max = precursor chargestate
		for my $chargestate ( $minchargestate .. $charge )
		{
			$verbose && print "searching for +$chargestate charged x-link-ions <br>\n";
			my @transformedpeaks = ();
			## makes a transformed peaklist where the heavy peaks are shifted to the left by the isodiff/chargestates
			## if these are compared with the light list they should match if a xlink is there
			foreach my $heavypeak (@peaklist_heavy)
			{
				push @transformedpeaks, [ $heavypeak->[0] - ( $isotopediff / $chargestate ), $heavypeak->[1] ];
			}
			push @xlinkpeaks_light, matchpeaks( \@peaklist_light, \@transformedpeaks, $defhash->{'xlink_ms2tolerance'}, $lower_ratio, $upper_ratio, $chargestate );
			$verbose && ( push @xlinkpeaks_heavy, matchpeaks( \@transformedpeaks, \@peaklist_light, $defhash->{'xlink_ms2tolerance'}, $lower_ratio, $upper_ratio, $chargestate ) );
		}
	}

	# END OF PEAKMATCHING
	if ($verbose)
	{

		# print the matched peaks
		print "\n#### Common PEAKS MATCHED LIGHT: ####\n";
		foreach my $peak (@commonpeaks_light)
		{
			printpeak($peak);
		}
		print "#### END OF MATCHED COMMON PEAKLIST LIGHT ####\n";
		print "\n#### Common PEAKS MATCHED HEAVY:####\n";
		foreach my $peak (@commonpeaks_heavy)
		{
			printpeak($peak);
		}
		print "#### END OF MATCHED COMMON PEAKLIST HEAVY ###\n";
	}
	my $commonmin = 0;
	my $xlinkmin  = 0;

	# sort out if nhighest param is used
	if ( scalar(@commonpeaks_light) > $usenhighest )
	{
		my @tmp = sort { intensity($b) <=> intensity($a) } @commonpeaks_light;
		$commonmin = intensity( $tmp[$usenhighest] );
	}
	if ( scalar(@xlinkpeaks_light) > $usenhighest )
	{
		my @tmp = sort { intensity($b) <=> intensity($a) } @xlinkpeaks_light;
		$xlinkmin = intensity( $tmp[$usenhighest] );
	}

	# AVERAGE the common peaks from the light and heavy spectrum
	if ($average_common_peaks)
	{
		for $i ( 0 .. $#commonpeaks_light )
		{
			$meancommon_mz[$i] = ( mz( $commonpeaks_light[$i] ) + mz( $commonpeaks_heavy[$i] ) ) / 2;
			$verbose && ( $meancommon_delta[$i] = sprintf( "%.4f", ( mz( $commonpeaks_light[$i] ) - mz( $commonpeaks_heavy[$i] ) ) ) );
			$meancommon_int[$i] = ( intensity( $commonpeaks_light[$i] ) + intensity( $commonpeaks_heavy[$i] ) ) / 2;
		}
	}
	$verbose && print "common peaks tolerance=$tolerance, required peakratio = $peakratio\n";
	$verbose && print "mz1\tintensity1\tmz2\tintensity2\tmeanmz\tmeanintensity\tdelta\n";
	my $ncommonpeaks = 0;
	## GENERATE THE FILES FOR THE COMMON AND XLINKER PEAKS
	print COMMONPEAKS "$id\n$mz\n$charge_orginal\n";
	print XLINKERPEAKS "$id\n$mz\n$charge_orginal\n";
	for ( $i = 0 ; $i <= $#commonpeaks_light ; $i += 1 )
	{
		if ( intensity( $commonpeaks_light[$i] ) >= $commonmin )
		{
			$ncommonpeaks++;
			$verbose && print mz( $commonpeaks_light[$i] ), "\t", intensity( $commonpeaks_light[$i] ), "\t\t", mz( $commonpeaks_heavy[$i] ), "\t", intensity( $commonpeaks_heavy[$i] ), "\t", $meancommon_mz[$i], "\t", $meancommon_int[$i], "\t", $meancommon_delta[$i], "\n";
			if ($average_common_peaks)
			{
				print COMMONPEAKS $meancommon_mz[$i], "\t", $meancommon_int[$i], "\t0", "\n";
			} else
			{
				print COMMONPEAKS mz( $commonpeaks_light[$i] ), "\t", intensity( $commonpeaks_light[$i] ), "\t0", "\n";
			}
		}
	}
	if ($verbose)
	{
		print "xlinker peaks\n";
		for ( $i = 0 ; $i <= $#xlinkpeaks_light ; $i += 1 )
		{
			print mz( $xlinkpeaks_light[$i] ), "\t", intensity( $xlinkpeaks_light[$i] ), "\tcharge: ", charge( $xlinkpeaks_light[$i] ), "\t";
			print mz( $xlinkpeaks_heavy[$i] ), "\t", intensity( $xlinkpeaks_heavy[$i] ), "\tcharge: ", charge( $xlinkpeaks_heavy[$i] ), "\t";
			print "\n";
		}
	}
	my $nxlinkpeaks = 0;
	

	
	foreach my $xlinkpeak (@xlinkpeaks_light)
	{
		if ( intensity($xlinkpeak) >= $xlinkmin )
		{
			$nxlinkpeaks++;
			#print XLINKERPEAKS mz($xlinkpeak), "\t", intensity($xlinkpeak), "\t", charge($xlinkpeak), "\n";		
			### Deconvolute if the option is selected
			if ($deconvolutexlinkions)
			{
				my $deconvolutedpeak = deconvolute($xlinkpeak);
				print XLINKERPEAKS mz($deconvolutedpeak) . "\t" . intensity($deconvolutedpeak) . "\t" . charge($deconvolutedpeak) . "\n";
			} else
			{
				print XLINKERPEAKS mz($xlinkpeak) . "\t" . intensity($xlinkpeak) . "\t" . charge($xlinkpeak) . "\n";
			}	
		}
	}
	### PRINT OVERLAY SPECTRUM ################################################
	if ($printspectra)
	{
		require specplot;
		my $specfile = join "", $outfile, "_specplot.png";
		foreach my $peakpair (@PL2_plot)
		{
			$peakpair->[1] *= -1;
		}
		foreach my $peakpair (@commonpeaks_heavy)
		{
			$peakpair->[1] *= -1;
		}
		foreach my $peakpair (@xlinkpeaks_heavy)
		{
			$peakpair->[1] *= -1;
		}
		my $specobj1 = specplot->new();

		#$specobj1->setcolor( "grey", "grey", "green", "red" );
		$specobj1->plotdata( $minionsize, $maxionsize, [ "grey", "grey", "green", "red" ], \@PL1_plot, \@PL2_plot, \@commonpeaks_light, \@xlinkpeaks_light );
		$specobj1->drawlegend( 200, 10, { basename($outfile), "black" } );
		my %colorhash = (
						  "positive: light spectrum peaks" => "grey",
						  "negative: heavy spectrum peaks" => "grey",
						  "matched common-spectrum peaks"  => "green",
						  "matched xlink-spectrum peaks"   => "red",
		);
		$specobj1->drawlegend( 600, 20, \%colorhash );
		$specobj1->printimage($specfile);
	}
	if ($weboutput)
	{
		print '<br><b>', basename($input1), ",", basename($input2), ': sorting peaks into common peaks and x-linker peaks:</b><br>';
		print 'number of peaks spectrum1  > ', $threshold, ': ',
		  scalar(@peaklist_light), '<br>number of peaks spectrum2 > ',
		  $threshold, ': ', scalar(@peaklist_heavy),
		  '<br>number of  common-peaks > ', sprintf( "%.2f", $commonmin ), ' ',
		  $ncommonpeaks, '<br>number of xlink-peaks > ',
		  sprintf( "%.2f", $xlinkmin ), ' ', $nxlinkpeaks, '<br>';
	} else
	{
		print "#peaks1 > ", $threshold, ": ", scalar(@peaklist_light), "|  #peaks2 > ", $threshold, ": ", scalar(@peaklist_heavy), "| #commonpeaks > ", sprintf( "%.2f", $commonmin ), "% ", $ncommonpeaks, "| #xlinkpeaks > ", sprintf( "%.2f", $xlinkmin ), "% ", $nxlinkpeaks, "<br>\n";
	}
	close(COMMONPEAKS);
	close(XLINKERPEAKS);

	#exit;
	return 1;
}

sub difference
{
	my $list1 = shift;
	my $list2 = shift;
	my ( @a, @b );
	foreach (@$list1)
	{
		push @a, mz($_);
	}
	foreach (@$list2)
	{
		push @b, mz($_);
	}
	my @union = my @isect = my @diff = ();
	my %union = my %isect = ();
	my %count = ();
	my ($e);
	foreach $e ( @a, @b ) { $count{$e}++ }
	foreach $e ( keys %count )
	{
		push( @union, $e );
		push @{ $count{$e} == 2 ? \@isect : \@diff }, $e;
	}
	return \@diff;
}

sub isect
{
	my $list1 = shift;
	my $list2 = shift;
	my @union = my @isect = my @diff = ();
	my %union = my %isect = ();
	my %count = ();
	my ($e);
	foreach $e ( @$list1, @$list2 ) { $union{ mz($e) }++ && $isect{ mz($e) }++ }
	@isect = keys %isect;
	return ( \@isect );
}

sub mz
{
	my $pairs = shift;
	return $pairs->[0];
}

sub intensity
{
	my $pairs = shift;
	return $pairs->[1];
}

sub printpeak
{
	my $pairs = shift;
	print $pairs->[0] . " " . $pairs->[1] . "\n";
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

sub charge
{
	my $pairs = shift;
	return $pairs->[2];
}

sub setcharge
{
	my $pairs  = shift;
	my $charge = shift;
	$pairs->[2] = $charge;
}

sub matchpeaks
{

	# Matchfunction for 2 peaklists: prerequesite is that the 2 peaklists are sorted ascending by mz
	# Changed on 03 / 09 / 2009 by TW, OR
	my $ions         = shift;
	my $peaks        = shift;
	my $ms2tolerance = shift;
	my $lower_ratio  = shift;
	my $upper_ratio  = shift;
	my $chargestate  = shift;
	my $i;
	my $j         = 0;
	my $nhits     = 0;
	my $lastindex = 0;
	my @matched   = ();
### CP TOLERANCES for PPM MATCHING
	my $tolerance  = $ms2tolerance;
	my $ppmmeasure = 0;

	if ( $tolerancemeasure =~ /^ppm/i )
	{
		$ppmmeasure = 1;

		#print "tolerancemeasure: ppm tolerance $ms2tolerance\n";
	}
	for $i ( 0 .. $#$ions )
	{
		if ( $lastindex == 0 )
		{
			$j = 0;
		} else
		{
			$j = $lastindex;
		}
		if ($ppmmeasure)
		{

			#print $ms2tolerance;
			$tolerance = $ms2tolerance * 1e-6 * mz( $ions->[$i] );    #ppm to amu measure
		}

		#print "Allowed tolerance for ion with ".mz( $ions->[$i])." m/z: $tolerance m/z\n";
		$verbose && print "ion:$i lastindex:$lastindex ";
		$verbose && print "ion: ", mz( $ions->[$i] ) . " ";
		$verbose && print "peak: $j: ", mz( $peaks->[$j] ), " delta: ", mz( $peaks->[$j] ) - mz( $ions->[$i] ), "\n";

		# 1. IF THE PEAK (LIST 2) is smaller than the ION plus tolerance (-->goes in the until loop if cond is false)
		#    THEN THE PEAK HAS A CHANCE TO MATCH --> go in the loop
		#    IF PEAK is is bigger than the ION, then no matching is possible, goes not in the loop and
		#    takes the next ION in the list and test the next ION
		until ( ( mz( $peaks->[$j] ) > ( ( mz( $ions->[$i] ) + $tolerance ) ) ) || $lastindex > $#$peaks )
		{
			$verbose && print "in until ion $i: ", mz( $ions->[$i] ) . " peak $j: ", mz( $peaks->[$j] ) . " delta: ", mz( $peaks->[$j] ) - mz( $ions->[$i] ), "\n";

			# 2. If the PEAK is smaller or equal than the ION-tolerance
			# then this peak will never match against a larger ION in the next round
			# therefore lastindex is set to this peaknumber (j)
			if ( ( mz( $peaks->[$j] ) <= ( mz( $ions->[$i] ) - $tolerance ) ) )
			{
				$lastindex = $j;
				$verbose && print "delta: ", mz( $peaks->[$j] ) - mz( $ions->[$i] ), " lastindex =j: " . $lastindex . "\n";
			}

			# 4. Test if the ION matches the PEAK within the tolerance
			#	 and test if the match is within the required intensity ratio
			#	 push matched ion into the @matched array
			if (
				 ( abs( mz( $ions->[$i] ) - mz( $peaks->[$j] ) ) <= $tolerance )
				 && ( intensity( $ions->[$i] ) / intensity( $peaks->[$j] ) >= $lower_ratio
					  && ( intensity( $ions->[$i] ) / intensity( $peaks->[$j] ) ) <= $upper_ratio )
			  )
			{
				$verbose && print "matched peak: " . mz( $ions->[$i] ), " ", mz( $peaks->[$j] ), " $chargestate\n";
				$ions->[$i]->[2] = $chargestate;
				$nhits++;
				push @matched, [ mz( $ions->[$i] ), intensity( $ions->[$i] ), $chargestate ];
			}

			# 5. Cout up the counter j to test the next PEAK with the same ION
			$j++;
		}
		$verbose && print "PEAK to BIG to MATCH->go for next ION:peak: $j: ", mz( $peaks->[$j] ), " ion:$i ", mz( $ions->[$i] ) . " delta: ", mz( $peaks->[$j] ) - mz( $ions->[$i] ), "\n";
	}
	return @matched;
}

sub matchpeaks_oldbug
{
	my $ions        = shift;
	my $peaks       = shift;
	my $tolerance   = shift;
	my $lower_ratio = shift;
	my $upper_ratio = shift;
	my $chargestate = shift;
	my $i;
	my $j         = 0;
	my $nhits     = 0;
	my $lastindex = 0;
	my @matched   = ();

	# Matchfunction by TW: prerequesite is that the 2 peaklists are sorted ascending by mz
	for $i ( 0 .. $#$ions )
	{
		$j = $lastindex;

		#print mz($ions->[$i])."\n";
		#print mz($peaks->[$i])."\n";
		# if the condition
		until ( ( ( mz( $peaks->[$j] ) - mz( $ions->[$i] ) ) >= $tolerance ) || $lastindex > $#$peaks )
		{
			if (
				 ( abs( mz( $ions->[$i] ) - mz( $peaks->[$j] ) ) <= $tolerance )
				 && ( intensity( $ions->[$i] ) / intensity( $peaks->[$j] ) >= $lower_ratio
					  && ( intensity( $ions->[$i] ) / intensity( $peaks->[$j] ) ) <= $upper_ratio )
			  )
			{
				$verbose && print mz( $ions->[$i] ), " ", mz( $peaks->[$j] ), " $chargestate\n";
				$ions->[$i]->[2] = $chargestate;
				$nhits++;
				push @matched, [ mz( $ions->[$i] ), intensity( $ions->[$i] ), $chargestate ];
			}
			if ( ( mz( $peaks->[$j] ) - mz( $ions->[$i] ) ) < $tolerance )
			{
				$lastindex++;
			}
			$j++;
		}
	}
	return @matched;
}

sub readtables
{
	my $xquestdef = shift;
	open DEF, "<$xquestdef"
	  or die "cannot open xquest definition file $xquestdef $!";

	#read in definitions
	#read in definitions
	my (%DIGEST);
	my $configdef     = undef;
	my $enzymedef     = undef;
	my $modifications = undef;
	while (<DEF>)
	{
		chomp;
		my $line = $_;
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
			my @splitted = split( " ", $line );
			$DIGEST{ $splitted[0] } = $splitted[1];

			#print $splitted[0], "\t", $DIGEST{ $splitted[0] }, "\n";
		}
	}
	return ( \%DIGEST );
}

sub print_affi_and_changelog
{
	my $scriptinfo = shift;
	print "version " . $scriptinfo->{'version'} . " written by ";
	print $scriptinfo->{'author'} . "\n";
	print "Affiliation: " . $scriptinfo->{'affi'} . "\n";
	print "In case of troubles mailto: ", $scriptinfo->{'mailto'} . "\n";
	print "Changelog:\n";
	foreach my $key ( sort { $a <=> $b } keys %{ $scriptinfo->{'clog'} } )
	{
		print "Version " . $key . ": " . $scriptinfo->{'clog'}->{$key};
	}
	return;
}

sub usage()
{
print "
	SOFTWARE: ", basename($0), " version $version
	
	AUTHOR: Thomas Walzthoeni modified from orginal script by Oliver Rinner.

	INFORMATION: Program to compare light and heavy MS/MS scan pairs.
 
 	USAGE: ", basename($0), " -Option [Parameter]

	REQUIRED OPTIONS: 
	-match [] isotopepair file; format: xmm.pl output file
	-dir [] directory where files are stored
	-def [] xquest.def filename  
	-resultdir [] directory for search results
 	
 	OTHER OPTIONS [defaults]:
 	-genxml [] mzXML filename where spectra are extracted from, creates a spec.xml file where the compared spectra are stored.

	EXAMPLE:
	$0 -match FN-XL_matched.txt -dir FN-XL_matcheddir -def xquest.def -genxml FN-XL.c.mzXML -resultdir FN-XL_matched

";
	exit;
}
