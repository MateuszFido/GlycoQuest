#!/usr/bin/perl
use strict;

#---------------------------------------------------------------------------
# runXquest.pl
## A wrapper to run the xQuest search pipeline. It generates (creates *.sh or *.bat files) or executes (optional) the sequence of commands to
## to run a xquest search, or multiple searches.
## Genreal workflow of this program
## A. Extract spectra with MzXml2Search (tool from the TPP)
## B. Search for isotopic scan pairs (xmm.pl).
## C. Comparision of light and heavy MS/MS scans by compare_peaks3.pl
## D. xQuest Search
## Not necessarily the whole sequence of commands is generated or executed:
## If a spectrum directory is present, step A is skipped.
## If *matched.txt file is present, steps A and B are skipped. (see -help for details)
#
# Execute runXquest.pl -help to display information and usage options.
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
use File::Basename;
use Getopt::Long;
use File::Spec;
use Cwd;
use Data::Dumper;
use File::Copy;
##########################################################
# Include modules dir as lib that is relative to the Script path
##########################################################
use FindBin;
use lib "$FindBin::Bin/../modules";
##########################################################
use Environment;
my ( $getdef, $sh, $xmmonly, $runonprottools, $verbose, $help, $xmlfiles, $copydb, $matchlist, $nocompare, $timetolive, $memory, $dev, $listfile, $cponly, $pseudosh, $xmmextracmds, $cp_extracmds, $xq_extracmds, $ssh, $xmlmode, $cpforce, $noxmm, $xmmforce, $brutus );
## Here the path to a Perl installation can be defined if perl is not defined in the PATH variable.
my $perlpath = "";    ## eg: add "/path/to/perl/perl-5 " ## this is prepeded to all commands if used
## Standard Parameters
my $xquestdef                   = "xquest.def";                     ## default file name of xquest definition name
my $xmmdef                      = "xmm.def";                        ## xmm definition file
my $MASTERMAP                   = "./MASTER_RUN/MASTER_RUN.txt";    ## default file name of MASTER_MAP file
my $force                       = 0;
my $resultdir                   = "XQUESTresults";                  ## default resultdirectory
my $converted_spectra_directory = "XQUESTdir";                      ## default converted spectra directory
my $email                       = 'walzthoeni@imsb.biol.ethz.ch';
my $parallel                    = 0;
my $ms2extract_cmds             = "-T15000";
$xq_extracmds = "-nidx";
my $scriptinfo = {};
$scriptinfo->{'version'} = "1.2";
$scriptinfo->{'author'}  = "Thomas Walzthoeni";
#### Create the Options Hash and initialize with the standard values
my %PARAMS = (                                                      ## the options that are passed as arguments ovverride the standard values if they are set
	## common settings
	'results'    => $resultdir,                                     ## Result directory name
	'xmlfiles'   => $xmlfiles,                                      ## mzXml filenames can be passed comma separeted
	'force'      => $force,
	'prottools'  => $runonprottools,
	'list'       => $listfile,
	'parallel'   => $parallel,
	'timetolive' => $timetolive,
	'memory'     => $memory,
	'dev'        => $dev,
	'ssh'        => $ssh,
	'verbose'    => $verbose,
	'help'       => $help,
	## Extraction of dta spectra from mzxml files
	## xmm parameters
	'xmm_master'    => $MASTERMAP,      ## Master Map file, where xmm will search isotopic pairs
	'xmmdef'        => $xmmdef,         ## xmm definition file
	'xmm_extracmds' => $xmmextracmds,
	'xmmonly'       => $xmmonly,        ## option to run only xmm.pl to for pairfinding
	'pseudoSH'      => $pseudosh,
	'xmmforce'      => $xmmforce,
	## compare_peaks2.pl parameters
	'matchlist'       => $matchlist,
	'cponly'          => $cponly,
	'dir'             => $converted_spectra_directory,    ## Directory name for converted spectra, if not specified the default name is $mastermapbasename+dir
	'cp_extracmds'    => $cp_extracmds,
	'cpforce'         => $cpforce,
	'xq_extracmds'    => $xq_extracmds,
	'ms2extract_cmds' => $ms2extract_cmds,
	'nocompare'       => $nocompare,
	## xquest parameters
	'xquestdef' => $xquestdef,
	'testcmd'   => undef,
	'xmlmode'   => $xmlmode,
	'noxmm'     => $noxmm,
	'euler'    => $brutus,
	'sh'        => $sh,
	'getdef'    => $getdef,
);
## HERE ALL PARAMS THAT COME FROM THE CMD LINE MUST BE DEFINED
GetOptions(
			\%PARAMS,         'ssh',            'xquestdef=s', 'matchlist=s', 'parallel=i', 'verbose',     'list=s',    'pseudoSH', 'xmm_master=s', 'xmm_extracmds=s',
			'cp_extracmds=s', 'xq_extracmds=s', 'force',       'dev=s',       'help',       'xquestdef=s', 'xmmdef=s',  'testcmd',  'xmmonly',      'cponly',
			'timetolive=i',   'xmlmode',        'cpforce',     'noxmm',       'xmmforce',   'brutus',      'nocompare', 'sh',       'getdef'
);
my $PARAMS = \%PARAMS;
unless ( $PARAMS->{'matchlist'} or $PARAMS{'list'} or $PARAMS{'getdef'} )
{
### Then check if there is a Master Map file
	my ( $status, $msg ) = check_file( $PARAMS->{'xmm_master'} );
	unless ( $status == 1 )
	{
		print "$msg";
		&usage();
	}
}
&usage() if ( $PARAMS->{'help'} );
my $OS = $^O;
my $ftype;
### Check which OS is used:
if ( $OS =~ m/MSWin/ )
{
	$ftype = "bat";
}
if ( $OS eq "linux" )
{
	$ftype = "sh";
}
## Set the rootdirectory where we are now
$PARAMS->{'rootdir'} = getcwd;
## Set the username
my $user_name = $ENV{'LOGNAME'};
chomp($user_name);
$PARAMS->{'username'} = $user_name;
##
my $env    = Environment->new;
my $server = $env->get_env;
print "Current Server: $server\n";
### Set the xquest paths
my ( $cp_path, $bin );
my $basedir;

#---------------------------------------------------------------------------
#  Set the Working directories of xQuest // Defined in Environment object
#---------------------------------------------------------------------------
if ( $PARAMS->{'dev'} )
{
	### Concat if a specific version is used
	$basedir = File::Spec->catfile( $env->get_path('xquest_dev'), $PARAMS->{'dev'} );
	$cp_path = File::Spec->catfile($basedir);
	$bin     = File::Spec->catfile( $basedir, "bin" );
} else
{
	### Use the currently defined stable version see ENV Obj
	$basedir = File::Spec->catfile( $env->get_path('xquest_stable') );
	$cp_path = File::Spec->catfile($basedir);
	$bin     = File::Spec->catfile( $basedir, "bin" );
}

#---------------------------------------------------------------------------
#  Searchlog
#---------------------------------------------------------------------------
use SimpleAuth;
my $authobj = SimpleAuth->new($basedir);
my $usagelogfile = File::Spec->catfile( $basedir, "logs", "usage.log" );
$authobj->set_logfile($usagelogfile);
my $user = getlogin();
my $msg  = "User: $user, xQuest Search Submitted from $PARAMS->{'rootdir'}.\n";
$authobj->writelog($msg);
$PARAMS->{'bin'}        = $bin;
$PARAMS->{'xquestpath'} = $cp_path;
$PARAMS->{'server'}     = $env->get_env;

#---------------------------------------------------------------------------
#  Check if the paths are set
#---------------------------------------------------------------------------
unless ( $PARAMS->{'bin'} and $PARAMS->{'xquestpath'} and $PARAMS->{'server'} )
{
	print "Error: Cannot find xquest path or server name!";
	&usage();
}

#---------------------------------------------------------------------------
#  Verbose params
#---------------------------------------------------------------------------
if ( $PARAMS->{'verbose'} )
{
	print_params($PARAMS);
}
## checks if a directory list is provided
my @directories = ();

#---------------------------------------------------------------------------
#  get deffiles and exit
#---------------------------------------------------------------------------
if ( $PARAMS->{'getdef'} )
{
	$basedir = File::Spec->catfile( $env->get_path('xquest_stable') );
	my $tmppathxquest0 = File::Spec->catfile( $basedir, 'deffiles', 'xQuest', 'xquest.def' );
	my $tmppathxmm0    = File::Spec->catfile( $basedir, 'deffiles', 'xmm',    'xmm.def' );
	if ( -e $tmppathxquest0 )
	{
		print "Copy xquest.def template to current folder\n";
		copy( $tmppathxquest0, $PARAMS->{'rootdir'} );
	} else
	{
		warn "Cannot find xquest.def template $tmppathxquest0\n";
	}
	if ( -e $tmppathxmm0 )
	{
		print "Copy xmm.def template to current folder\n";
		copy( $tmppathxmm0, $PARAMS->{'rootdir'} );
	} else
	{
		warn "Cannot find xmm.def template $tmppathxmm0\n";
	}
	exit;
}

#---------------------------------------------------------------------------
#  Check if deffiles exist
#---------------------------------------------------------------------------
if ( $PARAMS{'list'} )
{
	$listfile = $PARAMS{'list'};
	## open the input file
	open INPUT, "<$listfile" or die "Cannot open file $listfile";
	## read the lines
	while ( my $line = <INPUT> )
	{
		chomp $line;
		push @directories, $line;
	}
	close(INPUT);
	if ( $PARAMS{'verbose'} )
	{
		print "Directory list was specified, file contains the following lines:\n";
		foreach my $dir (@directories)
		{
			print "Checking $dir for definition files.\n";
		}
	}
	## check if the deffiles are in there
	foreach my $dir (@directories)
	{
		my $tmppathxquest = File::Spec->catfile( $PARAMS->{'rootdir'}, $dir, $PARAMS->{'xquestdef'} );
		my $tmppathxmm    = File::Spec->catfile( $PARAMS->{'rootdir'}, $dir, $PARAMS->{'xmmdef'} );
		unless ( -e $tmppathxquest )
		{
			print "Error: Cannot find xquest definition file $tmppathxquest $!<br>";
			&usage();
		}
		unless ( -e $tmppathxmm )
		{
			print "Error: Cannot find xmm definition file $tmppathxmm $!<br>";
			&usage();
		}
	}
} else
{
## if no list file was provided
## check if xquest definition file is readable
	if ( $PARAMS{'verbose'} )
	{
		print "Checking for definition files in working directory.\n";
	}
	unless ( -e $PARAMS->{'xquestdef'} )
	{
		die "Cannot find xquest definition file $PARAMS->{'xquestdef'} $!<br>";
	}
## check if xmm definition file is readable
	unless ( -e $PARAMS->{'xmmdef'} )
	{
		die "Cannot find xmm definition file  $PARAMS->{'xmmdef'} $!<br>";
	}
}

#---------------------------------------------------------------------------
#  Generate a new run.sh file in every subfolder (if a listfile is provided),
#  oherwise generate a run.sh file in the cwd
#---------------------------------------------------------------------------
my $rootdir = $PARAMS->{'rootdir'};
if ( $PARAMS->{'list'} )
{
	## delete all run.sh files
	foreach my $dir (@directories)
	{
		print "\n#### Delete existing sh files in subfolder $dir ####\n";
		my $shfilewc = File::Spec->catfile( $rootdir, $dir, "runxq*.$ftype" );
		my @shfiles = glob($shfilewc);
		## Deletes all sh files
		foreach my $file (@shfiles)
		{

			#print "Delete sh file $file, will create new\n";
			unlink($file);
		}
	}
} else
{
	print "\n#### Delete existing sh files in dir $rootdir ####\n";
	my $shfilewc = File::Spec->catfile( $rootdir, "runxq*.$ftype" );
	my @shfiles = glob($shfilewc);
	## Deletes all png files
	foreach my $file (@shfiles)
	{

		#print "Delete sh file $file, will create new\n";
		unlink($file);
	}
}

#---------------------------------------------------------------------------
#  Extraction of MS2 scans
#---------------------------------------------------------------------------
print "\n#### Extration of MS2 spectra from mzXML files #### \n";
my ( $extractms2, $extractms2path ) = _extract_ms2($PARAMS);
for ( my $i = 0 ; $i < scalar(@$extractms2) ; $i++ )
{
	print "Change directory to: $extractms2path->[$i]\n";
	print "Command: $extractms2->[$i]\n";
	unless ( $PARAMS->{'testcmd'} )
	{
		### Always execute directly
		chdir( $extractms2path->[$i] );
		system( $extractms2->[$i] );
	}
}

#---------------------------------------------------------------------------
#  xmm.pl, extract master map
#---------------------------------------------------------------------------
##########################################################
# MASTER_MAP.txt --> xmm.pl --> Matchlist (*.matched.txt)#
# matchlist option -iso matchlistname					 #
# xmm.def is the param file								 #
# xmm Searches for isotopic pairs in a feature list    	 #
##########################################################
print "\n#### Search for isotopic pairs with xmm.pl #### \n";
my ( $xmmcmds, $xmmlocations ) = _get_xmm_cmds( \%PARAMS );
( $PARAMS->{'noxmm'} ) && print "-noxmm selected: will skip xmm.pl\n";
for ( my $i = 0 ; $i < scalar(@$xmmcmds) ; $i++ )
{
	print "Change directory to: $xmmlocations->[$i]\n";
	print "Command: $xmmcmds->[$i]\n";
	unless ( $PARAMS->{'testcmd'} or $PARAMS->{'noxmm'} )
	{
		### Always execute xmm directly
		chdir( $xmmlocations->[$i] );
		system( $xmmcmds->[$i] );
	}
}
if ( $PARAMS->{'xmmonly'} )
{
	print "-xmmonly option selected, will exit now, bye!\n";
	exit;
}

#---------------------------------------------------------------------------
#  PARALELL processing, spliting of the matchlists
#---------------------------------------------------------------------------
if ( $PARAMS->{'parallel'} )
{
	print "\n#### Parallel processing #### \n";
	my $matchfilesarrayref = _parallel_processing($PARAMS);
}

#---------------------------------------------------------------------------
#  Compare Peaks Section
#---------------------------------------------------------------------------
## generated the inclusion list for the
## xquest search: e.g  M09-12171_matched_isotopepairs.txt
## creates also the xquest result dir
print "\n#### Prepare compare_peaks3.pl commands ####\n";
#### If Parallel processing is selected then for every matchlist a cpm must be performed
my ( $cp_cmds, $cp_locationarray ) = _get_cp_cmds( \%PARAMS );
my @cp_unique_array;
my @waitcmd;
for ( my $i = 0 ; $i < scalar(@$cp_cmds) ; $i++ )
{
	### check if we want to submit the job to brutus
	if ( $PARAMS->{'brutus'} )
	{
		my $unique_id = unique_id();
		my $waitcmd   = "-w " . '"ended (' . $unique_id . ')"';
		push @waitcmd, $waitcmd;

		#my $cmd = "bsub -R\"select[imsb]\" -J $unique_id $cp_cmds->[$i]";
		my $cmd = "bsub -J $unique_id $cp_cmds->[$i]";
		print "Compare_peaks command: $cmd\n";
		print "Change directory to: $cp_locationarray->[$i]\n";
		unless ( $PARAMS->{'testcmd'} || $PARAMS->{'nocompare'} )
		{
			my $runxqfile = File::Spec->catfile( $cp_locationarray->[$i], "runxq$i.$ftype" );
			open( MYFILE, '>>', $runxqfile );
			print MYFILE "cd $cp_locationarray->[$i]\n";
			print MYFILE $perlpath . $cmd . ">cp$i.log\n";
			close(MYFILE);
			unless ( $PARAMS->{'sh'} )
			{
				chdir( $cp_locationarray->[$i] );
				system( $cmd. ";" );
			}
		}
	} else
	{
		## run on the local machine
		print "Compare_peaks command: ";
		print "Change directory to: $cp_locationarray->[$i]\n";
		print "Command: $cp_cmds->[$i]\n";
		unless ( $PARAMS->{'testcmd'} || $PARAMS->{'nocompare'} )
		{
			my $runxqfile = File::Spec->catfile( $cp_locationarray->[$i], "runxq$i.$ftype" );
			open( MYFILE, '>>', $runxqfile );
			print MYFILE "cd $cp_locationarray->[$i]\n";
			print MYFILE $perlpath . $cp_cmds->[$i] . ">cp$i.log\n";
			close(MYFILE);
			unless ( $PARAMS->{'sh'} )
			{
				chdir( $cp_locationarray->[$i] );
				system( $cp_cmds->[$i] );
			}
		}
	}
}
if ( $PARAMS->{'nocompare'} )
{
	print "-nocompare option selected,cmds are not executed, bye!\n";
	@waitcmd = ();
}
if ( $PARAMS->{'cponly'} )
{
	print "-cponly option selected, will exit now, bye!\n";
	exit;
}

#---------------------------------------------------------------------------
#  xQuest Section
#---------------------------------------------------------------------------
print "\n#### Prepare xQuest commands ####\n";
my ( $xq_cmds, $xq_locationarray ) = _get_xquest_cmds( \%PARAMS );

#print Dumper($xq_cmds);
#print Dumper($xq_locationarray);
for ( my $i = 0 ; $i < scalar(@$xq_cmds) ; $i++ )
{
	if ( $PARAMS->{'brutus'} )
	{
		my $bsuboptions;
		## check if a time flag is set
		if ( ( $PARAMS->{'timetolive'} ) )
		{
			$bsuboptions .= "-W " . $PARAMS->{'timetolive'} . ":00";
		} else
		{
			$bsuboptions .= "-W 4:00";
		}
		#my $runXquest = "bsub -e errorfile $waitcmd[$i] $bsuboptions -R \"select[model==Opteron8380]\" $xq_cmds->[$i]";
		my $runXquest = "bsub -e errorfile $waitcmd[$i] $bsuboptions $xq_cmds->[$i]";
		print "Change directory to: $xq_locationarray->[$i]\n";
		print "xQuest command: $runXquest\n";
		unless ( $PARAMS->{'testcmd'} )
		{
			my $runxqfile = File::Spec->catfile( $xq_locationarray->[$i], "runxq$i.$ftype" );
			open( MYFILE, '>>', $runxqfile );
			print MYFILE "cd $xq_locationarray->[$i]\n";
			print MYFILE $perlpath . $xq_cmds->[$i] . ">xq$i.log\n";
			close(MYFILE);
			unless ( $PARAMS->{'sh'} )
			{
				chdir( $xq_locationarray->[$i] );
				system( $runXquest. ";" );
			}
		}
	} else
	{
		print "Run xQuest on localmachine:";
		print "Change directory to: $xq_locationarray->[$i]\n";
		print "Command: $xq_cmds->[$i]\n";
		unless ( $PARAMS->{'testcmd'} )
		{
			my $runxqfile = File::Spec->catfile( $xq_locationarray->[$i], "runxq$i.$ftype" );
			open( MYFILE, '>>', $runxqfile );
			print MYFILE "cd $xq_locationarray->[$i]\n";
			print MYFILE $perlpath . $xq_cmds->[$i] . ">xq$i.log\n";
			close(MYFILE);
			unless ( $PARAMS->{'sh'} )
			{
				chdir( $xq_locationarray->[$i] );
				system( $xq_cmds->[$i] );
			}
		}
	}
}

# only extracts spectra if no spectra folder is present or the -force option is used
# if a list of directories is used go into every directory
sub _extract_ms2
{
	my $PARAMS          = shift;
	my $bin             = $PARAMS->{'xquestpath'};
	my $verbose         = $PARAMS->{'verbose'};
	my $rootdir         = $PARAMS->{'rootdir'};
	my $ms2extract_cmds = $PARAMS->{'ms2extract_cmds'};
	my $force           = $PARAMS->{'force'};
	my @mzXML_files;
	my @extractms2;
	my @extractms2path;
	chdir($rootdir);

	if ( $PARAMS->{'xmlmode'} )
	{
		print "xmlmode selected, -won't extract ms/ms scans\n";
		return \@extractms2, \@extractms2path;
	}
	if ( $PARAMS->{'list'} )
	{
		foreach my $dir (@directories)
		{
			print "\n#### Extraction of dta files ####\n";
			my $tmppath = File::Spec->catfile( $rootdir, $dir );

			#$verbose && print "Change to directory $tmppath\n";
			@mzXML_files = glob("$tmppath/*.mzXML");

			#print Dumper(@mzXML_files);
			foreach my $mzxml (@mzXML_files)
			{

				#$verbose && print "Found mzXML file: $mzxml\n";
				my $specdir = basename( $mzxml, ".mzXML" );

				#$verbose && print "Spectrum directory will be: $specdir\n";
				my $specdirpath = File::Spec->catfile( $tmppath, $specdir );
				## Spectra are only extracted if there is not already a spectrum dirctory
				if ( -e $specdirpath && !$force )
				{
					print "Directory $specdir already exists, spectra are not extracted. (use the -force option to extract)\n";
					next;
				}
				my $mzxml2search = "MzXML2Search -dta $ms2extract_cmds $mzxml";
				print "Extract MS2 scans:-->$mzxml2search\n";
				push @extractms2,     $mzxml2search;
				push @extractms2path, $tmppath;
			}
		}
	} else
	{
## if no listfile was provided --> then the user should be in the targetfolder
		#print "\n#### Extraction of dta files ####\n";
		my $tmppath = File::Spec->catfile($rootdir);
		@mzXML_files = glob("*.mzXML");
		unless (@mzXML_files)
		{
			print "Warning: No mzXML files found in this folder, wont extract any scans.\n";
		}
		foreach my $mzxml (@mzXML_files)
		{
			my $specdir = basename( $mzxml, ".mzXML" );

			#$verbose && print "Found mzXML file: $mzxml\n";
			## Spectra are only extracted if there is nprotot already a spectrum dirctory
			if ( -e $specdir && !$force )
			{
				print "Directory $specdir already exists, spectra are not extracted. (use the -force option to extract)\n";
				next;
			}
			my $mzxml2search = "MzXML2Search -dta -T20000 $mzxml";
			push @extractms2,     $mzxml2search;
			push @extractms2path, $tmppath;
		}
	}
	return \@extractms2, \@extractms2path;
}

sub _get_xmm_cmds
{
	my $PARAMS    = shift;
	my $bin       = $PARAMS->{'bin'};
	my $listfile  = $PARAMS->{'list'};
	my $rootdir   = $PARAMS->{'rootdir'};
	my $verbose   = $PARAMS->{'verbose'};
	my $MASTERMAP = $PARAMS->{'xmm_master'};
	my $xmmdef    = $PARAMS->{'xmmdef'};
	my @xmmcmds;
	my @xmmlocations;
	my $tmppath;
	my $matchlist;
	my $basenamematchlist;

	#my $xmm_extracmds;
	if ( $PARAMS->{'matchlist'} )
	{
		print "Matchlist is specified, will skip xmm.pl.\n";
		return \@xmmcmds, \@xmmlocations;
	}
	my $xxquest;
	chdir($rootdir);    #change to the rootdir to find the list file
### Generate xmm cmd for all folders
	if ($listfile)
	{
### Get the directories
		my @directories = read_file($listfile);
		foreach my $dir (@directories)
		{
			my $xmm_extracmds;
			if ( $PARAMS->{'pseudoSH'} )
			{
				$xmm_extracmds .= " -pseudoSH";
			}
			$tmppath = File::Spec->catfile( $rootdir, $dir );
			$matchlist = $dir . "_matched.txt";
			## change directory
			chdir($tmppath);
			## Get all mzxmlfiles
			my @mzXML_files = glob("*.mzXML");
			### check if matchlist exists
			if ( -e $matchlist && !$PARAMS->{'xmmforce'} )
			{
				print "Matchlist $matchlist already exists, will skip xmm.pl. Use -xmmforce to re-extract\n";
				return \@xmmcmds, \@xmmlocations;
			}
			foreach my $mzxml (@mzXML_files)
			{
				$verbose && print "Searching for isotopic pairs in $MASTERMAP, parameters are defined in $xmmdef...\n";
				$xmm_extracmds .= " -mz \"$mzxml\"";
				$xxquest = "$bin\/xmm.pl -master $MASTERMAP -iso $matchlist -def $xmmdef $xmm_extracmds";
				push @xmmcmds,      $xxquest;
				push @xmmlocations, $tmppath;
			}
		}
	} else
	{
		### We are in the working directory
		$tmppath   = File::Spec->catfile($rootdir);
		$matchlist = basename($rootdir) . "_matched.txt";
		### check if matchlist exists
		if ( -e $matchlist && !$PARAMS->{'xmmforce'} )
		{
			print "Matchlist $matchlist already exists, will skip xmm.pl. Use -xmmforce to re-extract\n";
			return \@xmmcmds, \@xmmlocations;
		}
		my $xmm_extracmds;
		if ( $PARAMS->{'pseudoSH'} )
		{
			$xmm_extracmds .= " -pseudoSH";
		}
		my @mzXML_files = glob("*.mzXML");
		foreach my $mzxml (@mzXML_files)
		{
			$verbose && print "Searching for isotopic pairs in $MASTERMAP, parameters are defined in $xmmdef...\n";
			$xmm_extracmds .= " -mz \"$mzxml\"";
			$xxquest = "$bin\/xmm.pl -master $MASTERMAP -iso $matchlist -def $xmmdef $xmm_extracmds";
			push @xmmcmds,      $xxquest;
			push @xmmlocations, $tmppath;
		}
	}
	return \@xmmcmds, \@xmmlocations;
}

sub _parallel_processing
{
	my $PARAMS   = shift;
	my $listfile = $PARAMS->{'list'};
	my $verbose  = $PARAMS->{'verbose'};
	my $rootdir  = $PARAMS->{'rootdir'};
	my $matchfilesarrayref;
	chdir($rootdir);    ##change to the rootdir to find the list file
### Generate xmm cmd for all folders
	if ($listfile)
	{
		### Get the directories
		my @directories = read_file($listfile);
		foreach my $dir (@directories)
		{
			if( $dir ne "" )
			{
				my $tmppath = File::Spec->catfile( $rootdir, $dir );
				chdir($tmppath);
				my $basenamematchlist = $dir . "_matched";
				my $mlfullpath = File::Spec->catfile( $tmppath, $basenamematchlist . ".txt" );
				($verbose) && print "Parallelization: Split Match file ($basenamematchlist) into $PARAMS->{'parallel'} lines per file\n";
				$matchfilesarrayref = _split_matchlist( $basenamematchlist, $mlfullpath, $PARAMS );
			}
		}
	} else
	{
		my $tmppath = File::Spec->catfile($rootdir);
		my $basenamematchlist;
		if ( $PARAMS->{'matchlist'} )
		{
			## Then the basename is defined.
			$basenamematchlist = $PARAMS->{'matchlist'};
			#$basenamematchlist =~ s/\.\w+//;    ## chop the .txt #Doesn't work on paths with . in the name
			my ($basenamematchlist_file, $basenamematchlist_dir) = File::Basename::fileparse($basenamematchlist);
			$basenamematchlist = $basenamematchlist_dir . $basenamematchlist_file;
			
		} else
		{
			$basenamematchlist = basename($rootdir) . "_matched";    ## try if it is the same as the dirname
		}
		my $mlfullpath = File::Spec->catfile( $tmppath, $basenamematchlist . ".txt" );
		($verbose) && print "Parallelization: Split Match file ($basenamematchlist) into $PARAMS->{'parallel'} lines per file\n";
		$matchfilesarrayref = _split_matchlist( $basenamematchlist, $mlfullpath, $PARAMS );
	}
	return $matchfilesarrayref;
}

sub _get_cp_cmds
{
	my $PARAMS    = shift;
	my $bin       = $PARAMS->{'bin'};
	my $listfile  = $PARAMS->{'list'};
	my $rootdir   = $PARAMS->{'rootdir'};
	my $verbose   = $PARAMS->{'verbose'};
	my $xquestdef = $PARAMS->{'xquestdef'};
	my @cp_cmds;
	my @cp_locationarray;
	my $comparecmd;
	my $tmppath;
	chdir($rootdir);    #change to the rootdir to find the list file

	# Is cpforce selected?
	if ( $PARAMS->{'cpforce'} )
	{
		$PARAMS->{'cp_extracmds'} .= " -cpforce";
	}
### Generate xmm cmd for all folders
	if ($listfile)
	{
### Get the directories
		my @directories = read_file($listfile);
		foreach my $dir (@directories)
		{
			$tmppath = File::Spec->catfile( $rootdir, $dir );
			my $cp_extracmds = $PARAMS->{'cp_extracmds'};
			if ( $PARAMS->{'xmlmode'} )
			{
				## Get the mzXML filename
				my $mzxmlfn = _get_mzxmlfile($tmppath);
				$cp_extracmds .= " -genxml $mzxmlfn";
			}
			## Parallel processing selected
			if ( $PARAMS->{'parallel'} )
			{
				## read the matchfiles from file
				my $parallelfn = File::Spec->catfile( $tmppath, "parallel_matchfiles" );
				my @cpfilesarray = read_file($parallelfn);
				print "Parallel processing of @cpfilesarray\n";
				foreach my $matchlistbn (@cpfilesarray)
				{
					my $matchlist                   = $matchlistbn . ".txt";
					my $converted_spectra_directory = File::Spec->catfile( $tmppath, $dir . "_matcheddir" );
					my $resultdir                   = File::Spec->catfile( $tmppath, $matchlistbn );
					_create_dir( $converted_spectra_directory, $PARAMS, "converted spectra directory" );
					_create_dir( $resultdir,                   $PARAMS, "xQuest result directory" );
					## now change the resultdir and converted spectra dir vars to only the basenames
					$resultdir                   = $matchlistbn;
					$converted_spectra_directory = $dir . "_matcheddir";
					### for xml mode create for every ml a spec xml file
					if ( $PARAMS->{'xmlmode'} )
					{
						## Get the mzXML filename
						#my $mzxmlfn = _get_mzxmlfile($tmppath);
						my $specfilenameout = $matchlistbn . ".spec.xml";
						my $mzxmlfn         = _get_mzxmlfile($tmppath);
						$cp_extracmds = " -genxml $mzxmlfn -specfilenameout $specfilenameout " . $PARAMS->{'cp_extracmds'};
					}
					## generate the cp command
					$comparecmd = "$bin\/compare_peaks3.pl -match $matchlist -dir $converted_spectra_directory -def $xquestdef $cp_extracmds -resultdir $resultdir";
					push @cp_cmds,          $comparecmd;
					push @cp_locationarray, $tmppath;
				}
			} else
			{
				## Normal processing with filelist
				my $matchlist                   = $dir . "_matched.txt";
				my $converted_spectra_directory = File::Spec->catfile( $tmppath, $dir . "_matcheddir" );
				my $resultdir                   = File::Spec->catfile( $tmppath, $dir . "_matched" );
				_create_dir( $converted_spectra_directory, $PARAMS, "converted spectra directory" );
				_create_dir( $resultdir,                   $PARAMS, "xQuest result directory" );
				## change the conveterd specdir to only the basename
				$converted_spectra_directory = $dir . "_matcheddir";
				$resultdir                   = $dir . "_matched";
				$comparecmd                  = "$bin\/compare_peaks3.pl -match $matchlist -dir $converted_spectra_directory -def $xquestdef $cp_extracmds -resultdir $resultdir";
				push @cp_cmds,          $comparecmd;
				push @cp_locationarray, $tmppath;
			}
		}
	} else
	{
		### No listfile was used //matchhlist or xmm with MasterRun must be defined
		## Parallel processing selected
		$tmppath = File::Spec->catfile($rootdir);
		my $dir          = basename($rootdir);
		my $cp_extracmds = $PARAMS->{'cp_extracmds'};
		if ( $PARAMS->{'xmlmode'} )
		{
			## Get the mzXML filename
			my $mzxmlfn = _get_mzxmlfile($tmppath);
			$cp_extracmds .= " -genxml $mzxmlfn";
		}
		if ( $PARAMS->{'parallel'} )
		{
			## read the matchfiles from file
			my $parallelfn = File::Spec->catfile( $tmppath, "parallel_matchfiles" );
			my @cpfilesarray = read_file($parallelfn);
			print "Parallel processing of @cpfilesarray\n";
			foreach my $matchlistbn (@cpfilesarray)
			{
				my $matchlist = $matchlistbn . ".txt";
				### check if the matchlistfile exists
				### Then check if there is a Master Map file
				my ( $status, $msg ) = check_file($matchlist);
				unless ( $status == 1 )
				{
					print "$msg";
					exit;
				}
				### Generate the Directories
				my $converted_spectra_directory = File::Spec->catfile( $tmppath, $matchlistbn . "dir" );
				## converted spe dir basename
				my $converted_spectra_directory_bn = $matchlistbn . "dir";
				my $resultdir                      = File::Spec->catfile( $tmppath, $matchlistbn );
				my $resultdir_bn                   = $matchlistbn;
				_create_dir( $converted_spectra_directory, $PARAMS, "converted spectra directory" );
				_create_dir( $resultdir,                   $PARAMS, "xQuest result directory" );
				## now change the resultdir and converted spectra dir vars to only the basenames
				$comparecmd = "$bin\/compare_peaks3.pl -match $matchlist -dir $converted_spectra_directory_bn -def $xquestdef $cp_extracmds -resultdir $resultdir_bn";
				push @cp_cmds,          $comparecmd;
				push @cp_locationarray, $tmppath;
			}
		} else
		{
			## Normal processing with matchlist
			### check if a matchlist was provided as input parameters
			my $converted_spectra_directory;
			my $resultdir;
			my $matchlist;
			if ( $PARAMS->{'matchlist'} )
			{
				$matchlist = $PARAMS->{'matchlist'};
				my @filenamesplit = split( /\./, $PARAMS->{'matchlist'} );
				my $matchbn = $filenamesplit[0];

				#my $matchbn=basename($matchlist);
				#print ($matchbn);
				$converted_spectra_directory = File::Spec->catfile( $tmppath, $matchbn . "dir" );
				$resultdir                   = File::Spec->catfile( $tmppath, $matchbn );
			} else    ## else use the working dir bn
			{
				$matchlist                   = $dir . "_matched.txt";
				$converted_spectra_directory = File::Spec->catfile( $tmppath, $dir . "_matcheddir" );
				$resultdir                   = File::Spec->catfile( $tmppath, $dir . "_matched" );
			}
			_create_dir( $converted_spectra_directory, $PARAMS, "converted spectra directory" );
			_create_dir( $resultdir,                   $PARAMS, "xQuest result directory" );
			## change the conveterd specdir to only the basename
			#$converted_spectra_directory = $dir . "_matcheddir";
			#$resultdir                   = $dir . "_matched";
			$converted_spectra_directory = basename($converted_spectra_directory);
			$resultdir                   = basename($resultdir);
			$comparecmd                  = "$bin\/compare_peaks3.pl -match $matchlist -dir $converted_spectra_directory -def $xquestdef $cp_extracmds -resultdir $resultdir";
			push @cp_cmds,          $comparecmd;
			push @cp_locationarray, $tmppath;
		}
	}
	return \@cp_cmds, \@cp_locationarray;
}

#---------------------------------------------------------------------------
#  xQuest Section
#---------------------------------------------------------------------------
sub _get_xquest_cmds
{
	my $PARAMS   = shift;
	my $bin      = $PARAMS->{'bin'};
	my $listfile = $PARAMS->{'list'};
	my $rootdir  = $PARAMS->{'rootdir'};
	my $verbose  = $PARAMS->{'verbose'};
	my $tmppath;
	my @xq_cmds;
	my @xq_paths;
#### !! ADD THE  -moveresultdir where needed !!!
	chdir($rootdir);    #change to the rootdir to find the list file
	                    #my $listfilesummary="resultdirectories";

	if ($listfile)
	{
		### generate the listfilesummary
		my $listfilesummary       = File::Spec->catfile( $rootdir, "resultdirectories_fullpath" );
		my $listfilesummary_short = File::Spec->catfile( $rootdir, "resultdirectories" );
		open LISTFILESUMMARY,  ">$listfilesummary"       or die "cannot open $listfilesummary for writing! $!\n";
		open LISTFILESUMMARY2, ">$listfilesummary_short" or die "cannot open $listfilesummary for writing! $!\n";
		my @directories = read_file($listfile);
		foreach my $dir (@directories)
		{
			$tmppath = File::Spec->catfile( $rootdir, $dir );
			### Check if parallel processing is selected
			if ( $PARAMS->{'parallel'} )
			{
				## read the matchfiles (basenames) from file
				my $parallelfn = File::Spec->catfile( $tmppath, "parallel_matchfiles" );
				my @cpfilesarray = read_file($parallelfn);
				print "Parallel processing of @cpfilesarray\n";
				foreach my $matchlistbn (@cpfilesarray)
				{
					my $resultdir  = $matchlistbn;
					my $isopairs   = $matchlistbn . "_isotopepairs.txt";
					my $statusfile = $matchlistbn . ".stat";
					my $xqcmd      = _gen_xquest_cmd( $PARAMS, $resultdir, $isopairs, $statusfile );
					push @xq_cmds,  $xqcmd;
					push @xq_paths, $tmppath;
					## print to listfile
					print LISTFILESUMMARY "$dir/$resultdir\n";
					print LISTFILESUMMARY2 "$resultdir\n";
				}
			} else
			{
				$tmppath = File::Spec->catfile( $rootdir, $dir );
				my $resultdir  = $dir . "_matched";
				my $isopairs   = $resultdir . "_isotopepairs.txt";
				my $statusfile = $resultdir . ".stat";
				my $xqcmd      = _gen_xquest_cmd( $PARAMS, $resultdir, $isopairs, $statusfile );
				push @xq_cmds,  $xqcmd;
				push @xq_paths, $tmppath;
				## print to listfile
				print LISTFILESUMMARY "$dir/$resultdir\n";
				print LISTFILESUMMARY2 "$resultdir\n";
			}
		}
		close(LISTFILESUMMARY);
		close(LISTFILESUMMARY2);
	} else
	{
		## no listfile
		$tmppath = File::Spec->catfile($rootdir);
		my $dir = basename($rootdir);
		### Check if parallel processing is selected
		if ( $PARAMS->{'parallel'} )
		{
			### generate the listfilesummary
			my $listfilesummary       = File::Spec->catfile( $tmppath, "resultdirectories" );
			my $listfilesummary_short = File::Spec->catfile( $rootdir, "resultdirectories" );
			open LISTFILESUMMARY,  ">$listfilesummary"       or die "cannot open $listfilesummary for writing! $!\n";
			open LISTFILESUMMARY2, ">$listfilesummary_short" or die "cannot open $listfilesummary for writing! $!\n";
			## read the matchfiles (basenames) from file
			my $parallelfn = File::Spec->catfile( $tmppath, "parallel_matchfiles" );
			my @cpfilesarray = read_file($parallelfn);
			print "Parallel processing of @cpfilesarray\n";
			foreach my $matchlistbn (@cpfilesarray)
			{
				my $resultdir  = $matchlistbn;
				my $isopairs   = $matchlistbn . "_isotopepairs.txt";
				my $statusfile = $matchlistbn . ".stat";
				my $xqcmd      = _gen_xquest_cmd( $PARAMS, $resultdir, $isopairs, $statusfile );
				push @xq_cmds,  $xqcmd;
				push @xq_paths, $tmppath;
				## print to listfile
				print LISTFILESUMMARY "$dir/$resultdir\n";
				print LISTFILESUMMARY2 "$resultdir\n";
			}
			close(LISTFILESUMMARY);
			close(LISTFILESUMMARY2);
		} else
		{
			### Check if a matchlist was provided as input parameter
			my $resultdir;
			my $isopairs;
			my $statusfile;
			## check if a Master map was extracted
			#if ( $PARAMS->{'xmm_master'} )
			#{
			if ( $PARAMS->{'matchlist'} )
			{    ## check if a matchlist was provided by input param
				my $basenamematchlist = $PARAMS->{'matchlist'};
				#$basenamematchlist =~ s/\.\w+//;    ## chop the .txt #Doesn't work on paths with . in the name
				my ($basenamematchlist_file, $basenamematchlist_dir) = File::Basename::fileparse($basenamematchlist);
				$basenamematchlist = $basenamematchlist_dir . $basenamematchlist_file;
				$tmppath    = File::Spec->catfile($rootdir);
				$resultdir  = $basenamematchlist;
				$isopairs   = $basenamematchlist . "_isotopepairs.txt";
				$statusfile = $basenamematchlist . ".stat";
			} else
			{
				## processing started with a MM. take the dirname as bn
				$tmppath    = File::Spec->catfile($rootdir);
				$resultdir  = $dir . "_matched";
				$isopairs   = $resultdir . "_isotopepairs.txt";
				$statusfile = $resultdir . ".stat";
			}
			my $xqcmd = _gen_xquest_cmd( $PARAMS, $resultdir, $isopairs, $statusfile );
			push @xq_cmds,  $xqcmd;
			push @xq_paths, $tmppath;
		}
	}
	return \@xq_cmds, \@xq_paths;
}
### used for xmlmode to find the mzxml
sub _get_mzxmlfile
{
	my $tmppath = shift;
## Get the mzXML filename
	$tmppath = File::Spec->catfile( $tmppath, "*.mzXML" );

	#print "$tmppath\n";
	my @mzXML_files = glob($tmppath);
	my $mzxmlfn     = basename( $mzXML_files[0] );
	unless ($mzxmlfn)
	{
		die "xmlmode: No mzXML found in $tmppath\n";
	}
	if ( scalar(@mzXML_files) > 1 )
	{

		#print Dumper(@mzXML_files);
		die "xmlmode ERROR: There is more than one mzXML file present in $tmppath\n";
	}
	return $mzxmlfn;
}

sub _gen_xquest_cmd
{
	my $PARAMS     = shift;
	my $resultdir  = shift;
	my $isopairs   = shift;
	my $statusfile = shift;
	my $bin        = $PARAMS->{'bin'};
	my $extracmds  = $PARAMS->{'xq_extracmds'};
	my $xquestdef  = $PARAMS->{'xquestdef'};
	if ( $PARAMS->{'xmlmode'} )
	{
		$extracmds .= " -specxml $resultdir.spec.xml";
	}
	my $runXquest = "$bin\/xquest.pl $extracmds -def $xquestdef -xquestdir $bin -list $isopairs -resdir $resultdir";
	return $runXquest;
}

sub _create_dir
{
	my $dirpath = shift;
	my $PARAMS  = shift;
	my $msg     = shift;
	my $verbose = $PARAMS->{'verbose'};
	unless ( -e $dirpath )
	{
		$verbose && print "Create $msg directory: $dirpath\n";
		mkdir($dirpath);
	} else
	{
		$verbose && print "$msg: $dirpath exists\n";
	}
	return;
}

sub unique_id
{
	my $sessionId = "";
	my $length    = 16;
	for ( my $i = 0 ; $i < $length ; )
	{
		my $j = chr( int( rand(127) ) );
		if ( $j =~ /[a-zA-Z0-9]/ )
		{
			$sessionId .= $j;
			$i++;
		}
	}
	return $sessionId;
}
### Function to split the Matchfile of xmm.pl
sub _split_matchlist
{
	my $basename      = shift;
	my $mlfullpath    = shift;
	my $PARAMS        = shift;
	my @array         = read_file($mlfullpath);
	my $linesperfile  = $PARAMS->{'parallel'};
	my $verbose       = $PARAMS->{'verbose'};
	my $test          = $PARAMS->{'test'};
	my $numofelements = scalar(@array);
	my $numberoffiles = int( $numofelements / $linesperfile );
	my $remaining     = $numofelements % $linesperfile;

	#if ($remaining){$numberoffiles++;}
	my @matchfiles;
	my $parallelfilename = "parallel_matchfiles";
	($verbose) && print "Splitting Matchfile into pieces of $linesperfile lines per file\n";
	($verbose) && print "Total number of lines: $numofelements\n";
	($verbose) && print "Number of files: $numberoffiles\n";
	($verbose) && print "Remaining for last file $remaining\n";
	### Save the parts to files
	my $start = 0;
	my $end   = $start + $linesperfile;
	### Make a fresh parallel out file
	open MYOUTFILE, ">", "$parallelfilename" or die "Cannot open file $parallelfilename $!";
	close(MYOUTFILE);
	if ( $numberoffiles == 0 )
	{
		print "Parallel processing: Nothing to split\n";
		save_to_file( $parallelfilename, $basename . "\n", 1 );
		push( @matchfiles, $basename . ".txt" );
		return \@matchfiles;
	}
	if ($remaining)
	{
		$numberoffiles++;
	}
	for ( my $i = 0 ; $i < $numberoffiles ; $i++ )
	{
		my $filename = $basename . "$i.txt";
		## save the filename to the parallel file
		save_to_file( $parallelfilename, $basename . "$i" . "\n", 1 );
		push( @matchfiles, $filename );
		print "Saving partial Matchfile to: $filename start:$start end:$end-1\n";
		## CREATE A NEW matchfile for the partial list
		open MYOUTFILE, ">", "$filename" or die "Cannot open file $filename $!";
		close(MYOUTFILE);
		foreach my $line ( $start .. ( $end - 1 ) )
		{

			#print "index is $line filename is $filename\n";
			chomp($line);
			if ( $line eq "" ) { next }    ### jump over empty lines
			my $line = $array[$line] . "\n";
			save_to_file( $filename, $line, 1 );
		}

		#Reset Start and End
		$start = $end;
		$end   = $start + $linesperfile;
		if ( $end > $numofelements )
		{
			## then save from start +remaining
			$end = $start + $remaining;
		}
	}
	return \@matchfiles;
}

#===  FUNCTION  ================================================================
#  NAME:        read_file
#  PURPOSE:     read_file
#  DESCRIPTION: Read file line by line and put results into an array
#  PARAMETERS:  $filename (including path if script is not executed in the folder)
#  RETURNS:     array
#===============================================================================
sub read_file
{
	my $filename = shift;
	my $verbose  = shift;
	my @array;
	if ($verbose) { print "Reading from file $filename\n" }
	open FILE, $filename or die "File $filename not found $!";
	while ( my $line = <FILE> )
	{
		if ($verbose) { print "Reading line $line" }
		chomp($line);
		push( @array, $line );
	}
	($verbose) && print "\n";
	close FILE;
	return @array;
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
		open MYOUTFILE, ">>", "$filename" or die "cannot open file $filename $!";
	} else
	{
		open MYOUTFILE, ">", "$filename" or die "cannot open file $filename $!";
	}
	print MYOUTFILE $text;
	close MYOUTFILE;
	return 1;
}

sub print_params
{
	my $hashref = shift;
	foreach my $key ( sort keys %$hashref )
	{
		my $value = $hashref->{$key};
		if ( ref($value) )
		{
			$value = $$value;
		}
		unless ($value) { $value = "not defined" }
		my $length = length($key);
		if ( $length > 6 )
		{
			print "$key \t  =>  $value\n";
		} else
		{
			print "$key \t\t  =>  $value\n";
		}
	}
}

#===  FUNCTION  ================================================================
#  NAME:        check_file
#  PURPOSE:     Check file for existence and readability
#  DESCRIPTION: Check file for existence and readability
#  PARAMETERS:  $filename (including path if script is not executed in the folder)
#  RETURNS:     void
#===============================================================================
sub check_file
{
	my $filename = shift;
	my $verbose  = shift;
	my $msg;
	my $status;

	# -e is for exists, -r for readable
	unless ( ( -e $filename ) && ( -r $filename ) )
	{
		$msg    = "Error: Cannot find/read the file $filename.\n";
		$status = -1;
	} else
	{
		$msg    = "Checking file $filename ok.\n";
		$status = 1;
	}
	return ( $status, $msg );
}

sub print_affi_and_changelog
{
	print "\n############################################## ", basename($0), " ##############################################\n";
	my $scriptinfo = shift;
	print "version " . $scriptinfo->{'version'} . " written by ";
	print $scriptinfo->{'author'} . "\n";
	print "Affiliation: " . $scriptinfo->{'affi'} . "\n";
	print "In case of troubles mailto: ", $scriptinfo->{'mailto'} . "\n";
	print "Changelog:\n";
	foreach my $key ( sort { $a <=> $b } keys %{ $scriptinfo->{'clog'} } )
	{
		print "Version " . $key . ": " . $scriptinfo->{'clog'}->{$key} . "\n";
		print "###########################################################################################################\n";
		print '
### Subversion ###################################################
# $Revision: 200 $
# $LastChangedDate: 2012-12-03 22:17:53 +0100 (Mon, 03 Dec 2012) $
# $Author: walzthoeni $
# $LastChangedBy: walzthoeni $
##################################################################
';
	}
	return;
}

sub usage()
{
	print "
	SOFTWARE: ", basename($0), " version $scriptinfo->{'version'}
	
	AUTHOR: Thomas Walzthoeni

	INFORMATION: A wrapper to run the xQuest search pipeline. It generates (creates *.sh or *.bat files) or executes (optional) the sequence of commands to
	run a xQuest search, or multiple searches. The directory structure can be generated by pQuest.pl.
	The program can be used with a filelist, to search multiple mzXML files or in a single search directory. 	
 	
 	WORKFLOW:
	A. Extract spectra with MzXml2Search (tool from the TPP, must be installed)
	--> is skipped if a spectrumdirectory already exists or -xmlmode is used, can be forced with -force option
	B. Search for isotopic scan pairs (xmm.pl).
	--> is skipped if a *_matched.txt file already exists, can be forced with -xmmforce option
	C. Comparision of light and heavy MS/MS scans by compare_peaks3.pl 
	--> is skipped if a *_isotopepairs.txt file already exists, can be forced with -cpforce option
	D. xQuest Search
 	
 	USAGE: ", basename($0), " -Option [Parameter]
	
	REQUIRED OPTIONS: -list [] or -matchlist [] or -master [./MASTER_RUN/MASTER_RUN.txt]
	-list [] filename with a list of the folder name(s) to search, same list as used by pQuest.pl. Execute the command in the root folder.
	Info: This is the standard way to execute or prepare a search. For further options, using the -matchlist or -master option see below at FURTHER OPTIONS.
	
	THE following options are recommended:
	-xmlmode option, spectra will be saved in a xml file. 
	-pseudoSH, scannumbers are directly used from pseudosh.pl output file (MasterMap with scannumbers)
	-sh option generate runxq.sh (linux) or runxq.bat (Win) files for submission to a cluster or queuing system, otherwise commands are executed directly on the shell.
	
	SPECIFIC OPTIONS
	
	-getdef copy template definition files to the current folder.
		
	Step B:	xmm.pl options
	-xmmonly run only xmm.pl then exit
	-xmmforce force to run xmm
	-pseudoSH ->scannumbers are directily used from pseudosh output (MasterMap with scannumbers)
				
	Step C: compare_peaks3.pl options
	-cponly run comparepeaks only then exit
	-nocompare (if set) then compare_peaks3.pl is not run
	-cpforce force compare peaks to re-run

	EXECUTION options
	-sh option generate runxq.sh (linux) or runxq.bat (Win) files for submission to a cluster or queuing system. 
	Without the -sh option the series of commands is executed directly by the program.
	
	PARALLEL PROCESSING options
	-parallel [int] split the matchlist into equal parts for semi-parallel processing

	RESULT DIRECTORIES and spectrum directories are automatically generated
	-results [] result directory
	-dir [] directory to hold processed spectra
	                          
	DEFINITION files              
	-def [xquest.def] xquest.def file to use
	-xmmdef [xmm.def] xmm.def file to use
	    		
	DEBUG options
	-testcmd print commands without executing
	            	
	## BRUTUS cluster options	
	-timetolive option to set the timeflag (in h) for bsub
	-brutus usage on Brutus cluster
	
	FURTHER OPTIONS
	-master [./MASTER_RUN/MASTER_RUN.txt] MasterMap file, (if used in this mode you have to be in the direcory where the mzXML is), starts at step B
	-matchlist [] Matchlistfile (if used in this mode you have to be in the direcory where the mzXML is), starts at step C
	
	EXAMPLE: ", basename($0), " -list files -pseudoSH -xmlmode -sh
	";
	exit;
}
