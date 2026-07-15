#!/usr/bin/env perl
use strict;
use warnings;

#---------------------------------------------------------------------------
# xquest.pl
# A software to identify cross-linked peptides from MS/MS spectra.
# Execute xquest.pl -help to display information and usage options.
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

#System libraries
use File::Path;
use MIME::Base64;
use Getopt::Long;
use File::Basename;
use File::Spec;
use File::Copy;
use Data::Dumper;
use Cwd;

##########################################################
# include the directory for the xquest modules ###########
use File::Spec::Functions qw(rel2abs);
use File::Basename;
use FindBin;
use lib "$FindBin::Bin/../modules";
##########################################################

#use lib "$FindBin::Bin/../../perl5";
use Combinatorics;
use PepObj;
use Match;
use Spectrum;
use Index;
use Ionseries;
use Environment;
use Read_Params;
use XML::TreeBuilder;
use MIME::Base64;

# omit output buffering ## useful on cluster
$| = 1;
my $OS = $^O;    # get the current OS
my $xquestdir;
my $version = "xquest 2.1.7";
if ( defined( $ENV{'XQUEST_DIR'} ) )
{
	$xquestdir = $ENV{'XQUEST_DIR'};
} else
{

	#warn "environment variable XQUEST_DIR is not defined $!, setting it to .";
	$xquestdir = '.';
}
use vars qw($PARAMS $MSTAB $ENZ $WEBPARAMS );
my ( $pepobj, $verbose, $outfile, $help, $spectrum, $i, $j, $parention, $superverbose, @spectrum, @scantypes, @isotopeshifts, @rttimesscans, @mzscans, $newindex,, $onlyindex, $stdout, $basename, $speclist, $resultdir, @xlinkhits, $statusfilehandle, $spectrumdir, $benchmark, $readversion,
	 $enumeration, $enumerationindexname, $enumerationindexname_mono, $enumerationindexname_intra, $masslist, $indexinfo, $dbindex, $directory, $xquestdef, $masstable, $webconfig, $statusfile, $progressfile, $msg, $webversion, $moveresultdir, $resultdirfinal, $specxml, $owndb );
GetOptions(
			'def=s'           => \$xquestdef,
			'info'            => \$indexinfo,
			'out=s'           => \$outfile,
			'spec=s'          => \$spectrum,
			'xquestdir=s'     => \$xquestdir,
			'masstab=s'       => \$masstable,
			'masslist=s'      => \$masslist,
			'list=s'          => \$speclist,        ### is used by xQuest pipeline
			'dbindex'         => \$dbindex,
			'nidx'            => \$newindex,
			'enum'            => \$enumeration,
			'oidx'            => \$onlyindex,
			'resdir=s'        => \$resultdir,
			'dir=s'           => \$spectrumdir,
			'sv'              => \$superverbose,
			'version'         => \$readversion,
			'stdout'          => \$stdout,
			'bm'              => \$benchmark,
			'verbose'         => \$verbose,
			'stat=s'          => \$statusfile,
			'help'            => \$help,
			'progress=s'      => \$progressfile,
			'webversion'      => \$webversion,
			'moveresultdir=s' => \$moveresultdir,
			'specxml=s'       => \$specxml,
			'owndb'           => \$owndb,
);
&usage if $help;

## Use the Env object to get the Paths
my $env = Environment->new;
$webconfig = $env->get_path('web.config');


if ( !$onlyindex )
{
	&usage unless ($resultdir);
	&usage unless ( $speclist || $spectrum );
}
unless ($xquestdef)
{
	$xquestdef = File::Spec->catfile( $xquestdir, "xquest.def" );
}
unless ($masstable)
{
	#$masstable = File::Spec->catfile( $xquestdir, "mass_table.def" );
	$masstable = $env->get_path('mass.def');	
}

unless ($webconfig)
{
	$webconfig = File::Spec->catfile( $xquestdir, "web.config" );
}

# open masstable and definition table
# basename is the database path/basename is used for the db indices
( $MSTAB, $PARAMS, $ENZ, $basename, $WEBPARAMS ) = Read_Params::readtables( $masstable, $xquestdef, $webconfig, $masslist, $xquestdir, $verbose );

#$verbose=1;
if ($verbose)
{
	Read_Params::print_params($PARAMS);
}
### Define the Stat and progressfile locations (may be overwritten with the -moveresdir option)
# get the location of the current directory
my $dir = getcwd;
######
# moveresultdir option allows to place the results directly in the webserver resultdirectory
# Progress and stat files are directly written there
# resultfolder is copied then at the end of the search
if ($moveresultdir)
{

	# get the resultdirbase from the webparams
	my $resultdirbase = $WEBPARAMS->{'resultdirbase'};

	# cat the username to the resultdirectory
	my $user_name = get_username();

	#print $user_name;
	my $homedir;
	if ( $moveresultdir ne "" )
	{
		$homedir = File::Spec->catfile( $resultdirbase, $user_name, $moveresultdir );
	} else
	{
		$homedir = File::Spec->catfile( $resultdirbase, $user_name );
	}
	unless ( -e $homedir )
	{
		warn "directory $homedir does not exist, create new...\n";
		mkdir($homedir) or die $!;
	}

	# Specify the status- and progressfilename and if not defined by a param
	unless ($statusfile)
	{
		$statusfile = "$resultdir.stat";
	}
	unless ($progressfile)
	{
		$progressfile = "$resultdir.progress";
	}

	# Specifies the $resultdir folder where the results are stored
	$resultdirfinal = File::Spec->catfile( $homedir, $resultdir );
	## remove the resultdir if it already exists.
	if ( -e $resultdirfinal )
	{
		print "Resultdirectory already exists, will remove it \n";
		eval { rmtree( $resultdirfinal, 1, 0 ) };
		if ($@)
		{
			print "Couldn't delete $resultdirfinal: $@";
		}
	}
	unless ( -e $resultdirfinal )
	{
		print "Resultdirectory $resultdirfinal does not exist, create new...\n";

		#mkdir($resultdirfinal) or die $!;
		## must be writable to see the results
		eval { mkpath( $resultdirfinal, 1, 0775 ) };
		if ($@)
		{
			print "Couldn't create $resultdirfinal: $@";
		}
	}

	# specify the path for the statusfile and progressfile
	$statusfile   = File::Spec->catfile( $resultdirfinal, $statusfile );
	$progressfile = File::Spec->catfile( $resultdirfinal, $progressfile );
} else
{

	# standard file names and locations for stat and progress files are defined
	# usually the statfile name is defined as a param, if not define standard statfile name
	unless ($statusfile)
	{
		$statusfile = File::Spec->catfile( $dir, $resultdir, "$resultdir.stat" );
	}
	$progressfile = File::Spec->catfile( $dir, $resultdir, "$resultdir.progress" );
}
## EOF copydirect section
######
if ($webversion)
{
	unless ($progressfile)
	{
		my @filename = split( "\\\\", $resultdir );
		print $resultdir. "<br>";
		my $foldername = pop(@filename);
		$progressfile = File::Spec->catfile( $resultdir, "$foldername.progress" );
	}
}

# creates a standard file name for the progress file
#print "STAT file location: $statusfile\n";
#print "PROGRESS file location: $progressfile\n";
#exit;
printversion($version) if $readversion;
usage($version) unless $spectrum || $speclist || $onlyindex || $masslist || $indexinfo;
usage($version) if $help;
my $time1 = time;
if ($stdout)
{
	$statusfilehandle = *STDOUT;
} else
{
	open STATUS, ">", "$statusfile" or die "cannot open status file $statusfile $!";
	$statusfilehandle = *STATUS;
}
if ($resultdir)
{    #override param file
	$PARAMS->{'outputpath'} = $resultdir;
}
if ( $PARAMS->{'outputpath'} )
{
	$directory = $PARAMS->{'outputpath'};
	unless ( -e $directory )
	{
		print "Cannot find directory $directory ... make new one \n";
		mkdir $directory;
	}
	copy( $xquestdef, File::Spec->catfile( $directory, basename($xquestdef) ) )
	  or warn "copying of $xquestdef to $directory was not successfull: $!";
}
#### Generation of the database ######################################################
######### open protein database, search also in the current directory ################
my $indb;
my $databasefile;
my $indb_dc;
my $databasefile_dc;
if ( -e $PARAMS->{'database'} )
{
	$indb = Bio::SeqIO->new( -file => $PARAMS->{'database'}, -format => 'fasta' );
	$databasefile = $PARAMS->{'database'};
	if ( $PARAMS->{'RuntimeDecoys'} )
	{
		if ( -e $PARAMS->{'database_dc'} )
		{
			$indb_dc = Bio::SeqIO->new( -file => $PARAMS->{'database_dc'}, -format => 'fasta' );
			$databasefile_dc = $PARAMS->{'database_dc'};
		} else
		{
			die "Cannot open decoy database file", $PARAMS->{'database_dc'}, "$!";
		}
	}
} elsif ( -e File::Spec->catfile( $xquestdir, $PARAMS->{'database'} ) )
{
	$indb = Bio::SeqIO->new( -file => File::Spec->catfile( $xquestdir, $PARAMS->{'database'} ), -format => 'fasta' );
	$databasefile = File::Spec->catfile( $xquestdir, $PARAMS->{'database'} );
} else
{
	die "cannot open database file ", $PARAMS->{'database'}, "$!";
}
my $basename_dc;
if ( $dbindex || $PARAMS->{'copydb2resdir'} )
{
	copy( $databasefile, $directory )    #copy database to resultdirectory if required
	  or die "copying of $databasefile to $directory was not successfull: $!";
	#### ALSO COPY THE DB TO THE CWD/DB, MAKE NEW INDEX, AND CHANGE THE DB PATH
	my $dbdir = File::Spec->catfile( $dir, $resultdir, "db" );
	_create_dir( $dbdir, 1, "Database Directory" );
	copy( $databasefile, $dbdir );       #copy database to resultdirectory if required
	my $bndb = basename($databasefile);
	$databasefile = File::Spec->catfile( $dbdir, $bndb );    ### change the path of the db to cwd/db/db.fasta
	my $dbcopy = $databasefile;
	#$dbcopy =~ s/\.\w+//; #Doesn't work on paths with . in the name
	my ($dbcopy_file, $dbcopy_dir) = File::Basename::fileparse($dbcopy);
	$dbcopy = $dbcopy_dir . $dbcopy_file;
	### change the db basename-> is used for the db indexes
	$basename = $dbcopy;
	if ( $PARAMS->{'RuntimeDecoys'} )
	{
		copy( $databasefile_dc, $directory )                 #copy dc database to resultdirectory if required
		  or die "copying of $databasefile to $directory was not successfull: $!";
		copy( $databasefile_dc, $dbdir );
		my $bndb_dc = basename($databasefile_dc);
		$databasefile_dc = File::Spec->catfile( $dbdir, $bndb_dc );
		my $dbcopy_dc = $databasefile_dc;
		#$dbcopy_dc =~ s/\.\w+//; #Doesn't work on paths with . in the name
		my ($dbcopy_file, $dbcopy_dir) = File::Basename::fileparse($dbcopy);
		$dbcopy = $dbcopy_dir . $dbcopy_file;
		$basename_dc = $dbdir; #$dbcopy_dc;
		print "Target Database is: " . $dbcopy . "\n";
		print "Decoy database is: " . $dbcopy_dc . "\n";
	}
}
########################################################################################
#prepare index files names #############################################################
my $ionindexname = $basename . "_ion.db";
my $pepindexname = $basename . "_peps.db";
### check out the params which DB is required ##########################################
$enumeration = $PARAMS->{'enumerate'};
my $iontag = $PARAMS->{'Iontagmode'};
### generate additionalfilenames for enumeration DB ####################################
if ($enumeration)
{
	$enumerationindexname       = join "_", $basename, "xlink_enum.db";
	$enumerationindexname_mono  = join "_", $basename, "mono_enum.db";
	$enumerationindexname_intra = join "_", $basename, "intra_enum.db";
}
######## MAKE THE NEW DATABASE #########################################################
#my $index;
my $index = Index->new( $indb, $basename, $statusfilehandle, $PARAMS, $ENZ, $MSTAB, $verbose, $superverbose, $newindex, $onlyindex, $enumeration, $iontag, $indexinfo );
### GET THE DATABASES INTO VARS #####
my $IONINDEX  = $index->get_ionindex;
my $PEPINDEX  = $index->get_pepindex;
my $INFOINDEX = $index->get_infoindex;
######## MAKE INDEX OF THE DC DATABASE
my ( $IONINDEX_DC, $PEPINDEX_DC, $INFOINDEX_DC, $index_dc );
if ( $PARAMS->{'RuntimeDecoys'} )
{
	print "\nGenerating new index for dc db\n";
	$index_dc     = Index->new_dc( $indb_dc, $basename_dc, $statusfilehandle, $PARAMS, $ENZ, $MSTAB, $verbose, $superverbose, $newindex, $onlyindex, $enumeration, $iontag, $indexinfo );
	$IONINDEX_DC  = $index_dc->get_ionindex_dc;
	$PEPINDEX_DC  = $index_dc->get_pepindex_dc;
	$INFOINDEX_DC = $index_dc->get_infoindex_dc;
}

#print $INFOINDEX;
#print $statusfilehandle "indexed ", $index->getnprots, " proteins\nindexed ",
#  $index->getnpeps, " peptides\nindexed ", $index->getnions, " ions\n";
my $ENUMERATION_INTERXLINKS = $index->get_enumeration_inter_xlinks_index;
my $ENUMERATION_INTRAXLINKS = $index->get_enumeration_intra_xlinks_index;
my $ENUMERATION_INTRA       = $index->get_enumeration_intralinks_index;
my $ENUMERATION_MONO        = $index->get_enumeration_monolinks_index;

#search indicated spectra
my $xlinkoutfile;
if ($masslist)
{
	open MASSLIST, "<$masslist" or die "cannot open masslist list $masslist!\n";
	while (<MASSLIST>)
	{
		chomp;
		my @tmp = split;
		if ( $tmp[0] =~ /\d+\.\d+/ && $tmp[1] =~ /\d/ )
		{

			#print "masslist: mz",$tmp[0]," charge: ",$tmp[1],"Mr: ",$tmp[0]*$tmp[1]-$tmp[1]*$MSTAB->{'Hatom'}->{'native'},"\n";
			push @spectrum, $tmp[0] * $tmp[1] - $tmp[1] * $MSTAB->{'Hatom'}->{'native'};
		} elsif ( $tmp[0] =~ /\d+\.\d+/ )
		{
			push @spectrum, $tmp[0];
		} else
		{
			warn "not a number ", $tmp[0], "$0";
		}
	}
	close(MASSLIST);
} elsif ($spectrum)
{
	chomp($spectrum);
	push @spectrum,      File::Spec->catfile($spectrum);
	push @isotopeshifts, 0;
	$xlinkoutfile = join "", $spectrum, "_xlinks.txt";
} elsif ($speclist)
{
	$xlinkoutfile = $speclist;
	$xlinkoutfile =~ s/.txt/_xlinks.txt/;
	open SPEC, "<$speclist" or die "cannot open spectrum list $speclist!\n";
	while (<SPEC>)
	{
		chomp;
		my ( $specname, $isotopeshift, $scantype, $rttimes, $mzscan );
		( $specname, $isotopeshift, $scantype, $rttimes, $mzscan ) = split;
		unless ($isotopeshift)
		{
			$isotopeshift = 0;
		}
		my $specfile = File::Spec->catfile( ".", $spectrumdir, $specname ); # or die "cannot open $specname from $speclist\n";
		push @spectrum,      $specfile;
		push @isotopeshifts, $isotopeshift;
		push @scantypes,     $scantype;
		push @rttimesscans,  $rttimes;
		push @mzscans,       $mzscan;
	}
}
close(SPEC);
my $xmlfile;
#### Define the xml out filename for search results
unless ($outfile)
{
	$xmlfile = File::Spec->catfile( $PARAMS->{'outputpath'}, "xquest.xml" );
	open XLINKS, ">$xmlfile"
	  or die "cannot open $xmlfile for writing! $!\n";
} else
{
	$xmlfile = File::Spec->catfile( $PARAMS->{'outputpath'}, ( join "", $outfile, ".xml" ) );
	open XLINKS, ">$xmlfile" or die "cannot open $xmlfile for writing! $!\n";
}
openXMLheader( *XLINKS, $PARAMS, $WEBPARAMS );    #writes a xml
#### Define the xml out filename for the spec.xml file
#### ONLY GENERATE IF $specxml IS NOT DEFINED
#### OTHERWISE PARSE THE SPECTRA INTO A HASH
my $filenamespec;
my $xmlspectfile;
my $spechash = {};
unless ($specxml)
{
	my $xmlspectfile;                             ##full path
	my $specoutfile;                              ##filename
#### CREATING A SPEC.XML FILE OR PARSE A SPEC.XML FILE IF PROVIDED
	unless ($specoutfile)
	{
		## outputpath is already the resultdirectory
		$filenamespec = "$resultdir" . ".spec." . "xml";
		$xmlspectfile = File::Spec->catfile( $PARAMS->{'outputpath'}, $filenamespec );
		open XMLSPECFILE, ">$xmlspectfile" or die "cannot open $xmlspectfile for writing! $!\n";
	} else
	{
		$filenamespec = "$specoutfile" . "xml";
		$xmlspectfile = File::Spec->catfile( $PARAMS->{'outputpath'}, ( join "", $specoutfile, ".xml" ) );
		open XMLSPECFILE, ">$xmlspectfile" or die "cannot open $xmlspectfile for writing! $!\n";
	}
	openXMLspecHeader(*XMLSPECFILE);
} else
{
### OTHERWISE CHECK IF THE SPECFILE EXISTS AND PARSE THE SPECTRA INTO A HASH
	$PARAMS->{'specxml'} = $specxml;
	$filenamespec = $specxml;
	$xmlspectfile = File::Spec->catfile( $PARAMS->{'outputpath'}, $filenamespec );
	unless ( ( -e $xmlspectfile ) && ( -r $xmlspectfile ) )
	{
		die "Cannot read SpecXML file $xmlspectfile\n";
	} else
	{
		print "\nXML Mode: SpecXML file $xmlspectfile found and parsed.\n";
### PARSING
		$spechash = parse_spec_xml_file($xmlspectfile);
	}
}
$i = 0;
my $nspectra = $#spectrum + 1;
my $xcorr_findmaxpeak = $PARAMS->{'xcorr_tolerance_window'} || 0;
#### SEARCH EVERY SPECTRUM ###
for $i ( 0 .. $#spectrum )
{

	#print "xQuest searching $i of $nspectra...\n";
	if ($progressfile)
	{
		my $nthspec = $i + 1;
		$msg = "\n ### Searching spectrum $nthspec of $nspectra... ###";
		print $msg. "\n";
		printprogress( $progressfile, $msg );
	}
	my $spec          = $spectrum[$i];
	my $isotopicshift = $isotopeshifts[$i];
	my $scantype      = $scantypes[$i];
	my $rttimes       = $rttimesscans[$i];
	my $mzs           = $mzscans[$i];
	$PARAMS->{'addedmass'} = $isotopicshift;
	$isotopicshift && print $statusfilehandle "\nsearching ", basename($spec), " normalizing x-linker peaks with isotopic shift of $isotopicshift\n";
	###### INITIALIZE A SPECTRUM OBJECT #######
	my $spectrumObj;
	if ($masslist)
	{
		$spectrumObj = Spectrum->newmasslist( $PARAMS, $spec );
	} else
	{
		### Create a Spectrum object using the spechash (empty if no spec.xml was provided)
		$spectrumObj = Spectrum->new( $spec, $PARAMS, $scantype, $isotopicshift, $spechash );
	}
	print $statusfilehandle "Searching spectrum ", ( $i + 1 ), "\/$nspectra: ", $spectrumObj->getSpecname(), "\n";
	my $timestart = time;
	###### INITIALIZE A MATCHOBJECT ###########
	my $matchobj = Match->new( $spectrumObj, $IONINDEX, $PEPINDEX, $PARAMS, $MSTAB, $verbose, $statusfilehandle );
	### add the INFOINDEX to the matchobject
	### Add the infoindex to the Matchobject
	$matchobj->set_infoindex( "infoindex_target", $INFOINDEX );

	#---------------------------------------------------------------------------
	#  Generation of the dc matchobject
	#---------------------------------------------------------------------------
	my $matchobj_dc;
	if ( $PARAMS->{'RuntimeDecoys'} )
	{

		#print "Creating dc matchobject\n";
		$matchobj_dc = Match->new( $spectrumObj, $IONINDEX_DC, $PEPINDEX_DC, $PARAMS, $MSTAB, $verbose, $statusfilehandle );
		$matchobj_dc->set_infoindex( "infoindex_decoy", $INFOINDEX_DC );
		#### Generate only full decoy hits or also target-decoy hybrid hits ####
		unless ( $PARAMS->{'Fulldecoysonly'} )
		{
			print "Geneating hybrid decoys: Adding DC peps to target matchlist\n";
			$matchobj->add_dc_candidates_to_matchlists( $matchobj_dc->hitlist, $matchobj_dc->gethithash );
			### $matchobj->set_infoindex("infoindex_decoy", $INFOINDEX_DC);
			### for hybrids add the infoindex for apriori score calculations
			$matchobj->set_infoindex( "infoindex_decoy", $INFOINDEX_DC );
		}
	}
	$benchmark && ( print $statusfilehandle "matching took ", time - $timestart, "s\n" );
	my $timestart2 = time;
	if ( $PARAMS->{'printcandidatepeps'} )
	{
		$matchobj->printhits($statusfilehandle);
	}
	### GENERATION OF THE xl candidates (Search)
	my $xlinkhits;
	my $xlinkhits_dc;
	if ($enumeration)
	{
		$xlinkhits = $matchobj->makexlinks_enumerate( $PARAMS, $statusfilehandle, $verbose, $PEPINDEX, $ENUMERATION_INTERXLINKS, $ENUMERATION_INTRAXLINKS, $ENUMERATION_MONO, $ENUMERATION_INTRA, $MSTAB, $INFOINDEX );
	} else
	{
		$xlinkhits = $matchobj->makexlinks_iontag( $PARAMS, $statusfilehandle, $verbose );
		if ( $PARAMS->{'Fulldecoysonly'} )
		{
			$xlinkhits_dc = $matchobj_dc->makexlinks_iontag( $PARAMS, $statusfilehandle, $verbose );
			## Join the xlink hits
			#print "Number of target xlink candidates: ", scalar (@$xlinkhits),"\n";
			#print "Number of decoy xlink candidates: ", scalar (@$xlinkhits_dc),"\n";
			push @$xlinkhits, @$xlinkhits_dc;    ### Add the decoy candidates to the xlink candidates
			print "Number of target+decoy xlink candidates: ", scalar(@$xlinkhits), "\n";
		}
	}
	$benchmark && ( print $statusfilehandle "xlinksearch took ", time - $timestart2, "s\n" );
	$verbose && print "Generation of Xlinks took ", time - $timestart2, "s\n";
	my $timestart3 = time;
	openXMLresult( *XLINKS, basename($spec), $matchobj->getnhits, $spectrumObj, $rttimes, $mzs );

	#---------------------------------------------------------------------------
	#  Save Spectra to Spectrum XML File(base64 encoded) IF NO SPECXML FILE WAS SPECIFIED
	#---------------------------------------------------------------------------
	unless ($specxml)
	{
		my ( $lightspec, $heavyspec ) = map { File::Spec->catfile( $resultdir, $_ ) } split /,/, $spectrumObj->getspectrumheader;
		save_spectra_to_XML_base64enc( $lightspec, "light", *XMLSPECFILE );
		save_spectra_to_XML_base64enc( $heavyspec, "heavy", *XMLSPECFILE );

		#print "LightSpectrum is: $lightspec\n";
		#print "heavySpectrum is: $heavyspec\n";
		my $commonspectrum = $spectrumObj->{'commonspecname'};
		my $xlspectrum     = $spectrumObj->{'xlinkspecname'};
		save_spectra_to_XML_base64enc( $commonspectrum, "common",  *XMLSPECFILE );
		save_spectra_to_XML_base64enc( $xlspectrum,     "xlinker", *XMLSPECFILE );
	}
	$matchobj = ();
	if ( @$xlinkhits >= 1 )
	{
		my $prescorerank = 0;
		my @scoredhits   = ();
		unless ( $PARAMS->{'massmatchonly'} )
		{
			my %seen = ();
			print "Number of hits " . @$xlinkhits . "\n";
			foreach my $hit ( sort { $b->getprescore <=> $a->getprescore } @{$xlinkhits} )
			{
				if ( defined($hit) )
				{

					#print "HIT TYPE: ".$hit->get_target_decoy_label." HIT ID: ".$hit->getid." Prescore:".$hit->getprescore."\n";
					### HERE ALL THOSE HITS WITH THE SAME ID WHICH IS PEPA-PEPB TOPA TOPB are sorted out
					unless ( $seen{ $hit->getid } )
					{
						## Check if this is a possible xl topology (Param: possibleTopology)
						if ( $PARAMS->{'possibleTopology'} )
						{
							unless ( $hit->check_topology )
							{

								#print "Skipping hit\n";
								next;
							}
						}
						if ( $prescorerank < $PARAMS->{'usenprescores'} )
						{

							#print "HIT TYPE: ".$hit->get_target_decoy_label." HIT ID: ".$hit->getid." Prescore:".$hit->getprescore."\n";
							$hit->calcfullscore($xcorr_findmaxpeak);

							#$hit->print_subscores;
							#print "\n";
							$prescorerank++;
							push @scoredhits, $hit;
							$seen{ $hit->getid } = 1;
						}
					} else
					{
						$seen{ $hit->getid }++;
					}
				}
			}

			#print Dumper (\%seen);
			#exit;
			my $rank = 0;
			print "Sorting ", scalar(@scoredhits), " scored hits\n";
			foreach my $hit ( sort { $b->getscore <=> $a->getscore } @scoredhits )
			{
				$rank++;
				if ( $rank <= $PARAMS->{'reportnbesthits'} )
				{
					my $pics;

					#	print $statusfilehandle "$rank: ";
					$PARAMS->{'printflatresultfile'}
					  && $hit->printxlinks($statusfilehandle);

					#$verbose && $hit->printhittable;
					if ( $PARAMS->{'drawspectra'}
						 && !$PARAMS->{'massmatchonly'} )
					{
						$pics = $hit->drawxlinkspec( $PARAMS->{'drawlogscale'}, ( $PARAMS->{'waterloss'} || $PARAMS->{'nh3loss'} ), undef, undef, undef, 1, 1 );
					}
					$hit->printxlinksXML( *XLINKS, $rank, $pics, $WEBPARAMS, $PARAMS->{'printionmatches'} );
				} else
				{
					last;
				}

				#exit;
			}
		} else
		{
			foreach my $hit ( @{$xlinkhits} )
			{

				#$hit->printxlinks($statusfilehandle);
				$hit->printxlinksXML( *XLINKS, 0, 0, $WEBPARAMS, $PARAMS->{'printionmatches'} );
			}
		}
	}
	closeXMLresult(*XLINKS);
	$spectrumObj = ();
	$benchmark && ( print $statusfilehandle "xlink evaluation took ", time - $timestart3, "s\n" );
	$verbose && print "xlink evaluation took ", time - $timestart3, "s\n";

	#exit (0);
	#$pm->finish($i); # pass an exit code to finish
}

#$pm->wait_all_children; ## wait for the child processes
#exit;
closeXMLheader(*XLINKS);
unless ($specxml)
{
	closeXMLspecHeader(*XMLSPECFILE);
}
close(XLINKS);
close(XMLSPECFILE);
print $statusfilehandle "done, execution time: ", time - $time1, "s\n";
print "done, execution time: ", time - $time1, "s\n";
close(STATUS);
if ($progressfile)
{
	$msg = "Search finished: $nspectra spectra were searched.";
	printprogress( $progressfile, $msg );
	close(PROGRESS);
}

#print "move $statusfile";
print "Copy $statusfile to $statusfile.done\n";
copy( $statusfile, $statusfile . '.done' ) or warn "could not move  $statusfile to $statusfile.done $!";
#################################### end of main #################################
if ($moveresultdir)
{
	my $finaldestination;
	my $resultdirbase = $WEBPARAMS->{'resultdirbase'};
	if ( $moveresultdir ne "" )
	{
		$finaldestination = File::Spec->catfile( $resultdirbase, get_username(), $moveresultdir );

		#$homedir = File::Spec->catfile( $resultdirbase, $user_name, $moveresultdir );
	} else
	{
		$finaldestination = File::Spec->catfile( $resultdirbase, get_username() );
	}
	my $cmd;
	($verbose) && print "copy $xquestdef,$xmlfile,$databasefile to $resultdirfinal\n";
#### Copy xquest.xml,spec.xml, database, and xquest.def to the resultdirectory
	$cmd = "cp $xquestdef $resultdirfinal";
	system($cmd);
	$cmd = "cp $xmlfile $resultdirfinal";
	system($cmd);
	$cmd = "cp $xmlspectfile $resultdirfinal";
	system($cmd);
	$cmd = "cp $databasefile $resultdirfinal";
	system($cmd);
## chmod the resultdirectory to 775
	$cmd = "chmod 775 $resultdirfinal";
	system($cmd);
	print "Delete tmp resultfolder\n";
	eval { rmtree( $resultdir, 0, 0 ) };

	if ($@)
	{
		print "Couldn't delete $resultdir: $@";
	}

	#$cmd="rm -r $resultdir";
	#system ($cmd);
}
################################## delete the tmp folder if webversion param set ########
################################## END OF xQUEST ########################################
################################## FUNCTIONS ############################################
sub _create_dir
{
	my $dirpath = shift;
	my $verbose = shift;
	my $msg     = shift;

	#my $verbose = $PARAMS->{'verbose'};
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

#---------------------------------------------------------------------------
#  SUB parseXML for PARSING THE XML FILE
#---------------------------------------------------------------------------
sub parse_spec_xml_file
{
	my $xmlfilename = shift;
	my $spechash    = {};
	my $tree        = XML::TreeBuilder->new();
	print "Parsing mzXML  $xmlfilename\n";
	$tree->parse_file($xmlfilename);
	my @resultsheader = map { $_ } $tree->find('xquest_spectra');
	### index the spectra by filename
	foreach my $header (@resultsheader)
	{
		my @spectra = map { $_ } $header->find('spectrum');
		foreach my $spectrum (@spectra)
		{
			my $filename = $spectrum->attr('filename');
			my $content  = $spectrum->content();
			$spechash->{$filename} = $content->[0];
		}
	}
	$tree->delete();
	return $spechash;
}
## functions for saving the spectra in the xmlresult
## Save the spectra in a separate file
## print spectrum light/heavy/common/xlink into one tag base 64 encoded
sub printToXML
{
	my $xlinkfile = shift;
	my $content   = shift;
	print $xlinkfile $content;
}

sub read_spectrum
{
	my $path = shift;
	my @array;
	open FILE, "<", $path or die $!;
	while (<FILE>)
	{
		my $line = $_;
		chomp $line;
		push @array, ($line);
	}
	close(FILE);
	return @array;
}

sub slurp
{
	my $file = shift;
	print "Filename: $file\n";
	local ( $/, *FH );
	open( FH, $file ) or die $!;
	my $buffer = <FH>;
	return $buffer;
}

sub save_spectra_to_XML_base64enc
{
	my $filename      = shift;    ##specfilename
	my $typeattribute = shift;
	my $xmlfilename   = shift;    ## ref to filehandle *SPECFILE
	my $text;
	## load the file
	my $buffer = slurp($filename);
	## encode the content
	my $encoded_text = encode_base64($buffer);
	## create an xml element
	my $filenamebn = basename($filename);
	$text = "<spectrum filename=\"$filenamebn\" type=\"$typeattribute\">";
	$text .= $encoded_text;
	$text .= "</spectrum>\n";
	### Store the spectra in the xml file
	printToXML( $xmlfilename, $text );
	return;
}

#--------------------------------- start xml functions --------------------------#
sub openXMLheader
{
	my $xlinkfile = shift;
	my $PARAMS    = shift;
	my $webparams = shift;
	my $date      = localtime;
	my $database  = $PARAMS->{'database'};
	print $xlinkfile "<?xml version=\"1.0\" encoding=\"UTF-8\"?>
<?xml-stylesheet type=\"text\/xsl\" href=\"", $webparams->{'xslt_stylesheet'}, "\"?>
<xquest_results xquest_version=\"$version\" date=\"$date\" author=\"Thomas Walzthoeni,Oliver Rinner\" homepage=\"http://proteomics.ethz.ch\" deffile=\"", basename($xquestdef), "\" ";

	foreach ( keys %$PARAMS )
	{
		if ( defined( $PARAMS->{$_} ) )
		{
			unless (/^\#/)
			{
				print $xlinkfile $_, "=\"", $PARAMS->{$_}, "\" ";
			}
		}
	}
	print $xlinkfile ">\n";
}

sub openXMLspecHeader
{
	my $specoutfile = shift;
	my $date        = localtime;
	print $specoutfile "<?xml version=\"1.0\" encoding=\"UTF-8\"?>";
	print $specoutfile "<xquest_spectra xquest_version=\"$version\" author=\"Thomas Walzthoeni,Oliver Rinner\" homepage=\"http://proteomics.ethz.ch\" resultdir=\"", $resultdir, "\" deffile=\"", basename($xquestdef), "\" ";
	print $specoutfile ">\n";
}

sub closeXMLspecHeader
{
	my $specoutfile = shift;
	print $specoutfile "<\/xquest_spectra>";
}

sub closeXMLheader
{
	my $xlinkfile = shift;
	print $xlinkfile "<\/xquest_results>";
}

sub openXMLresult
{
	my $xlinkfile   = shift;
	my $spectrum    = shift;
	my $ncandidates = shift;
	my $specObj     = shift;
	my $rttimes     = shift;
	my $mzs         = shift;
	if ( !$ncandidates )
	{
		$ncandidates = "";
	}
	my $ncommon          = $specObj->getncommonions;
	my $nxlink           = $specObj->getnxlinkions;
	my $Mrprecursor      = $specObj->getms1mass;
	my $addedMass        = $specObj->get_addedmass;
	my $scantype         = $specObj->get_scantype;
	my $mz_precursor     = $specObj->getprecursorMz;
	my $charge_precursor = $specObj->getprecursorCharge;
	my $apriori_common   = $specObj->get_apriori_pcommon;
	my $apriori_xlink    = $specObj->get_apriori_pxlink;
	print $xlinkfile
	  "<spectrum_search spectrum=\"$spectrum\" mean_ionintensity=\"",
	  $specObj->get_average_ionintensity, "\" ionintensity_stdev=\"",
	  $specObj->get_ionintensity_stdev, "\" addedMass=\"", $addedMass,
	  "\" iontag_ncandidates=\"$ncandidates\"  apriori_pmatch_common=\"",
	  sprintf( "%.4f", $apriori_common ), "\" apriori_pmatch_xlink=\"",
	  sprintf( "%.4f", $apriori_xlink ),
	  "\" ncommonions=\"$ncommon\" nxlinkions=\"$nxlink\" mz_precursor=\"$mz_precursor\" scantype=\"$scantype\" charge_precursor=\"$charge_precursor\" Mr_precursor=\"$Mrprecursor\" rtsecscans=\"$rttimes\" mzscans=\"$mzs\" >\n";
}

sub closeXMLresult
{
	my $xlinkfile = shift;
	print $xlinkfile "<\/spectrum_search>\n";
}
################################# end xml functions ##########################
sub readdefaults
{
	my $defaultfdef = shift;
	my ( %PARAMS, %ENZ, %MOD );
	open DEF, "<$defaultfdef" or warn "cannot open xquest defaults definition file $defaultfdef $!";
	my $configdef     = undef;
	my $enzymedef     = undef;
	my $modifications = undef;
	while ( my $line = <DEF> )
	{
		chomp;

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
			my @results = split( ' ', $line );
			$PARAMS{ $results[0] } = $results[1];

			#print $_[0], "\t", $PARAMS{ $_[0] }, "\n";
		} elsif ( defined($enzymedef) )
		{
			my @results = split( ' ', $line );
			$ENZ{ $results[0] }->{'name'}     = $results[1];
			$ENZ{ $results[0] }->{'cutAA'}    = $results[3];
			$ENZ{ $results[0] }->{'notcutAA'} = $results[4];
		} elsif ( defined($modifications) )
		{
			my @results = split( ' ', $line );
			$MOD{ $results[0] } = $results[1];
		}
	}
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
			$verbose && print "Fixed Modification defined: ", $key, " ", $MSTAB{$key}->{'native'}, "+ ", $MOD{$key}, "\n";
			$MSTAB{$key}->{'native'} += $MOD{$key};

			#print $MSTAB{$_}->{'native'},"\n";
		}
	}

	#define variable modification X, U, B, J
	if ( $PARAMS{'variable_mod'} )
	{
	#	my ( $AA, $delta ) = split /,|:/, $PARAMS{'variable_mod'};
	#	$MSTAB{'X'}->{'native'} = $MSTAB{$AA}->{'native'} + $delta;
		my @AAlist = ('X', 'U', 'B', 'J');
		my @AAdelta = split /,|:/, $PARAMS{'variable_mod'};
		my $nmods = int( ( scalar @AAdelta ) / 2 );
		for my $i ( 0 .. $nmods - 1 ) {
			my $pseudo = $AAlist[$i];
			$MSTAB{$pseudo}->{'native'} = $MSTAB{$AAdelta[2*$i]}->{'native'} + $AAdelta[2*$i + 1];
		}
		
	}
	if ( $MSTAB{'X'}->{'native'} )
	{
	#	$verbose && print "modificaton X: ", $MSTAB{'X'}->{'native'}, "\n";
		my @AAlist = ('X', 'U', 'B', 'J');
 		my @AAdelta = split /,|:/, $PARAMS{'variable_mod'};
		for my $i (0..($#AAdelta / 2 - 1) ){
			$verbose && print "modification X: ", $MSTAB{$AAdelta[2*$i]}->{'native'}, "\n";
		}
	
	}
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
	unless ( $PARAMS{'usenprescores'} )
	{
		$PARAMS{'usenprescores'} = 100;
	}
	$verbose && ( print "xlink targets: ", $PARAMS{'AArequired'}, "\n" );
	my $dbname;
	if ( -e $PARAMS{'database'} )
	{
		$dbname = $PARAMS{'database'};
	} elsif ( -e File::Spec->catfile( $xquestdir, $PARAMS{'database'} ) )
	{
		$dbname = File::Spec->catfile( $xquestdir, $PARAMS{'database'} );
	} else
	{
		die "cannot open database file $dbname $!";
	}
	## Get the Db path/basename is used for the db indices
	#$dbname =~ s/\.\w+//; #Doesn't work on paths with . in the name
	my ($dbname_file, $dbname_dir) = File::Basename::fileparse($dbname);
	$dbname = $dbname_dir . $dbname_file;
	open WEBCONFIG, "<$webconfig"
	  or warn "could not open web config file $! ignoring";
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

sub printenumeration
{
	my $ENUMERATIONINDEX = shift;
	foreach ( keys %$ENUMERATIONINDEX )
	{
		my @combinations = @{ $ENUMERATIONINDEX->{$_} };
		print "$_ @combinations\n";
	}
}
## @method void printprogress()
# Print the progress of the search to a progressfile
sub printprogress
{
	my $progressfile = shift;    # the filename
	my $message      = shift;    # the message to print

	# overrides the old file, there should always be only one line in the file
	open PROGRESS, ">", "$progressfile" or warn "cannot open progress file $progressfile $!";
	unless (*PROGRESS)
	{
		return;
	}
### MAKE THE PROGRESS FILEHANDLE HOT OTHERWISE THE OUTPUT IS BUFFERED!
	select( ( select(PROGRESS), $| = 1 )[0] );
	print PROGRESS "$message";
	close(PROGRESS);
}

sub printversion
{
	my $version = shift;
	print "$version\n";
	exit 0;
}

sub printversion_changes
{
	my $version = shift;
	print "$version \nfixed bug in ion-tag mode that set the max pepsize to small missing  loop- and mon-links
\n";
	exit 0;
}

sub get_username
{
	my $OS = $^O;    # get the current OS
	my $user_name;
	if ( $OS eq "MSWin32" )
	{
		$user_name = "WindowsUser";    # username hardcoded because on prottools no Win32 package
	}
	if ( $OS eq "linux" )
	{
		$user_name = $ENV{'LOGNAME'};
	}
	chomp($user_name);
	unless ($user_name)
	{
		die "cannot read your username from $!";
	}
	return $user_name;
}

sub usage
{
	print "
	SOFTWARE: ", basename($0), " version $version
	
	AUTHOR: Thomas Walzthoeni based on original version by Oliver Rinner.

	INFORMATION: A software to identify cross-linked peptides from MS/MS spectra.
 
 	USAGE: ", basename($0), " -Option [Parameter]

	REQUIRED OPTIONS [defaults]:
	-def [] xquest definition file (xquest.def)
	-xquestdir [] xquest installation directory
	-list [] *_isotopepairs.txt file, filelist from compare_peaks program, holds the spectra to be searched
	-resdir [] the resultdirectory
	
	OTHER OPTIONS [defaults]:
	-specxml [] a spectrum xml file that holds the spectra
	-nidx [] create new database indicies
	-oidx [] make only index and exit
	-stat [auto] status filename, if not set filename is automatically generated during run
	-progress [auto] progress filename, if not set filename is automatically generated during run
	
	-h print this help
	
	EXAMPLE:
	$0 -nidx -specxml FN-XL_matched.spec.xml -def xquest.def -xquestdir /xquest/root/path -list FN-XL_matched_isotopepairs.txt -resdir FN-XL_matched
	
	
	";
	exit;
}
