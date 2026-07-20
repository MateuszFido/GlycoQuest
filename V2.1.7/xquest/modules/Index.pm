package Index;
use strict;
#---------------------------------------------------------------------------
# Module: Index.pm
# Author(s): Thomas Walzthoeni, xQuest specific modules are based on 
# original work by Oliver Rinner.
# Description: Module for indexing databases.
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

use FindBin;
#use lib "$FindBin::Bin/../../perl5";
use MLDBM qw(DB_File Storable);
use Bio::Perl;
use DB_File;
use Fcntl;
use Xquest_Digest;
use Data::Dumper;

# Berkeley DB 1.x can exhaust its finite pool of overflow pages when dense ion
# bins are stored with the default 4 KiB hash pages.  xQuest ion indexes contain
# large Storable arrays, so use the largest hash page size accepted by the
# bundled Berkeley DB 1.x runtime (64 KiB is rejected by this build).
# Each tie receives a fresh HASHINFO object because target and decoy indexes may
# be open at the same time.
sub ion_hash_options {
	my $hashinfo = DB_File::HASHINFO->new();
	$hashinfo->{'bsize'} = 32768;
	$hashinfo->{'cachesize'} = 64 * 1024 * 1024;
	return $hashinfo;
}

sub new {
	my $class = shift();
	my $self  = {};
	bless $self, $class;
	my $indb             = shift;
	my $basename         = shift;
	my $statusfilehandle = shift;
	my $PARAMS           = shift;
	my $ENZ              = shift;
	my $MSTAB            = shift;
	my $verbose          = shift;
	my $superverbose     = shift;
	my $newindex         = shift;
	my $onlyindex        = shift;
	my $enumeration      = shift;
	my $iontagmode       = shift;
	my $indexinfo        = shift;

	$self->{'iontagmode'}  = $iontagmode;
	$self->{'enumeration'} = $enumeration;

	my $openindex    = undef;
	my $makenewindex = undef;
	
	### define if a new db should be created or another index should be opened #####
	if ( $newindex || $onlyindex ) {
		$makenewindex = defined;
	}
	else {
		$openindex = defined;
	}

	my (
		$i,
		$j,
		$enumerationindexname_inter_xlinks,
		$enumerationindexname_intra_xlinks,
		$enumerationindexname_mono,
		$enumerationindexname_intra
	);
	
	use vars qw(%IONINDEX %PEPINDEX $WEBPARAMS %INFOINDEX %ENUMERATION_INTER_XLINKS %ENUMERATION_INTRA_XLINKS %ENUMERATION_MONO %ENUMERATION_INTRA);
	my $indexinfoname = $basename . "_info.db";
	my $ionindexname  = $basename . "_ion.db";
	my $pepindexname  = $basename . "_peps.db";
	my $infoindex     = $basename . "_info.db";
	my $indexstatus   = $basename . "_index.stat";

##### DEFINE PARAMS FROM THE PARAM ARRAY #########################################
	$enumeration = $PARAMS->{'enumerate'};
	my $decoy = $PARAMS->{'decoy'};

	unless ( $PARAMS->{'Iontag_writeaftern'} ) {
		$PARAMS->{'Iontag_writeaftern'} = int(1000000);
	}

	my $writetodisk_after_nporteins = $PARAMS->{'Iontag_writeaftern'};
	my $matchall                    = $PARAMS->{'matchall'};
	my $massmatchonly               = $PARAMS->{'massmatchonly'};

my $checkcode = join "::", $PARAMS->{'database'}, $PARAMS->{'enzyme_num'},
	  $PARAMS->{'missed_cleavages'},     $PARAMS->{'requiredmissed_cleavages'},
	  $PARAMS->{'tryptic_termini'},      $PARAMS->{'mindigestlength'},
	  $PARAMS->{'maxdigestlength'},      $PARAMS->{'nocutatxlink'},
	  $PARAMS->{'variable_mod'},         $PARAMS->{'nvariable_mod'},
	  $PARAMS->{'ionindexintprecision'}, $PARAMS->{'intprecision'},
	  $PARAMS->{'enumeration_index_mode'}, $PARAMS->{'AArequired'},
	  $PARAMS->{'xkinkerID'}, $PARAMS->{'xlinkermw'}, $PARAMS->{'monolinkmw'},;

	#print "checkcode $checkcode\n";

	if ($enumeration) {
		$enumerationindexname_inter_xlinks = join "_", $basename,"inter_xlink_enum.db";
		$enumerationindexname_intra_xlinks = join "_", $basename,"intra_xlink_enum.db";
		$enumerationindexname_mono  = join "_", $basename, "mono_enum.db";
		$enumerationindexname_intra = join "_", $basename, "intra_enum.db";
	}

	#make new ion and peptide index or open existing
	if ( defined($openindex) || defined($indexinfo) ) {
		print $statusfilehandle "open infoindex $infoindex...checking consistency\n";
		tie( %INFOINDEX, 'MLDBM', $infoindex ) or die "can't open tie to $infoindex: $!";
		#print "true";
		#print $INFOINDEX{'checkcode'}."\n";
		#print $checkcode;
		unless ( $INFOINDEX{'checkcode'} eq $checkcode ) {
		print "parameters have changed since last index build. You should rebuild the index (-option -nidx)\n";
		exit;
		}
		if ($indexinfo) {
			$self->{'INFOINDEX'} = \%INFOINDEX;

			$self->printindexinfo(*STDOUT);
			exit 1;
		}

		print $statusfilehandle "open peptideindex $pepindexname...\n";
		tie( %PEPINDEX, 'MLDBM', $pepindexname ) or die "can't open tie to $pepindexname: $!";

		if ($iontagmode) {
			print $statusfilehandle "open ionindex $ionindexname...\n";
			tie( %IONINDEX, 'MLDBM', $ionindexname, O_RDWR, 0666, ion_hash_options() ) or die "can't open tie to $ionindexname: $!";
		}

		elsif ($enumeration) {

			print $statusfilehandle
			  "open enumerationindexes $enumerationindexname_inter_xlinks...\n";
			tie( %ENUMERATION_INTER_XLINKS, 'MLDBM',
				$enumerationindexname_inter_xlinks )
			  or die "can't open tie to $enumerationindexname_inter_xlinks: $!";

			print $statusfilehandle
			  "open enumerationindexes $enumerationindexname_intra_xlinks...\n";
			tie( %ENUMERATION_INTRA_XLINKS, 'MLDBM',
				$enumerationindexname_intra_xlinks )
			  or die "can't open tie to $enumerationindexname_intra_xlinks: $!";

			print $statusfilehandle
			  "open enumerationindexes $enumerationindexname_mono...\n";
			tie( %ENUMERATION_MONO, 'MLDBM', $enumerationindexname_mono )
			  or die "can't open tie to $enumerationindexname_intra: $!";

			print $statusfilehandle
			  "open enumerationindexes $enumerationindexname_intra...\n";
			tie( %ENUMERATION_INTRA, 'MLDBM', $enumerationindexname_intra )
			  or die "can't open tie to $enumerationindexname_intra: $!";

		}

	}

	else {    #make new indexes delete old if present
		print $statusfilehandle "writing new peptide index $pepindexname\n";
		open STATUS, ">$indexstatus" or die "cannot open index status file $indexstatus $!";
		my @proteinids = ();
		if ( -e $pepindexname ) {
			unlink($pepindexname);
		}

		tie( %PEPINDEX, 'MLDBM', $pepindexname, O_CREAT | O_RDWR, 0666 )
		  or die "can't open tie to $pepindexname: $!";

		print $statusfilehandle "writing infoindex $infoindex...\n";

		if ( -e $infoindex ) {
			unlink($infoindex);
		}
		tie( %INFOINDEX, 'MLDBM', $infoindex, O_CREAT | O_RDWR, 0666 )
		  or die "can't open tie to $infoindex: $!";

		if ($iontagmode) {
			print $statusfilehandle "writing new ion index $ionindexname\n";
			if ( -e $ionindexname ) {
				unlink($ionindexname);
			}
			tie( %IONINDEX, 'MLDBM', $ionindexname, O_CREAT | O_RDWR, 0666, ion_hash_options() )
			  or die "can't open tie to $ionindexname: $!";
			for $i ( $PARAMS->{'minionsize'} ... $PARAMS->{'maxionsize'} ) {
				$IONINDEX{ int($i) } = [];
			}
		}
		elsif ($enumeration) {
			if ( $PARAMS->{'search_intercrosslinks'} ) {
				if ( -e $enumerationindexname_inter_xlinks ) {
					print $statusfilehandle
					  "deletete $enumerationindexname_inter_xlinks...\n";

					unlink($enumerationindexname_inter_xlinks);
				}
				tie(
					%ENUMERATION_INTER_XLINKS, 'MLDBM',
					$enumerationindexname_inter_xlinks,
					O_CREAT | O_RDWR, 0666
				  )
				  or die
				  "can't open tie to $enumerationindexname_inter_xlinks: $!";
				print $statusfilehandle
				  "create $enumerationindexname_inter_xlinks...\n";

			}
			if ( $PARAMS->{'search_intracrosslinks'} ) {
				if ( -e $enumerationindexname_intra_xlinks ) {
					unlink($enumerationindexname_intra_xlinks);
					print $statusfilehandle
					  "deletete $enumerationindexname_intra_xlinks...\n";

				}
				tie(
					%ENUMERATION_INTRA_XLINKS, 'MLDBM',
					$enumerationindexname_intra_xlinks,
					O_CREAT | O_RDWR, 0666
				  )
				  or die
				  "can't open tie to $enumerationindexname_intra_xlinks: $!";
				print $statusfilehandle
				  "create $enumerationindexname_intra_xlinks...\n";
			}
			if ( $PARAMS->{'search_monolinks'} ) {
				if ( -e $enumerationindexname_mono ) {
					unlink($enumerationindexname_mono);
					print $statusfilehandle
					  "deletete $enumerationindexname_mono...\n";

				}

				tie( %ENUMERATION_MONO, 'MLDBM', $enumerationindexname_mono,
					O_CREAT | O_RDWR, 0666 )
				  or die "can't open tie to $enumerationindexname_mono: $!";
				print $statusfilehandle
				  "create $enumerationindexname_mono...\n";
			}
			if ( $PARAMS->{'search_intralinks'} ) {
				if ( -e $enumerationindexname_intra ) {
					unlink($enumerationindexname_intra);
					print $statusfilehandle
					  "deletete $enumerationindexname_mono...\n";

				}
				tie( %ENUMERATION_INTRA, 'MLDBM', $enumerationindexname_intra,
					O_CREAT | O_RDWR, 0666 )
				  or die "can't open tie to $enumerationindexname_intra: $!";
				print $statusfilehandle
				  "create $enumerationindexname_intra...\n";
			}
		}
		my $nproteins = 0;
		my $npeptides = 0;
		my $nions     = 0;
		my ($pepobj);
		my @sequences = ();
		
		my $numpeptideshash={};
		
		while ( my $seq = $indb->next_seq() ) {
			$nproteins++;
			push @sequences, $seq;
#			if ($decoy) {
#				my $desc   = $seq->desc;
#				my $id     = $seq->id;
#				my $revseq = reverse $seq->seq;
#				my $newseq = Bio::Seq->new(
#					-seq  => $revseq,
#					-desc => "reverse of $desc",
#					-id   => "reverse_$id",
#				);
#				$nproteins++;
#				push @sequences, $newseq;
#			}
		}
		my %TEMPIONINDEX = ();

		my $proteincounter = 0;

		foreach my $protseq (@sequences) {
			my $proteinID = $protseq->id;
			push @proteinids, $proteinID;
			$proteincounter++;
			$verbose && print "Indexing ", $protseq->id, " \n";
			
			##### Xquest_digest generates PepObjects
			my $peptides = Xquest_Digest::getpeps( $MSTAB, $PARAMS, $ENZ, $protseq, $superverbose );
			$verbose && print $protseq->id, ": \#peptides matching requirements: ", $#$peptides + 1, "\n";
			
			### Store the number of peptides with the protein ID (is then stored in Infoindex)
			$numpeptideshash->{$protseq->id}=$#$peptides + 1;
			
			foreach $pepobj (@$peptides) {
			$npeptides +=$self->storeIndex( $pepobj, \%TEMPIONINDEX, \%PEPINDEX,$PARAMS, $superverbose, $proteinID );
			}
			$protseq->DESTROY;

			print STATUS "$proteincounter\/$nproteins: indexed $proteinID with ",$#$peptides + 1, " peptides  total unique peptides=$npeptides\n";
			if (   ( $proteincounter % $writetodisk_after_nporteins == 0 )	|| ( $proteincounter >= $nproteins ) )
			{
				if ($iontagmode) {
					my $ionentries = 0;
					print STATUS "-> writing to ionidex ... ";
					foreach my $ionmz ( keys %TEMPIONINDEX ) {
						#print "$ionmz\n";
						my $ionlist = $IONINDEX{ $ionmz }; #needs to be done that weired for complex hashes in MLDBM
						foreach my $peptid ( @{ $TEMPIONINDEX{$ionmz} } ) {
							push @$ionlist, $peptid;
							$ionentries++;
							$nions++;
						}
						$IONINDEX{$ionmz} = $ionlist;
					}
					print STATUS "... wrote $ionentries ion-peptide associations to index\n";
				}
				%TEMPIONINDEX = ();
			}

		}
		#print Dumper(\%PEPINDEX);
		#exit;
		print STATUS "-> done with peptide index ... \n";
		#print "nions: $nions\n";
		$INFOINDEX{'nions'}     = $nions;
		$INFOINDEX{'nprots'}    = $nproteins;
		$INFOINDEX{'npeps'}     = $npeptides;
		$INFOINDEX{'checkcode'} = $checkcode;
		$INFOINDEX{'protids'}   = \@proteinids;
		$INFOINDEX{'npeptideshash'}   = $numpeptideshash;
		
		if ($enumeration) {
			my ( $ncombinations, $nmonolinks, $nintralinks,
				$intraxlinkcombinations );
			if ( $PARAMS->{'enumeration_index_mode'} =~ /bigmem/ ) {
				print STATUS "making enumeration index in large mem mode\n";
				( $ncombinations, $nmonolinks, $nintralinks ) =
				  PepObj::storecombinationIndex_bigmem( \%PEPINDEX, $PARAMS,
					\%ENUMERATION_INTER_XLINKS, \%ENUMERATION_MONO,
					\%ENUMERATION_INTRA, $verbose );
			}
			elsif ( $PARAMS->{'enumeration_index_mode'} =~ /smarthash/ ) {
				print STATUS "making enumeration index in  smarthash mode\n";
				(
					$ncombinations, $intraxlinkcombinations, $nmonolinks,
					$nintralinks
				  )
				  = Index::storecombinationIndex_smarthash(
					\%PEPINDEX,                 $PARAMS,
					\%ENUMERATION_INTER_XLINKS, \%ENUMERATION_INTRA_XLINKS,
					\%ENUMERATION_MONO,         \%ENUMERATION_INTRA,
					\%INFOINDEX,                $verbose,
					*STATUS
				  );
			}
			elsif ( $PARAMS->{'enumeration_index_mode'} =~ /split_db/ ) {
				print STATUS "making enumeration index in split DB mode\n";
				( $ncombinations, $nmonolinks, $nintralinks ) =
				  Index::storecombinationIndex_smarthash_splitdb( \%PEPINDEX,
					$PARAMS, \%ENUMERATION_INTER_XLINKS, \%ENUMERATION_MONO,
					\%ENUMERATION_INTRA, $verbose, $basename );
			}

			else {
				print STATUS "making enumeration index in low mem mode\n";
				( $ncombinations, $nmonolinks, $nintralinks ) =
				  Index::storecombinationIndex( \%PEPINDEX, $PARAMS,
					\%ENUMERATION_INTER_XLINKS, \%ENUMERATION_MONO,
					\%ENUMERATION_INTRA, $verbose );

			}

			$INFOINDEX{'nmonolinks'}   = $nmonolinks;
			$INFOINDEX{'nintralinks'}  = $nintralinks;
			$INFOINDEX{'nintraxlinks'} = $intraxlinkcombinations;
			$INFOINDEX{'ninterxlinks'} = $ncombinations;

print "\nindexed ", $ncombinations, "xlink-combinations $intraxlinkcombinations intra-xlinkcombinations $nmonolinks monolinks and $nintralinks intralinks\n";

		}

		#		print $statusfilehandle "indexed ", $self->getnprots,
		#		  " proteins\nindexed ", $self->getnpeps, " peptides\nindexed ",
		#		  $self->getnions, " ions\n";
		$onlyindex && exit 0;

	}
	$self->{'INFOINDEX'} = \%INFOINDEX;

	$self->{'IONINDEX'} = \%IONINDEX;
	$self->{'PEPINDEX'} = \%PEPINDEX;

	$self->{'ENUMERATION_INTER_XLINKS'} = \%ENUMERATION_INTER_XLINKS;
	$self->{'ENUMERATION_INTRA_XLINKS'} = \%ENUMERATION_INTRA_XLINKS;
	$self->{'ENUMERATION_MONO'}         = \%ENUMERATION_MONO;
	$self->{'ENUMERATION_INTRA'}        = \%ENUMERATION_INTRA;

	$self->printindexinfo($statusfilehandle);


	if ( $PARAMS->{'printpeptides'} ) {

		my $pepfilename = $basename . "_peptides.txt";
		print "writing pepties into $pepfilename<br>";
		open PEPTIDES, ">$pepfilename" or die "cannot open peptidelistfile $pepfilename for writing$!";
		foreach my $pepid ( keys %PEPINDEX ) {
			( my $coor = $pepid ) =~ s/(^\w+::)//;
			print PEPTIDES ( join ",", @{ $PEPINDEX{$pepid}->{'proteinID'} } ),
			  "\t", $PEPINDEX{$pepid}->{'desc'}, "\t",
			  $PEPINDEX{$pepid}->{'seq'},"\t", 
			  $PEPINDEX{$pepid}->{'mw'}, 
			  "\n";
		}
	}
	close(STATUS);
	return $self;
}

sub new_dc {
	my $class = shift();
	my $self  = {};
	bless $self, $class;
	my $indb             = shift;
	my $basename         = shift;
	my $statusfilehandle = shift;
	my $PARAMS           = shift;
	my $ENZ              = shift;
	my $MSTAB            = shift;
	my $verbose          = shift;
	my $superverbose     = shift;
	my $newindex         = shift;
	my $onlyindex        = shift;
	my $enumeration      = shift;
	my $iontagmode       = shift;
	my $indexinfo        = shift;

	$self->{'iontagmode'}  = $iontagmode;
	$self->{'enumeration'} = $enumeration;

	my $openindex    = undef;
	my $makenewindex = undef;
	
	### define if a new db should be created or another index should be opened #####
	if ( $newindex || $onlyindex ) {
		$makenewindex = defined;
	}
	else {
		$openindex = defined;
	}

	my (
		$i,
		$j,
		$enumerationindexname_inter_xlinks,
		$enumerationindexname_intra_xlinks,
		$enumerationindexname_mono,
		$enumerationindexname_intra
	);
	
	use vars qw(%IONINDEX_DC %PEPINDEX_DC $WEBPARAMS %INFOINDEX_DC %ENUMERATION_INTER_XLINKS %ENUMERATION_INTRA_XLINKS %ENUMERATION_MONO %ENUMERATION_INTRA);
	my $indexinfoname = $basename . "_info.db";
	my $ionindexname  = $basename . "_ion.db";
	my $pepindexname  = $basename . "_peps.db";
	my $infoindex     = $basename . "_info.db";
	my $indexstatus   = $basename . "_index.stat";

##### DEFINE PARAMS FROM THE PARAM ARRAY #########################################
	$enumeration = $PARAMS->{'enumerate'};
	my $decoy = $PARAMS->{'decoy'};

	unless ( $PARAMS->{'Iontag_writeaftern'} ) {
		$PARAMS->{'Iontag_writeaftern'} = int(1000000);
	}

	my $writetodisk_after_nporteins = $PARAMS->{'Iontag_writeaftern'};
	my $matchall                    = $PARAMS->{'matchall'};
	my $massmatchonly               = $PARAMS->{'massmatchonly'};

my $checkcode = join "::", $PARAMS->{'database_db'}, $PARAMS->{'enzyme_num'},
	  $PARAMS->{'missed_cleavages'},     $PARAMS->{'requiredmissed_cleavages'},
	  $PARAMS->{'tryptic_termini'},      $PARAMS->{'mindigestlength'},
	  $PARAMS->{'maxdigestlength'},      $PARAMS->{'nocutatxlink'},
	  $PARAMS->{'variable_mod'},         $PARAMS->{'nvariable_mod'},
	  $PARAMS->{'ionindexintprecision'}, $PARAMS->{'intprecision'},
	  $PARAMS->{'enumeration_index_mode'}, $PARAMS->{'AArequired'},
	  $PARAMS->{'xkinkerID'}, $PARAMS->{'xlinkermw'}, $PARAMS->{'monolinkmw'},;

	#print "checkcode $checkcode\n";

	if ($enumeration) {
		$enumerationindexname_inter_xlinks = join "_", $basename,"inter_xlink_enum.db";
		$enumerationindexname_intra_xlinks = join "_", $basename,"intra_xlink_enum.db";
		$enumerationindexname_mono  = join "_", $basename, "mono_enum.db";
		$enumerationindexname_intra = join "_", $basename, "intra_enum.db";
	}

	#make new ion and peptide index or open existing
	if ( defined($openindex) || defined($indexinfo) ) {
		print $statusfilehandle "open infoindex $infoindex...checking consistency\n";
		tie( %INFOINDEX_DC, 'MLDBM', $infoindex ) or die "can't open tie to $infoindex: $!";
		#print "true";
		#print $INFOINDEX{'checkcode'}."\n";
		#print $checkcode;
		unless ( $INFOINDEX_DC{'checkcode'} eq $checkcode ) {
		print "parameters have changed since last index build. You should rebuild the index (-option -nidx)\n";
		exit;
		}
		if ($indexinfo) {
			$self->{'INFOINDEX_DC'} = \%INFOINDEX_DC;

			$self->printindexinfo(*STDOUT);
			exit 1;
		}

		print $statusfilehandle "open peptideindex $pepindexname...\n";
		tie( %PEPINDEX_DC, 'MLDBM', $pepindexname ) or die "can't open tie to $pepindexname: $!";

		if ($iontagmode) {
			print $statusfilehandle "open ionindex $ionindexname...\n";
			tie( %IONINDEX_DC, 'MLDBM', $ionindexname, O_RDWR, 0666, ion_hash_options() ) or die "can't open tie to $ionindexname: $!";
		}

		elsif ($enumeration) {

			print $statusfilehandle
			  "open enumerationindexes $enumerationindexname_inter_xlinks...\n";
			tie( %ENUMERATION_INTER_XLINKS, 'MLDBM',
				$enumerationindexname_inter_xlinks )
			  or die "can't open tie to $enumerationindexname_inter_xlinks: $!";

			print $statusfilehandle
			  "open enumerationindexes $enumerationindexname_intra_xlinks...\n";
			tie( %ENUMERATION_INTRA_XLINKS, 'MLDBM',
				$enumerationindexname_intra_xlinks )
			  or die "can't open tie to $enumerationindexname_intra_xlinks: $!";

			print $statusfilehandle
			  "open enumerationindexes $enumerationindexname_mono...\n";
			tie( %ENUMERATION_MONO, 'MLDBM', $enumerationindexname_mono )
			  or die "can't open tie to $enumerationindexname_intra: $!";

			print $statusfilehandle
			  "open enumerationindexes $enumerationindexname_intra...\n";
			tie( %ENUMERATION_INTRA, 'MLDBM', $enumerationindexname_intra )
			  or die "can't open tie to $enumerationindexname_intra: $!";

		}

	}

	else {    #make new indexes delete old if present
		print $statusfilehandle "writing new peptide index $pepindexname\n";
		open STATUS, ">$indexstatus" or die "cannot open index status file $indexstatus $!";
		my @proteinids = ();
		
		my $numpeptideshash={};
		
		if ( -e $pepindexname ) {
			unlink($pepindexname);
		}

		tie( %PEPINDEX_DC, 'MLDBM', $pepindexname, O_CREAT | O_RDWR, 0666 )
		  or die "can't open tie to $pepindexname: $!";

		print $statusfilehandle "writing infoindex $infoindex...\n";

		if ( -e $infoindex ) {
			unlink($infoindex);
		}
		tie( %INFOINDEX_DC, 'MLDBM', $infoindex, O_CREAT | O_RDWR, 0666 )
		  or die "can't open tie to $infoindex: $!";

		if ($iontagmode) {
			print $statusfilehandle "writing new ion index $ionindexname\n";
			if ( -e $ionindexname ) {
				unlink($ionindexname);
			}
			tie( %IONINDEX_DC, 'MLDBM', $ionindexname, O_CREAT | O_RDWR, 0666, ion_hash_options() )
			  or die "can't open tie to $ionindexname: $!";
			for $i ( $PARAMS->{'minionsize'} ... $PARAMS->{'maxionsize'} ) {
				$IONINDEX_DC{ int($i) } = [];
			}
		}
		elsif ($enumeration) {
			if ( $PARAMS->{'search_intercrosslinks'} ) {
				if ( -e $enumerationindexname_inter_xlinks ) {
					print $statusfilehandle
					  "deletete $enumerationindexname_inter_xlinks...\n";

					unlink($enumerationindexname_inter_xlinks);
				}
				tie(
					%ENUMERATION_INTER_XLINKS, 'MLDBM',
					$enumerationindexname_inter_xlinks,
					O_CREAT | O_RDWR, 0666
				  )
				  or die
				  "can't open tie to $enumerationindexname_inter_xlinks: $!";
				print $statusfilehandle
				  "create $enumerationindexname_inter_xlinks...\n";

			}
			if ( $PARAMS->{'search_intracrosslinks'} ) {
				if ( -e $enumerationindexname_intra_xlinks ) {
					unlink($enumerationindexname_intra_xlinks);
					print $statusfilehandle
					  "deletete $enumerationindexname_intra_xlinks...\n";

				}
				tie(
					%ENUMERATION_INTRA_XLINKS, 'MLDBM',
					$enumerationindexname_intra_xlinks,
					O_CREAT | O_RDWR, 0666
				  )
				  or die
				  "can't open tie to $enumerationindexname_intra_xlinks: $!";
				print $statusfilehandle
				  "create $enumerationindexname_intra_xlinks...\n";
			}
			if ( $PARAMS->{'search_monolinks'} ) {
				if ( -e $enumerationindexname_mono ) {
					unlink($enumerationindexname_mono);
					print $statusfilehandle
					  "deletete $enumerationindexname_mono...\n";

				}

				tie( %ENUMERATION_MONO, 'MLDBM', $enumerationindexname_mono,
					O_CREAT | O_RDWR, 0666 )
				  or die "can't open tie to $enumerationindexname_mono: $!";
				print $statusfilehandle
				  "create $enumerationindexname_mono...\n";
			}
			if ( $PARAMS->{'search_intralinks'} ) {
				if ( -e $enumerationindexname_intra ) {
					unlink($enumerationindexname_intra);
					print $statusfilehandle
					  "deletete $enumerationindexname_mono...\n";

				}
				tie( %ENUMERATION_INTRA, 'MLDBM', $enumerationindexname_intra,
					O_CREAT | O_RDWR, 0666 )
				  or die "can't open tie to $enumerationindexname_intra: $!";
				print $statusfilehandle
				  "create $enumerationindexname_intra...\n";
			}
		}
		my $nproteins = 0;
		my $npeptides = 0;
		my $nions     = 0;
		my ($pepobj);
		my @sequences = ();

		while ( my $seq = $indb->next_seq() ) {
			$nproteins++;
			push @sequences, $seq;
			
			if ($decoy) {
				my $desc   = $seq->desc;
				my $id     = $seq->id;
				my $revseq = reverse $seq->seq;
				my $newseq = Bio::Seq->new(
					-seq  => $revseq,
					-desc => "reverse of $desc",
					-id   => "reverse_$id",
				);
				$nproteins++;
				push @sequences, $newseq;
			}
		}
		my %TEMPIONINDEX_DC = ();

		my $proteincounter = 0;

		foreach my $protseq (@sequences) {
			my $proteinID = $protseq->id;
			push @proteinids, $proteinID;
			$proteincounter++;

			$verbose && print "indexing ", $protseq->id, " \n";
			##### Xquest_digest generates PepObjects
			my $peptides = Xquest_Digest::getpeps( $MSTAB, $PARAMS, $ENZ, $protseq, $superverbose );
			$verbose && print $protseq->id, ": \#peptides matching requirements: ", $#$peptides + 1, "\n";
			### Store the number of peptides with the protein ID (is then stored in Infoindex)
			$numpeptideshash->{$protseq->id}=$#$peptides + 1;
			foreach $pepobj (@$peptides) {
			$npeptides +=$self->storeIndex( $pepobj, \%TEMPIONINDEX_DC, \%PEPINDEX_DC,$PARAMS, $superverbose, $proteinID );
			}
			$protseq->DESTROY;
#			print Dumper(\%TEMPIONINDEX);
#			print Dumper(\%PEPINDEX);
#			exit;
			
			print STATUS "$proteincounter\/$nproteins: indexed $proteinID with ",$#$peptides + 1, " peptides  total unique peptides=$npeptides\n";
			if (   ( $proteincounter % $writetodisk_after_nporteins == 0 )	|| ( $proteincounter >= $nproteins ) )
			{

				if ($iontagmode) {
					my $ionentries = 0;
					print STATUS "-> writing to ionidex ... ";
					foreach my $ionmz ( keys %TEMPIONINDEX_DC ) {
						#print "$ionmz\n";
						my $ionlist = $IONINDEX_DC{ $ionmz }; #needs to be done that weired for complex hashes in MLDBM
						foreach my $peptid ( @{ $TEMPIONINDEX_DC{$ionmz} } ) {
							push @$ionlist, $peptid;
							$ionentries++;
							$nions++;
						}
						$IONINDEX_DC{$ionmz} = $ionlist;
					}
					print STATUS "... wrote $ionentries ion-peptide associations to index\n";
				}
				%TEMPIONINDEX_DC = ();
			}

		}
		#print Dumper(\%PEPINDEX);
		#exit;
		print STATUS "-> done with peptide index ... \n";

		#print "nions: $nions\n";
		$INFOINDEX_DC{'nions'}     = $nions;
		$INFOINDEX_DC{'nprots'}    = $nproteins;
		$INFOINDEX_DC{'npeps'}     = $npeptides;
		$INFOINDEX_DC{'checkcode'} = $checkcode;
		$INFOINDEX_DC{'protids'}   = \@proteinids;
		$INFOINDEX_DC{'npeptideshash'}   = $numpeptideshash;
		
		if ($enumeration) {
			my ( $ncombinations, $nmonolinks, $nintralinks,
				$intraxlinkcombinations );
			if ( $PARAMS->{'enumeration_index_mode'} =~ /bigmem/ ) {
				print STATUS "making enumeration index in large mem mode\n";
				( $ncombinations, $nmonolinks, $nintralinks ) =
				  PepObj::storecombinationIndex_bigmem( \%PEPINDEX_DC, $PARAMS,
					\%ENUMERATION_INTER_XLINKS, \%ENUMERATION_MONO,
					\%ENUMERATION_INTRA, $verbose );
			}
			elsif ( $PARAMS->{'enumeration_index_mode'} =~ /smarthash/ ) {
				print STATUS "making enumeration index in  smarthash mode\n";
				(
					$ncombinations, $intraxlinkcombinations, $nmonolinks,
					$nintralinks
				  )
				  = Index::storecombinationIndex_smarthash(
					\%PEPINDEX_DC,                 $PARAMS,
					\%ENUMERATION_INTER_XLINKS, \%ENUMERATION_INTRA_XLINKS,
					\%ENUMERATION_MONO,         \%ENUMERATION_INTRA,
					\%INFOINDEX,                $verbose,
					*STATUS
				  );
			}
			elsif ( $PARAMS->{'enumeration_index_mode'} =~ /split_db/ ) {
				print STATUS "making enumeration index in split DB mode\n";
				( $ncombinations, $nmonolinks, $nintralinks ) =
				  Index::storecombinationIndex_smarthash_splitdb( \%PEPINDEX_DC,
					$PARAMS, \%ENUMERATION_INTER_XLINKS, \%ENUMERATION_MONO,
					\%ENUMERATION_INTRA, $verbose, $basename );
			}

			else {
				print STATUS "making enumeration index in low mem mode\n";
				( $ncombinations, $nmonolinks, $nintralinks ) =
				  Index::storecombinationIndex( \%PEPINDEX_DC, $PARAMS,
					\%ENUMERATION_INTER_XLINKS, \%ENUMERATION_MONO,
					\%ENUMERATION_INTRA, $verbose );

			}

			$INFOINDEX_DC{'nmonolinks'}   = $nmonolinks;
			$INFOINDEX_DC{'nintralinks'}  = $nintralinks;
			$INFOINDEX_DC{'nintraxlinks'} = $intraxlinkcombinations;
			$INFOINDEX_DC{'ninterxlinks'} = $ncombinations;

#print $statusfilehandle "\nindexed ", $ncombinations, xlink-combinations $intraxlinkcombinations intra-xlinkcombinations $nmonolinks monolinks and $nintralinks intralinks\n";

		}

		#		print $statusfilehandle "indexed ", $self->getnprots,
		#		  " proteins\nindexed ", $self->getnpeps, " peptides\nindexed ",
		#		  $self->getnions, " ions\n";
		$onlyindex && exit 0;

	}
	$self->{'INFOINDEX_DC'} = \%INFOINDEX_DC;

	$self->{'IONINDEX_DC'} = \%IONINDEX_DC;
	$self->{'PEPINDEX_DC'} = \%PEPINDEX_DC;

	$self->{'ENUMERATION_INTER_XLINKS_DC'} = \%ENUMERATION_INTER_XLINKS;
	$self->{'ENUMERATION_INTRA_XLINKS_DC'} = \%ENUMERATION_INTRA_XLINKS;
	$self->{'ENUMERATION_MONO_DC'}         = \%ENUMERATION_MONO;
	$self->{'ENUMERATION_INTRA_DC'}        = \%ENUMERATION_INTRA;

	#$self->printindexinfo($statusfilehandle);

#print Dumper (\%PEPINDEX_DC);

	if ( $PARAMS->{'printpeptides'} ) {

		my $pepfilename = $basename . "_peptides.txt";
		print "writing pepties into $pepfilename<br>";
		open PEPTIDES, ">$pepfilename" or die "cannot open peptidelistfile $pepfilename for writing$!";
		foreach my $pepid ( keys %PEPINDEX_DC ) {
			( my $coor = $pepid ) =~ s/(^\w+::)//;
			#print $pepid."\n";
			#print @{ $PEPINDEX_DC{$pepid}->{'proteinID'} };
			print PEPTIDES ( join ",", @{ $PEPINDEX_DC{$pepid}->{'proteinID'} } ),"\t", $PEPINDEX_DC{$pepid}->{'desc'}, "\t",
			  $PEPINDEX_DC{$pepid}->{'seq'},"\t", 
			  $PEPINDEX_DC{$pepid}->{'mw'}, 
			  "\n";
		}
	}
	close(STATUS);
	return $self;
}


sub getnprots {
	my $self = shift;
	return $self->get_infoindex->{'nprots'};
}

sub getnpeps {
	my $self = shift;
	return $self->get_infoindex->{'npeps'};
}

sub getnions {
	my $self = shift;
	return $self->get_infoindex->{'nions'};

}

sub getnmonolinks {
	my $self = shift;
	if ( defined( $self->get_infoindex->{'nmonolinks'} ) ) {
		return $self->get_infoindex->{'nmonolinks'};
	}
	else {
		return 0;
	}
}

sub getnintralinks {
	my $self = shift;
	if ( defined( $self->get_infoindex->{'nintralinks'} ) ) {
		return $self->get_infoindex->{'nintralinks'};
	}
	else {
		return 0;
	}
}

sub getnintraxlinks {
	my $self = shift;
	if ( defined( $self->get_infoindex->{'nintraxlinks'} ) ) {
		return $self->get_infoindex->{'nintraxlinks'};
	}
	else {
		return 0;
	}
}

sub getninterxlinks {
	my $self = shift;
	if ( defined( $self->get_infoindex->{'ninterxlinks'} ) ) {
		return $self->get_infoindex->{'ninterxlinks'};
	}
	else {
		return 0;
	}
}

sub get_infoindex {
	my $self = shift;
	return $self->{'INFOINDEX'};
}

sub get_infoindex_dc {
	my $self = shift;
	return $self->{'INFOINDEX_DC'};
}

sub get_ionindex {
	my $self = shift;
	return $self->{'IONINDEX'};
}

sub get_ionindex_dc {
	my $self = shift;
	return $self->{'IONINDEX_DC'};
}

sub get_pepindex {
	my $self = shift;
	return $self->{'PEPINDEX'};
}

sub get_pepindex_dc {
	my $self = shift;
	return $self->{'PEPINDEX_DC'};
}

sub get_enumeration_inter_xlinks_index {
	my $self = shift;
	return $self->{'ENUMERATION_INTER_XLINKS'};
}

sub get_enumeration_intra_xlinks_index {
	my $self = shift;
	return $self->{'ENUMERATION_INTRA_XLINKS'};
}

sub get_enumeration_intralinks_index {
	my $self = shift;
	return $self->{'ENUMERATION_INTRA'};
}

sub get_enumeration_monolinks_index {
	my $self = shift;
	return $self->{'ENUMERATION_MONO'};
}

sub printindexinfo {
	my $self       = shift;
	my $filehandle = shift;

	print $filehandle "number of proteins: ", $self->getnprots, "\n";
	my $nrpeptides = $self->getnpeps;
	print "number of proteins: ", $self->getnprots, "<br>\n";

	my $nrpeptides = $self->getnpeps;

	print $filehandle "number of peptides: ", $nrpeptides, "\n";
	print "number of peptides: ", $nrpeptides, "<br>\n";

	print $filehandle "number of theoretical x-link combinations (n^2/2 + n/2): ",
	  ( $nrpeptides**2 + $nrpeptides ) / 2, "\n";
	print "number of theoretical x-link combinations (n^2/2 + n/2): ",
	  ( $nrpeptides**2 + $nrpeptides) / 2, "<br>\n";

	if ( $self->{'iontagmode'} ) {
		print $filehandle "number of indexed ions: ", $self->getnions, "\n";
		print "number of indexed ions: ", $self->getnions, "<br>\n";
	}
	elsif ( $self->{'enumeration'} ) {
		print $filehandle "number of monolinks: ", $self->getnmonolinks,
		  "<br>\n";
		print $filehandle "number of intralinks: ", $self->getnintralinks,
		  "<br>\n";
		print $filehandle "number of intra-protein cross-links: ",
		  $self->getnintraxlinks, "<br>\n";
		print $filehandle "number of inter-protein cross-links: ",
		  $self->getninterxlinks, "<br>\n";
	}

}

sub storecombinationIndex {
	my $PEPINDEX          = shift;
	my $PARAMS            = shift;
	my $ENUMERATION       = shift;
	my $ENUMERATION_MONO  = shift;
	my $ENUMERATION_INTRA = shift;
	my $verbose           = shift;
	my $xlinkermw         = $PARAMS->{'xlinkermw'};
	my $ms1intprecision   = $PARAMS->{'intprecision'};
	my ( $i, $j, %masshash );
	my $ncombinations = 0;
	my $nmonolinks    = 0;
	my $nintralinks   = 0;
	my @allids        = keys %$PEPINDEX;

	my $id2numberhash = $ENUMERATION->{'id2num'};
	for $i ( 0 .. $#allids ) {
		$id2numberhash->{$i} = $allids[$i];
		$masshash{$i} = $PEPINDEX->{ $allids[$i] }->{'mw'};
	}
	$ENUMERATION->{'id2num'}        = $id2numberhash;
	$ENUMERATION->{'indexprecision'} = $PARAMS->{'intprecision'};

	for $i ( 0 .. $#allids ) {

		#		print "searching combinations for ", $id2numberhash->{$i },
		#		  "\n",;
		my %bufferhash;
		for $j ( $i .. $#allids ) {
			my $newid = [
				$j, $i,
				sprintf( "%.5f",
					( $masshash{$i} + $masshash{$j} + $xlinkermw ) )
			];

			#my $newid = [ $i, $j ];

			#			my $xlinkmass =
			#			  $PEPINDEX->{ $allids[$i] }->{'mw'} +
			#			  $PEPINDEX->{ $allids[$j] }->{'mw'} + $PARAMS->{'xlinkermw'};

			my $xlinkmass = $masshash{$i} + $masshash{$j} + $xlinkermw;
			my $intmass   = int( $xlinkmass * $ms1intprecision );

			push @{ $bufferhash{$intmass} }, $newid;
			$ncombinations++;
		}

		#transfer from temporary buffer into DB hash
		foreach my $intmass ( keys %bufferhash ) {
			my $masscombinations = $ENUMERATION->{$intmass};
			push @{$masscombinations}, @{ $bufferhash{$intmass} };
			$ENUMERATION->{$intmass} = $masscombinations;
		}

		#enumerate monolinks

		my $monolinkmass = $masshash{$i} + $PARAMS->{'monolinkmw'};

		#		  $PEPINDEX->{ $allids[$i] }->{'mw'} + $PARAMS->{'monolinkmw'};
		my $intmass = int( $monolinkmass * $ms1intprecision );

		#get index entry
		my $masscombinations = $ENUMERATION_MONO->{$intmass};
		push @{$masscombinations}, $i;

		#writeindexentry
		$ENUMERATION_MONO->{$intmass} = $masscombinations;
		$nmonolinks++;

		#enumerate intralinks
		my $intralinkmass = $masshash{$i} + $PARAMS->{'xlinkermw'};

		#$PEPINDEX->{ $allids[$i] }->{'mw'} + $PARAMS->{'xlinkermw'};
		my $intmass = int( $intralinkmass * $ms1intprecision );

		#get index entry
		my $masscombinations = $ENUMERATION_INTRA->{$intmass};
		push @{$masscombinations}, $i;

		#writeindexentry
		$ENUMERATION_INTRA->{$intmass} = $masscombinations;
		$nintralinks++;
	}
	return ( $ncombinations, $nmonolinks, $nintralinks );
}

sub storecombinationIndex_smarthash_massbins {
	my $PEPINDEX          = shift;
	my $PARAMS            = shift;
	my $ENUMERATION       = shift;
	my $ENUMERATION_MONO  = shift;
	my $ENUMERATION_INTRA = shift;
	my $verbose           = shift;
	my $xlinkermw         = $PARAMS->{'xlinkermw'};
	my $ms1intprecision   = $PARAMS->{'intprecision'};
	my $peptideminsize    = $PARAMS->{'minpepmr'};
	my $tolerance         = $PARAMS->{'ms1tolerance'};
	my $massbins          = 1;

	if ( $PARAMS->{'tolerancemeasure'} =~ /^ppm/i ) {
		$tolerance =
		  $PARAMS->{'ms1tolerance'} * 1e-6 *
		  1000;    #ppm to amu measure for a 1000 Da peptide
	}
	my $nbins = int( ( $ms1intprecision * $tolerance ) * $massbins + 1 );

	#bins for integer mass within tolerance
	print "nbins tolerance $nbins\n";
	my ( $i, $j, %masshash );
	my $ncombinations = 0;
	my $nmonolinks    = 0;
	my $nintralinks   = 0;
	print "sorting peptides for size ";
	my @allids =
	  sort { $PEPINDEX->{$a}->{'mw'} <=> $PEPINDEX->{$b}->{'mw'} }
	  keys %$PEPINDEX;
	print "... done\n";

	my $number2idhash = $ENUMERATION->{'id2num'};
	for $i ( 0 .. $#allids ) {
		$number2idhash->{$i} = $allids[$i];
	}
	$ENUMERATION->{'id2num'}        = $number2idhash;
	$ENUMERATION->{'indexprecision'} = $ms1intprecision;

	for $i ( 0 .. $#allids ) {
		$masshash{$i} = $PEPINDEX->{ $allids[$i] }->{'mw'};

		#print "$i ", $masshash{$i}, "\n";
	}

	my ( %pepindexmassbins, $id1 );
	for $id1 ( 0 .. $#allids ) {
		push @{ $pepindexmassbins{
				int( $masshash{$id1} * $ms1intprecision * $massbins ) } }, $id1;
	}

	my $xlinkminsize =
	  int( $ms1intprecision * ( 2 * $PARAMS->{'minpepmr'} + $xlinkermw ) -
		  $ms1intprecision );
	my $xlinkmaxsize =
	  int( $ms1intprecision * ( 2 * $PARAMS->{'maxpepmr'} + $xlinkermw ) +
		  $ms1intprecision );

	my $targetxlinksize;
	my $sizerange    = $xlinkmaxsize - $xlinkminsize;
	my $percentrange = int( $sizerange / 100 );
	my $counter      = 0;

	for $targetxlinksize ( $xlinkminsize .. $xlinkmaxsize ) {
		$verbose && print "############# targetxlinksize= ",
		  $targetxlinksize / $ms1intprecision, " ##########\n";
		$counter++;
		my @possible_masscombinations = ();
		for $i ( 0 .. $#allids ) {

			#			if ( ( ( $masshash{$i} + $xlinkermw ) * $ms1intprecision ) >
			#				( $targetxlinksize - $peptideminsize*$ms1intprecision ) )
			if (
				( $masshash{$i} * $ms1intprecision ) >= (
					$ms1intprecision + $targetxlinksize - $xlinkermw *
					  $ms1intprecision
				) / 2
			  )
			{
				next;
			}
			my $target_pepsize = int(
				$massbins * (
					$targetxlinksize - $ms1intprecision *
					  ( $masshash{$i} + $xlinkermw )
				)
			);

			my $binmr;
			for $binmr (
				int( $target_pepsize - $nbins ) ..
				int( $target_pepsize + $nbins ) )
			{
				if ( $pepindexmassbins{$binmr} ) {

					my @pepsinrange = @{ $pepindexmassbins{$binmr} };
					foreach my $pep (@pepsinrange) {
						$ncombinations++;
						my $newid = [
							$pep, $i,
							sprintf(
								"%.5f",
								(
									$masshash{$i} + $masshash{$pep} + $xlinkermw
								)
							)
						];

		#	print sprintf("%.4f",($masshash{$i}+$masshash{$pep}+$xlinkermw)),"\n";
						push @possible_masscombinations, $newid;

						if ($verbose) {
							print " selected pep $i ", $number2idhash->{$i},
							  " ", $masshash{$i}, " targetpepsize= ",
							  $target_pepsize, " found fit: ";
							print "number: @pepsinrange ",
							  $number2idhash->{$pep}, " $pep ", $masshash{$pep},
							  " bin: $binmr final mass ",
							  sprintf(
								"%.5f ",
								(
									$masshash{$i} + $masshash{$pep} + $xlinkermw
								)
							  );
						}
					}

					$verbose && print " \n";
				}
			}
		}
		if ( scalar(@possible_masscombinations) > 0 ) {
			my $masscombinations = $ENUMERATION->{$targetxlinksize};
			push @{$masscombinations}, @possible_masscombinations;
			$ENUMERATION->{$targetxlinksize} = $masscombinations;
		}

	}

	for $i ( 0 .. $#allids ) {

		#enumerate monolinks

		my $monolinkmass =
		  $PEPINDEX->{ $allids[$i] }->{'mw'} + $PARAMS->{'monolinkmw'};
		my $intmass = int( $monolinkmass * $ms1intprecision );

		#get index entry
		my $masscombinations = $ENUMERATION_MONO->{$intmass};
		push @{$masscombinations}, $i;

		#writeindexentry
		$ENUMERATION_MONO->{$intmass} = $masscombinations;
		$nmonolinks++;

		#enumerate intralinks
		my $intralinkmass =
		  $PEPINDEX->{ $allids[$i] }->{'mw'} + $PARAMS->{'xlinkermw'};
		my $intmass = int( $intralinkmass * $ms1intprecision );

		#get index entry
		my $masscombinations = $ENUMERATION_INTRA->{$intmass};
		push @{$masscombinations}, $i;

		#writeindexentry
		$ENUMERATION_INTRA->{$intmass} = $masscombinations;
		$nintralinks++;
	}
	return ( $ncombinations, $nmonolinks, $nintralinks );
}

sub storecombinationIndex_bigmem {
	my $PEPINDEX          = shift;
	my $PARAMS            = shift;
	my $ENUMERATION       = shift;
	my $ENUMERATION_MONO  = shift;
	my $ENUMERATION_INTRA = shift;
	my $verbose           = shift;

	my $xlinkermw       = $PARAMS->{'xlinkermw'};
	my $ms1intprecision = $PARAMS->{'intprecision'};
	my ( $i, $j, %masshash );
	my $ncombinations = 0;
	my $nmonolinks    = 0;
	my $nintralinks   = 0;
	my @allids        = keys %$PEPINDEX;

	my $id2numberhash = $ENUMERATION->{'id2num'};
	for $i ( 0 .. $#allids ) {
		$id2numberhash->{$i} = $allids[$i];
		$masshash{$i} = $PEPINDEX->{ $allids[$i] }->{'mw'};
	}
	$ENUMERATION->{'id2num'}        = $id2numberhash;
	$ENUMERATION->{'indexprecision'} = $PARAMS->{'intprecision'};
	my %bufferhash;

	for $i ( 0 .. $#allids ) {
		for $j ( $i .. $#allids ) {
			my $newid = [
				$j, $i,
				sprintf( "%.5f",
					( $masshash{$i} + $masshash{$j} + $xlinkermw ) )
			];

			my $xlinkmass = $masshash{$i} + $masshash{$j} + $xlinkermw;
			my $intmass   = int( $xlinkmass * $ms1intprecision );

			push @{ $bufferhash{$intmass} }, $newid;
			$ncombinations++;
		}

		#enumerate monolinks

		my $monolinkmass = $masshash{$i} + $PARAMS->{'monolinkmw'};
		my $intmass      = int( $monolinkmass * $ms1intprecision );

		#get index entry
		my $masscombinations = $ENUMERATION_MONO->{$intmass};
		push @{$masscombinations}, $i;

		#writeindexentry
		$ENUMERATION_MONO->{$intmass} = $masscombinations;
		$nmonolinks++;

		#enumerate intralinks
		my $intralinkmass = $masshash{$i} + $PARAMS->{'xlinkermw'};

		#$PEPINDEX->{ $allids[$i] }->{'mw'} + $PARAMS->{'xlinkermw'};
		my $intmass = int( $intralinkmass * $ms1intprecision );

		#get index entry
		my $masscombinations = $ENUMERATION_INTRA->{$intmass};
		push @{$masscombinations}, $i;

		#writeindexentry
		$ENUMERATION_INTRA->{$intmass} = $masscombinations;
		$nintralinks++;
	}

	foreach my $intmass ( keys %bufferhash ) {
		if ( $bufferhash{$intmass} ) {
			my $masscombinations = $ENUMERATION->{$intmass};
			push @{$masscombinations}, @{ $bufferhash{$intmass} };
			$ENUMERATION->{$intmass} = $masscombinations;
		}
	}

	return ( $ncombinations, $nmonolinks, $nintralinks );
}

sub storecombinationIndex_smarthash_splitdb {
	my $PEPINDEX          = shift;
	my $PARAMS            = shift;
	my $ENUMERATION       = shift;
	my $ENUMERATION_MONO  = shift;
	my $ENUMERATION_INTRA = shift;
	my $verbose           = shift;
	my $basename          = shift;
	my $xlinkermw         = $PARAMS->{'xlinkermw'};
	my $ms1intprecision   = $PARAMS->{'intprecision'};
	my $peptideminsize    = $PARAMS->{'minpepmr'};
	open STATUS, ">index.stat" or die $!;
	my ( $i, $j, %masshash );
	my $ncombinations = 0;
	my $nmonolinks    = 0;
	my $nintralinks   = 0;
	print "sorting peptides for size ";
	my @allids =
	  sort { $PEPINDEX->{$a}->{'mw'} <=> $PEPINDEX->{$b}->{'mw'} }
	  keys %$PEPINDEX;
	print "... done\n";

	my $number2idhash = $ENUMERATION->{'id2num'};
	for $i ( 0 .. $#allids ) {
		$number2idhash->{$i} = $allids[$i];
	}
	$ENUMERATION->{'id2num'}        = $number2idhash;
	$ENUMERATION->{'indexprecision'} = $ms1intprecision;

	for $i ( 0 .. $#allids ) {
		$masshash{$i} = $PEPINDEX->{ $allids[$i] }->{'mw'};

		#print "$i ", $masshash{$i}, "\n";
	}

	my %pepindexmassbins;
	for $i ( 0 .. $#allids ) {
		push @{ $pepindexmassbins{ int( $masshash{$i} * $ms1intprecision ) } },
		  $i;

	#		push @{ $pepindexmassbins{ int( $masshash{$i} * $ms1intprecision +1) } },
	#		  $i;
	#		push @{ $pepindexmassbins{ int( $masshash{$i} * $ms1intprecision -1) } },
	#		  $i;
	}

	my $xlinkminsize =
	  int( $ms1intprecision * ( 2 * $PARAMS->{'minpepmr'} + $xlinkermw ) -
		  $ms1intprecision );
	my $xlinkmaxsize =
	  int( $ms1intprecision * ( 2 * $PARAMS->{'maxpepmr'} + $xlinkermw ) +
		  $ms1intprecision );

	my $targetxlinksize;
	my $sizerange    = $xlinkmaxsize - $xlinkminsize;
	my $percentrange = int( $sizerange / 100 );
	my $counter      = 0;
	print STATUS "scanning dtabase for $sizerange comination sizes\n";

	for $targetxlinksize ( $xlinkminsize .. $xlinkmaxsize ) {

		#		$verbose && print "############# targetxlinksize= ",
		#		  $targetxlinksize / $ms1intprecision, " ##########\n";
		$counter++;

		#		if ( $counter % $percentrange == 0 ) {
		#			print STATUS 100 * $counter / $sizerange, "\% calculated so far\n";
		#		}
		my @possible_masscombinations = ();
		for $i ( 0 .. $#allids ) {

			#			if ( ( ( $masshash{$i} + $xlinkermw ) * $ms1intprecision ) >
			#				( $targetxlinksize - $peptideminsize*$ms1intprecision ) )
			if (
				( $masshash{$i} * $ms1intprecision ) >= (
					$ms1intprecision + $targetxlinksize - $xlinkermw *
					  $ms1intprecision
				) / 2
			  )
			{
				next;
			}
			my $target_pepsize = int(
				(
					$targetxlinksize - $ms1intprecision *
					  ( $masshash{$i} + $xlinkermw )
				)
			);
			if ( $pepindexmassbins{$target_pepsize} ) {

				my @pepsinrange = @{ $pepindexmassbins{$target_pepsize} };
				foreach my $pep (@pepsinrange) {
					$ncombinations++;
					my $newid = [
						$pep, $i,
						sprintf( "%.5f",
							( $masshash{$i} + $masshash{$pep} + $xlinkermw ) )
					];

		#	print sprintf("%.4f",($masshash{$i}+$masshash{$pep}+$xlinkermw)),"\n";
					push @possible_masscombinations, $newid;

					#					if ($verbose) {
					#						print " selected pep ", $number2idhash->{$i}, " ",
					#						  $masshash{$i}, " targetpepsize= ", $target_pepsize,
					#						  " found fit:";
					#						print "number: ", $number2idhash->{$pep}, " $pep ",
					#						  $masshash{$pep};
					#					}
				}

				#				$verbose && print " \n";
			}

		}
		if ( scalar(@possible_masscombinations) > 0 ) {
			my $indexname = join "", $basename, "_", $targetxlinksize,
			  "_xlink_db";
			my %SPLITDB;
			tie( %SPLITDB, 'MLDBM', $indexname, O_CREAT | O_RDWR, 0666 )
			  or die "can't open tie to $indexname: $!";

			my $masscombinations = $SPLITDB{$targetxlinksize};
			push @{$masscombinations}, @possible_masscombinations;
			$SPLITDB{$targetxlinksize} = $masscombinations;
			$ENUMERATION->{$targetxlinksize} = $indexname;
			untie %SPLITDB;
		}

	}

	for $i ( 0 .. $#allids ) {

		#enumerate monolinks

		my $monolinkmass =
		  $PEPINDEX->{ $allids[$i] }->{'mw'} + $PARAMS->{'monolinkmw'};
		my $intmass = int( $monolinkmass * $ms1intprecision );

		#get index entry
		my $masscombinations = $ENUMERATION_MONO->{$intmass};
		push @{$masscombinations}, $i;

		#writeindexentry
		$ENUMERATION_MONO->{$intmass} = $masscombinations;
		$nmonolinks++;

		#enumerate intralinks
		my $intralinkmass =
		  $PEPINDEX->{ $allids[$i] }->{'mw'} + $PARAMS->{'xlinkermw'};
		my $intmass = int( $intralinkmass * $ms1intprecision );

		#get index entry
		my $masscombinations = $ENUMERATION_INTRA->{$intmass};
		push @{$masscombinations}, $i;

		#writeindexentry
		$ENUMERATION_INTRA->{$intmass} = $masscombinations;
		$nintralinks++;
	}
	close(STATUS);
	return ( $ncombinations, $nmonolinks, $nintralinks );
}

sub storecombinationIndex_smarthash {
	my $PEPINDEX                 = shift;
	my $PARAMS                   = shift;
	my $ENUMERATION_INTER_XLINKS = shift;
	my $ENUMERATION_INTRA_XLINKS = shift;
	my $ENUMERATION_MONO         = shift;
	my $ENUMERATION_INTRA        = shift;
	my $INFOINDEX                = shift;
	my $verbose                  = shift;
	my $statusfilehandle         = shift;

	my $xlinkermw         = $PARAMS->{'xlinkermw'};
	my $ms1intprecision   = $PARAMS->{'intprecision'};
	my $peptideminsize    = $PARAMS->{'minpepmr'};
	my $tolerance         = $PARAMS->{'ms1tolerance'};
	my $writetodiskaftern = $PARAMS->{'writetodiskaftern'};

	my %intraxlinkhash;

	if ( $PARAMS->{'tolerancemeasure'} =~ /^ppm/i ) {
		$tolerance =
		  $PARAMS->{'ms1tolerance'} * 1e-6 *
		  1000;    #ppm to amu measure for a 1000 Da peptide
	}
	my $xlinkminsize =
	  int( $ms1intprecision * ( 2 * $PARAMS->{'minpepmr'} + $xlinkermw ) -
		  $ms1intprecision );
	my $xlinkmaxsize =
	  int( $ms1intprecision * ( 2 * $PARAMS->{'maxpepmr'} + $xlinkermw ) +
		  $ms1intprecision );
	my $nkeys = ( $xlinkmaxsize - $xlinkminsize ) * $ms1intprecision;
	open STATUS, ">index.stat" or die $!;
	my ( $i, $j, %masshash );
	my $ncombinations          = 0;
	my $intraxlinkcombinations = 0;
	my $nmonolinks             = 0;
	my $nintralinks            = 0;
	my @allids;

	#critical
	#if ( $PARAMS->{'search_intercrosslinks'} ) {
	print $statusfilehandle "sorting peptides for size ";

	my @temp = map {[$PEPINDEX->{$_}->{'mw'},$_]} keys %$PEPINDEX;
	@temp=sort {$a->[0] <=> $b->[0]} @temp;
	@allids=map {$_->[1]} @temp;
	
	print $statusfilehandle "... done\n";

	my $number2idhash = $INFOINDEX->{'id2num'};
	for $i ( 0 .. $#allids ) {
		$number2idhash->{$i} = $allids[$i];
	}
	$INFOINDEX->{'id2num'}= $number2idhash; 


	for $i ( 0 .. $#allids ) {
		$masshash{$i} = $PEPINDEX->{ $allids[$i] }->{'mw'};
	}

	if ( $PARAMS->{'search_intracrosslinks'} ) {
		print $statusfilehandle
		  "creating intra-protein cross-link index in smarthash mode\n";

		print "creating intra-protein cross-link index in smarthash mode\n";
		my ( $protid, $k, @idsofprot );
		print $statusfilehandle "sorting pepids to proteins ...\n";
		for $i ( 0 .. $#allids ) {
			my @proteinsIDs = @{ $PEPINDEX->{ $allids[$i] }->{'proteinID'} };
			foreach my $protein (@proteinsIDs) {
				push @{ $intraxlinkhash{$protein} }, $i;
			}

		}
		my @protids         = keys %intraxlinkhash;
		my $numberofprotids = scalar(@protids);
		print "indexing intra-protein cross-links\n";
		my %bufferhash = ();
		keys(%bufferhash) = $nkeys;
		my $nprots               = 0;
		my $combinationsinmemory = 0;

		foreach $protid (@protids) {
			@idsofprot = @{ $intraxlinkhash{$protid} };
			print $statusfilehandle "indexing $protid\n";
			for $i ( 0 .. $#idsofprot ) {
				for $j ( $i .. $#idsofprot ) { # Shouldn't j = i be an inter-protein cross-link?? 
					my $newid = [
						$idsofprot[$j],
						$idsofprot[$i],
						sprintf(
							"%.5f",
							(
								$masshash{ $idsofprot[$i] } +
								  $masshash{ $idsofprot[$j] } + $xlinkermw
							)
						)
					];

					my $xlinkmass = $masshash{ $idsofprot[$i] } +
					  $masshash{ $idsofprot[$j] } + $xlinkermw;
					my $intmass = int( $xlinkmass * $ms1intprecision );
					if (   $intmass >= $xlinkminsize
						&& $intmass <= $xlinkmaxsize )
					{
						push @{ $bufferhash{$intmass} }, $newid;
						$intraxlinkcombinations++;
						$combinationsinmemory++;

					}

				}
			}
			$nprots++;
			if (   ( $nprots >= $numberofprotids )
				|| ( $combinationsinmemory > $writetodiskaftern ) )
			{
				print $statusfilehandle
"$nprots: temporary index contains $combinationsinmemory combinations \n->writing to disk after $protid... ";

				$combinationsinmemory = 0;

				foreach my $intmass ( keys %bufferhash ) {
					if ( $bufferhash{$intmass} ) {
						my $masscombinations =
						  $ENUMERATION_INTRA_XLINKS->{$intmass};
						push @{$masscombinations}, @{ $bufferhash{$intmass} };
						$ENUMERATION_INTRA_XLINKS->{$intmass} =
						  $masscombinations;
					}
				}
				print $statusfilehandle "... done\n";
				%bufferhash = ();
				keys(%bufferhash) = $nkeys;
			}

		}
		my %bufferhash;

	}

	if ( $PARAMS->{'search_intercrosslinks'} ) {
		print "indexing inter-protein cross-links....\n";
		print $statusfilehandle "indexing inter-protein cross-links....\n";
		print $statusfilehandle "--> indexing massbins...";
		my %pepindexmassbins;
		for $i ( 0 .. $#allids ) {
			foreach my $tolerant_intmass (
				int_with_tolerance(
					$masshash{$i}, $tolerance, $ms1intprecision
				)
			  )
			{
				push @{ $pepindexmassbins{$tolerant_intmass} }, $i;
			}
		}
		print $statusfilehandle "... done\n";


		my $targetxlinksize;
		$verbose = 0;
		my %seen = ();
		my $k    = 0;

		for $targetxlinksize ( $xlinkminsize .. $xlinkmaxsize ) {

			$seen{ $k - 2 } = ();
			$verbose && print "############# targetxlinksize= ", $targetxlinksize, " ####################\n";

			my @possible_masscombinations = ();
			 for $i ( 0 .. $#allids ) {
#
#				if (
#					( $masshash{$i} * $ms1intprecision ) >= (
#						$ms1intprecision + $targetxlinksize
#						+ $xlinkermw *  $ms1intprecision
#					) / 2
#				  )
#				{
#					if($verbose){
#					print "current mass = ",$masshash{$i} * $ms1intprecision," -->next target\n";
#					}
#				next target;
#				}
				my $target_pepsize = int(
					(
						$targetxlinksize - $ms1intprecision *
						  ( $masshash{$i} + $xlinkermw )
					)
				);
				if ( $pepindexmassbins{$target_pepsize} ) {
					my @pepsinrange = @{ $pepindexmassbins{$target_pepsize} };
					foreach my $pep (@pepsinrange) {

						my $idstring = join "::",
						  (	sort ( $number2idhash->{$i},$number2idhash->{$pep} ) );
								
								
						if (   !( $seen{ $k - 1 }->{$idstring}++ )
							&& !( $seen{$k}->{$idstring}++ ) )
						{
							$ncombinations++;
							my $newid = [
								$pep, $i,
								sprintf(
									"%.5f",
									(
										$masshash{$i} + $masshash{$pep} +
										  $xlinkermw
									)
								)
							];
							push @possible_masscombinations, $newid;

							if ($verbose) {
								print " selected pep ", $number2idhash->{$i},
								  " ", $masshash{$i}, " targetpepsize= ",
								  $target_pepsize, " found fit:";
								print "number: ", $number2idhash->{$pep},
								  " $pep ", $masshash{$pep}, " final size ",
								  sprintf(
									"%.5f",
									(
										$masshash{$i} + $masshash{$pep} +
										  $xlinkermw
									)
								  ),
								  " delta: ",
								  abs(
									$targetxlinksize / $ms1intprecision - (
										$masshash{$i} + $masshash{$pep} +
										  $xlinkermw
									)
								  ),
								  , "\n";
							}
						}
					}
				}

			}
			if ( scalar(@possible_masscombinations) > 0 ) {
				my $masscombinations = $ENUMERATION_INTER_XLINKS->{$targetxlinksize};
				push @{$masscombinations}, @possible_masscombinations;
				$ENUMERATION_INTER_XLINKS->{$targetxlinksize} = $masscombinations;
			}
			$k++;
		}
	}
	if ( $PARAMS->{'search_monolinks'} ) {
		print "indexing monolinks ...\n";
		print $statusfilehandle "indexing monolinks ...\n";
		my %bufferhash = ();
		my @monolinkmasses = split /,/, $PARAMS->{'monolinkmw'};
		for $i ( 0 .. $#allids ) {

			#enumerate monolinks

			foreach my $mass (@monolinkmasses) {
				my $monolinkmass = $PEPINDEX->{ $allids[$i] }->{'mw'} + $mass;
				my $intmass      = int( $monolinkmass * $ms1intprecision );

				#get index entry
				push @{ $bufferhash{$intmass} }, $i;
				$nmonolinks++;
			}
		}
		foreach my $intmass ( keys %bufferhash ) {
			if ( $bufferhash{$intmass} ) {
				my $masscombinations = $ENUMERATION_MONO->{$intmass};
				push @{$masscombinations}, @{ $bufferhash{$intmass} };
				$ENUMERATION_MONO->{$intmass} = $masscombinations;
			}
		}
	}

	if ( $PARAMS->{'search_intralinks'} ) { #also called loop links
		print "indexing intralinks ...\n";
		print $statusfilehandle "indexing intralinks ...\n";
		my %bufferhash = ();

		my @xlinkmasses = split /,/, $PARAMS->{'monolinkmw'};
		for $i ( 0 .. $#allids ) {

			#enumerate intralinks
			foreach my $mass (@xlinkmasses) {
				my $intralinkmass =
				$PEPINDEX->{ $allids[$i] }->{'mw'} + $mass; #$PARAMS->{'xlinkermw'};
				my $intmass = int( $intralinkmass * $ms1intprecision );
				push @{ $bufferhash{$intmass} }, $i;
				$nintralinks++;
			}
		}
		foreach my $intmass ( keys %bufferhash ) {
			if ( $bufferhash{$intmass} ) {
				my $masscombinations = $ENUMERATION_INTRA->{$intmass};
				push @{$masscombinations}, @{ $bufferhash{$intmass} };
				$ENUMERATION_INTRA->{$intmass} = $masscombinations;
			}
		}

	}
	close(STATUS);
	return ( $ncombinations, $intraxlinkcombinations, $nmonolinks,$nintralinks );
}

sub round {
	my ($number) = shift;
	return int( $number + .5 );
}

sub int_with_tolerance {
	my $number       = shift;
	my $tolerance    = shift;
	my $intprecision = shift;
	my %intnumbers;
	$intnumbers{ int( $number * $intprecision ) } = defined;
	$intnumbers{ int( ( $number + $tolerance ) * $intprecision ) } = defined;
	$intnumbers{ int( ( $number - $tolerance ) * $intprecision ) } = defined;

	my @intnumbers = keys %intnumbers;

	#	print $number+$tolerance," ",$number-$tolerance," $number @intnumbers \n";
	return @intnumbers;
}

sub storeIndex {
	my $self             = shift;
	my $pepobj           = shift;
	my $IONINDEX         = shift;
	my $PEPINDEX         = shift; ## hashref
	my $PARAMS           = shift;
	my $verbose          = shift;
	my $protid           = shift;
	my $statusfilehandle = shift;
	my ( @entry1, $ion, $intprecision );
	my $minionsize    = $PARAMS->{'minionsize'};
	my $maxionsize    = $PARAMS->{'maxionsize'};
	my $enumerate     = $PARAMS->{'enumerate'};
	my $matchall      = $PARAMS->{'matchall'};
	my $massmatchonly = $PARAMS->{'massmatchonly'};
	$intprecision = $PARAMS->{'ionindexintprecision'};

	my $redundant = undef;
	if ( defined( $PEPINDEX->{ $pepobj->id } ) ) {
		$redundant = defined;
		#print "TRUE";
	}

	#even if peptide is redundant peptide is stored again with addition of protein id. but ions are not again indexed
	my $peptideinfo = $PEPINDEX->{ $pepobj->id };
	## $PEPINDEX is a hash in hash first key is the peptidesequence, second keys are seq, desc, mw, and protid (arrayreference)
	$peptideinfo->{'seq'}  = $pepobj->seq;
	$peptideinfo->{'desc'} = $pepobj->desc;
	$peptideinfo->{'mw'}   = $pepobj->molweight;
	push @{ $peptideinfo->{'proteinID'} }, $protid;   #add protein id to peptide

	$PEPINDEX->{ $pepobj->id } = $peptideinfo;
#	print Dumper ($peptideinfo);
#	print "\nindexing ", $pepobj->id, "\n";
#	exit;
	$verbose && $pepobj->printtable;
	unless ( $enumerate	|| defined($redundant)	|| $massmatchonly	|| $matchall )
	{
		my $ions = $pepobj->iontag_getionsforindex;
		
		foreach $ion (@$ions) {
			my $intion = int( $ion * $intprecision );    #get 1/inpresciosn resolution
			#$IONINDEX->{$intion}=[];
			if ( $ion >= $minionsize && $ion <= $maxionsize ) {
				## deref $IONINDEX->{$intion} (which is an arrayref) and push in
				push @{ $IONINDEX->{$intion} }, $pepobj->id;
			}
		}
	}
	if ( defined($redundant) ) {
		return 0;
	}
	else {
		return 1;
	}
}

sub contains {
	my $ionindexentry = shift;
	my $id            = shift;
	my %seen          = ();
	foreach my $id (@$ionindexentry) {
		if ( $seen{$id}++ ) {

			#print "seen $id in ",join " ",@$ionindexentry,"\n";
			return 1;
		}
	}
	return 0;
}

sub printenumeration {
	my $ENUMERATIONINDEX = shift;
	foreach ( keys %$ENUMERATIONINDEX ) {
		my @combinations = @{ $ENUMERATIONINDEX->{$_} };
		print "$_ @combinations\n";
	}
}

1;
