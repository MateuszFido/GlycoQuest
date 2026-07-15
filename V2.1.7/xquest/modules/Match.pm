package Match;
use strict;
#---------------------------------------------------------------------------
# Module: Match.pm
# Author(s): Thomas Walzthoeni, xQuest specific modules are based on 
# original work by Oliver Rinner.
# Description: Search Module.
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
use Spectrum;
use LinkObj;
use Data::Dumper;

sub new
{
	my $class = shift();
	my $self  = {};
	bless $self, $class;
	my $spectrumobj = shift;
	my $IONINDEX    = shift;
	my $PEPINDEX    = shift;
	my $PARAMS      = shift;
	my $MSTAB       = shift;
	my $verbose     = shift;
	my $filehandle  = shift;

	#	$self->{'spectrum'} = basename($spectrum);
	$self->{'PARAMS'}   = $PARAMS;
	$self->{'spectrum'} = $spectrumobj;
	$self->{'MSTAB'}    = $MSTAB;
	
	if ( $PARAMS->{'enumerate'} )
	{    #enumeration mode1 slow index build fast search
		$self->enumerate;
	} elsif ( $PARAMS->{'matchall'} )
	{    #enumeration mode2 fast index build slow search
		( $self->{'hitlist'}, $self->{'sortedbymass'} ) = $self->matchall( $spectrumobj, $IONINDEX, $PEPINDEX, $PARAMS, $MSTAB, $verbose, $filehandle );
	} elsif ( $PARAMS->{'massmatchonly'} )
	{
		( $self->{'hitlist'}, $self->{'sortedbymass'} ) = $self->massmatch( $spectrumobj, $IONINDEX, $PEPINDEX, $PARAMS, $MSTAB, $verbose, $filehandle );
	} else
	{    # ionindex mode
		( $self->{'hitlist'}, $self->{'sortedbymass'} ) = $self->match_iontag( $spectrumobj, $IONINDEX, $PEPINDEX, $PARAMS, $MSTAB, $verbose, $filehandle );
	}
	return $self;
}

sub set_infoindex{
my $self=shift;
my $name=shift;
my $index=shift;
$self->{$name} = $index;
}

sub get_infoindex{
my $self=shift;
my $name=shift;
return $self->{$name};
}


sub new_matchobj_dummy
{
	my $class = shift();
	my $self  = {};
	bless $self, $class;
	my $spectrumobj = shift;
	$self->{'spectrum'} = $spectrumobj;
	
	return $self;
}

sub gethithash
{
	my $self = shift;
	return $self->{'sortedbymass'};
}

sub getMSTAB
{
	my $self = shift;
	return $self->{'MSTAB'};
}

sub getspectrumname
{
	my $self = shift;
	return $self->getSpecObj->getSpecname;
}

sub getms1mass
{
	my $self = shift;
	return $self->getSpecObj->getms1mass;
}

sub getms1Mz
{
	my $self = shift;
	return $self->getSpecObj->getprecursorMz;
}

sub getms1Charge
{
	my $self = shift;
	return $self->getSpecObj->getprecursorCharge;
}

sub printhits
{
	my $self             = shift;
	my $statusfilehandle = shift;
	my $hitlist          = $self->hitlist;
	print $statusfilehandle "<p>candidate peptides with > ", $self->getParams->{'minhits'}, " hits<br>\n";
	my $counter=1;
	foreach ( sort { $hitlist->{$b}->{'prescore'} <=> $hitlist->{$a}->{'prescore'} } ( keys %{$hitlist} ) )
	{
		print $statusfilehandle "hit $counter:", $_, "\t";
		print $statusfilehandle $hitlist->{$_}->{'header'}, "\t";
		print $statusfilehandle "spectrum peaks matched: ", $hitlist->{$_}->{'matchratio'},  "\t";
		print $statusfilehandle "nhits common ions: ",      $hitlist->{$_}->{'ncommonhits'}, "\t";
		print $statusfilehandle "nhits xlink ions: ",       $hitlist->{$_}->{'nxlinkhits'},  "\t";
		print $statusfilehandle "common ions tested: ",     $hitlist->{$_}->{'nions'},       "\t";
		print $statusfilehandle "xlink ions tested: ",      $hitlist->{$_}->{'nxlinkions'},  "\t";
		print $statusfilehandle "intsumtotal: ", sprintf( "%.2f", $hitlist->{$_}->{'intsumtotal'} ), "\t";
		print $statusfilehandle "prescore: ", sprintf( "%.2f", $hitlist->{$_}->{'prescore'} ), "\n";

	$counter++;
	}
}

sub hitlist
{
	my $self = shift;
	return $self->{'hitlist'};
}

sub match_iontag
{
	my $self       = shift;
	my $spectrum   = shift;
	my $ionindex   = shift;
	my $pepindex   = shift;
	my $PARAMS     = shift;
	my $MSTAB      = shift;
	my $verbose    = shift;
	my $filehandle = shift;
	my $ms1mass    = $self->getms1mass;
	my $maxpepsize = $ms1mass;
	
	# Define mass of H2O
	my $H2O=2 * $MSTAB->{'Hatom'}->{'native'} + $MSTAB->{'Oatom'}->{'native'};

	$verbose = 0;
	my ( $intprecision, $ion, $peak, $testion );
	my $minhits = $PARAMS->{'minhits'};
	unless ( $PARAMS->{'miniontaghits'} )
	{
		$PARAMS->{'miniontaghits'} = 1;
	}
	my $miniontaghits     = $PARAMS->{'miniontaghits'};
	my $MS2tolerance      = $PARAMS->{'ms2tolerance'};
	my $MS1tolerance      = $PARAMS->{'ms1tolerance'};
	my $MS2xlinktolerance = $PARAMS->{'xlink_ms2tolerance'};
	my $picktolerance     = $PARAMS->{'picktolerance'};
	unless ( $PARAMS->{'search_maxcandidate_peps'} )
	{
		$PARAMS->{'search_maxcandidate_peps'} = 10000;
	}
	my $searchmaxncandidates = $PARAMS->{'search_maxcandidate_peps'};
	if ( $PARAMS->{'tolerancemeasure'} =~ /^ppm/i )
	{
		$MS1tolerance = $MS1tolerance * 1e-6 * $ms1mass;    #ppm to amu measure
	}

## IF EDAC is used as cross-linker add the positive mass tolerance to the precursor mass plus H2O (for looplinks), otherwise peptides with a positive shift 
## are not considered, which misses monolink candidates and looplink candidates (+ H2O) since no additional mass is added by the cross-linker
if ($PARAMS->{'xkinkerID'} eq "EDAC" ){
$maxpepsize=$maxpepsize+$MS1tolerance+$H2O;
}

	
	if ( $PARAMS->{'picktolerance_measure'} =~ /^ppm/i )
	{
		$verbose && print "picktolerance adjusted dynamically for each ion: $picktolerance ppm\n";
		print $filehandle "picktolerance adjusted dynamically for each ion: $picktolerance ppm\n";
	} else
	{
		print $filehandle "picktolerance: $picktolerance m/z\n";
		$verbose && print "picktolerance: $picktolerance m/z\n";
	}
	if ( $PARAMS->{'tolerancemeasure_ms2'} =~ /^ppm/i )
	{
		$verbose && print "matchtolerance adjusted dynamically for each ion: $MS2tolerance ppm\n";
		print $filehandle "matchtolerance adjusted dynamically for each ion: $MS2tolerance ppm\n";
	} else
	{
		$verbose && print "matchtolerance: $MS2tolerance m/z\n";
		print $filehandle "matchtolerance: $MS2tolerance m/z\n";
	}
	if ( $PARAMS->{'ionindexintprecision'} )
	{
		$intprecision = $PARAMS->{'ionindexintprecision'};
	} else
	{
		$intprecision = 10;
	}
	## Todo: CHECK IF ITS NECCASARY TO do +1!
	my $inttolerance = int( $intprecision * $picktolerance ) + 1;

	#my $inttolerance = int( $intprecision * $picktolerance ) ;
	my ( %matchlist, %matchlistsorted );
	my %seenpeps    = {};
	my @commonpeaks = sort { $a <=> $b } @{ $spectrum->getcommonpeaks };
	my @testionlist = @{ $spectrum->gettestions };
	my $ionintensityhash = $spectrum->get_ionintensityhash;
	
	$verbose && print "peaks in spectrum: ", scalar(@commonpeaks), "\n";
	print $filehandle "peaks in spectrum: ", scalar(@commonpeaks), "\n";
	my @candidateids = ();
	my %seen         = ();
	
	#### SELECTION OF THE CANIDATE PEPTIDES ###
	my $ids_ions = {};
	$verbose && print "Match.pm: Testing " . scalar(@testionlist) . " common-ions\n";
	foreach $testion (@testionlist)
	{
		my $intionmz = int( $testion * $intprecision );    #get intprecision resolution
		if ( $PARAMS->{'picktolerance_measure'} =~ /^ppm/i )
		{
			$picktolerance = $PARAMS->{'picktolerance'} * 1e-6 * $testion;
			$inttolerance  = int( $intprecision * $picktolerance ) + 1;
		}

		#print "picktolerance for ion with $testion m\/z: $picktolerance m/z\n";
		#print "inttolerance: $inttolerance\n";
		# exit;
		# get all peptides that contain the specific candidate ion
		# $idinrange are peptide seqs
		for my $rangex ( $intionmz - $inttolerance ... $intionmz + $inttolerance )
		{
			foreach my $idinrange ( @{ $ionindex->{$rangex} } )
			{
				unless ( $seen{$idinrange}++ )
				{
					push @candidateids, $idinrange;
					## id in range is a peptideseq
					#print $idinrange."\n";
					push @{ $ids_ions->{$idinrange} }, $testion;
					
					#$ids_ions->{$idinrange} = $testion;
				} else
				{
					my $ionseen = $ids_ions->{$idinrange};
					## check if the seen ion is larger, otherwise put the smallest ion into the hash
					if ( $ionseen > $testion )
					{
						push @{ $ids_ions->{$idinrange} }, $testion;

						#$ids_ions->{$idinrange} = $testion;
					}
				}
			}
		}
	}
	
	## Get all candidates
#	my @indexkeys= keys %$ionindex;
#	foreach my $key (@indexkeys){
#	foreach my $idinrange ( @{ $ionindex->{$key} } ){
#		unless ( $seen{$idinrange}++ )
#				{
#				#push @candidateids, $idinrange;	
#				}
#	}
#	}
#my $numcandidates=	@candidateids;
#print "Number of candidates: $numcandidates\n";
	
	
	#exit;
	
	
	
	
	
	
	
	$verbose && print "First pass: Number of unique candidate peptides associated to ion tags: " . scalar(@candidateids) . "\n";
	print $filehandle "number of candidate peptides associated to ion tags: ", scalar(@candidateids), "\n";

	#print Dumper ($ids_ions);
	### Delete all candidates with smaller number of seenions than defined in $miniontaghits
	foreach my $id ( keys %seen )
	{

		#$verbose && print "$id was seen ", $seen{$id}, " times\n";
		## ONLY FOR DC validation
		my $maxiontaghits = $PARAMS->{'maxiontaghits'};
		if ($maxiontaghits)
		{

			#print "MaxIThits :$maxiontaghits, Number of seenions: $seen{$id}\n";
			unless ( ( $seen{$id} >= $miniontaghits ) && ( $seen{$id} <= $maxiontaghits ) )
			{

				#print "Will be deleted!\n";
				delete( $seen{$id} );
				delete( $ids_ions->{$id} );
			}
		} else
		{
			unless ( ( $seen{$id} >= $miniontaghits ) )
			{

				#print "Will be deleted!\n";
				delete( $seen{$id} );
				delete( $ids_ions->{$id} );
			}
		}
	}
	## Push all candidates in an array
	my @considered_candidateids = ();
	foreach my $candidateid (@candidateids)
	{
		if ( $seen{$candidateid} )
		{
			push @considered_candidateids, $candidateid;
		}
	}
	@candidateids = ();
	$verbose && print "Second pass: Number of candidate peptides associated to ion tags >=$miniontaghits miniontaghits: ", scalar(@considered_candidateids), "\n";
	print $filehandle "number of candidate peptides associated to ion tags >=$miniontaghits: ", scalar(@considered_candidateids), "\n";
	my @candidatepepObjs = ();
	my $ncandidates      = 0;
	my $discarded        = 0;
	## Generate Pepobjects for all candidate peptides that are smaller than the precursor
	foreach my $id (@considered_candidateids)
	{
		my $indexentry = $pepindex->{$id};
		my $desc       = $indexentry->{'desc'};
		my $seq        = $indexentry->{'seq'};
		my $molweight  = $indexentry->{'mw'};
		my $proteinID  = $indexentry->{'proteinID'};
		my $seenion    = $ids_ions->{$id};             ### is now an arrayreference

		#print "Seenion for peptide $seq was $seenion\n";
		if ( $molweight > $maxpepsize )
		{

			#$verbose && print "discarding $id ... too big\n";
			$discarded++;
			next;
		}

		#$verbose && print "creating new peptide object for $id\n";
		my $pepobj = PepObj->new( $seq, $id, $desc, $MSTAB, $PARAMS, $verbose, $proteinID, $seenion );
		$seenpeps{$id} = $pepobj;
		push @candidatepepObjs, $pepobj;
	}
	$verbose && print "Third Pass: Discarded $discarded pepobjects >> peptide too large (> Mr precursor)\n";
	$verbose && print "Generated " . scalar(@candidatepepObjs) . " canditate peptide objects\n";
	my $i = 0;
	
	## Matching of the canidate petides (single pep.)
	my $timestart = time;
	foreach my $candidate (@candidatepepObjs)
	{
		$i++;
		my $ntestions             = 0;
		my $nxlinktestions        = 0;
		my $npeaks                = 0;
		my $candidate_common_ions = getcandidate_common_ions( $candidate, $PARAMS );
		my $nhits                 = 0;
		my $ncommonhits           = 0;
		my $nxlinkhits            = 0;
		
		my $intsumcommon=0;
		my $intsumxlink=0;
		
		$ntestions += scalar(@$candidate_common_ions);
		$npeaks    += scalar(@commonpeaks);

my $sequence =  $candidate->seq;

#if ($sequence eq "LTKIENKVDK"){
#print "Sequence is KIHLKELITK\n";
#print "CAndidate common ions\n";
#print Dumper ($candidate_common_ions);
#my $deltaMr      = $self->getms1mass - $candidate->molweight;
#print "MS1 mass: $self->getms1mass, candidate mass: $candidate->molweight, deltamr: $deltaMr\n";
##print "ProtId: ";
##print Dumper($candidate->protid);
#print "\n";
##exit;	
##print Dumper ($ionintensityhash);
#print "Commonpeaks in spectrum: \n";
#print Dumper (\@commonpeaks);
##exit;
#}
		#print "Counting Hits\n";
		#$ncommonhits += $self->counthits( \@commonpeaks, $candidate_common_ions, $MS2tolerance, $ionintensityhash );
		($ncommonhits, $intsumcommon) =  $self->counthits( \@commonpeaks, $candidate_common_ions, $MS2tolerance, $ionintensityhash );
		
		if ( $PARAMS->{'iontag_match_xlinkions'} )
		{
			my $mincharge     = $self->minioncharge_xlinks;
			my $maxcharge     = $self->maxioncharge_xlinks;
			my $xlinkionshash = $self->getcandidate_xlink_ions( $candidate, $PARAMS );

			#print Dumper ($xlinkionshash);
			#	exit;
			for my $charge ( $mincharge .. $maxcharge )
			{
				my @xlinkpeaks = ( @{ $spectrum->getxlinkpeaks($charge) }, @{ $spectrum->getxlinkpeaks("0") } );

				#print Dumper (\@xlinkpeaks);
				my $xlinkions = $xlinkionshash->{$charge};
				$nxlinktestions += scalar(@$xlinkions);
				$npeaks         += scalar(@xlinkpeaks);
				#$nxlinkhits     += $self->counthits( \@xlinkpeaks, $xlinkions, $MS2xlinktolerance, $ionintensityhash );
			my ($nxltmp, $intsumxltmp)=$self->counthits( \@xlinkpeaks, $xlinkions, $MS2xlinktolerance, $ionintensityhash );
			$nxlinkhits+=$nxltmp;
			$intsumxlink+=$intsumxltmp;
			}

			#exit;
		}
		$nhits = $ncommonhits + $nxlinkhits;
		## print the canditate
		#$verbose && print "creating new peptide object for ", $candidate->seq," COMMONHITS: $ncommonhits XLhits: $nxlinkhits\n";
		if ( $nhits >= $minhits )
		{
			my $newid = $candidate->id;
			$matchlist{$newid}->{'pepObj'} = $candidate;
			### HERE THE PRESCORE FOR THE SINGLEPEPS IS CALCULATED
			$matchlist{$newid}->{'prescore'}    = $nhits / $candidate->getlength;
			$matchlist{$newid}->{'nxlinkhits'}  = $nxlinkhits;
			$matchlist{$newid}->{'ncommonhits'} = $ncommonhits;
			$matchlist{$newid}->{'nions'}       = $ntestions;                       #scalar(@candidate_common_ions);
			$matchlist{$newid}->{'nxlinkions'}  = $nxlinktestions;                  #scalar(@candidate_common_ions);
			$matchlist{$newid}->{'matchratio'}  = join "", $nhits, "/", $npeaks;    #scalar(@commonpeaks);
			$matchlist{$newid}->{'intsumtotal'}  = $intsumcommon + $intsumxlink; 
			$ncandidates++;
		}
	}
	my $matchingduration = time - $timestart;
	$verbose && print "Matched $i single peptide candidates. Matching took $matchingduration s\n";
	$verbose && print "Fourth pass: Filtered $ncandidates candidates with at least $minhits Iontag hits (COMMON and XL)\n";
	my (@candidatekeys);
	### Sort out the nbest candidates
	### Put them into the matchlistsorted which is indexed by ms1mass*100
	if ( $ncandidates > $searchmaxncandidates )
	{
		my @temp = map { [ $matchlist{$_}->{'prescore'}, $_ ] } keys %matchlist;
		@temp = sort { $b->[0] <=> $a->[0] } @temp;
		my @candidate_pepobs = map { $_->[1] } @temp;
		my %bestcandidatehash;
		@candidatekeys = @candidate_pepobs[ 0 .. $searchmaxncandidates - 1 ];
		my $counter=1;
		
		foreach my $pep (@candidatekeys)
		{
#			### Generate semitryptic peptides of these candidates // not used yet
#			my @semitryptic=generate_semitryptic_candidates($matchlist{$pep}->{'pepObj'}, $PARAMS, $MSTAB);
#
#			foreach my $semitryppepobj (@semitryptic){
#			print "Generate semitryptic candidate: $semitryppepobj->seq\n";
#			push @{ $matchlistsorted{ int( 100 * $semitryppepobj->molweight ) } }, $semitryppepobj;
#			$bestcandidatehash{$pep} = $matchlist{$pep};			
#			}
			### set the rank
			$matchlist{$pep}->{'pepObj'}->{'iontagrank'}=$counter;
			$counter++;
			push @{ $matchlistsorted{ int( 100 * $matchlist{$pep}->{'pepObj'}->molweight ) } }, $matchlist{$pep}->{'pepObj'};
			$bestcandidatehash{$pep} = $matchlist{$pep};
		}
		$verbose && print "Fifth pass: discarding ", scalar(@candidate_pepobs) - $searchmaxncandidates, " candidates with number of hits >= $minhits. Keeping best ", $searchmaxncandidates, " matches\n";
		print $filehandle "discarding ", scalar(@candidate_pepobs) - $searchmaxncandidates, " candidates with number of hits >= $minhits. Keeping best ", $searchmaxncandidates, " matches\n";
		$self->{'nhits'} = scalar(@candidatekeys);
		return \%bestcandidatehash, \%matchlistsorted;
	} else
	{
		@candidatekeys = keys(%matchlist);
		foreach my $pep (@candidatekeys)
		{
			push @{ $matchlistsorted{ int( 100 * $matchlist{$pep}->{'pepObj'}->molweight ) } }, $matchlist{$pep}->{'pepObj'};
		}
		print $filehandle "Number of candidates with  >", $minhits, " matches: ", scalar(@candidatekeys), "\n";
		$verbose && print "Number of candidates with  >", $minhits, " matches: ", scalar(@candidatekeys), "\n";
		$self->{'nhits'} = scalar(@candidatekeys);
		return \%matchlist, \%matchlistsorted;
	}
}

sub generate_semitryptic_candidates{
my $pepobject=shift;
my $PARAMS=shift;
my $MSTAB=shift;

my $verbose;
my @pepobjects;

my $sequence = $pepobject->seq;
#print "Sequence is: $sequence\n";

### Start depleting AAs from the end. Check if the peptide is still xlinkable.
my @seq = split //, $sequence;
#print Dumper (@seq);

my $length = @seq;

for my $i ( 0 .. $length-1 ) {
my $offset = $length -$i;
my $subsequence = substr $sequence, 0, $offset;
my $check = test_xlinkrequirements($subsequence,$PARAMS);
#print "Subseq: $subsequence :check $check\n";

if ($check){
## Generate pepObject
#print "Subseq: $subsequence :check $check\n";
my $desc       = $pepobject->{'desc'};
my $seq        = $subsequence;
my $proteinID  = $pepobject->{'proteinID'};
my @seenions;

my $pepobj=PepObj->new( $seq, 0, $desc, $MSTAB, $PARAMS, $verbose, $proteinID, \@seenions );
push @pepobjects, $pepobj;	
}else{
last;	
}
}	

return @pepobjects;

}




sub test_xlinkrequirements {
	my $sequence   = shift;
	my $PARAMS     = shift;
	my $xids       = shift;
	my @xlinksites = ();

	while ( $sequence =~ /$PARAMS->{'AArequired'}/gi ) {
		push @xlinksites, pos($sequence);
	}

	my $nrxlinksites = $#xlinksites + 1;

	if ( $nrxlinksites == 0 ) {

		#	$superverbose && print "$xids $sequence -> no xlinksite discarded\n";
		return 0;
	}
	elsif ( $nrxlinksites > 1 ) {

 #	$superverbose &&  print "$xids $sequence nr of xlinks = ",$nrxlinksites,"\n";
		return 1;
	}

	elsif (( $nrxlinksites == 1 )
		&& ( $xlinksites[0] < length($sequence) ) )
	{    #only one K but in the middle of peptide
		 #	print "$xids $sequence nr sites only ",$nrxlinksites," but in the middle at:",$xlinksites[0],"\n";
		return 1;
	}
	elsif ($PARAMS->{'nocutatxlink'}
		&& ( $nrxlinksites == 1 )
		&& ( $xlinksites[0] == length($sequence) ) )
	{    #only one K at the end of sequence
		 #	print "$xids $sequence nr sites only ",$nrxlinksites," at the end at ",$xlinksites[0]," discarded\n";
		return 0;
	}
	elsif (!( $PARAMS->{'nocutatxlink'} )
		&& ( $nrxlinksites == 1 )
		&& ( $xlinksites[0] == length($sequence) ) )
	{    #only one K but in the middle of peptide
		 #	print "$xids $sequence nr sites only ",$nrxlinksites," at the end at ",$xlinksites[0]," accepted\n";
		return 1;
	}
	else {
		warn "problem with $sequence $!";
		return 0;
	}

}



sub add_dc_candidates_to_matchlists
{
	my $self                 = shift;
	my $matchlisttoadd       = shift;
	my $matchlistsortedtoadd = shift;

	#print (keys %$matchlisttoadd);
### merge the matchlists (indexed by pepseq)
	my $matchlist      = $self->hitlist;
	my $nredundantpeps = 0;
	my $nadded         = 0;
	my $redundantpeps  = {};
	my @nummatchlist = keys %$matchlist;
	my $count = scalar (@nummatchlist);
	
### Merge the matchlists
	foreach my $pepkey ( keys %$matchlisttoadd )
	{
		if ( $matchlist->{$pepkey} )
		{
			$nredundantpeps++;
			$redundantpeps->{$pepkey} = 1;
		} else
		{
### Add it to the matchlist
			$matchlist->{$pepkey} = $matchlisttoadd->{$pepkey};
			$nadded++;
		}
	}
	#print "Merged matchlists: #of target peptides in matchlist $count # of decoy peptides added: $nadded \n";
### Merge the sorted matchlist (key is bin, value arrayref with pepobs)
	my $matchlistsorted = $self->gethithash;
	$nredundantpeps = 0;
	$nadded         = 0;
	foreach my $bin ( keys %$matchlistsortedtoadd )
	{
### is an arrayref
		my $pepobjref = $matchlistsortedtoadd->{$bin};
## get the pepobjects out
		foreach my $pepobj (@$pepobjref)
		{

			#print "Seq is ". $pepobj->seq."\n";
			my $pepseq = $pepobj->seq;
### Check if this seq is redundant in the target peplist
			if ( $redundantpeps->{$pepseq} )
			{
				$nredundantpeps++;
			} else
			{
### Add it to the matchlistsorted
				push @{ $matchlistsorted->{$bin} }, $pepobj;
				$nadded++;
			}
		}
	}
	#print "Merged sorted matchlists: # of redundant peptides: $nredundantpeps, # of peptides added: $nadded\n";
	return;
}

sub enumerate
{
	my $self = shift;
}

sub matchall
{
	my $self       = shift;
	my $spectrum   = shift;
	my $ionindex   = shift;
	my $pepindex   = shift;
	my $PARAMS     = shift;
	my $MSTAB      = shift;
	my $verbose    = shift;
	my $filehandle = shift;
	my $ms1mass    = $self->getms1mass;
	my ( @peaks, $ion, $peak, $testion );
	my $minhits       = $PARAMS->{'minhits'};
	my $tolerance     = $PARAMS->{'ms2tolerance'};
	my $picktolerance = $PARAMS->{'picktolerance'};
	my ( %matchlist, %matchlistsorted );
	my %seenpeps = {};
	$verbose && print "Spectrum peaks: ";
	$verbose && print "@peaks\n";
	$verbose && print "ions tested: ";

	#print "assumed parentmass = $testion\n";
	my ( $candidateion, $candidateids );
	my $range = int( $testion * 10 );    #get 0.1 resolution
	my ( $candidates, $id, $rangex, %seenids );

	#get all peptides that contain the specific candidate ion
	foreach $id ( keys %$pepindex )
	{
		my $indexentry = $pepindex->{$id};
		my $desc       = $indexentry->{'desc'};
		my $seq        = $indexentry->{'seq'};
		my $molweight  = $indexentry->{'mw'};
		my $proteinID  = $indexentry->{'proteinID'};
		my $pepobj;

		#take all ions into account that match $testionmass within 1 Da
		unless ( $molweight >= $ms1mass )
		{
			my $pepobj = PepObj->new( $seq, $id, $desc, $MSTAB, $PARAMS, $verbose, $proteinID );
			push @$candidates, $pepobj;
		}
	}
	my $i = 0;
	foreach my $candidate (@$candidates)
	{
		$i++;

		#	print "$i: testing ", $candidate->id, " ", $candidate->seq,		" M: ", $candidate->molweight, "\n";
		#my @candidate_ions = getcandidateions( $candidate, $PARAMS );
		#my $nhits = counthits( \@candidate_ions, \@peaks, $tolerance );
		my $newid = $candidate->id;
		$matchlist{$newid}->{'pepObj'}   = $candidate;
		$matchlist{$newid}->{'matchObj'} = $self;
	}
	my @keys = keys(%matchlist);
	foreach my $pep (@keys)
	{
		push @{ $matchlistsorted{ int( 100 * $matchlist{$pep}->{'pepObj'}->molweight ) } }, $matchlist{$pep}->{'pepObj'};
	}
	return \%matchlist, \%matchlistsorted;
}

sub massmatch
{
	my $self       = shift;
	my $spectrum   = shift;
	my $ionindex   = shift;
	my $pepindex   = shift;
	my $PARAMS     = shift;
	my $MSTAB      = shift;
	my $verbose    = shift;
	my $filehandle = shift;
	my $ms1mass    = $self->getms1mass;
	my ( @peaks, $ion, $peak, $testion );
	my $minhits       = $PARAMS->{'minhits'};
	my $tolerance     = $PARAMS->{'ms2tolerance'};
	my $picktolerance = $PARAMS->{'picktolerance'};
	my ( %matchlist, %matchlistsorted );
	my %seenpeps = {};

	#print "assumed parentmass = $testion\n";
	my ( $candidateion, $candidateids );
	my $range = int( $testion * 10 );    #get 0.1 resolution
	my ( $candidates, $id, $rangex, %seenids );

	#get all peptides that contain the specific candidate ion
	foreach $id ( keys %$pepindex )
	{
		my $indexentry = $pepindex->{$id};
		my $desc       = $indexentry->{'desc'};
		my $seq        = $indexentry->{'seq'};
		my $molweight  = $indexentry->{'mw'};
		my $proteinID  = $indexentry->{'proteinID'};
		my $pepobj;

		#take all ions into account that match $testionmass within 1 Da
		unless ( $molweight >= $ms1mass )
		{
			my $pepobj = PepObj->new( $seq, $id, $desc, $MSTAB, $PARAMS, $verbose, $proteinID );
			push @$candidates, $pepobj;
		}
	}
	my $i = 0;
	print $filehandle "number of candidate peptides to combine in matchall mode: ", $#$candidates + 1, "\n";
	foreach my $candidate (@$candidates)
	{
		$i++;
		my $newid = $candidate->id;
		$matchlist{$newid}->{'pepObj'}   = $candidate;
		$matchlist{$newid}->{'matchObj'} = $self;
	}
	my @keys = keys(%matchlist);
	foreach my $pep (@keys)
	{
		push @{ $matchlistsorted{ int( 100 * $matchlist{$pep}->{'pepObj'}->molweight ) } }, $matchlist{$pep}->{'pepObj'};
	}
	return \%matchlist, \%matchlistsorted;
}

sub DESTROY
{
	my $this = shift;
}

sub getcandidate_common_ions
{
	my $candidate  = shift;
	my $candidates = $candidate->iontag_getpossiblecommonions;
	return $candidates;
}

sub getcandidate_xlink_ions
{
	my $self         = shift;
	my $candidate    = shift;
	my $PARAMS       = shift;
	my $Hatom        = $self->getMSTAB->{'Hatom'}->{'native'};
	my $minioncharge = $self->minioncharge_xlinks;
	my $maxioncharge = $self->maxioncharge_xlinks;
	my $deltaMr      = $self->getms1mass - $candidate->molweight;
	my $xlinkionhash = $candidate->iontag_getpossiblexlinkions( $deltaMr, $minioncharge, $maxioncharge );
	return $xlinkionhash;
}

sub getnhits
{
	my $self = shift;
	return $self->{'nhits'};
}

sub nonred
{
	my $array = shift;
	my %seen  = ();
	my $entry;
	my @unique = ();
	foreach $entry (@$array)
	{
		if ( defined($entry) )
		{
			my $xid = join "", $entry->getuniqueid,;
			push( @unique, $entry ) unless $seen{$xid}++;
		}
	}
	@$array = @unique;
}

sub nonred_list
{
	my $array = shift;
	my %seen  = ();
	my $entry;
	my @unique = ();
	foreach $entry (@$array)
	{
		if ( defined($entry) )
		{
			push( @unique, $entry ) unless $seen{$entry}++;
		}
	}
	@$array = @unique;
}

sub makexlinks_iontag
{
	my $self             = shift;
	my $PARAMS           = shift;
	my $filehandle       = shift;
	my $verbose          = shift;
	my $ms1parentmass    = $self->getms1mass;
	my $xlinkermass      = $PARAMS->{'xlinkermw'};
	my $monolinkmasslist = $PARAMS->{'monolinkmw'};
	my @monolinkmasses   = split /,/, $monolinkmasslist;
	my $ms1tolerance     = $PARAMS->{'ms1tolerance'};
	print $filehandle "MS1_mass: ", $ms1parentmass, "\n";
	print $filehandle "MS1_Mz: ",     $self->getms1Mz,     "\n";
	print $filehandle "MS1_Charge: ", $self->getms1Charge, "\n";
	my ( @hits, $hit, @monolinks, @intralinkmasses, @XLINKS );
	my $hitlist = $self->hitlist;    #get all matching peptides
	print $filehandle "hitlist size ", scalar( keys %{$hitlist} ), "\n";
	print "Generation of iontag xlinks:\n";
	print "Number of peptides in hitlist: " . scalar( keys %{$hitlist} ) . "\n";

	# my $ms1parent = $ms1mass;
	if ( $PARAMS->{'tolerancemeasure'} =~ /^ppm/i )
	{
		$ms1tolerance = $PARAMS->{'ms1tolerance'} * 1e-6 * $ms1parentmass;    #ppm to amu measure
	}
	my ( $minborder, $maxborder );
	if ( $PARAMS->{'ms1tol_minborder'} && $PARAMS->{'ms1tol_maxborder'} )
	{
		if ( $PARAMS->{'tolerancemeasure'} =~ /^ppm/i )
		{
			$minborder = $PARAMS->{'ms1tol_minborder'} * 1e-6 * $ms1parentmass;
			$maxborder = $PARAMS->{'ms1tol_maxborder'} * 1e-6 * $ms1parentmass;
		} else
		{
			$minborder = $PARAMS->{'ms1tol_minborder'};
			$maxborder = $PARAMS->{'ms1tol_maxborder'};
		}
	}
	print $filehandle "precursor masstolerance: $ms1tolerance\n";

	#	print "<br>precursor mass: $ms1parentmass, ms1 tolerance: $ms1tolerance";
	if ( $PARAMS->{'search_intralinks'} )
	{
		#search for loop-links
		foreach my $hitlistmember ( keys %{$hitlist} )
		{    #cycle over all peptides

			#		my $ms1parent = $ms1mass;
			my $hit = getpepObj( $hitlist->{$hitlistmember} );    #retrieve peptide object

			#print "looplink: parentmass=$ms1parentmass hit: ",$hit->seq, " mass: ",$hit->molweight + $xlinkermass,"ms1tolerance: $ms1tolerance empircal delta: ",$hit->molweight + $xlinkermass - $ms1parentmass ,"\n";
			if ( abs( $hit->molweight + $xlinkermass - $ms1parentmass ) <= $ms1tolerance )
			{
				if ( $PARAMS->{'ms1tol_minborder'} && $PARAMS->{'ms1tol_maxborder'} )
				{
					my $errorrel = ( $ms1parentmass - ($hit->molweight + $xlinkermass) );
					my $deltappm = 1e6 * ( $errorrel / $ms1parentmass );

					#print "Relativ precursor error is: $errorrel Da and $deltappm ppm minborder/maxborder is $minborder, $maxborder\n";
					if ( $errorrel < $minborder || $errorrel > $maxborder )
					{

						#print "Will Skip this hit\n";
						next;
					}
				}
				my $topology = $self->getIntraCombinations($hit);
				foreach my $combination ( @{$topology} )
				{
					$verbose && print "<br>evaluating intralinks for ", $hit->id, " topology @$combination\n";
					push @intralinkmasses, LinkObj->new( 'intralink', [$hit], $xlinkermass, $combination, $self, $PARAMS, $self->getMSTAB );
				}
			}
		}
	}
	if ( $PARAMS->{'search_monolinks'} )
	{
		#search for monolinks
		foreach my $hitlistmember ( keys %{$hitlist} )
		{    #cycle over all peptides

			#			my $ms1parent = $ms1mass;
			my $hit = getpepObj( $hitlist->{$hitlistmember} );    #retrieve peptide object
			foreach my $monolinkmass (@monolinkmasses)
			{
				if ( abs( $hit->molweight + $monolinkmass - $ms1parentmass ) <= $ms1tolerance )
				{
				
				if ( $PARAMS->{'ms1tol_minborder'} && $PARAMS->{'ms1tol_maxborder'} )
				{
					my $errorrel = ( $ms1parentmass - ($hit->molweight + $monolinkmass) );
					my $deltappm = 1e6 * ( $errorrel / $ms1parentmass );

					#print "Relativ precursor error is: $errorrel Da and $deltappm ppm minborder/maxborder is $minborder, $maxborder\n";
					if ( $errorrel < $minborder || $errorrel > $maxborder )
					{
						#print "Will Skip this hit\n";
						next;
					}
				}
				
					$verbose && print $hit->id, " with Mr: ", $hit->molweight + $monolinkmass, " matches MS1 mass = $ms1parentmass\n";
					my $topology = $self->getMonoCombinations($hit);
					foreach my $position ( @{$topology} )
					{
						$verbose && print "evaluating monolinks for ", $hit->id, " topology $position\n";
						push @monolinks, LinkObj->new( 'monolink', [$hit], $monolinkmass, [$position], $self, $PARAMS, $self->getMSTAB );
					}
				}
			}
		}
	}
	### Search xlinks
	if ( $PARAMS->{'search_intracrosslinks'} || $PARAMS->{'search_intercrosslinks'} )
	{
		#search intercrosslinks
		my ( $hit1, $hit2, $hitobj1, $hitobj2, $i, $j );
		my @candidatepeps = keys %{$hitlist};
		my $nxlinks       = 0;
		my $ndecoys       = 0;
		my $ntargets      =0;
		my $seen          = {};
		$verbose && print "Number of candidate peptides: $#candidatepeps\n";
		## Iterate over all candidate peptides
		for ( $i = 0 ; $i < $#candidatepeps ; $i++ )
		{
			$hitobj1 = $candidatepeps[$i];
			$hit1    = getpepObj( $hitlist->{$hitobj1} );
			my $hit1mw = $hit1->molweight;
			## SELECT ALL CANDIDATES FROM THE CANDIDATES THAT FIT THE MISSING PRECURSOR MASS
			my $candidatepepmatch = $self->selectxlinkcandidates( $hit1, $ms1tolerance );
			
			## COMBINE THESE WITH THE HIT1
			for ( $j = 0 ; $j <= $#$candidatepepmatch ; $j++ )
			{
				$hit2 = $candidatepepmatch->[$j];
				my $hit2mw       = $hit2->molweight;
				my $xlinkpepmass = $hit1mw + $hit2mw + $xlinkermass;
				
				#print "Candidate1 is: ".$hit1->id."\n";
				#print "Candidate2 is: ".$hit2->id."\n";
				if ( abs( $xlinkpepmass - $ms1parentmass ) <= $ms1tolerance )
				{
					if ( $PARAMS->{'ms1tol_minborder'} && $PARAMS->{'ms1tol_maxborder'} )
					{
						my $errorrel = ( $ms1parentmass - $xlinkpepmass );
						my $deltappm = 1e6 * ( $errorrel / $ms1parentmass );

						#print "Relativ precursor error is: $errorrel Da and $deltappm ppm minborder/maxborder is $minborder, $maxborder\n";
						if ( $errorrel < $minborder || $errorrel > $maxborder )
						{
							#print "Will Skip this hit\n";
							next;
						}
					}
					### Sorting of the hits!
					### If the length is equal then sort alphabetically, otherwise by length
					my @hits;
					if ( ( $hit1->getlength ) == ( $hit2->getlength ) )
					{
						@hits = sort { $a->seq cmp $b->seq } ( $hit1, $hit2 );
					} else
					{
						## sort by length descending
						@hits = sort { $b->getlength <=> $a->getlength } ( $hit1, $hit2 );
					}
					
					### Check if this hit was already seen
					my $id = $hits[0]->seq . "-" . $hits[1]->seq;
																
					if ( $seen->{$id} )
					{
						next;
					} else
					{
						$seen->{$id} = 1;
					}


					####################
					# GENERATION OF THE XLINK HITS
					####################
					my $topology = $self->getXlinkCombinations( $hits[0], $hits[1] );
					### all topologies are put in an array [3,7]
					### filter out those combinations that have the terminal AA as topology if nocutatxlink is on
					foreach my $combination (@$topology)
					{
						#print "Combination: ", $combination->[0], " ", $combination->[1], " Peptides: ",$hits[0]->seq, " ", $hits[1]->seq, "\n";
						## check if one of the two combos has the eastern AA as topo, then check if its allowed to xlink there, otherwise donot make an xlink
						## Also check if the AArequired contains : which indicates pairs of AAs; then check if properly matched
						my $booleanHetero = 0;
						if( scalar(split /:/, $PARAMS->{'AArequired2'}) > 1 ) {
							my $AApairs = [split /,|:/, $PARAMS->{'AArequired2'}];
							for (my $inc = 0; $inc < scalar @{$AApairs}; $inc += 2) {	
								$booleanHetero = $booleanHetero ||
								( (substr $hits[0]->seq, $combination->[0] - 1,1) =~ $AApairs->[$inc] && (substr $hits[1]->seq, $combination->[1] - 1,1) =~ $AApairs->[$inc+1]) ||
								( (substr $hits[1]->seq, $combination->[1] - 1,1) =~ $AApairs->[$inc] && (substr $hits[0]->seq, $combination->[0] - 1,1) =~ $AApairs->[$inc+1]);
							}
						} else {
							$booleanHetero = 1;
						}
						unless ( (( $PARAMS->{'nocutatxlink'} ) && ( ( $hits[0]->getlength == $combination->[0] ) || ( $hits[1]->getlength == $combination->[1] ) ))  || !$booleanHetero )
						{    #test for xlinks on eastern K in first half of unless condition
							my $linkobj = LinkObj->new( 'xlink', \@hits, $xlinkermass, $combination, $self, $PARAMS, $self->getMSTAB );
							$linkobj->set_target_dc_label("target");
							
							my $protstring = join (",", $hits[0]->protidstring, $hits[1]->protidstring);
							
							if ($protstring =~ m/decoy/){
								$ndecoys++;
							}else{
								$ntargets++;
							}
							push @XLINKS, $linkobj;
							$nxlinks++;
							#exit;
						}
					}
				}
			}
		}
		#print "Generated $ntargets target xlinks and $ndecoys decoy xlinks\n";
	}
##################
##################
	#nonred( \@XLINKS );
	my @results = ( @XLINKS, @monolinks, @intralinkmasses );
	print $filehandle "found ", $#XLINKS + 1, " xlinks\n";
	print "found ", $#XLINKS + 1,          " xlinks\n";
	print "found ", $#intralinkmasses + 1, " intralinks\n";
	print "found ", $#monolinks + 1,       " monolinks\n";
	return \@results;
}

sub makexlinks_enumerate
{
	my $self                         = shift;
	my $PARAMS                       = shift;
	my $filehandle                   = shift;
	my $verbose                      = shift;
	my $PEPINDEX                     = shift;
	my $ENUMERATION_INTERXLINKS      = shift;
	my $ENUMERATIONINDEX_INTRAXLINKS = shift;
	my $ENUMERATIONINDEX_MONO        = shift;
	my $ENUMERATIONINDEX_INTRA       = shift;
	my $MSTAB                        = shift;
	my $INFOINDEX                    = shift;
	my $ms1mass                      = $self->getms1mass;
	my $xlinkermass                  = $PARAMS->{'xlinkermw'};
	my $monolinkmass                 = $PARAMS->{'monolinkmw'};
	my @monolinkmasses               = split /,/, $monolinkmass;
	my $tolerance                    = $PARAMS->{'ms1tolerance'};
	my $ms1intprecision              = $PARAMS->{'intprecision'};
	print $filehandle "MS1_mass: ", $ms1mass, "\n";
	print $filehandle "MS1_Mz: ",     $self->getms1Mz,     "\n";
	print $filehandle "MS1_Charge: ", $self->getms1Charge, "\n";
	my ( @candidatepeps_xlink, @candidatepeps_mono, @candidatepeps_intra, @MONOLINKS, @INTRALINKS, @XLINKS );
	my $hitlist   = $self->hitlist;    #get all matching peptides
	my $ms1parent = $ms1mass;

	if ( $PARAMS->{'tolerancemeasure'} =~ /^ppm/i )
	{
		$tolerance = $PARAMS->{'ms1tolerance'} * 1e-6 * $ms1parent;    #ppm to amu measure
	}
	print $filehandle "precursor masstolerance: $tolerance\n";
	my $intmass       = int( $ms1intprecision * $ms1mass );
	my $nbins         = int( $ms1intprecision * $tolerance ) + 1;
	my $id2numberhash = undef;

	#	if($PARAMS->{'search_intercrosslinks'} ){
	#	  $id2numberhash = $ENUMERATION_INTERXLINKS->{'id2num'};    #get number to peptide-id hash
	#	}elsif($PARAMS->{'search_intracrosslinks'} ){
	#		$id2numberhash = $ENUMERATION_INTERXLINKS->{'id2num'};
	#	}
	$id2numberhash = $INFOINDEX->{'id2num'};

	#print Dumper($INFOINDEX);
	#exit;
	unless ($id2numberhash)
	{
		die "Id 2 Number Hash is empty. Rebuild the enum index.\n";
	}
	my ( $hit1, $hit2, $hitobj1, $hitobj2, $i, $j, %seenpeps );
	my @bins = ();
	if ( $PARAMS->{'search_monolinks'} )
	{
		for $i ( $intmass - $nbins .. $intmass + $nbins )
		{
			if ( $ENUMERATIONINDEX_MONO->{$i} )
			{
				push @candidatepeps_mono, @{ $ENUMERATIONINDEX_MONO->{$i} };
			}
		}
		print $filehandle "number of peptide monolinks within ", $PARAMS->{'ms1tolerance'}, " ", $PARAMS->{'tolerancemeasure'}, " = ", $#candidatepeps_mono + 1, "\n";
	}
	if ( $PARAMS->{'search_intralinks'} )
	{
		for $i ( $intmass - $nbins .. $intmass + $nbins )
		{
			if ( $ENUMERATIONINDEX_INTRA->{$i} )
			{
				push @candidatepeps_intra, @{ $ENUMERATIONINDEX_INTRA->{$i} };
			}
		}
		print $filehandle "number of peptide intralinks within ", $PARAMS->{'ms1tolerance'}, " ", $PARAMS->{'tolerancemeasure'}, " = ", $#candidatepeps_intra + 1, "\n";
	}

	#search monolinks
	my %seen = ();
	foreach my $candidate_monolink (@candidatepeps_mono)
	{
		my $id = $id2numberhash->{$candidate_monolink};    #convert number from combinationhash into peptide-id
		unless ( $seen{$id}++ )
		{
			my $indexentry = $PEPINDEX->{$id};
			my $desc       = $indexentry->{'desc'};
			my $seq        = $indexentry->{'seq'};
			my $molweight  = $indexentry->{'mw'};
			my $proteinID  = $indexentry->{'proteinID'};
			my $pepobj;
			if ( defined( $seenpeps{$id} ) )
			{
				$pepobj = $seenpeps{$id};
			} else
			{
				$pepobj = PepObj->new( $seq, $id, $desc, $MSTAB, $PARAMS, $verbose, $proteinID );
				$seenpeps{$id} = $pepobj;
			}
			if ($verbose)
			{
				my $error = 1e6 * abs( $ms1parent - ( $pepobj->molweight + $monolinkmass ) ) / $ms1parent;
				print $pepobj->id, " ", int($error), " ppm\n";
			}
			foreach my $mass (@monolinkmasses)
			{
				my $candidate_mw = $pepobj->molweight + $mass;
				if (    ( $candidate_mw >= $ms1parent - $tolerance )
					 && ( $candidate_mw <= $ms1parent + $tolerance ) )
				{
					my $topology = $self->getMonoCombinations($pepobj);
					foreach my $position ( @{$topology} )
					{
						push @MONOLINKS, LinkObj->new( 'monolink', [$pepobj], $mass, [$position], $self, $PARAMS, $self->getMSTAB );
					}
				}
			}
		}
	}

	#serch intralinks
	%seen = ();
	foreach my $candidate_intralinks (@candidatepeps_intra)
	{
		my $id = $id2numberhash->{$candidate_intralinks};    #convert number from combinationhash into peptide-id
		unless ( $seen{$id}++ )
		{
			my $indexentry = $PEPINDEX->{$id};
			my $desc       = $indexentry->{'desc'};
			my $seq        = $indexentry->{'seq'};
			my $molweight  = $indexentry->{'mw'};
			my $proteinID  = $indexentry->{'proteinID'};
			my $pepobj;
			if ( defined( $seenpeps{$id} ) )
			{
				$pepobj = $seenpeps{$id};
			} else
			{
				$pepobj = PepObj->new( $seq, $id, $desc, $MSTAB, $PARAMS, $verbose, $proteinID );
				$seenpeps{$id} = $pepobj;
			}
			if ($verbose)
			{
				my $error = 1e6 * abs( $ms1parent - ( $pepobj->molweight + $xlinkermass ) ) / $ms1parent;
				print $pepobj->id, " ", int($error), " ppm\n";
			}
			my $candidate_mw = $pepobj->molweight + $xlinkermass;
			if (    ( $candidate_mw >= $ms1parent - $tolerance )
				 && ( $candidate_mw <= $ms1parent + $tolerance ) )
			{
				my $topology = $self->getIntraCombinations($pepobj);
				foreach my $combination ( @{$topology} )
				{
					push @INTRALINKS, LinkObj->new( 'intralink', [$pepobj], $xlinkermass, $combination, $self, $PARAMS, $self->getMSTAB );
				}
			}
		}
	}
	if ( $PARAMS->{'search_intercrosslinks'} )
	{
		for $i ( $intmass - $nbins .. $intmass + $nbins )
		{
			push @bins, $i;
			if ( $ENUMERATION_INTERXLINKS->{$i} )
			{
				push @candidatepeps_xlink, @{ $ENUMERATION_INTERXLINKS->{$i} };
			}
		}
	}
	if ( $PARAMS->{'search_intracrosslinks'} )
	{
		@bins = ();
		for $i ( $intmass - $nbins .. $intmass + $nbins )
		{
			push @bins, $i;
			if ( $ENUMERATIONINDEX_INTRAXLINKS->{$i} )
			{
				push @candidatepeps_xlink, @{ $ENUMERATIONINDEX_INTRAXLINKS->{$i} };
			}
		}
	}
	print $filehandle "number of peptide combinations within hash for bins ";
	foreach my $bin (@bins)
	{
		print $filehandle $bin / $ms1intprecision, " ";
	}
	print $filehandle " = ", scalar(@candidatepeps_xlink), "\n";

	#search xlinks prefilter for mass from index
	my @candidatepeps_xlink_massfiltered = ();
	foreach my $candidatecombination (@candidatepeps_xlink)
	{
		my $hit1 = $candidatecombination->[0];
		my $hit2 = $candidatecombination->[1];
		my $hitmw;

		#print $candidatecombination->[3];
		#exit;
		if ( $candidatecombination->[2] )
		{
			$hitmw = pop @$candidatecombination;

			#print "candidate: ",$hitmw,"\n";
			## check if the hit is within the tolerance
			if ( abs( $ms1parent - $hitmw ) <= $tolerance )
			{
				push @candidatepeps_xlink_massfiltered, $candidatecombination;
			}
		} else
		{
			push @candidatepeps_xlink_massfiltered, $candidatecombination;
		}
	}

	#search only on filtered entries
	print $filehandle "number of filtered candidates with mass difference <= $tolerance amu = ", scalar(@candidatepeps_xlink_massfiltered), "\n";
	foreach my $candidatecombination (@candidatepeps_xlink_massfiltered)
	{
		my $hit1 = $candidatecombination->[0];    ### the index numbers
		my $hit2 = $candidatecombination->[1];
		my @hits = ();

		#print Dumper (@$candidatecombination);
		#exit;
		foreach my $id (@$candidatecombination)
		{
			$id = $id2numberhash->{$id};          #convert number from combinationhash into peptide-id
			my $indexentry = $PEPINDEX->{$id};

			#print Dumper ($PEPINDEX);
			my $desc      = $indexentry->{'desc'};
			my $seq       = $indexentry->{'seq'};
			my $proteinID = $indexentry->{'proteinID'};

			#	my $molweight  = $indexentry->{'mw'};
			my $pepobj;
			if ( defined( $seenpeps{$id} ) )
			{
				$pepobj = $seenpeps{$id};
				push @hits, $pepobj;
			} else
			{
				$pepobj = PepObj->new( $seq, $id, $desc, $MSTAB, $PARAMS, $verbose, $proteinID );
				$seenpeps{$id} = $pepobj;
				push @hits, $pepobj;
			}
		}

		#$verbose=1;
		if ($verbose)
		{
			my $error = 1e6 * abs( $ms1parent - ( $hits[0]->molweight + $hits[1]->molweight + $xlinkermass ) ) / $ms1parent;
			print $hits[0]->id, " ", $hits[1]->id, " ", int($error), " ppm\n";
		}

		#TODO:
		# -Generation of ENUM decoys
		# -Filtering for seenpeptides so that not redundant decoys come in (see IT mode)
		#@hits = sort { $a->seq cmp $b->seq } @hits;
		### Sorting of the hits!
		### If the length is equal then sort alphabetically, otherwise by length
		if ( ( $hits[0]->getlength ) == ( $hits[1]->getlength ) )
		{

			#print "TRUE\n";
			@hits = sort { $a->seq cmp $b->seq } (@hits);
		} else
		{
			## sort by length descending
			@hits = sort { $b->getlength <=> $a->getlength } (@hits);
		}
		## sort by length descending
		#@hits = sort { $b->getlength cmp $a->getlength } @hits;
		my $topology = $self->getXlinkCombinations( $hits[0], $hits[1] );
		foreach my $combination (@$topology)
		{
			my $booleanHetero = 0;
			if( scalar(split /:/, $PARAMS->{'AArequired2'}) > 1 ) {
				my $AApairs = [split /,|:/, $PARAMS->{'AArequired2'}];
				for (my $inc = 0; $inc < scalar @{$AApairs}; $inc += 2) {	
					$booleanHetero = $booleanHetero ||
					( (substr $hits[0]->seq, $combination->[0] - 1,1) =~ $AApairs->[$inc] && (substr $hits[1]->seq, $combination->[1] - 1,1) =~ $AApairs->[$inc+1]) ||
						( (substr $hits[1]->seq, $combination->[1] - 1,1) =~ $AApairs->[$inc] && (substr $hits[0]->seq, $combination->[0] - 1,1) =~ $AApairs->[$inc+1]);
				}
			} else {
				$booleanHetero = 1;
			}
			
			unless ( (( $PARAMS->{'nocutatxlink'} ) && ( ( $hits[0]->getlength == $combination->[0] ) || ( $hits[1]->getlength == $combination->[1] ) ))  || !$booleanHetero )
			{    #test for xlinks on eastern K
				push @XLINKS, LinkObj->new( 'xlink', \@hits, $xlinkermass, $combination, $self, $PARAMS, $MSTAB );
			}
		}
	}
	my @results = ( @XLINKS, @MONOLINKS, @INTRALINKS );
	print $filehandle "found ", scalar(@XLINKS), " possible xlink topologies, ", scalar(@MONOLINKS), " possible mononlinks, ", scalar(@INTRALINKS), " possible intralinks ", "within ", $PARAMS->{'ms1tolerance'}, " ", $PARAMS->{'tolerancemeasure'}, "\n";
	return \@results;
}

sub getParams
{
	my $self = shift;
	return $self->{'PARAMS'};
}

sub verbose
{
	my $self = shift;
	return $self->getParams->{'verbose'};
}

sub selectxlinkcandidates
{
	my $self        = shift;
	my $hit1        = shift;
	my $tolerance   = shift;
	my $hithash     = $self->gethithash;
	my $PARAMS      = $self->getParams;
	my $ms1mass     = $self->getms1mass;
	my $xlinkermass = $PARAMS->{'xlinkermw'};
	my $targetmass  = int( 100 * ( $ms1mass - ( $hit1->molweight + $xlinkermass ) ) );
	my $i;
	my @targets = ();
	my $nbins   = 1 + int( $tolerance * 100 );

	for $i ( ( $targetmass - $nbins ) .. ( $targetmass + $nbins ) )
	{
		if ( defined( $hithash->{$i} ) )
		{
			push @targets, @{ $hithash->{$i} };
		}
	}
	return \@targets;
}

sub getMs1parentmass
{
	my $self    = shift;
	my $hit     = shift;
	my $ms1mass = $self->{$hit}->{'precursorMz'} * $self->{'hit'}->{'precursorcharge'} - $self->{$hit}->{'precursorcharge'} * 1.007825032;
	return $ms1mass;
}

sub getpepObj
{
	my $self = shift;
	return $self->{'pepObj'};
}

sub getXlinkCombinations
{
#Determines the possible combinations for a given pair of peptides (i.e. if there are multiples Ks in each, etc.)
	my $self    = shift;
	my $seq1Obj = shift;
	my $seq2Obj = shift;
	my $seq1    = $seq1Obj->seq;
	my $seq2    = $seq2Obj->seq;
	my ( @hit1pos, @hit2pos, @combinations );
	my $xlinkAA = $self->getParams->{'AArequired'};
	@hit1pos      = ();
	@hit2pos      = ();
	@combinations = ();

	while ( $seq1 =~ /($xlinkAA)/gi )
	{
		push @hit1pos, pos($seq1);

		#@hit1pos=searchposition($seq1,$xlinkAA);
	}
	while ( $seq2 =~ /($xlinkAA)/gi )
	{
		push @hit2pos, pos($seq2);

		#@hit2pos=searchposition($seq2,$xlinkAA);
	}

	#print "$seq1 @hit1pos\n";
	#print "$seq2 @hit2pos\n";
	foreach my $x1 (@hit1pos)
	{
		foreach my $x2 (@hit2pos)
		{
			push @combinations, [ $x1, $x2 ];
		}
	}
	return \@combinations;
}

sub getMonoCombinations
{
	my $self   = shift;
	my $seqObj = shift;
	my $seq    = $seqObj->seq;
	my $length = $seqObj->length;
	my ( @hitpos, @combinations, $x );
	my $xlinkAA       = $self->getParams->{'AArequired'};
	my $nocutatxlinks = $self->nocutatxlinks;
	while ( $seq =~ /($xlinkAA)/gi )
	{

		unless ( $nocutatxlinks && ( $length == pos($seq) ) )
		{
			push @hitpos, pos($seq);
		}
	}
	return \@hitpos;
}

sub searchposition
{
	my $seq        = shift;
	my $AArequired = shift;
	my @positions  = ();
	my @sequence   = split //, $seq;
	my $i          = 0;
	for $i ( 0 .. $#sequence )
	{
		if ( $sequence[$i] =~ /$AArequired/i )
		{
			push @positions, $i + 1;
		}
	}
	return @positions;
}

sub getIntraCombinations
{
	my $self          = shift;
	my $seqObj        = shift;
	my $seq           = $seqObj->seq;
	my $nocutatxlinks = $self->nocutatxlinks;
	my $length        = $seqObj->length;
	my ( @hitpos, @combinations, $x1, $x2 );
	my $xlinkAA = $self->getParams->{'AArequired'};
	while ( $seq =~ /($xlinkAA)/gi )
	{

		unless ( $nocutatxlinks && ( $length == pos($seq) ) )
		{
			push @hitpos, pos($seq);
		}
	}
	for $x1 ( 0 .. $#hitpos )
	{
		for $x2 ( $x1 .. $#hitpos )
		{
			if ( $hitpos[$x2] > $hitpos[$x1] )
			{
				push @combinations, [ $hitpos[$x1], $hitpos[$x2] ];
			}
		}
	}
	return \@combinations;
}
### No multiple matching of peaks to ions was allowed
sub counthits_old
{
	my $self         = shift;
	my $liste1       = shift;
	my $liste2       = shift;
	my $ms2tolerance = shift;              ### Tolerance for common or xlink ms2 peaks // Da or ppm
	my $PARAMS       = $self->getParams;
	## initially set the ms2 tolerance
	## if ppm is selected the tolerance is then recalculated
	my $tolerance = $ms2tolerance;
	my $ppm       = undef;
	if ( $PARAMS->{'tolerancemeasure_ms2'} =~ /^ppm/i )
	{
		$ppm = defined;
	}

	#print "MS2 tolerance: $ms2tolerance\n";
	#push @$liste2,1002.901;
	#my @liste2     = sort { $a <=> $b } @$liste2;
	#print "Liste1:\n";
	#$liste2=\@liste2;
	#print Dumper ($liste1);
	#print "\nListe2:\n";
	#print Dumper ($liste2);
	#print "Tolerance: $tolerance\n";
	#exit;
	my ( $i, $j );

	#nonred_list($liste1);
	#nonred_list($liste2);
	my @list1 = sort { $a <=> $b } @$liste1;
	my @list2 = sort { $a <=> $b } @$liste2;
	$j = 0;
	my $nhits = 0;
	for $i ( 0 .. $#list1 )
	{
		## Adjust tolerance dynamically if ppm measure is selected
		## Tolerance is calculated from the peak perspective
		if ( defined($ppm) )
		{
			$tolerance = $list1[$i] * 1e-6 * $ms2tolerance;
		}
		##    1. Aslong the difference between ION (List2)-PEAK (List1) is not larger than
		##    the tolerance then the ion has a chance to match the peak.
		##    Then the until loop is executed aslong as the condition is not true.
		until ( ( ( $list2[$j] - $list1[$i] ) >= $tolerance ) || $j > $#list2 )
		{

			#print "Tolerance $tolerance ,DELTA ION: ".$list2[$j]." - PEAK: ".$list1[$i].":".($list2[$j]-$list1[$i])."\n";
			## 1.b. Recalculate the error if ppm measure is selected
			#if (defined($ppm)){
			#$tolerance = $list2[$j] * 1e-6 * $ms2tolerance;
			#}
			## 2. Check if the peak matches the ion
			if ( abs( $list1[$i] - $list2[$j] ) <= $tolerance )
			{
				## 2.b if yes count hit
				print "peak ", $list1[$i], " hit ", $list2[$j], " delta " . abs( $list1[$i] - $list2[$j] ) . " tolerance allowed: $tolerance\n";
				$nhits++;
			}
			## 3. Count j up, this is the counter where the list 2 will start
			##    in the next iteration. There it may start at least at the next position of list 2 because the ion in list 1 will be at least equal or larger.
			$j++;
		}
	}
	exit;
	return $nhits;
}
### Added by TW
### Multiple matching of peaks to ions, or ions to peaks is allowed
### added intensity summing
sub counthits
{
	my $self         = shift;
	my $liste1       = shift;              # Peaks
	my $liste2       = shift;              # IONs theoretical
	my $ms2tolerance = shift;              ### Tolerance for common or xlink ms2 peaks // Da or ppm
	my $ionintensityhash=shift;				### ionintensithash mz--> int
	my $verbose = shift;
	my $PARAMS       = $self->getParams;
	## initially set the ms2 tolerance
	## if ppm is selected the tolerance is then recalculated
	my $tolerance = $ms2tolerance;
	my $ppm       = undef;
	if ( $PARAMS->{'tolerancemeasure_ms2'} =~ /^ppm/i )
	{
		$ppm = defined;
	}

	#print "MS2 tolerance: $ms2tolerance\n";
	#push @$liste2,1002.901;
	#my @liste2     = sort { $a <=> $b } @$liste2;
	#print "Liste1:\n";
	#$liste2=\@liste2;
	#print Dumper ($liste1);
	#print "\nListe2:\n";
	#print Dumper ($liste2);
	#print "Tolerance: $tolerance\n";
	#exit;
	my ( $i, $j );

	#nonred_list($liste1);
	#nonred_list($liste2);
	my @list1 = sort { $a <=> $b } @$liste1;
	my @list2 = sort { $a <=> $b } @$liste2;
	$j = 0;
	my $nhits     = 0;
	my $lastindex = 0;
	
	my $intsum=0;
	for $i ( 0 .. $#list2 )    ## Go through the ions in the outer loop
	{
		## Adjust tolerance dynamically to the ion if ppm measure is selected
		## Tolerance is calculated from the ion perspective
		if ( defined($ppm) )
		{
			$tolerance = $list2[$i] * 1e-6 * $ms2tolerance;
		}
		### This is the counter for the inner loop of the peaks
		### For a new ion the counter is set to the lastindex
		if ( $lastindex == 0 )
		{
			$j = 0;
		} else
		{
			$j = $lastindex;
		}
		##    1. Aslong the peaks is not larger than the ion+tol
		##    it has a chance to match and it goes in the loop (or the peak was the last in the list)
		##    Then the until loop is executed aslong as the condition is not true.
		## $verbose &&    print "TESTING Theoretical ION:".$list2[$i]."starting at j: $j\n";
		until ( ( $list1[$j] > ( $list2[$i] + $tolerance ) ) || $j > $#list1 )
		{

			# 2a. If the PEAK is smaller or equal than the ION-tolerance
			# then this peak will never match against a larger ION in the next round
			# therefore lastindex is set to this peaknumber (j), where matching will start in the next round
			# $verbose && print "Testing Peak: $list1[$j]\n";
			if ( ( $list1[$j] <= ( $list2[$i] - $tolerance ) ) )
			{
				$lastindex = $j;

				#print "Tolerance $tolerance,DELTA ION: ".$list2[$i]." - PEAK: ".$list1[$j].":".($list2[$j]-$list1[$i])."\n";
				#$verbose && print "delta: ", mz( $peaks->[$j] ) - mz( $ions->[$i] ), " lastindex =j: " . $lastindex . "\n";
			}
			## 2b. Check if the peak matches the ion
			if ( abs( $list2[$i] - $list1[$j] ) <= $tolerance )
			{
				## 2.b if yes count hit
				# $verbose && print "peak matches: Peak ",$list1[$j]," --> ION: ",$list2[$i], " delta ".abs( $list2[$i] - $list1[$j] )." tolerance allowed: $tolerance\n";
				$nhits++;
				## add the intensity to the intensity sum
				$intsum+=$ionintensityhash->{$list1[$j]};
				#$verbose && print "Matched peak in spectrum: $list1[$j] intensity: $ionintensityhash->{$list1[$j]}\n";
			
#			unless ($ionintensityhash->{$list1[$j]}){
#				warn "NO PEAK FOUND IN INTHASH!\n";
#			}
			
			
			}
			## 2c. Count j up, to test the next peak with the same ion
			##    in the next iteration.
			$j++;
		}
	}

	#exit;
	return $nhits, $intsum;
}

sub counthits2
{
	my $liste1    = shift;
	my $liste2    = shift;
	my $tolerance = shift;
	my ( $i, $j );
	my @list1 = sort { $a <=> $b } @$liste1;
	my @list2 = sort { $a <=> $b } @$liste2;
	$j = 0;
	my $nhits = 0;

	for $i ( 0 .. $#list1 )
	{
		until ( ( ( $list2[$j] - $list1[$i] ) > $tolerance ) || $j >= $#list2 )
		{
			if ( abs( $list1[$i] - $list2[$j] ) <= $tolerance )
			{
				print "$i peak ", $list1[$i], " hit $j ", $list2[$j], "\n";
				$nhits++;
			}
			$j++;
		}
	}
	return $nhits;
}
### Changed by TW
sub matchpeaks
{
	my $ionsref      = shift;
	my $peaksref     = shift;
	my $PARAMS       = shift;
	my $ms2tolerance = shift;
	my $matchhash    = shift;
	my $delta        = shift;    #used for matching second isotope
	my $errorhashppm = shift;
	my ( $i, @matchedpeaks );
	my $ppmmeasure = 0;

	if ( $PARAMS->{'tolerancemeasure_ms2'} =~ /^ppm/i )
	{
		$ppmmeasure = 1;
	}
	my $tolerance = $ms2tolerance;
	my %seen      = ();

	#	push @$peaksref,968.51;
	my @peaks = sort { $a <=> $b } @$peaksref;
	my @ions  = sort { $a <=> $b } @$ionsref;
	my $j     = 0;
	my $nhits = 0;
	my @matched   = ();
	my $lastindex = 0;
	for $i ( 0 .. $#peaks )
	{
		$j = $lastindex;
		### recal the tolerance to amu units!
		if ($ppmmeasure)
		{

			#print $ms2tolerance;
			$tolerance = $ms2tolerance * 1e-6 * $peaks[$i];    #ppm to amu measure
		}

		#print "Tolerance: $tolerance<br>\n";
		until ( ( $ions[$j] - ( $peaks[$i] ) > $tolerance ) || $j > $#ions )
		{
			if ( abs( $ions[$j] - $peaks[$i] ) <= $tolerance )
			{
				my $delta2 = $peaks[$i] - $ions[$j];
				## calc the error in ppm
				my $deltappm = 1e6 * ( $delta2 / $ions[$j] );

				#print "Matched ion-peak:  $ions[$j] - $peaks[$i] ,tolerance: $tolerance, delta: $delta2 deltappm $deltappm<br>\n";
				#print "Index Peak is: $i, Index ion is: $j<br>";
				### Check if the ion was already seen
				### Then dont report it twice
				### Modified by TW
				### 1. Checks if this ion has already been matched by a peak
				### 2.a. If not, then the peak is added
				### 2.b. If the ion has already a match, then it is checked if this peak has already been assigned
				###      to this ion (this can happen if alpha and beta chain have same ions)
				###		 If not then the peak is added
				unless ( $seen{ $ions[$j] }++ )
				{
					$nhits++;
					push @matched,      $ions[$j];
					push @matchedpeaks, $peaks[$i];
					push @{ $matchhash->{ $ions[$j] } },    $peaks[$i];
					push @{ $errorhashppm->{ $ions[$j] } }, $deltappm;
				} else
				{
					## check if this peak is already seen or if it is another peak!
					unless ( { map { $_ => 1 } @{ $matchhash->{ $ions[$j] } } }->{ $peaks[$i] } )
					{
						push @matchedpeaks, $peaks[$i];
						push @{ $matchhash->{ $ions[$j] } },    $peaks[$i];
						push @{ $errorhashppm->{ $ions[$j] } }, $deltappm;
					}
				}
			}
			if ( ( $ions[$j] <= ( $peaks[$i] - $tolerance ) ) )
			{
				$lastindex = $j;
			}

			#
			#			if ( ( $ions[$j] - ( $peaks[$i] ) < $tolerance ) ) {
			#				$lastindex++;
			#			}
			$j++;
		}
	}
	### DELTA IS USED FOR ONLY SECOND ISOTOPEPAIR MATCHING
	if ($delta)
	{
		$j         = 0;
		$lastindex = 0;
		for $i ( 0 .. $#peaks )
		{
			$j = $lastindex;
			if ($ppmmeasure)
			{
				$tolerance = $ms2tolerance * 1e-6 * $peaks[$i];    #ppm to amu measure
			}
			until ( ( ( $ions[$j] + $delta - $peaks[$i] ) > $tolerance ) || $j > $#ions )
			{
				if ( abs( $ions[$j] + $delta - $peaks[$i] ) <= $tolerance )
				{
					my $delta2 = $ions[$j] - $peaks[$i] + $delta;
					my $deltappm = $delta2 / ( $ions[$j] * 1e-6 );

					#print "Matched 2nd isotope peak: delta: $delta, ion-peak:  $ions[$j] - $peaks[$i] ,tolerance: $tolerance, delta: $delta2<br>\n";
					## 1. Check the ion (the original ion) has already a hit
					## 2.a. If not then assign the peak to this ion
					## 2.b. If it has already a peak assigned check if it it the same peak, only add if not
					unless ( $seen{ $ions[$j] }++ )
					{
						$nhits++;
						push @matched,      $ions[$j];
						push @matchedpeaks, $peaks[$i];
						push @{ $matchhash->{ $ions[$j] } },    $peaks[$i];
						push @{ $errorhashppm->{ $ions[$j] } }, $deltappm;

						#print "Matched 2nd isotope peak: delta: $delta, ion-peak:  $ions[$j] - $peaks[$i] ,tolerance: $tolerance, delta: $delta2<br>\n";
					} else
					{
						### check if this peak has already seen for this ion (from the orginal peak matching)
						unless ( { map { $_ => 1 } @{ $matchhash->{ $ions[$j] } } }->{ $peaks[$i] } )
						{
							push @matchedpeaks, $peaks[$i];
							push @{ $matchhash->{ $ions[$j] } },    $peaks[$i];
							push @{ $errorhashppm->{ $ions[$j] } }, $deltappm;

							#print "THE $peaks[$i] PEAK has been seen\n";
							#next;
						}
					}
				}
				if ( ( $ions[$j] + $delta ) <= ( $peaks[$i] - $tolerance ) )
				{
					$lastindex = $j;
				}
				$j++;
			}
		}
	}
	return \@matched, \@matchedpeaks, $matchhash;
}

sub matchpeaks_modif
{
	my $ionsref      = shift;
	my $peaksref     = shift;
	my $PARAMS       = shift;
	my $ms2tolerance = shift;
	my $matchhash    = shift;
	my $delta        = shift;
	my $specobj      = shift;
	my ( $i, @matchedpeaks );
	my $ppmmeasure = 0;

	if ( $PARAMS->{'tolerancemeasure_ms2'} =~ /^ppm/i )
	{
		$ppmmeasure = 1;
	}
	my $tolerance = $ms2tolerance;
	my %seen      = ();
	my @peaks     = sort { $a <=> $b } @$peaksref;
	my @ions      = sort { $a <=> $b } @$ionsref;
	my $j         = 0;
	my $nhits     = 0;
	my @matched   = ();
	my $lastindex = 0;
	for $i ( 0 .. $#peaks )
	{

		#print $peaks[$i];
		my $ionintensity = $specobj->get_ionintensity( $peaks[$i] );
		$j = $lastindex;
		if ($ppmmeasure)
		{
			$tolerance = $ms2tolerance * 1e-6 * $peaks[$i];    #ppm to amu measure
		}
		until ( ( $ions[$j] - ( $peaks[$i] ) > $tolerance ) || $j > $#ions )
		{
			if ( abs( $ions[$j] - $peaks[$i] ) <= $tolerance )
			{
				if ( $seen{ $ions[$j] } )
				{
					if ( $seen{ $ions[$j] } < $ionintensity )
					{
						pop @matched;
						pop @matchedpeaks;
						push @matched,      $ions[$j];
						push @matchedpeaks, $peaks[$i];
						$matchhash->{ $ions[$j] } = $peaks[$i];
						$seen{ $ions[$j] } = $ionintensity;
					}
				} else
				{
					$nhits++;
					push @matched,      $ions[$j];
					push @matchedpeaks, $peaks[$i];
					$matchhash->{ $ions[$j] } = $peaks[$i];
					$seen{ $ions[$j] } = $ionintensity;
				}
			}
			if ( ( $ions[$j] - ( $peaks[$i] ) < $tolerance ) )
			{
				$lastindex++;
			}
			$j++;
		}
	}
	if ($delta)
	{
		$j         = 0;
		$lastindex = 0;
		for $i ( 0 .. $#peaks )
		{
			$j = $lastindex;
			if ($ppmmeasure)
			{
				$tolerance = $ms2tolerance * 1e-6 * $peaks[$i];    #ppm to amu measure
			}
			until ( ( ( $ions[$j] + $delta - $peaks[$i] ) > $tolerance ) || $j >= $#ions )
			{
				if ( abs( $ions[$j] + $delta - $peaks[$i] ) <= $tolerance )
				{
					unless ( $seen{ $ions[$j] }++ )
					{
						$nhits++;
						push @matched,      $ions[$j];
						push @matchedpeaks, $peaks[$i];
						$matchhash->{ $ions[$j] } = $peaks[$i];
					}
				}
				if ( ( $ions[$j] + $delta - ( $peaks[$i] ) < $tolerance ) )
				{
					$lastindex++;
				}
				$j++;
			}
		}
	}
	return \@matched, \@matchedpeaks, $matchhash;
}

sub matchpeaks_2ndisotopewewe
{
	my $ionsref            = shift;
	my $peaksref           = shift;
	my $PARAMS             = shift;
	my $ms2tolerance       = shift;
	my $matchhash          = shift;
	my $charge             = shift;
	my $secondisotopeshift = 1 / $charge;
	my ( $i, @matchedpeaks );
	my $ppmmeasure = 0;
	my %seen       = ();

	if ( $PARAMS->{'tolerancemeasure_ms2'} =~ /^ppm/i )
	{
		$ppmmeasure = 1;
	}
	my $tolerance = $ms2tolerance;
	my @ions      = sort { $a <=> $b } @$ionsref;
	my @peaks     = sort { $a <=> $b } @$peaksref;
	my $j         = 0;
	my $nhits     = 0;
	my @matched   = ();
	for $i ( 0 .. $#ions )
	{

		if ($ppmmeasure)
		{
			$tolerance = $ms2tolerance * 1e-6 * $ions[$i];    #ppm to amu measure
		}
		until ( ( ( $peaks[$j] - $ions[$i] ) > ( $tolerance + $secondisotopeshift ) ) || $j >= $#peaks )
		{
			if ( ( abs( $ions[$i] - $peaks[$j] ) <= $tolerance ) )
			{
				unless ( $seen{ $ions[$i] }++ )
				{
					$nhits++;
					push @matched,      $ions[$i];
					push @matchedpeaks, $peaks[$j];
					$matchhash->{ $ions[$i] } = $peaks[$j];
				}
			} elsif ( abs( $ions[$i] - ( $peaks[$j] - $secondisotopeshift ) ) <= $tolerance )
			{
				unless ( $seen{ $ions[$i] }++ )
				{
					$nhits++;
					push @matched,      $ions[$i];
					push @matchedpeaks, $peaks[$j];
					$matchhash->{ $ions[$i] } = $peaks[$j];
				}
			}
			$j++;
		}
	}
	return \@matched, \@matchedpeaks, $matchhash;
}

sub calcxcorr
{
	my $series1  = shift;
	my $series2  = shift;
	my $maxdelay = shift;
	my $PARAMS   = shift;
	my $dotproduct;
	if ( $PARAMS->{'xcorrdelay'} == 0 )
	{
		$dotproduct = defined;
	}
	my @x = @{$series1};
	my @y = @{$series2};
	my $n = $#$series1;
	my ( $i, $j, $delay, $mx, $my, $sx, $sy, $sxy, $denom, @r );

	#   /* Calculate the mean of the two series x[], y[] */
	$my = 0;
	$my = 0;
	for ( $i = 0 ; $i <= $n ; $i++ )
	{
		$mx += $x[$i];
		$my += $y[$i];
	}

	#print "n $n mean x $mx mean y =$my\n";
	$mx /= ($n);
	$my /= ($n);

	#   /* Calculate the denominator */
	$sx = 0;
	$sy = 0;
	for ( $i = 0 ; $i <= $n ; $i++ )
	{
		$sx += ( $x[$i] - $mx ) * ( $x[$i] - $mx );
		$sy += ( $y[$i] - $my ) * ( $y[$i] - $my );
	}
	$denom = sqrt( $sx * $sy );

	#   /* Calculate the correlation series */
	for ( $delay = -$maxdelay ; $delay <= $maxdelay ; $delay++ )
	{
		$sxy = 0;
		for ( $i = 0 ; $i < $n ; $i++ )
		{
			$j = $i + $delay;
			unless ( $j < 0 || $j >= $n )
			{
				$sxy += ( $x[$i] - $mx ) * ( $y[$j] - $my );
			} else
			{
				next;
			}

			#         /* Or should it be (?)
			#			if ( $j < 0 || $j >= $n ) { next; }
			#			else {
			#				$sxy += ( $x[$i] - $mx ) * ( $y[$j] - $my );
			#
			#				#$sxy += ($x[$i] - $mx) * (-$my);
			#				#$sxy += ( $x[$i] - $mx ) * ( $y[$j] - $my );
			#			}
			#*/;
		}
		if ( $denom != 0 )
		{
			push @r, [ $delay, ( $sxy / $denom ) ]    #      /* r is the correlation coefficient at "delay" */
			                                          #print "$delay ",$sxy / $denom,"\n";
		} else
		{
			push @r, [ $delay, 0 ]                    #      /* r is the correlation coefficient at "delay" */
			                                          #print "$delay ",$sxy / $denom,"\n";
		}
	}
	return \@r;
}

sub round
{
	my ($number) = shift;
	return int( $number + .5 );
}

sub maxioncharge_xlinks
{
	my $self                = shift;
	my $definedmaxioncharge = $self->getParams->{'ioncharge_xlink'}->[-1];
	my $spectrumcharge      = $self->getSpecObj->getprecursorCharge;
	if ( $definedmaxioncharge > $self->getSpecObj->getprecursorCharge )
	{
		return $spectrumcharge;
	} else
	{
		return $definedmaxioncharge;
	}
}

sub minioncharge_xlinks
{
	my $self = shift;
	return $self->getParams->{'ioncharge_xlink'}->[0];
}

sub maxioncharge_common
{
	my $self                = shift;
	my $definedmaxioncharge = $self->getParams->{'ioncharge_common'}->[-1];
	my $spectrumcharge      = $self->getSpecObj->getprecursorCharge;
	if ( $definedmaxioncharge > $self->getSpecObj->getprecursorCharge )
	{
		return $spectrumcharge;
	} else
	{
		return $definedmaxioncharge;
	}
}

sub minioncharge_common
{
	my $self = shift;
	return $self->getParams->{'ioncharge_common'}->[0];
}

sub nocutatxlinks
{
	my $self = shift;
	return $self->getParams->{'nocutatxlinks'};
}

sub getSpecObj
{
	my $self = shift;
	return $self->{'spectrum'};
}

sub getspectrum
{
	my $self = shift;
	return $self->{'spectrum'};
}
1;
