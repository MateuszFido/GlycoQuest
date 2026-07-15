#!/usr/bin/env perl
use strict;

#---------------------------------------------------------------------------
# xdecoy.pl
# A software/script to generate decoy sequences.
# Execute xdecoy.pl -help to display information and usage options.
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
use Bio::Perl;
use Math::Random;
use File::Basename;
use Getopt::Long;
use File::Spec;
use Data::Dumper;

my ( $verbose,$shuffle,$descstring, $sv, $help, $input, $reverse, $reverse_internal_pep, $targetseq, $outfile, $dcstring, $shuffleids, $decoy_mode, $binsize, $foreigndb );

my %scramblefunctions = (
						  'nterm'    => \&_scrambleprot_keep_nterm,
						  'cterm'    => \&_scrambleprot_keep_cterm,
						  'ctermpep' => \&_reversepep_keep_cterm,
						  'ntermpep' => \&_reversepep_keep_nterm,
						  'reverse'  => \&_reverse_sequence,
);

my $scriptinfo = {};
$scriptinfo->{'version'} = "1.3";
$scriptinfo->{'author'}  = "Thomas Walzthoeni based on orginal work by Oliver Rinner";
my $enzymenum = 1;


#---------------------------------------------------------------------------
# Standard settings (are overwritten if used as params)  
#---------------------------------------------------------------------------
$decoy_mode = "reverse";
#$dcstring   = "decoy";
$binsize    = 50;

my %PARAMS = (
			   'db'          => $input,
			   'out'         => $outfile,
			   'enz'         => $enzymenum,
			   'decoy'       => $decoy_mode,
			   'verbose'     => $verbose,
			   'help'        => $help,
			   'targetseq'   => $targetseq,
			   'decoystring' => $dcstring,
			   'shuffleids'  => $shuffleids,
			   'binsize'     => $binsize,
			   'sverbose'    => $sv,
			   'foreigndb'   => $foreigndb,
			   'shuffle'=>$shuffle,
			   'descstring'=>$descstring,
);
## HERE ALL PARAMS THAT COME FROM THE CMD LINE MUST BE DEFINED
GetOptions( \%PARAMS, 'db=s', 'out=s', 'decoy=s', 'verbose', 'help', 'targetseq', 'decoystring=s', 'shuffleids', 'binsize=i', 'sverbose', 'foreigndb=s', 'shuffle', 'descstring=s' );
my $PARAMS = \%PARAMS;

if ( $PARAMS->{'verbose'} || $PARAMS->{'sverbose'} )
{
	print_params($PARAMS);
}
&usage() unless ( $PARAMS->{'db'} );
&usage() if $PARAMS->{'help'};

#---------------------------------------------------------------------------
#  Set some global variables
#---------------------------------------------------------------------------
$decoy_mode = $PARAMS->{'decoy'};
$dcstring = $PARAMS->{'decoystring'};
$descstring = $PARAMS->{'descstring'};
$outfile = $PARAMS->{'out'};

unless ($dcstring){
if ($decoy_mode eq "reverse"){
$dcstring = "reverse";	
}else{
$dcstring = "decoy";		
}
}

unless ($descstring){
if ($decoy_mode eq "reverse"){
$descstring = "Reverse";	
}else{
$descstring = "Decoy";		
}
} 


( $PARAMS->{'verbose'} || $PARAMS->{'sverbose'} ) && print "Decoymode: $decoy_mode\n";
&usage() unless ( $decoy_mode eq "cterm" || $decoy_mode eq "nterm" || $decoy_mode eq "ctermpep" || $decoy_mode eq "ntermpep" || $decoy_mode eq "reverse" );
my $bndb = $PARAMS->{'db'};
#$bndb =~ s/\.\w+//; #Doesn't work on paths with . in the name
my ($bndb_file, $bndb_dir) = File::Basename::fileparse($bndb);
$bndb = $bndb_dir . $bndb_file;
unless ($outfile)
{
	$outfile = $bndb . "_decoy.fasta";
}
$PARAMS->{'db_bn'}=$bndb;


print "outfilename: $outfile\n";
#exit;
my $indb  = Bio::SeqIO->new( -file => "$PARAMS->{'db'}", -format => 'fasta' );
my $outdb = Bio::SeqIO->new( -file => ">$outfile",       -format => 'fasta' );

#---------------------------------------------------------------------------
#  GENERATION of a Forein DB if this option is selected
#  A forein database is used to mimic the actual target database
#---------------------------------------------------------------------------
my ( $foreigndbin, $foreignseqs, $foreigngrouped, $targetgrouped );
my $foreigndbinhash;

if ( $PARAMS->{'foreigndb'} )
{
	( $PARAMS->{'verbose'} || $PARAMS->{'sverbose'} ) && print "Reading foreign database:  $PARAMS->{'foreigndb'}\n";
	$foreigndbin = Bio::SeqIO->new( -file => "$PARAMS->{'foreigndb'}", -format => 'fasta' );
	$foreignseqs = read_db_into_hash($foreigndbin);
	my $longest = get_longest_prot($foreignseqs);
	$foreigngrouped = gen_group_hash( $foreignseqs, $longest );
	( $PARAMS->{'verbose'} || $PARAMS->{'sverbose'} ) && print "Reading template database:  $PARAMS->{'db'}\n";
	my $targetdbseqs = read_db_into_hash($indb);
	my $longest2     = get_longest_prot($targetdbseqs);
	$targetgrouped = gen_group_hash( $targetdbseqs, $longest2 );
	my $foreigndecoys=create_db_based_on_templtegrouphash( $targetdbseqs, $targetgrouped, $foreignseqs, $foreigngrouped );
	my $outfilename=$PARAMS->{'db_bn'}."_foreign_decoys.fasta";
	store_db($foreigndecoys,$outfilename);
	exit 0;
	
}
my $ruid      = 0;
my $sequences = {};

#---------------------------------------------------------------------------
#  Generation of a normal decoy database
#  Options are nterm, cterm, ntermpep,ctermpep, [reverse]
#---------------------------------------------------------------------------
while ( my $seqobj = $indb->next_seq() )
{
	my $id       = $seqobj->id();
	my $desc     = $seqobj->desc();
	my $sequence = $seqobj->seq;
	( $sequence = $sequence ) =~ s/\*//g;
	( $PARAMS->{'sverbose'} ) && print "$id $desc\n";
	$PARAMS->{'sverbose'} && print $sequence, "\n";
	my $cuts         = undef;
	my $scrambledseq = "";

	unless ( $decoy_mode eq "reverse" )
	{
		$cuts = getpeps( $sequence, $enzymenum );
		$scrambledseq = &{ $scramblefunctions{$decoy_mode} }( $seqobj, $cuts );
	} else
	{
		$PARAMS->{'sverbose'} && print "Generate reverse Sequence of $id\n";
		$scrambledseq = reverse($sequence);
	}
	( $PARAMS->{'sverbose'} ) && print ">$dcstring" . "_$id $descstring of $desc\n";
	$PARAMS->{'sverbose'} && print "$scrambledseq\n";
	my $newseq = Bio::Seq->new(
								-seq  => $scrambledseq,
								-desc => "$descstring of $desc",
								-id   => "$dcstring" . "_$id",
	);
	if ($targetseq)
	{
		print "Adding also traget proteins to the database\n";
		$sequences->{$ruid}->{'seqobject'} = $seqobj;
		$sequences->{$ruid}->{'length'}    = length($sequence);
		$id++;
	}
	$sequences->{$ruid}->{'seqobject'} = $newseq;
	$sequences->{$ruid}->{'length'}    = length($scrambledseq);
	$ruid++;
}
( $PARAMS->{'verbose'} || $PARAMS->{'sverbose'} ) && print "Indexed ", $ruid , " protein(s)\n";

#---------------------------------------------------------------------------
#  Generate shuffled ID Decoys if selected
#  If selected this approach uses the sequences generated by the decoy approach
#  shifts the identifiers
#---------------------------------------------------------------------------
my $shifteddecoys = {};
if ( $PARAMS->{'shuffleids'} )
{
### get the longest sequence
	my $longest = get_longest_prot($sequences);

	#$PARAMS->{'longest'} = $longest;
	( $PARAMS->{'verbose'} || $PARAMS->{'sverbose'} ) && print "Largest protein contains $longest AAs\n";
### Generate a grouped hash
	my $grouphash = generate_grouped_hash_join_single_bins( $sequences, $longest );
### Shift the ids of the proteins within the bins
	my $newseqhash = {};
	( $PARAMS->{'verbose'} ) && print "Shifting ProteinIds within BINS\n";

	#print Dumper ($grouphash);
	my $ruid = 0;
	foreach my $bin ( keys %$grouphash )
	{
## Get all protein Ids
		my $proteinids = $grouphash->{$bin};

		#print "BIN $bin:",join (",",@$proteinids),"\n";
## Create a shifted array where the last element is the first
		my @idarray      = @{ $grouphash->{$bin} };
		my @shiftedarray = @{ $grouphash->{$bin} };
		unshift( @shiftedarray, pop(@shiftedarray) );

		#print "Shifted array: ",join (",",@shiftedarray),"\n";
### These prot ids are now assigned to the orinal ids
## modify the description to see which protein was shifted with which
		for ( my $i = 0 ; $i < scalar(@idarray) ; $i++ )
		{

			#print "CREATE DC\n";
### take the id out and modify the id
			my $id       = $idarray[$i];
			my $pseudoid = $shiftedarray[$i];
			my $idshift  = $shiftedarray[$i];
			my $seqobj1  = $sequences->{$id}->{'seqobject'};
			my $seqobj2  = $sequences->{$idshift}->{'seqobject'};
### Generate the new seqobject with the sequence of id1 and the identifier of id2
			my $sequence  = $seqobj1;
			my $id1       = $seqobj1->id();
			my $desc1     = $seqobj1->desc();
			my $sequence1 = $seqobj1->seq;
			my $id2       = $seqobj2->id();
			my $desc2     = $seqobj2->desc();
### Newseqobject
			my $desc       = "$desc2 with sequence of $desc1";
			my $identifier = "$id2";
			my $newseq = Bio::Seq->new(
										-seq  => $sequence1,
										-desc => $desc,
										-id   => $identifier,
			);
			( $PARAMS->{'sverbose'} ) && print ">$identifier" . " $desc\n";
#### Add the newsequence to the dchash
			$shifteddecoys->{$ruid}->{'seqobject'} = $newseq;
			$ruid++;
		}
	}
	( $PARAMS->{'verbose'} || $PARAMS->{'sverbose'} ) && print "Created  ", $ruid, " decoys with shifted identifiers\n";
}

#---------------------------------------------------------------------------
#  Write the results for the normal decoy approach
#---------------------------------------------------------------------------
my $ruids = 0;
foreach my $id ( keys %$sequences )
{
	my $seqobj = $sequences->{$id}->{'seqobject'};
	$outdb->write_seq($seqobj);
	$ruids++;
}
print "Wrote $ruids proteins to $outfile\n";

#---------------------------------------------------------------------------
#  Write the resultdatabase if shuffeld ids were used
#---------------------------------------------------------------------------
if ( $PARAMS->{'shuffleids'} )
{
	$outfile = $bndb . "_xdc_shifted.fasta";
	my $outdb = Bio::SeqIO->new( -file => ">$outfile", -format => 'fasta' );
	my $ruids = 0;
	foreach my $id ( keys %$shifteddecoys )
	{
		my $seqobj = $shifteddecoys->{$id}->{'seqobject'};
		$outdb->write_seq($seqobj);
		$ruids++;
	}
	print "Wrote $ruids proteins to $outfile\n";
}



#---------------------------------------------------------------------------
#  FUNCTIONS
#---------------------------------------------------------------------------

sub store_db{
my $seq=shift;
my $dbfilename=shift;
my $outdb = Bio::SeqIO->new( -file => ">$dbfilename", -format => 'fasta' );
my $ruids = 1;
foreach my $id ( keys %$seq )
{
	my $seqobj = $seq->{$id}->{'seqobject'};
	$outdb->write_seq($seqobj);
	$ruids++;
}
print "Wrote $ruids proteins to $dbfilename\n";
}


### foreigndb: create a database sampled from the foreign database
### based on the grouphash of the original database
sub create_db_based_on_templtegrouphash
{
	my $templatedbseqs     = shift;
	my $templategrouphash  = shift;
	my $foreigndbseqs      = shift;
	my $foreigndbgrouphash = shift;
	my $foreigndecoys;
	my $ruid = 0;
	( $PARAMS->{'verbose'} ) && print "Generate foreign decoys\n";
	( $PARAMS->{'verbose'} && $PARAMS->{'shuffle'} ) && print "Will shuffle the foreign bins\n";
	
	foreach my $bin ( sort { $a <=> $b } keys %$templategrouphash )
	{	
	#	print "BIN IS: $bin\n";
		($PARAMS->{'sverbose'} ) && print "Generating foreing decoys for BIN:  $bin\n";
	
		my $templateseqarray = $templategrouphash->{$bin};
### here we can shuffle if we want
		
		my $foreignseqarray = $foreigndbgrouphash->{$bin};
		#print Dumper ($foreignseqarray);
		if ( $PARAMS->{'shuffle'}){
		fisher_yates_shuffle( $foreignseqarray );  
		}  
		#print Dumper ($foreignseqarray);
		
		unless ($foreignseqarray)
		{
			die "Bin $bin does not exist in foreign database\n";
		}
		if ( scalar(@$templateseqarray) > scalar(@$foreignseqarray) )
		{
			die "The bin $bin contains too few sequences in comparision to the template database\n";
		}

		#print Dumper ($templateseqarray);
		for ( my $i = 0 ; $i < scalar(@$templateseqarray) ; $i++ )
		{
### get the id of the templatearray
			my $tempid    = $templateseqarray->[$i];
			my $foreignid = $foreignseqarray->[$i];
			#print "sample from database: tempid $tempid , forId $foreignid\n";
### get the corresponding array from the foreign db
			my $templateseqobj = $templatedbseqs->{$tempid}->{'seqobject'};
			my $foreignseqobj  = $foreigndbseqs->{$foreignid}->{'seqobject'};
### Generate the id and desc of the template protein
			my $id1   = $templateseqobj->id();
			my $desc1 = $templateseqobj->desc();
### Generate the id and desc of the foreign protein
			my $sequence2 = $foreignseqobj->seq;

			#my $id2       = $foreignseqobj->id();
			my $desc2 = $foreignseqobj->desc();
### Newseqobject
			my $desc       = "foreign decoy of $desc1 with sequence of $desc2";
			my $identifier = "$dcstring" . "_$id1";
			my $newseq = Bio::Seq->new(
										-seq  => $sequence2,
										-desc => $desc,
										-id   => $identifier,
			);
			( $PARAMS->{'sverbose'} ) && print ">$identifier" . " $desc\n";
#### Add the newsequence to the dchash
			$foreigndecoys->{$ruid}->{'seqobject'} = $newseq;
			$ruid++;
		}
	}
	return $foreigndecoys;
}










sub get_longest_prot
{
	my $sequencehash = shift;
	my $longest;
	foreach my $id ( keys %$sequencehash )
	{
		$longest = $sequencehash->{$id}->{'length'} if $sequencehash->{$id}->{'length'} > $longest;
	}
	return $longest;
}

sub read_db_into_hash
{
	my $indb      = shift;
	my $ruid      = 0;
	my $sequences = {};
	while ( my $seqobj = $indb->next_seq() )
	{
		$sequences->{$ruid}->{'seqobject'} = $seqobj;
		my $sequence = $seqobj->seq;
		$sequences->{$ruid}->{'length'} = length($sequence);
		$ruid++;
	}
	( $PARAMS->{'verbose'} || $PARAMS->{'sverbose'} ) && print "Indexed ", $ruid + 1, " protein(s)\n";
	return $sequences;
}

sub gen_group_hash
{
	my $seqhash        = shift;
	my $longestprotein = shift;
	my $grouphash      = {};
	my @groupkeys;
### Genererate the BINS
### The array for the score bins
	my $numbins = int( $longestprotein / $PARAMS->{'binsize'} );
	my $rest    = $longestprotein % $PARAMS->{'binsize'};
	if ($rest) { $numbins++ }
	( $PARAMS->{'verbose'} ) && print "Generating grouped hash with $numbins bins, binsize $PARAMS->{'binsize'} AA\n";
### Generate an array with the bins
	my @groupsarray;
	my $indexedprots = 0;

	foreach my $i ( 0 .. $numbins )
	{
		push @groupsarray, $i * $PARAMS->{'binsize'};
		push @groupkeys,   $i;
	}
### Iter through the bins and assign the protein ids to the grouphash
	for ( my $i = 0 ; $i < scalar(@groupsarray) ; $i++ )
	{
		foreach my $id ( keys %$seqhash )
		{
			my $length = $seqhash->{$id}->{'length'};
			if ( ( $length >= $groupsarray[$i] ) && ( $length < $groupsarray[ $i + 1 ] ) )
			{
				push @{ $grouphash->{ $groupsarray[$i] } }, $id;
				$indexedprots++;
			}
		}
	}
	( $PARAMS->{'verbose'} ) && print "Indexed ".($indexedprots+1)." Proteins into ProteingroupHash\n";
	### Print Histogram of the grouphash
	( $PARAMS->{'verbose'} ) && print "Orginal Grouphash (BIN::# Elements)\n";
	foreach my $bin ( sort { $a <=> $b } keys %$grouphash )
	{
		( $PARAMS->{'verbose'} ) && print $bin. ": ", @{ $grouphash->{$bin} } . "\n";
	}
	return $grouphash;
}

sub generate_grouped_hash_join_single_bins
{
	my $seqhash   = shift;
	my $longest   = shift;
	my $grouphash = {};
	my @groupkeys;
### Genererate the BINS
### The array for the score bins
	my $numbins = int( $longest / $PARAMS->{'binsize'} );
	my $rest    = $longest % $PARAMS->{'binsize'};
	if ($rest) { $numbins++ }
	( $PARAMS->{'verbose'} ) && print "Genrating grouped hash with $numbins bins, binsize $PARAMS->{'binsize'} AA\n";
### Generate an array with the bins
	my @groupsarray;
	my $indexedprots = 0;

	foreach my $i ( 0 .. $numbins )
	{
		push @groupsarray, $i * $PARAMS->{'binsize'};
		push @groupkeys,   $i;
	}
### Iter through the bins and assign the protein ids to the grouphash
	for ( my $i = 0 ; $i < scalar(@groupsarray) ; $i++ )
	{
		foreach my $id ( keys %$seqhash )
		{
			my $length = $seqhash->{$id}->{'length'};
			if ( ( $length >= $groupsarray[$i] ) && ( $length < $groupsarray[ $i + 1 ] ) )
			{
				push @{ $grouphash->{ $groupsarray[$i] } }, $id;
				$indexedprots++;
			}
		}
	}
	( $PARAMS->{'verbose'} ) && print "Indexed $indexedprots Proteins into ProteingroupHash\n";
	### Print Histogram of the grouphash
	( $PARAMS->{'verbose'} ) && print "Orginal Grouphash (BIN::# Elements)\n";
	foreach my $bin ( sort { $a <=> $b } keys %$grouphash )
	{
		( $PARAMS->{'verbose'} ) && print $bin. ": ", @{ $grouphash->{$bin} } . "\n";
	}
	### Sort out those bins that have only a single protein in there and sort them to the nearest neighbour
	( $PARAMS->{'verbose'} ) && print "Sorting grouphash bins with only 1 protein into nearest neighbour bin\n";
	my @keys         = sort { $a <=> $b } keys %$grouphash;
	my @keysbackward = sort { $b <=> $a } keys %$grouphash;
	my @binswithone;
	## Sort from largest to smallest
	for ( my $i = 0 ; $i < scalar(@keysbackward) ; $i++ )
	{
		my $bin            = $keysbackward[$i];
		my $numelementsbin = @{ $grouphash->{$bin} };
		my $nearestback;
		my $nearestfwd;
		if ( $numelementsbin == 1 )
		{

			#print "BIN: $bin";
			#print " is 1\n";
			### Search for the nearest neighbour, calc delta
			## go backwards until 0
			my $deltaback;
			my $deltaforward;
			### Search the nearest fwd if the bins would be sorted asc.
			### (attention we search from backto front)
			## Search FWD (in acs order this is the + direction)
			my $startindex = @keys - $i;
			for my $binnum ( $startindex .. $#keys )
			{
				my $binindex    = $keys[$binnum];
				my $numelements = @{ $grouphash->{$binindex} };

				#my $binAA = $grouphash->{$binindex};
				if ($nearestfwd)
				{
					next;
				} else
				{

					#print "Checking bins fwd: i is $binnum bin $keys[$binnum]\n";
					if ($numelements)
					{
						$nearestfwd   = $binindex;
						$deltaforward = abs( $bin - $binindex );
					}
				}
			}

			#print "Nearest fwd bin is $nearestfwd delta is $deltaforward\n";
			## Search Backwards (in acs order this is the - direction)
			$startindex = $i + 1;
			for my $binnum ( $startindex .. $#keysbackward )
			{
				my $binindex    = $keysbackward[$binnum];
				my $numelements = @{ $grouphash->{$binindex} };
				if ($nearestback)
				{
					next;
				} else
				{

					#print "Checking bins backward: i is $binnum bin $keysbackward[$binnum]\n";
					if ($numelements)
					{
						$nearestback = $binindex;
						$deltaback   = abs( $bin - $binindex );
					}
				}
			}

			#print "Nearest backward bin is $nearestback delta is $deltaback\n";
			### REMOVE FROM THE OLD BIN AND PUT INTO BIN WITH SMALLEST DELTA
			if ( !( $deltaforward || $deltaback ) )
			{
				print "ERROR: nor fwd or backwd bin is set!\n";
				exit;
			}
			my $bintopush;
			if ( $deltaforward && $deltaback )
			{
				### if delta is equal then put into smaller bin
				if ( $deltaforward < $deltaback )
				{

					#print "put into forward bin\n";
					$bintopush = $nearestfwd;
				} else
				{

					#print "put into backward bin\n";
					$bintopush = $nearestback;
				}
			} elsif ($deltaforward)
			{

				#print "no backward bin availaible put into forward bin\n";
				$bintopush = $nearestfwd;
			} elsif ($deltaback)
			{
				$bintopush = $nearestback;

				#print "no fwd bin availaible put into backward bin\n";
			}
			push @{ $grouphash->{$bintopush} }, @{ $grouphash->{$bin} };
			my @emtpy;
			$grouphash->{$bin} = \@emtpy;
		}
	}
	### Print Histogram of the grouphash
	( $PARAMS->{'verbose'} ) && print "Sorted Grouphash (BIN::# Elements)\n";
	foreach my $bin ( sort { $a <=> $b } keys %$grouphash )
	{
		( $PARAMS->{'verbose'} ) && print $bin. ": ", @{ $grouphash->{$bin} } . "\n";
	}
	return $grouphash;
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

   sub fisher_yates_shuffle {
        my $deck = shift;  # $deck is a reference to an array
        my $i = @$deck;
        while ($i--) {
            my $j = int rand ($i+1);
            @$deck[$i,$j] = @$deck[$j,$i];
        }
    }



sub usage()
{
print "
	SOFTWARE: ", basename($0), " version $scriptinfo->{'version'}
	
	AUTHOR: $scriptinfo->{'author'}

	INFORMATION: A software/script to generate decoy sequences.
 	
 	USAGE: ", basename($0), " -Option [Parameter]

	Per default a reverse database from the original input sequences (target sequences are not added) is created
	
	REQUIRED OPTIONS:
	
	-db [] input fasta database
	
	OTHER OPTIONS [default]:
	
	-out [database_decoy.fasta] output database name
	-decoy [reverse|cterm|nterm|ctermpep|ntermpep] 
	
	[reverse]: Total protein sequence is simply reversed;
		ctermpep: c-terminus of each peptide is fixed and the remaining sequence reversed LHIK -> IHLK
		ntermpep: analogous to ctermpep for c-terminus DHIL -> DLIH
		cterm: c-terminus of each peptide is fixed and the remaining sequence randomly scrambled, e.g LHIK->HLIK 
		nterm: analogous to cterm, e.g. DHIL -> DHLI                       				
	-enz [1]
		1: Trypsin
		2: Chymotrypsin
		3: TrypsR
		4: LysC
	
	-v print control information
	-sv print more control information
	-targetseq write also target proteins into the database
	-decoystring [decoy|reverse] the string that is used to label the decoy proteins 
	Please note: in the standard \"reverse\" mode the string \"reverse\" is used otherwise \"decoy\".

	-h print this help

	EXAMPLES:

	Reverse sequences:
	", basename($0)," -db flybase.fasta -out flybase_decoy.fasta 
	Reverse and shuffle sequences:
	", basename($0)," -db flybase.fasta -out flybase_decoy.fasta 
	", basename($0)," -db flybase_decoy.fasta -out flybase_rev_shuffle.fasta -decoy cterm
";
	exit;
}

sub getpeps
{
	my $sequence  = shift;
	my $enzymenum = shift;
	my $verbose   = shift;
	my ( $i, $j, @sequences, $mincuts, @cuts );
	@cuts = ();
	if ( $enzymenum == 1 )
	{
		push @cuts, _digestTryps($sequence);
	} elsif ( $enzymenum == 2 )
	{
		push @cuts, _digestChymoTryps($sequence);
	} elsif ( $enzymenum == 9 )
	{
		push @cuts, _digestTrypsR($sequence);
	} elsif ( $enzymenum == 14 )
	{
		push @cuts, _digestLysC($sequence);
	} elsif ( $enzymenum == 10 )
	{
		push @cuts, _digestAspN($sequence);
	} elsif ( $enzymenum == 15 )
	{
		push @cuts, _digestTrypsAspN($sequence);
	} elsif ( $enzymenum == 16 )
	{
		push @cuts, _digestChymoTryps_lowspec($sequence);
	} elsif ( $enzymenum == 17 )
	{
		push @cuts, _digestTryps_lowspec($sequence);
	} else
	{
		die "no digest method defined for enzyme ", "$!";
	}
	push @cuts, length($sequence) + 1;
	push @cuts, 1;
	my %seen = ();
	@cuts = grep { !$seen{$_}++ } @cuts;
	@cuts = sort { $a <=> $b } @cuts;
	$verbose && print "cuts for $sequence:\n@cuts\n";
	return \@cuts;
}

sub _scrambleprot_keep_cterm
{
	my $seqobj = shift;
	my $cuts   = shift;
	my @cuts   = @$cuts;
	my ($i);
	my $scrambledprot = "";
	for $i ( 0 .. $#cuts - 1 )
	{
		my $seq = $seqobj->subseq( $cuts[$i], $cuts[ $i + 1 ] - 1 );

		# 	print $cuts[$i]," ",$cuts[ $i + 1 ]," ",$seq,"\n";
		my @subseq = split //, $seq;
		my $newseq;
		my $ntrials = 0;
		do
		{
			$newseq = "";
			$ntrials++;
			my @index = random_permuted_index( length($seq) - 1 );

			# 	print "@index\n";
			foreach my $index (@index)
			{
				$newseq .= $subseq[$index];
			}
			$newseq .= $subseq[-1];
			$verbose && print "trial$ntrials; random vector: @index; $seq -> $newseq\n";
		} while ( ( $newseq eq $seq ) && $ntrials <= 3 && length($seq) > 2 );
		$scrambledprot .= $newseq;
	}

	#print $scrambledprot,"\n";
	return $scrambledprot;
}

sub _reversepep_keep_cterm
{
	my $seqobj = shift;
	my $cuts   = shift;
	my @cuts   = @$cuts;
	my ($i);
	my $scrambledprot = "";
	for $i ( 0 .. $#cuts - 1 )
	{
		my $fullpep = $seqobj->subseq( $cuts[$i], $cuts[ $i + 1 ] - 1 );
		my $nontermpep = "";
		if ( ( $cuts[ $i + 1 ] - 2 ) < $cuts[$i] )
		{
			$nontermpep = "";
		}

		#		elsif($cuts[ $i + 1 ] - 2) == $cuts[$i]){
		#			$nontermpep=
		#		}
		else
		{
			$nontermpep = $seqobj->subseq( $cuts[$i], $cuts[ $i + 1 ] - 2 );
		}

		# 	print $cuts[$i]," ",$cuts[ $i + 1 ]," ",$seq,"\n";
		my @subseq  = split //, $fullpep;
		my $ctermAA = $subseq[-1];
		#print $ctermAA."\n";
		my $newseq  = reverse($nontermpep) . $ctermAA;
		$verbose && print "$fullpep -> $newseq\n";
		$scrambledprot .= $newseq;
	}

	#print $scrambledprot,"\n";
	return $scrambledprot;
}

sub _scrambleprot_keep_nterm
{
	my $seqobj = shift;
	my $cuts   = shift;
	my @cuts   = @$cuts;
	my ($i);
	my $scrambledprot = "";
	for $i ( 0 .. $#cuts - 1 )
	{
		my $seq = $seqobj->subseq( $cuts[$i], $cuts[ $i + 1 ] - 1 );

		#print $cuts[$i]," ",$cuts[ $i + 1 ]," ",$seq,"\n";
		my @subseq = split //, $seq;
		my @index = random_permuted_index( length($seq) - 1 );

		# 	print "@index\n";
		my $newseq = "";
		$newseq .= $subseq[0];
		foreach my $index (@index)
		{
			$newseq .= $subseq[ $index + 1 ];
		}
		$verbose && print "@index $seq -> $newseq\n";
		$scrambledprot .= $newseq;
	}
	return $scrambledprot;
}

sub _reverse_sequence
{
	my $sequence = shift;
	return reverse($sequence);
}

sub _digestTryps
{
	my $seq      = shift;
	my @seq      = split //, $seq;
	my @cutsites = ();
	while ( $seq =~ /(K|R)/gi )
	{
		unless ( $seq[ pos($seq) ] =~ /P/i )
		{
			push @cutsites, pos($seq) + 1;
		}
	}
	return @cutsites;
}

sub _digestTryps_lowspec
{
	my $seq      = shift;
	my @seq      = split //, $seq;
	my @cutsites = ();
	while ( $seq =~ /(K|R)/gi )
	{
		push @cutsites, pos($seq) + 1;
	}
	return @cutsites;
}

sub _digestChymoTryps
{
	my $seq = shift;
	my @seq = split //, $seq;

	#	print "$seq\n";
	my @cutsites = ();
	while ( $seq =~ /(F|Y|W)/gi )
	{
		unless ( $seq[ pos($seq) ] =~ /P/i )
		{
			push @cutsites, pos($seq) + 1;
		}
	}
	return @cutsites;
}

sub _digestChymoTryps_lowspec
{
	my $seq = shift;
	my $ENZ = shift;
	my @seq = split //, $seq;

	#	print "$seq\n";
	my @cutsites = ();
	my $cutAA    = $ENZ->{'cutAA'};
	my $notcutAA = $ENZ->{'notcutAA'};

	#	while ( $seq =~ /(K[^P]|R[^P])/gi ) {
	while ( $seq =~ /(F|Y|W|M|L)/gi )
	{
		unless ( $seq[ pos($seq) ] =~ /P/i )
		{
			push @cutsites, pos($seq) + 1;
		}

		#push @cutsites, pos($seq);
	}

	#print "@cutsites\n";
	return @cutsites;
}

sub _digestTrypsAspN
{
	my $seq      = shift;
	my $ENZ      = shift;
	my @cutsites = ();
	my $cutAA    = $ENZ->{'cutAA'};
	my $notcutAA = $ENZ->{'notcutAA'};
	while ( $seq =~ /(D|K[^P]|R[^P])/gi )
	{
		push @cutsites, pos($seq);
	}
	return @cutsites;
}

sub _digestTrypsR
{
	my $seq = shift;
	my $ENZ = shift;
	my @seq = split //, $seq;

	#	print "$seq\n";
	my @cutsites = ();
	my $cutAA    = $ENZ->{'cutAA'};
	my $notcutAA = $ENZ->{'notcutAA'};

	#	while ( $seq =~ /(K[^P]|R[^P])/gi ) {
	while ( $seq =~ /(R)/gi )
	{
		unless ( $seq[ pos($seq) ] =~ /P/i )
		{
			push @cutsites, pos($seq) + 1;
		}

		#push @cutsites, pos($seq);
	}

	#print "@cutsites\n";
	return @cutsites;
}

sub _digestAspN
{
	my $seq      = shift;
	my $ENZ      = shift;
	my @cutsites = ();
	while ( $seq =~ /(D)/gi )
	{
		push @cutsites, pos($seq);
	}
	return @cutsites;
}
###
#
###
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
