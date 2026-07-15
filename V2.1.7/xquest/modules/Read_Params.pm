package Read_Params;
use strict;
#---------------------------------------------------------------------------
# Module: Read_Params.pm
# Author(s): Thomas Walzthoeni, xQuest specific modules are based on 
# original work by Oliver Rinner.
# Description: Module for parameter handling.
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
use warnings;
use File::Basename;
use Data::Dumper;
use Cwd;

sub readtables
{
	my $table     = shift;
	my $xquestdef = shift;
	my $webconfig = shift;
	my $masslist  = shift;
	
	my $xquestdir=shift;
	my $verbose =shift;
	my $sessionid=shift; 		## usually empty only for web scripts used
	my $resultdirbase=shift;
	#print "$xquestdef\n";
	
	
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
			if ($line ne "")
			{
				#print $line."\n";
				my @results = split(' ', $line );
				$PARAMS{ $results[0] } = $results[1];
				#print "key:".$results[0]."value ".$results[1]."\n";
			}

			#print $_[0], "\t", $PARAMS{ $_[0] }, "\n";
		} elsif ( defined($enzymedef) )
		{
			if ($line ne "")
			{
				my @results = split( ' ', $line );
				$ENZ{ $results[0] }->{'name'}     = $results[1];
				$ENZ{ $results[0] }->{'cutAA'}    = $results[3];
				$ENZ{ $results[0] }->{'notcutAA'} = $results[4];
			}
		} elsif ( defined($modifications) )
		{
			if ($line ne "")
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
			$MSTAB{$key}->{'average'} += $MOD{$key};
			#print $MSTAB{$_}->{'native'},"\n";
		}
	}

	#define variable modification X
	if ( $PARAMS{'variable_mod'} )
	{
	#	my ( $AA, $delta ) = split /,|:/, $PARAMS{'variable_mod'};
	#	$MSTAB{'X'}->{'native'} = $MSTAB{$AA}->{'native'} + $delta;
	#	$MSTAB{'X'}->{'average'} = $MSTAB{$AA}->{'average'} + $delta;
		my @AAlist = ('X', 'U', 'B', 'J');
		my @AAdelta = split /,|:/, $PARAMS{'variable_mod'};
		my $nmods = int( ( scalar @AAdelta ) / 2 );
		for my $i ( 0 .. $nmods - 1 ) {
			my $pseudo = $AAlist[$i];
			$MSTAB{$pseudo}->{'native'} = $MSTAB{$AAdelta[2*$i]}->{'native'} + $AAdelta[2*$i + 1];
			$MSTAB{$pseudo}->{'average'} = $MSTAB{$AAdelta[2*$i]}->{'average'} + $AAdelta[2*$i + 1];
			$verbose && print "modification $pseudo ($AAdelta[2*$i]): ", $MSTAB{$pseudo}->{'native'}, "\n";
		}
		# Maxes out at 4 variable modifications specifiable, could also use more but running out of space (O for C-terminus? Z for N-terminus coded)
	}
	#if ($MSTAB{'X'}->{'native'} ){
	#$verbose && print "modificaton X: ", $MSTAB{'X'}->{'native'}, "\n";
	#}
	#if ( $PARAMS{'AArequired'} )
	{
		$PARAMS{'AArequired2'} = $PARAMS{'AArequired'};
		chomp( $PARAMS{'AArequired'} );
		my @xlinktargets = split /,|:|\|/, $PARAMS{'AArequired'};
		my $AAstring;
		foreach my $AA (@xlinktargets)
		{
			$AAstring = join "|", @xlinktargets;
		}
		$PARAMS{'AArequired'} = $AAstring;
	}
	
	if ( $PARAMS{'possibleTopology'} )
	{
		chomp( $PARAMS{'possibleTopology'} );
		my @possibletopologies = split /,/, $PARAMS{'possibleTopology'};
		my $topos={};
		foreach my $topo (@possibletopologies)
		{
			$topos->{$topo} = 1;
			# also index the reverse
			$topos->{reverse($topo)} = 1;
		}
		$PARAMS{'possibleTopology'} = $topos;
	}	
	
	
	
	
	unless ( $PARAMS{'usenprescores'} )
	{
		$PARAMS{'usenprescores'}=100;
	}
	$verbose && ( print "xlink targets: ", $PARAMS{'AArequired'}, "\n" );
	
	my $dbname;
	my $dbbasename=basename($PARAMS{'database'});
	#print "DB:".$PARAMS{'database'}."\n";

	if ( -e $PARAMS{'database'} )
	{
	$dbname = $PARAMS{'database'};
	} elsif ( -e File::Spec->catfile( $xquestdir, $PARAMS{'database'} ) )
	{
	$dbname = File::Spec->catfile( $xquestdir, $PARAMS{'database'} );
	}elsif ( -e File::Spec->catfile( $resultdirbase, $sessionid, $dbbasename ) ){
	$dbname=File::Spec->catfile( $resultdirbase, $sessionid, $dbbasename );
	}elsif(-e File::Spec->catfile($dbbasename)){
	### Then the db file is in the cwd
	my $dir = getcwd;
	$dbname=File::Spec->catfile( $dir, $dbbasename	);
	}
	 else
	{
		die "cannot open database file $PARAMS{'database'} $!";
	}
	## Get the Db path/basename is used for the db indices
	$dbname =~ s/\.\w+//;

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
	unless ( $PARAMS{'writetodiskaftern'} )
	{
		$PARAMS{'writetodiskaftern'} = 100;
	}
	unless ( $PARAMS{'tolerancemeasure_ms2'} )
	{
		$PARAMS{'tolerancemeasure_ms2'} = "Da";    #$PARAMS->{'tolerancemeasure'};
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

sub print_params
{
	my $hashref = shift;
	foreach my $key ( sort keys %$hashref )
	{
		my $value = $hashref->{$key};
		#print ref($value)."\n";

    	if ( (ref($value) eq "HASH") || (ref($value) eq "ARRAY") )
		{
		print "$key ==> ".ref($value)."-Reference\t\n";
		#print "key: ".$key."\t\n";
		
		if ((ref($value) eq "HASH")){
		foreach my $key (keys %$value){
		print "\t".$key."=>".$value->{$key}."\n";	
		}	
		}
		if ((ref($value) eq "ARRAY")){		
		for (my $i=0;$i<scalar(@$value);$i++){
		print "\t";
		print $i."=>".$value->[$i];
		print "\n";	
		}		
		}
		if ((ref($value) eq "SCALAR")){		
		print "$key \t  =>  $$value\n";	
		}
		}else{
		unless ($value) { $value = "not defined" }
		
		if ($key=~ m/\#/){
		#print "TRUE";
		#exit;
			next;
		}
		my $length = length($key);
		
		if ( $length <6  )
		{
			print "$key \t\t\t  =>  $value\n";
		} elsif ($length <=14)
		{
			print "$key \t\t  =>  $value\n";
		}else{
		print "$key \t  =>  $value\n";	
		}
				
		}
		

	}
}


1;