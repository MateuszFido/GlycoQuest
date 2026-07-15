package Common;
#---------------------------------------------------------------------------
# Module: Common.pm
# Author(s): Thomas Walzthoeni, xQuest specific modules are based on 
# original work by Oliver Rinner.
# Description: Some frequently used functions.
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

#---------------------------------------------------------------------------
# prettytime retuns a nicely formatted time string
#---------------------------------------------------------------------------
sub prettytime {
	my @months   = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
	my @weekDays = qw(Sun Mon Tue Wed Thu Fri Sat Sun);
	my (
		$second,     $minute,    $hour,
		$dayOfMonth, $month,     $yearOffset,
		$dayOfWeek,  $dayOfYear, $daylightSavings
	  )
	  = localtime();
	my $year    = 1900 + $yearOffset;
	my $theTime = "$hour:$minute:$second, $weekDays[$dayOfWeek] $months[$month] $dayOfMonth, $year";
	return $theTime;
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
	open FILE, $filename or die $!;
	while ( my $line = <FILE> )
	{
		chomp($line);
		unless ($line){
			next;
		}
		if ($verbose) { print "Reading line $line\n" }
		
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

#===  FUNCTION  ================================================================
#  NAME:        save_array_to_file
#  PURPOSE:     save_array_to_file
#  DESCRIPTION: Save an array to a file (every arrayelement is one line, array elements should contain \n)
#  PARAMETERS:  $filename, $arrayreference, $append (1 if string should be appended)
#  RETURNS:     1 if done correctly
#===============================================================================
sub save_array_to_file
{
	my $filename = shift;
	my $arraref  = shift;
	my $lineend  = shift;
	my $append   = shift;
	my $verbose =shift;
	
	if ($append)
	{
		open MYOUTFILE, ">>", "$filename" or die "cannot open file $filename $!";
	} else
	{
		open MYOUTFILE, ">", "$filename" or die "cannot open file $filename $!";
	}
	
	foreach my $line (@$arraref){
	#print "writing $line\n";
	print MYOUTFILE $line;
	if ($lineend){
	print MYOUTFILE "\n";
	}
	}
	close MYOUTFILE;
	if ($verbose){
	print "Wrote ". @$arraref. " lines to file $filename\n";
	}
	
	return 1;
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
	if ($verbose)
	{
		print "Checking file: " . $filename . "...";
	}

	# -e is for exists, -r for readable
	unless ( ( -e $filename ) && ( -r $filename ) )
	{
		print "Error: Cannot find/read the file $filename.\n";
		exit 0;
	} else
	{
		if ($verbose) { print " ok.\n" }
	}
	return;
}


sub print_affi_and_changelog
{
my $scriptinfo = shift;
print "\nVersion " . $scriptinfo->{'version'} . " written by ";
print $scriptinfo->{'author'} . "\n";
print "Affiliation: " . $scriptinfo->{'affi'} . "\n";
print "mailto: ", $scriptinfo->{'mailto'} . "\n";
print "Changelog:\n";
foreach my $key ( sort { $a <=> $b } keys %{ $scriptinfo->{'clog'} } )
{
print "Version " . $key . ": " . $scriptinfo->{'clog'}->{$key};
}
return;
}

#===  FUNCTION  ================================================================
#  NAME:        print_params
#  PURPOSE:     Printing a parameter hash
#  DESCRIPTION: Sorts a parameter Hash alphabetically and prints it
#  PARAMETERS:  $PARAMS (Hashreference)
#  RETURNS:     void
#===============================================================================
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





1;