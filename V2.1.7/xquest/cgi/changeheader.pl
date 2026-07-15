#!/usr/bin/perl
use strict;
#---------------------------------------------------------------------------
# changeheader.pl
# A software/script to change cgi headers
# Execute changeheader.pl -help to display information and usage options.
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
use Getopt::Long;
use File::Copy;
my $OS = $^O; # get the current OS
my $newdb;
my $olddb;

if($OS eq "MSWin32") {
	$newdb ='#!C:/Perl/bin/perl.exe';
	$olddb='#!C:/WINDOWS/Perl/bin/perl.exe';
	}
if($OS eq "linux") {
	$newdb ='#!/usr/bin/perl';
	$olddb='#!C:/Perl/bin/perl.exe';
	}

my $cgifiles='*.cgi';
my $help;

GetOptions(
        'xml=s' =>\$cgifiles,
        'old=s'       => \$olddb,
        'new=s'        => \$newdb,
        'help'        => \$help,

  );

unless ($newdb || $olddb){
print "newdb or olddb not defined";
&usage();
} 
&usage if $help;


my @xmlfiles=glob($cgifiles);
foreach my $xmlfile ( @xmlfiles){

# get fh to the actual file
open XMLFILE, "$xmlfile" or die $!;
# def filename of the modified file
my $xmlmod=join ".",$xmlfile,"mod";
#get fh the mod file
open XMLOUT, ">$xmlmod" or die $!;
chomp;
my $newname;
($newname=$xmlfile)=~s/.xml/_mod.xml/;
print "$xmlfile -> $newname\n";

my $outfile = join "",$xmlfile,".mod";

## swith the olddb to the newdb and write into mod file
while(<XMLFILE>){
s/$olddb/$newdb/;
print XMLOUT $_;
}
 close(XMLFILE);
 close(XMLOUT);
# 
copy ($xmlmod,$xmlfile);
## delete the mod file
unlink ($xmlmod);
 }


sub usage {
    print "Changes the local Perl path tag in cgi files.\n";
    print "usage: $0\n";
    print "options
        -old (default $olddb)
        -new (default $newdb)\n
        -help print this help.
";
    exit;
}
