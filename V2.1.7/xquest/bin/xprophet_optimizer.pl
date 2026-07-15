#!/usr/bin/env perl
use strict;
#---------------------------------------------------------------------------
# xprophet_optimizer.pl
# Optimize a parameter by caling xprophet, and writing to output file, gen optimal xprophet file.
# Execute xprophet.pl -help to display information and usage options.
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
##########################################################
# Include modules dir as lib that is relative to the Script path
##########################################################
use FindBin;
use lib "$FindBin::Bin/../modules";
##########################################################
use Environment;
use Read_Params;
use Data::Dumper;
use Storable;
use Cwd;
use File::Path;
use XML::TreeBuilder;
use Getopt::Long;
use File::Basename;

#---------------------------------------------------------------------------
#  Variables
#---------------------------------------------------------------------------
my ( $input, $outfilename, $deffile, $help, $verbose, $version, $opti, $forceindex, $fdr );

#---------------------------------------------------------------------------
#  Default values
#---------------------------------------------------------------------------
$version     = "1.0";
$input       = "xquest.xml";
$outfilename = "optimizer_results.xls";
$deffile     = "xproph.opt.def";
$fdr=0.05;
my $author      = "Thomas Walzthoeni";
my $affiliation = "ETH Zurich, Institute of Molecular Systems Biology, Wolfgang Pauli-Str. 16, CH-8093 Zurich";
my $mailto      = "walzthoeni\@imsb.biol.ethz.ch";
my $cwd         = getcwd();

#---------------------------------------------------------------------------
#  Options that are passed by arguments
#---------------------------------------------------------------------------
GetOptions(
			'in=s'   => \$input,          ## input filename"
			'out=s'  => \$outfilename,    ## output filename of the optimizer results
			'def=s'  => \$deffile,        ## deffilename
			'help'   => \$help,
			'v'      => \$verbose,
			'opti=s' => \$opti,
			'forceindex'=>\$forceindex,
			'fdrcutoff=s'=>\$fdr,
);

#---------------------------------------------------------------------------
#  Define paths to programs
#---------------------------------------------------------------------------
my $env       = Environment->new;
my $xqbasedir = File::Spec->catfile( $env->get_path('xquest_stable') );
my $xqbin     = File::Spec->catfile( $xqbasedir, "bin" );

# Define paths to programs
my $mergeprogram          = File::Spec->catfile( $xqbin, "mergexml.pl" );
my $annotateprogram       = File::Spec->catfile( $xqbin, "annotatexml.pl" );
my $indexspecfilesprogram = File::Spec->catfile( $xqbin, "indexspecfiles.pl" );
my $xprophet              = File::Spec->catfile( $xqbin, "xprophet.pl" );

#---------------------------------------------------------------------------
#  Check the input
#---------------------------------------------------------------------------
&usage() if $help;
unless ( -e $input )
{
	print "ERROR: Can't read input file $input\n";
	exit;
}

#---------------------------------------------------------------------------
#  Generate the file to store the optimizer results, generate new
#---------------------------------------------------------------------------
my $optiresultsfile = File::Spec->catfile( $cwd, "$outfilename" );
my $headerline= "minionsmatched \t FDR \t intra-protein xls \t inter-protein xls \t mono and looplinks\n";
open RES, ">", $outfilename or die "cannot open file $outfilename $!";
print RES $headerline;
close (RES);

#---------------------------------------------------------------------------
#  1. Check if there is already a template def file in the folder;
#---------------------------------------------------------------------------
if ( -e $deffile )
{
print "Deffile is availaible\n";
} else
{
## Run xproph and generate deffile
	my $cmd = "$xprophet -in $input -def $deffile -getdef -v";
	my $res = system($cmd);
	die "xprophet returned an error; error code: $res\n" if ($res);
}

#---------------------------------------------------------------------------
#  2. Run optimizer
#---------------------------------------------------------------------------
for my $i (0 .. 10){
print "N: $i\n";
my $cmd = "$xprophet -in $input -def $deffile -v -optifile \"$optiresultsfile\" -opti $fdr -minionsmatched $i";
print "CMD: $cmd\n";
my $res = system($cmd);
}


sub usage()
{
	print "
	SOFTWARE: ", basename($0), " version $version
	
	AUTHOR: Thomas Walzthoeni

	INFORMATION: A software/script optimize xprphet parameters. Currently only the minionsmatched parameter.
 	
 	USAGE: ", basename($0), " -Option [Parameter]

	REQUIRED OPTIONS [defaults]: 
	-in [$input] xquest.xml input filename
	-out [$outfilename] output filename.
 	-def [xproph.def] xProphet definition file, see DEFINTION FILE for further information.
    -fdrcutoff[$fdr] FDR cutoff used by xProphet to filter the results.
	OTHER OPTIONS:
 	-help print this information.
 	

	EXAMPLE: ", basename($0), " -in merged_xquest.xml -out $outfilename
 		
	";
	exit;
}
