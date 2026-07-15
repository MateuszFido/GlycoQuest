#!C:/Perl/bin/perl.exe
use strict;
#---------------------------------------------------------------------------
# xannotate.cgi
# A software/script to to show protein sequence.
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

use CGI qw/:standard :html3/;
use CGI::Carp 'fatalsToBrowser';
use CGI::FastTemplate;
use Socket;
use File::Basename;
use File::Spec;
use File::Path;
use Bio::Perl;
use Bio::Index::Fasta;
# include the directory for the xquest modules ###########
use File::Spec::Functions qw(rel2abs);
use File::Basename;
use FindBin;
use lib "$FindBin::Bin/../modules";
##########################################################
use Environment;
my $env=Environment->new();
my $webconfig=$env->get_path('web.config');

my $form              = new CGI;
my $myself            = self_url();
my $time              = localtime();

my $sessionid         = $form->param('id');
my $hitid             = $form->param('hitid');
my $protein1        = $form->param('protein1');
my $protein2        = $form->param('protein2');
my $topology           =$form->param('topology');
my $xmlfilename       = $form->param('resultxml');
my $seq1=$form->param('seq1');
my $seq2=$form->param('seq2');
my $urlbase           = $form->param('urlbase');
my $spectrumfilename  = $form->param('spectrum');
my $database          = $form->param('database');
my $xquestdeffile     = $form->param('xquestdef');
my $WEBPARAMS         = readtables($webconfig);
my $resultdirbase     = $WEBPARAMS->{'resultdirbase'};
my $resulturlbase     = $WEBPARAMS->{'resulturlbase'};
my @xlinkpositions=split /,/,$form->param('xlinkposition');


my $JSCRIPT = <<EOF;
<!--
function bookmark()
{
var bookmarkurl="$myself"
var bookmarktitle="xQuest results $time";
if(document.all)
{
window.external.AddFavorite(bookmarkurl,bookmarktitle);
}
}
// -->
EOF


print_header();

print '<h1>xQuest results</h1><h3>', $hitid,      '</h3>';
print '<h3>',(join "-",$protein1,$protein2), '</h3><hr>';

print $form->start_form;
print '<p><input type="button" value="Close" onclick="self.close()">';

setparams();
print submit( -name => 'reload' );


print checkbox( -name => 'show sequence',-label =>'show protein sequences');

my $idx;
#print $database;
unless (-e $database){
	
  print '<div class="error">database ',$database,' does not exist</div><p>';
  print "run xquest with the -db option <i>or</i> copy ", basename($database)," to resultsfolder $resultdirbase\/$sessionid<br>";
  print "for questions write to ",$WEBPARAMS->{'email'};

print $form->end_form;
print $form->end_html;
exit 1;
}

unless(-e $database.".idx"){
 #    print "creating index <br>";
     $idx = Bio::Index::Fasta->new(-filename   => $database.".idx",
            -write_flag => 1);
          $idx->make_index($database);
}
else{
  #   print "index exists <br>";
     $idx = Bio::Index::Fasta->new($database.".idx");
}

my @proteins1=split /,/,$protein1;
my @proteins2=split /,/,$protein2;


if(@proteins2){
print "<h2>alpha chain</h2>"  ;
}
print '<table border="1">';
print '<tr class = "tableheading">';
print "<td>";
print "proteinID";
print "</td>";
print "<td>";
print "protein annotation";
print "</td>";
print "<td>";
print "xlink-position";
print "</td>";
print "</tr>";

foreach my $protein (@proteins1){
print "<tr>";
#print $protein."\n";
my $seqobj=$idx->fetch($protein);

print "<td>";

print $seqobj->id;
print "</td>";
print "<td>";
print $seqobj->desc;
print "</td>";
absposition($seqobj,$seq1,$xlinkpositions[0]);

print "</tr>";
}
print '</table>';

if(@proteins2){
print "<h2>beta chain</h2>"  ;

print '<table border="1">';
print '<tr class = "tableheading">';
print "<td>";
print "proteinID";
print "</td>";
print "<td>";
print "protein annotation";
print "</td>";
print "<td>";
print "xlink-position";
print "</td>";
print "</tr>";

foreach my $protein (@proteins2){
print "<tr>";
my $seqobj=$idx->fetch($protein);


print "<td>";

print $seqobj->id;
print "</td>";
print "<td>";
print $seqobj->desc;
print "</td>";
absposition($seqobj,$seq2,$xlinkpositions[1]);

print "</tr>";
}
print '</table>';


}



print '<p><input type="button" value="Close" onclick="self.close()">';
print $form->end_form;

print $form->end_html;



sub print_header {
        print $form->header('text/html');
        print $form->start_html(
                -title   => 'xQuest Annotation',
                -author  => 'rinner@xquest.org',
                  -base    => '/cgi-bin/xquest/',
                -expires => '+30s',
                -script  => $JSCRIPT,
                -meta    => {
                        'keywords'  => 'ms cross-link ',
                        'copyright' => 'copyright 2006 Oliver Rinner'
                },
                -style => { 'src' => $WEBPARAMS->{'css_stylesheet'}, },
        );

}


sub absposition{
my $seqobj=shift;
my $xlinkedseq=shift;
my $xlinksite=shift;
my @xlinksites=();
my $sequence=$seqobj->seq;
my $sitestring="";
chomp($sequence);
        while ( $sequence =~ /$xlinkedseq/gi ) {
                push @xlinksites, pos($sequence);
        }
   foreach my $site(@xlinksites){
   #  $sitestring.=$site-length($xlinkedseq)+$xlinksite;
print "<td>";
print $site-length($xlinkedseq)+$xlinksite;
print "</td>";
     if(param('show sequence')){
     print "<tr>";
print "<td>";
print "</td>";
print "<td>";
    #print $site;
     printformatseq($sequence,[$site-length($xlinkedseq),$site-1],$site-length($xlinkedseq)+$xlinksite);
    print "</td>";

      print "</tr>";
    }

   }
  return   $sitestring;
}


sub printformatseq{
my $seq=shift;
my $highlight=shift;
my $xlinkposition=shift;
my @seq=split //,$seq;
my $i;
print '<div class="sequence">';
for $i (0..$#seq+1){
 if ($i%50 == 0 && $i>0){
  print "<br>";
  }
  if($i>=$highlight->[0] && $i<=$highlight->[1]){
  if ($i==$xlinkposition-1){
     print '<b>',$seq[$i],'</b>';
  }else{
   print '<font color="#FF8000">',$seq[$i],'</font>';
                               }
   }
   else{
   print $seq[$i];
   }
}
print "</div>"
}



sub readtables {
        my $webparam = shift;
        my %WEBPARAMS;
        open WEBCONFIG, "<$webparam"
          or warn "could not open web config file $! ignoring";
        while (<WEBCONFIG>) {
                chomp;
                my @keyvalue = split /::/;
                $WEBPARAMS{ $keyvalue[0] } = $keyvalue[1];
        }
        return \%WEBPARAMS;
}



sub setparams {

        print $form->hidden(
                -name  => 'id',
                -value => param('id'),
        );
        print $form->hidden(
                -name  => 'protein1',
                -value => param('protein1'),
        );


              print $form->hidden(
                -name  => 'protein2',
                -value => param('protein2'),
        );
        print $form->hidden(
                -name  => 'hitid',
                -value => param('hitid'),
        );

        print $form->hidden(
                -name  => 'topology',
                -value => param('topology'),
        );
        print $form->hidden(
                -name  => 'seq1',
                -value => param('seq1'),
        );
        print $form->hidden(
                -name  => 'seq2',
                -value => param('seq2'),
        );

        print $form->hidden(
                -name  => 'xlinkposition',
                -value => param('xlinkposition'),
        );

        print $form->hidden(
                -name  => 'database',
                -value => param('database'),
        );





}
