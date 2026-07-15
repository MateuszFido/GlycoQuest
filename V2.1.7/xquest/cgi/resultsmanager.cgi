#!C:/Perl/bin/perl.exe
use strict;

#---------------------------------------------------------------------------
# resultsmanager.cgi
# A software/script to manage search results.
# Author: Thomas Walzthoeni
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
use strict;
use CGI qw/:standard :html3/;
use CGI::Carp qw(fatalsToBrowser);
use CGI::FastTemplate;
use File::Basename;
use File::Spec;
use File::Path;
use Data::Dumper;
use Cwd;
use File::stat;
use Time::localtime;
use Time::Piece;
use Storable;

# include the directory for the xquest modules ###########
use File::Spec::Functions qw(rel2abs);
use File::Basename;
use FindBin;
use lib "$FindBin::Bin/../modules";
##########################################################
use Environment;
use SimpleAuth;

#---------------------------------------------------------------------------
#  Version info
#---------------------------------------------------------------------------
my $version = "1.1";

#---------------------------------------------------------------------------
#  Create the cgi object and initialize parameters
#---------------------------------------------------------------------------
my $form = new CGI;

#---------------------------------------------------------------------------
#  Read the webparams: defined in the web.config file, path must be defined in Environment.pm
#---------------------------------------------------------------------------
my $env       = Environment->new();
my $webconfig = $env->get_path('web.config');
my $WEBPARAMS = readtables($webconfig);

#---------------------------------------------------------------------------
#  Create the path to the resultdirectory that should be displayed
#---------------------------------------------------------------------------
my $resultdirectory;

# The resutdirectory base e.g. /home/xquest/
my $resultdirbase = $WEBPARAMS->{'resultdirbase'};    # defined in the webparams

# Check if a subfolder should be displayed, useful if more people use it: eg. /username
my $resultdirname;

# $resultdirname=$form->param('resultdir'); ## You can use this if you want to pass the param by URL, just uncomment
# Otherwise use the username, therfore you have to add a .htaccess file to the cgi directory and require authentication
unless ($resultdirname)
{
	$resultdirname = $ENV{'REMOTE_USER'};
}

# Check if a subfolder is defined; this is a specific folder that can be used e.g. /xquestresults
my $resultsubfolder = $WEBPARAMS->{'resultsubfolder'};

# Create relative path for the URL construction
my $relpathtobase;
unless ($resultsubfolder)
{
	$relpathtobase = $resultdirname;
} else
{
	$relpathtobase = $resultdirname . "\/" . $resultsubfolder;
}
## Concatenate the paths, the path is then e.g. /home/xquest/username/xquestresults
$resultdirectory = File::Spec->catfile( $resultdirbase, $resultdirname, $resultsubfolder );

#---------------------------------------------------------------------------
#  Create or load the indexfilename where the structure is stored
#---------------------------------------------------------------------------
my $indexfilename = File::Spec->catfile( $resultdirectory, "index.db" );
my $make_new_index;
my $INDEX = {};
if ( ( -e $indexfilename ) && ( !$make_new_index ) )
{
## reads the Hash from the DB
	$INDEX = retrieve($indexfilename);
} else
{

	# Create new
	store $INDEX, $indexfilename;
}

#---------------------------------------------------------------------------
#  Print the HTML header
#---------------------------------------------------------------------------
print_header();
print "<h1>xQuest/xProphet results manager</h1>";
my $debugparam = 0;
if ($debugparam)
{
	debug_param($form);
}
my $myself  = $form->url();
my $infourl = $myself . "?info=1";
## Print Info if selected
print '<a href="' . $infourl . '">How to use</a>';
print "<hr>";

#---------------------------------------------------------------------------
#  Print info if selected
#---------------------------------------------------------------------------
if ( $form->param('info') )
{
	print "<div>";
	print "<h2>Usage</h2>";
	print "<div class=\"doc\">";
	print "With the xQuest/xProphet manager, search results can be organized.<br>";
	print "Availaible folders in the result rootdirectory (see below) are displayed, newly added folders are added to the root project folder.<br>";
	print "Create or delete project folders using the textfield, project folders are then availaible in the dropdown lists.<br>";
	print "To move search results into the project folders, select them using the checkboxes, and then select the project folder from the dropdown list (Move selected result(s) to:), then use update.<br>";
	print "To delete a project folder, type the name of the folder into the textbox and select delete, then update.<br>";
	print "If a project folder is deleted, search results are moved back to the root folder<br>";
	print "To display a project folder, select the folder from the dropdown list, then use update.<br>";
	print "Note: The manager uses a file (index.db) to store the structure, it never deletes any resultfolders.<br>";
	print "If you want to reset the structure, delete the index.db file in your root resultdirectory.<br>";
	print "Folder can be deleted by moving them into the Trash folder. In the Trash folder selected folders can be deleted.\n";
	print "</div>";
	print "<br><div>";
	print '<a href="' . $myself . '">Back</a>';
	print "</div>";
	print $form->end_form;
	print "<br>";

	#print "For debugging:<br>";
	print '<div class="code">';
	print "Folders used with this settings: <br>";
	print "Resultdirbase: $resultdirbase<br>";
	print "Resultdirname: $resultdirname<br>";
	print "Resultsubfolder: $resultsubfolder<br>";
	print "Result root directory: $resultdirectory<br>";
	print "</div>";
	print_footer();
	print $form->end_html;
	exit;
}

#---------------------------------------------------------------------------
#  For debuging
#---------------------------------------------------------------------------
#debug_param($form);
#---------------------------------------------------------------------------
#  Statusvariable, is printed below the options
#---------------------------------------------------------------------------
my $statmsg;

#---------------------------------------------------------------------------
#  Set the Params
#---------------------------------------------------------------------------
my $newpjxfoldername = $form->param('newpjxfoldername');
unless ($newpjxfoldername) { $newpjxfoldername = "" }
my $dispfolder = $form->param('dispfolder');
unless ($dispfolder) { $dispfolder = 'root' }
my $sortby = $form->param('sortby');
unless ($sortby) { $sortby = 'date' }
my $faction = $form->param('faction');
unless ($faction) { $faction = 'create' }
my $deleteselected = $form->param('del_form');
my $jsvalue        = $form->param('jsvalue');

#---------------------------------------------------------------------------
#  Check if the action is to delete folders from the Trahs
#---------------------------------------------------------------------------
if ($deleteselected)
{
	$statmsg .= delete_folders( $form, $INDEX, 0 );
	#print "True\n";
}

#---------------------------------------------------------------------------
#  Read the directories into a hash, also generate the links
#---------------------------------------------------------------------------
my $directoryhash = read_dirs();


#---------------------------------------------------------------------------
#  Update the index, where the folderstructure is stored, check for new folders, remove deleted folders
#---------------------------------------------------------------------------
$statmsg .= update_index( $directoryhash, $INDEX );


#---------------------------------------------------------------------------
#  Add a new projectfolder if selected or delete
#---------------------------------------------------------------------------
if ($newpjxfoldername)
{
## Add to the index
	$statmsg .= add_project_folder( $newpjxfoldername, $faction, $INDEX );
## unset the params
	$form->param( 'newpjxfoldername', "" );
	$dispfolder = "root";
	$form->param( 'faction', "create" );
}

#---------------------------------------------------------------------------
#  Move a folder if selected
#---------------------------------------------------------------------------
unless($deleteselected){
$statmsg .= movefolders( $form, $INDEX );
}
#---------------------------------------------------------------------------
#  Print the options form
#---------------------------------------------------------------------------
print_options_form();
print "<div class=\"highlight\">$statmsg<\/div>";

#---------------------------------------------------------------------------
#  Get the folders that should be displayed
#---------------------------------------------------------------------------
my $folderstbd = get_folders_tobedisplayed( $dispfolder, $sortby, $directoryhash, $INDEX );
print "<h2>Results for project folder: $dispfolder</h2>";

#---------------------------------------------------------------------------
#  Display the selected folder
#---------------------------------------------------------------------------
display_folders( $folderstbd, $directoryhash );

#---------------------------------------------------------------------------
#  END
#---------------------------------------------------------------------------
print $form->end_form;
print_footer();
print $form->end_html;

#---------------------------------------------------------------------------
#  Functions
#---------------------------------------------------------------------------
sub delete_folders
{
	my $form  = shift;
	my $INDEX = shift;
	my $test=shift;
	my $error = 0;
## First check if the selected folders are in the Trash folder
	my @parameters = $form->param();
	my $statusmsg;
	my @folderstobedeleted;
	for ( my $i = 0 ; $i < @parameters ; $i++ )
	{
		my $name = $parameters[$i];
		if ( $name =~ m/selected/ )
		{
			my @validationarray = split( /:_:/, $name );
			my $foldername = $validationarray[1];

			# my $targetfolder    = $form->param('targetfolder');
			push @folderstobedeleted, $foldername
			  ## Change Index
		}
	}
## Check if these folders are in the Trash folder
	foreach my $folder (@folderstobedeleted)
	{
		unless ( $INDEX->{'resdirs'}->{$folder} eq "Trash" )
		{
			$error     = 1;
			$statusmsg = "Folder is not in the Trash folder<br>";
		}
	}
## Check if jsvalue is 1
	if ($jsvalue==0)
	{
		$error     = 1;
		$statusmsg = "Action canceled<br>";
	}
## Delete function
if ($error == 0){	
	foreach my $folder (@folderstobedeleted)
	{
		my $foldertodelete = File::Spec->catfile( $resultdirbase, $relpathtobase, $folder );
## Remove the folder
		#print "Folder to delete: $foldertodelete\n";
		unless ( -e $foldertodelete )
		{
			$statusmsg .= "No folder to delete!<br>";
		} else
		{
			unless ($test){
			rmtree( [$foldertodelete] );
			}
			$statusmsg .= "Folder $foldertodelete deleted!<br>";
		}
	}
}
	return $statusmsg;
}


sub print_footer
{
	my $htmltable;
	$htmltable .= "<div style=\"position:relative;border:0px solid #00ff00;bottom:-25px; text-align: right\" >";
	$htmltable .= "xQuest/xProphet results manager, version $version, Thomas Walzthoeni / ETH Zurich";
	$htmltable .= '</div>';
	print $htmltable;
}

sub movefolders
{
	my $form       = shift;
	my $INDEX      = shift;
	my @parameters = $form->param();
	my $statusmsg;
	for ( my $i = 0 ; $i < @parameters ; $i++ )
	{
		my $name = $parameters[$i];
		if ( $name =~ m/selected/ )
		{
			my @validationarray = split( /:_:/, $name );
			my $foldername      = $validationarray[1];
			my $targetfolder    = $form->param('targetfolder');
			$statusmsg .= "Moved results folder: <b>$foldername</b> to <b>$targetfolder</b><br>";
			$INDEX->{'resdirs'}->{$foldername} = $targetfolder;
			## Change Index
		}
	}

	# store index
	store $INDEX, $indexfilename;
	return $statusmsg;
}

sub display_folders
{
	my $folderstbd    = shift;
	my $directoryhash = shift;
	my $htmltable;
	$htmltable .= '<div class="tabletext" width="100%">' . "\n";
	$htmltable .= '<table border="1"  width="100%">';
	$htmltable .= '<tr class="tableheading">';
	$htmltable .= '<tr class="tableheading">' . "\n";
	$htmltable .= '<td>Select</td>';
	$htmltable .= '<td>Name</td>';
	$htmltable .= '<td>Date (last modified)</td>';
	$htmltable .= '<td>Results</td>';
	$htmltable .= '<td>xQuest Deffile</td>';
	$htmltable .= '<td>xProphet Deffile</td>';
	$htmltable .= '</tr>' . "\n";
	foreach my $folder (@$folderstbd)
	{
## Generate a checkbox for every folder
		my $movecheckname = "selected:_:$folder";
		$htmltable .= '<tr>' . "\n";
		$htmltable .= '<td>' . $form->checkbox( -name => $movecheckname, -value => 1, -checked => 0, -label => '' ) . '</td>';
## Get url
		my $url = $directoryhash->{'dirs'}->{$folder};
		$htmltable .= '<td>' . $folder . '</td>';
## Get date
		my $datelmod = $directoryhash->{'time'}->{$folder};
		$htmltable .= '<td>' . $datelmod . '</td>';
## URL
		my $link = "<a target=\"_blank\" href=$url>view</a><br>";
		$htmltable .= '<td>' . $link . '</td>';
		my $xqdeffileurl = $directoryhash->{'xqdeffileurl'}->{$folder};
		my $xqdeffillink = "<a target=\"_blank\" href=$xqdeffileurl>view</a><br>";
		$htmltable .= '<td>' . $xqdeffillink . '</td>';
		my $xpdeffileurl = $directoryhash->{'xpdeffileurl'}->{$folder};
		my $xpdeffillink = "<a target=\"_blank\" href=$xpdeffileurl>view</a><br>";
		$htmltable .= '<td>' . $xpdeffillink . '</td>';
		$htmltable .= '</tr>' . "\n";
	}
	$htmltable .= '</table>' . "\n";
	print $htmltable;
}

sub get_folders_tobedisplayed
{
	my $dispfolder    = shift;    ## The folder to be displayed
	my $sortby        = shift;
	my $directoryhash = shift;
	my $INDEX         = shift;
	my $folderstodisp = {};

	#print "SEL F: $dispfolder<br>";
## get all directories
	my @allresdirs = keys %{ $INDEX->{'resdirs'} };
	foreach my $resdir (@allresdirs)
	{
## check if in the displayed folder
		my $projectfolder = $INDEX->{'resdirs'}->{$resdir};

		#
		if ( $dispfolder eq $projectfolder )
		{

			#print "$folder";
			$folderstodisp->{$resdir} = 1;
		}
	}
## ITERATE over all folders
	my $foldersbysort;
	if ( $sortby eq "name" )
	{
		$foldersbysort = $directoryhash->{'resdirectories'};
	}
	if ( $sortby eq "date" )
	{
		$foldersbysort = $directoryhash->{'resdirectoriesdesctime'};
	}
	my @folders;
	foreach my $folder (@$foldersbysort)
	{

		#print "$folder<br>";
		if ( $folderstodisp->{$folder} )
		{
			push @folders, $folder;
		}
	}
	return \@folders;
}

#---------------------------------------------------------------------------
#  Print the directories of the selected project folder
#---------------------------------------------------------------------------
sub add_project_folder
{
	my $newpjxfoldername = shift;
	my $faction          = shift;
	my $INDEX            = shift;
	my $statmsg;
	if ( $newpjxfoldername eq "root" )
	{
		$statmsg .= "Folder <b>root</b> cannot be added/deleted.<br>";
		return $statmsg;
	}
	if ( $newpjxfoldername eq "Trash" )
	{
		$statmsg .= "Folder <b>root</b> cannot be added/deleted.<br>";
		return $statmsg;
	}
	if ( $faction eq "create" )
	{
		unless ( $INDEX->{'projectfolders'}->{$newpjxfoldername} )
		{
			$INDEX->{'projectfolders'}->{$newpjxfoldername} = 1;
			$statmsg .= "Added project folder <b>$newpjxfoldername<\/b><br>";
		} else
		{
			$statmsg = "Projectfolder <b>$newpjxfoldername</b> already exists.<br>";
		}
	}
	if ( $faction eq "delete" )
	{
		if ( $INDEX->{'projectfolders'}->{$newpjxfoldername} )
		{
			## move results to the root folder
			my @folders = keys %{ $INDEX->{'resdirs'} };
			foreach my $f (@folders)
			{
				my $pjfolder = $INDEX->{'resdirs'}->{$f};
				if ( $pjfolder eq $newpjxfoldername )
				{
					## move back to root folder
					$INDEX->{'resdirs'}->{$f} = 'root';
					$statmsg .= "Moved result folder <b>$f</b> back to <b>root</b>.<br>";
				}
			}
			## delete the folder
			delete $INDEX->{'projectfolders'}->{$newpjxfoldername};
			$statmsg .= "Folder <b>$newpjxfoldername</b> deleted.<br>";
		} else
		{
			$statmsg .= "Folder <b>$newpjxfoldername</b> does not exist.<br>";
		}
	}
## Store index after update
	store $INDEX, $indexfilename;
	return $statmsg;
}

sub update_index
{
	my $directoryhash = shift;
	my $INDEX         = shift;
	my $statmsg;
## First check if new folders were added
	my $currentdirectories = $directoryhash->{'resdirectories'};
	my $c                  = 0;
	foreach my $rsdir (@$currentdirectories)
	{
		unless ( $INDEX->{'resdirs'}->{$rsdir} )
		{
## Add a directory to the index, set root as projectfolder
			$INDEX->{'resdirs'}->{$rsdir} = "root";
			$c++;
		}
	}
	if ( $c > 0 )
	{
		if ( $c == 1 )
		{
			$statmsg .= "$c result folder was added to the index.<br>";
		} else
		{
			$statmsg .= "$c result folders were added to the index.<br>";
		}
	}
## Check if some folders have been removed, also remove them from the index
	my @indexedfolders = keys %{ $INDEX->{'resdirs'} };
	$c = 0;
	foreach my $rsdir (@indexedfolders)
	{
		unless ( $directoryhash->{'dirs'}->{$rsdir} )
		{
## remove from INDEX
			delete( $INDEX->{'resdirs'}->{$rsdir} );
			$c++;
		}
	}
	if ( $c > 0 )
	{
		if ( $c == 1 )
		{
			$statmsg .= "$c result folder was removed from the index.<br>";
		} else
		{
			$statmsg .= "$c result folders were removed from the index.<br>";
		}
	}

	#print "Removed $c directories from index.<br>";
	@indexedfolders = keys %{ $INDEX->{'resdirs'} };

	# print "Currently indexed folders: ", scalar(@indexedfolders), "<br>";
## ADD the root folder if not set yet
	unless ( $INDEX->{'projectfolders'}->{'root'} )
	{
		$INDEX->{'projectfolders'}->{'root'} = 1;
	}
## ADD the TRASH folder if it does not exists
	unless ( $INDEX->{'projectfolders'}->{'Trash'} )
	{
		$INDEX->{'projectfolders'}->{'Trash'} = 1;
	}
## Store index after update
	store $INDEX, $indexfilename;
	return $statmsg;
}

sub read_dirs
{
	my $directoryhash;
	my $datetranslationhash = {};
	opendir( DIR, $resultdirectory );
	my @directories = readdir(DIR);
	closedir(DIR);
	foreach my $resultdir ( sort @directories )
	{
		unless ( $resultdir =~ /\./ || $resultdir =~ /\.\./ )
		{
			## Concat the full path
			my $fullpath     = File::Spec->catfile( $resultdirectory, $resultdir );
			my $statobject   = ( stat $fullpath );                                    # gives back a statobject
			my $datemodified = $statobject->mtime;                                    # is the epoch time
			my ( $sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst ) = localtime($datemodified);    # convert epoch time
			## Create a number that can be sorted
			my $y = $year + 1900;
			my ( $m, $d );
			if ( ( $mon + 1 ) >= 10 )
			{
				$m = $mon + 1;
			} else
			{
				$m = "0" . ( $mon + 1 );
			}
			if ( ($mday) >= 10 )
			{
				$d = $mday;
			} else
			{
				$d = "0" . ($mday);
			}
			my $sortabledate = $y . $m . $d;

			# Generate a human readable format
			my $humanreadabledate = $d . "-" . $m . "-" . $y;
			$datetranslationhash->{$sortabledate} = $humanreadabledate;

			# Generate URL
			my $url = $WEBPARAMS->{'xmlparser'} . "?id=" . "$relpathtobase/$resultdir";
			$directoryhash->{'dirs'}->{$resultdir} = $url;

			# Deffile URLS
			my $xqdeffileurl = $WEBPARAMS->{'webdir'} . "dispdef.cgi?id=" . "$relpathtobase/$resultdir/xquest.def";
			$directoryhash->{'xqdeffileurl'}->{$resultdir} = $xqdeffileurl;
			my $xpdeffileurl = $WEBPARAMS->{'webdir'} . "dispdef.cgi?id=" . "$relpathtobase/$resultdir/xproph.def";
			$directoryhash->{'xpdeffileurl'}->{$resultdir} = $xpdeffileurl;
			push @{ $directoryhash->{'resdirectories'} }, $resultdir;    ## already sorted by name
			push @{ $directoryhash->{'timemod'}->{$sortabledate} }, $resultdir;    ##
			$directoryhash->{'time'}->{$resultdir} = $humanreadabledate;
		}
	}
## Sorting folders by time
## get all the times
	my @timesdesc = sort { $b <=> $a } ( keys %$datetranslationhash );
	my @resultfoldersdesc;
## Get all resdirs
	foreach my $date (@timesdesc)
	{
## get all resultfolders with this date
		my $resfolderarray = $directoryhash->{'timemod'}->{$date};
		foreach my $folder (@$resfolderarray)
		{
			push @resultfoldersdesc, $folder;
		}
	}
	$directoryhash->{'resdirectoriesdesctime'} = \@resultfoldersdesc;

	#foreach my $resfolder (@resultfoldersdesc){
	#print "$resfolder last mod: ".$directoryhash->{'time'}->{$resfolder}."<br>";
	#}
	return $directoryhash;
}

sub print_options_form
{

	#---------------------------------------------------------------------------
	#  GENERATE THE FORM FOR FILTERING AND SORTING
	#---------------------------------------------------------------------------
	print start_multipart_form( -name => 'cgivalues' );
	print '<TABLE BORDER=0 CELLSPACING=1 CELLPADDING=3>';
	print '<TR BGCOLOR=#EEEEFF NOWRAP>';
	print '<th style="color: #0000FF; font-size: 12pt" align="left" colspan="4">Options</td>';
	print '</TR>';

	#---------------------------------------------------------------------------
	#  ROW 1
	#---------------------------------------------------------------------------
	print '<TR BGCOLOR=#EEEEFF NOWRAP>';
	print '<TD>';
	print "Create/delete a project folder:</TD><TD> ",
	  $form->textfield(
						-name      => 'newpjxfoldername',
						-value     => $newpjxfoldername,
						-size      => 20,
						-maxlength => 20
	  );

	#print $form->br;
	#print '</TD>';
	my @values = [ 'create', 'delete' ];
	print "",
	  $form->popup_menu(
						 -name    => 'faction',
						 -values  => @values,
						 -default => 'create',
	  );
	print '</TD>';

	#print '</TR>';
	#---------------------------------------------------------------------------
	#  ROW 2
	#---------------------------------------------------------------------------
	#print '<TR BGCOLOR=#EEEEFF NOWRAP>';
	print '<TD>';
	## GET THE project folders from the INDEX
	@values = keys %{ $INDEX->{'projectfolders'} };
	## sort by name
	@values = sort(@values);
	print "Move selected result(s) to:</TD><TD> ",
	  $form->popup_menu(
						 -name    => 'targetfolder',
						 -values  => \@values,
						 -default => 'root',
	  );
	print '</TD>';
	print '</TR>';

	#  ROW 2
	#---------------------------------------------------------------------------
	print '<TR BGCOLOR=#EEEEFF NOWRAP>';

	#---------------------------------------------------------------------------
	#  ROW 2
	#---------------------------------------------------------------------------
	print '<TD>';
	## GET THE project folders from the INDEX
	@values = keys %{ $INDEX->{'projectfolders'} };
	## sort by name
	@values = sort(@values);
	print "Display project folder:</TD><TD> ",
	  $form->popup_menu(
						 -name    => 'dispfolder',
						 -values  => \@values,
						 -default => 'root',
	  );
	print '</TD>';
	print '<TD>';
	## Sortby
	@values = [ 'date', 'name' ];
	print "Sort results by:</TD><TD> ",
	  $form->popup_menu(
						 -name    => 'sortby',
						 -values  => @values,
						 -default => 'date',
	  );
	print '</TD>';
	print '</TR>';

	#---------------------------------------------------------------------------
	#  ROW 3
	#---------------------------------------------------------------------------
	print '<TR BGCOLOR=#EEEEFF NOWRAP>';
	print '<TD>';
	print "Refresh:</TD><TD> ",
	  $form->submit( -name  => 'submit_form',
					 -value => 'Update', );
	print '</TD>';
	print '<TD>';
	if ( $dispfolder eq "Trash" )
	{
		print $form->submit(
							 -name    => 'del_form',
							 -value   => 'Delete selected',
							 -onClick => "checkBoxDelete()",
		);
	}
	print '</TD>';
	print '<TD>';
	print '</TD>';
	print '</TR>';

	#---------------------------------------------------------------------------
	#  END OF TABLE
	#---------------------------------------------------------------------------
	print '</TABLE>';

	#---------------------------------------------------------------------------
	#  Hidden field
	#---------------------------------------------------------------------------
	print $form->hidden( 'jsvalue', 0 );

	#---------------------------------------------------------------------------
	#  Hidden fields ->resultxmlfilename
	#---------------------------------------------------------------------------
	#print $form->hidden( 'resultxml', basename($xmlfilename) );
	#---------------------------------------------------------------------------
	#  END THE FORM
	#---------------------------------------------------------------------------
}

sub readtables
{
	my $webparam = shift;
	my %WEBPARAMS;
	open WEBCONFIG, "<$webparam" or warn "could not open web config file $! ignoring";
	while (<WEBCONFIG>)
	{
		chomp;
		my @keyvalue = split /::/;
		$WEBPARAMS{ $keyvalue[0] } = $keyvalue[1];
	}
	return \%WEBPARAMS;
}

sub print_header
{
	print $form->header('text/html');
## JS Code
	my $jscode = qq|
function checkBoxDelete() {
if (confirm('Are you sure that you want to delete the selected folders, this action is irreversible.')) {
oFormObject = document.forms['cgivalues'];
oFormElement = oFormObject.elements["jsvalue"];
oFormObject.elements["jsvalue"].value = '1';
} else {
oFormObject = document.forms['cgivalues'];
oFormElement = oFormObject.elements["jsvalue"];
oFormObject.elements["jsvalue"].value = '0';}
}|;
	print $form->start_html(
							 -title   => 'xQuest/xProphet results manager',
							 -author  => 'Thomas Walzthoeni',
							 -expires => '+30s',
							 -meta    => { 'keywords' => 'MS cross-link xQuest xProhet', },
							 -style   => { 'src' => $WEBPARAMS->{'css_stylesheet'}, },
							 -script  => $jscode,
	);
}

sub debug_param
{

	#my $var=self;
	my $var = shift;

	# for debuging of param()
	my @parameters = $var->param();
	@parameters = sort (@parameters);
	print "<table border=\"1\">";
	print "<tr><td>name</td><td>value</td></tr>";
	for ( my $i = 0 ; $i < @parameters ; $i++ )
	{
		print "<tr><td>";
		print "$parameters[$i] </td>";
		print "<td>" . $var->param( $parameters[$i] ) . "</td></tr>";
	}
	print "</table>";
	return;
}
