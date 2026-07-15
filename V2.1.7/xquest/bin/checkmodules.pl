#!/usr/bin/env perl
use strict;

my @modules = qw (
Bio::Index::Fasta
MLDBM
CGI::Fast
Data::Dumper
HTML::Template
File::Copy
GD::Graph
IO::Compress::Base
MIME::Base64
CGI::Session
HTML::TreeBuilder
XML::Writer 
XML::Parser
XML::TreeBuilder
CGI::FastTemplate
HTML::PageIndex
Statistics::Descriptive
CGI::FormBuilder
Math::Random
);

foreach my $module (@modules){
(my $fn="$module.pm")=~s|::|/|g; # Foo::Bar::Baz => Foo/Bar/Baz.pm
if  (eval {require $fn;1;}) {
	#module loaded
	print "Module $fn installed\n";
}else{
	print "Module $fn -->not availaible\n";
}	
}
