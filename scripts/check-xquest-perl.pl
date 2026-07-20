#!/usr/bin/env perl
# Fail-fast check that the Perl environment can compile the xQuest search path
# used by GlycoQuest (compare_peaks3.pl → xquest.pl → Index.pm).
#
# Usage:
#   scripts/check-xquest-perl.pl --xquest-root /path/to/V2.1.7/xquest
#   scripts/check-xquest-perl.pl   # defaults to <repo>/V2.1.7/xquest
#
# Honors PERL5OPT / PERL5LIB already set by the caller (e.g. scripts/run.sh),
# but strips any .../1209/lib64/... entries (wrong Perl ABI).
# Exit 0 on success; non-zero with install hints on failure.

use strict;
use warnings;
use File::Basename qw(dirname);
use File::Spec;
use Cwd qw(abs_path);

my $script_dir = abs_path(dirname(__FILE__));
my $repo_root  = abs_path(File::Spec->catdir($script_dir, '..'));
my $xquest_root;

while (@ARGV) {
    my $arg = shift @ARGV;
    if ( $arg eq '--xquest-root' ) {
        $xquest_root = shift @ARGV
          or die "error: --xquest-root requires a path\n";
    }
    elsif ( $arg eq '-h' || $arg eq '--help' ) {
        print <<"EOF";
Usage: check-xquest-perl.pl [--xquest-root DIR]

Compile-checks GlycoQuest's xQuest search scripts and required modules.
EOF
        exit 0;
    }
    else {
        die "error: unknown argument: $arg\n";
    }
}

$xquest_root //= File::Spec->catdir( $repo_root, 'V2.1.7', 'xquest' );
$xquest_root = abs_path($xquest_root)
  or die "error: xQuest root not found: $xquest_root\n";

my $compare = File::Spec->catfile( $xquest_root, 'bin', 'compare_peaks3.pl' );
my $xquest  = File::Spec->catfile( $xquest_root, 'bin', 'xquest.pl' );
for my $path ( $compare, $xquest ) {
    -f $path or die "error: missing xQuest script: $path\n";
}

# Never use 1209/lib64 — legacy XS .so files crash modern Perl
# (e.g. undefined symbol: Perl_Gthr_key_ptr).
sub strip_lib64 {
    my ($path) = @_;
    return () unless defined $path && length $path;
    return grep { length && $_ !~ m{/1209/lib64(?:/|$)} } split /:/, $path;
}

# Match job PERL5LIB (src/xquest/perl_deps.rs::xquest_perl5lib).
my @lib_dirs = (
    File::Spec->catdir( $xquest_root, '1209', 'lib',   'perl5' ),
    File::Spec->catdir( $xquest_root, '1209', 'share', 'perl5' ),
    File::Spec->catdir( $xquest_root, 'modules' ),
);
{
    require Config;
    my $arch = $Config::Config{archname} // '';
    if ($arch) {
        my $arch_lib =
          File::Spec->catdir( $xquest_root, '1209', 'lib', 'perl5', $arch );
        unshift @lib_dirs, $arch_lib if -d $arch_lib;
    }
}
if ( defined $ENV{HOME} && -d "$ENV{HOME}/perl5/lib/perl5" ) {
    unshift @lib_dirs, "$ENV{HOME}/perl5/lib/perl5";
}

my @existing_libs = grep { -d $_ } @lib_dirs;
require lib;
lib->import(@existing_libs);

my $perl5lib = join( ':', @existing_libs, strip_lib64( $ENV{PERL5LIB} ) );
$ENV{PERL5LIB} = $perl5lib;

# Strip -I.../1209/lib64 from PERL5OPT if a previous run added them.
if ( defined $ENV{PERL5OPT} ) {
    my @opts = grep { $_ !~ m{^-I.*/1209/lib64} } split /\s+/, $ENV{PERL5OPT};
    $ENV{PERL5OPT} = join( ' ', @opts );
    delete $ENV{PERL5OPT} unless length $ENV{PERL5OPT};
}

# Required for GlycoQuest batch search. GD is optional (soft-loaded; drawspectra 0).
my @critical = qw(
  DB_File
  MLDBM
  Bio::Perl
  XML::TreeBuilder
  XML::Parser
  XML::Element
  HTML::Tagset
  HTML::Entities
  MIME::Base64
  Storable
  Statistics::Descriptive
);

my @missing;
my @load_errors;
for my $mod (@critical) {
    ( my $fn = "$mod.pm" ) =~ s{::}{/}g;
    eval {
        require $fn;
        1;
    } or do {
        my $err = $@ // 'unknown error';
        chomp $err;
        push @missing, $mod;
        push @load_errors, "$mod: $err";
    };
}

my @compile_errors;
for my $script ( $compare, $xquest ) {
    my $out = qx{env PERL5LIB=\Q$ENV{PERL5LIB}\E perl -c \Q$script\E 2>&1};
    my $rc  = $? >> 8;
    if ( $rc != 0 || $out !~ /syntax OK/ ) {
        chomp $out;
        push @compile_errors, "$script\n$out";
    }
}

if ( !@missing && !@compile_errors ) {
    print "xQuest Perl search path OK ($xquest_root)\n";
    print "  modules: " . join( ', ', @critical ) . "\n";
    print "  compile: compare_peaks3.pl, xquest.pl\n";
    print "  note: GD optional (spectrum PNGs); GlycoQuest uses drawspectra 0\n";
    exit 0;
}

print STDERR "error: xQuest Perl environment is incomplete (failing before job launch).\n";
print STDERR "  xquest-root: $xquest_root\n";
print STDERR "  PERL5LIB=$ENV{PERL5LIB}\n";
print STDERR "  PERL5OPT=" . ( $ENV{PERL5OPT} // '<unset>' ) . "\n";
print STDERR "  LD_LIBRARY_PATH=" . ( $ENV{LD_LIBRARY_PATH} // '<unset>' ) . "\n";

if (@load_errors) {
    print STDERR "\nModule load failures:\n";
    print STDERR "  $_\n" for @load_errors;
}

if (@compile_errors) {
    print STDERR "\nCompile check failures:\n";
    print STDERR "$_\n\n" for @compile_errors;
}

print STDERR <<"EOF";

Install once (login node), then re-run this check:

  scripts/bootstrap-euler-perl.sh
  # IMPORTANT: do NOT add V2.1.7/xquest/1209/lib64 to PERL5LIB
  # (legacy XS .so files break modern Perl).
  scripts/check-xquest-perl.pl --xquest-root $xquest_root

XS runtime libs: libdb (DB_File), libexpat (XML::Parser).
GD is optional for GlycoQuest (drawspectra 0).

Fedora/RHEL: perl-DB_File perl-XML-Parser
Debian/Ubuntu: libdb-file-perl libxml-parser-perl
EOF

exit 1;
