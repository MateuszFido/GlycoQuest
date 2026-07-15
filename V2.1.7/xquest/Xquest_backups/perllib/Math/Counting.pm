# $Id$
package Math::Counting;
use strict;
use warnings;
use base 'Exporter';
use vars qw( @EXPORT_OK @EXPORT %EXPORT_TAGS $VERSION );
use Carp;
use Math::BigInt; #try => 'GMP';

$VERSION = 0.0704;
@EXPORT_OK = ();
@EXPORT = qw(
    f  factorial    bfact
    P  permutation  bperm
    C  combination  bcomb
);
%EXPORT_TAGS = (
    short => [qw( f P C )],
    long  => [qw( factorial permutation combination )],
    big   => [qw( bfact bperm bcomb )],
);

sub f { factorial(@_) }
sub factorial {
    my( $n ) = @_;
    my $product = 1;
    while( $n > 0 ) {
        $product *= $n--;
    }
    return $product;
}
# The right way to do it:
sub bfact {
    my $n = Math::BigInt->new( shift );
    return $n->bfac();
}

sub P { permutation(@_) }
sub permutation {
    my( $n, $r ) = @_;
    my $product = 1;
    while( $r > 0 ) {
        $product *= $n--;
        $r--;
    }
    return $product;
}
sub bperm {
    my( $n, $r ) = @_;
    $n = Math::BigInt->new( $n );
    $r = Math::BigInt->new( $r );
    $r = $n - $r;
    return $n->bfac() / $r->bfac();
}

sub C { combination(@_) }
sub combination {
    my( $n, $r ) = @_;
    my $product = 1;
    while( $r > 0 ) {
        $product *= $n--;
        $product /= $r--;
    }
    return $product;
}
sub bcomb {
    my( $n, $k ) = @_;
    $n = Math::BigInt->new( $n );
    $k = Math::BigInt->new( $k );
    my $r = $n - $k;
    return $n->bfac() / ($k->bfac() * $r->bfac());
}

1;

__END__

=head1 NAME

Math::Counting - Combinatorial counting operations

=head1 SYNOPSIS

  ## The algorithm-student versions:

  use Math::Counting ':long';
  printf "Given n=%d and r=%d:\nFact=%d\nPerm=%d\nComb=%d\n",
    $n, $r, factorial($n), permutation($n, $r), combination($n, $r);

  use Math::Counting ':short';
  printf "Given n=%d and r=%d:\nF=%d\nP=%d\nC=%d\n",
    $n, $r, f($n), P($n, $r), C($n, $r);

  ## The "Right way to do it":

  use Math::Counting ':big';
  printf "n=%d, r=%d:\nBig F=%d\n Big P=%d\nBig C=%d\n",
    $n, $r, bfact($n), bperm($n, $r), bcomb($n, $r);

=head1 DESCRIPTION

This package computes the factorial, number of (n,r) permutations and
(n,r) combinations using either L<Math::BigInt> and the technique of
"tail call elimination" (TCE).

No functions are exported by default.

=head1 FUNCTIONS

=head2 f, factorial

  $f = f($n);
  $f = factorial($n);

Return the number of arrangements of B<n> with TCE.

=head2 bfact

  $f = bfact($n);

Return the value of the L<Math::BigInt/bfact> function, which is the
"Right Way To Do It."

=head2 P, permutation

  $p = P($n, $r);
  $p = permutation($n, $r);

Return the number of arrangements of B<r> elements drawn from a set of
B<n> elements with TCE.  (B<nPn> is the same as B<n!>.)

=head2 bperm

  $p = bperm($n, $r);

Return the C<Math::BigInt> computation: B<n!/(n-r)!>.

=head2 C, choose, combination

  $c = C($n, $r);
  $c = choose($n, $r);
  $c = combination($n, $r);

Return the number of ways to choose B<r> elements from a set of B<n>
elements with TCE.

=head2 bcomb

  $c = bcomb($n, $r);

Return the C<Math::BigInt> computation: B<n!/r!(n-r)!>.

=head1 TO DO

Figure out how to allow the use of different C<BigInt> variations,
like C<GMP>.

Provide the gamma function for the factorial of non-integer numbers?

=head1 SEE ALSO

B<Higher Order Perl> by Mark Jason Dominus
(L<http://hop.perl.plover.com/>).

B<Mastering Algorithms with Perl> by Orwant, Hietaniemi & Macdonald
(L<http://www.oreilly.com/catalog/maperl/>).

The B<bfac> function of the C<Math::BigInt> module.

L<Algorithm::Combinatorics>

L<String::OrderedCombination> for B<nPk> list generation.

L<Math::Combinatorics> for B<nPn> and B<nCr> list generation.

L<http://en.wikipedia.org/wiki/Factorial>

L<http://en.wikipedia.org/wiki/Permutation>

L<http://en.wikipedia.org/wiki/Combination>

=head1 AUTHOR

Not me. I am but a pebble on the beach...

=head1 COPYRIGHT

Copyright 2007 Gene Boggs All Rights Reserved

You may use this package under the same terms as Perl itself.

=cut
