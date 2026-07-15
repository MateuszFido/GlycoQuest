# $File: //member/autrijus/Acme-Hello/lib/Acme/Hello.pm $ 
# $Revision: #2 $ $Change: 2772 $ $DateTime: 2002/12/19 03:19:17 $

package Acme::Hello;
$Acme::Hello::VERSION = '0.02';

use strict;
use Acme::Hello::I18N;

use Exporter;
use base 'Exporter';
use vars '@EXPORT';

@EXPORT = 'hello';

=head1 NAME

Acme::Hello - Print a greeting message

=head1 VERSION

This document describes version 0.02 of B<Acme::Hello>.

=head1 SYNOPSIS

    use Acme::Hello;	# exports hello() by default
    hello();		# procedure call interface

    my $obj = Acme::Hello->new;
    $obj->hello;	# object-oriented interface

=cut

sub new {
    my ($class, %args) = @_;
    $class = ref($class) if (ref $class);

    $args{lh} ||= Acme::Hello::I18N->get_handle($args{language})
        or die "Cannot find handle for language: $args{language}.\n";

    return bless(\%args, $class);
}

sub hello {
    my $self = ref($_[0]) ? $_[0] : __PACKAGE__->new;

    print $self->loc("Hello, world!"), "\n";
}

sub lh {
    my $self = shift;
    $self->{lh} = shift if @_;
    return $self->{lh};
}

sub loc {
    my $self = shift;
    return $self->lh->maketext(@_);
}


1;

__END__

=head1 AUTHORS

Autrijus Tang E<lt>autrijus@autrijus.orgE<gt>

=head1 COPYRIGHT

Copyright 2002 by Autrijus Tang E<lt>autrijus@autrijus.orgE<gt>.

This program is free software; you can redistribute it and/or 
modify it under the same terms as Perl itself.

See L<http://www.perl.com/perl/misc/Artistic.html>

=cut

# Local variables:
# c-indentation-style: bsd
# c-basic-offset: 4
# indent-tabs-mode: nil
# End:
# vim: expandtab shiftwidth=4:
