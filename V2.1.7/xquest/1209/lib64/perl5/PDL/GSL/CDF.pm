
#
# GENERATED WITH PDL::PP! Don't modify!
#
package PDL::GSL::CDF;

@EXPORT_OK  = qw( PDL::PP gsl_cdf_ugaussian_P PDL::PP gsl_cdf_ugaussian_Pinv PDL::PP gsl_cdf_gaussian_P PDL::PP gsl_cdf_gaussian_Pinv PDL::PP gsl_cdf_gamma_P PDL::PP gsl_cdf_gamma_Pinv PDL::PP gsl_cdf_cauchy_P PDL::PP gsl_cdf_cauchy_Pinv PDL::PP gsl_cdf_laplace_P PDL::PP gsl_cdf_laplace_Pinv PDL::PP gsl_cdf_rayleigh_P PDL::PP gsl_cdf_rayleigh_Pinv PDL::PP gsl_cdf_chisq_P PDL::PP gsl_cdf_chisq_Pinv PDL::PP gsl_cdf_exponential_P PDL::PP gsl_cdf_exponential_Pinv PDL::PP gsl_cdf_exppow_P PDL::PP gsl_cdf_tdist_P PDL::PP gsl_cdf_tdist_Pinv PDL::PP gsl_cdf_fdist_P PDL::PP gsl_cdf_fdist_Pinv PDL::PP gsl_cdf_beta_P PDL::PP gsl_cdf_beta_Pinv PDL::PP gsl_cdf_flat_P PDL::PP gsl_cdf_flat_Pinv PDL::PP gsl_cdf_lognormal_P PDL::PP gsl_cdf_lognormal_Pinv PDL::PP gsl_cdf_gumbel1_P PDL::PP gsl_cdf_gumbel1_Pinv PDL::PP gsl_cdf_gumbel2_P PDL::PP gsl_cdf_gumbel2_Pinv PDL::PP gsl_cdf_weibull_P PDL::PP gsl_cdf_weibull_Pinv PDL::PP gsl_cdf_pareto_P PDL::PP gsl_cdf_pareto_Pinv PDL::PP gsl_cdf_logistic_P PDL::PP gsl_cdf_logistic_Pinv PDL::PP gsl_cdf_binomial_P PDL::PP gsl_cdf_poisson_P PDL::PP gsl_cdf_geometric_P PDL::PP gsl_cdf_negative_binomial_P PDL::PP gsl_cdf_pascal_P PDL::PP gsl_cdf_hypergeometric_P );
%EXPORT_TAGS = (Func=>[@EXPORT_OK]);

use PDL::Core;
use PDL::Exporter;
use DynaLoader;



   
   @ISA    = ( 'PDL::Exporter','DynaLoader' );
   push @PDL::Core::PP, __PACKAGE__;
   bootstrap PDL::GSL::CDF ;





$PDL::onlinedoc->scan(__FILE__) if $PDL::onlinedoc;

=head1 NAME

PDL::GSL::CDF - PDL interface to GSL Cumulative Distribution Functions

=head1 DESCRIPTION

This is an interface to the Cumulative Distribution Function package present in the GNU Scientific Library. 

=head1 SYNOPSIS

    use PDL;
    use PDL::GSL::CDF;

    my $p = gsl_cdf_tdist_P( $t, $df );

    my $t = gsl_cdf_tdist_Pinv( $p, $df );

=cut







=head1 FUNCTIONS



=cut






=head2 gsl_cdf_ugaussian_P

=for sig

  Signature: (double x();  [o]out())

=for ref



=for bad

gsl_cdf_ugaussian_P does handle bad values.
It will set the bad-value flag of all output piddles if the flag is set for any of the input piddles.


=cut






*gsl_cdf_ugaussian_P = \&PDL::gsl_cdf_ugaussian_P;





=head2 gsl_cdf_ugaussian_Pinv

=for sig

  Signature: (double p();  [o]out())

=for ref



=for bad

gsl_cdf_ugaussian_Pinv does handle bad values.
It will set the bad-value flag of all output piddles if the flag is set for any of the input piddles.


=cut






*gsl_cdf_ugaussian_Pinv = \&PDL::gsl_cdf_ugaussian_Pinv;





=head2 gsl_cdf_gaussian_P

=for sig

  Signature: (double x(); double sigma();  [o]out())

=for ref



=for bad

gsl_cdf_gaussian_P does handle bad values.
It will set the bad-value flag of all output piddles if the flag is set for any of the input piddles.


=cut






*gsl_cdf_gaussian_P = \&PDL::gsl_cdf_gaussian_P;





=head2 gsl_cdf_gaussian_Pinv

=for sig

  Signature: (double p(); double sigma();  [o]out())

=for ref



=for bad

gsl_cdf_gaussian_Pinv does handle bad values.
It will set the bad-value flag of all output piddles if the flag is set for any of the input piddles.


=cut






*gsl_cdf_gaussian_Pinv = \&PDL::gsl_cdf_gaussian_Pinv;





=head2 gsl_cdf_gamma_P

=for sig

  Signature: (double x(); double a(); double b();  [o]out())

=for ref



=for bad

gsl_cdf_gamma_P does handle bad values.
It will set the bad-value flag of all output piddles if the flag is set for any of the input piddles.


=cut






*gsl_cdf_gamma_P = \&PDL::gsl_cdf_gamma_P;





=head2 gsl_cdf_gamma_Pinv

=for sig

  Signature: (double p(); double a(); double b();  [o]out())

=for ref



=for bad

gsl_cdf_gamma_Pinv does handle bad values.
It will set the bad-value flag of all output piddles if the flag is set for any of the input piddles.


=cut






*gsl_cdf_gamma_Pinv = \&PDL::gsl_cdf_gamma_Pinv;





=head2 gsl_cdf_cauchy_P

=for sig

  Signature: (double x(); double a();  [o]out())

=for ref



=for bad

gsl_cdf_cauchy_P does handle bad values.
It will set the bad-value flag of all output piddles if the flag is set for any of the input piddles.


=cut






*gsl_cdf_cauchy_P = \&PDL::gsl_cdf_cauchy_P;





=head2 gsl_cdf_cauchy_Pinv

=for sig

  Signature: (double p(); double a();  [o]out())

=for ref



=for bad

gsl_cdf_cauchy_Pinv does handle bad values.
It will set the bad-value flag of all output piddles if the flag is set for any of the input piddles.


=cut






*gsl_cdf_cauchy_Pinv = \&PDL::gsl_cdf_cauchy_Pinv;





=head2 gsl_cdf_laplace_P

=for sig

  Signature: (double x(); double a();  [o]out())

=for ref



=for bad

gsl_cdf_laplace_P does handle bad values.
It will set the bad-value flag of all output piddles if the flag is set for any of the input piddles.


=cut






*gsl_cdf_laplace_P = \&PDL::gsl_cdf_laplace_P;





=head2 gsl_cdf_laplace_Pinv

=for sig

  Signature: (double p(); double a();  [o]out())

=for ref



=for bad

gsl_cdf_laplace_Pinv does handle bad values.
It will set the bad-value flag of all output piddles if the flag is set for any of the input piddles.


=cut






*gsl_cdf_laplace_Pinv = \&PDL::gsl_cdf_laplace_Pinv;





=head2 gsl_cdf_rayleigh_P

=for sig

  Signature: (double x(); double sigma();  [o]out())

=for ref



=for bad

gsl_cdf_rayleigh_P does handle bad values.
It will set the bad-value flag of all output piddles if the flag is set for any of the input piddles.


=cut






*gsl_cdf_rayleigh_P = \&PDL::gsl_cdf_rayleigh_P;





=head2 gsl_cdf_rayleigh_Pinv

=for sig

  Signature: (double p(); double sigma();  [o]out())

=for ref



=for bad

gsl_cdf_rayleigh_Pinv does handle bad values.
It will set the bad-value flag of all output piddles if the flag is set for any of the input piddles.


=cut






*gsl_cdf_rayleigh_Pinv = \&PDL::gsl_cdf_rayleigh_Pinv;





=head2 gsl_cdf_chisq_P

=for sig

  Signature: (double x(); double nu();  [o]out())

=for ref



=for bad

gsl_cdf_chisq_P does handle bad values.
It will set the bad-value flag of all output piddles if the flag is set for any of the input piddles.


=cut






*gsl_cdf_chisq_P = \&PDL::gsl_cdf_chisq_P;





=head2 gsl_cdf_chisq_Pinv

=for sig

  Signature: (double p(); double nu();  [o]out())

=for ref



=for bad

gsl_cdf_chisq_Pinv does handle bad values.
It will set the bad-value flag of all output piddles if the flag is set for any of the input piddles.


=cut






*gsl_cdf_chisq_Pinv = \&PDL::gsl_cdf_chisq_Pinv;





=head2 gsl_cdf_exponential_P

=for sig

  Signature: (double x(); double mu();  [o]out())

=for ref



=for bad

gsl_cdf_exponential_P does handle bad values.
It will set the bad-value flag of all output piddles if the flag is set for any of the input piddles.


=cut






*gsl_cdf_exponential_P = \&PDL::gsl_cdf_exponential_P;





=head2 gsl_cdf_exponential_Pinv

=for sig

  Signature: (double p(); double mu();  [o]out())

=for ref



=for bad

gsl_cdf_exponential_Pinv does handle bad values.
It will set the bad-value flag of all output piddles if the flag is set for any of the input piddles.


=cut






*gsl_cdf_exponential_Pinv = \&PDL::gsl_cdf_exponential_Pinv;





=head2 gsl_cdf_exppow_P

=for sig

  Signature: (double x(); double a(); double b();  [o]out())

=for ref



=for bad

gsl_cdf_exppow_P does handle bad values.
It will set the bad-value flag of all output piddles if the flag is set for any of the input piddles.


=cut






*gsl_cdf_exppow_P = \&PDL::gsl_cdf_exppow_P;





=head2 gsl_cdf_tdist_P

=for sig

  Signature: (double x(); double nu();  [o]out())

=for ref



=for bad

gsl_cdf_tdist_P does handle bad values.
It will set the bad-value flag of all output piddles if the flag is set for any of the input piddles.


=cut






*gsl_cdf_tdist_P = \&PDL::gsl_cdf_tdist_P;





=head2 gsl_cdf_tdist_Pinv

=for sig

  Signature: (double p(); double nu();  [o]out())

=for ref



=for bad

gsl_cdf_tdist_Pinv does handle bad values.
It will set the bad-value flag of all output piddles if the flag is set for any of the input piddles.


=cut






*gsl_cdf_tdist_Pinv = \&PDL::gsl_cdf_tdist_Pinv;





=head2 gsl_cdf_fdist_P

=for sig

  Signature: (double x(); double nua(); double nub();  [o]out())

=for ref



=for bad

gsl_cdf_fdist_P does handle bad values.
It will set the bad-value flag of all output piddles if the flag is set for any of the input piddles.


=cut






*gsl_cdf_fdist_P = \&PDL::gsl_cdf_fdist_P;





=head2 gsl_cdf_fdist_Pinv

=for sig

  Signature: (double p(); double nua(); double nub();  [o]out())

=for ref



=for bad

gsl_cdf_fdist_Pinv does handle bad values.
It will set the bad-value flag of all output piddles if the flag is set for any of the input piddles.


=cut






*gsl_cdf_fdist_Pinv = \&PDL::gsl_cdf_fdist_Pinv;





=head2 gsl_cdf_beta_P

=for sig

  Signature: (double x(); double a(); double b();  [o]out())

=for ref



=for bad

gsl_cdf_beta_P does handle bad values.
It will set the bad-value flag of all output piddles if the flag is set for any of the input piddles.


=cut






*gsl_cdf_beta_P = \&PDL::gsl_cdf_beta_P;





=head2 gsl_cdf_beta_Pinv

=for sig

  Signature: (double p(); double a(); double b();  [o]out())

=for ref



=for bad

gsl_cdf_beta_Pinv does handle bad values.
It will set the bad-value flag of all output piddles if the flag is set for any of the input piddles.


=cut






*gsl_cdf_beta_Pinv = \&PDL::gsl_cdf_beta_Pinv;





=head2 gsl_cdf_flat_P

=for sig

  Signature: (double x(); double a(); double b();  [o]out())

=for ref



=for bad

gsl_cdf_flat_P does handle bad values.
It will set the bad-value flag of all output piddles if the flag is set for any of the input piddles.


=cut






*gsl_cdf_flat_P = \&PDL::gsl_cdf_flat_P;





=head2 gsl_cdf_flat_Pinv

=for sig

  Signature: (double p(); double a(); double b();  [o]out())

=for ref



=for bad

gsl_cdf_flat_Pinv does handle bad values.
It will set the bad-value flag of all output piddles if the flag is set for any of the input piddles.


=cut






*gsl_cdf_flat_Pinv = \&PDL::gsl_cdf_flat_Pinv;





=head2 gsl_cdf_lognormal_P

=for sig

  Signature: (double x(); double zeta(); double sigma();  [o]out())

=for ref



=for bad

gsl_cdf_lognormal_P does handle bad values.
It will set the bad-value flag of all output piddles if the flag is set for any of the input piddles.


=cut






*gsl_cdf_lognormal_P = \&PDL::gsl_cdf_lognormal_P;





=head2 gsl_cdf_lognormal_Pinv

=for sig

  Signature: (double p(); double zeta(); double sigma();  [o]out())

=for ref



=for bad

gsl_cdf_lognormal_Pinv does handle bad values.
It will set the bad-value flag of all output piddles if the flag is set for any of the input piddles.


=cut






*gsl_cdf_lognormal_Pinv = \&PDL::gsl_cdf_lognormal_Pinv;





=head2 gsl_cdf_gumbel1_P

=for sig

  Signature: (double x(); double a(); double b();  [o]out())

=for ref



=for bad

gsl_cdf_gumbel1_P does handle bad values.
It will set the bad-value flag of all output piddles if the flag is set for any of the input piddles.


=cut






*gsl_cdf_gumbel1_P = \&PDL::gsl_cdf_gumbel1_P;





=head2 gsl_cdf_gumbel1_Pinv

=for sig

  Signature: (double p(); double a(); double b();  [o]out())

=for ref



=for bad

gsl_cdf_gumbel1_Pinv does handle bad values.
It will set the bad-value flag of all output piddles if the flag is set for any of the input piddles.


=cut






*gsl_cdf_gumbel1_Pinv = \&PDL::gsl_cdf_gumbel1_Pinv;





=head2 gsl_cdf_gumbel2_P

=for sig

  Signature: (double x(); double a(); double b();  [o]out())

=for ref



=for bad

gsl_cdf_gumbel2_P does handle bad values.
It will set the bad-value flag of all output piddles if the flag is set for any of the input piddles.


=cut






*gsl_cdf_gumbel2_P = \&PDL::gsl_cdf_gumbel2_P;





=head2 gsl_cdf_gumbel2_Pinv

=for sig

  Signature: (double p(); double a(); double b();  [o]out())

=for ref



=for bad

gsl_cdf_gumbel2_Pinv does handle bad values.
It will set the bad-value flag of all output piddles if the flag is set for any of the input piddles.


=cut






*gsl_cdf_gumbel2_Pinv = \&PDL::gsl_cdf_gumbel2_Pinv;





=head2 gsl_cdf_weibull_P

=for sig

  Signature: (double x(); double a(); double b();  [o]out())

=for ref



=for bad

gsl_cdf_weibull_P does handle bad values.
It will set the bad-value flag of all output piddles if the flag is set for any of the input piddles.


=cut






*gsl_cdf_weibull_P = \&PDL::gsl_cdf_weibull_P;





=head2 gsl_cdf_weibull_Pinv

=for sig

  Signature: (double p(); double a(); double b();  [o]out())

=for ref



=for bad

gsl_cdf_weibull_Pinv does handle bad values.
It will set the bad-value flag of all output piddles if the flag is set for any of the input piddles.


=cut






*gsl_cdf_weibull_Pinv = \&PDL::gsl_cdf_weibull_Pinv;





=head2 gsl_cdf_pareto_P

=for sig

  Signature: (double x(); double a(); double b();  [o]out())

=for ref



=for bad

gsl_cdf_pareto_P does handle bad values.
It will set the bad-value flag of all output piddles if the flag is set for any of the input piddles.


=cut






*gsl_cdf_pareto_P = \&PDL::gsl_cdf_pareto_P;





=head2 gsl_cdf_pareto_Pinv

=for sig

  Signature: (double p(); double a(); double b();  [o]out())

=for ref



=for bad

gsl_cdf_pareto_Pinv does handle bad values.
It will set the bad-value flag of all output piddles if the flag is set for any of the input piddles.


=cut






*gsl_cdf_pareto_Pinv = \&PDL::gsl_cdf_pareto_Pinv;





=head2 gsl_cdf_logistic_P

=for sig

  Signature: (double x(); double a();  [o]out())

=for ref



=for bad

gsl_cdf_logistic_P does handle bad values.
It will set the bad-value flag of all output piddles if the flag is set for any of the input piddles.


=cut






*gsl_cdf_logistic_P = \&PDL::gsl_cdf_logistic_P;





=head2 gsl_cdf_logistic_Pinv

=for sig

  Signature: (double p(); double a();  [o]out())

=for ref



=for bad

gsl_cdf_logistic_Pinv does handle bad values.
It will set the bad-value flag of all output piddles if the flag is set for any of the input piddles.


=cut






*gsl_cdf_logistic_Pinv = \&PDL::gsl_cdf_logistic_Pinv;





=head2 gsl_cdf_binomial_P

=for sig

  Signature: (ushort k(); double p(); ushort n();  [o]out())

=for ref



=for bad

gsl_cdf_binomial_P does handle bad values.
It will set the bad-value flag of all output piddles if the flag is set for any of the input piddles.


=cut






*gsl_cdf_binomial_P = \&PDL::gsl_cdf_binomial_P;





=head2 gsl_cdf_poisson_P

=for sig

  Signature: (ushort k(); double mu();  [o]out())

=for ref



=for bad

gsl_cdf_poisson_P does handle bad values.
It will set the bad-value flag of all output piddles if the flag is set for any of the input piddles.


=cut






*gsl_cdf_poisson_P = \&PDL::gsl_cdf_poisson_P;





=head2 gsl_cdf_geometric_P

=for sig

  Signature: (ushort k(); double p();  [o]out())

=for ref



=for bad

gsl_cdf_geometric_P does handle bad values.
It will set the bad-value flag of all output piddles if the flag is set for any of the input piddles.


=cut






*gsl_cdf_geometric_P = \&PDL::gsl_cdf_geometric_P;





=head2 gsl_cdf_negative_binomial_P

=for sig

  Signature: (ushort k(); double p(); double n();  [o]out())

=for ref



=for bad

gsl_cdf_negative_binomial_P does handle bad values.
It will set the bad-value flag of all output piddles if the flag is set for any of the input piddles.


=cut






*gsl_cdf_negative_binomial_P = \&PDL::gsl_cdf_negative_binomial_P;





=head2 gsl_cdf_pascal_P

=for sig

  Signature: (ushort k(); double p(); ushort n();  [o]out())

=for ref



=for bad

gsl_cdf_pascal_P does handle bad values.
It will set the bad-value flag of all output piddles if the flag is set for any of the input piddles.


=cut






*gsl_cdf_pascal_P = \&PDL::gsl_cdf_pascal_P;





=head2 gsl_cdf_hypergeometric_P

=for sig

  Signature: (ushort k(); ushort na(); ushort nb(); ushort t();  [o]out())

=for ref



=for bad

gsl_cdf_hypergeometric_P does handle bad values.
It will set the bad-value flag of all output piddles if the flag is set for any of the input piddles.


=cut






*gsl_cdf_hypergeometric_P = \&PDL::gsl_cdf_hypergeometric_P;



;


=head1 AUTHOR

Copyright (C) 2009 Maggie J. Xiong <maggiexyz users.sourceforge.net>

The GSL CDF module was written by J. Stover.

All rights reserved. There is no warranty. You are allowed to redistribute this software / documentation as described in the file COPYING in the PDL distribution.

=cut





# Exit with OK status

1;

		   