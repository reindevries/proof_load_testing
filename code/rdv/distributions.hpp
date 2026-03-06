// distributions.hpp
//
// Implements common statistical distribution functions. Some of these require the Boost C++
// library. In these cases the required include file is provided in the comments.
//
// Author: Rein de Vries
// Date: 21 August 2025

#ifndef __RDV_DISTRIBUTIONS__
#define __RDV_DISTRIBUTIONS__

#include "general.hpp"
#include "erf_inv.hpp"

namespace rdv {

// Gamma distribution
// https://en.wikipedia.org/wiki/Gamma_distribution
// #include <boost/math/special_functions/gamma.hpp>

#ifdef BOOST_MATH_SF_GAMMA_HPP
inline double gamma_pdf(double x, double k, double theta) {
	ASSERT(x >= 0.0);
	ASSERT(k > 0.0 && theta > 0.0);
	return boost::math::gamma_p_derivative(k, x / theta) / theta;
}

inline double gamma_cdf(double x, double k, double theta) {
	ASSERT(x >= 0.0);
	ASSERT(k > 0.0 && theta > 0.0);
	return boost::math::gamma_p(k, x / theta);
}

inline double gamma_sf(double x, double k, double theta) {
	ASSERT(x >= 0.0);
	ASSERT(k > 0.0 && theta > 0.0);
	return boost::math::gamma_q(k, x / theta);
}
	
inline double gamma_cdf_inv(double F, double k, double theta) {
	ASSERT(F >= 0.0 && F <= 1.0);
	ASSERT(k > 0.0 && theta > 0.0);
	if (F == 1.0)
		return M_INF;
	return theta * boost::math::gamma_p_inv(k, F);
}

inline double gamma_sf_inv(double S, double k, double theta) {
	ASSERT(S >= 0.0 && S <= 1.0);
	ASSERT(k > 0.0 && theta > 0.0);
	if (S == 0.0)
		return M_INF;
	return theta * boost::math::gamma_q_inv(k, S);
}

inline double gamma_quantile(double F, double k, double theta) {
	return gamma_cdf_inv(F, k, theta);
}
#endif

inline void gamma_ms(double mean, double stddev, double& k, double& theta) {
	ASSERT(mean > 0.0 && stddev > 0.0);
	k = sq(mean / stddev);
	theta = sq(stddev) / mean;
}

inline void gamma_mv(double mean, double coeff_of_var, double& k, double& theta) {
	ASSERT(coeff_of_var >= 0.0);
	gamma_ms(mean, abs(coeff_of_var * mean), k, theta);
}

inline double gamma_mean(double k, double theta) {
	ASSERT(k > 0.0 && theta > 0.0);
	return k * theta;
}

inline double gamma_variance(double k, double theta) {
	ASSERT(k > 0.0 && theta > 0.0);
	return k * sq(theta);
}

inline double gamma_stddev(double k, double theta) {
	return sqrt(gamma_variance(k, theta));
}

inline double gamma_cov(double k, double theta) {
	return gamma_stddev(k, theta) / gamma_mean(k, theta);
}

// Gumbel maximum distribution (extreme value type I)
// https://en.wikipedia.org/wiki/Gumbel_distribution

inline double gumbel_max_pdf(double x, double mu, double beta) {
	ASSERT(x == x && beta > 0.0);
	double z = (x - mu) / beta;
	return exp(-(z + exp(-z))) / beta;
}

inline double gumbel_max_log_pdf(double x, double mu, double beta) {
	ASSERT(x == x && beta > 0.0);
	double z = (x - mu) / beta;
	return -log(beta) - z - exp(-z);
}

inline double gumbel_max_cdf(double x, double mu, double beta) {
	ASSERT(x == x && beta > 0.0);
	double z = (x - mu) / beta;
	return exp(-exp(-z));
}

inline double gumbel_max_log_cdf(double x, double mu, double beta) {
	ASSERT(x == x && beta > 0.0);
	double z = (x - mu) / beta;
	return -exp(-z);
}

inline double gumbel_max_sf(double x, double mu, double beta) {
	ASSERT(x == x && beta > 0.0);
	double z = (x - mu) / beta;
	return -expm1(-exp(-z));
}
	
inline double gumbel_max_log_sf(double x, double mu, double beta) {
	ASSERT(x == x && beta > 0.0);
	double z = (x - mu) / beta;
	return log1p(-exp(-exp(-z)));
}

inline double gumbel_max_cdf_inv(double F, double mu, double beta) {
	ASSERT(F >= 0.0 && F <= 1.0 && beta > 0.0);
	return mu - beta * log(-log(F));
}

inline double gumbel_max_sf_inv(double S, double mu, double beta) {
	ASSERT(S >= 0.0 && S <= 1.0 && beta > 0.0);
	return mu - beta * log(-log1p(-S));
}

inline double gumbel_max_quantile(double F, double mu, double beta) {
	return gumbel_max_cdf_inv(F, mu, beta);
}

inline void gumbel_max_ms(double mean, double stddev, double& mu, double& beta) {
	ASSERT(stddev > 0.0);
	beta = sqrt(6.0 * sq(stddev) / sq(M_PI));
	mu = mean - beta * M_EM;
}

inline void gumbel_max_mv(double mean, double coeff_of_var, double& mu, double& beta) {
	ASSERT(coeff_of_var > 0.0);
	gumbel_max_ms(mean, abs(coeff_of_var * mean), mu, beta);
}

inline double gumbel_max_mean(double mu, double beta) {
	ASSERT(beta > 0.0);
	return mu + beta * M_EM;
}

inline double gumbel_max_variance(double mu, double beta) {
	ASSERT(beta > 0.0);
	return sq(M_PI) / 6.0 * sq(beta);
}

inline double gumbel_max_stddev(double mu, double beta) {
	return sqrt(gumbel_max_variance(mu, beta));
}

inline double gumbel_max_cov(double mu, double beta) {
	return gumbel_max_stddev(mu, beta) / abs(gumbel_max_mean(mu, beta));
}

// Lognormal distribution
// https://en.wikipedia.org/wiki/Lognormal_distribution

inline double lognormal_pdf(double x, double mu, double sigma) {
	ASSERT(x == x);
	return M_1_SQRT2PI * exp(-sq(log(x) - mu) / (2.0 * sq(sigma))) / (x * sigma);
}

inline double lognormal_cdf(double x, double mu, double sigma) {
	ASSERT(x == x);
	return 0.5 * erfc(-M_SQRT1_2 * (log(x) - mu) / sigma);
}

inline double lognormal_sf(double x, double mu, double sigma) {
	ASSERT(x == x);
	return 0.5 * erfc(M_SQRT1_2 * (log(x) - mu) / sigma);
}

inline double lognormal_cdf_inv(double F, double mu, double sigma) {
	ASSERT(F >= 0.0 && F <= 1.0);
	return exp(mu + M_SQRT2 * sigma * erf_inv(2.0 * F - 1.0));
}

inline double lognormal_sf_inv(double S, double mu, double sigma) {
	ASSERT(S >= 0.0 && S <= 1.0);
	return exp(mu + M_SQRT2 * sigma * erf_inv(1.0 - 2.0 * S));
}

inline double lognormal_quantile(double F, double mu, double sigma) {
	return lognormal_cdf_inv(F, mu, sigma);
}

inline double lognormal_log_pdf(double x, double mu, double sigma) {
	return -log(x * sigma) - 0.5 * M_LOG_2PI - 0.5 * sq((log(x) - mu) / sigma);
}

inline void lognormal_ms(double mean, double stddev, double& mu, double& sigma) {
	ASSERT(mean >= 0.0);
	ASSERT(stddev >= 0.0);
	double x = sqrt(sq(mean) + sq(stddev));
	mu = log(sq(mean) / x);
	sigma = sqrt(2.0 * log(x / mean));
}

inline void lognormal_mv(double mean, double coeff_of_var, double& mu, double& sigma) {
	ASSERT(coeff_of_var >= 0.0); lognormal_ms(mean, abs(coeff_of_var * mean), mu, sigma);
}

inline double lognormal_mean(double mu, double sigma) {
	ASSERT(sigma >= 0.0);
	return exp(mu + 0.5 * sq(sigma));
}

inline double lognormal_variance(double mu, double sigma) {
	ASSERT(sigma >= 0.0);
	return (exp(sq(sigma)) - 1.0) * exp(2.0 * mu + sq(sigma));
}

inline double lognormal_stddev(double mu, double sigma) {
	return sqrt(lognormal_variance(mu, sigma));
}

inline double lognormal_cov(double mu, double sigma) {
	ASSERT(sigma >= 0.0);
	return sqrt(exp(sq(sigma)) - 1.0);
}

// Normal distribution
// https://en.wikipedia.org/wiki/Normal_distribution

inline double normal_pdf(double x, double mu = 0.0, double sigma = 1.0) {
	ASSERT(x == x);
	return M_1_SQRT2PI * exp(-0.5 * sq((x - mu) / sigma)) / sigma;
}

inline double normal_cdf(double x, double mu = 0.0, double sigma = 1.0) {
	ASSERT(x == x);
	return 0.5 * erfc(-M_SQRT1_2 * (x - mu) / sigma);
}

inline double normal_sf(double x, double mu = 0.0, double sigma = 1.0) {
	ASSERT(x == x);
	return 0.5 * erfc(M_SQRT1_2 * (x - mu) / sigma);
}

inline double normal_cdf_inv(double F, double mu = 0.0, double sigma = 1.0) {
	ASSERT(F >= 0.0 && F <= 1.0);
	return mu + M_SQRT2 * sigma * erf_inv(2.0 * F - 1.0);
}

inline double normal_sf_inv(double S, double mu = 0.0, double sigma = 1.0) {
	ASSERT(S >= 0.0 && S <= 1.0);
	return mu + M_SQRT2 * sigma * erf_inv(1.0 - 2.0 * S);
}

inline double normal_quantile(double F, double mu = 0.0, double sigma = 1.0) {
	return normal_cdf_inv(F, mu, sigma);
}

inline double normal_log_pdf(double x, double mu, double sigma) {
	return -0.5 * M_LOG_2PI - log(sigma) - 0.5 * sq((x - mu) / sigma);
}

inline void normal_ms(double mean, double stddev, double& mu, double& sigma) {
	ASSERT(stddev > 0.0); mu = mean; sigma = stddev;
}

inline void normal_mv(double mean, double coeff_of_var, double& mu, double& sigma) {
	ASSERT(coeff_of_var > 0.0); normal_ms(mean, abs(coeff_of_var * mean), mu, sigma);
}

inline double normal_mean(double mu, double sigma) {
	ASSERT(sigma >= 0.0);
	return mu;
}

inline double normal_variance(double mu, double sigma) {
	ASSERT(sigma >= 0.0);
	return sq(sigma);
}

inline double normal_stddev(double mu, double sigma) {
	ASSERT(sigma >= 0.0);
	return sigma;
}
	
inline double normal_cov(double mu, double sigma) {
	return normal_stddev(mu, sigma) / normal_mean(mu, sigma);
}

// Student's t-distribution
// https://en.wikipedia.org/wiki/Student%27s_t-distribution
// #include <boost/math/distributions/students_t.hpp>

#ifdef BOOST_STATS_STUDENTS_T_HPP
inline double students_t_pdf(double x, double nu, double mu = 0.0, double s = 1.0) {
	ASSERT(x == x);
	ASSERT(nu > 0.0 && s > 0.0);
	boost::math::students_t_distribution dist(nu);
	return pdf(dist, (x - mu) / s) / s;
}

inline double students_t_cdf(double x, double nu, double mu = 0.0, double s = 1.0) {
	ASSERT(x == x);
	ASSERT(nu > 0.0 && s > 0.0);
	boost::math::students_t_distribution dist(nu);
	return cdf(dist, (x - mu) / s);
}

inline double students_t_sf(double x, double nu, double mu = 0.0, double s = 1.0) {
	return students_t_cdf(-(x - mu) / s, nu);
}

inline double students_t_cdf_inv(double F, double nu, double mu = 0.0, double s = 1.0) {
	ASSERT(F >= 0.0 && F <= 1.0);
	ASSERT(nu > 0.0 && s > 0.0);
	if (F == 0.0) return -M_INF;
	if (F == 1.0) return M_INF;
	boost::math::students_t_distribution dist(nu);
	return mu + s * quantile(dist, F);
}

inline double students_t_quantile(double F, double nu, double mu = 0.0, double s = 1.0) {
	return students_t_cdf_inv(F, nu, mu, s);
}

inline double students_t_sf_inv(double S, double nu, double mu = 0.0, double s = 1.0) {
	return mu - s * students_t_cdf_inv(S, nu);
}
#else
inline double students_t_pdf(double x, double nu, double mu = 0.0, double s = 1.0) {
	ASSERT(x == x);
	ASSERT(nu > 0.0 && s > 0.0);
	double c = 0.5 * (nu + 1.0);
	double z = (x - mu) / s;
	return tgamma(c) / (s * sqrt(nu * M_PI) * tgamma(0.5 * nu)) * pow(1.0 + sq(z) / nu, -c);
}
#endif

inline double students_t_mean(double nu, double mu = 0.0, double s = 1.0) {
	ASSERT(nu > 0.0 && s > 0.0);
	if (nu > 1.0)
		return mu;
	return M_NAN;
}

inline double students_t_variance(double nu, double mu = 0.0, double s = 1.0) {
	ASSERT(nu > 0.0 && s > 0.0);
	if (nu > 2.0)
		return sq(s) * nu / (nu - 2.0);
	else if (nu > 1.0)
		return M_INF;
	return M_NAN;
}

inline double students_t_stddev(double nu, double mu = 0.0, double s = 1.0) {
	return sqrt(students_t_variance(nu, mu, s));
}

inline double students_t_cov(double nu, double mu, double s) {
	return students_t_stddev(nu, mu, s) / students_t_mean(nu, mu, s);
}

// Triangular distribution
// https://en.wikipedia.org/wiki/Triangular_distribution

inline double triangular_pdf(double x, double a, double b, double c) {
	ASSERT(x >= a && x <= b);
	ASSERT(a < b && c >= a && c <= b);
	return x < c ? 2.0 * (x - a) / ((b - a) * (c - a)) :
			       2.0 * (b - x) / ((b - a) * (b - c));
}

inline double triangular_cdf(double x, double a, double b, double c) {
	ASSERT(x >= a && x <= b);
	ASSERT(a < b && c >= a && c <= b);
	return x < c ? sq(x - a) / ((b - a) * (c - a)) :
	               1.0 - sq(b - x) / ((b - a) * (b - c));
}

inline double triangular_sf(double x, double a, double b, double c) {
	return 1.0 - triangular_cdf(x, a, b, c);
}

inline double triangular_cdf_inv(double F, double a, double b, double c) {
	ASSERT(F >= 0.0 && F <= 1.0);
	ASSERT(a < b && c >= a && c <= b);
	double F_c = (c - a) / (b - a);
	return F < F_c ? a + sqrt(F * (b - a) * (c - a)) :
	                 b - sqrt((1.0 - F) * (b - a) * (b - c));
}

inline double triangular_sf_inv(double S, double a, double b, double c) {
	return triangular_cdf_inv(1.0 - S, a, b, c);
}

inline double triangular_quantile(double F, double a, double b, double c) {
	return triangular_cdf_inv(F, a, b, c);
}

inline double triangular_mean(double a, double b, double c) {
	ASSERT(a < b && c >= a && c <= b);
	return (a + b + c) / 3.0;
}

inline double triangular_variance(double a, double b, double c) {
	ASSERT(a < b && c >= a && c <= b);
	return (a*a + b*b + c*c - a*b - a*c - b*c) / 18.0;
}

inline double triangular_stddev(double a, double b, double c) {
	return sqrt(triangular_variance(a, b, c));
}

inline double triangular_cov(double a, double b, double c) {
	return triangular_stddev(a, b, c) / triangular_mean(a, b, c);
}

// Uniform distribution
// https://en.wikipedia.org/wiki/Uniform_distribution_(continuous)

inline double uniform_pdf(double x, double a, double b) {
	ASSERT(x >= a && x <= b);
	return 1.0 / (b - a);
}

inline double uniform_cdf(double x, double a, double b) {
	ASSERT(x >= a && x <= b);
	return (x - a) / (b - a);
}

inline double uniform_sf(double x, double a, double b) {
	ASSERT(x >= a && x <= b);
	return (b - x) / (b - a);
}

inline double uniform_cdf_inv(double F, double a, double b) {
	ASSERT(F >= 0.0 && F <= 1.0 && a < b);
	return a + F * (b - a);
}

inline double uniform_sf_inv(double S, double a, double b) {
	ASSERT(S >= 0.0 && S <= 1.0 && a < b);
	return b - S * (b - a);
}

inline double uniform_quantile(double F, double a, double b) {
	return uniform_cdf_inv(F, a, b);
}

inline void uniform_ms(double mean, double stddev, double& a, double& b) {
	ASSERT(stddev > 0.0);
	a = mean - M_SQRT3 * stddev;
	b = mean + M_SQRT3 * stddev;
}

inline void uniform_mv(double mean, double coeff_of_var, double& a, double& b) {
	ASSERT(coeff_of_var > 0.0);
	uniform_ms(mean, abs(coeff_of_var * mean), a, b);
}

inline double uniform_mean(double a, double b) {
	ASSERT(a < b);
	return 0.5 * (a + b);
}

inline double uniform_variance(double a, double b) {
	ASSERT(a < b);
	return sq(b - a) / 12.0;
}

inline double uniform_stddev(double a, double b) {
	return sqrt(uniform_variance(a, b));
}

inline double uniform_cov(double a, double b) {
	return uniform_stddev(a, b) / uniform_mean(a, b);
}

}

#endif  // __RDV_DISTRIBUTIONS__
