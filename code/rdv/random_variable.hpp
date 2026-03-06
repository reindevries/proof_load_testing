// random_variable.hpp
// 
// In this file the abstract base class random_variable is defined along with macros to create
// derived random variable classes with 2 or 3 parameters.
// 
// Author: Rein de Vries
// Date: 21 November 2023

#ifndef __RDV_RANDOM_VARIABLE__
#define __RDV_RANDOM_VARIABLE__

#include "general.hpp"
#include "distributions.hpp"
#include "random.hpp"

namespace rdv
{

// Class random_variable

class random_variable
{
public:
	random_variable()
		{}
	
	~random_variable()
		{}
	
	virtual double pdf(double x) const = 0;
	virtual double cdf(double x) const = 0;
	virtual double sf(double x) const = 0;
	virtual double cdf_inv(double F) const = 0;
	virtual double sf_inv(double S) const = 0;
	virtual double quantile(double F) const = 0;
	virtual double mean() const = 0;
	virtual double stddev() const = 0;
	virtual double cov() const = 0;
	
	template <class G>
	double random(G& gen) const
		{return cdf_inv(rdv::random_uniform(gen));}
	
	double random() const  // not thread-safe
	{
		static std::mt19937 gen(std::random_device{}());
		return random(gen);
	}
	
	template <class G>
	double operator ()(G& gen) const
		{return random(gen);}
	
	virtual double from_std_norm(double u) const
	{
		if (u < 0.0)
			return cdf_inv(normal_cdf(u));
		return sf_inv(normal_sf(u));
	}
	
	virtual double to_std_norm(double x) const
	{
		if (x < mean())
			return normal_cdf_inv(cdf(x));
		return normal_sf_inv(sf(x));
	}
	
	double interval(double x_1, double x_2) const
	{
		if (x_1 < mean())
			return cdf(x_2) - cdf(x_1);
		return sf(x_1) - sf(x_2);
	}
};

// Macro to create random variable with 2 parameters

#define RANDOM_VARIABLE_2(name, param1, param2, ...)                                  \
                                                                                      \
class rv_##name : public random_variable                                              \
{                                                                                     \
public:                                                                               \
	rv_##name() : m_##param1(M_NAN), m_##param2(M_NAN)                                \
		{}                                                                            \
	                                                                                  \
	rv_##name(double param1, double param2) : m_##param1(param1), m_##param2(param2)  \
		{}                                                                            \
	                                                                                  \
	static rv_##name ms(double mean, double stddev)                                   \
	{                                                                                 \
		double param1, param2;                                                        \
		name##_ms(mean, stddev, param1, param2);                                      \
		return rv_##name(param1, param2);                                             \
	}                                                                                 \
	                                                                                  \
	static rv_##name mv(double mean, double coeff_of_var)                             \
	{                                                                                 \
		double param1, param2;                                                        \
		name##_mv(mean, coeff_of_var, param1, param2);                                \
		return rv_##name(param1, param2);                                             \
	}                                                                                 \
	                                                                                  \
	static rv_##name kv(double F, double char_value, double coeff_of_var)             \
	{                                                                                 \
		double param1, param2;                                                        \
		name##_mv(1.0, coeff_of_var, param1, param2);                                 \
		double mean = char_value / name##_cdf_inv(F, param1, param2);                 \
		name##_mv(mean, coeff_of_var, param1, param2);                                \
		return rv_##name(param1, param2);                                             \
	}                                                                                 \
	                                                                                  \
	void set_ms(double mean, double stddev)                                           \
		{name##_ms(mean, stddev, m_##param1, m_##param2);}                            \
	                                                                                  \
	void set_mv(double mean, double coeff_of_var)                                     \
		{name##_mv(mean, coeff_of_var, m_##param1, m_##param2);}                      \
	                                                                                  \
	double pdf(double x) const override                                               \
		{return name##_pdf(x, m_##param1, m_##param2);}                               \
	                                                                                  \
	double cdf(double x) const override                                               \
		{return name##_cdf(x, m_##param1, m_##param2);}                               \
	                                                                                  \
	double sf(double x) const override                                                \
		{return name##_sf(x, m_##param1, m_##param2);}                                \
	                                                                                  \
	double cdf_inv(double F) const override                                           \
		{return name##_cdf_inv(F, m_##param1, m_##param2);}                           \
	                                                                                  \
	double sf_inv(double S) const override                                            \
		{return name##_sf_inv(S, m_##param1, m_##param2);}                            \
	                                                                                  \
	double quantile(double F) const override                                          \
		{return name##_quantile(F, m_##param1, m_##param2);}                          \
	                                                                                  \
	double mean() const override                                                      \
		{return name##_mean(m_##param1, m_##param2);}                                 \
	                                                                                  \
	double stddev() const override                                                    \
		{return name##_stddev(m_##param1, m_##param2);}                               \
	                                                                                  \
	double cov() const override                                                       \
		{return name##_cov(m_##param1, m_##param2);}                                  \
	                                                                                  \
	double m_##param1, m_##param2;                                                    \
	                                                                                  \
	__VA_ARGS__                                                                       \
};

// Macro to create random variable with 3 parameters

#define RANDOM_VARIABLE_3(name, param1, param2, param3, ...)               \
                                                                           \
class rv_##name : public random_variable                                   \
{                                                                          \
public:                                                                    \
	rv_##name() : m_##param1(M_NAN), m_##param2(M_NAN), m_##param3(M_NAN)  \
		{}                                                                 \
	                                                                       \
	rv_##name(double param1, double param2, double param3) :               \
			m_##param1(param1), m_##param2(param2), m_##param3(param3)     \
		{}                                                                 \
	                                                                       \
	double pdf(double x) const override                                    \
		{return name##_pdf(x, m_##param1, m_##param2, m_##param3);}        \
	                                                                       \
	double cdf(double x) const override                                    \
		{return name##_cdf(x, m_##param1, m_##param2, m_##param3);}        \
	                                                                       \
	double sf(double x) const override                                     \
		{return name##_sf(x, m_##param1, m_##param2, m_##param3);}         \
	                                                                       \
	double cdf_inv(double F) const override                                \
		{return name##_cdf_inv(F, m_##param1, m_##param2, m_##param3);}    \
	                                                                       \
	double sf_inv(double S) const override                                 \
		{return name##_sf_inv(S, m_##param1, m_##param2, m_##param3);}     \
	                                                                       \
	double quantile(double F) const override                               \
		{return name##_quantile(F, m_##param1, m_##param2, m_##param3);}   \
	                                                                       \
	double mean() const override                                           \
		{return name##_mean(m_##param1, m_##param2, m_##param3);}          \
	                                                                       \
	double stddev() const override                                         \
		{return name##_stddev(m_##param1, m_##param2, m_##param3);}        \
	                                                                       \
	double cov() const override                                            \
		{return name##_cov(m_##param1, m_##param2, m_##param3);}           \
	                                                                       \
	double m_##param1, m_##param2, m_##param3;                             \
	                                                                       \
	__VA_ARGS__                                                            \
};

// Random variable list

#ifdef BOOST_MATH_SF_GAMMA_HPP
RANDOM_VARIABLE_2(gamma, k, theta)
#endif
RANDOM_VARIABLE_2(gumbel_max, mu, beta,
	rv_gumbel_max shifted(double factor) const
		{return rv_gumbel_max(m_mu + m_beta * ln(factor), m_beta);}
)
RANDOM_VARIABLE_2(lognormal, mu, sigma,
	double from_std_norm(double u) const override
		{return exp(m_mu + m_sigma * u);}
	double to_std_norm(double x) const override
		{return (ln(x) - m_mu) / m_sigma;}
)
RANDOM_VARIABLE_2(normal, mu, sigma, 
	double from_std_norm(double u) const override
		{return m_mu + m_sigma * u;}
	double to_std_norm(double x) const override
		{return (x - m_mu) / m_sigma;}
)
RANDOM_VARIABLE_2(uniform, a, b)
#ifdef BOOST_STATS_STUDENTS_T_HPP
RANDOM_VARIABLE_3(students_t, nu, mu, s)
#endif
RANDOM_VARIABLE_3(triangular, a, b, c)

}

#endif  // __RDV_RANDOM_VARIABLE__
