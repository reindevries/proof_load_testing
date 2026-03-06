// reliability.hpp
//
// Implements the first-order reliability method (FORM).
//
// Author: Rein de Vries
// Date: 7 December 2023

#ifndef __RDV_RELIABILITY__
#define __RDV_RELIABILITY__

#include "general.hpp"
#include "distributions.hpp"
#include <cstring>  // for memset

namespace rdv {

template <typename T>
inline size_t FORM(
		__in size_t n,                       // number of random variables
		__in T z_func,                       // limit state function
		__out double& beta,                  // reliability index
		__out_opt double* alpha = nullptr,   // internal storage of influence coefficients
		__inout_opt double* u = nullptr,     // input value is used to set starting point
		__out_opt double* dzdu = nullptr,    // internal storage of derivatives
		__out_opt double* u_prev = nullptr,  // internal storage of previous point
		__in_opt size_t max_iter = 100,      // maximum number of iterations
		__in_opt double du = 1e-3,           // use 0.1 or larger with FEM calculations
		__in_opt double beta_tol = 1e-6,     // use 1e-3 as a practical value
		__in_opt double max_step = 1.0) {    // maximum step size between iterations
	// Constants for relaxation
	constexpr size_t i_relax = 8;    // start relaxation after this many iterations
	constexpr double c_relax = 0.6;  // relaxation coefficient between 0 and 1
	
	// Allocate memory if needed
	size_t length = n * sizeof(double);
	bool free_u = !u;
	bool free_alpha = !alpha;
	bool free_dzdu = !dzdu;
	bool free_u_prev = !u_prev;
	if (free_u) u = (double*)malloc(length);
	if (free_alpha) alpha = (double*)malloc(length);
	if (free_dzdu) dzdu = (double*)malloc(length);
	if (free_u_prev) u_prev = (double*)malloc(length);
	
	// Initialize vectors
	if (free_u) memset(u, 0, length);
	memset(alpha, 0, length);
	memset(dzdu, 0, length);
	memcpy(u_prev, u, length);
	
	// Start iterations
	double half_du = 0.5 * du;
	double beta_prev = 0.0;
	
	for (size_t iteration = 1; iteration <= max_iter; ++iteration) {
		// Calculate value in current point
		double z = z_func((const double*)u);
		
		// Calculate derivatives in current point
		for (size_t i = 0; i != n; ++i) {
			// Save previous derivatives
			double dzdu_prev = dzdu[i];
			
			// One-sided finite difference
			//du = copysign(du, -u[i]);  // towards origin, failure region may be noisy
			//
			//u[i] += du;
			//dzdu[i] = (z_func((const double*)u) - z) / du;
			//u[i] -= du;
			
			// Central finite difference
			u[i] += half_du;
			double z_1 = z_func((const double*)u);
			u[i] -= half_du;
			u[i] -= half_du;
			double z_2 = z_func((const double*)u);
			u[i] += half_du;
			
			dzdu[i] = (z_1 - z_2) / du;
			
			// Counteract oscillations with relaxation
			if (iteration > i_relax)
				dzdu[i] = (1.0 - c_relax) * dzdu[i] + c_relax * dzdu_prev;
		}
		
		// Calculate mean and stddev of linearized limit state function
		//
		// Random variables are standard normally distributed
		// and linearized in the point (U_1*, U_2*, ...)
		// Z = Z(U_1, U_2, ...) + dZ/dU_1 * (U_1 - U_1*) + dZ/dU_2 * (U_2 - U_2*) + ...
		double mu_Z = z;
		for (size_t i = 0; i != n; ++i)
			mu_Z -= dzdu[i] * u[i];
		
		double sigma_Z = 0.0;
		for (size_t i = 0; i != n; ++i)
			sigma_Z += sq(dzdu[i]);
		sigma_Z = sqrt(sigma_Z);
		
		// Calculate reliability index and influence coefficients
		beta = mu_Z / sigma_Z;  // = -U_Z
		for (size_t i = 0; i != n; ++i)
			alpha[i] = dzdu[i] / sigma_Z;
		
		// Calculate new point
		for (size_t i = 0; i != n; ++i)
			u[i] = -alpha[i] * beta;
		
		// Calculate step size
		double step = 0.0;
		for (size_t i = 0; i != n; ++i)
			step += sq(u[i] - u_prev[i]);
		step = sqrt(step);
		
		// Check if converged
		if (iteration > 1 && step < max_step && abs(beta - beta_prev) < beta_tol) {
			if (free_u) free(u);
			if (free_alpha) free(alpha);
			if (free_dzdu) free(dzdu);
			if (free_u_prev) free(u_prev);
			return iteration;
		}
		
		// Limit step size
		if (step > max_step) {
			double c = max_step / step;
			for (size_t i = 0; i != n; ++i)
				u[i] = u_prev[i] + c * (u[i] - u_prev[i]);
		}
		
		// Counteract oscillations with relaxation
		if (iteration > i_relax) {
			for (size_t i = 0; i != n; ++i)
				u[i] = (1.0 - c_relax) * u[i] + c_relax * u_prev[i];
		}
		
		// Save result for next iteration
		beta_prev = beta;
		memcpy(u_prev, u, length);
	}
	
	// We end up here in case of no convergence
	if (free_u) free(u);
	if (free_alpha) free(alpha);
	if (free_dzdu) free(dzdu);
	if (free_u_prev) free(u_prev);
	
	ASSERT(false);
	return 0;
}

}

#endif  // __RDV_RELIABILITY__
