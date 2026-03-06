// golden_section.hpp
//
// This file contains an implementation of the Golden-section search to find the minimum of a
// one-dimensional function. The tolerance and maximum number of iterations can be specified.
//
// Author: Rein de Vries
// Date: 3 September 2024

#ifndef __RDV_GOLDEN_SECTION__
#define __RDV_GOLDEN_SECTION__

#include "general.hpp"

namespace rdv {

template <typename F>
inline size_t golden_section_minimize(F&& f, double x_min, 
		double x_max, double& x, double x_tol = 1e-6, size_t max_iter = 100) {
	// Set initial bounds
	double a = x_min;
	double b = x_max;
	
	for (size_t i = 0; i < max_iter; i++) {
		// Determine points to evaluate
		double x_1 = b - (b - a) / M_PHI;
		double x_2 = a + (b - a) / M_PHI;
		
		// Evaluate function to determine which bound to shift
		if (f(x_1) < f(x_2))
			b = x_2;
		else
			a = x_1;
		
		// Check for convergence
		if (b - a <= x_tol) {
			x = 0.5 * (a + b);
			return i + 1;
		}
	}
	
	// Failed to converge within the specified number of iterations
	ASSERT(false);
	x = M_NAN;
	return 0;
}

}

#endif  // __RDV_GOLDEN_SECTION__
