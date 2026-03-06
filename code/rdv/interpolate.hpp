// interpolate.hpp
//
// This file provides an interpolation routine that operates on two vectors. One vector contains the
// x-values and the other vector contains the corresponding y-values. Given a query x-value, the 
// routine returns the interpolated y-value. The routine supports different extrapolation methods
// for query x-values that fall outside the range of the provided x-values.
//
// Author: Rein de Vries
// Date: 3 September 2024

#ifndef __RDV_INTERPOLATE__
#define __RDV_INTERPOLATE__

#include "general.hpp"

namespace rdv {

enum class extrapolation {
	none,
	constant,
	linear
};

template<typename V, typename T>
inline T interpolate(const V& x, const V& y, T xq, extrapolation extrapol = extrapolation::none) {
	// Check vector size
	size_t n = size(x);
	if (n < 1 || size(y) != n)
		{ASSERT(false); return std::numeric_limits<T>::quiet_NaN();}
	
	// Handle endpoints
	if (xq == x[0])
		return y[0];
	
	if (xq == x[n-1])
		return y[n-1];
	
	// Extrapolate left
	if (xq < x[0]) {
		if (extrapol == extrapolation::constant)
			return y[0];
		else if (extrapol == extrapolation::linear && n >= 2) {
			T dx = x[1] - x[0];
			T dy = y[1] - y[0];
			ASSERT(dx > 0.0);
			return y[0] + dy / dx * (xq - x[0]);
		}
		
		ASSERT(false); 
		return std::numeric_limits<T>::quiet_NaN();
	}
	
	// Extrapolate right
	if (xq > x[n-1]) {
		if (extrapol == extrapolation::constant)
			return y[n-1];
		else if (extrapol == extrapolation::linear && n >= 2) {
			T dx = x[n-1] - x[n-2];
			T dy = y[n-1] - y[n-2];
			ASSERT(dx > 0.0);
			return y[n-1] + dy / dx * (xq - x[n-1]);
		}
		
		ASSERT(false); 
		return std::numeric_limits<T>::quiet_NaN();
	}
	
	// Find the index of the first element with the same or a larger x-value
	ASSERT(n >= 2);
	size_t i = std::distance(begin(x), std::lower_bound(begin(x), end(x), xq));
	ASSERT(i > 0 && i < n);
	
	// Perform the lineair interpolation
	T dx = x[i] - x[i-1];
	ASSERT(dx > 0.0);
	double f = double(xq - x[i-1]) / dx;
	return (1.0 - f) * y[i-1] + f * y[i];
}

}

#endif  // __RDV_INTERPOLATE__
