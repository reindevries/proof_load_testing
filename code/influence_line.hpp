// influence_line.hpp
//
// This header file contains influence line functions for single-span and two-span bridges. They
// have been derived analytically based on the fourth-order beam differential equation.
//
// Author: Rein de Vries <rein.devries@tno.nl>
// Date: 10 November 2025

#ifndef __INFLUENCE_LINE__
#define __INFLUENCE_LINE__

#include "rdv/general.hpp"

inline double single_span_M(double x_M, double x_P, double L, double P){
	ASSERT(x_M >= 0.0 && x_M <= L);  // point of interest not on bridge
	if (x_P < 0.0 || x_P > L)
		return 0.0;
	if (x_M < x_P)
		return P * x_M * (L - x_P) / L;
	return P * x_P * (L - x_M) / L;
}

inline double single_span_V(double x_V, double x_P, double L, double P) {
	ASSERT(x_V >= 0.0 && x_V <= L);  // point of interest not on bridge
	if (x_P < 0.0 || x_P > L)
		return 0.0;
	if (x_V < x_P)
		return P * (L - x_P) / L;
	return -P * x_P / L;
}

inline double two_spans_M(double x_M, double x_P, double L_1, double L_2, double P) {
	double L = L_1 + L_2;
	ASSERT(x_M >= 0.0 && x_M <= L);  // point of interest not on bridge
	if (x_P < 0.0 || x_P > L)
		return 0.0;
	if (x_P > L_1)
		return two_spans_M(L - x_M, L - x_P, L_2, L_1, P);  // mirror
	if (x_M < x_P)
		return 0.5 * P * (2.0 * L_1*L_1*L_1 + 2.0 * L_1*L_1 * L_2 - 3.0 * L_1*L_1 * x_P - 
				2.0 * L_1 * L_2 * x_P + x_P*x_P*x_P) * x_M / (L_1*L_1 * L);
	if (x_M < L_1)
		return 0.5 * P * (2.0 * L_1*L_1*L_1 + 2.0 * L_1*L_1 * L_2 - 3.0 * L_1*L_1 * x_M -
				2.0 * L_1 * L_2 * x_M + x_M * x_P*x_P) * x_P / (L_1*L_1 * L);
	return -0.5 * P * x_P * (L_1*L_1 - x_P*x_P) * (L - x_M) / (L_1 * L_2 * L);
}

inline double two_spans_V(double x_V, double x_P, double L_1, double L_2, double P) {
	double L = L_1 + L_2;
	ASSERT(x_V >= 0.0 && x_V <= L);  // point of interest not on bridge
	if (x_P < 0.0 || x_P > L)
		return 0.0;
	if (x_P > L_1)
		return -two_spans_V(L - x_V, L - x_P, L_2, L_1, P);  // mirror
	if (x_V < x_P)
		return 0.5 * P * (2.0 * L_1*L_1*L_1 + 2.0 * L_1*L_1 * L_2 - 3.0 * L_1*L_1 * x_P - 
				2.0 * L_1 * L_2 * x_P + x_P*x_P*x_P) / (L_1*L_1 * L);
	if (x_V < L_1)
		return -0.5 * P * x_P * (3.0 * L_1*L_1 + 2.0 * L_1 * L_2 - x_P*x_P) / (L_1*L_1 * L);
	return 0.5 * P * x_P * (L_1*L_1 - x_P*x_P) / (L_1 * L_2 * L);
}

#endif  // __INFLUENCE_LINE__
