// load_effect.hpp
//
// This header file contains functions to calculate the load effect on single-span bridges given
// customary load models LM1, HS20 and HL93.
//
// Author: Rein de Vries <rein.devries@tno.nl>
// Date: 10 November 2025

#ifndef __LOAD_EFFECT__
#define __LOAD_EFFECT__

#include <vector>
#include "influence_line.hpp"

inline void load_effect(double q, double L, const std::vector<double>& axle_x, 
		const std::vector<double>& axle_P, double& M, double& V) {
	// Load effect locations
	const double x_M = L / 2.0;
	const double x_V = 1.0;
	
	// Load effects from distributed load
	double M_dist = 1.0/8.0 * q * L * L;
	double V_dist = q * (L / 2.0 - x_V);
	
	// Load effects from axles
	double M_axles = 0.0;
	double V_axles = 0.0;
	
	for (double x = -axle_x.back(); x < L; x += 0.01) {
		double M = 0.0;
		double V = 0.0;
		for (size_t i = 0; i < size(axle_x); i++) {
			M += single_span_M(x_M, x + axle_x[i], L, axle_P[i]);
			V += single_span_V(x_V, x + axle_x[i], L, axle_P[i]);
		}
		M_axles = std::max(M, M_axles);
		V_axles = std::max(V, V_axles);
	}
	
	// Result via superposition
	M = M_dist + M_axles;
	V = V_dist + V_axles;
}

inline void load_effect(double q, const std::vector<double>& L, const std::vector<double>& axle_x,
		const std::vector<double>& axle_P, std::vector<double>& M_Q, std::vector<double>& V_Q) {
	// Prepare result vectors
	M_Q.resize(size(L));
	V_Q.resize(size(L));
	
	// Calculate load effect for every length
	for (size_t i = 0; i < size(L); i++)
		load_effect(q, L[i], axle_x, axle_P, M_Q[i], V_Q[i]);
}

template <typename T>
inline void load_effect_LM1(const T& L, T& M_Q, T& V_Q) {
	load_effect(3.0 * 9e3, L, {0.0, 1.2}, {300e3, 300e3}, M_Q, V_Q);
}

template <typename T>
inline void load_effect_HS20(const T& L, T& M_Q, T& V_Q) {
	load_effect(0.0, L, {0.0, 4.27, 8.54}, {142.3e3, 142.3e3, 35.6e3}, M_Q, V_Q);
}

template <typename T>
inline void load_effect_HL93(const T& L, T& M_Q, T& V_Q) {
	// Truck and distributed load
	std::vector<double> M_Q_truck, V_Q_truck;
	load_effect(9.34e3, L, {0.0, 4.27, 8.54}, {142.3e3, 142.3e3, 35.6e3}, M_Q_truck, V_Q_truck);
	
	// Tandem and distributed load
	std::vector<double> M_Q_tandem, V_Q_tandem;
	load_effect(9.34e3, L, {0.0, 1.22}, {111.2e3, 111.2e3}, M_Q_tandem, V_Q_tandem);
	
	// Maximum of both
	M_Q.resize(size(L));
	V_Q.resize(size(L));
	for (size_t i = 0; i < size(L); i++) {
		M_Q[i] = std::max(M_Q_truck[i], M_Q_tandem[i]);
		V_Q[i] = std::max(V_Q_truck[i], V_Q_tandem[i]);
	}
}

#endif  // __LOAD_EFFECT__
