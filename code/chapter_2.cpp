// chapter_2.cpp
//
// This file contains the Monte Carlo simulation procedure to produce graphs displaying the 
// evolution of reliability over time. Importance sampling is used on variables f_y, C_0Q and θ_E to 
// increase the accuracy and ultimately reduce the run time. The output is displayed in Figures 2.2 
// and 2.3 of the dissertation.
//
// Author: Rein de Vries <rein.devries@tno.nl>
// Date: 10 November 2025

#include "rdv/msvcr.hpp"
#include <boost/math/special_functions/gamma.hpp>
#include "rdv/random_variable.hpp"
#include "rdv/vector_ext.hpp"
#include <omp.h>
#include "influence_line.hpp"

int main(int argc, char *argv[]) {
	// Initialisation
	rdv::tic();
	std::cout.setf(std::ios::left);
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
	
	// Constants
	constexpr double L = 10.0;
	constexpr double h = 0.8;
	constexpr double b = 3.0;
	constexpr double g = 9.81;
	
	// Maximum moment permanent load (G)
	double g_DL = 2450.0 * g * h * b;
	double M_DL = 1.0/8.0 * g_DL * L * L;
	double g_SDL = (1.2e3 + 1.5e3) * b;
	double M_SDL = 1.0/8.0 * g_SDL * L * L;
	double M_G = M_DL + M_SDL;
	
	std::cout << "M_DL = " << M_DL / 1e3 << " kNm\n";
	std::cout << "M_SDL = " << M_SDL / 1e3 << " kNm\n";
	std::cout << "M_G = " << M_G / 1e3 << " kNm\n";
	
	// Maximum moment variable load (Q)
	double q = 400.0 * g * b;
	double M_dist = 1.0/8.0 * q * L * L;
	std::vector<double> axle_x = {0.0,      1.5,      7.5};
	std::vector<double> axle_P = {20e3 * g, 20e3 * g, 20e3 * g};
	
	double M_tandem = 0.0;
	for (double x = -axle_x.back(); x < L; x += 0.01) {
		double M = 0.0;
		for (size_t i = 0; i < axle_x.size(); i++)
			M += single_span_M(L / 2.0, x + axle_x[i], L, axle_P[i]);
		M_tandem = std::max(M, M_tandem);
	}
	
	double S = 1.0 + 40.0 / (100.0 + L);
	double M_Q = (M_dist + M_tandem) * S;
	double M_Ed = M_G + M_Q;
	
	std::cout << "M_dist = " << M_dist / 1e3 << " kNm\n";
	std::cout << "M_tandem = " << M_tandem / 1e3 << " kNm\n";
	std::cout << "S = " << S << "\n";
	std::cout << "M_Q = " << M_Q / 1e3 << " kNm\n";
	std::cout << "M_Ed = M_G + M_Q = " << M_Ed / 1e3 << " kNm\n";
	
	// Geometry
	// A_s1 * d_1 + A_s2 * d_2 = A_s * d
	// d = (A_s1 * d_1 + A_s2 * d_2) / A_s
	double A_s1 = M_PI * rdv::sq(0.028 / 2.0) * b / 0.207 * 2.0;
	double A_s2 = M_PI * rdv::sq(0.022 / 2.0) * b / 0.207;
	double d_1 = h - 0.030 - 0.028 / 2.0;
	double d_2 = h - 0.030 - 0.028 - 0.028 - 0.022 / 2.0;
	double A_s = A_s1 + A_s2;
	double d = (A_s1 * d_1 + A_s2 * d_2) / A_s;
	double ρ = A_s / (b * d);
	
	std::cout << "A_s = " << A_s << " m^2\n";
	std::cout << "d = " << d << "\n";
	std::cout << "ρ = " << ρ * 100.0 << "%\n";
	
	// Design resistance
	double f_yd = 147e6;
	double f_cd = 5.9e6;
	double F_s = A_s * f_yd;
	double M_Rd = F_s * (d - F_s / (2.0 * 0.85 * f_cd * b));
	
	std::cout << "M_Rd = " << M_Rd / 1e3 << " kNm\n";
	std::cout << "UC = M_Ed / M_Rd = " << M_Ed / M_Rd << "\n";
	
	// Random variables
	auto rv_h = rdv::rv_normal::mv(h, 0.02);              // [m]
	auto rv_a = rdv::rv_gamma::mv(h - d, 0.17);           // [m]
	auto rv_f_c = rdv::rv_lognormal::mv(21.1e6, 0.20);    // [N/m^2] B250
	auto rv_f_y = rdv::rv_lognormal::mv(261e6, 0.05);     // [N/m^2] QR24
	auto rv_θ_R = rdv::rv_lognormal::mv(1.0, 0.05);       // [-]
	auto rv_G_DL = rdv::rv_normal::mv(M_DL, 0.05);        // [Nm]
	auto rv_G_SDL = rdv::rv_normal::mv(M_SDL, 0.10);      // [Nm]
	auto rv_Q = rdv::rv_gumbel_max::mv(1150e3, 0.025);    // [Nm]
	auto rv_C_0Q = rdv::rv_lognormal::mv(1.1, 0.10);      // [-]
	auto rv_θ_E = rdv::rv_lognormal::mv(1.0, 0.11);       // [-]
	auto rv_c_Q0 = rdv::rv_lognormal::mv(0.78, 0.10);     // [-]
	auto rv_Δ_c_Q = rdv::rv_lognormal::mv(0.004, 0.10);   // increase per year
	auto rv_c_R_t0 = rdv::rv_lognormal::mv(20.0, 0.10);   // [-]
	auto rv_Δ_c_R = rdv::rv_lognormal::mv(0.0025, 0.10);  // decrease per year
	auto rv_Q_PL_1 = rdv::rv_normal::mv(1800e3, 0.01);    // [Nm]
	auto rv_Q_PL_2 = rdv::rv_normal::mv(2000e3, 0.01);    // [Nm]
	
	std::cout << "m_a = " << rv_a.mean() << " m\n";
	std::cout << "m_f_c (60 years) = " << exp(0.25 * (1.0 - sqrt(28.0 / 
			(60.0 * 365.25)))) * rv_f_c.mean() / 1e6 << " MPa (not used)\n";
	std::cout << "f_ck = " << rv_f_c.quantile(0.05) / 1e6 << " MPa\n";
	std::cout << "f_yk = " << rv_f_y.quantile(0.05) / 1e6 << " MPa\n";
	std::cout << "Q_k = " << rv_Q.quantile(0.999) / 1e3 << " kNm\n";
	
	// Helper random variables for importance sampling
	auto rv_f_y_h = rv_f_y.mv(rv_f_y.mean() / 1.1, rv_f_y.cov() * 1.2);
	auto rv_C_0Q_h = rv_C_0Q.mv(rv_C_0Q.mean() * 1.1, rv_C_0Q.cov() * 1.2);
	auto rv_θ_E_h = rv_θ_E.mv(rv_θ_E.mean() * 1.1, rv_θ_E.cov() * 1.2);
	
	// Prepare random sampling
	std::vector<uint32_t> seeds = rdv::generate_seeds(omp_get_max_threads());
	size_t samples = size_t(1e8);
	size_t years = 101;
	std::vector<double> failures(years, 0.0);
	double total_w = 0.0;
	std::cout << "Performing Monte Carlo simulation...\n";
	
	// Start parallel threads
	#pragma omp parallel reduction(+: failures, total_w)
	{
		// Thread-local random number generator
		std::mt19937 gen(seeds[omp_get_thread_num()]);
		
		#pragma omp for
		for (size_t i_sample = 0; i_sample < samples; i_sample++) {
			// Sample time-independent variables
			double h = rv_h(gen);
			double a = rv_a(gen);
			double f_c = rv_f_c(gen);
			double f_y = rv_f_y_h(gen);
			double θ_R = rv_θ_R(gen);
			double G_DL = rv_G_DL(gen);
			double G_SDL = rv_G_SDL(gen);
			double C_0Q = rv_C_0Q_h(gen);
			double θ_E = rv_θ_E_h(gen);
			double c_Q0 = rv_c_Q0(gen);
			double Δ_c_Q = rv_Δ_c_Q(gen);
			double c_R_t0 = rv_c_R_t0(gen);
			double Δ_c_R = rv_Δ_c_R(gen);
			
			// Enforce physical constraints
			if (h < 0.0) h = 0.0;
			if (G_DL < 0.0) G_DL = 0.0;
			if (G_SDL < 0.0) G_SDL = 0.0;
			
			// Importance sampling weight
			double w = 1.0;
			w *= rv_f_y.pdf(f_y) / rv_f_y_h.pdf(f_y);
			w *= rv_C_0Q.pdf(C_0Q) / rv_C_0Q_h.pdf(C_0Q);
			w *= rv_θ_E.pdf(θ_E) / rv_θ_E_h.pdf(θ_E);
			total_w += w;
			
			// Calculate resistance
			double d = h - a;
			double F_s = A_s * f_y;
			double R = F_s * (d - F_s / (2.0 * 0.85 * f_c * b));
			
			// Loop over the years
			for (size_t i_year = 0; i_year < years; i_year++) {
				// Traffic load trend
				double c_Q = c_Q0 + Δ_c_Q * i_year;
				
				// Traffic load effect, including proof load effect
				double Q = c_Q * C_0Q * rv_Q(gen);
				if (i_year == 60)
					Q = std::max(Q, rv_Q_PL_1(gen));
				else if (i_year == 70)
					Q = std::max(Q, rv_Q_PL_2(gen));
				
				// Resistance deterioration
				double c_R;
				if (i_year < c_R_t0)
					c_R = 1.0;
				else
					c_R = std::max(1.0 - Δ_c_R * (i_year - c_R_t0), 0.0);
				
				// Evaluate limit state function
				double Z = θ_R * c_R * R - θ_E * (G_DL + G_SDL + Q);
				if (Z < 0.0) {
					failures[i_year] += w;
					break;
				}
			}
		}
	}
	
	// Report relative importance sampling weight
	std::cout << "total_w / samples = " << total_w / samples << "\n";
	
	// Print results
	std::cout << std::setw(10) << "Year" << "Reliability\n";
	for (size_t i_year = 0; i_year < years; i_year++) {
		double P_f_i = failures[i_year] / (samples - sum(head(failures, i_year)));
		std::cout << std::setw(10) << 1960 + i_year << -rdv::normal_quantile(P_f_i) << "\n";
	}
	
	// Print elapsed time
	rdv::toc();
	return 0;
}
