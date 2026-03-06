// chapter_4_time-dependent.cpp
//
// This code is quite similar to chapter_2.cpp as it also performs a time-dependent reliability
// analysis. However, in this case a low-informative prior distribution is assumed for the resistance
// with its mean value directly calculated from the mean values of the load variables. The output of
// the analysis is provided in Figure 4.4 of the dissertation.
//
// Author: Rein de Vries <rein.devries@tno.nl>
// Date: 10 November 2025

#include "rdv/msvcr.hpp"
#include "rdv/random_variable.hpp"
#include "rdv/vector_ext.hpp"
#include <omp.h>

int main(int argc, char *argv[]) {
	// Initialisation
	rdv::tic();
	std::cout.setf(std::ios::left);
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
	
	// Random variables
	auto rv_θ_R = rdv::rv_lognormal::mv(1.0, 0.05);
	auto rv_G_DL = rdv::rv_normal::mv(721e3, 0.05);
	auto rv_G_SDL = rdv::rv_normal::mv(101e3, 0.1);
	auto rv_Q = rdv::rv_gumbel_max::mv(1150e3, 0.025);
	auto rv_C_0Q = rdv::rv_lognormal::mv(1.1, 0.1);
	auto rv_θ_E = rdv::rv_lognormal::mv(1.0	, 0.11);
	auto rv_c_Q0 = rdv::rv_lognormal::mv(0.78, 0.1);
	auto rv_Δ_c_Q = rdv::rv_lognormal::mv(0.004, 0.1);
	auto rv_c_R_t0 = rdv::rv_lognormal::mv(20.0, 0.1);
	auto rv_Δ_c_R = rdv::rv_lognormal::mv(0.0025, 0.1);
	auto rv_Q_PL_1 = rdv::rv_normal::mv(1800e3, 0.01);  // now needs 2000
	auto rv_Q_PL_2 = rdv::rv_normal::mv(2000e3, 0.01);  // now needs 2100
	
	// Low-informative prior for the resistance
	double m_R = 1.5 * (rv_G_DL.mean() + rv_G_SDL.mean() + rv_Q.mean());
	auto rv_R = rdv::rv_lognormal::mv(m_R, 0.5);
	std::cout << "m_R = " << m_R / 1e3 << " kNm\n";
	std::cout << "s_R = " << rv_R.stddev() / 1e3 << " kNm\n";
	
	// Helper distributions for importance sampling
	auto rv_R_h = rv_R.mv(rv_R.mean() / 1.1, rv_R.cov() * 1.2);
	auto rv_C_0Q_h = rv_C_0Q.mv(rv_C_0Q.mean() * 1.1, rv_C_0Q.cov() * 1.2);
	auto rv_θ_E_h = rv_θ_E.mv(rv_θ_E.mean() * 1.1, rv_θ_E.cov() * 1.2);
	
	// Prepare random sampling
	std::vector<uint32_t> seeds = rdv::generate_seeds(omp_get_max_threads());
	size_t samples = size_t(2e8);
	size_t years = 101;
	std::vector<double> failures(years, 0.0);
	double total_w = 0.0;
	std::cout << "Performing Bayesian Monte Carlo simulation...\n";
	
	// Start parallel threads
	#pragma omp parallel reduction(+:failures, total_w)
	{
		// Thread-local random number generator
		std::mt19937 gen(seeds[omp_get_thread_num()]);
		
		#pragma omp for
		for (size_t i_sample = 0; i_sample < samples; i_sample++) {
			// Sample time-independent variables
			double R = rv_R_h(gen);
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
			if (G_DL < 0.0) G_DL = 0.0;
			if (G_SDL < 0.0) G_SDL = 0.0;
			
			// Importance sampling weight
			double w = 1.0;
			w *= rv_R.pdf(R) / rv_R_h.pdf(R);
			w *= rv_C_0Q.pdf(C_0Q) / rv_C_0Q_h.pdf(C_0Q);
			w *= rv_θ_E.pdf(θ_E) / rv_θ_E_h.pdf(θ_E);
			total_w += w;
			
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
