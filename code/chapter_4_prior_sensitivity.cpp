// chapter_4_prior_sensitivity.cpp
//
// Using this code, the sensitivity to various choices for the prior distribution of the resistance
// is explored. The distribution is changed by uncommenting the relevant rv_R_hat definition. The
// result of the calculation is the posterior distribution of the resistance, assuming survival of
// the proof load effect. The results of these calculations are displayed in Figures 4.2 and 4.3 of
// the dissertation.
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
	auto rv_R_hat = rdv::rv_normal::mv(1150e3, 0.5);
	//auto rv_R_hat = rdv::rv_normal::mv(1725e3, 0.5);
	//auto rv_R_hat = rdv::rv_uniform(0.0, 3450e3);
	//auto rv_R_hat = rdv::rv_triangular(0.0, 3450e3, 0.0);
	auto rv_Q = rdv::rv_gumbel_max::mv(1150e3, 0.025);
	auto rv_C_0Q = rdv::rv_lognormal::mv(1.1, 0.1);
	auto rv_Q_PL = rdv::rv_normal::mv(1800e3, 0.01);
	
	// Set up bins for posterior distribution plot
	size_t bin_count = 200;
	double bin_width = 3500e3 / bin_count;
	std::vector<double> bins(bin_count, 0.0);
	
	// Prepare random sampling
	std::vector<uint32_t> seeds = rdv::generate_seeds(omp_get_max_threads());
	size_t samples = size_t(1e9);
	size_t failures = 0;
	size_t survived_PL = 0;
	std::cout << "Performing Bayesian Monte Carlo simulation...\n";
	
	// Start parallel threads
	#pragma omp parallel reduction(+:failures, survived_PL, bins)
	{
		// Thread-local random number generator
		std::mt19937 gen(seeds[omp_get_thread_num()]);
		
		#pragma omp for
		for (size_t i_sample = 0; i_sample < samples; i_sample++) {
			// Sample variables
			double R_hat = rv_R_hat(gen);
			double Q = rv_Q(gen);
			double C_0Q = rv_C_0Q(gen);
			double Q_PL = rv_Q_PL(gen);
			
			// Check if proof load test is survived
			double Z_PL = R_hat - Q_PL;
			if (Z_PL >= 0.0) {
				// Count survivals for conditional failure probability
				survived_PL++;
				
				// Collect posterior distribution samples
				if (R_hat >= 0.0) {
					size_t i_bin = size_t(R_hat / bin_width);
					if (i_bin < bin_count)
						bins[i_bin] += 1.0;
				}
				
				// Evaluate regular traffic load situation
				double Z = R_hat - C_0Q * Q;
				if (Z < 0.0)
					failures++;
			}
		}
	}
	
	// Output posterior distribution of residual resistance
	double dx = bin_width / 1e3;
	std::cout << std::setw(15) << "R_hat" << "PDF\n";
	for (size_t i_bin = 0; i_bin < bin_count; i_bin++) {
		double x = (i_bin + 0.5) * dx;
		double y = (bins[i_bin] / survived_PL) / dx;
		std::cout << std::setw(15) << x << y << "\n";
	}
	
	// Output reliability after proof load survival
	double P_f = double(failures) / survived_PL;
	double beta = -rdv::normal_quantile(P_f);
	std::cout << "beta = " << beta << "\n";
	
	// Print elapsed time
	rdv::toc();
	return 0;
}
