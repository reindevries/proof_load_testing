// chapter_5_case_study_2_reliability_updating.cpp
//
// This code is similar to chapter_5_case_study_1.cpp, but in this case the resistance is predicted
// using the numerically simulated data. The resulting reliability indices are provided in Table 5.2
// of the dissertation.
//
// Author: Rein de Vries <rein.devries@tno.nl>
// Date: 10 November 2025

#include "rdv/msvcr.hpp"
#include "rdv/random_variable.hpp"
#include <vector>
#include <omp.h>

int main(int argc, char* argv[]) {
	// Initialisation
	rdv::tic();
	std::cout.setf(std::ios::left);
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
	
	// Define random variables
	auto rv_θ_R = rdv::rv_lognormal::mv(1.0, 0.15);
	auto rv_θ_E = rdv::rv_lognormal::mv(1.0, 0.10);
	auto rv_G_DL = rdv::rv_normal::mv(604.0, 0.05);
	auto rv_G_SDL = rdv::rv_normal::mv(143.0, 0.10);
	auto rv_C_0Q = rdv::rv_lognormal::mv(1.1, 0.10);
	auto rv_Q = rdv::rv_gumbel_max::mv(1301.0, 0.058);
	auto rv_θ_E_PL = rdv::rv_lognormal::mv(1.0, 0.10);
	auto rv_U = rdv::rv_normal(0.0, 1.0);
	
	// Characteristic traffic load
	double Q_k_LM1 = 2466.0;
	double Q_k_WIM = rv_Q.quantile(0.999);
	std::cout << "Q_k_LM1 = " << Q_k_LM1 << " kNm\n";
	std::cout << "Q_k_WIM = " << Q_k_WIM << " kNm\n";
	
	// Correlation between load effect calculations
	double ρ_θ_E = 0.7;
	double ρ_θ_E_compl = sqrt(1.0 - rdv::sq(ρ_θ_E));
	
	// Test loads and corresponding strength ratio statistics
	std::vector<double> m_Q_PL = {1125.0, 1875.0, 2625.0, 3188.0, 3281.0};
	std::vector<double> m_X    = {  2.42,   1.76,   1.36,   1.17,   1.15};
	std::vector<double> V_X    = { 0.090,  0.082,  0.068,  0.056,  0.054};
	
	// Assess reliability during and after each proof load test
	std::cout << "Evaluation of reliability via MCMC...\n";
	std::cout << std::setw(15) << "factor_WIM" << 
	             std::setw(15) << "factor_LM1" <<
	             std::setw(15) << "beta_during" << 
	             std::setw(15) << "beta_after" << 
	             std::setw(15) << "beta_lower" << 
	             std::setw(15) << "accept_rate" << "\n";
	std::cout.precision(3);
	
	// Loop over target loads
	for (size_t i = 0; i < size(m_Q_PL); i++) {
		// Random variables specific to the target load
		auto rv_Q_PL = rdv::rv_normal::mv(m_Q_PL[i], 0.02);
		auto rv_X = rdv::rv_normal::mv(m_X[i], V_X[i]);
		
		// MCMC: Metropolis-Hastings algorithm parameters
		std::vector<uint32_t> seeds = rdv::generate_seeds(omp_get_max_threads());
		size_t samples = size_t(1e9);  // total
		size_t burn_in = 1000;         // per thread
		double θ_R_width = 0.25 * rv_θ_R.mean();
		double X_width = 0.25 * rv_X.mean();
		
		// Keep track of failures
		size_t accept_count = 0, post_samples = 0;
		size_t failures_during = 0, failures_after = 0, failures_lower = 0;
		
		#pragma omp parallel reduction(+:accept_count, post_samples, \
				failures_during, failures_after, failures_lower)
		{
			// Thread-local random number generator and parameters
			std::mt19937 gen(seeds[omp_get_thread_num()]);
			size_t thread_samples = samples / omp_get_num_threads();
			
			double θ_R_curr = rv_θ_R.mean();
			double X_curr = rv_X.mean();
			double p_curr = rv_θ_R.pdf(θ_R_curr) * rv_X.pdf(X_curr);
			
			// Loop over thread samples
			for (size_t i_sample = 0; i_sample < thread_samples; i_sample++) {
				// Get new point from proposal (jumping) distribution
				double θ_R = θ_R_curr + rv_U(gen) * θ_R_width;
				double X = X_curr + rv_U(gen) * X_width;
				
				// Survival condition (likelihood)
				if (θ_R > 0.0 && X > 0.0 && θ_R * X > 1.0) {
					// Evaluate prior
					double p = rv_θ_R.pdf(θ_R) * rv_X.pdf(X);
					
					// Accept or reject jump
					double p_accept = std::min(1.0, p / p_curr);
					double p_rand = std::generate_canonical<double, size_t(-1)>(gen);
					
					if (p_rand < p_accept) {
						θ_R_curr = θ_R;
						X_curr = X;
						p_curr = p;
						accept_count++;
					}
				}
				
				// Make use of the posterior sample (θ_R_curr, X_curr)
				if (i_sample > burn_in) {
					post_samples++;
					
					// Standard normals for correlated model uncertainty of the load effect
					double U_θ_E = rv_U(gen);
					double U_θ_E_PL = ρ_θ_E * U_θ_E + ρ_θ_E_compl * rv_U(gen);
					
					// Realizations for each variable
					double θ_R = θ_R_curr;
					double X = X_curr;
					double θ_E = rv_θ_E.from_std_norm(U_θ_E);
					double G_DL = rv_G_DL(gen);
					double G_SDL = rv_G_SDL(gen);
					double C_0Q = rv_C_0Q(gen);
					double Q = rv_Q(gen);
					double θ_E_PL = rv_θ_E_PL.from_std_norm(U_θ_E_PL);
					double Q_PL = rv_Q_PL(gen);
					
					// Evaluate limit state function after a succesful test
					double E_PL = θ_E * (G_DL + G_SDL) + θ_E_PL * Q_PL;
					double E = θ_E * (G_DL + G_SDL + C_0Q * Q);
					double Z_after = θ_R * X * E_PL - E;
					if (Z_after < 0.0)
						failures_after++;
					
					// Evaluate limit state function during testing (plain MCS)
					θ_R = rv_θ_R(gen);
					X = rv_X(gen);
					double Z_during = θ_R * X * E_PL - E_PL;
					if (Z_during < 0.0)
						failures_during++;
						
					// Lower-bound approach (plain MCS)
					double Z_lower = θ_E_PL * Q_PL - θ_E * C_0Q * Q;
					if (Z_lower < 0.0)
						failures_lower++;
				}
			}
		}
		
		// Calculate the failure probability and reliability
		double P_f_during = double(failures_during) / post_samples;
		double P_f_after = double(failures_after) / post_samples;
		double P_f_lower = double(failures_lower) / post_samples;
		double β_during = -rdv::normal_quantile(P_f_during);
		double β_after = -rdv::normal_quantile(P_f_after);
		double β_lower = -rdv::normal_quantile(P_f_lower);
		
		std::cout << std::setw(15) << m_Q_PL[i] / Q_k_WIM <<
		             std::setw(15) << m_Q_PL[i] / Q_k_LM1 <<
		             std::setw(15) << β_during <<
		             std::setw(15) << β_after <<
		             std::setw(15) << β_lower <<
		             std::setw(15) << double(accept_count) / samples << "\n";
	}
	
	// Print elapsed time
	rdv::toc();
	return 0;
}
