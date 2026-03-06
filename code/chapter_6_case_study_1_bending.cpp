// chapter_6_case_study_1_bending.cpp
//
// This code calculates transfer factors for case study 1, considering bending. The flag
// 'correlated_Q' switches between the correlated and uncorrelated situations. Using interpolation, 
// the target load corresponding to beta = 4 may be determined. The ratio between the load required 
// in each configuration and the n = N case results in the transfer factors. The results for bending 
// are provided in Table 6.3 of the dissertation.
//
// Author: Rein de Vries <rein.devries@tno.nl>
// Date: 10 November 2025

#include "rdv/msvcr.hpp"
#include "rdv/random_variable.hpp"
#include "rdv/vector_ext.hpp"
#include <Eigen/Dense>
#include <omp.h>

inline void calculate_reliability(size_t N, size_t n, double m_Q_PL, double& β) {
	// Settings
	constexpr bool correlated_Q = false;
	
	// Random variables
	auto rv_θ_R = rdv::rv_lognormal::mv(1.0, 0.05);    // value for bending
	auto rv_m_R = rdv::rv_lognormal::mv(2.55, 0.5);    // 1.5 times the mean load effect
	auto rv_V_R = rdv::rv_lognormal::mv(0.05, 0.1);    // value for bending
	auto rv_θ_E = rdv::rv_lognormal::mv(1.0, 0.1);     // model uncertainty load effect
	auto rv_G = rdv::rv_normal::mv(0.7, 0.07);         // value for bending
	auto rv_C_0Q = rdv::rv_lognormal::mv(1.0, 0.1);    // stat uncertainty traffic load, no DAF
	auto rv_Q = rdv::rv_gumbel_max::mv(1.0, 0.025);    // value for bending, small spans
	auto rv_θ_E_PL = rdv::rv_lognormal::mv(1.0, 0.1);  // model uncert, plus correlation 0.7
	auto rv_Q_PL = rdv::rv_normal::mv(m_Q_PL, 0.02);   // proof load effect
	auto rv_U = rdv::rv_normal(0.0, 1.0);
	
	// Prepare random sampling
	std::vector<uint32_t> seeds = rdv::generate_seeds(omp_get_max_threads());
	size_t samples = size_t(1e8);
	size_t post_samples = 0, failures = 0;
	
	// Cholesky decomposition of permanent load correlation matrix
	Eigen::MatrixXd R_G = Eigen::MatrixXd::Constant(N, N, 0.7);
	R_G.diagonal().setOnes();
	Eigen::LLT<Eigen::MatrixXd> llt_G(R_G);
	auto L_G = llt_G.matrixL();
	
	#pragma omp parallel reduction(+:post_samples, failures)
	{
		// Thread objects
		std::mt19937 gen(seeds[omp_get_thread_num()]);
		Eigen::VectorXd R(N), G(N), u_G(N), v_G(N);
		
		// Sampling procedure
		#pragma omp for
		for (size_t i_sample = 0; i_sample < samples; i_sample++) {
			// Fully correlated between components
			double θ_R = rv_θ_R(gen);
			double m_R = rv_m_R(gen);
			double V_R = rv_V_R(gen);
			double C_0Q = rv_C_0Q(gen);
			
			// Model uncertainty (proof) load effect
			double v_θ_E = rv_U(gen);
			double v_θ_E_PL = 0.7 * v_θ_E + 0.714143 * rv_U(gen);
			double θ_E = rv_θ_E.from_std_norm(v_θ_E);
			double θ_E_PL = rv_θ_E_PL.from_std_norm(v_θ_E_PL);
			
			// Correlated standard normals for permanent loads
			for (size_t i = 0; i < N; i++)
				u_G[i] = rv_U(gen);
			v_G.noalias() = L_G * u_G;
			
			// Per-component realizations of resistance and permanent loads
			auto rv_R = rdv::rv_lognormal::mv(m_R, V_R);
			for (size_t i = 0; i < N; i++) {
				R[i] = rv_R(gen);
				G[i] = rv_G.from_std_norm(v_G[i]);
			}
			
			// Subject all N components to a year of traffic load (in-service proven strength)
			bool failed = false;
			double Q = rv_Q(gen);
			for (size_t i = 0; i < N; i++) {
				Q = (correlated_Q ? Q : rv_Q(gen));
				double Z = θ_R * R[i] - θ_E * (G[i] + C_0Q * Q);
				if (Z < 0.0) {
					failed = true;
					break;
				}
			}
			if (failed)
				continue;
			
			// Subject n components to proof load test
			failed = false;
			for (size_t i = 0; i < n; i++) {
				double Q_PL = rv_Q_PL(gen);  // uncorrelated, describes operation variability
				double Z = θ_R * R[i] - θ_E * G[i] - θ_E_PL * Q_PL;
				if (Z < 0.0) {
					failed = true;
					break;
				}
			}
			if (failed)
				continue;
			
			// Subject all N components to a year of traffic load (future)
			failed = false;
			post_samples++;
			Q = rv_Q(gen);
			for (size_t i = 0; i < N; i++) {
				Q = (correlated_Q ? Q : rv_Q(gen));
				double Z = θ_R * R[i] - θ_E * (G[i] + C_0Q * Q);
				if (Z < 0.0) {
					failed = true;
					break;
				}
			}
			if (failed)
				failures++;
		}
	}
	
	// Calculate reliability of the object containing N components
	double P_f = double(failures) / post_samples;
	β = (P_f > 0.0 ? -rdv::normal_quantile(P_f) : 6.0);
	
	// Display result
	std::cout << std::setw(10) << m_Q_PL << β << "\n";
}

int main(int argc, char *argv[]) {
	// Initialisation
	rdv::tic();
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
	std::cout.setf(std::ios::left);
	std::cout.precision(3);
	
	// Sample size and population sizes
	std::vector<double> ns = {1, 2, 3, 4, 5, 6, 8};
	std::vector<double> Ns = {1, 2, 3, 4, 5, 6, 8, 10, 20};
	
	// Consider the various configurations
	std::cout << "m_Q_PL    Reliability\n";
	for (size_t n : ns) {
		for (size_t N : Ns) {
			if (N < n)
				continue;
			
			std::cout << "n=" << n << " N=" << N << "\n";
			
			// Target loads depend on the number of components
			double start = 1.4;
			double step = (N <= 3 ? 0.1 : 0.2);
			
			// Calculate reliability for various proof load values
			double β = -M_INF;
			for (size_t i_PL = 0; i_PL < 5; i_PL++) {
				double m_Q_PL = start + step * i_PL;
				if (β < 4.0)
					calculate_reliability(N, n, m_Q_PL, β);
				else
					std::cout << std::setw(10) << m_Q_PL << 6.0 << "\n";
			}
		}
		
		std::cout << "\n";
	}
	
	// Print elapsed time
	rdv::toc();
	return 0;
}
