// chapter_6_load_correlation_decay.cpp
//
// In this code a Monte Carlo simulation is performed to determine the decay of the correlation
// between load effects when the reference period increases. Gumbel distributed random variables 
// with given initial correlation and block sizes are defined using vectors. The output is provided 
// in Figure 6.5 of the dissertation.
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
	SetConsoleOutputCP(CP_UTF8);
	
	// Random variables
	rdv::rv_normal U(0.0, 1.0);
	rdv::rv_gumbel_max E(2.0, 0.3);  // any parameter combination will do
	
	// Prepare Monte Carlo simulation
	std::vector<uint32_t> seeds = rdv::generate_seeds(omp_get_max_threads());
	size_t samples = size_t(1e8);
	std::vector<size_t> block_sizes = {3, 10, 30, 100, 300, 1000, 3000, 10000, 30000, 100000};
	std::vector<double> initial_correlations = {0.95, 0.9, 0.8, 0.7, 0.6, 0.5};
	
	// Print header
	std::cout << std::setw(15) << "n" << std::setw(15) << "rho = 1";
	for (double rho : initial_correlations) {
		std::stringstream ss;
		ss << "rho = " << rho;
		std::cout << std::setw(15) << ss.str();
	}
	std::cout << "\n";
	
	// Initial block correlations (n = 1)
	std::cout << std::setw(15) << 1 << std::setw(15) << 1;
	for (double rho : initial_correlations)
		std::cout << std::setw(15) << rho;
	std::cout << "\n";
	
	// Loop over block sizes
	for (size_t n : block_sizes) {
		std::cout << std::setw(15) << n << std::setw(15) << 1.0;
		
		// Random variable with distribution describing the n-block maximum
		auto E_n = E.shifted(n);
		
		// Loop over initial correlations
		for (double rho : initial_correlations) {
			double ohr = sqrt(1.0 - rho*rho);
			
			// Determine the number of blocks to sample
			size_t blocks = std::min(samples / n, size_t(1e7));
			
			// Loop over blocks and calculate the correlation
			double rho_n = 0.0;
			
			#pragma omp parallel
			{
				// Thread objects
				std::mt19937 gen(seeds[omp_get_thread_num()]);
				
				#pragma omp for reduction(+:rho_n)
				for (size_t i_block = 0; i_block < blocks; i_block++) {
					// Obtain samples and save max
					double e_1_max = -M_INF;
					double e_2_max = -M_INF;
					for (size_t i = 0; i < n; i++) {
						double u_1 = U(gen);
						double u_2 = rho * u_1 + ohr * U(gen);
						double e_1 = E.from_std_norm(u_1);
						double e_2 = E.from_std_norm(u_2);
						e_1_max = std::max(e_1, e_1_max);
						e_2_max = std::max(e_2, e_2_max);
					}
					
					// Standard normals of maxima
					double u_1_n = E_n.to_std_norm(e_1_max);
					double u_2_n = E_n.to_std_norm(e_2_max);
					
					// Calculate correlation
					rho_n += u_1_n * u_2_n;
				}
			}
			
			// Output correlation
			rho_n /= blocks;
			std::cout << std::setw(15) << rho_n;
		}
		
		std::cout << "\n";
	}
	
	// Print elapsed time
	rdv::toc();
	return 0;
}
