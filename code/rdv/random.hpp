// random.hpp
//
// Utilities for random number generation, including seed generation and standard uniform 
// distribution sampling from both std::mt19937 and std::mt19937_64 generators.
//
// Author: Rein de Vries
// Date: 10 November 2025

#ifndef __RDV_RANDOM__
#define __RDV_RANDOM__

#include "general.hpp"
#include <random>
#include <vector>

namespace rdv {

inline std::vector<uint32_t> generate_seeds(size_t n, uint32_t start = 0) {
	std::vector<uint32_t> seeds(n);
	if (start == 0) {
		std::random_device rd;
		std::seed_seq seq{rd(), rd()};
		seq.generate(begin(seeds), end(seeds));
	}
	else {
		std::seed_seq seq{start};
		seq.generate(begin(seeds), end(seeds));
	}
	return seeds;
}

inline double random_uniform(std::mt19937& gen) {
	return ((((uint64_t{gen()} << 32) | gen()) >> 12) + 0.5) / (1ULL << 52);
}

}

#endif  // __RDV_RANDOM__
