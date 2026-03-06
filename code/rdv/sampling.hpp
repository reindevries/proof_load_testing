// sampling.hpp
// 
// This file provides a Latin hypercube sampling routine for generating samples in a multi-
// dimensional space. The result is stored in a vector of vectors, where each inner vector
// represents a random sample with values for each dimension.
// 
// Author: Rein de Vries
// Date: 8 January 2024

#ifndef __RDV_SAMPLING__
#define __RDV_SAMPLING__

#include <random>
#include <vector>

namespace rdv {

// Latin hypercube sampling
template <class V>
inline void latin_hyp_sampling(size_t dimension, size_t samples, std::vector<V>& result) {
	// Prepare random generator
	std::mt19937 gen(std::random_device{}());
	std::uniform_real_distribution<> dist;
	
	// Generate indexes and shuffle
	std::vector<std::vector<size_t>> indexes(dimension);
	for (size_t i = 0; i < dimension; i++) {
		indexes[i].resize(samples);
		for (size_t j = 0; j < samples; j++)
			indexes[i][j] = j;
		shuffle(begin(indexes[i]), end(indexes[i]), gen);
	}
	
	// Generate samples
	result.resize(samples);
	for (size_t j = 0; j < samples; j++) {
		result[j].resize(dimension);
		for (size_t i = 0; i < dimension; i++)
			result[j][i] = (indexes[i][j] + dist(gen)) / samples;
	}
}

}

#endif  // __RDV_SAMPLING__
