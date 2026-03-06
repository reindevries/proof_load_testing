// vector_ext.hpp
//
// Implements vector extensions such as summation, minimum, etc. It commits the ultimate sin by
// defining extra functions in the std namespace, but this is done consciously to provide a seamless
// experience when working with std::vector.
//
// Author: Rein de Vries
// Date: 20 August 2024

#ifndef __RDV_VECTOR_EXT__
#define __RDV_VECTOR_EXT__

#define VECTOR_NPOS SIZE_MAX

#include "general.hpp"
#include <type_traits>
#include <vector>

namespace std {

template <typename T>
inline size_t size(const vector<T>& vec) {
	return vec.size();
}

// Numerical properties

template <typename T>
inline T sum(const vector<T>& vec) {
	T res = (T)0;
	for (size_t i = 0; i != vec.size(); ++i)
		res += vec[i];
	return res;
}

template <typename T>
inline T min(const vector<T>& vec) {
	if (vec.size() > 0)	{
		T res = vec[0];
		for (size_t i = 1; i != vec.size(); ++i)
			if (vec[i] < res)
				res = vec[i];
		return res;
	}
	
	ASSERT(false);
	return std::numeric_limits<T>::quiet_NaN();
}

template <typename T>
inline T max(const vector<T>& vec) {
	if (vec.size() > 0)	{
		T res = vec[0];
		for (size_t i = 1; i != vec.size(); ++i)
			if (vec[i] > res)
				res = vec[i];
		return res;
	}
	
	ASSERT(false);
	return std::numeric_limits<T>::quiet_NaN();
}

template <typename T>
inline T norm(const vector<T>& vec) {
	T res = (T)0;
	for (size_t i = 0; i != vec.size(); ++i)
		res += vec[i]*vec[i];
	return sqrt(res);
}

// Element-wise math operations

template <typename T, enable_if_t<std::is_arithmetic_v<T>, int> = 0>
inline vector<T> log(const vector<T>& vec) {
	vector<T> result(size(vec));
	for (size_t i = 0; i != size(vec); ++i)
		result[i] = ::log(vec[i]);
	return result;
}

// Finding and sorting

template <typename T, typename T2>
inline size_t find(const vector<T>& vec, T2&& val) {
	for (size_t i = 0; i != vec.size(); ++i)
		if (vec[i] == val)
			return i;
	return VECTOR_NPOS;
}

template <typename T, typename T2>
inline size_t find_first_of(const vector<T>& vec, T2&& values) {
	for (auto&& val : values) {
		size_t res = find(vec, val);
		if (res != VECTOR_NPOS)
			return res;
	}
	return VECTOR_NPOS;
}

template <typename T>
inline size_t min_index(const vector<T>& vec) {
	T res = (T)0;
	size_t res_i = 0;
	if (vec.size() > 0)	{
		res = vec[0];
		res_i = 0;
		
		for (size_t i = 1; i != vec.size(); ++i) {
			if (vec[i] < res) {
				res = vec[i];
				res_i = i;
			}
		}
	}
	return res_i;
}

template <typename T>
inline size_t max_index(const vector<T>& vec) {
	T res = (T)0;
	size_t res_i = 0;
	if (vec.size() > 0)	{
		res = vec[0];
		res_i = 0;
		
		for (size_t i = 1; i != vec.size(); ++i) {
			if (vec[i] > res) {
				res = vec[i];
				res_i = i;
			}
		}
	}
	return res_i;
}

template <typename T>
inline void sort(vector<T>& vec, bool small_to_large = true) {
	if (small_to_large)
		sort(vec.begin(), vec.end(), less<T>());
	else
		sort(vec.begin(), vec.end(), greater<T>());
}

// Statistical properties

template <typename T>
inline T mean(const vector<T>& vec) {
	return sum(vec) / vec.size();
}

template <typename T>
inline T variance_p(const vector<T>& vec, T m) {
	if (vec.size() > 0)	{
		T res = (T)0;
		for (size_t i = 0; i != vec.size(); ++i)
			res += rdv::sq(vec[i] - m);
		res /= vec.size();
		return res;
	}
	
	ASSERT(false);
	return std::numeric_limits<T>::quiet_NaN();
}

template <typename T>
inline T variance_p(const vector<T>& vec) {
	return variance_p(vec, mean(vec));
}

template <typename T>
inline T variance_s(const vector<T>& vec, T m) {
	if (vec.size() > 1) {
		T res = (T)0;
		for (size_t i = 0; i != vec.size(); ++i)
			res += rdv::sq(vec[i] - m);
		res /= (vec.size() - 1);
		return res;
	}
	
	ASSERT(false);
	return std::numeric_limits<T>::quiet_NaN();
}

template <typename T>
inline T variance_s(const vector<T>& vec) {
	return variance_s(vec, mean(vec));
}

template <typename T>
inline T stddev_p(const vector<T>& vec, T m) {
	return sqrt(variance_p(vec, m));
}

template <typename T>
inline T stddev_p(const vector<T>& vec) {
	return stddev_p(vec, mean(vec));
}

template <typename T>
inline T stddev_s(const vector<T>& vec, T m) {
	return sqrt(variance_s(vec, m));
}

template <typename T>
inline T stddev_s(const vector<T>& vec) {
	return stddev_s(vec, mean(vec));
}

template <typename T>
inline T cov_p(const vector<T>& vec, T m) {
	return stddev_p(vec, m) / abs(m);
}

template <typename T>
inline T cov_p(const vector<T>& vec) {
	return cov_p(vec, mean(vec));
}

template <typename T>
inline T cov_s(const vector<T>& vec, T m) {
	return stddev_s(vec, m) / abs(m);
}

template <typename T>
inline T cov_s(const vector<T>& vec) {
	return cov_s(vec, mean(vec));
}

// Correlation, assuming normal distribution
// https://en.wikipedia.org/wiki/Pearson_correlation_coefficient

template <typename T>
inline T correlation_p(const vector<T>& vec1, const vector<T>& vec2) {
	if (vec1.size() != vec2.size())
		{ASSERT(false); return std::numeric_limits<T>::quiet_NaN();}

	size_t n = vec1.size();

	if (n < 1)
		{ASSERT(false); return std::numeric_limits<T>::quiet_NaN();}

	T m1 = mean(vec1);
	T m2 = mean(vec2);
	
	T s1 = stddev_p(vec1, m1);
	T s2 = stddev_p(vec2, m2);

	T res = T(0);
	for (size_t i = 0; i != n; ++i)
		res += (vec1[i] - m1) * (vec2[i] - m2);
	res /= n * s1 * s2;

	return res;
}

template <typename T>
inline T correlation_s(const vector<T>& vec1, const vector<T>& vec2) {
	if (vec1.size() != vec2.size())
		{ASSERT(false); return std::numeric_limits<T>::quiet_NaN();}
	size_t n = vec1.size();
	if (n < 2)
		{ASSERT(false); return std::numeric_limits<T>::quiet_NaN();}
	return correlation_p(vec1, vec2) * n / (n - 1);
}

// Resizing of multidimensional vectors
// Here std::vector may contain any type of vector that has the member function 'resize'

template <typename V>
void resize(vector<V>& vec, size_t n1, size_t n2) {
	vec.resize(n1);
	for (auto&& x : vec)
		x.resize(n2);
}

template <typename V, typename T>
void resize(vector<V>& vec, size_t n1, size_t n2, T val) {
	vec.resize(n1);
	for (auto&& x : vec)
		x.resize(n2, val);
}

template <typename V>
void resize(vector<vector<V>>& vec, size_t n1, size_t n2, size_t n3) {
	vec.resize(n1);
	for (auto&& x : vec)
		resize(x, n2, n3);
}

template <typename V, typename T>
void resize(vector<vector<V>>& vec, size_t n1, size_t n2, size_t n3, T val) {
	vec.resize(n1);
	for (auto&& x : vec)
		resize(x, n2, n3, val);
}

template <typename V>
void resize(vector<vector<vector<V>>>& vec, size_t n1, size_t n2, size_t n3, size_t n4) {
	vec.resize(n1);
	for (auto&& x : vec)
		resize(x, n2, n3, n4);
}

template <typename V, typename T>
void resize(vector<vector<vector<V>>>& vec, size_t n1, size_t n2, size_t n3, size_t n4, T val) {
	vec.resize(n1);
	for (auto&& x : vec)
		resize(x, n2, n3, n4, val);
}

// Multiplication and division by a scalar

template <typename T, typename T2>
vector<T> operator *(const vector<T>& vec, T2 val) {
	vector<T> res(vec.size());
	for (size_t i = 0; i != vec.size(); ++i)
		res[i] = vec[i] * val;
	return res;
}

template <typename T, typename T2>
vector<T> operator *(T2 val, const vector<T>& vec) {
	vector<T> res(vec.size());
	for (size_t i = 0; i != vec.size(); ++i)
		res[i] = vec[i] * val;
	return res;
}

template <typename T, typename T2>
vector<T> operator /(const vector<T>& vec, T2 val) {
	vector<T> res(vec.size());
	for (size_t i = 0; i != vec.size(); ++i)
		res[i] = vec[i] / val;
	return res;
}

// Addition and subtraction of vectors (of the same length)

template <typename T>
vector<T> operator +(const vector<T>& vec1, const vector<T>& vec2) {
	ASSERT(vec1.size() == vec2.size());
	vector<T> res(min(vec1.size(), vec2.size()));
	for (size_t i = 0; i != res.size(); ++i)
		res[i] = vec1[i] + vec2[i];
	return res;
}

template <typename T>
vector<T> operator -(const vector<T>& vec1, const vector<T>& vec2) {
	ASSERT(vec1.size() == vec2.size());
	vector<T> res(min(vec1.size(), vec2.size()));
	for (size_t i = 0; i != res.size(); ++i)
		res[i] = vec1[i] - vec2[i];
	return res;
}

template <typename T>
vector<T> operator -(const vector<T>& vec) {
	vector<T> res(vec.size());
	for (size_t i = 0; i != res.size(); ++i)
		res[i] = -vec[i];
	return res;	
}

// Get a part of the vector

template <typename T>
vector<T> part(const vector<T>& vec, size_t offset, size_t size) {
	vector<T> result;
	if (offset + size > vec.size())
		{ASSERT(false); return result;}
	result.reserve(size);
	size_t end = offset + size;
	for (size_t i = offset; i != end; ++i)
		result.push_back(vec[i]);
	return result;
}

template <typename T>
vector<T> head(const vector<T>& vec, size_t size) {
	vector<T> result;
	if (size > vec.size())
		{ASSERT(false); return result;}
	result.reserve(size);
	for (size_t i = 0; i != size; ++i)
		result.push_back(vec[i]);
	return result;
}

template <typename T>
vector<T> tail(const vector<T>& vec, size_t size) {
	vector<T> result;
	if (size > vec.size())
		{ASSERT(false); return result;}
	result.reserve(size);
	size_t offset = vec.size() - size;
	for (size_t i = offset; i != vec.size(); ++i)
		result.push_back(vec[i]);
	return result;
}

// Stream output

template <typename T>
ostream& operator <<(ostream& os, const vector<T>& v) {
	os << "{";
	for (auto&& it = v.begin(); it != v.end(); ++it) 	{
		if (it != v.begin())
			os << ", ";
		os << *it;
	}
	os << "}";

	return os;
}

// Append another vector

template <typename T>
void append(vector<T>& a, const vector<T>& b) {
	a.reserve(a.size() + b.size());
	for (auto&& it = b.begin(); it != b.end(); ++it)
		a.push_back(*it);
}

}

// Custom OpenMP reductions for plus (+) operation

#ifdef __GNUG__
  #define DO_PRAGMA(x) _Pragma(#x)
  
  #define STD_VEC_PLUS(type) \
    DO_PRAGMA(omp declare reduction(+ : std::vector<type> : \
      std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), \
      omp_out.begin(), std::plus<type>())) initializer(omp_priv = omp_orig))
  
  STD_VEC_PLUS(double)
  STD_VEC_PLUS(float)
  STD_VEC_PLUS(int)
  STD_VEC_PLUS(size_t)
#endif

#endif // __RDV_VECTOR_EXT__
