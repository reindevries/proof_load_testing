// general.hpp
//
// In this file commonly used functions are provided. Most of them are small and will be inlined by 
// the compiler such that no performance loss occurs.
//
// Author: Rein de Vries
// Date: 3 September 2024

#ifndef __RDV_GENERAL__
#define __RDV_GENERAL__

#define _USE_MATH_DEFINES  // for cmath

#include <cassert>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <type_traits>

// Debugging assertions

#define ASSERT assert

#if _DEBUG
  #define ASSERTRETURN(x) {auto y = x; ASSERT(y); return y;}
  #define VERIFY(x) ASSERT(x)
#else
  #define ASSERTRETURN(x) {return x;}
  #define VERIFY(x) x
#endif

// Parameter decorations

#ifndef __in
  #define __in
  #define __out
  #define __inout
  #define __in_opt
  #define __out_opt
  #define __inout_opt
  #define __callback
#endif

// Mathematical constants

#define M_1_SQRT2PI  0.398942280401432677940  // 1 / sqrt(2 * pi)
#define M_EM         0.577215664901532860607  // Euler–Mascheroni constant
#define M_SQRT3      1.73205080756887729353   // sqrt(3)
#define M_PHI        1.61803398874989484820   // (1 + sqrt(5)) / 2, golden ratio
#define M_2PI        6.28318530717958647693   // 2 * pi
#define M_LOG_2PI    1.83787706640934548356   // log(2 * pi)

#define M_INF        std::numeric_limits<double>::infinity()
#define M_NAN        std::numeric_limits<double>::quiet_NaN()

// Avoiding the C integer abs version
using std::abs;

namespace rdv {
	
// Use of standard functions
using std::abs;
using std::min;
using std::max;
using std::round;
using std::sqrt;

// Use of iterators
using std::size;
using std::begin;
using std::end;
using std::size_t;

// Mathematical functions

template <typename T, std::enable_if_t<std::is_arithmetic_v<T>, int> = 0>
inline T sq(T x) {
	return x*x;
}

template <typename T, std::enable_if_t<std::is_arithmetic_v<T>, int> = 0>
inline T ln(T x) {
	return log(x);
}

template <typename T>
inline T ln1p(T x) {
	return log1p(x);
}

template <typename T>
inline T pow10(T x) {
	return pow(10, x);
}

template <typename T>
inline int sgn(T x)
	{return int(T(0) < x) - int(x < T(0));}

template <typename T>
inline T max_abs(T a, T b) {
	T abs_a = abs(a);
	T abs_b = abs(b);
	return (abs_a > abs_b ? abs_a : abs_b);
}

template <typename T>
inline T min_abs(T a, T b) {
	T abs_a = abs(a);
	T abs_b = abs(b);
	return (abs_a < abs_b ? abs_a : abs_b);
}

template <typename T>
inline T min(T a, T b, T c) {
	return min(min(a, b), c);
}

template <typename T>
inline T min(T a, T b, T c, T d) {
	return min(min(a, b), min(c, d));
}

template <typename T>
inline T max(T a, T b, T c) {
	return max(max(a, b), c);
}

template <typename T>
inline T max(T a, T b, T c, T d) {
	return max(max(a, b), max(c, d));
}

inline double deg_to_rad(double x) {
	return x / 180.0 * M_PI;
}
	
inline double rad_to_deg(double x) {
	return x / M_PI * 180.0;
}

inline double round(double x, int n) {
	double m = pow10((double)n);
	return round(x * m) / m;
}

// Enable generalisations for complex numbers

template <typename T>
using complex_part_t = decltype(abs(T(0)));

#define REAL_FUNCS(type)   \
inline type conj(type x) { \
	return x;              \
}                          \
                           \
inline type real(type x) { \
	return x;              \
}                          \
                           \
inline type imag(type x) { \
	return 0;              \
}

REAL_FUNCS(float)
REAL_FUNCS(double)
REAL_FUNCS(short)
REAL_FUNCS(int)
REAL_FUNCS(int64_t)

// Timing functions

inline void _tictoc(bool start) {
	static std::chrono::system_clock::time_point start_time;
	if (start)
		start_time = std::chrono::system_clock::now();
	else {
		std::chrono::milliseconds dt = std::chrono::duration_cast<
				std::chrono::milliseconds>(std::chrono::system_clock::now() - start_time);
		std::cout << "Time elapsed: " << dt.count() / 1000.0 << " seconds\n";
	}
}

inline void tic() {
	_tictoc(true);
}

inline void toc() {
	_tictoc(false);
}

}

#endif  // __RDV_GENERAL__
