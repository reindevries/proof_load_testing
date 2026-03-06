// msvcr.hpp
// 
// In this file common mathematical functions from the Microsoft Visual C++ Runtime library, if 
// installed, are made available. These functions provide better performance than those supplied
// with GNU C++. If the Microsoft compiler is used, then making use of this module is unnecessary. 
// This file is meant to be used on the Windows platform. Only double precision is supported for 
// the function parameters.
// 
// From C++0x onwards static function variable initialization is thread-safe. Therefore these
// functions may be safely used in a multi-threaded environment.
// 
// Author: Rein de Vries
// Date: 6 November 2025

#ifndef __RDV_MSVCR__
#define __RDV_MSVCR__

#include "windows.hpp"
#include "general.hpp"

#ifndef _MSC_VER

// Loading of module
inline HMODULE msvcr_init_handle() {
	HMODULE handle = LoadLibrary(_T("msvcr130.dll"));
	if (handle == NULL) handle = LoadLibrary(_T("msvcr120.dll"));
	if (handle == NULL) handle = LoadLibrary(_T("msvcr110.dll"));
	if (handle == NULL) handle = LoadLibrary(_T("msvcr100.dll"));
	if (handle == NULL) handle = LoadLibrary(_T("msvcr90.dll"));
	if (handle == NULL)
		std::cerr << "Error: Could not load msvcr130.dll or a substitute.\n";
	return handle;
}

// Get module handle
inline HMODULE msvcr_handle() {
	static HMODULE handle = msvcr_init_handle();
	return handle;
}

// Macro to define functions with 1 parameter
#define MSVCR_FUNC_1(name)                                                          \
inline double name(double x) {                                                      \
	typedef double(*Func)(double);                                                  \
	static Func p = reinterpret_cast<Func>(GetProcAddress(msvcr_handle(), #name));  \
	return p ? p(x) : M_NAN;                                                        \
}

// Macro to define functions with 2 parameters
#define MSVCR_FUNC_2(name)                                                          \
inline double name(double x, double y) {                                            \
	typedef double(*Func)(double, double);                                          \
	static Func p = reinterpret_cast<Func>(GetProcAddress(msvcr_handle(), #name));  \
	return p ? p(x, y) : M_NAN;                                                     \
}

// Function list
MSVCR_FUNC_1(atan)
MSVCR_FUNC_2(atan2)
MSVCR_FUNC_1(cos)
MSVCR_FUNC_1(cosh)
MSVCR_FUNC_1(erf)
MSVCR_FUNC_1(erfc)
MSVCR_FUNC_1(exp)
MSVCR_FUNC_1(exp2)
MSVCR_FUNC_1(expm1)
MSVCR_FUNC_1(log)
MSVCR_FUNC_1(log10)
MSVCR_FUNC_1(log1p)
MSVCR_FUNC_1(log2)
MSVCR_FUNC_2(pow)
MSVCR_FUNC_1(sin)
MSVCR_FUNC_1(sinh)
MSVCR_FUNC_1(sqrt)
MSVCR_FUNC_1(tan)
MSVCR_FUNC_1(tanh)
MSVCR_FUNC_1(tgamma)

#endif  // _MSC_VER

#endif  // __RDV_MSVCR__
