#include <cxxabi.h>

#ifndef demangle
#define demangle(res) abi::__cxa_demangle(res.get_typeinfo().name(), 0, 0, nullptr)

#define demanglep(res) abi::__cxa_demangle(res->get_typeinfo().name(), 0, 0, nullptr)
#endif

bool isTestMode ();
