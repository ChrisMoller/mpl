#ifndef MAIN_H
#define MAIN_H
#include <cxxabi.h>
#include <antlr4-runtime.h>

#pragma once

#include "SymbolTable.h"


#ifndef demangle
#define demangle(res)  abi::__cxa_demangle(res.get_typeinfo().name(), \
					   0, 0, nullptr)
#define demanglep(res) abi::__cxa_demangle(res->get_typeinfo().name(), \
					   0, 0, nullptr)
#define demangled(res) abi::__cxa_demangle(res, 0, 0, nullptr)
#endif

bool isTestMode ();

typedef enum
  {
   SOURCE_NONE,
   SOURCE_CMDLINE,
   SOURCE_FILE
  } source_e;

source_e get_source ();
bool     isFromFile ();
bool     isFromCmdLine ();

void push_symbol_table (SymbolTable *st);
void pop_symbol_table ();
void insert_symbol_value (std::string sym, antlrcpp::Any &val);
antlrcpp::Any get_symbol_value (std::string sym);

#endif // MAIN_H
