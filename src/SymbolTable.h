#include "antlr4-runtime.h"

#pragma once

class SymbolTable {
public:
  SymbolTable ();

  ~SymbolTable ();

  SymbolTable *push (std::map<std::string, antlrcpp::Any>*symb);
  void insert (std::string sym, antlrcpp::Any &val);
  antlrcpp::Any lookup (std::string sym);

private:
  SymbolTable *parent = nullptr;
  std::map<std::string, antlrcpp::Any>*symbols = nullptr;
};

