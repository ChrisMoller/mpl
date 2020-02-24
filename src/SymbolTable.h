#include "antlr4-runtime.h"

class SymbolTable {
public:
  SymbolTable () {}

  ~SymbolTable () {}

  void insert (std::string sym, antlrcpp::Any &val);
  antlrcpp::Any lookup (std::string sym);

private:
  std::map<std::string, antlrcpp::Any>symbols;
};

SymbolTable* get_global_symtab ();
