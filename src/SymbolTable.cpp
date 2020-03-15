#include <antlr4-runtime.h>
#include "SymbolTable.h"
#include "main.h"

using namespace antlr4;

SymbolTable::SymbolTable ()
{
  symbols = new std::map<std::string, antlrcpp::Any>;
  parent = nullptr;
}

SymbolTable:: ~SymbolTable ()
{
  if (parent) {
    delete parent;
    parent = nullptr;
  }
  if (symbols) {
    delete symbols;
    symbols = nullptr;
  }
}

void
SymbolTable::insert (std::string sym, antlrcpp::Any &val)
{
  antlrcpp::Any current = lookup (sym);
  if (current.isNotNull ()) symbols->erase (sym);
  antlrcpp::Any nv = antlrcpp::Any (val);
  symbols->insert ({sym, val});
}

antlrcpp::Any
SymbolTable::lookup (std::string sym)
{
  antlrcpp::Any rv;
  if (!symbols->empty ()) {
    auto itr = symbols->find (sym);
    if (itr != symbols->end ()) {
      rv = antlrcpp::Any (itr->second);
    }
  }
  return rv;
}

SymbolTable global_symtab;

SymbolTable* get_global_symtab () {return &global_symtab; }

