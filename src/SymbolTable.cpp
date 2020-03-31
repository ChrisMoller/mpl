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
  if (val.get_typeinfo() == typeid(antlr4::tree::ParseTree*)) {
    std::cout << "insert " << sym 
	      << "(" << demangle (val) << ")"
	      << std::endl;
    antlr4::tree::ParseTree *tree =
	    val.as<antlr4::tree::ParseTree *>();
    printf ("val = 0x%p tree = %p\n", val, tree);
  }
  antlrcpp::Any current = lookup (sym);
  if (current.isNotNull ()) symbols->erase (sym);
  antlrcpp::Any nv = antlrcpp::Any (val);
  symbols->insert ({sym, val});
}

antlrcpp::Any
SymbolTable::lookup (std::string sym)
{
  antlrcpp::Any rv (nullptr);
  if (!symbols->empty ()) {
    auto itr = symbols->find (sym);
    if (itr != symbols->end ()) {
      if (itr->second.get_typeinfo() == typeid(antlr4::tree::ParseTree*)) {
	std::cout << "lookup " << sym
		  << "(" << demangle (itr->second) << ")"
		  << std::endl;
	antlr4::tree::ParseTree *tree =
	  itr->second.as<antlr4::tree::ParseTree *>();
	printf ("itr->second = 0x%p tree = %p\n", itr->second, tree);
      }
      return itr->second;
      //      rv = antlrcpp::Any (itr->second);
    }
  }
  return rv;
}

SymbolTable global_symtab;

SymbolTable* get_global_symtab () {return &global_symtab; }

