#include <iostream>
#include <cxxabi.h>

#include <antlr4-runtime.h>
#include "main.h"
#include "MPLLexer.h"
#include "MPLParser.h"
#include "MPLParserBaseVisitor.h"
#include "MPLVisitor.h"
#include "Error.h"
#include "Matrix.h"

#include "SymbolTable.h"
#include "DyadicFcns.h"
#include "MonadicFcns.h"
#include "Print.h"

/*** indices ***/

antlrcpp::Any
MPLVisitor::visitMPLIndex(MPLParser::MPLIndexContext *ctx)
{
  antlrcpp::Any rc;
#if 0
  std::cout << "n = "
	    <<  ctx->children.size ()
	    << std::endl;

  for (size_t i = 0; i < ctx->children.size (); i++) {
    antlrcpp::Any res = ctx->children[i]->accept(this);
    std::cout << "res " << i << " expr \"" << lexpr << "\""
	      << " type \n\t"
	      << "(" << demangle (res) << ")"
	      << std::endl;
  }
#endif

  antlrcpp::Any value = ctx->children[0]->accept(this);
  
#if 0
  std::cout << "value type \n\t"
	    << "(" << demangle (value) << ")"
	    << std::endl;
#endif
  
  antlrcpp::Any index = ctx->children[2]->accept(this);
  
#if 0
  std::cout << "index type \n\t"
	    << "(" << demangle (index) << ")"
	    << std::endl;
#endif

  if (value.get_typeinfo() == typeid(std::string)) {
    std::string ss = value.as<std::string>();
    value = get_global_symtab ()->lookup (ss);
  }
  if (index.get_typeinfo() == typeid(std::string)) {
    std::string ss = index.as<std::string>();
    index = get_global_symtab ()->lookup (ss);
  }


  // fixme-- what if any idx val = nan
  // fixme-- what if any idx is bool
  
  if (value.get_typeinfo() == typeid(std::vector<double>*)) {
    std::vector<double> *vals = value.as<std::vector<double>*>();
    if (index.get_typeinfo() == typeid(double)) {
      // (1 2 3)[ 2 ]
      // 0 <= ix < rho-v
      // res = double, i.e., res = typeid index
      double idx = index.as<double>();
      if (idx >= 0 && idx < vals->size ())
	rc = (*vals)[size_t (idx)];
      else {
	rc = Error(Error::ERROR_OUT_OF_RANGE, ", index");
      }
    }
    else if (index.get_typeinfo() == typeid(std::vector<double>*)) {
      // (1 2 3)[ 0 2 ]
      // foreach ix,  0 <= ix < rho-v
      // res = vector double, i.e., res = typeid index,
      // rho-res = rho idx
      std::vector<double> *ix = index.as<std::vector<double>*>();
      size_t rho_v = vals->size ();
      for (size_t i = 0; i < ix->size (); i++) {
	if ((*ix)[i] >= rho_v) {
	  rc = Error(Error::ERROR_OUT_OF_RANGE, ", index");
	  break;
	}
      }
      std::vector<double> *res = new std::vector<double>;
      res->resize (ix->size ());
      for (size_t i = 0; i < res->size (); i++)
	(*res)[i] = (*vals)[(*ix)[i]];
      rc = res;
    }
    else if (index.get_typeinfo() == typeid(Matrix*)) {
      // (2 3 4)[ 2 3# 0 2 1... ]
      // foreach ix,  0 <= ix < rho-v
      // res = matrix, i.e., res = typeid index
      // rhorho res = rhorho idx, rho-res = rho idx
      Matrix *idx = index.as<Matrix*>();
      size_t rho_v = vals->size ();
      std::vector<double> *idx_vals = idx->get_data ();
      for (size_t i = 0; i < idx_vals->size (); i++) {
	if ((*idx_vals)[i] >= rho_v) {
	  rc = Error(Error::ERROR_OUT_OF_RANGE, ", index");
	  break;
	}
      }
      std::vector<size_t> *idx_dims = idx->get_dims ();
      std::vector<size_t> *ndims =
	new std::vector<size_t>(idx_dims->size ());
      for (size_t i = 0; i < idx_dims->size (); i++)
	(*ndims)[i] = (*idx_dims)[i];
      std::vector<double_t> *ndata = new std::vector<double>;
      ndata->resize(idx_vals->size ());
      for (size_t i = 0; i < idx_vals->size (); i++)
	(*ndata)[i] = (*vals)[int((*idx_vals)[i])];
      Matrix *res = new Matrix (ndims, ndata);
      rc = res;
    }
    else {
      rc = Error(Error::ERROR_UNKNOWN_DATA_TYPE, ", index");
    }
  }
  else if (value.get_typeinfo() == typeid(Matrix *)) {
    Matrix *mtx = value.as<Matrix *>();
    if (index.get_typeinfo() == typeid(std::vector<double>)) {
      //  (2 3#::6)[ 1 2]
      //  rho idx = rhorho v
      // res = double
    }
    if (index.get_typeinfo() == typeid(Matrix *)) {
      //  (2 3#::6)[ 1 2]
      //  rho idx = rhorho v
      // res = double
      
    }
    else {
      rc = Error(Error::ERROR_UNKNOWN_DATA_TYPE, ", index");
    }
  }
  else {
  std::cout << "value type \n\t"
	      << "(" << demangle (value) << ")"
	      << std::endl;
    rc = Error(Error::ERROR_UNKNOWN_DATA_TYPE, ", argument");
  }

  return rc;
}

/*** vectors ***/

antlrcpp::Any
MPLVisitor::visitMPLVector(MPLParser::MPLVectorContext *ctx)
{
  std::vector<antlr4::tree::TerminalNode *> vector = ctx->Number ();
  size_t n = vector.size ();
  std::vector<double> *array = new std::vector<double>;
  for (size_t i = 0; i < n; i++) {
    tree::TerminalNode *ety = ctx->Number (i);
    Token *token = ety->getSymbol();
    
    int offset = 0;
    if ('~' == token->getText()[0]) offset = 1;
    std::string str = token->getText ();
    double val = std::stod(str.substr(offset));
    if (offset == 1) val = -val;
    array->push_back (val);
    //std::cout << "ety " << i << " = " << val << std::endl;
  }

  antlrcpp::Any value = array;
  return value;
}


/*** identifiers ***/

antlrcpp::Any
MPLVisitor::visitMPLIdentifier(MPLParser::MPLIdentifierContext *ctx)
{
  tree::TerminalNode *id = ctx->identifier->ID ();
  Token *token = id->getSymbol();
#ifdef SHOW_TRACE
  std::string indent = std::string (3 * depth, ' ');
  std::cout << indent << "ID " << token->getText () << std::endl;
#endif
  antlrcpp::Any value = token->getText ();

  return value;
}


/*** numbers ***/

antlrcpp::Any
MPLVisitor::visitMPLNumber(MPLParser::MPLNumberContext *ctx)
{
  tree::TerminalNode *number = ctx->Number ();
  Token *token = number->getSymbol();
#ifdef SHOW_TRACE
  std::string indent = std::string (3 * depth, ' ');
  std::cout << indent << "Nr " << token->getText () << std::endl;
#endif
  int offset = 0;
  if ('~' == token->getText()[0]) offset = 1;
  std::string str = token->getText ();
  double val = std::stod(str.substr(offset));
  if (offset == 1) val = -val;
  antlrcpp::Any value = val;

  return value;
}

/*** parens ***/
    
antlrcpp::Any
MPLVisitor::visitMPLParen(MPLParser::MPLParenContext *ctx)
{
  antlrcpp::Any rc;	// just to provide a null
  
  size_t n = ctx->children.size();
  
  for (size_t i = 0; i < n; i++) {
    antlrcpp::Any result = ctx->children[i]->accept(this);
    
    /********
	     FIXME
	     
	     this will just return the first evaluated expr
	     but is would be possible to have it return a vector
	       
    *******/
    if (result.isNotNull()) return result;
  }

  return rc;
}

/*** monadics ***/

antlrcpp::Any
MPLVisitor::visitMPLMonadic(MPLParser::MPLMonadicContext *ctx)
{
  Token *token = ctx->op ()->getStart ();
#ifdef SHOW_TRACE
  std::string indent = std::string (3 * depth, ' ');
  std::cout << indent << "Mo " << token->getText () << std::endl;
#endif
  auto token_idx = token->getType();
  antlrcpp::Any result;
  
  size_t n = ctx->children.size();
  
  if (n != 2) {
    result = Error(Error::ERROR_MALFORMED_EXPRESSION);
    return result;
  }

#ifdef SHOW_TRACE
  depth++;
#endif
  antlrcpp::Any right_any = ctx->children[1]->accept(this);
#ifdef SHOW_TRACE
  --depth;
#endif
  if (right_any.get_typeinfo() == typeid(std::string)) {
    std::string ss = right_any.as<std::string>();
    right_any = get_global_symtab ()->lookup (ss);
  }

  if (right_any.isNotNull ()) {
    mfunc fcn = monadic_get_func (token_idx);
    if (fcn) (*fcn)(result, right_any);
    else {
      result = Error(Error::ERROR_MALFORMED_EXPRESSION);
    }
  }
  else {
    result = Error(Error::ERROR_MALFORMED_EXPRESSION);
  }

  return result;
}


  /*** dyadics ***/
  
antlrcpp::Any
MPLVisitor::visitMPLDyadic(MPLParser::MPLDyadicContext *ctx)
{
  Token *token = ctx->op ()->getStart ();
#ifdef SHOW_TRACE
  std::string indent = std::string (3 * depth, ' ');
  std::cout << indent << "Dy " << token->getText ()
	    << " depth = " << depth << std::endl;
#endif
  auto token_idx = token->getType();
  antlrcpp::Any result;
  std::string *left_copy = nullptr;

  size_t n = ctx->children.size();

  if (n != 3) {
    result = Error(Error::ERROR_MALFORMED_EXPRESSION);
    return result;
  }

#ifdef SHOW_TRACE
  depth++;
#endif
  antlrcpp::Any left_any = ctx->children[0]->accept(this);
#ifdef SHOW_TRACE
  depth--;
#endif
  if (token_idx == MPLLexer::OpEqual) {
    if (left_any.get_typeinfo() == typeid(std::string)) {
      std::string ss = left_any.as<std::string>();
      left_copy = new std::string (ss);
      left_any = left_copy;
    }
  }
  else {
    if (left_any.get_typeinfo() == typeid(std::string)) {
      std::string ss = left_any.as<std::string>();
      left_any = get_global_symtab ()->lookup (ss);
    }
  }

#ifdef SHOW_TRACE
  depth++;
#endif
  antlrcpp::Any right_any = ctx->children[2]->accept(this);
#ifdef SHOW_TRACE
  depth--;
#endif
  if (right_any.get_typeinfo() == typeid(std::string)) {
    std::string ss = right_any.as<std::string>();
    right_any = get_global_symtab ()->lookup (ss);
  }

  if (left_any.isNotNull () && right_any.isNotNull ()) {
    dfunc fcn = dyadic_get_func (token_idx);
    if (fcn) (*fcn)(result, left_any, right_any);
    else {
      result = Error(Error::ERROR_MALFORMED_EXPRESSION);
      return result;
    }
  }
  else {
    result = Error(Error::ERROR_MALFORMED_EXPRESSION);
    return result;
  }

  if (left_copy) {
    left_copy->clear ();
    delete left_copy;
  }

  last_token = token_idx;
  return result;
}

/*** statement ***/
  
antlrcpp::Any
MPLVisitor::visitMPLStatement(MPLParser::MPLStatementContext *ctx)
{
  size_t n = ctx->children.size();
  for (size_t i = 0; i < n; i++) {
    last_token = MPLLexer::DUMMY;
    
#ifdef SHOW_TRACE
    depth++;
#endif
    antlrcpp::Any res  = ctx->children[i]->accept(this);
#ifdef SHOW_TYPES
    std::cout << "res " << i << " expr \"" << lexpr << "\""
	      << " type \n\t"
	      << "(" << demangle (res) << ")"
	      << std::endl;
#endif
#ifdef SHOW_TRACE
    --depth;
#endif
      
    if (last_token != MPLLexer::OpEqual) {
      if (res.isNotNull ()) {
	if (res.get_typeinfo() == typeid(std::string)) {
	  // FIXME : create a lookup that takes an Any
	  std::string ss = res.as<std::string>();
	  antlrcpp::Any newres = get_global_symtab ()->lookup (ss);
	    
#if SHOW_TYPES
	  std::cout << "newres " << i << " expr \"" << lexpr << "\""
		    << " type \n\t"
		    << "(" << demangle (newres) << ")"
		    << std::endl;
#endif
	  
	  if (newres.isNull ())
	    print_str (show, lexpr,ss);
	  else 
	    print_val (show, lexpr, newres);
	}	
	else 
	  print_val (show, lexpr, res);
      }
    }
  }
  antlrcpp::Any rc = 0; // some sort of okay
  return rc;
}

