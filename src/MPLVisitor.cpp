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
#include "Programme.h"



/*** programs ***/

// ./mpl -e 'pgm = (ggga=(7-9); f; g;){7+6; 8-2; 2*4;}; 8-777; pgm'

//  ./mpl -e 'pgm = (gg = 6;){gg;}; pgm;' 

typedef enum
  {
   STATE_WAITING_FOR_PARAMS,	      
   STATE_IN_PARAMS,	      
   STATE_WAITING_FOR_INITIALISER,
   STATE_WAITING_FOR_STATEMENTS,
   STATE_IN_STATEMENTS
  } state_e;

antlrcpp::Any
MPLVisitor::visitMPLProgramme(MPLParser::MPLProgrammeContext *ctx)
{
  std::string latest_param;
  std::string latest_init;
  antlr4::tree::ParseTree *latest_pt = nullptr;
    
  antlrcpp::Any rc = Error (Error::ERROR_INTERNAL, " visitMPLProgramme");

  Programme *programme = new Programme (this);

  state_e state = STATE_WAITING_FOR_PARAMS;
  size_t nr_kids = ctx->children.size ();
  for (size_t i = 0; i < nr_kids; i++) {
    antlr4::tree::ParseTree *pt =  ctx->children[i];
#if 0
    Token *token = pt->start;
    //    antlrcpp::Any result = ctx->children[i]->accept(this);
    //    antlrcpp::Any result = ctx->children[i];
    std::cout << "toString " << i << " \"" << pt->toString () << "\"\n";;
    std::cout << "toTree " << pt->toStringTree (false) << std::endl;
    std::cout << "getText " << pt->getText () << std::endl;
    std::cout << "\n";
#if 0
    std::cout << "child " << i << " "
	      << " type \n\t"
	      << "(" << demanglep (pt) << ")"
	      << std::endl;
#endif
#endif

    //    std::string str = pt->getText ();
    //std::cout << i << " " << str << std::endl;
    
    switch(state) {
    case STATE_WAITING_FOR_PARAMS:
      {
	std::string ss =  pt->getText ();
	if ('(' == ss.front ()) {
	  state = STATE_IN_PARAMS;
	}
      }
      break;
    case STATE_IN_PARAMS:
      {
	std::string ss =  pt->getText ();
	switch (ss.front ()) {
	case '=': state = STATE_WAITING_FOR_INITIALISER; break;
	case ';':
	  programme->update_symtab (latest_param, latest_init, latest_pt);
	  state = STATE_IN_PARAMS;
	  break;
	case ')':
	  programme->update_symtab (latest_param, latest_init, latest_pt);
	  state = STATE_WAITING_FOR_STATEMENTS;
	  break;
	default:
	  latest_param = ss;
	  break;
	}
      }
      break;
    case STATE_WAITING_FOR_INITIALISER:
      {
	std::string ss =  pt->getText ();
	switch (ss.front ()) {
	case ';':
	  programme->update_symtab (latest_param, latest_init, latest_pt);
	  state = STATE_IN_PARAMS;
	  break;
	case ')':
	  programme->update_symtab (latest_param, latest_init, latest_pt);
	  state = STATE_WAITING_FOR_STATEMENTS;
	  break;
	default: latest_init = ss; break;
	}
      }
      break;
    case STATE_WAITING_FOR_STATEMENTS:
      {
	std::string ss =  pt->getText ();
	if ('{' == ss.front ()) {
	  state = STATE_IN_STATEMENTS;
	}
      }
      break;
    case STATE_IN_STATEMENTS:
      {
	std::string ss =  pt->getText ();
	switch (ss.front ()) {
	case ';': state = STATE_WAITING_FOR_PARAMS; break;
	case '}': state = STATE_WAITING_FOR_PARAMS; break;
	default: programme->push_statement (pt); break;
      }
      break;
    }
      
    
    }
  }

  std::vector<antlr4::tree::ParseTree *>*stmts =
    programme->get_stmts ();
  size_t n = stmts->size ();
  for (size_t i = i; i < n; i++) {
    antlr4::tree::ParseTree *px = (*stmts)[i];
    antlr4::tree::ParseTreeVisitor *that = programme->get_this ();
    //    antlrcpp::Any rc = px->accept (that);
  }

  rc = programme;
  
  return rc;
}

/*** indices ***/

antlrcpp::Any
MPLVisitor::visitMPLIndex(MPLParser::MPLIndexContext *ctx)
{
  antlrcpp::Any rc = Error(Error::ERROR_INTERNAL, " visitMPLIndex");
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
    value = get_symbol_value (ss);
    if  (value.get_typeinfo() == typeid(Programme *)) {
      Programme *programme = value.as<Programme *>();
      value = programme->run ();
    }
  }
  if (index.get_typeinfo() == typeid(std::string)) {
    std::string ss = index.as<std::string>();
    index = get_symbol_value (ss);
    if  (index.get_typeinfo() == typeid(Programme *)) {
      Programme *programme = value.as<Programme *>();
      index = programme->run ();
    }
  }


  // fixme-- what if any idx val = nan
  // fixme-- what if any idx is bool
  
  if (value.get_typeinfo() == typeid(std::vector<double>*)) {
    std::vector<double> *vals = value.as<std::vector<double>*>();
    if (index.get_typeinfo() == typeid(double)) {
      // vector [ scalar ]
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
      // vector [ vector ]
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
      // vector [ matrix ]
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
      rc = Error(Error::ERROR_UNKNOWN_DATA_TYPE, ", vector index");
    }
  }
  else if (value.get_typeinfo() == typeid(Matrix *)) {
    Matrix *mtx = value.as<Matrix *>();
    if (index.get_typeinfo() == typeid(std::vector<double>*)) {
      // matrix [ vector ]
      //  (2 3#::6)[ 1 2]
      //  rho idx = rhorho v
      // res = double
      std::vector<double> *idx = index.as<std::vector<double>*>();
      if (mtx->get_rhorho () == idx->size ()) {
	std::vector<size_t>*dims = mtx->get_dims ();
	if (mtx->get_rhorho () == dims->size ()) {
	  size_t i;
	  for (i = 0; i < idx->size (); i++) {
	    if ((*idx)[i] >= (*dims)[i]) break;
	  }
	  if (i == idx->size ()) {
	    double val = mtx->get_value (idx);
	    if (std::isnan (val)) {
	      rc = Error(Error::ERROR_DIMENSION_ERROR);
	      std::string em = mtx->get_errmsg ();
	      std::cout << em << std::endl;
	    }
	    rc = val;
	  }
	  else {
	    rc = Error (Error::ERROR_DIMENSION_MISMATCH,
			" Index element value exceeds matrix axis bounds.");
	  }
	}
	else {
	  rc = Error (Error::ERROR_INTERNAL,
		      "  Matrix rhorho != Matrix rho vector length.");
        }
      }
      else {
	rc = Error (Error::ERROR_DIMENSION_MISMATCH,
		    "  The length of the index vector must equal the order of the matrix.");
      }
    }
    else if (index.get_typeinfo() == typeid(Matrix *)) {
      // matrix [ matrix ]
      //  (2 3#::6)[ 1 2]
      //  rho idx = rhorho v
      // res = double
      
    }
    else {
      rc = Error(Error::ERROR_UNKNOWN_DATA_TYPE, ", matrix index");
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

    std::string str = token->getText ();
    std::replace (str.begin (), str.end (), '~', '-');
    double val = std::stod(str);
    array->push_back (val);
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

  std::string str = token->getText ();
  std::replace (str.begin (), str.end (), '~', '-');
  double val = std::stod(str);
  antlrcpp::Any value = val;

  return value;
}

/*** parens ***/
    
antlrcpp::Any
MPLVisitor::visitMPLParen(MPLParser::MPLParenContext *ctx)
{
  antlrcpp::Any rc = Error(Error::ERROR_INTERNAL, " visitMPLParen");
  
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
  auto token_idx = token->getType();
  antlrcpp::Any result = Error(Error::ERROR_INTERNAL, " visitMPLMonadic");
  
  size_t n = ctx->children.size();

  antlrcpp::Any qual_any;
  antlrcpp::Any right_any;
  if (n == 2) {
    qual_any  = nullptr;
    right_any = ctx->children[1]->accept(this);
  }
  else if (n == 5) {
    qual_any  = ctx->children[2]->accept(this);
    right_any = ctx->children[4]->accept(this);
  }
  else {
    result = Error(Error::ERROR_MALFORMED_EXPRESSION,
		   ", invalid parameter count");
    return result;
  }

  if (right_any.get_typeinfo() == typeid(std::string)) {
    std::string ss = right_any.as<std::string>();
    right_any = get_symbol_value (ss);
    if  (right_any.get_typeinfo() == typeid(Programme *)) {
      Programme *programme = right_any.as<Programme *>();
      right_any = programme->run ();
    }
  }

  if (right_any.isNotNull ()) {
    mfunc fcn = monadic_get_func (token_idx);
    if (fcn) (*fcn)(result, right_any, qual_any);
    else {
      result = Error(Error::ERROR_MALFORMED_EXPRESSION, ", invalid operation");
    }
  }
  else {
    result = Error(Error::ERROR_MALFORMED_EXPRESSION, ", null monadic operand");
  }

  return result;
}


  /*** dyadics ***/
  
antlrcpp::Any
MPLVisitor::visitMPLDyadic(MPLParser::MPLDyadicContext *ctx)
{
  Token *token = ctx->op ()->getStart ();
  auto token_idx = token->getType();
  antlrcpp::Any result = Error(Error::ERROR_INTERNAL, " visitMPLDyadic");
  std::string *left_copy = nullptr;

  antlrcpp::Any qual_any;
  antlrcpp::Any right_any;
  
  size_t n = ctx->children.size();

  if (n == 3) {
    qual_any  = nullptr;
    right_any = ctx->children[2]->accept(this);
  }
  else if (n == 6) {
    qual_any  = ctx->children[3]->accept(this);
    right_any = ctx->children[5]->accept(this);
  }
  else {
    result = Error(Error::ERROR_MALFORMED_EXPRESSION);
    return result;
  }

  antlrcpp::Any left_any = ctx->children[0]->accept(this);

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
      left_any = get_symbol_value (ss);
      if  (left_any.get_typeinfo() == typeid(Programme *)) {
	Programme *programme = left_any.as<Programme *>();
	left_any = programme->run ();
      }
    }
  }

  if (right_any.get_typeinfo() == typeid(std::string)) {
    std::string ss = right_any.as<std::string>();
    right_any = get_symbol_value (ss);
    if  (right_any.get_typeinfo() == typeid(Programme *)) {
      Programme *programme = right_any.as<Programme *>();
      right_any = programme->run ();
    }
  }

  if (left_any.isNotNull () && right_any.isNotNull ()) {
    dfunc fcn = dyadic_get_func (token_idx);
    if (fcn) (*fcn)(result, left_any, right_any, qual_any);
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
  antlrcpp::Any res;
  antlrcpp::Any newres;
  antlrcpp::Any rc = Error(Error::ERROR_NONE, ", statement");
  size_t n = ctx->children.size();
  //  std::cout << "n = " << n << std::endl;
  for (size_t i = 0; i < n; i++) {
    last_token = MPLLexer::DUMMY;
    
#ifdef SHOW_TRACE
    depth++;
#endif
    res  = ctx->children[i]->accept(this);
#if 0
    std::cout << "res " << i << " expr \"" << lexpr << "\""
	      << " type \n\t"
	      << "(" << demangle (res) << ")"
	      << std::endl;
#endif
#ifdef SHOW_TRACE
    --depth;
#endif
      
    //    if (last_token != MPLLexer::OpEqual) {
    if (res.isNotNull ()) {
      if (res.get_typeinfo() == typeid(std::string)) {
	// FIXME : create a lookup that takes an Any
	std::string ss = res.as<std::string>();
	newres = get_symbol_value (ss);
	    
#if 0
	std::cout << "newres " << i << " expr \"" << lexpr << "\""
		  << " type \n\t"
		  << "(" << demangle (newres) << ")"
		  << "at " << __FILE__ ", " << __LINE__ << std::endl;
#endif
	  
	if (newres.isNull ()) {
	  if (last_token != MPLLexer::OpEqual)
	    print_str (show, lexpr,ss);
	  rc = newres;
	}
#if 0
	else if (newres.get_typeinfo() == typeid(antlr4::tree::ParseTree*)) {
	  antlr4::tree::ParseTree *tree =
	    newres.as<antlr4::tree::ParseTree *>();
	  newres = tree->accept (this);
	}
#endif
	else {
	  if  (newres.get_typeinfo() == typeid(Programme *)) {
	    Programme *programme = newres.as<Programme *>();
	    newres = programme->run ();
	  }
	  if (last_token != MPLLexer::OpEqual)
	    print_val (show, lexpr, newres);
	  rc = newres;
	}
      }		
      else {
	if  (res.get_typeinfo() == typeid(Programme *)) {
	  Programme *programme = newres.as<Programme *>();
	  newres = programme->run ();
	  if (last_token != MPLLexer::OpEqual)
	    print_val (show, lexpr, newres);
	}
	else {
	  if (last_token != MPLLexer::OpEqual)
	    print_val (show, lexpr, res);
	}
      }
    }
    //    }
  }
  return rc;
}

  
