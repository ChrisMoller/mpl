#include <iostream>
#include <cxxabi.h>

#include <antlr4-runtime.h>
#include "MPLLexer.h"
#include "MPLParser.h"
#include "MPLParserBaseVisitor.h"

#include "DyadicFcns.h"
#include "MonadicFcns.h"

void set_params (std::string p, std::string e);

using namespace mpl;
using namespace antlr4;


/*** identifiers ***/

class  MPLVisitor : public MPLParserBaseVisitor {
public:
  MPLVisitor () {}
  MPLVisitor (bool v, std::string&str) {show = v; lexpr = str;}

  antlrcpp::Any
  visitMPLProgramme(MPLParser::MPLProgrammeContext *ctx) override;

  antlrcpp::Any
  visitMPLIndex(MPLParser::MPLIndexContext *ctx) override;
  
  antlrcpp::Any
  visitMPLIdentifier(MPLParser::MPLIdentifierContext *ctx) override;

  antlrcpp::Any
  visitMPLNumber(MPLParser::MPLNumberContext *ctx) override;

  antlrcpp::Any
  visitMPLVector(MPLParser::MPLVectorContext *ctx) override;

  antlrcpp::Any
  visitMPLParen(MPLParser::MPLParenContext *ctx) override;

  antlrcpp::Any
  visitMPLMonadic(MPLParser::MPLMonadicContext *ctx) override;

  antlrcpp::Any
  visitMPLDyadic(MPLParser::MPLDyadicContext *ctx) override;

  antlrcpp::Any
  visitMPLStatement(MPLParser::MPLStatementContext *ctx) override;

private:
  std::string lexpr;
  bool show = false;
#ifdef SHOW_TRACE
  int depth = 0;
#endif
  int last_token = MPLLexer::DUMMY;
};
