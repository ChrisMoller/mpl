parser grammar MPLParser;

options {
	tokenVocab = MPLLexer;
}

// These are all supported parser sections:

// Parser file header. Appears at the top in all parser related files.
// Use e.g. for copyrights.

@parser::header {/* parser/listener/visitor header section */}

// Appears before any #include in h + cpp files.

@parser::preinclude {/* parser precinclude section */}

// Follows directly after the standard #includes in h + cpp files.

@parser::postinclude {
/* parser postinclude section */
#include "main.h"
#pragma GCC diagnostic ignored "-Wunused-parameter"
}

// Directly preceeds the parser class declaration in the h file
//(e.g. for additional types etc.).

@parser::context {/* parser context section */}

// Appears in the private part of the parser in the h file.
// The function bodies could also appear in the definitions section,
// but I want to maximize
// Java compatibility, so we can also create a Java parser from this grammar.
// Still, some tweaking is necessary after the Java file generation
// (e.g. bool -> boolean).

@parser::members {
/* public parser declarations/members section */
//void myAction(antlr4::Token *tok)
//{
//  std::cout << "myAction\n";
//}
//bool myAction() { return true; }
//bool doesItBlend() { return true; }
//void cleanUp() {}
//void doInit() {}
//void doAfter() {}
}

// Appears in the public part of the parser in the h file.

@parser::declarations {
/* private parser declarations section */
}

// Appears in line with the other class member definitions in the cpp file.

@parser::definitions {
/* parser definitions section */
}

// Additionally there are similar sections for (base)listener
// and (base)visitor files.
@parser::listenerpreinclude {/* listener preinclude section */}
@parser::listenerpostinclude {/* listener postinclude section */}
@parser::listenerdeclarations {
  /* listener public declarations/members section */
}
@parser::listenermembers {/* listener private declarations/members section */}
@parser::listenerdefinitions {/* listener definitions section */}

@parser::baselistenerpreinclude {/* base listener preinclude section */}
@parser::baselistenerpostinclude {/* base listener postinclude section */}
@parser::baselistenerdeclarations {
  /* base listener public declarations/members section */
}
@parser::baselistenermembers {
  /* base listener private declarations/members section */
}
@parser::baselistenerdefinitions {/* base listener definitions section */}

@parser::visitorpreinclude {/* visitor preinclude section */}
@parser::visitorpostinclude {/* visitor postinclude section */}
@parser::visitordeclarations {
/* visitor public declarations/members section */
}
@parser::visitormembers {
/* visitor private declarations/members section */
}
@parser::visitordefinitions {/* visitor definitions section */}

@parser::basevisitorpreinclude {/* base visitor preinclude section */}
@parser::basevisitorpostinclude {/* base visitor postinclude section */}
@parser::basevisitordeclarations {
  /* base visitor public declarations/members section */
}
@parser::basevisitormembers {
  /* base visitor private declarations/members section */
}
@parser::basevisitordefinitions {/* base visitor definitions section */}

// Actual grammar start.

main	  : stat+ EOF
	  | Semicolon
	  ;

stat: (expr)? eos      				# MPLStatement
      ;

eos	: {isFromCmdLine()}? Semicolon | EOL
	| {isFromFile()}? Semicolon
	;


expr	: <assoc = right> expr op expr		# MPLDyadic
      	| <assoc = right> op expr 		# MPLMonadic
      	| <assoc = right> op
	    OpenSquare op_or_expr CloseSquare expr 	# MPLQualMono
	| OpenPar expr ClosePar			# MPLParen
	| expr OpenSquare expr CloseSquare	# MPLIndex 
    	| identifier = id			# MPLIdentifier
  	| Number     				# MPLNumber
  	| Number Number+		   	# MPLVector
    	| String				# MPLString
	| OpenPar (id (OpEqual expr)? eos)* ClosePar
	  OpenCurly stat+ CloseCurly            #MPLProgramme
//	| OpenPar (expr Semicolon)* ClosePar
//	  OpenCurly stat+ CloseCurly            #MPLProgramme
	;


//param   : id (OpEqual expr)?;
//param   : p=id {set_param($p.text);} (OpEqual e=expr {set_init ($e.text);})?;

// e.start:  type antlr4::Token

// MPLProgrammeContext : public ExprContext
//          std::vector<StatContext *> stat();
//          StatContext* stat(size_t i);
//          std::vector<ExprContext *> expr();
//          ExprContext* expr(size_t i);
// ExprContext : public antlr4::ParserRuleContext
//          nothing useful
// ParserRuleContext : public RuleContext
//          Token *start;
//          Token *stop;
// RuleContext : public tree::ParseTree
//           virtual std::string getText() override;
//           virtual std::string toStringTree(bool pretty = false) override;
//           virtual std::string toString() override;
//           std::vector<ParseTree *> children;
// 

op_or_expr : op | expr ;
	
op	: OpStar
	| OpSlash
	| OpPlus
	| OpMinus
	| OpDollar
	| OpPercent
	| OpHat
	| OpLn
	| OpLog
	| OpExp
	| OpRoot
	| OpSin
	| OpCos
	| OpTan
	| OpAsin
	| OpAcos
	| OpAtan
	| OpSind
	| OpCosd
	| OpTand
	| OpAsind
	| OpAcosd
	| OpAtand
	| OpBar
	| OpBang
	| OpCeil
	| OpFloor
	| OpRound
	| OpColonColon
	| OpEqual
	| OpLeftAngle
	| OpRightAngle
	| OpQMark
	| OpQEqual
	| OpBar
	| OpPound
	| OpComma
	| OpColon
	| OpQBangLeftAngleEqual
	| OpQBangRightAngleEqual
	| OpQLeftAngleEqual
	| OpQRightAngleEqual
	| OpQBangEqual
	| OpQBangLeftAngle
	| OpQBangRightAngle
	| OpQLeftAngle
	| OpQRightAngle
	| OpBSStar
	| OpDet
	| OpBSSlash
	| OpSlashPlus
	| OpSlashStar
	| OpBSI
	;

id	: ID;

//Array	: OpenCurly el += Number (Comma el += Number)* CloseCurly;
//idarray : OpenCurly element += id (Comma element += id)* CloseCurly;
//any	: t = .;
