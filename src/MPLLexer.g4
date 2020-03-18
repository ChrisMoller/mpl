lexer grammar MPLLexer;

// These are all supported lexer sections:

// Lexer file header. Appears at the top of h + cpp files.
// Use e.g. for copyrights.

@lexer::header {/* lexer header section */}

// Appears before any #include in h + cpp files.

@lexer::preinclude {/* lexer precinclude section */}

// Follows directly after the standard #includes in h + cpp files.

@lexer::postinclude {
/* lexer postinclude section */
#include "main.h"
#pragma GCC diagnostic ignored "-Wunused-parameter"
}

// Directly preceds the lexer class declaration in the h file
// (e.g. for additional types etc.).

@lexer::context {/* lexer context section */}

// Appears in the public part of the lexer in the h file.
@lexer::members {/* public lexer declarations section */
//bool canTestFoo() { return true; }
//bool isItFoo() { return true; }
//bool isItBar() { return true; }

//void myFooLexerAction() { /* do something*/ };
//void myBarLexerAction() { /* do something*/ };
}

// Appears in the private part of the lexer in the h file.
@lexer::declarations {
/* private lexer declarations/members section */
}

// Appears in line with the other class member definitions in the cpp file.
@lexer::definitions {/* lexer definitions section */}

channels { CommentsChannel, DirectiveChannel }

tokens {
	DUMMY
}

//Return: 'return';
//Continue: 'continue';


Number	          : '~'?(Real | Imag | Complex | NaN);
fragment Complex  : Real[+\-]Imag ;
fragment Imag	  : Real[ij] ;
fragment Real	  : '.'Digit+ | Digit+('.'Digit*)?([Ee][~+\-]?Digit+)? ;
fragment Digit	  : [0-9];
fragment NaN	  : 'NAN' | 'NaN' | 'nan' | 'inf';

OpStar		  	: '*';	    	// nullptr             dyadicStar
OpSlash		  	: '/';		// monadicSlash	       dyadicSlash
OpPlus		  	: '+';		// nullptr	       dyadicPlus
OpMinus		  	: '-';		// monadicMinus	       dyadicMinus
OpDollar		: '$';		// monadicTranspose    dyadicTranspose
OpPercent		: '%';		// nullptr	       nullptr
OpHat		  	: '^';		// nullptr	       dyadicHat
OpLn			: 'ln';		// monadicLn	       nullptr
OpLog			: 'log';	// monadicLog	       dyadicLog
OpExp			: 'exp';	// monadicExp	       nullptr
OpRoot			: 'root';	// monadicRoot	       dyadicRoot
OpSin			: 'sin';	// monadicSin	       nullptr
OpCos			: 'cos';	// monadicCos	       nullptr
OpTan			: 'tan';	// monadicTan	       nullptr
OpAsin			: 'asin';	// monadicAsin	       nullptr
OpAcos			: 'acos';	// monadicAcos	       nullptr
OpAtan			: 'atan';	// monadicAtan	       dyadicAtan
OpSind			: 'sind';	// monadicSind	       nullptr
OpCosd			: 'cosd';	// monadicCosd	       nullptr
OpTand			: 'tand';	// monadicTand	       nullptr
OpAsind			: 'asind';	// monadicAsind	       nullptr
OpAcosd			: 'acosd';	// monadicAcosd	       nullptr
OpAtand			: 'atand';	// monadicAtand	       dyadicAtand
OpBar			: '|';		// monadicBar	       nullptr
OpBang			: '!';		// monadicBang	       nullptr
OpCeil			: 'ceil';	// monadicCeil	       nullptr
OpFloor			: 'floor';	// monadicFloor	       nullptr
OpRound			: 'round';	// monadicRound	       nullptr
OpLeftAngle 	  	: '<';		// monadicGradeDown    nullptr
OpRightAngle	  	: '>';		// monadicGradeUp      nullptr
OpQMark			: '?';		// monadicRand	       nullptr
OpEqual			: '=';		// nullptr	       dyadicEqual
OpPound			: '#';		// monadicShape	       dyadicShape
OpComma			: ',';		// nullptr	       nullptr
OpColon			: ':';		// nullptr	       nullptr
OpBSStar		: '\\*';	// nullptr	       dyadicMatMult
OpColonColon		: '::';		// monadicRange	       dyadicRange
OpQEqual		: '?=';		// nullptr	       dyadicTestEq
OpQLeftAngle	  	: '?<';		// nullptr	       dyadicTestLT
OpQRightAngle  		: '?>';		// nullptr	       dyadicTestGT
OpQLeftAngleEqual	: '?<=';	// nullptr	       dyadicTestLE
OpQRightAngleEqual 	: '?>=';	// nullptr	       dyadicTestGE
OpQBangEqual	  	: '?!=';	// nullptr	       dyadicTestNE
OpQBangLeftAngle	: '?!<';	// nullptr	       dyadicTestGE
OpQBangRightAngle	: '?!>';	// nullptr	       dyadicTestLE
OpQBangLeftAngleEqual	: '?!<=';	// nullptr	       dyadicTestGT
OpQBangRightAngleEqual	: '?!>=';	// nullptr	       dyadicTestLT
OpDet			: 'det';	// monadicDeterminant  nullptr
OpBSSlash		: '\\/';	// monadicInverse      dyadicMatSolve
OpSlashPlus		: '/+';		// monadicSum	       nullptr
OpSlashStar		: '/*';		// monadicProduct      nullptr
OpBSI			: '\\I';	// monadicIdentity     dyadicIdentity

And		: 'and';
Semicolon	: ';';
OpenPar		: '(';
ClosePar	: ')';
OpenSquare	: '[';
CloseSquare	: ']';
OpenCurly	: '{';
CloseCurly	: '}';
//Ampersand	: '&' -> type(DUMMY);

//Vector  : Number Number+ ;

ID: LETTER (LETTER | Digit)*;
fragment LETTER : [_a-zA-Z];

String		: '"' .*? '"';

LASTTOKEN	: '><';
Comment 	: '<!--' .*? '-->' -> skip;

EOL		: {isFromCmdLine()}? [\n\r]+ ;

WS		: ({isFromCmdLine()}?  [ \t]+
		| {isFromFile()}?   [ \n\t]+) -> skip
		;

//WS		: [ \t]+ -> channel(99);

//EOL		: {isFromCmdLine()}? [\n\r]+ ;
//EOL		: [\n\r]+ ;

//WS		: ({isFromCmdLine()}?  [ \t]+
//		|  EOF
//		| {isFromFile()}?   [ \n\t]+) -> skip
//		;
//
//WS		: ({isFromCmdLine()}? (EOF | [ \t]+)
//		| {isFromFile()}?   (EOF | [ \n\t]+)) -> channel(99)
//		;

