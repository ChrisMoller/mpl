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
@lexer::declarations {/* private lexer declarations/members section */}

// Appears in line with the other class member definitions in the cpp file.
@lexer::definitions {/* lexer definitions section */}

channels { CommentsChannel, DirectiveChannel }

tokens {
	DUMMY
}

//Return: 'return';
//Continue: 'continue';


Number	  : '~'?(Real | Imag | Complex | NaN);
fragment Complex  : Real[+\-]Imag ;
fragment Imag	  : Real[ij] ;
fragment Real	  : '.'Digit+ | Digit+('.'Digit*)?([Ee][+\-]?Digit+)? ;
fragment Digit	  : [0-9];
fragment NaN	  : 'NAN' | 'NaN' | 'nan' | 'inf';


OpStar		  	: '*';
OpSlash		  	: '/';
OpPlus		  	: '+';
OpMinus		  	: '-';
OpDollar		: '$';
OpPercent		: '%';
OpHat		  	: '^';
OpLn			: 'ln';
OpLog			: 'log';
OpExp			: 'exp';
OpRoot			: 'root';
OpSin			: 'sin';
OpCos			: 'cos';
OpTan			: 'tan';
OpAsin			: 'asin';
OpAcos			: 'acos';
OpAtan			: 'atan';
OpSind			: 'sind';
OpCosd			: 'cosd';
OpTand			: 'tand';
OpAsind			: 'asind';
OpAcosd			: 'acosd';
OpAtand			: 'atand';
OpBar			: '|';
OpBang			: '!';
OpCeil			: 'ceil';
OpFloor			: 'floor';
OpRound			: 'round';
OpLeftAngle 	  	: '<';
OpRightAngle	  	: '>';
OpQMark			: '?';
OpEqual			: '=';
OpPound			: '#';
OpComma			: ',';
OpColon			: ':';
OpColonColon		: '::';
OpQEqual		: '?=';
OpQLeftAngle	  	: '?<';
OpQRightAngle  		: '?>';
OpQLeftAngleEqual	: '?<=';
OpQRightAngleEqual 	: '?>=';
OpQBangEqual	  	: '?!=';
OpQBangLeftAngle	: '?!<';
OpQBangRightAngle	: '?!>';
OpQBangLeftAngleEqual	: '?!<=';
OpQBangRightAngleEqual	: '?!>=';

And		: 'and';
Semicolon	: ';';
OpenPar		: '(';
ClosePar	: ')';
OpenSquare	: '[';
CloseSquare	: ']';
//Ampersand	: '&' -> type(DUMMY);

//Vector  : Number Number+ ;

ID: LETTER (LETTER | '0'..'9')*;
fragment LETTER : [_a-zA-Z];

String		: '"' .*? '"';

LASTTOKEN	: '><';
Comment 	: '<!--' .*? '-->' -> skip;
EOL		: [\n\r]+ ;
WS		: [ \t]+ -> channel(99);

