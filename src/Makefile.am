bin_PROGRAMS = mpl
mpl_SOURCES = \
               DyadicFcns.cpp DyadicFcns.h \
               Error.cpp Error.h \
               main.cpp main.h \
               Matrix.cpp Matrix.h \
               MonadicFcns.cpp MonadicFcns.h \
               MPLVisitor.cpp MPLVisitor.h \
               Print.cpp Print.h \
               SymbolTable.cpp SymbolTable.h

nodist_mpl_SOURCES =  MPLLexer.cpp MPLParser.cpp
BUILT_SOURCES = MPLLexer.cpp MPLParser.cpp MPLParserBaseVisitor.h
CLEANFILES = MPLLexer.cpp MPLParser.cpp MPLParserBaseVisitor.h \
             MPLLexer.interp MPLLexer.tokens MPLParserBaseListener.cpp \
             MPLParserBaseListener.h MPLParserBaseVisitor.cpp \
             MPLParser.h MPLParser.interp MPLParserListener.cpp \
             MPLParserListener.h MPLParser.tokens MPLParserVisitor.cpp \
             MPLParserVisitor.h MPLLexer.h 

DISTCLEANFILES = .deps/*.P

AUTOMAKE_OPTIONS = dejagnu

#$(bin_PROGRAMS)-eval.o: eval.c opsdata.h printext.h node.h


AM_CPPFLAGS = -std=c++11 -Wno-attributes \
   -I/usr/local/antlr4/runtime/Cpp/run/usr/local/include/antlr4-runtime 

# AM_CPPFLAGS += --coverage

mpl_LDFLAGS = -L/usr/local/antlr4/runtime/Cpp/run/usr/local/lib -lantlr4-runtime

ANTLRFLAGS = -visitor -listener

ANTLRJAR = /home/moller/javastuff/antlr-4.8-complete.jar

MPLLexer.h MPLLexer.cpp  : MPLLexer.g4
	java  -jar $(ANTLRJAR) -Dlanguage=Cpp $(ANTLRFLAGS) \
           -package mpl $<

MPLParser.h MPLParser.cpp  : MPLParser.g4
	java  -jar $(ANTLRJAR) -Dlanguage=Cpp $(ANTLRFLAGS) \
          -package mpl $<




# -leak-check=full
#  algrind --tool=memcheck --leak-check=full --leak-resolution=high -
#                -num-callers=20 --log-file=vgdump ./gfig
#  
#VALGRIND_FLAGS = --gen-suppressions=all --run-libc-freeres=yes \
#    --suppressions=valgrind.supp
#VALGRIND_FLAGS = --gen-suppressions=all --run-libc-freeres=yes

valgrind.log : clc
	valgrind $(VALGRIND_FLAGS) --log-file=$@ ./clc

save:
	gitpush

clean-local: clean-local-check
.PHONY: clean-local-check
clean-local-check:
	rm -rf clc-parse.c clc-parse.h clc-lex.c
