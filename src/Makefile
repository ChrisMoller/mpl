INCLUDES = -I/usr/local/antlr4/runtime/Cpp/run/usr/local/include/antlr4-runtime 
    LIBS = -L/usr/local/antlr4/runtime/Cpp/run/usr/local/lib -lantlr4-runtime
   SHARE = /usr/local/antlr4/runtime/Cpp/run/usr/local/share
ANTLRJAR =  /home/moller/javastuff/antlr-4.8-complete.jar
JAVAOPTS = -visitor -listener

  CFLAGS = -g -std=c++11
 LDFLAGS = -std=c++11



# -listener
%.cpp : %.g4
	java  -jar $(ANTLRJAR) -Dlanguage=Cpp $(JAVAOPTS) -package mpl $<

%.o : %.cpp
	g++ -c $(CFLAGS) $(INCLUDES) -Wno-attributes $<

all::mpl

MPLLexer.h MPLLexer.cpp  : MPLLexer.g4

MPLLexer.o    : MPLLexer.cpp

MPLParser.h MPLParser.cpp : MPLParser.g4

MPLParser.o   : MPLParser.cpp

DyadicFcns.o  : DyadicFcns.cpp DyadicFcns.h SymbolTable.h Matrix.h
main.o        : main.cpp main.h MPLLexer.h MPLParser.h MPLParserBaseVisitor.h \
                SymbolTable.h DyadicFcns.h MonadicFcns.h MPLVisitor.h \
                Print.h Matrix.h Error.h
MonadicFcns.o : MonadicFcns.cpp SymbolTable.h MonadicFcns.h main.h Matrix.h
MPLVisitor.o  : MPLVisitor.cpp main.h MPLLexer.h MPLParser.h \
                MPLParserBaseVisitor.h MPLVisitor.h SymbolTable.h \
                DyadicFcns.h MonadicFcns.h Print.h
Print.o       : Print.cpp main.h Print.h Matrix.h
SymbolTable.o : SymbolTable.cpp SymbolTable.h main.h
Matrix.o      : Matrix.cpp Matrix.h
Error. o      : Error.cpp Error.h

mpl         : MPLLexer.o MPLParser.o \
              MPLVisitor.o MonadicFcns.o DyadicFcns.o SymbolTable.o \
              Print.o  Matrix.o Error.o main.o
	g++ -o $@ $(LDFLAGS) $(LIBS) $^

clean       :
	rm -f MPLLexer.interp MPLLexer.tokens MPLLexer.cpp MPLLexer.h \
              MPLParser.cpp MPLParserBaseListener.cpp MPLParserBaseVisitor.cpp \
              MPLParser.interp MPLParserListener.h MPLParserVisitor.cpp  \
              MPLParserBaseListener.h MPLParserBaseVisitor.h MPLParser.h \
              MPLParserListener.cpp MPLParser.tokens MPLParserVisitor.h \
              *.o

save:
	cvs ci -m "checkpoint"

veryclean : clean
	rm -f mpl

#VG_TEST = 1000#::6
VG_TEST = "(::20)[2 3#4 7 3 9 1 11]"

#	valgrind $(VALGRIND_FLAGS) --log-file=$@ ./mpl -e 'a = 2 3#::6; >a'  
VALGRIND_FLAGS = --gen-suppressions=all --run-libc-freeres=yes --leak-check=full
valgrind.log : FORCE
	valgrind $(VALGRIND_FLAGS) --log-file=$@ ./mpl -e 'a=1000#6'

testvg: FORCE
	./mpl -e $(VG_TEST)

FORCE:
