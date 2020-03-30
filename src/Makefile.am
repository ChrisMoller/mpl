bin_PROGRAMS = mpl
mpl_SOURCES = \
               DyadicFcns.cpp DyadicFcns.h \
               Error.cpp Error.h \
               main.cpp main.h \
               Matrix.cpp Matrix.h \
               MonadicFcns.cpp MonadicFcns.h \
               MPLVisitor.cpp MPLVisitor.h \
               Print.cpp Print.h \
               SymbolTable.cpp SymbolTable.h \
               Programme.cpp Programme.h

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


AM_CPPFLAGS = -std=c++11 -Wno-attributes $(GSL_CFLAGS) \
   -I/usr/local/antlr4/runtime/Cpp/run/usr/local/include/antlr4-runtime 

mpl_CFLAGS = -std=c++11 $(GSL_CFLAGS)

mpl_LDFLAGS = -L/usr/local/antlr4/runtime/Cpp/run/usr/local/lib \
      -lantlr4-runtime $(GSL_LIBS)

if ENABLE_GCOV
AM_CPPFLAGS += -fprofile-arcs -ftest-coverage
mpl_CFLAGS += -fprofile-arcs -ftest-coverage
mpl_LDFLAGS += -fprofile-arcs
AM_LDFLAGS = -fprofile-arcs
endif

ANTLRFLAGS = -visitor -listener

ANTLRJAR = /home/moller/javastuff/antlr-4.8-complete.jar

install-exec-hook:
	(cd ${bindir}; $(LN_S) mpl$(EXEEXT) mple$(EXEEXT))

uninstall-hook:
	(cd ${bindir}; rm -f  mple$(EXEEXT))

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

valgrind.log : mpl
	valgrind $(VALGRIND_FLAGS) --log-file=$@ \
             ./mpl  -e '(2 0 1) $ (2 3 4#1::24);' 

save:
	gitpush

clean-local: clean-local-check
.PHONY: clean-local-check
clean-local-check:
	rm -rf clc-parse.c clc-parse.h clc-lex.c *.gcno


if ENABLE_GCOV
cov-reset:
	rm -fr coverage
	find . -name "*.gcda" -exec rm {} \;
	lcov --directory . --zerocounters

cov-report:
	mkdir -p coverage
	lcov --compat-libtool --directory . --capture \
                 --output-file coverage/app.info
	lcov --remove coverage/app.info -o coverage/real.info \
           '*MPLLexer*' \
           '*MPLParser*' \
           '/usr/include/*' \
           '*antlr4-runtime*'
	genhtml --sort -o coverage/ coverage/real.info

cov:
	make cov-report
	google-chrome coverage/index.html &

#clean-local:
#	@make cov-reset

#check:
#	@make cov@find . -name "*.gcda" -exec rm {} \;
#	@lcov --directory . --zerocounters


#check:
#	@make cov

endif
