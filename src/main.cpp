#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <getopt.h>
#include <cxxabi.h>

#include <antlr4-runtime.h>
#include "main.h"
#include "MPLLexer.h"
#include "MPLParser.h"
#include "MPLParserBaseVisitor.h"

#include "SymbolTable.h"
#include "DyadicFcns.h"
#include "MonadicFcns.h"
#include "MPLVisitor.h"
#include "Print.h"

#include "Matrix.h"

using namespace mpl;
using namespace antlr4;

antlrcpp::Any::~Any ()
{
  /****

       The purpose of all this is to free things newly created
       with the pointer stored in an Any.  the "this->clear ()"
       does the actual freeing and also nulls out the pointer,
       just in case the enclosing Any gets freed more than once
       or the same pointer gets stored in multiple Any vars.
       (Both are unlikely, but you  never know...)
       
   ****/
  
  if (this->get_typeinfo() == typeid(std::vector<double>*) || 
      this->get_typeinfo() == typeid(std::vector<bool>*)   ||
      this->get_typeinfo() == typeid(std::vector<size_t>*) ||
      this->get_typeinfo() == typeid(Matrix*)	           ||
      this->get_typeinfo() == typeid(std::string*))
    this->clear ();
}

static void
do_eval (bool show, std::string str)
{
  //  std::cout << "\n\n\nParsing \"" << str << "\"" << std::endl;
  
  ANTLRInputStream input(str);
  
  MPLLexer lexer(&input);
  CommonTokenStream tokens(&lexer);
  tokens.fill();

#ifdef SHOW_TOKENS
  for (auto token : tokens.getTokens()) {
    std::cout << token->toString() << std::endl;
  }
#endif

  MPLParser parser(&tokens);
  tree::ParseTree* tree = parser.main();
  std::vector<tree::ParseTree *> children = tree->children;

#ifdef SHOW_TREE
  std::cout << tree->toStringTree(&parser) << std::endl << std::endl;
#endif

  MPLVisitor visitor (show, str);
  visitor.visit (tree);
}

int
main (int ac, char *av[])
{
  srand(static_cast<unsigned int>(clock()));
  int opt;
  
  {
    antlrcpp::Any pi = M_PI;
    antlrcpp::Any e  = M_E;
    get_global_symtab ()->insert ("pi", pi);
    get_global_symtab ()->insert ("e",  e);
  }

  bool show_exp = false;
  struct option options[] =
    {
     {"eval", required_argument, 0, 'e'},
     {"help", no_argument, 0, 'h'},
     {"show-expression", no_argument, 0, 's'},
     {0, 0, 0, 0}
    };

  while ( (opt = getopt_long(ac, av, "se:h", options, NULL)) != -1 ) { 
    switch ( opt ) {
    case 'h':
      std::cout << "Help!!\n";
      break;
    case 'e':
      do_eval (false, optarg);
      break;
    case 's':
      show_exp = true;
      break;
    case '?':  // unknown option...
      std::cerr << "Unknown option: '" << char(optopt) << "'!" << std::endl;
      break;
    }
  }

  if (optind < ac) {
    for (int i = optind; i < ac; i++) {
      //      std::cout << "reading \"" << av[i] << "\"";
      std::ifstream file(av[i]);
      if (file.is_open ()) {
	std::string str;
	 while (! file.eof() ) {
	   getline (file, str);
	   if (!str.empty ()) do_eval (show_exp, str);
	 }
	 file.close();
      }
      else std::cerr << "open failed\n";
    }
  }

  
#if 0
  //  do_eval ("a $ b + -c");
  //do_eval ("a$b+-c");
  do_eval ("4 - (11 - 7)");	// sb 0
  do_eval ("(4 - 11) - 7");	// sb -14
  do_eval ("4 - 11 - 7");
  do_eval ("4 + 11");
  do_eval ("4 - 11");
  do_eval ("4 - -11");
  do_eval ("1/4 * 11");
  do_eval ("4 / 11");
  do_eval ("4 ^ 3");
  do_eval ("-4");
  do_eval ("1/4");
  do_eval ("/4");		// error--no such op
  //do_eval ("4 $ 11 + -7");
  //do_eval ("4$11+-7");
  //  do_eval ("4j+5 - 6i * 99 ; 8+4i - 66;");
#endif
  return 0;
}
