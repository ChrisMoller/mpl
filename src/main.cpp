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

static source_e current_source;

source_e get_source () { return current_source; }
bool isFromFile ()     { return current_source == SOURCE_FILE; }
bool isFromCmdLine ()  { return current_source == SOURCE_CMDLINE; }
void set_source (source_e src) {current_source = src;}

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
do_eval (source_e src, bool show, std::string str)
{
  //  std::cout << "\n\n\nParsing \"" << str << "\"" << std::endl;

  set_source (src);
  
  ANTLRInputStream input(str);
  
  MPLLexer lexer(&input);
  CommonTokenStream tokens(&lexer);
  tokens.fill();

#if 0
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
  set_source (SOURCE_NONE);
}


static int test_mode = 0;

bool isTestMode () {bool rc = test_mode ? true : false; return rc; }

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
     {"test-mode", no_argument, &test_mode, 1},
     {0, 0, 0, 0}
    };

  while ( (opt = getopt_long(ac, av, "se:h", options, NULL)) != -1 ) { 
    switch ( opt ) {
    case 'h':
      std::cout << "Help!!\n";
      break;
    case 'e':
      {
	if (optarg[0] == '\'' &&
	    optarg[strlen (optarg) -1] == '\'') {
	  optarg[strlen (optarg) -1] = 0;
	  do_eval (SOURCE_CMDLINE, false, &optarg[1]);
	}
        else do_eval (SOURCE_CMDLINE, false, optarg);
      }
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
	std::stringstream buffer;
	buffer << file.rdbuf();
	std::string str = buffer.str ();
	if (!str.empty ()) {
	  do_eval (SOURCE_FILE, show_exp, str);
	}
      }
      else std::cerr << "open failed\n";


#if 0
      if (file.is_open ()) {
	std::string str;
	int i = 0;
	 while (! file.eof() ) {
	   getline (file, str);
	   if (!str.empty ()) {
	     std::cout << "line " << i++ << ": " << str <<  std::endl;
	     do_eval (SOURCE_FILE, show_exp, str);
	   }
	 }
	 file.close();
      }
      else std::cerr << "open failed\n";
#endif
    }
  }

  
  return 0;
}
