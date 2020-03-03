#include <iostream>

#include <antlr4-runtime.h>
#include "main.h"
#include "Matrix.h"
#include "Error.h"

void
print_str (bool show, std::string lexpr, std::string str)
{
  if (isTestMode ()) {
    std::cout << str << std::endl;
  }
  else {
    if (show) std::cout << lexpr << " ";
    std::cout << ">> " << str << std::endl;
  }
}

void
print_val (bool show, std::string lexpr, antlrcpp::Any &res)
{
  if (isTestMode ()) {
    if (res.get_typeinfo () == typeid(double)) {
      double res_val = res.as<double>();
      std::cout << res_val << std::endl;
    }
    else  if (res.get_typeinfo() == typeid(std::string)) {
      std::string str = res.as<std::string>();
      print_str (show, lexpr, str);
    }
    else  if (res.get_typeinfo() == typeid(bool)) {
      bool test = res.as<bool>();
      print_str (show, lexpr, test ? "true" : "false");
    }
    else  if (res.get_typeinfo() == typeid(std::vector<int>*)) {
      std::vector<int> *vec = res.as<std::vector<int>*>();
      for (auto dp = vec->begin (); dp != vec->end (); dp++)
	std::cout << *dp << " ";
      std::cout << std::endl;
    }
    else  if (res.get_typeinfo() == typeid(std::vector<double>*)) {
      std::vector<double> *vec = res.as<std::vector<double>*>();
      for (auto dp = vec->begin (); dp != vec->end (); dp++)
	std::cout << *dp << " ";
      std::cout << std::endl;
    }
    else  if (res.get_typeinfo() == typeid(std::vector<bool>*)) {
      std::vector<bool> *vec = res.as<std::vector<bool>*>();
      for (auto dp = vec->begin (); dp != vec->end (); dp++)
	std::cout << ((*dp) ? "true" : "false")  << " ";
      std::cout << std::endl;
    }
    else  if (res.get_typeinfo() == typeid(Matrix*)) {
      Matrix *mtx = res.as<Matrix *>();
      mtx->print ();
    }
    else  if (res.get_typeinfo() == typeid(Error)) {
      Error err = res.as<Error>();
      err.print ();
      //    err.print (": dimensions of the operands must match");
    }
    else  if (res.isNull ()) {//   if (res.get_typeinfo() == typeid(nullptr)) {
      // do nothing
    }	
    else {
      std::cout << "Unknown type in print fcn "
		<< "(" << demangle (res) << ")"
		<< std::endl;
    }
  }
  else {
    if (show) std::cout << lexpr;
    if (res.get_typeinfo () == typeid(double)) {
      double res_val = res.as<double>();
      if (show) std::cout << lexpr << "\n\t";
      std::cout << ">> " << res_val << std::endl;
    }
    else  if (res.get_typeinfo() == typeid(std::string)) {
      std::string str = res.as<std::string>();
      print_str (show, lexpr, str);
    }
    else  if (res.get_typeinfo() == typeid(bool)) {
      bool test = res.as<bool>();
      print_str (show, lexpr, test ? "true" : "false");
    }
    else  if (res.get_typeinfo() == typeid(std::vector<int>*)) {
      std::vector<int> *vec = res.as<std::vector<int>*>();
      std::cout << ">> ";
      for (auto dp = vec->begin (); dp != vec->end (); dp++)
	std::cout << *dp << " ";
      std::cout << std::endl;
    }
    else  if (res.get_typeinfo() == typeid(std::vector<double>*)) {
      std::vector<double> *vec = res.as<std::vector<double>*>();
      std::cout << ">> ";
      for (auto dp = vec->begin (); dp != vec->end (); dp++)
	std::cout << *dp << " ";
      std::cout << std::endl;
    }
    else  if (res.get_typeinfo() == typeid(std::vector<bool>*)) {
      std::vector<bool> *vec = res.as<std::vector<bool>*>();
      std::cout << ">> ";
      for (auto dp = vec->begin (); dp != vec->end (); dp++)
	std::cout << ((*dp) ? "true" : "false")  << " ";
      std::cout << std::endl;
    }
    else  if (res.get_typeinfo() == typeid(Matrix*)) {
      Matrix *mtx = res.as<Matrix *>();
      mtx->print ();
    }
    else  if (res.get_typeinfo() == typeid(Error)) {
      Error err = res.as<Error>();
      err.print ();
      //    err.print (": dimensions of the operands must match");
    }
    else  if (res.isNull ()) {//   if (res.get_typeinfo() == typeid(nullptr)) {
      // do nothing
    }	
    else {
      std::cout << "stop here\n";
      std::cout << "Unknown type in print fcn "
		<< "(" << demangle (res) << ")"
		<< std::endl;
    }
  }
}
