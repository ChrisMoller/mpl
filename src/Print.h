#ifndef PRINT_H
#define PRINT_H

#include <iostream>

#include <antlr4-runtime.h>

void print_str (bool show, std::string lexpr, std::string str);
void print_val (bool show, std::string lexpr, antlrcpp::Any &res);
#endif // PRINT_H
