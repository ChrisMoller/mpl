#include <vector>
#include <iostream>
#include <string>
#include "Error.h"

// fixme: generate this and enum from common source
static std::vector<std::string> descriptions =
  {
   "No error",			// ERROR_NONE
   "Unknown data type",		// ERROR_UNKNOWN_DATA_TYPE
   "Dimension mismatch",	// ERROR_DIMENSION_MISMATCH
   "Dimension error",		// ERROR_DIMENSION_ERROR
   "Mal-formed expression",	// ERROR_MALFORMED_EXPRESSION
   "Mal-formed transpose",	// ERROR_MALFORMED_TRANSPOSE
   "Failed transpose",		// ERROR_FAILED_TRANSPOSE
   "Index out of range",	// ERROR_OUT_OF_RANGE
   "Lvalue can't be a constant." // ERROR_CONSTANT_LVALUE
   "Internal error",		// ERROR_INTERNAL
  };

Error::Error (int idx)
{
  error = idx;
}

Error::Error (int idx, std::string str)
{
  error = idx;
  cmt = str;
}

void
Error::print ()
{
  if (error > 0 && error < LAST_ERROR)
    std::cout << descriptions[error] << cmt
	      << std::endl;
}

void
Error::print (std::string cmt)
{
  if (error > 0 && error < LAST_ERROR)
    std::cout << descriptions[error] << cmt
	      << std::endl;
}

