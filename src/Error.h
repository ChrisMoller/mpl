#include <string>

class  Error {

public:
  // fixme: generate this and strings from common source
  enum
    {
     ERROR_NONE,
     ERROR_UNKNOWN_DATA_TYPE,
     ERROR_DIMENSION_MISMATCH,
     ERROR_DIMENSION_ERROR,
     ERROR_MALFORMED_EXPRESSION,
     ERROR_MALFORMED_TRANSPOSE,
     ERROR_FAILED_TRANSPOSE,
     ERROR_OUT_OF_RANGE,
     LAST_ERROR
    };
  
  Error (int idx);
  Error (int idx, std::string cmt);

  void print ();
  void print (std::string cmt);

private:
  int error;
  std::string cmt;
};
