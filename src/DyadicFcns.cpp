#include <cmath>
#include <cstring>
#include <antlr4-runtime.h>
//#include <gsl/gsl_blas.h>

#include "SymbolTable.h"
#include "DyadicFcns.h"
#include "Matrix.h"
#include "Error.h"
#include "main.h"

static double fragmentProduct (double lv, double rv){
  return lv * rv;
}
static double fragmentRatio   (double lv, double rv){return lv / rv;}
static double fragmentSum     (double lv, double rv){return lv + rv;}
static double fragmentDiff    (double lv, double rv){return lv - rv;}
static double fragmentPow     (double lv, double rv){return pow (lv, rv);}
static double fragmentLog   (double lv, double rv){return log (rv) / log (lv);}
static double fragmentRoot  (double lv, double rv){return exp (log (rv) / lv);}
static double fragmentAtan    (double lv, double rv){return atan2 (rv, lv);}
static double fragmentAtand (double lv, double rv){return RtD (atan2 (rv, lv));}
static bool   fragmentTestEq  (double lv, double rv){return (lv == rv);}
static bool   fragmentTestNE  (double lv, double rv){return (lv != rv);}
static bool   fragmentTestGE  (double lv, double rv){return (lv >= rv);}
static bool   fragmentTestLE  (double lv, double rv){return (lv <= rv);}
static bool   fragmentTestGT  (double lv, double rv){return (lv >  rv);}
static bool   fragmentTestLT  (double lv, double rv){return (lv <  rv);}

static void
do_dyadic  (antlrcpp::Any &rc, dyadic_op op,
	    antlrcpp::Any &left, antlrcpp::Any &right)
{
  if (right.get_typeinfo() == typeid(double) &&
      left.get_typeinfo()  == typeid(double)) {
    double lv = left.as<double>();
    double rv = right.as<double>();
    double res = (*op)(lv, rv);
    rc = res;
  }
  else if (right.get_typeinfo() == typeid(Matrix*) &&
	   left.get_typeinfo()  == typeid(double)) {
    double lv  = left.as<double>();
    Matrix *rv = right.as<Matrix *>();
    Matrix *mtx = new Matrix (op, lv, rv);
    rc = mtx;
  }
  else if (left.get_typeinfo() == typeid(Matrix*) &&
	   right.get_typeinfo()  == typeid(double)) {
    double rv  = right.as<double>();
    Matrix *lv = left.as<Matrix *>();
    Matrix *mtx = new Matrix (op, lv, rv);
    rc = mtx;
  }
  else if (left.get_typeinfo() == typeid(Matrix*) &&
	   right.get_typeinfo() == typeid(Matrix*)) {
    Matrix *lv = left.as<Matrix *>();
    Matrix *rv = right.as<Matrix *>();
    if (lv->isomorphic (rv)) {
      Matrix *mtx = new Matrix (op, lv, rv);
      rc = mtx;
    }
    else {
      rc = Error(Error::ERROR_DIMENSION_MISMATCH);
    }
  }
  else if (right.get_typeinfo() == typeid(std::vector<double>*) &&
	   left.get_typeinfo()  == typeid(double)) {
    std::vector<double> *right_vec = right.as<std::vector<double>*>();
    double lv = left.as<double>();
    size_t n = right_vec->size ();
    std::vector<double> *res = new std::vector<double>;
    res->resize (n);
    for (size_t i = 0;  i < n; i++)
      (*res)[i] = (*op)(lv, (*right_vec)[i]);
    rc = res;
  }
  else if (left.get_typeinfo() == typeid(std::vector<double>*) &&
	   right.get_typeinfo()  == typeid(double)) {
    std::vector<double> *left_vec = left.as<std::vector<double>*>();
    double rv = right.as<double>();
    size_t n = left_vec->size ();
    std::vector<double> *res = new std::vector<double>;
    res->resize (n);
    for (size_t i = 0;  i < n; i++)
      (*res)[i] = (*op)((*left_vec)[i], rv);
    rc = res;
  }
  else if (left.get_typeinfo() == typeid(std::vector<double>*) &&
	   right.get_typeinfo()  == typeid(std::vector<double>*)) {
    std::vector<double> *left_vec  = left.as<std::vector<double>*>();
    std::vector<double> *right_vec = right.as<std::vector<double>*>();
    size_t n = left_vec->size ();
    size_t m = right_vec->size ();
    std::vector<double> *res = new std::vector<double>;
    res->resize (n);
    if (n == m) {
      for (size_t i = 0;  i < n; i++) 
	(*res)[i] = (*op)((*left_vec)[i], (*right_vec)[i]);
      rc = res;
    }
    else {
      delete res;
      rc = Error(Error::ERROR_DIMENSION_MISMATCH,
		 ": Vector length mismatch");
    }
  }
  else {
    rc = Error(Error::ERROR_UNKNOWN_DATA_TYPE, ", dyadic operation."); 
  }
}

static void
dyadicStar (antlrcpp::Any &rc, antlrcpp::Any &left,
	    antlrcpp::Any &right, antlrcpp::Any &qual)
{
  do_dyadic (rc, fragmentProduct, left, right);
}

static void
dyadicSlash (antlrcpp::Any &rc, antlrcpp::Any &left,
	     antlrcpp::Any &right, antlrcpp::Any &qual)
{
  do_dyadic (rc, fragmentRatio, left, right);
}

static void
dyadicPlus  (antlrcpp::Any &rc, antlrcpp::Any &left,
	     antlrcpp::Any &right, antlrcpp::Any &qual)
{
  do_dyadic (rc, fragmentSum, left, right);
}

static void
dyadicMinus (antlrcpp::Any &rc, antlrcpp::Any &left,
	     antlrcpp::Any &right, antlrcpp::Any &qual)
{
   do_dyadic (rc, fragmentDiff, left, right);
}

static void
dyadicHat   (antlrcpp::Any &rc, antlrcpp::Any &left,
	     antlrcpp::Any &right, antlrcpp::Any &qual)
{
   do_dyadic (rc, fragmentPow, left, right);
}

static void
dyadicLog   (antlrcpp::Any &rc, antlrcpp::Any &left,
	     antlrcpp::Any &right, antlrcpp::Any &qual)
{
   do_dyadic (rc, fragmentLog, left, right);
}

static void
dyadicRoot  (antlrcpp::Any &rc, antlrcpp::Any &left,
	     antlrcpp::Any &right, antlrcpp::Any &qual)
{
   do_dyadic (rc, fragmentRoot, left, right);
}

static void
dyadicAtan  (antlrcpp::Any &rc, antlrcpp::Any &left,
	     antlrcpp::Any &right, antlrcpp::Any &qual)
{
   do_dyadic (rc, fragmentAtan, left, right);
}

static void
dyadicAtand (antlrcpp::Any &rc, antlrcpp::Any &left,
	     antlrcpp::Any &right, antlrcpp::Any &qual)
{
   do_dyadic (rc, fragmentAtand, left, right);
}

static void
dyadicEqual (antlrcpp::Any &rc, antlrcpp::Any &left,
	     antlrcpp::Any &right, antlrcpp::Any &qual)
{
  if (left.get_typeinfo() == typeid(std::string *)) {
    std::string *sym = left.as<std::string *>();
    insert_symbol_value (*sym, right);
    //    get_global_symtab ()->insert (*sym, right);
  }
  else {
    std::cout << "Lvalue can't be a constant.\n";
    rc = Error (Error::ERROR_CONSTANT_LVALUE);
  }
}

static void
dyadicRange (antlrcpp::Any &rc, antlrcpp::Any &left,
	     antlrcpp::Any &right, antlrcpp::Any &qual)
{
  if (right.get_typeinfo() == typeid(double) &&
      left.get_typeinfo()  == typeid(double)) {
    double lv = left.as<double>();
    double rv = right.as<double>();
    std::vector<double> *array = new std::vector<double>;
    int b = int (lv);
    int n = int (rv);
    int incr = (n >= b) ? 1 : -1;
    n += incr;

    for (int i = b; i != n; i+=incr)
      array->push_back (double(i));
    rc = array;
  }
  else {
    rc = Error(Error::ERROR_UNKNOWN_DATA_TYPE, ", range operation.");
  }
}

static void
do_test (antlrcpp::Any &rc, dyadic_test op,
	 antlrcpp::Any &left, antlrcpp::Any &right)
{
  /*
    fixme array of t/f, or one massive t/f
    maybe scan array of t/f
  */

  if (right.get_typeinfo() == typeid(double) &&
      left.get_typeinfo()  == typeid(double)) {
    double lv = left.as<double>();
    double rv = right.as<double>();
    bool res =  (*op)(lv, rv);
    rc = res;
  }
  else if (right.get_typeinfo() == typeid(std::vector<double>*) &&
	   left.get_typeinfo()  == typeid(double)) {
    std::vector<double> *right_vec = right.as<std::vector<double>*>();
    double lv = left.as<double>();
    size_t n = right_vec->size ();
    std::vector<bool> *array = new std::vector<bool>(n);
    for (size_t i = 0;  i < n; i++)
      (*array)[i] = (*op)(lv, (*right_vec)[i]);
    rc = array;
  }
  else if (left.get_typeinfo() == typeid(std::vector<double>*) &&
	   right.get_typeinfo()  == typeid(double)) {
    std::vector<double> *left_vec = left.as<std::vector<double>*>();
    double rv = right.as<double>();
    size_t n = left_vec->size ();
    std::vector<bool> *array = new std::vector<bool>(n);
    for (size_t i = 0;  i < n; i++)
      (*array)[i] = (*op)((*left_vec)[i], rv);
    rc = array;
  }
  else if (left.get_typeinfo() == typeid(std::vector<double>*) &&
	   right.get_typeinfo()  == typeid(std::vector<double>*)) {
    std::vector<double> *left_vec  = left.as<std::vector<double>*>();
    std::vector<double> *right_vec = right.as<std::vector<double>*>();
    size_t n = left_vec->size ();
    size_t m = right_vec->size ();
    std::vector<bool> *array = new std::vector<bool>(n);
    if (n == m) {
      for (size_t i = 0;  i < n; i++) 
	(*array)[i] = (*op)((*left_vec)[i], (*right_vec)[i]);
    }
    else {
      rc = Error(Error::ERROR_DIMENSION_MISMATCH,
		 ": Vector length mismatch");
    }
    rc = array;
  }
  else {
    rc = Error(Error::ERROR_UNKNOWN_DATA_TYPE, ", comparisin operation.");
  }
}

static void
dyadicTestEq (antlrcpp::Any &rc, antlrcpp::Any &left,
	      antlrcpp::Any &right, antlrcpp::Any &qual)
{
  do_test (rc, fragmentTestEq, left, right);
}

static void
dyadicTestLT (antlrcpp::Any &rc, antlrcpp::Any &left,
	      antlrcpp::Any &right, antlrcpp::Any &qual)
{
  do_test (rc, fragmentTestLT, left, right);
}

static void
dyadicTestGT (antlrcpp::Any &rc, antlrcpp::Any &left,
	      antlrcpp::Any &right, antlrcpp::Any &qual)
{
   do_test (rc, fragmentTestGT, left, right);
}

static void
dyadicTestLE (antlrcpp::Any &rc, antlrcpp::Any &left,
	      antlrcpp::Any &right, antlrcpp::Any &qual)
{
  return do_test (rc, fragmentTestLE, left, right);
}

static void
dyadicTestGE (antlrcpp::Any &rc, antlrcpp::Any &left,
	      antlrcpp::Any &right, antlrcpp::Any &qual)
{
  do_test (rc, fragmentTestGE, left, right);
}

static void
dyadicTestNE (antlrcpp::Any &rc, antlrcpp::Any &left,
	      antlrcpp::Any &right, antlrcpp::Any &qual)
{
  do_test (rc, fragmentTestNE, left, right);
}

static void
dyadicShape (antlrcpp::Any &rc, antlrcpp::Any &left,
	     antlrcpp::Any &right, antlrcpp::Any &qual)
{
  /*
    fixme array of t/f, or one massive t/f
    maybe scan array of t/f
  */

  if (right.get_typeinfo() == typeid(double) &&
      left.get_typeinfo()  == typeid(double)) {
    double lv = left.as<double>();
    double rv = right.as<double>();
    int n = int(lv);
    std::vector<double> *array = new std::vector<double>(n);
    if (n > 0) {
      for (size_t i = 0;  i < n; i++)
	(*array)[i] = rv;
    }
    rc = array;
  }
  else if (right.get_typeinfo() == typeid(std::vector<double>*) &&
	   left.get_typeinfo()  == typeid(double)) {
    std::vector<double> *right_vec = right.as<std::vector<double>*>();
    double lv = left.as<double>();
    int n = int(lv);
    size_t m = right_vec->size ();
    std::vector<double> *array = new std::vector<double>(n);
    if (n > 0 && m > 0) {
      for (size_t i = 0;  i < n; i++)
	(*array)[i] = (*right_vec)[i%m];
    }
    else {
      delete array;
      rc = Error(Error::ERROR_DIMENSION_ERROR,
		 ": Can't create a zero-dimension vector");
    }
    rc = array;
  }
  else if (left.get_typeinfo() == typeid(std::vector<double>*) &&
	   right.get_typeinfo()  == typeid(double)) {
    std::vector<double> *left_vec = left.as<std::vector<double>*>();
    if (left_vec->size () > 0) {
      Matrix *mtx = new Matrix (left);
      double rv = right.as<double>();
      mtx->fill (rv);
      rc = mtx;
    }
    else {
      rc = Error(Error::ERROR_DIMENSION_ERROR,
		 ": Can't create a zero-dimension matrix");
    }
  }
  else if (left.get_typeinfo() == typeid(std::vector<double>*) &&
	   right.get_typeinfo()  == typeid(std::vector<double>*)) {
    std::vector<double> *left_vec  = left.as<std::vector<double>*>();
    std::vector<double> *right_vec = right.as<std::vector<double>*>();
    if (left_vec->size () > 0 && right_vec->size () > 0) {
      Matrix *mtx = new Matrix (left);
      mtx->copy (right);
      rc = mtx;
    }
    else {
      rc = Error(Error::ERROR_DIMENSION_MISMATCH,
		 ": Matrix dimension mismatch");
    }
  }
  else if (left.get_typeinfo() == typeid(std::vector<double>*) &&
	   right.get_typeinfo()  == typeid(Matrix *)) {
    std::vector<double> *left_vec  = left.as<std::vector<double>*>();
    size_t cnt = 1;
    for (size_t i = 0; i < left_vec->size (); i++)
      cnt *= (*left_vec)[i];
    Matrix *right_mtx = right.as<Matrix*>();
    if (cnt == right_mtx->size ()) {
      Matrix *mtx = new Matrix (left);
      mtx->reshape (right);
      rc = mtx;
    }
    else {
      rc = Error(Error::ERROR_DIMENSION_MISMATCH,
		 ": Matrix dimension mismatch");
    }
  }
  else if (left.get_typeinfo() == typeid(double) &&
	   right.get_typeinfo()  == typeid(Matrix *)) {
    double cnt  = left.as<double>();
    Matrix *right_mtx = right.as<Matrix*>();
    if (cnt < 0.0 || cnt == right_mtx->size ()) {
      std::vector<double>*vec =
	new std::vector<double>(right_mtx->size ());
      std::memmove (vec->data (), right_mtx->get_data ()->data (),
		    right_mtx->size () * sizeof(double));  
      rc = vec;
    }
    else {
      rc = Error(Error::ERROR_DIMENSION_MISMATCH,
		 ": Matrix dimension mismatch");
    }
  }
  else {
    // errors
  }
}

static void
dyadicMatMult (antlrcpp::Any &rc, antlrcpp::Any &left,
	       antlrcpp::Any &right, antlrcpp::Any &qual)
{
  if (left.get_typeinfo() == typeid(Matrix*) &&
      right.get_typeinfo()  == typeid(Matrix*)) {
    Matrix *lv = left.as<Matrix *>();
    Matrix *rv = right.as<Matrix *>();
    Matrix *mtx = lv->multiply (rv);
    if (mtx) rc = mtx;
    else {
      rc = Error(Error::ERROR_FAILED_MTX_MULTIPLY);
      std::string em = rv->get_errmsg ();
      std::cout << em << std::endl;
    }
  }
  else if (left.get_typeinfo() == typeid(Matrix*) &&
      right.get_typeinfo()  == typeid(std::vector<double>*)) {
    Matrix *lv = left.as<Matrix *>();
    std::vector<double> *rv = right.as<std::vector<double> *>();
    std::vector<double> *vec = lv->multiply (Matrix::VM_VEC_RIGHT, rv);
    if (vec) rc = vec;
    else {
      rc = Error(Error::ERROR_FAILED_MTX_MULTIPLY);
      std::string em = lv->get_errmsg ();
      std::cout << em << std::endl;
    }
  }
  else if (right.get_typeinfo() == typeid(Matrix*) &&
      left.get_typeinfo()  == typeid(std::vector<double>*)) {
    Matrix *rv = right.as<Matrix *>();
    std::vector<double> *lv = left.as<std::vector<double> *>();
    std::vector<double> *vec = rv->multiply (Matrix::VM_VEC_LEFT, lv);
    if (vec) rc = vec;
    else {
      rc = Error(Error::ERROR_FAILED_MTX_MULTIPLY);
      std::string em = rv->get_errmsg ();
      std::cout << em << std::endl;
    }
  }
  else {
    rc = Error(Error::ERROR_UNKNOWN_DATA_TYPE, "Matrix multipy");
  }
}

static void
dyadicMatSolve (antlrcpp::Any &rc, antlrcpp::Any &left,
		antlrcpp::Any &right, antlrcpp::Any &qual)
{
  if (left.get_typeinfo() == typeid(Matrix*) &&
      right.get_typeinfo()  == typeid(Matrix*)) {
    Matrix *lv = left.as<Matrix *>();
    Matrix *rv = right.as<Matrix *>();
    Matrix *mtx = lv->solve (rv);
    if (mtx) rc = mtx;
    else {
      rc = Error(Error::ERROR_FAILED_MTX_MULTIPLY);
      std::string em = rv->get_errmsg ();
      std::cout << em << std::endl;
    }
  }
  else if (left.get_typeinfo() == typeid(Matrix*) &&
	   right.get_typeinfo()  == typeid(std::vector<double>*)) {
    Matrix *lv = left.as<Matrix *>();
    std::vector<double> *rv = right.as<std::vector<double> *>();
    std::vector<double> *vec = lv->solve (Matrix::VM_VEC_LEFT, rv);
    if (vec) rc = vec;
    else {
      rc = Error(Error::ERROR_FAILED_MTX_MULTIPLY);
      std::string em = lv->get_errmsg ();
      std::cout << em << std::endl;
    }
  }
  else if (right.get_typeinfo() == typeid(Matrix*) &&
	   left.get_typeinfo()  == typeid(std::vector<double>*)) {
    Matrix *rv = right.as<Matrix *>();
    std::vector<double> *lv = left.as<std::vector<double> *>();
    std::vector<double> *vec = rv->solve (Matrix::VM_VEC_RIGHT, lv);
    if (vec) rc = vec;
    else {
      rc = Error(Error::ERROR_FAILED_MTX_MULTIPLY);
      std::string em = rv->get_errmsg ();
      std::cout << em << std::endl;
    }
  }
  else {
    rc = Error(Error::ERROR_UNKNOWN_DATA_TYPE, "Matrix solve");
  }
}

static void
dyadicIdentity (antlrcpp::Any &rc, antlrcpp::Any &left,
		antlrcpp::Any &right, antlrcpp::Any &qual)
{
  if (right.get_typeinfo() == typeid(double) &&
      left.get_typeinfo()  == typeid(double)) {
    double rank = left.as<double>();
    double dim = right.as<double>();
    if (rank >= 2.0 && rank <= dim) {
      Matrix *rv = new Matrix (static_cast<size_t>(rank),
			       static_cast<size_t>(dim));
      rc = rv;
    }
    else {
      rc = Error(Error::ERROR_OUT_OF_RANGE, ", dyadic identity");
    }
  }
  else {
    rc = Error(Error::ERROR_UNKNOWN_DATA_TYPE, ", dyadic identity");
  }
}

void
do_vector_shift (antlrcpp::Any &rc, double shift,
		 antlrcpp::Any &right, antlrcpp::Any &qual)
{
  std::vector<double> *right_vec = right.as<std::vector<double>*>();
  std::vector<double> *rv = new std::vector<double> (right_vec->size ());
  shift = fmod (shift, static_cast<double>(right_vec->size ()));
  if (shift < 0.0) shift += static_cast<double>(right_vec->size ());
  int lv = static_cast<int>(shift);
  double *fm = right_vec->data ();
  double *to = rv->data ();
  std::memmove (to, fm + lv, (right_vec->size () - lv) * sizeof(double));
  std::memmove (to + (right_vec->size () - lv) , fm, lv * sizeof(double));
  rc = rv;
}

static void
dyadicLeft (antlrcpp::Any &rc, antlrcpp::Any &left,
	    antlrcpp::Any &right, antlrcpp::Any &qual)
{
  if (left.get_typeinfo() == typeid(double)) {
    double shift = left.as<double>();
    if (right.get_typeinfo() == typeid(double)) {
      double rv = right.as<double>();
      rc = rv;
    }
    else if (right.get_typeinfo()  == typeid(std::vector<double>*)) {
      do_vector_shift (rc, shift, right, qual);
    }
    else if (right.get_typeinfo() == typeid(Matrix*)) {
      Matrix *rv = right.as<Matrix *>();
      Matrix *mtx =rv->shift (shift, qual);
      if (mtx) rc = mtx;
      else {
	rc = Error(Error::ERROR_DIMENSION_MISMATCH, ", axes");
	std::string em = mtx->get_errmsg ();
	std::cout << em << std::endl;
      }
    }
    else {
      rc = Error(Error::ERROR_UNKNOWN_DATA_TYPE, ", dyadic left operation");
    }
  }
  else {
    rc = Error(Error::ERROR_UNKNOWN_DATA_TYPE, ", dyadic left operation");
  }
}

static void
dyadicRight (antlrcpp::Any &rc, antlrcpp::Any &left,
	     antlrcpp::Any &right, antlrcpp::Any &qual)
{
  if (left.get_typeinfo() == typeid(double)) {
    double shift = left.as<double>();
    if (right.get_typeinfo() == typeid(double)) {
      double rv = right.as<double>();
      rc = rv;
    }
    else if (right.get_typeinfo()  == typeid(std::vector<double>*)) {
      do_vector_shift (rc, -shift, right, qual);
    }
    else if (right.get_typeinfo() == typeid(Matrix*)) {
    }
    else {
      rc = Error(Error::ERROR_UNKNOWN_DATA_TYPE, ", dyadic left operation");
    }
  }
  else {
    rc = Error(Error::ERROR_UNKNOWN_DATA_TYPE, ", dyadic left operation");
  }
}

static dfunc dfuncs[] =
{
 nullptr,	// empty	 0
 nullptr,	// DUMMY	 1
 nullptr,	// Number	 2
 dyadicStar,	// OpStar	 3
 dyadicSlash,	// OpSlash	 4
 dyadicPlus,	// OpPlus	 5
 dyadicMinus,	// OpMinus	 6
 nullptr,	// OpDollar	 7
 nullptr,	// OpPercent	 8
 dyadicHat,	// OpHat	 9
 nullptr,	// OpLn		10
 dyadicLog,	// OpLog	11
 nullptr,	// OpExp	12
 dyadicRoot,	// OpRoot	13
 nullptr,	// OpSin	14
 nullptr,	// OpCos	15
 nullptr,	// OpTan	16
 nullptr,	// OpAsin	17
 nullptr,	// OpAcos	18
 dyadicAtan,	// OpAtan	19
 nullptr,	// OpSind	20
 nullptr,	// OpCosd	21
 nullptr,	// OpTand	22
 nullptr,	// OpAsind	23
 nullptr,	// OpAcosd	24
 dyadicAtand,	// OpAtand	25
 nullptr,	// OpBar	26
 nullptr,	// OpBang	27
 nullptr,	// OpCeil	28
 nullptr,	// OpFloor	29
 nullptr,	// OpRound	30
 nullptr,	// OpLeftAngle	31
 nullptr,	// OpRightAngle	32
 nullptr,	// OpQMark	33
 dyadicEqual,	// OpEqual	34
 dyadicShape,	// OpPound	35
 nullptr,	// OpComma	36
 nullptr,	// OpColon	37
 dyadicMatMult,	// OpBSStar	38
 dyadicRange,	// OpColonColon	39
 dyadicTestEq,	// OpQEqual		40
 dyadicTestLT,	// OpQLeftAngle		41
 dyadicTestGT,	// OpQRightAngle	42
 dyadicTestLE,	// OpQLeftAngleEqual	43
 dyadicTestGE,	// OpQRightAngleEqual	44
 dyadicTestNE,	// OpQBangEqual		45
 dyadicTestGE,	// OpQBangLeftAngle		46
 dyadicTestLE,	// OpQBangRightAngle		47
 dyadicTestGT,	// OpQBangLeftAngleEqual	48
 dyadicTestLT,	// OpQBangRightAngleEqual	49
 nullptr,	// OpDet			50
 dyadicMatSolve,	// OpBSSlash		51
 nullptr,		// OpSlashPlus		52
 nullptr,		// OpSlashStar		53
 dyadicIdentity,	// OpBSI		54
 dyadicLeft,		// OpLALA		55
 dyadicRight,		// OpRARA		56
};

dfunc
dyadic_get_func (int idx)
{
  dfunc res = nullptr;
  if (idx >= 0 && idx < sizeof(dfuncs) / sizeof(dfunc)) {
    res = dfuncs[idx];
  }

  return res;
}
