#include <cmath>
#include <cstdlib>


#include "Error.h"

//https://stackoverflow.com/questions/686353/random-float-number-generation

#include <antlr4-runtime.h>
#include "SymbolTable.h"
#include "MonadicFcns.h"
#include "Matrix.h"
#include "main.h"

#define DtR(d) ((M_PI / 180.0) * (d))
#define RtD(r) ((180.0 / M_PI) * (r))

static double fragmentInvert (double rv){return 1.0 / rv;}
static double fragmentNegate (double rv){return -rv;}
static double fragmentLn     (double rv){return log (rv);}
static double fragmentLog10  (double rv){return log10 (rv);}
static double fragmentExp    (double rv){return exp (rv);}
static double fragmentRoot   (double rv){return sqrt (rv);}
static double fragmentSin    (double rv){return sin (rv);}
static double fragmentCos    (double rv){return cos (rv);}
static double fragmentTan    (double rv){return tan (rv);}
static double fragmentAsin   (double rv){return asin (rv);}
static double fragmentAcos   (double rv){return acos (rv);}
static double fragmentAtan   (double rv){return atan (rv);}
static double fragmentSind   (double rv){return sin (DtR (rv));}
static double fragmentCosd   (double rv){return cos (DtR (rv));}
static double fragmentTand   (double rv){return tan (DtR (rv));}
static double fragmentAsind  (double rv){return RtD (asin (rv));}
static double fragmentAcosd  (double rv){return RtD (acos (rv));}
static double fragmentAtand  (double rv){return RtD (atan (rv));}
static double fragmentBar    (double rv){return fabs (rv);}
static double fragmentBang   (double rv){return tgamma (1.0 + rv);}
static double fragmentCeil   (double rv){return ceil (rv);}
static double fragmentFloor  (double rv){return floor (rv);}
static double fragmentRound  (double rv){return round (rv);}
static double fragmentRand   (double rv)
{
  return rv * double(rand()) / (double(RAND_MAX) + 1.0);
}

static void
do_monadic (antlrcpp::Any &rc, monadic_op op, antlrcpp::Any &right)
{
  if (right.get_typeinfo() == typeid(double)) {
    double rv = right.as<double>();
    double res = (*op)(rv);
    rc = res;
  }
  else if (right.get_typeinfo()  == typeid(std::vector<double>*)) {
    std::vector<double> *right_vec = right.as<std::vector<double>*>();
    size_t n = right_vec->size ();
    std::vector<double> *res = new std::vector<double>;
    res->resize (n);
    for (size_t i = 0;  i < n; i++)
      (*res)[i] = (*op)((*right_vec)[i]);
    rc = res;
  }
  else if (right.get_typeinfo() == typeid(Matrix*)) {
    Matrix *rv = right.as<Matrix *>();
    Matrix *mtx = new Matrix (op, rv);
    rc = mtx;
  }
  else {
    rc = Error(Error::ERROR_UNKNOWN_DATA_TYPE, ", monadic operation");
  }
}

static void
monadicSlash (antlrcpp::Any &rc, antlrcpp::Any &right)
{
  do_monadic (rc, fragmentInvert, right);
}

static void
monadicMinus (antlrcpp::Any &rc, antlrcpp::Any &right)
{
  do_monadic (rc, fragmentNegate, right);
}

static void
monadicLn (antlrcpp::Any &rc, antlrcpp::Any &right)
{
  do_monadic (rc, fragmentLn, right);
}

static void
monadicLog (antlrcpp::Any &rc, antlrcpp::Any &right)
{
  do_monadic (rc, fragmentLog10, right);
}

static void
monadicExp (antlrcpp::Any &rc, antlrcpp::Any &right)
{
  do_monadic (rc, fragmentExp, right);
}

static void
monadicRoot (antlrcpp::Any &rc, antlrcpp::Any &right)
{
  do_monadic (rc, fragmentRoot, right);
}

static void
monadicSin (antlrcpp::Any &rc, antlrcpp::Any &right)
{
  do_monadic (rc, fragmentSin, right);
}

static void
monadicCos (antlrcpp::Any &rc, antlrcpp::Any &right)
{
  do_monadic (rc, fragmentCos, right);
}

static void
monadicTan (antlrcpp::Any &rc, antlrcpp::Any &right)
{
  do_monadic (rc, fragmentTan, right);
}

static void
monadicAsin (antlrcpp::Any &rc, antlrcpp::Any &right)
{
  do_monadic (rc, fragmentAsin, right);
}

static void
monadicAcos (antlrcpp::Any &rc, antlrcpp::Any &right)
{
  do_monadic (rc, fragmentAcos, right);
}

static void
monadicAtan (antlrcpp::Any &rc, antlrcpp::Any &right)
{
  do_monadic (rc, fragmentAtan, right);
}

static void
monadicSind (antlrcpp::Any &rc, antlrcpp::Any &right)
{
  do_monadic (rc, fragmentSind, right);
}

static void
monadicCosd (antlrcpp::Any &rc, antlrcpp::Any &right)
{
  do_monadic (rc, fragmentCosd, right);
}

static void
monadicTand (antlrcpp::Any &rc, antlrcpp::Any &right)
{
  do_monadic (rc, fragmentTand, right);
}

static void
monadicAsind (antlrcpp::Any &rc, antlrcpp::Any &right)
{
  do_monadic (rc, fragmentAsind, right);
}

static void
monadicAcosd (antlrcpp::Any &rc, antlrcpp::Any &right)
{
  do_monadic (rc, fragmentAcosd, right);
}

static void
monadicAtand (antlrcpp::Any &rc, antlrcpp::Any &right)
{
  do_monadic (rc, fragmentAtand, right);
}

static void
monadicBar (antlrcpp::Any &rc, antlrcpp::Any &right)
{
  do_monadic (rc, fragmentBar, right);
}

static void
monadicBang (antlrcpp::Any &rc, antlrcpp::Any &right)
{
  do_monadic (rc, fragmentBang, right);
}

static void
monadicCeil (antlrcpp::Any &rc, antlrcpp::Any &right)
{
  do_monadic (rc, fragmentCeil, right);
}

static void
monadicFloor (antlrcpp::Any &rc, antlrcpp::Any &right)
{
  do_monadic (rc, fragmentFloor, right);
}

static void
monadicRound (antlrcpp::Any &rc, antlrcpp::Any &right)
{
  do_monadic (rc, fragmentRound, right);
}

static void
monadicRand (antlrcpp::Any &rc, antlrcpp::Any &right)
{
  do_monadic (rc, fragmentRand, right);
}

static void
monadicShape (antlrcpp::Any &rc, antlrcpp::Any &right)
{
  if (right.get_typeinfo() == typeid(double)) {
    rc = 0.0;
  }
  else if (right.get_typeinfo()  == typeid(std::vector<double>*)) {
    std::vector<double> *right_vec = right.as<std::vector<double>*>();
    size_t n = right_vec->size ();
    rc = double(n);
  }
  else if (right.get_typeinfo()  == typeid(Matrix*)) {
    Matrix *mtx = right.as<Matrix*>();
    
    rc = mtx->rho ();
  }
  else {
    rc = Error(Error::ERROR_UNKNOWN_DATA_TYPE, ", shape operation");
  }
}

static void
monadicTranspose (antlrcpp::Any &rc, antlrcpp::Any &right)
{
  if (right.get_typeinfo() == typeid(double)) {
    double rv = right.as<double>();
    rc = rv;
  }
  else if (right.get_typeinfo()  == typeid(std::vector<double>*)) {
    std::vector<double> *right_vec = right.as<std::vector<double>*>();
    rc = right_vec;
  }
  else if (right.get_typeinfo() == typeid(Matrix*)) {
    Matrix *rv = right.as<Matrix *>();
      Matrix *xmtx = rv->transpose ();
      if (xmtx) rc = xmtx;
      else {
	rc = Error(Error::ERROR_FAILED_TRANSPOSE);
	std::string em = rv->get_errmsg ();
	std::cout << em << std::endl;
      }
  }
  else {
    rc = Error(Error::ERROR_UNKNOWN_DATA_TYPE, ", transpose operation");
  }
}

static void
monadicRange (antlrcpp::Any &rc, antlrcpp::Any &right)
{
  if (right.get_typeinfo() == typeid(double)) {
    double rv = right.as<double>();
    std::vector<double> *array = new std::vector<double>;
    int n = int (rv);

    /****
	 fixme  do a quad-io index origin
     ****/
    for (int i = 0; i < n; i++)
      array->push_back (double(i));
    rc = array;
  }
  else {
    rc = Error(Error::ERROR_UNKNOWN_DATA_TYPE, ", monadic range operation");
  }
}

typedef struct {
  double  d;
  double  i;
} sort_ety;

typedef bool (*ety_sort_func) (sort_ety i, sort_ety j);

static bool ety_sort_up (sort_ety i, sort_ety j){  return (i.d < j.d); }
static bool ety_sort_dn (sort_ety i, sort_ety j){  return (i.d > j.d); }

static void
do_monadic_grade (antlrcpp::Any &rc, ety_sort_func func, antlrcpp::Any &right)
{
  if (right.get_typeinfo() == typeid(double)) {
    double rv = right.as<double>();
    rc = rv;
  }
  else if (right.get_typeinfo()  == typeid(std::vector<double>*)) {
    std::vector<double> *right_vec = right.as<std::vector<double>*>();
    size_t vlen = right_vec->size ();
    std::vector<sort_ety> *svec = new std::vector<sort_ety>;
    svec->resize (vlen);
    for (size_t p = 0; p < svec->size (); p++) {
      (*svec)[p].d = (*right_vec)[p];
      (*svec)[p].i = double(p);
    }
    std::sort (svec->begin (), svec->end (), func);
    std::vector<double> *res = new std::vector<double>;
    res->resize (vlen);
    for (size_t p = 0; p < svec->size (); p++) {
      sort_ety se =  (*svec)[p];
      (*res)[p] = se.i;
    }
    rc = res;
  }
  else if (right.get_typeinfo() == typeid(Matrix*)) {
    Matrix *rmtx = right.as<Matrix *>();
    std::vector<double> *right_data = rmtx->get_data ();
    size_t vlen = right_data->size ();
    std::vector<sort_ety> *svec = new std::vector<sort_ety>;
    svec->resize (vlen);
    for (size_t p = 0; p < svec->size (); p++) {
      (*svec)[p].d = (*right_data)[p];
      (*svec)[p].i = double(p);
    }
    std::sort (svec->begin (), svec->end (), func);
    std::vector<double> *ndata = new std::vector<double>;
    ndata->resize (vlen);
    for (size_t p = 0; p < svec->size (); p++) {
      sort_ety se =  (*svec)[p];
      (*ndata)[p] = se.i;
    }
    std::vector<size_t> *right_idx = rmtx->get_dims ();
    std::vector<size_t> *ndims = new std::vector<size_t>;
    ndims->resize (right_idx->size ());
    for (size_t p = 0; p < ndims->size (); p++)
      (*ndims)[p] = (*right_idx)[p];
    Matrix *mtx = new Matrix (ndims, ndata);
    rc = mtx;
  }
  else {
    rc = Error(Error::ERROR_UNKNOWN_DATA_TYPE, ", grade operation");
  }
}

static void
monadicGradeUp (antlrcpp::Any &rc, antlrcpp::Any &right)
{
  do_monadic_grade (rc, ety_sort_up, right);
}

static void
monadicGradeDown (antlrcpp::Any &rc, antlrcpp::Any &right)
{
  do_monadic_grade (rc, ety_sort_dn, right);
}

static void
monadicDeterminant (antlrcpp::Any &rc, antlrcpp::Any &right)
{
  Matrix *rv = right.as<Matrix *>();
  double det = rv->determinant ();

  if (std::isnan (det)) {
    rc = Error(Error::ERROR_FAILED_DETERMINANT);
    std::string em = rv->get_errmsg ();
    std::cout << em << std::endl;
  }
  else rc = det;
}

static void
monadicInverse (antlrcpp::Any &rc, antlrcpp::Any &right)
{
  Matrix *rv  = right.as<Matrix *>();
  Matrix *mtx = rv->inverse ();

  if (mtx) rc = mtx;
  else {
    rc = Error(Error::ERROR_FAILED_INVERSE);
    std::string em = rv->get_errmsg ();
    std::cout << em << std::endl;
  }
}

static mfunc mfuncs[] =
{
 nullptr,	// empty	 0
 nullptr,	// DUMMY	 1
 nullptr,	// Number	 2
 nullptr,	// OpStar	 3
 monadicSlash,	// OpSlash	 4
 nullptr,	// OpPlus	 5
 monadicMinus,	// OpMinus	 6
 monadicTranspose,	// OpDollar	 7
 nullptr,	// OpPercent	 8
 nullptr,	// OpHat	 9
 monadicLn,	// OpLn		10
 monadicLog,	// OpLog	11
 monadicExp,	// OpExp	12
 monadicRoot,	// OpRoot	13
 monadicSin,	// OpSin	14
 monadicCos,	// OpCos	15
 monadicTan,	// OpTan	16
 monadicAsin,	// OpAsin	17
 monadicAcos,	// OpAcos	18
 monadicAtan,	// OpAtan	19
 monadicSind,	// OpSind	20
 monadicCosd,	// OpCosd	21
 monadicTand,	// OpTand	22
 monadicAsind,	// OpAsind	23
 monadicAcosd,	// OpAcosd	24
 monadicAtand,	// OpAtand	25
 monadicBar,	// OpBar	26
 monadicBang,	// OpBang	27
 monadicCeil,	// OpCeil	28
 monadicFloor,	// OpFloor	29
 monadicRound,	// OpRound	30
 monadicGradeDown,	// OpLeftAngle	31
 monadicGradeUp,	// OpRightAngle	32
 monadicRand,	// OpQmark	33
 nullptr,	// OpEqual	34
 monadicShape,	// OpPound	35
 nullptr,	// OpComma	36
 nullptr,	// OpColon	37
 nullptr,	// OpBSStar	38
 monadicRange,	// OpColonColon	39
 nullptr,	// OpQEqual	40
 nullptr,	// OpQLeftAngle		41
 nullptr,	// OpQRightAngle	42
 nullptr,	// OpQLeftAngleEqual	43
 nullptr,	// OpQRightAngleEqual	44
 nullptr,	// OpQBangEqual		45
 nullptr,	// OpQBangLeftAngle		46
 nullptr,	// OpQBangRightAngle		47
 nullptr,	// OpQBangLeftAngleEqual	48
 nullptr,	// OpQBangRightAngleEqual	49
 monadicDeterminant,	// OpDet		50
 monadicInverse,	// OpBSSlash		51
};

mfunc
monadic_get_func (int idx)
{
  mfunc res = nullptr;
  if (idx >= 0 && idx < sizeof(mfuncs) / sizeof(mfunc)) {
    res = mfuncs[idx];
  }

  return res;
}
