#ifndef DYADICFCNS_H
#define DYNADICFCNS_H
#include <cmath>

#include <antlr4-runtime.h>


#define RtD(r) ((180.0 / M_PI) * (r))
typedef double (*dyadic_op)(double left, double right);
typedef bool   (*dyadic_test)(double left, double right);
      
typedef void (*dfunc)(antlrcpp::Any &rc,
		      antlrcpp::Any &left, antlrcpp::Any &right);

static void
dyadicStar  (antlrcpp::Any &rc, antlrcpp::Any &left, antlrcpp::Any &right);
static void
dyadicSlash (antlrcpp::Any &rc, antlrcpp::Any &left, antlrcpp::Any &right);
static void
dyadicPlus  (antlrcpp::Any &rc, antlrcpp::Any &left, antlrcpp::Any &right);
static void
dyadicMinus (antlrcpp::Any &rc, antlrcpp::Any &left, antlrcpp::Any &right);
static void
dyadicHat   (antlrcpp::Any &rc, antlrcpp::Any &left, antlrcpp::Any &right);
static void
dyadicLog   (antlrcpp::Any &rc, antlrcpp::Any &left, antlrcpp::Any &right);
static void
dyadicRoot  (antlrcpp::Any &rc, antlrcpp::Any &left, antlrcpp::Any &right);
static void
dyadicAtan  (antlrcpp::Any &rc, antlrcpp::Any &left, antlrcpp::Any &right);
static void
dyadicAtand (antlrcpp::Any &rc, antlrcpp::Any &left, antlrcpp::Any &right);
static void
dyadicEqual (antlrcpp::Any &rc, antlrcpp::Any &left, antlrcpp::Any &right);

void
do_vector_shift (antlrcpp::Any &rc, double shift, antlrcpp::Any &right);

dfunc dyadic_get_func (int idx);

#endif  // DYNADICFCNS_H
