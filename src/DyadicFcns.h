#ifndef DYADICFCNS_H
#define DYNADICFCNS_H
#include <cmath>

#include <antlr4-runtime.h>


#define RtD(r) ((180.0 / M_PI) * (r))
typedef double (*dyadic_op)(double left, double right);
typedef bool   (*dyadic_test)(double left, double right);
      
typedef void (*dfunc)(antlrcpp::Any &rc,
		      antlrcpp::Any &left,
		      antlrcpp::Any &right,
		      antlrcpp::Any &qual);

static void
dyadicStar  (antlrcpp::Any &rc, antlrcpp::Any &left,
	     antlrcpp::Any &right, antlrcpp::Any &qual);
static void
dyadicSlash (antlrcpp::Any &rc, antlrcpp::Any &left,
	     antlrcpp::Any &right, antlrcpp::Any &qual);
static void
dyadicPlus  (antlrcpp::Any &rc, antlrcpp::Any &left,
	     antlrcpp::Any &right, antlrcpp::Any &qual);
static void
dyadicMinus (antlrcpp::Any &rc, antlrcpp::Any &left,
	     antlrcpp::Any &right, antlrcpp::Any &qual);
static void
dyadicHat   (antlrcpp::Any &rc, antlrcpp::Any &left,
	     antlrcpp::Any &right, antlrcpp::Any &qual);
static void
dyadicLog   (antlrcpp::Any &rc, antlrcpp::Any &left,
	     antlrcpp::Any &right, antlrcpp::Any &qual);
static void
dyadicRoot  (antlrcpp::Any &rc, antlrcpp::Any &left,
	     antlrcpp::Any &right, antlrcpp::Any &qual);
static void
dyadicAtan  (antlrcpp::Any &rc, antlrcpp::Any &left,
	     antlrcpp::Any &right, antlrcpp::Any &qual);
static void
dyadicAtand (antlrcpp::Any &rc, antlrcpp::Any &left,
	     antlrcpp::Any &right, antlrcpp::Any &qual);
static void
dyadicEqual (antlrcpp::Any &rc, antlrcpp::Any &left,
	     antlrcpp::Any &right, antlrcpp::Any &qual);

void
do_vector_shift (antlrcpp::Any &rc, double shift,
		 antlrcpp::Any &right, antlrcpp::Any &qual);

dfunc dyadic_get_func (int idx);

#endif  // DYNADICFCNS_H
