#ifndef MONADICFCNS_H
#define MONADICFCNS_H
#include <cmath>

#include <antlr4-runtime.h>

typedef double (*monadic_op)(double right);
typedef void (*mfunc)(antlrcpp::Any &rc,
		      antlrcpp::Any &right,
		      antlrcpp::Any &qual);

static void
monadicMinus (antlrcpp::Any &rc, antlrcpp::Any &right);
static void
monadicLn    (antlrcpp::Any &rc, antlrcpp::Any &right);
static void
monadicLog   (antlrcpp::Any &rc, antlrcpp::Any &right);
static void
monadicExp   (antlrcpp::Any &rc, antlrcpp::Any &right);
static void
monadicRoot  (antlrcpp::Any &rc, antlrcpp::Any &right);
static void
monadicSin   (antlrcpp::Any &rc, antlrcpp::Any &right);
static void
monadicCos   (antlrcpp::Any &rc, antlrcpp::Any &right);
static void
monadicTan   (antlrcpp::Any &rc, antlrcpp::Any &right);
static void
monadicAsin  (antlrcpp::Any &rc, antlrcpp::Any &right);
static void
monadicAcos  (antlrcpp::Any &rc, antlrcpp::Any &right);
static void
monadicAtan  (antlrcpp::Any &rc, antlrcpp::Any &right);
static void
monadicSind  (antlrcpp::Any &rc, antlrcpp::Any &right);
static void
monadicCosd  (antlrcpp::Any &rc, antlrcpp::Any &right);
static void
monadicTand  (antlrcpp::Any &rc, antlrcpp::Any &right);
static void
monadicAsind (antlrcpp::Any &rc, antlrcpp::Any &right);
static void
monadicAcosd (antlrcpp::Any &rc, antlrcpp::Any &right);
static void
monadicAtand (antlrcpp::Any &rc, antlrcpp::Any &right);
static void
monadicBar   (antlrcpp::Any &rc, antlrcpp::Any &right);
static void
monadicBang  (antlrcpp::Any &rc, antlrcpp::Any &right);
static void
monadicCeil  (antlrcpp::Any &rc, antlrcpp::Any &right);
static void
monadicFloor (antlrcpp::Any &rc, antlrcpp::Any &right);
static void
monadicRound (antlrcpp::Any &rc, antlrcpp::Any &right);

mfunc monadic_get_func (int idx);

#endif // MONADICFCNS_H
