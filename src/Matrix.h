//#include <vector>
#include <antlr4-runtime.h>
#include "DyadicFcns.h"
#include "MonadicFcns.h"



class  Matrix {
public:
  Matrix ();
  Matrix (std::vector<size_t> *ndims, std::vector<double> *ndata);
  Matrix (monadic_op op, Matrix *rv);
  Matrix (antlrcpp::Any &dr);
  Matrix (dyadic_op op, double lv, Matrix *rv);
  Matrix (dyadic_op op, Matrix *lv, double rv);
  Matrix (dyadic_op op, Matrix *lv, Matrix *rv);
  Matrix (Matrix *rv);
  ~Matrix ();

  bool isomorphic (Matrix *rv);
  bool transpose ();
  bool transpose (antlrcpp::Any &permutation);
  void fill (double v);
  void copy (antlrcpp::Any &v);
  double &operator[](size_t i);
  antlrcpp::Any rho ();
  antlrcpp::Any values ();
  void print ();

  // diagnostics
  
  void print_rhorho ();

  void print_rho ();

  // internals
  
  double get_rhorho ();
  std::vector<size_t>* get_rho ();
  size_t size ();
  void clear ();
  void clear (bool clear_dims, bool clear_data);
  void set_data (std::vector<double> *ndata);
  void set_dims (std::vector<size_t> *ndims);
  std::vector<double>* get_data ();
  std::vector<size_t>* get_dims ();
  void set_errmsg (std::string str);
  std::string get_errmsg ();
  

private:
  std::vector<double> *data;
  std::vector<size_t> *dims;
  std::string errmsg;
};
