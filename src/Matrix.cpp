#include <cmath>
#include <cstring>
#include <antlr4-runtime.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector_double.h>

#include "DyadicFcns.h"
#include "MonadicFcns.h"
#include "Matrix.h"
#include "main.h"

// https://www.gnu.org/software/gsl/doc/html/blas.html

Matrix::Matrix ()
{
  data = nullptr;
  dims = nullptr;
}

// identity
Matrix::Matrix (size_t rank, size_t dim)
{
  size_t cnt =
    static_cast<size_t>(std::pow (static_cast<double>(dim),
				  static_cast<double>(rank)));
  size_t stride = 0;
  size_t scale = 1;
  for (size_t i = 0; i < rank; i++) {
    stride += scale;
    scale *= dim;
  }
  data = new std::vector<double>(cnt, 0.0);
  for (size_t i = 0, off = 0; i < dim; i++, off+=stride) {
    (*data)[off] = 1.0;
  }
  dims = new std::vector<size_t>(rank, dim);
}

Matrix::Matrix (std::vector<size_t> *ndims, std::vector<double> *ndata)
{
  dims = ndims;
  data = ndata;
}

Matrix::Matrix (antlrcpp::Any &dr)
{
  if (dr.get_typeinfo() == typeid(std::vector<double>*)) {
    std::vector<double> *dr_vec = dr.as<std::vector<double>*>();
    size_t rhorho = dr_vec->size ();
    if (rhorho > 0) {
      dims = new std::vector<size_t>(rhorho);
      int count = 1;
      for (int i = 0; i < rhorho; i++) {
	int d = int((*dr_vec)[i]);
	(*dims)[i] = int((*dr_vec)[i]);
	count *= d;
      }
      data = new std::vector<double>(count);
    }
    else {
      data = nullptr;
      dims = nullptr;
    }
  }
  else {
    data = nullptr;
    dims = nullptr;
  }
}

Matrix::Matrix (monadic_op op, Matrix *rv)
{
  data = rv->data;
  dims = rv->dims;
  for (auto rp = (*data).begin (); rp != (*data).end (); rp++)
    (*rp) = (*op)(*rp);
}

Matrix::Matrix (Matrix *rv)
{
  data = new std::vector<double>(rv->data->size ());
  std::memmove (data->data (), rv->data->data (),
		rv->data->size () * sizeof(double));
  dims = new std::vector<size_t>(rv->dims->size ());
  std::memmove (dims->data (), rv->dims->data (),
		rv->dims->size () * sizeof(size_t));
}

Matrix::Matrix (dyadic_op op, double lv, Matrix *rv)
{
  data = rv->data;
  dims = rv->dims;
  for (auto rp = (*data).begin (); rp != (*data).end (); rp++)
    (*rp) = (*op)(lv, (*rp));
}

Matrix::Matrix (dyadic_op op, Matrix *lv, double rv)
{
  data = lv->data;
  dims = lv->dims;
  for (auto rp = (*data).begin (); rp != (*data).end (); rp++)
    (*rp) = (*op)((*rp), rv);
}

Matrix::Matrix (dyadic_op op, Matrix *lv, Matrix *rv)
{
  // fixme assumes isomorphic
  dims = lv->dims;
  data = new  std::vector<double>(lv->data->size ());
  std::vector<double> *ld = lv->data;
  std::vector<double> *rd = rv->data;
  for (int i = 0; i < lv->data->size (); i++)
    (*data)[i] = (*op)((*ld)[i], (*rd)[i]);
}

Matrix::~Matrix ()
{
  //  std::cout << "deleting matrix\n";
  delete dims;
  delete data;
}

static size_t
fix_shift (double dshift, shift_direction_e direction, size_t dim)
{
  size_t rc;
  if (direction == SHIFT_REVERSE) dshift = -dshift;
  if  (dshift < 0.0) dshift += static_cast<double>(dim);
  return static_cast<size_t>(dshift);
}

Matrix *
Matrix::do_shift (Matrix *tp, size_t rhorho, std::vector<size_t>*incrs)
{
  std::vector<size_t> *new_dims = new std::vector<size_t>(rhorho);
  std::vector<double> *new_data = new std::vector<double>(tp->data->size ());
  std::memmove (new_dims->data (), dims->data (),
		rhorho * sizeof (double));
  auto ctrs = std::vector<size_t>(rhorho, 0);
  auto mpys = std::vector<size_t>(rhorho);
  mpys[rhorho - 1] = 1;
  for (int i = rhorho - 2; i >= 0; i--) mpys[i] = mpys[i+1] * (*dims)[i+1];
  bool run = true;
  size_t dest = 0;
  while (run) {
    size_t carry = 1;
    size_t offset = 0;
    for (int i = rhorho - 1; i >= 0; i--) {
      offset += ((ctrs[i] + (*incrs)[i]) % (*dims)[i]) * mpys[i];
      ctrs[i] += carry;
      if (ctrs[i] >= (*dims)[i]) {
	ctrs[i] = 0;
	carry = 1;
      }
      else carry = 0;
    }
    (*new_data)[dest++] = (*data)[offset];
    if (carry) run = false;
  }
  Matrix *mtx = new Matrix (new_dims, new_data);
  return mtx;
}

Matrix *
Matrix::shift (shift_direction_e direction, antlrcpp::Any &axes)
{
  size_t rhorho = dims->size ();
  auto incrs = std::vector<size_t>(rhorho, 0);
  if (axes.get_typeinfo() == typeid(std::vector<double>*)) {
    /***
	<<[a b c]2 3#::6
    ***/
    std::vector<double> *ixs = axes.as<std::vector<double>*>();
    for (size_t i = 0; i < ixs->size (); i++) {
      if ((*ixs)[i] < rhorho) {
	size_t ix = static_cast<size_t>((*ixs)[i]);
	incrs[ix] = fix_shift (1.0, direction, (*dims)[ix]);
      }
      else {
	set_errmsg (std::string ("Invalid axes rank."));
	return nullptr;
      }
    }
  }
  else if (axes.get_typeinfo() == typeid(double)) {
    /***
	<<[ a ]2 3#::6
    ***/
    size_t ix = static_cast<size_t>(axes.as<double>());
    if (ix < rhorho) {
      for (size_t i = 0; i < dims->size (); i++) {
	incrs[ix] = fix_shift (1.0, direction, (*dims)[ix]);
      }
    }
    else {
      set_errmsg (std::string ("Invalid axes rank."));
      return nullptr;
    }
  }
  else if (axes.get_typeinfo() == typeid(nullptr)) {
    /***
	x<<2 3#::6			incrs[rhorho - 1] = shift
    ***/
    incrs[rhorho - 1] = fix_shift (1.0, direction, (*dims)[rhorho - 1]);
  }
  else {
    set_errmsg (std::string ("Unknown axes datatype."));
    return nullptr;
  }

  return do_shift (this, rhorho, &incrs);
}
  
Matrix *
Matrix::shift (shift_direction_e direction, antlrcpp::Any &shift,
	       antlrcpp::Any &axes)
{
  size_t rhorho = dims->size ();
  auto incrs = std::vector<size_t>(rhorho, 0);
  if (axes.get_typeinfo() == typeid(std::vector<double>*)) {
    std::vector<double> *ixs = axes.as<std::vector<double>*>();
    for (size_t i = 0; i < ixs->size (); i++) {
      if ((*ixs)[i] < rhorho) {
	size_t ix = static_cast<size_t>((*ixs)[i]);
	if (shift.get_typeinfo() == typeid(double)) {
	  
	  // x<<[a b c]2 3#::6

	  double dshift = fabs (shift.as<double>());
	  dshift = fmod (dshift, static_cast<double>((*dims)[ix]));
	  if (direction == SHIFT_REVERSE) dshift = -dshift;
	  if  (dshift < 0.0) dshift += static_cast<double>((*dims)[ix]);
	  incrs[ix] = dshift;
	}
#if 0
	else if (shift.get_typeinfo() == typeid(std::vector<double>*)) {
	  
	  // (x y z)<<[a b c]2 3#::6	// unnecessary

	}
#endif
	else {
	  set_errmsg (std::string ("Unknown shift datatype."));
	  return nullptr;
	}
      }
      else {
	set_errmsg (std::string ("Invalid axes rank."));
	return nullptr;
      }
    }
  }
  else if (axes.get_typeinfo() == typeid(double)) {
    size_t ix = static_cast<size_t>(axes.as<double>());
    if (ix < rhorho) {
      if (shift.get_typeinfo() == typeid(double)) {
	  
	// x<<[ a ]2 3#::6

	double dshift = fabs (shift.as<double>());
	for (size_t i = 0; i < dims->size (); i++) {
	  dshift = fmod (dshift, static_cast<double>((*dims)[ix]));
	  if (direction == SHIFT_REVERSE) dshift = -dshift;
	  if (dshift < 0.0) dshift += static_cast<double>((*dims)[ix]);
	  incrs[ix] = dshift;
	}
      }
#if 0
      else if (shift.get_typeinfo() == typeid(std::vector<double>*)) {
	
	// (x y z)<<[ a ]2 3#::6		// makes no sense

      }
#endif
      else {
	set_errmsg (std::string ("Unknown shift datatype."));
	return nullptr;
      }
    }
    else {
      set_errmsg (std::string ("Invalid axes rank."));
      return nullptr;
    }
  }
  else if (axes.get_typeinfo() == typeid(nullptr)) {
    if (shift.get_typeinfo() == typeid(double)) {
      
      // x<<2 3#::6			incrs[rhorho - 1] = shift
      
      double dshift = shift.as<double>();
      dshift = fmod (dshift, static_cast<double>((*dims)[rhorho - 1]));
      if  (dshift < 0.0) dshift += static_cast<double>((*dims)[rhorho - 1]);
      incrs[rhorho - 1] = static_cast<size_t>(dshift);
    }
    else if (shift.get_typeinfo() == typeid(std::vector<double>*)) {

      // (x y z)<<2 3#::6
      
      if (axes.get_typeinfo() == typeid(nullptr)) {
	auto vshift = shift.as<std::vector<double>*>();
	if (rhorho == vshift->size ()) {
	  for (size_t i = 0; i < rhorho; i++) {
	    double dshift = fabs ((*vshift)[i]);
	    dshift = fmod (dshift, static_cast<double>((*dims)[i]));
	    if (direction == SHIFT_REVERSE) dshift = -dshift;
	    if (dshift < 0.0) dshift += static_cast<double>((*dims)[i]);
	    incrs[i] = dshift;
	  }
	}
	else {
	  set_errmsg (std::string ("Invalid shift rank."));
	  return nullptr;
	}
      }
      else {
	set_errmsg (std::string ("Unknown axes datatype."));
	return nullptr;
      }
    }
    else {
      set_errmsg (std::string ("Unknown shift datatype."));
      return nullptr;
    }
  }
  else {
    set_errmsg (std::string ("Unknown axes datatype."));
    return nullptr;
  }

  return do_shift (this, rhorho, &incrs);
}


class DimCounter
{
public:
  DimCounter (std::vector<size_t>* rho,		// 2 3 4
	      std::vector<size_t> *irho,	// 2 4 3
	      std::vector<size_t> *iperm)	// 2 0 1
  {
    rhorho = rho->size ();
    ctrs.resize (rhorho, 0);

    nperm.resize (rhorho, 0);
    std::memmove (nperm.data (), iperm->data (), rhorho * sizeof(double));
    
    incr.resize (rhorho, 0);
    incr[(*iperm)[rhorho - 1]] = 1;
    for (int i =  rhorho - 2; i >= 0; i--) 
      incr[(*iperm)[i]] = incr[(*iperm)[i+1]] * (*rho)[i+1];

    trip.resize (rhorho, 0);
    for (int i = 0, j = rhorho - 1; i < rhorho; i++, j--)
      trip[i] = incr[i] * (*irho)[j];
  }
  
  size_t get_next_index ()
  {
    size_t sum = 0;
    bool carry = 1;
    for (int i = rhorho - 1; i >= 0; i--) {
      sum += ctrs[i];
      if (carry) {
	ctrs[i] += incr[i];
	if (ctrs[i] >= trip[i]) {
	  ctrs[i] = 0;
	  carry = true;
	}
	else carry = false;
      }
    }
    return sum;
  }
  
private:
  size_t rhorho;
  std::vector<size_t> nperm;
  std::vector<size_t> incr;
  std::vector<size_t> trip;
  std::vector<size_t> ctrs;
};
  
//static bool dbl_sort (size_t i, size_t j) {return (i < j);}

static Matrix *
do_transpose (std::vector<size_t>*perm, Matrix *mtx)
{
  if (mtx->get_rhorho () >= 2) {	// no point in "transposing" vec,scalar
    size_t rhorho = mtx->get_rhorho ();
    std::vector<size_t>* iperm;
    if (perm) {
      if (rhorho == perm->size ()) {
	iperm = new std::vector<size_t>(rhorho);
	std::memmove (iperm->data (), perm->data (), rhorho * sizeof(double));
      }
      else {
	delete iperm;
	mtx->set_errmsg (std::string ("Permutation vector \
incompatible with given matrix"));
	return nullptr;
      }
    }
    else {
      iperm = new std::vector<size_t>(rhorho);	
      for (size_t i = 0; i < rhorho; i++)
	(*iperm)[i] = i;
      std::swap ((*iperm)[rhorho - 1], (*iperm)[rhorho - 2]);
    }

    std::vector<size_t> cperm (iperm->size ());
    std::memmove (cperm.data (), iperm->data (),
		  cperm.size () * sizeof(double));
    std::sort (cperm.begin (), cperm.end ());
    for (size_t jj = 0; jj < cperm.size (); jj++) {
      if (cperm[jj] != jj) {
	delete iperm;
  mtx->set_errmsg (std::string ("Permutation contains a non-unique element."));
	return nullptr;
      }
    }
    
    std::vector<size_t>*new_rho =  new std::vector<size_t>(rhorho);
    std::vector<size_t>* rho = mtx->get_rho ();
    for (size_t i = 0; i < rhorho; i++) {
      size_t t = (rhorho - 1) - (*iperm)[i];
      if (t < rhorho) 
	(*new_rho)[t] = (*rho)[i];
      else {
	delete iperm;
	delete new_rho;
	mtx->set_errmsg (std::string ("Permutation index out of range."));
	return nullptr;
      }
    }

    std::vector<double>*new_data =
      new std::vector<double>(mtx->size ());

    DimCounter tc = DimCounter (rho, new_rho, iperm);

    std::vector<double> *data = mtx->get_data ();
    for (int i = 0; i < new_data->size (); i++) {
      int to = tc.get_next_index ();
      (*new_data)[i] = (*data)[to];
    }

    std::vector<size_t>*c_rho =  new std::vector<size_t>(rhorho);
    for (size_t i = 0; i < rhorho; i++)
      (*c_rho)[i] = (*new_rho)[(rhorho - 1) - i];
    Matrix *nm = new Matrix (c_rho, new_data);

    delete iperm;
    delete new_rho; 
    return nm; 
  }
  return nullptr; 
}

std::vector<double> *
Matrix::solve (int direction, std::vector<double> *vec)
{
  Matrix *inv = this->inverse ();
  std::vector<double> *res = nullptr;
  if (inv) {
    res = this->multiply (direction, vec);
    delete inv;
  }
  return res;
}

Matrix *
Matrix::solve (Matrix *mtx)
{
#if 1
  //  ./mpl -e '(2 2#4::7) \/ 2 2#1::4'
  //  sb -0.5  1.5
  //     -1.5 2.5
  Matrix *inv = mtx->inverse ();
  Matrix *res = nullptr;
  if (inv) {
    res = this->multiply (inv);
    delete inv;
  }
  return res;

#else

  // needs work

  
  Matrix *nmtx = nullptr;
  gsl_matrix *L = nullptr;
  gsl_matrix *R = nullptr;
  
  size_t l_rows = (*dims)[0];
  size_t l_cols = (*dims)[1];
  std::vector<size_t> *r_dims = mtx->get_dims ();
  size_t r_rows = (*r_dims)[0];
  size_t r_cols = (*r_dims)[1];
  if (dims->size () == 2 &&
      r_dims->size () == 2 &&
      l_rows == l_cols &&
      r_rows == r_cols &&
      l_rows == r_rows) {

#if 0
    std::vector<double> *l_data = data;
    std::vector<double> *r_data = mtx->get_data ();
#else
    std::vector<double> *r_data = data;
    std::vector<double> *l_data = mtx->get_data ();
#endif
    
    L = gsl_matrix_alloc(l_rows, l_cols);
    for (size_t i = 0; i < data->size (); i++)
      L->data[i] = (*l_data)[i];
    
    R = gsl_matrix_alloc(r_rows, r_cols);
    for (size_t i = 0; i < r_data->size (); i++)
      R->data[i] = (*r_data)[i];

    int rc;
    int signum;
    gsl_permutation *P = gsl_permutation_calloc (l_rows);

    rc = gsl_linalg_LU_decomp(L, P, &signum);

    //  ./mpl -e '(2 2#4::7) \/ 2 2#1::4'
    //  sb -0.5  1.5
    //     -1.5 2.5
    //
    // CblasRight, CblasUpper,  CblasNoTrans, CblasNonUnit 
    // is 0.166667 2.5
    //    0.5      1.5       

    if (rc == 0) {
      // CBLAS_SIDE_t:       CblasLeft,              CblasRight￼
      // ￼                   B = \alpha op(inv(A))B  B = \alpha B op(inv(A))
      // CBLAS_UPLO_t:       CblasUpper, CblasLower
      // CBLAS_TRANSPOSE_t:  CblasNoTrans, CblasTrans, CblasConjTrans
      // CBLAS_DIAG_t:       CblasNonUnit, CblasUnit

#if 1
      gsl_vector *xv = nullptr;
      for (size_t j = 0; j < L->size2; j++) {
	gsl_vector_view x = gsl_matrix_row (L, j);
	xv = &(x.vector);
	gsl_linalg_LU_svx (R, P, xv);
	rc = gsl_linalg_LU_svx (R, P, xv);
	for (size_t k = 0; k < xv->size; k++)
	  std::cout << xv->data[k] << " ";
	std::cout << std::endl;
      }
#else
      rc = gsl_blas_dtrsm(CblasLeft, CblasUpper,  CblasNoTrans,
			  CblasNonUnit, 1.0, L, R);
#endif
      if (rc == 0) {	
	std::vector<size_t> *ndims = new std::vector<size_t>(2);
	ndims->resize (2);
	(*ndims)[0] = l_rows;
	(*ndims)[1] = l_cols;
	std::vector<double> *ndata = new std::vector<double>(l_rows * l_cols);
	ndata->resize (l_rows * l_cols);
	for (size_t o = 0; o < (l_rows * l_rows); o++)
	  (*ndata)[o] = R->data[o];
    
	nmtx = new Matrix (ndims, ndata);
      }
      else {
      }
    }
    else {
    }
    
  }
  else {
    set_errmsg (std::string ("Only two-dimensional matrices supported for matrix-matrix solve"));
  }
  if (L) gsl_matrix_free (L);
  if (R) gsl_matrix_free (R);

  return nmtx;
#endif
}

std::vector<double> *
Matrix::multiply (int direction, std::vector<double> *vec)
{
  std::vector<double> *res =  nullptr;

  gsl_matrix *M = nullptr;
  gsl_vector *V = nullptr;
  gsl_vector *Y = nullptr;
  
  if (dims->size () == 2) {
    size_t rows = (*dims)[0];
    size_t cols = (*dims)[1];
    if (direction == VM_VEC_RIGHT) {
      if (cols == vec->size ()) {
	M = gsl_matrix_alloc(rows, cols);
	for (size_t i = 0; i < data->size (); i++)
	  M->data[i] = (*data)[i];
	V = gsl_vector_alloc (vec->size ());
	std::memmove (V->data, vec->data (), vec->size () * sizeof(double));
	Y = gsl_vector_calloc (rows);
	int rc = gsl_blas_dgemv(CblasNoTrans, 1.0, M, V, 0.0, Y);
	if (0 == rc) {
	  res = new std::vector<double>(rows);
	  res->resize (rows);
	  for (size_t i = 0; i < rows; i++)
	    (*res)[i] = Y->data[i];
	}
	else {
	  set_errmsg (std::string ("Failure in matrix-vector multiplication."));
	}
      }
      else {
	set_errmsg (std::string ("Dimension mismatch in matrix-vector multiplication."));
      }
    }
    else {
      if (rows == vec->size ()) {
	gsl_matrix *M = gsl_matrix_alloc(rows, cols);
	std::memmove (M->data, data->data (), data->size () * sizeof(double));
	gsl_vector *V = gsl_vector_alloc (vec->size ());
	std::memmove (V->data, vec->data (), vec->size () * sizeof(double));
	gsl_vector *Y = gsl_vector_calloc (cols);
	int rc = gsl_blas_dgemv(CblasTrans, 1.0, M, V, 0.0, Y);
	if (0 == rc) {
	  res = new std::vector<double>(cols);
	  res->resize (cols);
	  for (size_t i = 0; i < cols; i++)
	    (*res)[i] = Y->data[i];
	}
	else {
	  set_errmsg (std::string ("Failure in matrix-vector multiplication."));
	}
      }
      else {
	set_errmsg (std::string ("Dimension mismatch in vector-matrix multiplication."));
      }
    }
  }
  else {
    set_errmsg (std::string ("Only two-dimensional matrices supported for vector-matrix multiplication"));
  }

  if (M) gsl_matrix_free (M);
  if (V) gsl_vector_free (V);
  if (Y) gsl_vector_free (Y);
  return res;
}

Matrix *
Matrix::multiply (Matrix *rv)
{
  Matrix *mtx = nullptr;
  std::vector<size_t> *l_rho = dims;
  std::vector<size_t> *r_rho = rv->get_dims ();
  size_t l_rows = (*l_rho)[0];
  size_t l_cols = (*l_rho)[1];
  size_t r_rows = (*r_rho)[0];
  size_t r_cols = (*r_rho)[1];

  if (dims->size () == 2 &&
      rv->get_rhorho () == 2) {
    if (l_cols == r_rows) {
      size_t c_rows = l_rows;
      size_t c_cols = r_cols;

      std::vector<double>* ld = data;
      std::vector<double>* rd = rv->get_data ();

      double *l_data = new double[l_rows * l_cols];
      double *r_data = new double[r_rows * r_cols];
      double *c_data = new double[c_rows * c_cols]();

      std::memmove (l_data, ld->data (), l_rows * l_cols * sizeof(double));
      std::memmove (r_data, rd->data (), r_rows * r_cols * sizeof(double));

      gsl_matrix_view L = gsl_matrix_view_array (l_data, l_rows, l_cols);
      gsl_matrix_view R = gsl_matrix_view_array (r_data, r_rows, r_cols);
      gsl_matrix_view C = gsl_matrix_view_array (c_data, c_rows, c_cols);
    
      gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
		      1.0, &L.matrix, &R.matrix,
		      0.0, &C.matrix);

      std::vector<size_t> *ndims = new std::vector<size_t>(2);
      ndims->resize (2);
      (*ndims)[0] = c_rows;
      (*ndims)[1] = c_cols;
      std::vector<double> *ndata = new std::vector<double>(c_rows * c_cols);
      ndata->resize (c_rows * c_cols);
      std::memmove (ndata->data (), c_data, c_rows * c_cols * sizeof(double));
    
      mtx = new Matrix (ndims, ndata);
    
      delete l_data;
      delete r_data;
      delete c_data;
    }
    else {
      set_errmsg (std::string ("Dimension mismatch in matrix multiplication."));
    }
  }
  else {
    set_errmsg (std::string ("Only two-dimensional matrices supported for matrix multiplication"));
  }
  return mtx;
}

Matrix *
Matrix::transpose (antlrcpp::Any &permutation)
{
  Matrix *rc = nullptr;
  if (permutation.get_typeinfo() != typeid(nullptr)) {
    std::vector<double> *perm = permutation.as<std::vector<double>*>();
    std::vector<size_t> *iperm = new std::vector<size_t>(perm->size ());
    for (size_t i = 0; i < perm->size (); i++)
      (*iperm)[i] = int((*perm)[i]);
    rc = do_transpose (iperm, this);
    delete iperm;
  }
  else rc = do_transpose (nullptr, this);
  return rc;
}

double
Matrix::sum ()
{
  double rc = NAN;
  if (data->size () > 0) {
    rc = 0.0;
    for (auto rp = (*data).begin (); rp != (*data).end (); rp++)
      rc += (*rp);
  }
  else {
    set_errmsg (std::string ("Empty matrix"));
  }
  return rc;
}

double
Matrix::product ()
{
  double rc = NAN;
  if (data->size () > 0) {
    rc = 1.0;
    for (auto rp = (*data).begin (); rp != (*data).end (); rp++)
      rc *= (*rp);
  }
  else {
    set_errmsg (std::string ("Empty matrix"));
  }
  return rc;
}

double
Matrix::determinant ()
{
  double det = nan ("");

  size_t rows = (*dims)[0];
  size_t cols = (*dims)[1];

  if (rows == cols) {
    gsl_matrix *C = gsl_matrix_alloc ((*dims)[1], (*dims)[0]);
    std::memmove (C->data, data->data (), data->size () * sizeof(double));
    gsl_permutation *P = gsl_permutation_calloc (rows);

    int signum;
    int rc =  gsl_linalg_LU_decomp(C, P, &signum);
    if (rc == 0) 
      det = gsl_linalg_LU_det(C, signum);
    else {
      set_errmsg (std::string ("Determinant failed in decompoosition."));
    }
    if (C) gsl_matrix_free (C);
    if (P)  gsl_permutation_free (P);
  }
  else {
    set_errmsg (std::string ("Determinant requires a square matrix."));
  }
  return det;
}

Matrix *
Matrix::inverse ()
{
  Matrix *mtx = nullptr;
  
  size_t rows = (*dims)[0];
  size_t cols = (*dims)[1];

  if (rows == cols) {
    gsl_matrix *C = gsl_matrix_alloc ((*dims)[1], (*dims)[0]);
    std::memmove (C->data, data->data (), data->size () * sizeof(double));
    gsl_permutation *P = gsl_permutation_calloc (rows);
    gsl_matrix *I = gsl_matrix_alloc (rows, cols);

    int signum;
    int rc =  gsl_linalg_LU_decomp(C, P, &signum);
    if (rc == 0) {
      rc = gsl_linalg_LU_invert(C, P, I);
      if (0 == rc) {
	std::vector<size_t> *ndims = new std::vector<size_t>(2);
	ndims->resize (2);
	(*ndims)[0] = rows;
	(*ndims)[1] = cols;
	std::vector<double> *ndata =
	  new std::vector<double>(data->size ());
	ndata->resize (data->size ());
	for (size_t i = 0; i < data->size (); i++)
	  (*ndata)[i] = I->data[i];
	mtx = new Matrix (ndims, ndata);
      }
    }
    else {
      set_errmsg (std::string ("Determinant failed in decompoosition."));
    }
    if (C) gsl_matrix_free (C);
    if (P) gsl_permutation_free (P);
    if (I) gsl_matrix_free (I);
  }
  else {
    set_errmsg (std::string ("Determinant requires a square matrix."));
  }	

  return mtx;
}

 bool
 Matrix::isomorphic (Matrix *rv)
 {
   bool rc = false;
   if (dims->size () == rv->dims->size () &&
       data->size () == rv->data->size ()) {
     size_t i;
     for (i = 0; i < dims->size (); i++) {
       if (dims[i] != rv->dims[i]) break;
     }
     if (i == (dims->size () - 1))
       rc = true;
   }
   return rc;
 }

void
Matrix::fill (double v)
{
  for (auto rp = (*data).begin (); rp != (*data).end (); rp++) (*rp) = v;
}

void
Matrix::copy (antlrcpp::Any &vec)
{
  std::vector<double> *vals = vec.as<std::vector<double>*>();
  size_t nr_vals = vals->size ();
  for (int i = 0; i < (*data).size (); i++)
    (*data)[i] = (*vals)[i%nr_vals];
}

void
Matrix::reshape (antlrcpp::Any &vec)
{
  Matrix *omtx = vec.as<Matrix*>();
  size_t nr_vals = data->size ();
  std::memmove (data->data (), omtx->data->data (), nr_vals * sizeof(double));
}

double &
Matrix::operator[](size_t i)
{
  return (*data)[i];
}
  

antlrcpp::Any
Matrix::rho ()
{
  antlrcpp::Any rc;
  std::vector<double> *rdims =  new std::vector<double>(dims->size ());
  for (int i = 0; i < dims->size (); i++)
    (*rdims)[i] = double((*dims)[i]);
  rc = rdims;
  return rc;
}
  
antlrcpp::Any
Matrix::values ()
{
  antlrcpp::Any rc = data;
  return rc;
}


// diagnostics
  
void
Matrix::print_rhorho ()
{
  std::cout << "rhorho " << (*dims).size () << std::endl;
}

void
Matrix::print_rho ()
{
  for (auto rp = (*dims).begin (); rp != (*dims).end (); rp++)
    std::cout << *rp << " ";
  std::cout << std::endl;
}

void
Matrix::print ()
{
  size_t rhorho = dims->size ();
  switch(rhorho) {
  case 0:			// empty
    std::cout << "''";
    break;
  case 1:			// vector
    {
      if (!isTestMode ()) std::cout << ">> ";
      for (auto rp = (*data).begin (); rp != (*data).end (); rp++)
	std::cout << *rp << " ";
      std::cout << std::endl;
    }
    break;
  case 2:			// 2d matrix
    {
      for (int i = 0; i < (*data).size (); i++) {
	if (0 == i % (*dims)[1]) {
	  std::cout << std::endl;
	  if (!isTestMode ()) std::cout << ">> ";
	}
	std::cout << (*data)[i] << " ";
      }
      std::cout << std::endl;
    }
    break;
  default:			// higher order
    {
      int bp1 = (*dims)[rhorho - 1];
      int bp2 = bp1 * (*dims)[rhorho - 2];
      std::vector<int> bp;

      int acc = bp2;
      for (int kk = 0; kk < rhorho - 2; kk++) {
	bp.push_back (acc);
	acc *= (*dims)[(rhorho - 3) -kk];
      }

      for (int i = 0; i < (*data).size (); i++) {
	if (i > 0 && 0 == i % bp1) std::cout <<  "\n  ";
	if (0 == i % bp2) {
	  if (!isTestMode ()) {
	    if (i > 0) std::cout << "\n\n";
	    std::cout << ">> [";
	    int j, k;
	    for (k = 0, j = rhorho - 3; j >= 0; j--, k++) 
	      std::cout << (i / (bp[j])) % (*dims)[k]    << " ";
	    std::cout << " * *]:\n  ";
	  }
	}
	std::cout << (*data)[i] << " ";
      }
      std::cout << "\n\n";
    }
    break;
  }
}


  // internals
  
double
Matrix::get_rhorho ()
{
  return double((*dims).size ());
}

std::vector<size_t>* 
Matrix::get_rho ()
{
  return dims;
}

std::vector<double>*
Matrix::get_data ()
{
  return data;
}

std::vector<size_t>*
Matrix::get_dims ()
{
  return dims;
}

size_t
Matrix::size ()
{
  return (*data).size ();
}

void
Matrix::clear ()
{
  delete dims;
  delete data;
}

void
Matrix::clear (bool clear_dims, bool clear_data)
{
  if (clear_dims) delete dims;
  if (clear_data) delete data;
}

void
Matrix::set_data (std::vector<double> *ndata)
{
  data = ndata;
}
void
Matrix::set_dims (std::vector<size_t> *ndims)
{
  dims = ndims;
}

void
Matrix::set_errmsg (std::string str)
{
  errmsg = str;
  //  errmsg = str.assign (str);
}

std::string
Matrix::get_errmsg ()
{
  return errmsg;
}

double
Matrix::get_value (std::vector<double>* idx)
{
  double rc = nan ("");
  if (dims->size () == idx->size ()) {
    double mpy = 1;
    double val = 0;
    for (int i = idx->size () - 1; i >= 0; i--) {
      val += mpy * (*idx)[i];
      mpy *= (*dims)[i];
    }

    size_t ix = size_t (val);
    if (ix < data->size ()) rc = (*data)[ix];
    else {
      set_errmsg (std::string ("Index vector out of range of given matrix"));
    }
  }
  else {
    set_errmsg (std::string ("Index vector incompatible with given matrix"));
  }
  return rc;
}
