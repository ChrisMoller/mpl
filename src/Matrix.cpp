//#include <vector>
#include <cmath>
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
  data = rv->data;
  dims = rv->dims;
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
  //  std::cout << "deleting\n";
  delete dims;
  delete data;
}

class DimCounter
{
public:
  DimCounter (std::vector<size_t> *irho, std::vector<size_t> *iperm)
  {
    rho  = irho;
    perm = iperm;
    size_t rhorho = rho->size ();
    ctr.resize (rhorho);
    stride.resize (rhorho);
    size_t mpy = 1;
    size_t p = (*perm)[rhorho - 1];
    stride[p] = 1;
    for (int ix = rhorho - 2; ix >=0; ix--) {
      p = (*perm)[ix];
      mpy = stride[p] = mpy * (*rho)[ix + 1];
    }
  }

  size_t get_next_index ()
  {
    int rhorho = ctr.size ();
    size_t idx = 0;

    for (int i = rhorho -1; i >= 0; i--) {
      idx += ctr[i] * stride[i];
    }

    int carry = 1;
    for (int i = rhorho -1; i >= 0; i--) {
      ctr[i] += carry;
      size_t p = (*perm)[i];
      if (ctr[i] >= (*rho)[p]) {
	ctr[i] = 0;
	carry = 1;
      }
      else carry = 0;
    }
    
    return idx;
  }

private:
  std::vector<size_t> *rho;
  std::vector<size_t> *perm;
  std::vector<size_t> ctr;
  std::vector<size_t> stride;
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
	for (size_t i = 0; i < rhorho; i++)
	  (*iperm)[i] = (*perm)[i];
      }
      else {
	delete iperm;
	mtx->set_errmsg (std::string ("Permutation vector incompatible with given matrix"));
	return nullptr;
      }
    }
    else {
      iperm = new std::vector<size_t>(rhorho);	
      for (size_t i = 0; i < rhorho; i++)
	(*iperm)[i] = i;
      std::swap ((*iperm)[rhorho - 2], (*iperm)[rhorho - 1]);
    }

    std::vector<size_t> cperm (iperm->size ());
    for (size_t jj = 0; jj < cperm.size (); jj++)
      cperm[jj] = (*iperm)[jj];
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
      size_t t = (*iperm)[i];
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

    DimCounter tc = DimCounter (new_rho, iperm);

    std::vector<double> *data = mtx->get_data ();
    for (int i = 0; i < new_data->size (); i++) {
      int to = tc.get_next_index ();
      (*new_data)[to] = (*data)[i];
    }

    Matrix *nm = new Matrix (new_rho, new_data);

    delete iperm;
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
	for (size_t i = 0; i < vec->size (); i++)
	  V->data[i] = (*vec)[i];
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
	for (size_t i = 0; i < data->size (); i++)
	  M->data[i] = (*data)[i];
	gsl_vector *V = gsl_vector_alloc (vec->size ());
	for (size_t i = 0; i < vec->size (); i++)
	  V->data[i] = (*vec)[i];
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
      double *c_data = new double[c_rows * c_cols];
	
      for (size_t o = 0; o < (l_rows * l_cols); o++)
	l_data[o] = (*ld)[o];
      for (size_t o = 0; o < (r_rows * r_cols); o++)
	r_data[o] = (*rd)[o];
      for (size_t o = 0; o < (c_rows * c_cols); o++)
	c_data[o] = 0.0;

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
      for (size_t o = 0; o < (c_rows * c_rows); o++)
	(*ndata)[o] = c_data[o];
    
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
  std::vector<double> *perm = permutation.as<std::vector<double>*>();
  std::vector<size_t> *iperm = new std::vector<size_t>(perm->size ());
  for (size_t i = 0; i < perm->size (); i++)
    (*iperm)[i] = int((*perm)[i]);
  Matrix * rc = do_transpose (iperm, this);
  delete iperm;
  return rc;
}

Matrix *
Matrix::transpose ()
{
  return do_transpose (nullptr, this);
}

double
Matrix::determinant ()
{
  double det = nan ("");

  size_t rows = (*dims)[0];
  size_t cols = (*dims)[1];

  if (rows == cols) {
    gsl_matrix *C = gsl_matrix_alloc ((*dims)[1], (*dims)[0]);
    for (size_t i = 0; i < data->size (); i++)
      C->data[i] = (*data)[i];
    gsl_permutation *P = gsl_permutation_calloc (rows);

    int signum;
    int rc =  gsl_linalg_LU_decomp(C, P, &signum);
    if (rc == 0) {
      det = gsl_linalg_LU_det(C, signum);
    }
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
    for (size_t i = 0; i < data->size (); i++)
      C->data[i] = (*data)[i];
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
   return 
     dims->size () == rv->dims->size () &&
     data->size () == rv->data->size () &&
     dims == rv->dims;
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
	if (0 == i % bp2) {
	  if (!isTestMode ()) {
	    if (i > 0) std::cout << "\n\n";
	    std::cout << ">> [";
	    int j, k;
	    for (k = 0, j = rhorho - 3; j >= 0; j--, k++) 
	      std::cout << (i / (bp[j])) % (*dims)[k]    << " ";
	    std::cout << " * *]:\n  ";
	  }
	  else if (i > 0 && 0 == i % bp1) std::cout <<  "\n  ";
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
