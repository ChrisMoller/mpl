//#include <vector>
#include <cmath>
#include <antlr4-runtime.h>
#include "DyadicFcns.h"
#include "MonadicFcns.h"
#include "Matrix.h"

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
      // fixme error now do i report that?
    }
  }
  else {
    // error
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

static bool dbl_sort (size_t i, size_t j) {return (i < j);}

static bool
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
	return false;
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
	return false;
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
	return false;
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

    mtx->clear (true, true);
    mtx->set_data (new_data);
    mtx->set_dims (new_rho);
    delete iperm;
  }
  return true;
}

bool
Matrix::transpose ()
{
  return do_transpose (nullptr, this);
}

bool
Matrix::transpose (antlrcpp::Any &permutation)
{
  std::vector<double> *perm = permutation.as<std::vector<double>*>();
  std::vector<size_t> *iperm = new std::vector<size_t>(perm->size ());
  for (size_t i = 0; i < perm->size (); i++)
    (*iperm)[i] = int((*perm)[i]);
  bool rc = do_transpose (iperm, this);
  delete iperm;
  return rc;
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
      std::cout << ">> ";
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
	  std::cout << ">> ";
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
	  if (i > 0) std::cout << "\n\n";
	  std::cout << ">> [";
	  int j, k;
	  for (k = 0, j = rhorho - 3; j >= 0; j--, k++) 
	    std::cout << (i / (bp[j])) % (*dims)[k]    << " ";
	  std::cout << " * *]:\n  ";
	}
	else if (i > 0 && 0 == i % bp1) std::cout <<  "\n  ";
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
      set_errmsg (std::string ("index vector out of range of given matrix"));
    }
  }
  else {
    set_errmsg (std::string ("index vector incompatible with given matrix"));
  }
  return rc;
}
