#ifndef EUCLIDEAN_MATRIX_H
#define EUCLIDEAN_MATRIX_H

#include <iostream>
#include <fstream>
#include <assert.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>

#include "euclidean_vector.h"

template <int D, typename T=double>
class euclidean_matrix {
  
  public:
  T x[D][D];
  
  public:
  
  euclidean_matrix();
  euclidean_matrix(const T[D][D]);           
  euclidean_matrix(const euclidean_matrix&);
  euclidean_matrix(const gsl_matrix*);
  ~euclidean_matrix();

  euclidean_matrix<D, T>& operator= (T val);
  euclidean_matrix<D, T>& operator= (euclidean_matrix<D,T>& m);
  euclidean_matrix<D, T>& operator= (gsl_matrix* m);
  euclidean_matrix<D, T>& operator+=(const euclidean_matrix<D, T>&);
  euclidean_matrix<D, T>& operator-=(const euclidean_matrix<D, T>&);
  euclidean_matrix<D, T>& operator*=(const T);
  euclidean_matrix<D, T>& operator/=(const T);

  euclidean_matrix<D, T>& invert();
  euclidean_matrix<D, T>& transpose();           // Transpose the matrix (not "get THE transpose")
  euclidean_matrix<D, T>  get_transpose() const; // Get the transpose (don't alter *this)

  euclidean_matrix<D, T> operator+(const euclidean_matrix<D, T>&) const;
  euclidean_matrix<D, T> operator-(const euclidean_matrix<D, T>&) const;
  euclidean_matrix<D, T> operator*(const euclidean_matrix<D, T>&) const; //A*B (w/ matrix)
  euclidean_vector<D, T> operator*(const euclidean_vector<D, T>&) const; //A*x (w/ vector)
  euclidean_matrix<D, T> operator/(const T) const;
  euclidean_matrix<D, T> operator%(const T) const;
  bool operator==(const euclidean_matrix<D, T> &a) const;


  euclidean_matrix       <D, int>    Integer() const;
  euclidean_matrix       <D, double> Double () const;

  static euclidean_matrix<D, int>    Integer(const euclidean_matrix<D, T>&); 
  static euclidean_matrix<D, double> Double (const euclidean_matrix<D, T>&); 
  
  T* operator[](const unsigned int);
 
  double        dot(const euclidean_matrix<D, T>&) const;
  static double dot(const euclidean_matrix<D, T>&, const euclidean_matrix<D, T>&);
  
  double        norm_sq() const;
  static double norm_sq(const euclidean_matrix<D, T>&);

  T        determinant() const;
  static T determinant(const euclidean_matrix<D, T>&);

  void read(std::ifstream&);
  void write(std::ofstream&) const;
};


template <int D, typename T>
std::ostream& operator<<(std::ostream&, const euclidean_matrix<D, T>&);




//==================================================================
// Implementation



// constructor
// ~~~~~~~~~~~
template <int D, typename T>
euclidean_matrix<D, T>::euclidean_matrix()
{
  for   (int i=0; i<D; i++)
    for (int j=0; j<D; j++)
      x[i][j] = 0;
}



template <int D, typename T>
euclidean_matrix<D, T>::euclidean_matrix(const T x_i[D][D])
{
  for   (int i=0; i<D; i++)
    for (int j=0; j<D; j++)
      x[i][j] = x_i[i][j];
}


template <int D, typename T>
euclidean_matrix<D, T>::euclidean_matrix(const euclidean_matrix<D, T> &m)
{
  for   (int i=0; i<D; i++)
    for (int j=0; j<D; j++)
      x[i][j] = m.x[i][j];
}

template <int D, typename T>
euclidean_matrix<D, T>::euclidean_matrix(const gsl_matrix* m)
{
  assert(m->size1 == D && m->size2==D);

  for   (int i=0; i<D; i++)
    for (int j=0; j<D; j++)
      x[i][j] = gsl_matrix_get(m,i,j);
}


// destructor
// ~~~~~~~~~~
template <int D, typename T>
euclidean_matrix<D, T>::~euclidean_matrix()
{
}


// =
// ~~
template <int D, typename T>
inline euclidean_matrix<D, T>& euclidean_matrix<D, T>::operator=(T val)
{
  for   (int i=0; i<D; i++)
    for (int j=0; j<D; j++)
      x[i][j] = val;

  return *this;
}

template <int D, typename T>
inline euclidean_matrix<D, T>& euclidean_matrix<D, T>::operator=(euclidean_matrix<D, T>& m)
{
  for   (int i=0; i<D; i++)
    for (int j=0; j<D; j++)
      x[i][j] = m[i][j];

  return *this;
}

template <int D, typename T>
inline euclidean_matrix<D, T>& euclidean_matrix<D, T>::operator=(gsl_matrix* m)
{
  assert(m->size1==DIM && m->size2==DIM);
  for   (int i=0; i<D; i++)
    for (int j=0; j<D; j++)
      x[i][j] = gsl_matrix_get(m,i,j);

  return *this;
}



// +=
// ~~
template <int D, typename T>
inline euclidean_matrix<D, T>& euclidean_matrix<D, T>::operator+=(const euclidean_matrix<D, T> &m)
{
  for   (int i=0; i<D; i++)
    for (int j=0; j<D; j++)
      x[i][j] += m.x[i][j];

  return *this;
}



// -=
// ~~
template <int D, typename T>
inline euclidean_matrix<D, T>& euclidean_matrix<D, T>::operator-=(const euclidean_matrix<D, T> &m)
{
  for   (int i=0; i<D; i++)
    for (int j=0; j<D; j++)
      x[i][j] -= m.x[i][j];

  return *this;
}



// *=
// ~~
template <int D, typename T>
inline euclidean_matrix<D, T>& euclidean_matrix<D, T>::operator*=(const T s)
{
  for   (int i=0; i<D; i++)
    for (int j=0; j<D; j++)
      x[i][j] *= s;

  return *this;
}



// /=
// ~~
template <int D, typename T>
inline euclidean_matrix<D, T>& euclidean_matrix<D, T>::operator/=(const T s)
{
  for   (int i=0; i<D; i++)
    for (int j=0; j<D; j++)
      x[i][j] /= s;

  return *this;
}



// Invert the matrix
// ~~~~~~~~~~~~~~~~~~~~~~
// (Use GSL libraries)
template <int D,typename T>
inline euclidean_matrix<D, T>& euclidean_matrix<D, T>::invert()
{
  //Load into GSL
  gsl_matrix *matrixToInvert = gsl_matrix_alloc(D,D);
  for   (int i=0; i<D; i++)
    for (int j=0; j<D; j++)
      gsl_matrix_set(matrixToInvert,i,j,x[i][j]);


  //Computations
  // variables needed
  int matrixSize = matrixToInvert->size1;
  int permInt = 0;                                             		     // necessary for gsl inverse calculation
  int successLU = 0;                                           		     // returns the success of the LU decomposition and inverse calc
  gsl_permutation *inverseP = gsl_permutation_calloc(matrixSize);            // necessary for gsl inverse calculation
  gsl_matrix *LUdecomp = gsl_matrix_calloc(matrixSize,matrixSize);           // necessary for the gsl inverse calculation
  gsl_matrix *inverseMatrix = gsl_matrix_calloc(matrixSize,matrixSize);      // the inverse matrix

  // calculate inverse
  gsl_matrix_memcpy(LUdecomp,matrixToInvert);

  successLU = gsl_linalg_LU_decomp(LUdecomp,inverseP,&permInt);
  if (successLU == 1) {
    std::cerr << "Inverse calculation failed at LU decomposition\n";
    for   (unsigned int i=0; i<matrixToInvert->size1; i++) {
      for (unsigned int j=0; j<matrixToInvert->size2; j++) {
	std::cerr<<"  "<<gsl_matrix_get(matrixToInvert,i,j);
      }
      std::cerr<<std::endl;
    }
    throw 1;
  }
  successLU = gsl_linalg_LU_invert(LUdecomp,inverseP,inverseMatrix);
  if (successLU == 1) {
    std::cerr << "Inverse calculation of the lambdas matrix failed at inverse calculation\n";
    for   (unsigned int i=0; i<matrixToInvert->size1; i++) {
      for (unsigned int j=0; j<matrixToInvert->size2; j++) {
	std::cerr<<"  "<<gsl_matrix_get(matrixToInvert,i,j);
      }
      std::cerr<<std::endl;
    }
    throw 1;
  }

  // memory free
  gsl_matrix_free(LUdecomp);
  gsl_permutation_free(inverseP);
  gsl_matrix_free(matrixToInvert);

  //Export from GSL
  for   (int i=0; i<D; i++)
    for (int j=0; j<D; j++)
      x[i][j] = gsl_matrix_get(inverseMatrix,i,j);
  gsl_matrix_free(inverseMatrix);

  return *this;
}



// transpose()
// ~~~~~~~~~~~
template <int D,typename T>
inline euclidean_matrix<D, T>& euclidean_matrix<D, T>::transpose()
{
  for   (int i=0  ; i<D; i++) {
    for (int j=i+1; j<D; j++) {
      T temp = x[i][j];
      x[i][j] = x[j][i];
      x[j][i] = temp;
    }
  }

  return *this;
}



// get_transpose()
// ~~~~~~~~~~~~~~~
template <int D,typename T>
inline euclidean_matrix<D, T> euclidean_matrix<D, T>::get_transpose() const
{
  euclidean_matrix<D,T> m(*this);
  m.transpose();

  return m;
}



// +
// ~ 
template <int D, typename T>
inline euclidean_matrix<D, T> euclidean_matrix<D, T>::operator+(const euclidean_matrix<D, T> &a) const
{
  euclidean_matrix<D, T> c;
  
  for   (int i=0; i<D; i++)
    for (int j=0; j<D; j++)
    c.x[i][j] = x[i][j] + a.x[i][j];

  return c;
}



// -
// ~
template <int D, typename T>
inline euclidean_matrix<D, T> euclidean_matrix<D, T>::operator-(const euclidean_matrix<D, T> &a) const
{
  euclidean_matrix<D, T> c;

  for   (int i=0; i<D; i++)
    for (int j=0; j<D; j++)
    c.x[i][j] = x[i][j] - a.x[i][j];

  return c;
}



// *
// ~
// This is matrix-matrix multiplication (not entry-by-entry)!
template <int D, typename T>
inline euclidean_matrix<D, T> euclidean_matrix<D, T>::operator*(const euclidean_matrix<D, T> &a) const
{
  euclidean_matrix<D, T> c;
  
  for     (int i=0; i<D; i++)
    for   (int j=0; j<D; j++)
      for (int k=0; k<D; k++)
	c[i][j] += x[i][k] * a.x[k][j];

  return c;
}



// *
// ~
// This is matrix multiplication with a vector
//  (matrix on left, vector on right)
template <int D, typename T>
inline euclidean_vector<D, T> euclidean_matrix<D, T>::operator*(const euclidean_vector<D, T> &a) const
{
  euclidean_vector<D, T> c;
  
  for   (int i=0; i<D; i++)
    for (int k=0; k<D; k++)
      c.x[i] += x[i][k] * a.x[k];

  return c;
}



// /
// ~
template <int D, typename T>
inline euclidean_matrix<D, T> euclidean_matrix<D, T>::operator/(const T s) const
{
  euclidean_matrix<D, T> c;

  for   (int i=0; i<D; i++)
    for (int j=0; j<D; j++)
    c[i][j] = x[i][j] / s;

  return c;
}



// ==
// ~
template <int D, typename T>
inline bool euclidean_matrix<D, T>::operator==(const euclidean_matrix<D, T> &a) const
{
  for   (int i=0; i<D; i++)
    for (int j=0; j<D; j++)
    {
      if (!(x[i][j]==a.x[i][j]))
	return false;
    }
  return true;
}



// %
// ~
template <int D, typename T>
inline euclidean_matrix<D, T> euclidean_matrix<D, T>::operator%(const T s) const
{
  euclidean_matrix<D, T> c;

  for   (int i=0; i<D; i++)
    for (int j=0; j<D; j++)
    c[i][j] = x[i][j] % s;

  return c;
}



// integer
// ~~~~~~~
template <int D, typename T>
inline euclidean_matrix<D, int> euclidean_matrix<D, T>::Integer() const
{
  euclidean_matrix<D, int> c;

  for   (int i=0; i<D; i++)
    for (int j=0; j<D; j++)
    c[i][j] = (int)x[i][j];

  return c;
}



template <int D, typename T>
inline euclidean_matrix<D, int> euclidean_matrix<D, T>::Integer(const euclidean_matrix<D, T>& m)
{
  return m.Integer();
}



// double
// ~~~~~~~
template <int D, typename T>
inline euclidean_matrix<D, double> euclidean_matrix<D, T>::Double() const
{
  euclidean_matrix<D, double> c;

  for   (int i=0; i<D; i++)
    for (int j=0; j<D; j++)
    c[i][j] = (double)x[i][j];

  return c;
}



template <int D, typename T>
inline euclidean_matrix<D, double> euclidean_matrix<D, T>::Double(const euclidean_matrix<D, T>& m)
{
  return m.Double();
}



// []
// ~~
template <int D, typename T>
inline T* euclidean_matrix<D, T>::operator[](const unsigned int i)
{
  return x[i];
}



// dot
// ~~~
template <int D, typename T>
inline double euclidean_matrix<D, T>::dot(const euclidean_matrix<D, T> &a) const
{
  double d=0;

  for   (int i=0; i<D; i++)
    for (int j=0; j<D; j++)
      d += x[i][j] * a.x[i][j];

  return d;
}



template <int D, typename T>
inline double euclidean_matrix<D, T>::dot(const euclidean_matrix<D, T> &a, const euclidean_matrix<D, T> &b)
{
  return a.dot(b);
}



// norm_squared
// ~~~~~~~~~~~
template <int D, typename T>
inline double euclidean_matrix<D, T>::norm_sq() const
{
  return dot(*this, *this);
}



template <int D, typename T>
inline double euclidean_matrix<D, T>::norm_sq(const euclidean_matrix<D, T>& m)
{
  return m.norm_sq();
}



// determinant
// ~~~~~~~~~~~

//Note: This should be done without GSL in the future if possible to prevent heap fragmentation if this is called frequently

template <int D, typename T>
inline T euclidean_matrix<D, T>::determinant() const
{
  switch(D) {
  case 2 : 
    return x[0][0]*x[1][1] - x[0][1]*x[1][0];
    break;
  case 3 :
    return 
        x[0][0]*(x[1][1]*x[2][2] - x[1][2]*x[2][1])
      - x[0][1]*(x[1][0]*x[2][2] - x[1][2]*x[2][0])
      + x[0][2]*(x[1][0]*x[2][1] - x[1][1]*x[2][0]);
    break;
  default : //Long (but slow) version for any size
    //Check squareness of the matrix:
    gsl_matrix* A = gsl_matrix_alloc(D,D);

    for (int i=0; i<D; i++)
      for (int j=0; j<D; j++)
	gsl_matrix_set(A,i,j,x[i][j]);

    //LU-decompose the copy
    gsl_permutation *perm = gsl_permutation_calloc(D);
    int permInt = 0; //Not sure why I need this, but whatever.

    gsl_linalg_LU_decomp(A,perm,&permInt);
    const double det = gsl_linalg_LU_det(A,1);

    //Deallocate:
    gsl_matrix_free(A);
    gsl_permutation_free(perm);

    return det;
  }
}



template <int D, typename T>
inline T euclidean_matrix<D, T>::determinant(const euclidean_matrix<D, T>& m)
{
  return m.determinant();
}



// read
// ~~~~
template <int D, typename T>
void euclidean_matrix<D, T>::read(std::ifstream& in)
{
  in.read((char*)x, sizeof(T)*D);
}



// write
// ~~~~~
template <int D, typename T>
void euclidean_matrix<D, T>::write(std::ofstream& out) const
{
  out.write((const char*)x, sizeof(T)*D);
}



// Insertion
// ~~~~~~~~~
template <int D, typename T>
std::ostream& operator<<(std::ostream& os, const euclidean_matrix<D, T>& v)
{
  for (int i=0; i<D; i++) {
    os << "[";
    for(int j=0; j<D-1; j++)
      os << v.x[i][j] << ", ";
    os << v.x[i][D-1] << "]";

    if (i<D-1)
      os<<"\n";
  }
  
  return os;
}

#endif 
