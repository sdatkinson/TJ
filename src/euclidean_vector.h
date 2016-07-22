#ifndef EUCLIDEAN_VECTOR_H
#define EUCLIDEAN_VECTOR_H

#include <iostream>
#include <fstream>
#include <assert.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <cmath>

template <int D, typename T=double>
class euclidean_vector {
  
  public:
  T x[D];
  
  public:
  
  euclidean_vector();
  euclidean_vector(const T); //Set all entries equal to this value
  euclidean_vector(const T[D]);           
  euclidean_vector(const euclidean_vector&);
  euclidean_vector(const gsl_vector*);
  euclidean_vector(const gsl_matrix*);
  ~euclidean_vector();

  euclidean_vector<D, T>& operator= (const euclidean_vector<D, T>&);
  euclidean_vector<D, T>& operator= (const gsl_matrix*);
  euclidean_vector<D, T>& operator+=(const euclidean_vector<D, T>&);
  euclidean_vector<D, T>& operator-=(const euclidean_vector<D, T>&);
  euclidean_vector<D, T>& operator*=(const T);
  euclidean_vector<D, T>& operator/=(const T);
  euclidean_vector<D, T>& operator%=(const double);

  euclidean_vector<D, T>& Increment(const T,  //min
				   const T); //max
  euclidean_vector<D, T>& Increment(const T,  //min
				   const T,  //inc
				   const T); //max
  euclidean_vector<D, T>& Increment(const euclidean_vector<D, T>&,  //min
				   const euclidean_vector<D, T>&,  //inc
				   const euclidean_vector<D, T>&); //max

  euclidean_vector<D, T> operator+(const euclidean_vector<D, T>&) const;
  euclidean_vector<D, T> operator+(const T) const;
  euclidean_vector<D, T> operator-(const euclidean_vector<D, T>&) const;
  euclidean_vector<D, T> operator-(const T) const;
  euclidean_vector<D, T> operator*(const T) const;
  //euclidean_vector<D, T> operator*(const euclidean_matrix<D, T>&) const; //x*A (w/ matrix)
  euclidean_vector<D, T> operator/(const T) const;
  euclidean_vector<D, T> operator%(const int   ) const;
  euclidean_vector<D, T> operator%(const double) const;
  euclidean_vector<D, T> operator%(const euclidean_vector<D,int   >&) const;
  euclidean_vector<D, T> operator%(const euclidean_vector<D,double>&) const;
  bool operator==(const euclidean_vector<D, T> &a) const;

  euclidean_vector       <D, int>    Integer() const;
  static euclidean_vector<D, int>    Integer(const euclidean_vector<D, T>&); 
  euclidean_vector       <D, int>    Integer_Round() const;
  static euclidean_vector<D, int>    Integer_Round(const euclidean_vector<D, T>&); 
  euclidean_vector       <D, double> Double () const;
  static euclidean_vector<D, double> Double (const euclidean_vector<D, T>&); 
  euclidean_vector       <D, T>      Abs    () const;
  static euclidean_vector<D, T>      Abs    (const euclidean_vector<D, T>&); 
  T                                  Max    () const;
  static T                           Max    (const euclidean_vector<D, T>&);
  static euclidean_vector<D, double> rand_sphere(); //Create a random vector within the unit sphere
  
  T& operator[](const unsigned int);

  void CopyFrom(const euclidean_vector<D, T>&);

  void get_minimum_image();
  void get_minimum_image(euclidean_vector<D, int>&);
 
  double        dot(const euclidean_vector<D, T>&) const;
  static double dot(const euclidean_vector<D, T>&, const euclidean_vector<D, T>&);
  
  double        norm_sq() const;
  static double norm_sq(const euclidean_vector<D, T>&);

  double        norm() const;
  static double norm(const euclidean_vector<D, T>&);

  euclidean_vector<D, T>& normalize();

  void read(std::ifstream&);
  void write(std::ofstream&) const;
};

template <int D, typename T>
std::ostream& operator<<(std::ostream&, const euclidean_vector<D, T>&);




//==================================================================
// Implementation



// constructor
// ~~~~~~~~~~~
template <int D, typename T>
euclidean_vector<D, T>::euclidean_vector()
{
  for(int k=0; k<D; k++)
    x[k] = 0;
}



template <int D, typename T>
euclidean_vector<D, T>::euclidean_vector(const T x_i)
{
  for (int k=0; k<D; k++)
    x[k]=x_i;
}



template <int D, typename T>
euclidean_vector<D, T>::euclidean_vector(const T x_i[D])
{
  for(int k=0; k<D; k++)
    x[k] = x_i[k];
}



template <int D, typename T>
euclidean_vector<D, T>::euclidean_vector(const euclidean_vector<D, T> &v)
{
  for(int k=0; k<D; k++)
    x[k] = v.x[k];
}



template <int D, typename T>
euclidean_vector<D, T>::euclidean_vector(const gsl_vector* v)
{
  assert(v->size==D);

  for(int k=0; k<D; k++)
    x[k] = gsl_vector_get(v,k);
}



template <int D, typename T>
euclidean_vector<D, T>::euclidean_vector(const gsl_matrix* v)
{
  assert(v->size1==D && v->size2==1);

  for(int k=0; k<D; k++)
    x[k] = gsl_matrix_get(v,k,0);
}



// destructor
// ~~~~~~~~~~
template <int D, typename T>
euclidean_vector<D, T>::~euclidean_vector()
{
}



// =
// ~~
template <int D, typename T>
inline euclidean_vector<D, T>& euclidean_vector<D, T>::operator=(const euclidean_vector<D, T> &v)
{
  for(int k=0; k<D; k++)
    x[k] = v.x[k];

  return *this;
}



template <int D, typename T>
inline euclidean_vector<D, T>& euclidean_vector<D, T>::operator=(const gsl_matrix* v)
{
  assert(v->size1==D && v->size2==1);
  for(int k=0; k<D; k++)
    x[k] = (T) gsl_matrix_get(v,k,0);

  return *this;
}



// +=
// ~~
template <int D, typename T>
inline euclidean_vector<D, T>& euclidean_vector<D, T>::operator+=(const euclidean_vector<D, T> &v)
{
  for(int k=0; k<D; k++)
    x[k] += v.x[k];

  return *this;
}



// -=
// ~~
template <int D, typename T>
inline euclidean_vector<D, T>& euclidean_vector<D, T>::operator-=(const euclidean_vector<D, T> &v)
{
  for(int k=0; k<D; k++)
    x[k] -= v.x[k];

  return *this;
}



// *=
// ~~
template <int D, typename T>
inline euclidean_vector<D, T>& euclidean_vector<D, T>::operator*=(const T s)
{
  for(int k=0; k<D; k++)
    x[k] *= s;

  return *this;
}



// /=
// ~~
template <int D, typename T>
inline euclidean_vector<D, T>& euclidean_vector<D, T>::operator/=(const T s)
{
  for(int k=0; k<D; k++)
    x[k] /= s;

  return *this;
}



// %=
// ~~
template <int D, typename T>
inline euclidean_vector<D, T>& euclidean_vector<D, T>::operator%=(const double s)
{
  for(int k=0; k<D; k++)
    x[k] = std::fmod(x[k],s);

  return *this;
}



// Increment
// ~~
template <int D, typename T>
inline euclidean_vector<D, T>& euclidean_vector<D, T>::Increment(const T minVal,
							       const T maxVal)
{
  this->Increment(minVal,(T) 1,maxVal);
  return *this;
}

template <int D, typename T>
inline euclidean_vector<D, T>& euclidean_vector<D, T>::Increment(const T minVal,
							       const T incVal,
							       const T maxVal)
{
  for (int k=0; k<D; k++) {
    x[k]+=incVal;
    if (x[k]<=maxVal)
      break;
    x[k]=minVal;
  }

  return *this;
}



template <int D, typename T>
inline euclidean_vector<D, T>& euclidean_vector<D, T>::Increment(const euclidean_vector<D, T> &minVal,
							       const euclidean_vector<D, T> &incVal,
							       const euclidean_vector<D, T> &maxVal)
{
  for (int k=0; k<D; k++) {
    x[k]+=incVal.x[k];
    if (x[k]<=maxVal.x[k])
      break;
    x[k]=minVal.x[k];
  }

  return *this;
}



// +
// ~ 
template <int D, typename T>
inline euclidean_vector<D, T> euclidean_vector<D, T>::operator+(const euclidean_vector<D, T> &a) const
{
  euclidean_vector<D, T> c;
  
  for(int k=0; k<D; k++)
    c.x[k] = x[k] + a.x[k];

  return c;
}



template <int D, typename T>
inline euclidean_vector<D, T> euclidean_vector<D, T>::operator+(const T s) const
{
  euclidean_vector<D, T> c;
  
  for(int k=0; k<D; k++)
    c.x[k] = x[k] + s;

  return c;
}



// -
// ~
template <int D, typename T>
inline euclidean_vector<D, T> euclidean_vector<D, T>::operator-(const euclidean_vector<D, T> &a) const
{
  euclidean_vector<D, T> c;

  for(int k=0; k<D; k++)
    c[k] = x[k] - a.x[k];

  return c;
}



template <int D, typename T>
inline euclidean_vector<D, T> euclidean_vector<D, T>::operator-(const T s) const
{
  euclidean_vector<D, T> c;
  
  for(int k=0; k<D; k++)
    c.x[k] = x[k] - s;

  return c;
}



// *
// ~
template <int D, typename T>
inline euclidean_vector<D, T> euclidean_vector<D, T>::operator*(const T s) const
{
  euclidean_vector<D, T> c;
  
  for(int k=0; k<D; k++)
    c[k] = x[k] * s;

  return c;
}



// /
// ~
template <int D, typename T>
inline euclidean_vector<D, T> euclidean_vector<D, T>::operator/(const T s) const
{
  euclidean_vector<D, T> c;

  for(int k=0; k<D; k++)
    c[k] = x[k] / s;

  return c;
}



// ==
// ~
template <int D, typename T>
inline bool euclidean_vector<D, T>::operator==(const euclidean_vector<D, T> &a) const
{
  for(int k=0; k<D; k++)
    {
      if (!(x[k]==a.x[k]))
	return false;
    }
  return true;
}



// %
// ~
template <int D, typename T>
inline euclidean_vector<D, T> euclidean_vector<D, T>::operator%(const int s) const
{
  euclidean_vector<D, T> c;

  for(int k=0; k<D; k++)
    c[k] = (int)x[k] % s;

  return c;
}



template <int D, typename T>
inline euclidean_vector<D, T> euclidean_vector<D, T>::operator%(const double s) const
{
  euclidean_vector<D, T> c;

  for(int k=0; k<D; k++)
    c[k] = fmod(x[k],s);

  return c;
}



template <int D, typename T>
  inline euclidean_vector<D, T> euclidean_vector<D, T>::operator%(const euclidean_vector<D,int>& a) const
{
  euclidean_vector<D, T> c;
  euclidean_vector<D,int> b(a);

  for(int k=0; k<D; k++)
    c[k] = (int)x[k] % b[k];

  return c;
}



template <int D, typename T>
  inline euclidean_vector<D, T> euclidean_vector<D, T>::operator%(const euclidean_vector<D,double>& a) const
{
  euclidean_vector<D, T> c;
  euclidean_vector<D,int> b(a);

  for(int k=0; k<D; k++)
    c[k] = fmod(x[k],b[k]);

  return c;
}



// Integer
// ~~~~~~~
template <int D, typename T>
inline euclidean_vector<D, int> euclidean_vector<D, T>::Integer() const
{
  euclidean_vector<D, int> c;

  for(int k=0; k<D; k++)
    c[k] = (int)x[k];

  return c;
}



template <int D, typename T>
inline euclidean_vector<D, int> euclidean_vector<D, T>::Integer(const euclidean_vector<D, T>& v)
{
  return v.Integer();
}



// Integer_Round
// ~~~~~~~~~~~~~
template <int D, typename T>
inline euclidean_vector<D, int> euclidean_vector<D, T>::Integer_Round() const
{
  euclidean_vector<D, int> c;

  for(int k=0; k<D; k++)
    c[k] = int(round(x[k]));

  return c;
}



template <int D, typename T>
inline euclidean_vector<D, int> euclidean_vector<D, T>::Integer_Round(const euclidean_vector<D, T>& v)
{
  return v.Integer_Round();
}



// Double
// ~~~~~~~
template <int D, typename T>
inline euclidean_vector<D, double> euclidean_vector<D, T>::Double() const
{
  euclidean_vector<D, double> c;

  for(int k=0; k<D; k++)
    c[k] = (double)x[k];

  return c;
}



template <int D, typename T>
inline euclidean_vector<D, double> euclidean_vector<D, T>::Double(const euclidean_vector<D, T>& v)
{
  return v.Double();
}



// abs
// ~~~~~~~
template <int D, typename T>
inline euclidean_vector<D, T> euclidean_vector<D, T>::Abs() const
{
  euclidean_vector<D, T> c;

  for(int k=0; k<D; k++)
    c[k] = std::abs(x[k]);

  return c;
}



template <int D, typename T>
inline euclidean_vector<D, T> euclidean_vector<D, T>::Abs(const euclidean_vector<D, T>& v)
{
  return v.Abs();
}


// max
// ~~~~~~~
template <int D, typename T>
inline T euclidean_vector<D, T>::Max() const
{
  T maxVal=x[0];

  for(int k=1; k<D; k++)
    if (x[k]>maxVal)
      maxVal=x[k];

  return maxVal;
}

template <int D, typename T>
inline T euclidean_vector<D, T>::Max(const euclidean_vector<D, T>& v)
{
  return v.Max();
}



// []
// ~~
template <int D, typename T>
inline T& euclidean_vector<D, T>::operator[](const unsigned int i)
{
  return x[i];
}



// CopyFrom()
//  Assignment in which "this" doesn't change
// ~~~~~~~~~~
template <int D, typename T>
  inline void euclidean_vector<D, T>::CopyFrom(const euclidean_vector<D, T>& v)
{
  for (int d=0; d<DIM; d++)
    x[d] = v.x[d];
}


// Get Minimum Image (L/2 method)
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template <int D, typename T>
inline void euclidean_vector<D, T>::get_minimum_image()
{
  for (int d=0; d<D; d++)
    x[d] = fmod(x[d]+1.5,1.0)-0.5;
}



template <int D, typename T>
inline void euclidean_vector<D, T>::get_minimum_image(euclidean_vector<D,int>& image)
// Goes a little slower since we need to keep track of what changes actually happened
{
  euclidean_vector<D,T> v0(*this);
  for (int d=0; d<D; d++)
    x[d] = fmod(x[d]+1.5,1.0)-0.5;
  image = euclidean_vector<D,T>::Integer_Round((*this) - v0);
}



// dot
// ~~~
template <int D, typename T>
inline double euclidean_vector<D, T>::dot(const euclidean_vector<D, T> &a) const
{
  double d=0;

  for(int k=0; k<D; k++)
    d += x[k] * a.x[k];

  return d;
}



template <int D, typename T>
inline double euclidean_vector<D, T>::dot(const euclidean_vector<D, T> &a, const euclidean_vector<D, T> &b)
{
  return a.dot(b);
}


// norm Squared
// ~~~~~~~~~~~~
template <int D, typename T>
inline double euclidean_vector<D, T>::norm_sq() const
{
  return dot(*this, *this);
}

template <int D, typename T>
inline double euclidean_vector<D, T>::norm_sq(const euclidean_vector<D, T>& v)
{
  return v.norm_sq();
}



// norm
// ~~~~
template <int D, typename T>
inline double euclidean_vector<D, T>::norm() const
{

  return sqrt(norm_sq());
}



template <int D, typename T>
inline double euclidean_vector<D, T>::norm(const euclidean_vector<D, T>& v)
{
  return v.norm();
}


// normalize
// ~~~~~~~~~
template <int D, typename T>
inline euclidean_vector<D, T>& euclidean_vector<D, T>::normalize()
{
  this/=this.norm();
  return *this;
}



// read
// ~~~~
template <int D, typename T>
void euclidean_vector<D, T>::read(std::ifstream& in)
{
  in.read((char*)x, sizeof(T)*D);
}



// write
// ~~~~~
template <int D, typename T>
void euclidean_vector<D, T>::write(std::ofstream& out) const
{
  out.write((const char*)x, sizeof(T)*D);
}



// Insertion
// ~~~~~~~~~
template <int D, typename T>
std::ostream& operator<<(std::ostream& os, const euclidean_vector<D, T>& v)
{
  os << "(";

  for(int k=0; k<D-1; k++)
    os << v.x[k] << ", ";

  os << v.x[D-1] << ")";

  return os;
}


#endif 
