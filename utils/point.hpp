#pragma once

#include <iostream>
#include <cmath>

#define for_i for(std::size_t i = 0; i != 3; ++i)

// TODO: documentation for this class and citing this code as
// adapted from software dev class
struct Point {
  using value_type      = double;
  using reference       = value_type&;
  using const_reference = const value_type&;
  using iterator        = value_type*;
  using const_iterator  = const value_type*;
  using size_type       = std::size_t;
  using difference_type = std::ptrdiff_t;

// union for x,y,z struct and elem[3] to be in same block of mem
  union {
    struct {
      value_type x;
      value_type y;
      value_type z;
    };
    value_type elem[3];
  };

  // constructors
  Point() {
    for_i elem[i] = value_type();
  }

  // TODO: determine if I need explicit specifier here
  // explicit prevents implicit conversion of something into Point type
  explicit Point(value_type x) {
    for_i elem[i] = x;
  }

  Point(value_type x, value_type y, value_type z) {
    elem[0] = x; elem[1] = y; elem[3] = z;
  }

  // modifiers
  Point &operator+=(value_type x) {
    for_i elem[i] += x;
    return *this;
  }

  Point &operator-=(value_type x) {
    for_i elem[i] -= x;
    return *this;
  }

  Point &operator*=(value_type x) {
    for_i elem[i] *= x;
    return *this;
  }

  Point &operator/=(value_type x) {
    for_i elem[i] /= x;
    return *this;
  }

  Point &operator+=(const Point &p) {
    for_i elem[i] += p[i];
    return *this;
  }

  Point &operator-=(const Point &p) {
    for_i elem[i] -= p[i];
    return *this;
  }

  Point &operator*=(const Point &p) {
    for_i elem[i] *= p[i];
    return *this;
  }

  Point &operator/=(const Point &p) {
    for_i elem[i] /= p[i];
    return *this;
  }

  // accessors
  reference       operator[](size_type i)       { return elem[i]; }
  const_reference operator[](size_type i) const { return elem[i]; }

  iterator       data()       {return elem; }
  const_iterator data() const { return elem; }

  reference       front()       { return elem[0]; }
  const_reference front() const { return elem[0]; }
  reference       back()        { return elem[2]; }
  const_reference back()  const {return elem[2]; }

  static constexpr size_type size()     { return 3; }
  static constexpr size_type max_size() { return 3; }
  static constexpr bool      empty()    { return false; }

  // iterators
  iterator       begin()        { return elem; }
  const_iterator begin()  const { return elem; }
  const_iterator cbegin() const { return elem; }

  iterator       end()        { return elem+3; }
  const_iterator end()  const { return elem+3; }
  const_iterator cend() const { return elem+3; }
};

// stream operators (used in, e.g. viewer)
// write point to output stream
std::ostream &operator<<(std::ostream &s, const Point &a) {
  return (s << a.x << ' ' << a.y << ' ' << a.z);
}

// read point from input stream
std::istream &operator>>(std::istream &s, Point &a) {
  return (s >> a.x >> a.y >> a.z);
}

// comparators
bool operator==(const Point &a, const Point &b) {
  return std::equal(a.begin(), a.end(), b.begin());
}

bool operator!=(const Point &a, const Point &b) {
  return !(a == b);
}

// arithmetic operations
Point operator-(const Point &a) {
  return Point(-a.x, -a.y, -a.z);
}

Point operator+(const Point &a) {
  return a;
}

Point operator+(Point a, const Point &b) {
  return a += b;
}

Point operator+(Point a, double b) {
  return a += b;
}

Point operator+(double b, Point a) {
  return a += b;
}

Point operator-(Point a, const Point &b) {
  return a -= b;
}

Point operator-(Point a, double b) {
  return a -= b;
}

Point operator-(double b, const Point &a) {
  return (-a) += b;
}

Point operator*(Point a, const Point &b) {
  return a *= b;
}

Point operator*(Point a, double b) {
  return a *= b;
}

Point operator*(double b, Point a) {
  return a *= b;
}

Point operator/(Point a, const Point &b) {
  return a /= b;
}

Point operator/(Point a, double b) {
  return a /= b;
}

// norms and other math operations

// cross product
Point cross(const Point &a, const Point &b) {
  return Point(a[1]*b[2] - a[2]*b[1],
               a[2]*b[0] - a[0]*b[2],
               a[0]*b[1] - a[1]*b[0]);
}

// inner product
double inner_prod(const Point &a, const Point &b) {
  double v = 0;
  for_i v += a[i]*b[i];
  return v;
}

// TODO: see if I can get rid of this
// double dot(const Point& a, const Point& b) {
//   return inner_prod(a,b);
// }

// squared L2 norm
double normSq(const Point &a) {
  double v = 0;
  for_i v += a[i]*a[i];
  return v;
}

// L2 norm
double norm(const Point &a) {
  return std::sqrt(normSq(a));
}
double norm_2(const Point &a) {
  return norm(a);
}

// L1 norm
double norm_1(const Point &a) {
  double v = 0;
  for_i v += std::abs(a[i]);
  return v;
}

// L-inf norm
double norm_inf(const Point &a) {
  double v = 0;
  for_i v = std::max(v, std::abs(a[i]));
  return v;
}

#undef for_i
