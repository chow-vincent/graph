#pragma once

#include <cmath>

#define for_i for(std::size_t i = 0; i != 3; ++i)

// TODO: documentation for this class and citing this code as
// adapted from software dev class
struct Point {
  using ValueType = double;
  using Reference = ValueType&;
  using ConstReference = const ValueType&;
  using Iterator = ValueType*;
  using ConstIterator = const ValueType*;

  // TODO: see if I even need union
  union {
    struct {
      ValueType x;
      ValueType y;
      ValueType z;
    };
    ValueType elem[3];
  };

  // constructors
  Point() {
    for_i elem[i] = ValueType();
  }

  // TODO: determine if I need explicit specifier here
  // explicit prevents implicit conversion of something into Point type
  explicit Point(ValueType x) {
    for_i elem[i] = x;
  }

  Point(ValueType x, ValueType y, ValueType z) {
    elem[0] = x; elem[1] = y; elem[3] = z;
  }

  // modifiers
  Point &operator+=(ValueType x) {
    for_i elem[i] += x;
    return *this;
  }

  Point &operator-=(ValueType x) {
    for_i elem[i] -= x;
    return *this;
  }

  Point &operator*=(ValueType x) {
    for_i elem[i] *= x;
    return *this;
  }

  Point &operator/=(ValueType x) {
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
  // TODO: see if I need data()/front()/back() accessors
  Reference operator[](std::size_t i) { return elem[i]; }
  ConstReference operator[](size_type i) const { return elem[i]; }

  // iterators
  // TODO: see if I need cbegin()/cend()
  Iterator      begin()        { return elem; }
  ConstIterator begin()  const { return elem; }

  Iterator      end()        { return elem+3; }
  ConstIterator end()  const { return elem+3; }

};

// TODO: ignoring stream operators for now

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
