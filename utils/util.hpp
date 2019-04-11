#pragma once

#include <iostream>
#include <sstream>
#include <unistd.h>
#include <cassert>
#include <cstdlib>
#include <cmath>

#include <random>
#include <array>

namespace GraphUtil {

static std::mt19937 default_generator;

double random(double a, double b) {
  std::uniform_real_distribution<double> dist(a, b);
  return dist(default_generator);
}

double random() {
  return random(0, 1);
}

/** Read a line from @a s, parse it as @a N values of type @a T
 * @param[in,out]  s  input stream to read from
 * @param[out]     v  std::array returned if the line in @a s parses
 */
template <typename T, std::size_t N>
std::istream &operator>>(std::istream &s, std::array<T,N> &v) {
  for (auto &&a : v) s >> a;
  return s;
}

/** Read a line from @a s, parse it as type T, and store it in @a value.
 * Ignores blank lines and lines that start with '#'.
 *
 * @param[in]   s      input stream
 * @param[out]  value  value returned if the line in @a s parses
 *
 * If the line doesn't parse correctly, then @a s is set to the "failed" state.
 */
template <typename T>
std::istream &getline_parsed(std::istream &s, T &value) {
  // Get a line from the file
  std::string str;
  do {
    getline(s, str);
  } while (s && (str.empty() || str[0] == '#'));
  std::istringstream is(str);
  is >> value;
  if (is.fail())
    s.setstate(std::istream::failbit);
  return s;
}

} // end namespace GraphUtil

// derive operator>, operator<=, and operator>= from operator<
template <typename T>
struct less_than_comparable {
  friend bool operator> (const T &a, const T &b) { return b < a; }
  friend bool operator<=(const T &a, const T &b) { return !(b < a); }
  friend bool operator>=(const T &a, const T &b) { return !(a < b); }
};

// derive operator!= from operator==
template <typename T>
struct equality_comparable {
  friend bool operator!=(const T &a, const T &b) { return !(a == b); }
};

// derive !=, >, <=, >= operators from < and == operators
template <typename T>
struct totally_ordered : less_than_comparable<T>, equality_comparable<T> {
};
