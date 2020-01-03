#ifndef TILE_HH
#define TILE_HH

#include <cassert>
#include <cmath>
#include <vector>

#include "auxiliary.hh"
#include "geometry3.hh"

using namespace std;


class Tile {

  double width, height; // for the default tile
  double shift;         // if row 0 has no x shift, row 1 is shifted right by 'shift'. 0 <= shift < 1. Typical values are 0, 0.333, 0.5 etc.
  structure3d struct3d;
  vector<size_t> indices_default;
  vector<size_t> indices_first;
  vector<size_t> indices_last;

public:
  Tile(double w, double h, double s, const structure3d& str, const vector<size_t>& _id, const vector<size_t>& _if, const vector<size_t>& _il): width(w), height(h), shift(s), struct3d(str), indices_default(_id), indices_first(_if), indices_last(_il) {}

  double get_height() const { return height; }
  double get_width() const { return width; }
  double get_shift() const { return shift; }

  size_t get_n_default() const { return indices_default.size(); }
  size_t get_n_first() const { return indices_first.size(); }
  size_t get_n_last() const { return indices_last.size(); }

  structure3d get_struct_default() const;
  structure3d get_struct_first() const;
  structure3d get_struct_last() const;

  //   friend ostream& operator<<(ostream& S, const structure2d &s2) {
  //     S << s2.dat;
  //     return S;
  //   }
};

#endif
