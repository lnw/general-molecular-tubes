
#include <cassert>
#include <cmath>
#include <vector>

#include "auxiliary.hh"
#include "geometry3.hh"
#include "tile.hh"

using namespace std;


structure3d Tile::get_struct_default() const {
  structure3d S;
  for (size_t i = 0; i < indices_default.size(); i++) {
    const size_t index = indices_default[i];
    S.push_back(struct3d.dat[index]);
    S.push_back_atom(struct3d.atomtypes[index]);
  }
  return S;
}

structure3d Tile::get_struct_first() const {
  structure3d S;
  for (size_t i = 0; i < indices_first.size(); i++) {
    const size_t index = indices_first[i];
    S.push_back(struct3d.dat[index]);
    S.push_back_atom(struct3d.atomtypes[index]);
  }
  return S;
}

structure3d Tile::get_struct_last() const {
  structure3d S;
  for (size_t i = 0; i < indices_last.size(); i++) {
    const size_t index = indices_last[i];
    S.push_back(struct3d.dat[index]);
    S.push_back_atom(struct3d.atomtypes[index]);
  }
  return S;
}
