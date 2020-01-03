#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>

#include "auxiliary.hh"

using namespace std;

// read atom selection from file, one line, starting at 0
vector<int> read_atom_selection(const char* input){
  ifstream file(input);
  string indices_string;
  getline(file, indices_string);
  stringstream indices_stream(indices_string);
  int index;
  vector<int> selection;
  while(indices_stream >> index){
    selection.push_back(index);
  }
  file.close();
  return selection;
}


