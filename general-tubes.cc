
#include <cassert>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>

#include "auxiliary.hh"
// #include "geometry2.hh"
#include "geometry3.hh"
#include "tile.hh"
#include "shape-gen.hh"
#include "tile-gen.hh"

using namespace std;



int main(int ac, char **av) {

//  cerr << ac << endl;
  if (ac != 5) {
    cout << "general-tube-gen" << endl;
    cout << "usage: " << av[0] << " <n_db> <n_rings> <shape> <tile>" << endl;
    cout << "  <n_db>: the number of double bonds" << endl;
    cout << "  <n_rings>: ..." << endl;
    cout << "  <shape>: ..." << endl;
    cout << "  <tile>: ..." << endl;
    abort();
  }
  const int n_db = stol(av[1], 0, 0);
  const int n_rings = stol(av[2], 0, 0);
  const int c = stol(av[3], 0, 0);
  const int t = stol(av[4], 0, 0);

  cout << "number of double bonds: " << n_db << endl;
  cout << "number of rings: " << n_rings << endl;
  cout << "shape type: " << c << endl;
  cout << "tile type: " << t << endl;

  TileGenerator tile_gen;
  auto tile = tile_gen.get_tile(t);

  const structure3d str_f(tile.get_struct_first());
  const structure3d str_d(tile.get_struct_default());
  const structure3d str_l(tile.get_struct_last());
  // cout << str_f << endl;
  // cout << str_d << endl;
  // cout << str_l << endl;

  // get r and R
  const double unit_cell_x = tile.get_width();
  const double unit_cell_y = tile.get_height();
  const double r = unit_cell_x*n_db / (2*M_PI);
  const double R = unit_cell_y*n_rings / (2*M_PI);
  cout << "R (major radius, length): " << R << endl;
  cout << "r (minor radius): " << r << endl;

  ShapeGenerator shape_gen;
  auto shape = shape_gen.get_shape(c);
  auto shape_tan = shape_gen.get_shape_tan(c);
  bool is_closed_curve = shape_gen.is_closed_curve(c);

  structure3d S2;
  for (int ring=0; ring<n_rings; ring++){
    structure3d str_;
    if(is_closed_curve)
      str_ = str_d;
    else if(ring==0)
      str_ = str_f;
    else if(ring==n_rings-1)
      str_ = str_l;
    else
      str_ = str_d;
     
    double ignore_me;
    const double shift = modf( tile.get_shift() * ring, &ignore_me);

    for (int db=0; db<n_db; db++){
      for(size_t atom=0; atom<str_.size(); atom++){
        S2.push_back(str_[atom] + coord3d(unit_cell_x*(db + shift), unit_cell_y*ring, 0.0));
        S2.push_back_atom(str_.atomtypes[atom]);
      }
    }
  }
  // cout << S2 << endl;

  const int n_steps(1000);
  double length_tot;
  const vector<double> dists_inv_acc = shape_gen.get_inverse_mapping(c, n_steps, length_tot);
  // cout << dists_inv_acc << endl;

  const double factor = 2* M_PI * R / length_tot;
  // cout << "fac: " << factor << endl;

  structure3d S3;
  for(size_t atom=0; atom<S2.size(); atom++){
    const double phi  ( (S2[atom][0])/r );
    const double theta( (S2[atom][1])/R );
    // cout << "phi, theta: " << theta << ", " << phi << endl;
    // modulate theta such that a tube is spaced equidistantly along its axis 
    double theta_prime;
    if(theta < 0 || theta > 2*M_PI){
      theta_prime = theta;
    }
    else {
      if ( abs(theta*n_steps/(2*M_PI) - round( theta*n_steps/(2*M_PI) )) < 1.e-3){ // don't interpolate
        const int index = int(theta*n_steps/(2*M_PI) + 0.5);
        // cout << "a " << index << endl;
        theta_prime = dists_inv_acc[index];
      }
      else { // interpolate
        const int index1 = int(theta*n_steps/(2*M_PI));
        const double weight2 = theta*n_steps/(2*M_PI) - index1;
        const int index2 = index1 + 1;
        const double weight1 = index2 - theta*n_steps/(2*M_PI);
        // cout << index1 << ", " << index2 << endl;
        // cout << weight1 << ", " << weight2 << endl;
        // cout << "b1 " << index1 << endl;
        // cout << "b2 " << index2 << endl;
        theta_prime = weight1 * dists_inv_acc[index1] + weight2 * dists_inv_acc[index2];
// cout << "t, t': {" << theta << ", " << theta_prime << "}," << endl;
      }
    }

    const coord3d nn(shape(theta_prime) * factor); // position vector back bone
    const coord3d tangent(shape_tan(theta_prime)); // tangent vector back bone
    const coord3d ref(0,0,1);
    const coord3d va( (tangent.cross(ref)).normalised() );
    const coord3d vb( (tangent.cross(va)).normalised() );
    const double rad = r+S2[atom][2];
    const coord3d c3d( nn + va*rad*cos(phi) + vb*rad*sin(phi) );
    // cout << c3d << endl;
    S3.push_back(c3d);
    S3.push_back_atom(S2.atomtypes[atom]);
  }

// cout << S3 << endl;


  ofstream xyz(("whatever-" + to_string(n_db) +"-"+ to_string(n_rings) + "-" + to_string(c) + "-" + to_string(t) + ".xyz").c_str());
  ofstream turbo(("whatever-" + to_string(n_db) +"-"+ to_string(n_rings) + "-" + to_string(c) + "-" + to_string(t) + ".coord").c_str());
  xyz << S3.to_xyz();
  turbo << S3.to_turbomole();
  xyz.close();
  turbo.close();

  return 0;
}


