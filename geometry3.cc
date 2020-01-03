
#include <iomanip>
#include <limits>
#include <cmath>

// #include "geometry2.hh"
#include "geometry3.hh"


coord3d coord3d::operator*(const matrix3d& m) const {
  return coord3d(x[0]*m(0,0)+x[1]*m(1,0)+x[2]*m(2,0),  x[0]*m(0,1)+x[1]*m(1,1)+x[2]*m(2,1),  x[0]*m(0,2)+x[1]*m(1,2)+x[2]*m(2,2));
}

// calculation of the angle beta at b(0,0,0)
double coord3d::angle(const coord3d& a, const coord3d& c) {
  const double L2 = a.dot(a);
  const double R2 = c.dot(c);
  const double M2 = (c-a).dot(c-a);
  //if (abs(M2)<1.e-10) return 0;
  const double den = 2.0*sqrt(L2 * R2);
  double arg = (L2+R2-M2)/den;
  if(arg > 1)  return 0;
  if(arg < -1) return M_PI;
  return acos(arg);
}

// calculation of the derivative of angle beta at b(0,0,0) according to coordinates a and c with fixed b
void coord3d::dangle(const coord3d& a, const coord3d& c, coord3d& da, coord3d& dc) {
  const double L2 = a.dot(a);
  const double R2 = c.dot(c);
  const double M2 = (c-a).dot(c-a);
  const double den = 2.0*sqrt(L2 * R2);
  const double arg = (L2+R2-M2)/den;

  const coord3d dM2__da = (a-c)*2.0;
  const coord3d dL2__da = a*2.0;
  const coord3d dden__da = dL2__da * R2/sqrt(L2*R2);
  const coord3d darg__da = dL2__da * 1.0/den - dM2__da * 1.0/den - dden__da * (L2+R2-M2)/(den*den);

  const coord3d dM2__dc = (c-a)*2.0;
  const coord3d dR2__dc = c*2.0;
  const coord3d dden__dc = dR2__dc * L2/sqrt(L2*R2);
  const coord3d darg__dc = dR2__dc * 1.0/den - dM2__dc * 1.0/den - dden__dc * (L2+R2-M2)/(den*den);

  da = -darg__da * 1.0/sqrt(1.0-arg*arg);
  dc = -darg__dc * 1.0/sqrt(1.0-arg*arg);
}

double coord3d::dihedral(const coord3d& b, const coord3d& c, const coord3d& d) {
  const coord3d ab = b; // a=0
  const coord3d bc = c-b;
  const coord3d cd = d-c;

  const coord3d abc = ab.cross(bc);
  const coord3d bcd = bc.cross(cd);

  const coord3d bc1 = bc/bc.norm();
  const coord3d abc1 = abc/abc.norm();
  const coord3d bcd1 = bcd/bcd.norm();
  const coord3d aux = abc1.cross(bc1);

  const double x = abc1.dot(bcd1);
  const double y = aux.dot(bcd1);

  return atan2(y,x);
}

// calculation of the derivative of dihedral angle theta at a(0,0,0), b, c and d  according to coordinates b, c and d with fixed a
void coord3d::ddihedral(const coord3d& b, const coord3d& c, const coord3d& d, coord3d& db, coord3d& dc, coord3d& dd) {
  const coord3d ab = b; // a=0
  const coord3d bc = c-b;
  const coord3d cd = d-c;

  const double bc_length_inv = 1.0/bc.norm();
  const coord3d bc1 = bc*bc_length_inv;

  const coord3d abc = ab.cross(bc);
  const coord3d bcd = bc.cross(cd);

  const double abc_length_inv = 1.0/abc.norm();
  const double bcd_length_inv = 1.0/bcd.norm();
  const coord3d abc1 = abc * abc_length_inv;
  const coord3d bcd1 = bcd * bcd_length_inv;

  const coord3d aux = abc1.cross(bc1);

  const double x = abc1.dot(bcd1);
  const double y = aux.dot(bcd1);

//  const double dihedral_abcd = atan2(y,x);
//  cout << "D: "<< dihedral_abcd<<endl;

  const matrix3d dab__db = matrix3d::unit_matrix();
  const matrix3d dbc__db = -matrix3d::unit_matrix();
  const matrix3d dbc__dc = matrix3d::unit_matrix();
  const matrix3d dcd__dc = -matrix3d::unit_matrix();
  const matrix3d dcd__dd = matrix3d::unit_matrix();

  // bc_length_inv=1/dsqrt(bc_x**2 + bc_y**2 + bc_z**2)
  const coord3d dbc_length_inv__dbc = -bc * pow(bc_length_inv, 3);

  // bc_length_inv=1/dsqrt(bc_x**2 + bc_y**2 + bc_z**2)
  // vec = vec * mtx
  const coord3d dbc_length_inv__db = dbc_length_inv__dbc * dbc__db;
  const coord3d dbc_length_inv__dc = dbc_length_inv__dbc * dbc__dc;

  const matrix3d dbc1__dbc = matrix3d::unit_matrix() * bc_length_inv;

  // bc1_x=bc_x*bc_length_inv
  // bc1_y=bc_y*bc_length_inv
  // bc1_z=bc_z*bc_length_inv
  // mtx = mtx * mtx + vec outer vec
  const matrix3d dbc1__db = dbc1__dbc*dbc__db + bc.outer(dbc_length_inv__db);
  const matrix3d dbc1__dc = dbc1__dbc*dbc__dc + bc.outer(dbc_length_inv__dc);

  // abc_x=ab_y*bc_z - ab_z*bc_y
  // abc_y=ab_z*bc_x - ab_x*bc_z
  // abc_z=ab_x*bc_y - ab_y*bc_x
  //FIXME is there a more elegant way of doing this?
  const matrix3d dabc__dab = matrix3d(0,bc[2],-bc[1], -bc[2],0,bc[0], bc[1],-bc[0],0);
  const matrix3d dabc__dbc = matrix3d(0,-ab[2],ab[1], ab[2],0,-ab[0], -ab[1],ab[0],0);

  // bcd_x=bc_y*cd_z - bc_z*cd_y
  // bcd_y=bc_z*cd_x - bc_x*cd_z
  // bcd_z=bc_x*cd_y - bc_y*cd_x
  //FIXME is there a more elegant way of doing this?
  const matrix3d dbcd__dbc = matrix3d(0,cd[2],-cd[1], -cd[2],0,cd[0], cd[1],-cd[0],0);
  const matrix3d dbcd__dcd = matrix3d(0,-bc[2],bc[1], bc[2],0,-bc[0], -bc[1],bc[0],0);

  // abc_x=-ab_y*bc_z + ab_z*bc_y
  // abc_y=-ab_z*bc_x + ab_x*bc_z
  // abc_z=-ab_x*bc_y + ab_y*bc_x
  // mtx = mtx * mtx + mtx * mtx
  const matrix3d dabc__db = dabc__dab*dab__db + dabc__dbc*dbc__db;
  const matrix3d dabc__dc =                     dabc__dbc*dbc__dc;

  // bcd_x=-bc_y*cd_z + bc_z*cd_y
  // bcd_y=-bc_z*cd_x + bc_x*cd_z
  // bcd_z=-bc_x*cd_y + bc_y*cd_x
  // mtx = mtx * mtx + mtx * mtx
  const matrix3d dbcd__db = dbcd__dbc*dbc__db;
  const matrix3d dbcd__dc = dbcd__dbc*dbc__dc + dbcd__dcd*dcd__dc;
  const matrix3d dbcd__dd =                     dbcd__dcd*dcd__dd;

  // abc_length_inv=1/dsqrt(abc_x**2 + abc_y**2 + abc_z**2)
  // bcd_length_inv=1/dsqrt(bcd_x**2 + bcd_y**2 + bcd_z**2)
  const coord3d dabc_length_inv__dabc = -abc*pow(abc_length_inv,3);
  const coord3d dbcd_length_inv__dbcd = -bcd*pow(bcd_length_inv,3);

  // abc_length_inv=1/dsqrt(abc_x**2 + abc_y**2 + abc_z**2)
  // vec = vec * mtx
  const coord3d dabc_length_inv__db = dabc_length_inv__dabc*dabc__db;
  const coord3d dabc_length_inv__dc = dabc_length_inv__dabc*dabc__dc;

  // bcd_length_inv=1/dsqrt(bcd_x**2 + bcd_y**2 + bcd_z**2)
  // vec = vec * mtx
  const coord3d dbcd_length_inv__db = dbcd_length_inv__dbcd * dbcd__db;
  const coord3d dbcd_length_inv__dc = dbcd_length_inv__dbcd * dbcd__dc;
  const coord3d dbcd_length_inv__dd = dbcd_length_inv__dbcd * dbcd__dd;

  // abc1_x=abc_x*abc_length_inv
  // abc1_y=abc_y*abc_length_inv
  // abc1_z=abc_z*abc_length_inv
  const matrix3d dabc1__dabc = matrix3d::unit_matrix() * abc_length_inv;

  // abc1_x=abc_x*abc_length_inv
  // abc1_y=abc_y*abc_length_inv
  // abc1_z=abc_z*abc_length_inv
  // mtx = mtx * mtx + vec outer vec
  const matrix3d dabc1__db = dabc1__dabc*dabc__db + abc.outer(dabc_length_inv__db);
  const matrix3d dabc1__dc = dabc1__dabc*dabc__dc + abc.outer(dabc_length_inv__dc);

  // bcd1_x=bcd_x*bcd_length_inv
  // bcd1_y=bcd_y*bcd_length_inv
  // bcd1_z=bcd_z*bcd_length_inv
  const matrix3d dbcd1__dbcd = matrix3d::unit_matrix() * bcd_length_inv;

  // bcd1_x=bcd_x*bcd_length_inv
  // bcd1_y=bcd_y*bcd_length_inv
  // bcd1_z=bcd_z*bcd_length_inv
  // mtx = mtx*mtx + vec outer vec
  const matrix3d dbcd1__db = dbcd1__dbcd * dbcd__db + bcd.outer(dbcd_length_inv__db);
  const matrix3d dbcd1__dc = dbcd1__dbcd * dbcd__dc + bcd.outer(dbcd_length_inv__dc);
  const matrix3d dbcd1__dd = dbcd1__dbcd * dbcd__dd + bcd.outer(dbcd_length_inv__dd);

  // aux_x=abc1_y*bc1_z-bc1_y*abc1_z
  // aux_y=abc1_z*bc1_x-bc1_z*abc1_x
  // aux_z=abc1_x*bc1_y-bc1_x*abc1_y
  //FIXME is there a more elegant way of doing this?
  const matrix3d daux__dabc1 = matrix3d(0,bc1[2],-bc1[1], -bc1[2],0,bc1[0], bc1[1],-bc1[0],0);
  const matrix3d daux__dbc1 = matrix3d(0,-abc1[2],abc1[1], abc1[2],0,-abc1[0], -abc1[1],abc1[0],0);

  // aux_x=abc1_y*bc1_z-bc1_y*abc1_z
  // aux_y=abc1_z*bc1_x-bc1_z*abc1_x
  // aux_z=abc1_x*bc1_y-bc1_x*abc1_y
  // mtx = mtx*mtx + mtx*mtx
  const matrix3d daux__db = daux__dabc1 * dabc1__db + daux__dbc1 * dbc1__db;
  const matrix3d daux__dc = daux__dabc1 * dabc1__dc + daux__dbc1 * dbc1__dc;

  // y=aux_x*bcd1_x + aux_y*bcd1_y + aux_z*bcd1_z
  // vec = vec * mtx
  const coord3d dy__db = bcd1 * daux__db + aux * dbcd1__db;
  const coord3d dy__dc = bcd1 * daux__dc + aux * dbcd1__dc;
  const coord3d dy__dd =                   aux * dbcd1__dd;

  // x=abc1_x*bcd1_x + abc1_y*bcd1_y + abc1_z*bcd1_z
  // vec = vec * mtx
  const coord3d dx__db = bcd1 * dabc1__db + abc1 * dbcd1__db;
  const coord3d dx__dc = bcd1 * dabc1__dc + abc1 * dbcd1__dc;
  const coord3d dx__dd =                    abc1 * dbcd1__dd;

  // df__dx=-y/(x**2 + y**2)
  // df__dy=x/(x**2 + y**2)
  const double df__dx = -y/(x*x + y*y);
  const double df__dy =  x/(x*x + y*y);

  // f=atan2(y,x)
  // vec = vec*sca + vec*sca
  db = dx__db*df__dx + dy__db*df__dy;
  dc = dx__dc*df__dx + dy__dc*df__dy;
  dd = dx__dd*df__dx + dy__dd*df__dy;
}

// coordinates of d-a in the dihedral abcd, where c2=b-a and c3=c-a
// so the call is coord = foo(b-a, c-a, ...) + a
coord3d coord3d::internal2cart(const coord3d& c1, const coord3d& c2, const double r, const double a, const double d){
  //ofstream debug("debug", ofstream::out | ofstream::app);

  //debug << "i2c: " << r << ", " << a << ", " << d << endl;
  //debug << "i2c: " << r << ", " << a*180/M_PI<< ", " << d*180/M_PI << endl;
// c0 is at (0,0,0)
  coord3d c0_work(coord3d(0,0,0)), c1_work(c1), c2_work(c2);

// orient c0, c1, c2: move by -c2
  c0_work -= c2;
  c1_work -= c2;
  c2_work -= c2;
  //debug << "geom (after shift) " << c0_work << ", " << c1_work << ", " << c2_work << endl;

// orient c0, c1, c2: rotate c1 (and c2) such that c1 lies on (-x,0,0)
  double theta_0 = coord3d::angle(c1_work, coord3d(-1,0,0));
  coord3d rot_axis_0 = c1_work.cross(coord3d(-1,0,0));
  //debug << theta_0 << ", " << abs(theta_0 - M_PI) << endl;
  if (abs(theta_0 - M_PI) < 1.e-7 || abs(theta_0) < 1.e-7) {
    rot_axis_0 = coord3d(0,0,1);
    //debug << "overwriting" << endl;
  }
  //debug << "rot axis " << rot_axis_0 << endl;
  rot_axis_0 = (rot_axis_0 / rot_axis_0.norm());
  double u=rot_axis_0[0], v=rot_axis_0[1], w=rot_axis_0[2];
  //debug << "theta_0 " << theta_0 << endl;
  matrix3d rot0 = matrix3d(u*u+(1-u*u)*cos(theta_0),  u*v*(1-cos(theta_0))-w*sin(theta_0),  u*w*(1-cos(theta_0))+v*sin(theta_0),
                           u*v*(1-cos(theta_0))+w*sin(theta_0),  v*v+(1-v*v)*cos(theta_0),  v*w*(1-cos(theta_0))-u*sin(theta_0),
                           u*w*(1-cos(theta_0))-v*sin(theta_0),  v*w*(1-cos(theta_0))+u*sin(theta_0),  w*w+(1-w*w)*cos(theta_0));
  c0_work = rot0 * c0_work;
  c1_work = rot0 * c1_work;
  //debug << "geom (after rot1) " << c0_work << ", " << c1_work << ", " << c2_work << endl;

// orient c0, c1, c2: rotate c0 to lie on (-x,y,0)
  coord3d rot_axis_1 = coord3d(1,0,0);
  //debug << "rot axis " << rot_axis_1 << endl;
  u=rot_axis_1[0], v=rot_axis_1[1], w=rot_axis_1[2];
  double theta_1 = coord3d::angle(coord3d(0,c0_work[1],c0_work[2]), coord3d(0,1,0));
  if(c0_work[2] > 0) theta_1 *= -1;
  //debug << "theta_1 " << theta_1 << endl;
  matrix3d rot1 = matrix3d(u*u+(1-u*u)*cos(theta_1),  u*v*(1-cos(theta_1))-w*sin(theta_1),  u*w*(1-cos(theta_1))+v*sin(theta_1),
                           u*v*(1-cos(theta_1))+w*sin(theta_1),  v*v+(1-v*v)*cos(theta_1),  v*w*(1-cos(theta_1))-u*sin(theta_1),
                           u*w*(1-cos(theta_1))-v*sin(theta_1),  v*w*(1-cos(theta_1))+u*sin(theta_1),  w*w+(1-w*w)*cos(theta_1));
  c0_work = rot1 * c0_work;
  //debug << "geom (after rot2) " << c0_work << ", " << c1_work << ", " << c2_work << endl;

// place c3
  coord3d c3_work(-r,0,0);

// rotate c3 by a around 0/0/z
  coord3d rot_axis_2 = coord3d(0,0,1);
  //debug << "rot axis " << rot_axis_2 << endl;
  u=rot_axis_2[0], v=rot_axis_2[1], w=rot_axis_2[2];
  //debug << "a " << a << endl;
  matrix3d rot2 = matrix3d(u*u+(1-u*u)*cos(-a),  u*v*(1-cos(-a))-w*sin(-a),  u*w*(1-cos(-a))+v*sin(-a),
                           u*v*(1-cos(-a))+w*sin(-a),  v*v+(1-v*v)*cos(-a),  v*w*(1-cos(-a))-u*sin(-a),
                           u*w*(1-cos(-a))-v*sin(-a),  v*w*(1-cos(a))+u*sin(-a),  w*w+(1-w*w)*cos(-a));
  c3_work = rot2 * c3_work;
  //debug << "geom " << c3_work << endl;

// rotate c3 by d around 0/0/z
  coord3d rot_axis_3 = coord3d(1,0,0);
  //debug << "rot axis " << rot_axis_3 << endl;
  u=rot_axis_3[0], v=rot_axis_3[1], w=rot_axis_3[2];
  //debug << "d " << d << endl;
  matrix3d rot3 = matrix3d(u*u+(1-u*u)*cos(-d),  u*v*(1-cos(-d))-w*sin(-d),  u*w*(1-cos(-d))+v*sin(-d),
                           u*v*(1-cos(-d))+w*sin(-d),  v*v+(1-v*v)*cos(-d),  v*w*(1-cos(-d))-u*sin(-d),
                           u*w*(1-cos(-d))-v*sin(-d),  v*w*(1-cos(-d))+u*sin(-d),  w*w+(1-w*w)*cos(-d));
  c3_work = rot3 * c3_work;
  //debug << "geom " << c3_work << endl;

// rotate everything back:
  u=rot_axis_1[0], v=rot_axis_1[1], w=rot_axis_1[2];
  matrix3d rot4 = matrix3d(u*u+(1-u*u)*cos(-theta_1),  u*v*(1-cos(-theta_1))-w*sin(-theta_1),  u*w*(1-cos(-theta_1))+v*sin(-theta_1),
                           u*v*(1-cos(-theta_1))+w*sin(-theta_1),  v*v+(1-v*v)*cos(-theta_1),  v*w*(1-cos(-theta_1))-u*sin(-theta_1),
                           u*w*(1-cos(-theta_1))-v*sin(-theta_1),  v*w*(1-cos(-theta_1))+u*sin(-theta_1),  w*w+(1-w*w)*cos(-theta_1));
  c3_work = rot4 * c3_work;

// rotate everything back:
  u=rot_axis_0[0], v=rot_axis_0[1], w=rot_axis_0[2];
  matrix3d rot5 = matrix3d(u*u+(1-u*u)*cos(-theta_0),  u*v*(1-cos(-theta_0))-w*sin(-theta_0),  u*w*(1-cos(-theta_0))+v*sin(-theta_0),
                           u*v*(1-cos(-theta_0))+w*sin(-theta_0),  v*v+(1-v*v)*cos(-theta_0),  v*w*(1-cos(-theta_0))-u*sin(-theta_0),
                           u*w*(1-cos(-theta_0))-v*sin(-theta_0),  v*w*(1-cos(-theta_0))+u*sin(-theta_0),  w*w+(1-w*w)*cos(-theta_0));
  c3_work = rot5 * c3_work;
  //debug << "geom " << c3_work << endl;

  c3_work += c2;
  //debug << "geom " << c3_work << endl;

  //debug.close();

  return c3_work;
}

void coord3d::scale(const double f){
  for(int i=0; i<3; i++){
    x[i] *= f;
  }
}

// following Klenin-2000
// returns writhe in multiples of 2pi
double structure3d::writhe() const {
  const unsigned int N = this->size();
  double writhe = 0.0;
  for (size_t i=0; i<N; i++){
    if ( atomtypes[i] == to_string('H')) continue;
    for (size_t j=i+2; j<N; j++){
      if (j-i == N-1) continue; // neighbouring segments lead to nan
      // if ( atomtypes[i] == to_string('C') && atomtypes[j] == to_string('C') && j-i == N-3) continue;
// cout << i << ", " << j << endl;
      const coord3d r12 = (*this)[(i+1)%N] - (*this)[i];
      const coord3d r13 = (*this)[j] - (*this)[i];
      const coord3d r14 = (*this)[(j+1)%N] - (*this)[i];
      const coord3d r23 = (*this)[j] - (*this)[(i+1)%N];
      const coord3d r24 = (*this)[(j+1)%N] - (*this)[(i+1)%N];
      const coord3d r34 = (*this)[(j+1)%N] - (*this)[j];
      const coord3d n1=(r13.cross(r14))/(r13.cross(r14)).norm();
      const coord3d n2=(r14.cross(r24))/(r14.cross(r24)).norm();
      const coord3d n3=(r24.cross(r23))/(r24.cross(r23)).norm();
      const coord3d n4=(r23.cross(r13))/(r23.cross(r13)).norm();
// cout << "n1-4 " << n1 << ", "  << n2 << ", " << n3 << ", " << n4 << endl;
      const double sign = ((r34.cross(r12)).dot(r13) > 0) - ((r34.cross(r12)).dot(r13) < 0);
// cout << "s " << sign << endl;
      const double omega = (asin(n1.dot(n2)) + asin(n2.dot(n3)) + asin(n3.dot(n4)) + asin(n4.dot(n1)))*sign;
// cout << "o " << omega << endl;
      writhe += omega;
    }
  }
// cout << "w " << writhe/(2*M_PI) << endl;
  return writhe/(2*M_PI);
}

// assumes that one branchless ring remains
double structure3d::writhe_stripped() const {
  return stripped_backbone().writhe();
}

// removes all atoms but C
structure3d structure3d::stripped_backbone() const {
  structure3d s = *this;
  // remove all but C
  for(int i=s.size()-1; i>=0; i--){
    if(s.atomtypes[i] != "C" && s.atomtypes[i] != "c"){
      s.dat.erase(s.dat.begin()+i);
      s.atomtypes.erase(s.atomtypes.begin()+i);
    }
  }
  return s;
}

structure3d structure3d::get_selection(const vector<int>& selection) const {
  const structure3d& full = *this;
  structure3d s;
  for(int i: selection){
    if(i >= full.size()){
      cerr << "invalid index, bailing out ..." << endl;
      abort();
    }
    s.push_back(full[i]);
    s.push_back_atom(full.atomtypes[i]);
  }

  return s;
}


bool structure3d::is_chain() const {
  const structure3d& s(*this);
  for(size_t i=0; i<size(); i++){
    const double dist = coord3d::dist(s[i], s[(i+1)%s.size()]);
    if(dist<1.2 || dist>1.7){
      cerr << ">>> unexpected distance of " << dist << " between atoms " << i << " and " << (i+1)%s.size() << ". You probably don't want this. <<<" << endl;
      return false;
    }
  }
  return true;
}


double structure3d::homa_stripped() const {
  const structure3d s3d =  stripped_backbone();

  vector<double> lengths;
  for(size_t i=0; i<s3d.size(); i++){
    for(size_t j=i; j<s3d.size(); j++){
      const double dist((s3d.dat[i]-s3d.dat[j]).norm());
      if(dist > 1.0 && dist < 2.0) lengths.push_back(dist);
    }
  }
// cout << "lengths" << lengths << endl;
// cout << "N:" << lengths.size() << endl ;

  double av = 0;
  for(double length: lengths) av += length;
  av /= lengths.size();
// cout << av << endl;

  double sqdiff = 0;
  for(double length: lengths) sqdiff += pow(1 - length/av, 2);
// cout << sqdiff << endl;

  const double alpha = 225;
  const double homa = 1.0-(alpha/lengths.size())*sqdiff;
  return homa;
}


string structure3d::to_turbomole() const {
  const double aa2bohr = 1.889716164632;
  ostringstream s;
  s << setprecision(8) << fixed;
  s << "$coord" << endl;
  for(size_t i=0; i<size(); ++i){
    s << setw(12) << (*this)[i][0]*aa2bohr << "  "<< setw(12) << (*this)[i][1]*aa2bohr << "  " << setw(12) << (*this)[i][2]*aa2bohr << "  " << atomtypes[i] << endl;
  }
  s << "$end" << endl;

  return s.str();
}

string structure3d::to_xyz() const {
  ostringstream s;
  s << setprecision(15) << fixed;
  s << size() << endl;
  s << "i could write something here" << endl;
  for(size_t i=0; i<size(); ++i){
    s << atomtypes[i] << "  " << setw(17) << (*this)[i][0] << "  " << setw(17) << (*this)[i][1] << "  " << setw(17) << (*this)[i][2] << endl;
  }
  return s.str();
}

structure3d structure3d::from_xyz(const string& filename)
{
  ifstream file(filename.c_str());
  unsigned int N;
  string Nstring, comment, element,line;
  structure3d s;

  getline(file, Nstring);
  getline(file, comment);

  N = stoi(Nstring);

  //  cout << "N = " << Nstring << "; comment = " << comment << endl;

  for(size_t i=0; i < N && getline(file,line); i++){
    stringstream l(line);
    coord3d x;

    l >> element;
    for(int j=0;j<3 && l.good(); j++){
      l >> x[j];
    }

    s.push_back(x);
    s.atomtypes.push_back(element);
    //    cout << i << ": " << x << endl;
  }
  file.close();

  assert(s.size() == N);

  return s;
}

void structure3d::scale(const double f){
  for (coord3d &c3d: dat) c3d.scale(f);
}

void structure3d::invert(){
  scale(-1);
}

void structure3d::move(const coord3d m){
  for (coord3d &c3d: dat) c3d+=m;
}

void structure3d::centre_at_origin(){
  coord3d barycentre({0,0,0});
  for (coord3d c3d: dat) barycentre += c3d;
  barycentre /= size();
  move(-barycentre);
}

void structure3d::clear(){
  dat.clear();
  atomtypes.clear();
}

void structure3d::resize(const size_t n){
  dat.resize(n);
  atomtypes.resize(n);
}

void structure3d::erase(const size_t n){
  dat.erase(dat.begin() + n);
  atomtypes.erase(atomtypes.begin() + n);
}

matrix3d coord3d::outer(const coord3d& y) const {
  return matrix3d(x[0]*y[0],x[0]*y[1],x[0]*y[2],  x[1]*y[0],x[1]*y[1],x[1]*y[2],  x[2]*y[0],x[2]*y[1],x[2]*y[2]);
}

void structure3d::push_back(const structure3d& s3d){
  for(coord3d c3d: s3d.dat) push_back(c3d);
  for(string a: s3d.atomtypes) push_back_atom(a);
}

matrix3d matrix3d::inverse() const {

  const double den = (-(values[2]*values[4]*values[6]) + values[1]*values[5]*values[6] + values[2]*values[3]*values[7] - values[0]*values[5]*values[7] - values[1]*values[3]*values[8] + values[0]*values[4]*values[8]);
  const double m1 = (-(values[5]*values[7]) + values[4]*values[8])/den;
  const double m2 = (values[2]*values[7] - values[1]*values[8])/den;
  const double m3 = (-(values[2]*values[4]) + values[1]*values[5])/den;
  const double m4 = (values[5]*values[6] - values[3]*values[8])/den;
  const double m5 = (-(values[2]*values[6]) + values[0]*values[8])/den;
  const double m6 = (values[2]*values[3] - values[0]*values[5])/den;
  const double m7 = (-(values[4]*values[6]) + values[3]*values[7])/den;
  const double m8 = (values[1]*values[6] - values[0]*values[7])/den;
  const double m9 = (-(values[1]*values[3]) + values[0]*values[4])/den;

  return matrix3d(m1, m2, m3, m4, m5, m6, m7, m8, m9);
}


