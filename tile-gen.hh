#ifndef TILES_HH
#define TILES_HH

#include <fstream>
#include <string>
#include <vector>

#include "auxiliary.hh"
#include "geometry3.hh"
#include "tile.hh"


using namespace std;


class TileGenerator {

  const double rCC_s = 1.54,
               rCC_sesqui = 1.39,
               rCC_t = 1.20,
               rCH = 1.09;


// annulene cycles with CC triple
Tile tile_1(){
  structure3d str;
  str.push_back(coord3d(0,                      0                                  , 0 )), str.push_back_atom("C");
  str.push_back(coord3d(rCC_sesqui*cos(M_PI/6), rCC_sesqui*sin(M_PI/6)             , 0 )), str.push_back_atom("C");
  str.push_back(coord3d(rCC_sesqui*cos(M_PI/6), rCC_sesqui*sin(M_PI/6)+rCC_s       , 0 )), str.push_back_atom("C");
  str.push_back(coord3d(rCC_sesqui*cos(M_PI/6), rCC_sesqui*sin(M_PI/6)+rCC_s+rCC_t , 0 )), str.push_back_atom("C");
  str.push_back(coord3d(0,                      -rCH                               , 0 )), str.push_back_atom("H");
  str.push_back(coord3d(rCC_sesqui*cos(M_PI/6), rCC_sesqui*sin(M_PI/6)+rCH         , 0 )), str.push_back_atom("H");
  Tile T(rCC_sesqui*2*cos(M_PI/6), // width
         rCC_sesqui*sin(M_PI/6) + 2*rCC_s + rCC_t, // height
         0.5, // shift
         str,
         vector<size_t>({0, 1, 2, 3}),
         vector<size_t>({0, 1, 2, 3, 4}),
         vector<size_t>({0, 1,          5})
         );
  return T;
};


// annulene cycles with gold and nitrogen
Tile tile_2(){
  structure3d str;
  str.push_back(coord3d(0,                      0                                      , 0   )), str.push_back_atom("C");
  str.push_back(coord3d(rCC_sesqui*cos(M_PI/6), rCC_sesqui*sin(M_PI/6)                 , 0   )), str.push_back_atom("N");
  str.push_back(coord3d(rCC_sesqui*cos(M_PI/6), rCC_sesqui*sin(M_PI/6)+rCC_s+0.5*rCC_t , 0.2 )), str.push_back_atom("Au");
  str.push_back(coord3d(0,                      -rCH                                   , 0   )), str.push_back_atom("H");
  Tile T(rCC_sesqui*2*cos(M_PI/6),
         rCC_sesqui*sin(M_PI/6) + 2*rCC_s + rCC_t,
         0.5, // shift
         str,
         vector<size_t>({0, 1, 2}),
         vector<size_t>({0, 1, 2, 3}),
         vector<size_t>({0, 1})
         );
  return T;
};


// annulene cycles with CC triple
Tile tile_3(){
  structure3d str;
  str.push_back(coord3d(0,                      rCC_sesqui*sin(M_PI/6)                   , 0 )), str.push_back_atom("C");
  str.push_back(coord3d(0,                      rCC_sesqui*(sin(M_PI/6)+1)               , 0 )), str.push_back_atom("C");
  str.push_back(coord3d(rCC_sesqui*cos(M_PI/6), 0                                        , 0 )), str.push_back_atom("C");
  str.push_back(coord3d(rCC_sesqui*cos(M_PI/6), rCC_sesqui*(2*sin(M_PI/6)+1)             , 0 )), str.push_back_atom("C");
  str.push_back(coord3d(rCC_sesqui*cos(M_PI/6), rCC_sesqui*(2*sin(M_PI/6)+1)+rCC_s       , 0 )), str.push_back_atom("C");
  str.push_back(coord3d(rCC_sesqui*cos(M_PI/6), rCC_sesqui*(2*sin(M_PI/6)+1)+rCC_s+rCC_t , 0 )), str.push_back_atom("C");
  str.push_back(coord3d(rCC_sesqui*cos(M_PI/6), -rCH                                     , 0 )), str.push_back_atom("H");
  str.push_back(coord3d(rCC_sesqui*cos(M_PI/6), rCC_sesqui*(2*sin(M_PI/6)+1)+rCH         , 0 )), str.push_back_atom("H");
  Tile T(rCC_sesqui*2*cos(M_PI/6),
         rCC_sesqui*(2*sin(M_PI/6)+1) + 2*rCC_s + rCC_t,
         0.0, // shift
         str,
         vector<size_t>({0, 1, 2, 3, 4, 5}),
         vector<size_t>({0, 1, 2, 3, 4, 5, 6}),
         vector<size_t>({0, 1, 2, 3,          7})
         );
  return T;
};


// annulene cycles with CC triple
Tile tile_4(){
  structure3d str;
  str.push_back(coord3d(0,                      rCC_sesqui*sin(M_PI/6)                       , 0 )), str.push_back_atom("C");
  str.push_back(coord3d(0,                      rCC_sesqui*(sin(M_PI/6)+1)                   , 0 )), str.push_back_atom("C");
  str.push_back(coord3d(rCC_sesqui*cos(M_PI/6), 0                                            , 0 )), str.push_back_atom("C");
  str.push_back(coord3d(rCC_sesqui*cos(M_PI/6), rCC_sesqui*(2*sin(M_PI/6)+1)                 , 0 )), str.push_back_atom("N");
  str.push_back(coord3d(rCC_sesqui*cos(M_PI/6), rCC_sesqui*(2*sin(M_PI/6)+1)+rCC_s+0.5*rCC_t , 0 )), str.push_back_atom("Au");
  str.push_back(coord3d(rCC_sesqui*cos(M_PI/6), -rCH                                         , 0 )), str.push_back_atom("H");
  Tile T(rCC_sesqui*2*cos(M_PI/6),
         rCC_sesqui*(2*sin(M_PI/6)+1) + 2*rCC_s + rCC_t,
         0.0, // shift
         str,
         vector<size_t>({0, 1, 2, 3, 4}),
         vector<size_t>({0, 1, 2, 3, 4, 5}),
         vector<size_t>({0, 1, 2, 3})
         );
  return T;
};


//
Tile tile_5(){
  const double spacer_height = sqrt( pow(2*rCC_s+rCC_t, 2) - pow(rCC_sesqui*sin(M_PI/6), 2));
  const double ratio_1 = (rCC_s)/(2*rCC_s+rCC_t);
  const double ratio_2 = (rCC_s+rCC_t)/(2*rCC_s+rCC_t);
  structure3d str;
  str.push_back(coord3d(0,                                                           0                                              , 0)), str.push_back_atom("C");
  str.push_back(coord3d(rCC_sesqui*sin(M_PI/6),                                      rCC_sesqui*cos(M_PI/6)                         , 0)), str.push_back_atom("C");
  str.push_back(coord3d(0,                                                           rCC_sesqui*2*cos(M_PI/6)                       , 0)), str.push_back_atom("C");
  str.push_back(coord3d(rCC_sesqui*sin(M_PI/6),                                      rCC_sesqui*3*cos(M_PI/6)                       , 0)), str.push_back_atom("C");
  str.push_back(coord3d(rCC_sesqui*(2*sin(M_PI/6)+1),                                0                                              , 0)), str.push_back_atom("C");
  str.push_back(coord3d(rCC_sesqui*(sin(M_PI/6)+1),                                  rCC_sesqui*cos(M_PI/6)                         , 0)), str.push_back_atom("C");
  str.push_back(coord3d(rCC_sesqui*(2*sin(M_PI/6)+1),                                rCC_sesqui*2*cos(M_PI/6)                       , 0)), str.push_back_atom("C");
  str.push_back(coord3d(rCC_sesqui*(sin(M_PI/6)+1),                                  rCC_sesqui*3*cos(M_PI/6)                       , 0)), str.push_back_atom("C");
  str.push_back(coord3d(ratio_2*rCC_sesqui*sin(M_PI/6),                              rCC_sesqui*3*cos(M_PI/6)+ratio_1*spacer_height , 0)), str.push_back_atom("C");
  str.push_back(coord3d(rCC_sesqui*(sin(M_PI/6)+1) + ratio_1*rCC_sesqui*sin(M_PI/6), rCC_sesqui*3*cos(M_PI/6)+ratio_1*spacer_height , 0)), str.push_back_atom("C");
  str.push_back(coord3d(ratio_1*rCC_sesqui*sin(M_PI/6),                              rCC_sesqui*3*cos(M_PI/6)+ratio_2*spacer_height , 0)), str.push_back_atom("C");
  str.push_back(coord3d(rCC_sesqui*(sin(M_PI/6)+1) + ratio_2*rCC_sesqui*sin(M_PI/6), rCC_sesqui*3*cos(M_PI/6)+ratio_2*spacer_height , 0)), str.push_back_atom("C");
  str.push_back(coord3d(rCC_sesqui*sin(M_PI/6),                                      -rCH*cos(M_PI/6)                               , 0)), str.push_back_atom("H");
  str.push_back(coord3d(rCC_sesqui*(sin(M_PI/6)+1),                                  -rCH*cos(M_PI/6)                               , 0)), str.push_back_atom("H");
  str.push_back(coord3d(0,                                                           rCC_sesqui*3*cos(M_PI/6)+rCH*cos(M_PI/6)       , 0)), str.push_back_atom("H");
  str.push_back(coord3d(rCC_sesqui*(2*sin(M_PI/6)+1),                                rCC_sesqui*3*cos(M_PI/6)+rCH*cos(M_PI/6)       , 0)), str.push_back_atom("H");
  Tile T(rCC_sesqui*(2*sin(M_PI/6)+2),
         spacer_height + rCC_sesqui*3*cos(M_PI/6),
         0.0, // shift
         str,
         vector<size_t>({0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11}), // default
         vector<size_t>({0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13}), // first
         vector<size_t>({0, 1, 2, 3, 4, 5, 6, 7,                       14, 15}) // last
         );
  return T;
};


//
Tile tile_6(){
  const double spacer_height = sqrt( pow(2*rCC_s+rCC_t, 2) - pow(rCC_sesqui*sin(M_PI/6), 2));
  structure3d str;
  str.push_back(coord3d(0,                                                       0                                          , 0)), str.push_back_atom("N");
  str.push_back(coord3d(rCC_sesqui*sin(M_PI/6),                                  rCC_sesqui*cos(M_PI/6)                     , 0)), str.push_back_atom("C");
  str.push_back(coord3d(0,                                                       rCC_sesqui*2*cos(M_PI/6)                   , 0)), str.push_back_atom("C");
  str.push_back(coord3d(rCC_sesqui*sin(M_PI/6),                                  rCC_sesqui*3*cos(M_PI/6)                   , 0)), str.push_back_atom("C");
  str.push_back(coord3d(rCC_sesqui*(2*sin(M_PI/6)+1),                            0                                          , 0)), str.push_back_atom("C");
  str.push_back(coord3d(rCC_sesqui*(sin(M_PI/6)+1),                              rCC_sesqui*cos(M_PI/6)                     , 0)), str.push_back_atom("C");
  str.push_back(coord3d(rCC_sesqui*(2*sin(M_PI/6)+1),                            rCC_sesqui*2*cos(M_PI/6)                   , 0)), str.push_back_atom("C");
  str.push_back(coord3d(rCC_sesqui*(sin(M_PI/6)+1),                              rCC_sesqui*3*cos(M_PI/6)                   , 0)), str.push_back_atom("N");
  str.push_back(coord3d(0.5*rCC_sesqui*sin(M_PI/6),                              rCC_sesqui*3*cos(M_PI/6)+0.5*spacer_height , 0)), str.push_back_atom("Au");
  str.push_back(coord3d(rCC_sesqui*(sin(M_PI/6)+1) + 0.5*rCC_sesqui*sin(M_PI/6), rCC_sesqui*3*cos(M_PI/6)+0.5*spacer_height , 0)), str.push_back_atom("Au");
  str.push_back(coord3d(rCC_sesqui*(sin(M_PI/6)+1),                              -rCH*cos(M_PI/6)                           , 0)), str.push_back_atom("H");
  str.push_back(coord3d(0,                                                       rCC_sesqui*3*cos(M_PI/6)+rCH*cos(M_PI/6)   , 0)), str.push_back_atom("H");
  Tile T(rCC_sesqui*(2*sin(M_PI/6)+2),
         spacer_height + rCC_sesqui*3*cos(M_PI/6),
         0.0, // shift
         str,
         vector<size_t>({0, 1, 2, 3, 4, 5, 6, 7, 8, 9}), // default
         vector<size_t>({0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10}), // first
         vector<size_t>({0, 1, 2, 3, 4, 5, 6, 7,           11}) // last
         );
  return T;
};


//
Tile tile_7(){
  const double spacer_height = sqrt( pow(2*rCC_s+rCC_t, 2) - pow(rCC_sesqui*sin(M_PI/6), 2));
  structure3d str;
  str.push_back(coord3d(0,                                                       0                                         , 0 )), str.push_back_atom("C");
  str.push_back(coord3d(rCC_sesqui*sin(M_PI/6),                                  rCC_sesqui*cos(M_PI/6)                    , 0 )), str.push_back_atom("C");
  str.push_back(coord3d(0,                                                       rCC_sesqui*2*cos(M_PI/6)                  , 0 )), str.push_back_atom("C");
  str.push_back(coord3d(rCC_sesqui*sin(M_PI/6),                                  rCC_sesqui*3*cos(M_PI/6)                  , 0 )), str.push_back_atom("N");
  str.push_back(coord3d(rCC_sesqui*(2*sin(M_PI/6)+1),                            0                                         , 0 )), str.push_back_atom("C");
  str.push_back(coord3d(rCC_sesqui*(sin(M_PI/6)+1),                              rCC_sesqui*cos(M_PI/6)                    , 0 )), str.push_back_atom("C");
  str.push_back(coord3d(rCC_sesqui*(2*sin(M_PI/6)+1),                            rCC_sesqui*2*cos(M_PI/6)                  , 0 )), str.push_back_atom("C");
  str.push_back(coord3d(rCC_sesqui*(sin(M_PI/6)+1),                              rCC_sesqui*3*cos(M_PI/6)                  , 0 )), str.push_back_atom("N");
  str.push_back(coord3d(0.5*rCC_sesqui*sin(M_PI/6),                              rCC_sesqui*3*cos(M_PI/6)+0.5*spacer_height, 0 )), str.push_back_atom("Au");
  str.push_back(coord3d(rCC_sesqui*(sin(M_PI/6)+1) + 0.5*rCC_sesqui*sin(M_PI/6), rCC_sesqui*3*cos(M_PI/6)+0.5*spacer_height, 0 )), str.push_back_atom("Au");
  str.push_back(coord3d(rCC_sesqui*sin(M_PI/6),                                  -rCH*cos(M_PI/6)                          , 0 )), str.push_back_atom("H");
  str.push_back(coord3d(rCC_sesqui*(sin(M_PI/6)+1),                              -rCH*cos(M_PI/6)                          , 0 )), str.push_back_atom("H");
  Tile T(rCC_sesqui*(2*sin(M_PI/6)+2),
         spacer_height + rCC_sesqui*3*cos(M_PI/6),
         0.0, // shift
         str,
         vector<size_t>({0, 1, 2, 3, 4, 5, 6, 7, 8, 9}), // default
         vector<size_t>({0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11}), // first
         vector<size_t>({0, 1, 2, 3, 4, 5, 6, 7}) // last
         );
  return T;
};


//
Tile tile_8(){
  structure3d str;
  str.push_back(coord3d(0,                                                       rCC_sesqui*sin(M_PI/6)                              , 0 )), str.push_back_atom("C");
  str.push_back(coord3d(0,                                                       rCC_sesqui*(sin(M_PI/6)+1)                          , 0 )), str.push_back_atom("C");
  str.push_back(coord3d(rCC_sesqui*cos(M_PI/6),                                  0                                                   , 0 )), str.push_back_atom("C");
  str.push_back(coord3d(rCC_sesqui*cos(M_PI/6),                                  rCC_sesqui*(2*sin(M_PI/6)+1)                        , 0 )), str.push_back_atom("C");
  str.push_back(coord3d(rCC_sesqui*2*cos(M_PI/6),                                rCC_sesqui*sin(M_PI/6)                              , 0 )), str.push_back_atom("C");
  str.push_back(coord3d(rCC_sesqui*2*cos(M_PI/6),                                rCC_sesqui*(sin(M_PI/6)+1)                          , 0 )), str.push_back_atom("C");
  str.push_back(coord3d(rCC_sesqui*3*cos(M_PI/6)+(2*rCC_s+rCC_t)*cos(M_PI/6),    rCC_sesqui*(sin(M_PI/6)+0.5)-rCC_t/2                , 0 )), str.push_back_atom("C"); // vert
  str.push_back(coord3d(rCC_sesqui*3*cos(M_PI/6)+(2*rCC_s+rCC_t)*cos(M_PI/6),    rCC_sesqui*(sin(M_PI/6)+0.5)+rCC_t/2                , 0 )), str.push_back_atom("C"); // vert
  str.push_back(coord3d(rCC_sesqui*2*cos(M_PI/6)+rCC_s*cos(M_PI/6),              rCC_sesqui*(3*sin(M_PI/6))+rCC_s*sin(M_PI/6)        , 0 )), str.push_back_atom("C"); // up
  str.push_back(coord3d(rCC_sesqui*2*cos(M_PI/6)+(rCC_s+rCC_t)*cos(M_PI/6),      rCC_sesqui*(3*sin(M_PI/6))+(rCC_s+rCC_t)*sin(M_PI/6), 0 )), str.push_back_atom("C"); // up
  str.push_back(coord3d(rCC_sesqui*4*cos(M_PI/6)+(3*rCC_s+2*rCC_t)*cos(M_PI/6),  rCC_sesqui*(3*sin(M_PI/6))+rCC_s*sin(M_PI/6)        , 0 )), str.push_back_atom("C"); // down
  str.push_back(coord3d(rCC_sesqui*4*cos(M_PI/6)+(3*rCC_s+rCC_t)*cos(M_PI/6),    rCC_sesqui*(3*sin(M_PI/6))+(rCC_s+rCC_t)*sin(M_PI/6), 0 )), str.push_back_atom("C"); // down
  str.push_back(coord3d(rCC_sesqui*cos(M_PI/6),                                  -rCH                                                , 0 )), str.push_back_atom("H");
  str.push_back(coord3d(rCC_sesqui*3*cos(M_PI/6)+(2*rCC_s+rCC_t)*cos(M_PI/6),    (2*rCC_sesqui+2*rCC_s+rCC_t)/2-rCH                  , 0 )), str.push_back_atom("H");
  str.push_back(coord3d(rCC_sesqui*2*cos(M_PI/6)+cos(M_PI/6)*rCH,                (rCC_sesqui-rCH)*sin(M_PI/6)                        , 0 )), str.push_back_atom("H");
  str.push_back(coord3d(rCC_sesqui*4*cos(M_PI/6)+(4*rCC_s+2*rCC_t)*cos(M_PI/6)-cos(M_PI/6)*rCH,  (rCC_sesqui-rCH)*sin(M_PI/6)        , 0 )), str.push_back_atom("H");
  str.push_back(coord3d(rCC_sesqui*cos(M_PI/6),                                  2*rCC_sesqui+rCH                                    , 0 )), str.push_back_atom("H");
  str.push_back(coord3d((rCC_sesqui*2+rCH)*cos(M_PI/6),                          rCC_sesqui*(sin(M_PI/6)+1)+rCH*sin(M_PI/6)          , 0 )), str.push_back_atom("H");
  str.push_back(coord3d(rCC_sesqui*4*cos(M_PI/6)+(4*rCC_s+2*rCC_t-rCH)*cos(M_PI/6),      rCC_sesqui*(sin(M_PI/6)+1)+rCH*sin(M_PI/6)  , 0 )), str.push_back_atom("H");
  str.push_back(coord3d(rCC_sesqui*3*cos(M_PI/6)+(2*rCC_s+rCC_t)*cos(M_PI/6),    0.34                                                , 0 )), str.push_back_atom("H"); // vert
  Tile T(rCC_sesqui*4*cos(M_PI/6)+(4*rCC_s+2*rCC_t)*cos(M_PI/6), // x
         (2*rCC_sesqui+2*rCC_s+rCC_t)/2, // y
         0.5, // shift
         str,
         vector<size_t>({0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11}), // default
         vector<size_t>({0, 1, 2, 3, 4, 5,       8, 9, 10, 11, 12, 13, 14, 15}), // first
         vector<size_t>({0, 1, 2, 3, 4, 5,                                     16, 17, 18, 19}) // last
         );
  return T;
};


//
Tile tile_9(){
  structure3d str;
  str.push_back(coord3d(0,                                                          rCC_sesqui*sin(M_PI/6)                                 , 0)), str.push_back_atom("C");
  str.push_back(coord3d(0,                                                          rCC_sesqui*(sin(M_PI/6)+1)                             , 0)), str.push_back_atom("N");
  str.push_back(coord3d(rCC_sesqui*cos(M_PI/6),                                     0                                                      , 0)), str.push_back_atom("N");
  str.push_back(coord3d(rCC_sesqui*cos(M_PI/6),                                     rCC_sesqui*(2*sin(M_PI/6)+1)                           , 0)), str.push_back_atom("C");
  str.push_back(coord3d(rCC_sesqui*2*cos(M_PI/6),                                   rCC_sesqui*sin(M_PI/6)                                 , 0)), str.push_back_atom("C");
  str.push_back(coord3d(rCC_sesqui*2*cos(M_PI/6),                                   rCC_sesqui*(sin(M_PI/6)+1)                             , 0)), str.push_back_atom("N");
  str.push_back(coord3d(rCC_sesqui*3*cos(M_PI/6)+(2*rCC_s+rCC_t)*cos(M_PI/6),       rCC_sesqui*(sin(M_PI/6)+0.5)                           , 0)), str.push_back_atom("Au"); // vert
  str.push_back(coord3d(rCC_sesqui*2*cos(M_PI/6)+(rCC_s+rCC_t/2)*cos(M_PI/6),       rCC_sesqui*(3*sin(M_PI/6))+(rCC_s+rCC_t/2)*sin(M_PI/6) , 0)), str.push_back_atom("Au"); // up
  str.push_back(coord3d(rCC_sesqui*4*cos(M_PI/6)+(3*rCC_s+1.5*rCC_t)*cos(M_PI/6),   rCC_sesqui*(3*sin(M_PI/6))+(rCC_s+rCC_t/2)*sin(M_PI/6) , 0)), str.push_back_atom("Au"); // down
  str.push_back(coord3d(rCC_sesqui*2*cos(M_PI/6)+cos(M_PI/6)*rCH,                   (rCC_sesqui-rCH)*sin(M_PI/6)                           , 0)), str.push_back_atom("H");
  str.push_back(coord3d(rCC_sesqui*4*cos(M_PI/6)+(4*rCC_s+2*rCC_t-rCH)*cos(M_PI/6), (rCC_sesqui-rCH)*sin(M_PI/6)                           , 0)), str.push_back_atom("H");
  str.push_back(coord3d(rCC_sesqui*cos(M_PI/6),                                     2*rCC_sesqui+rCH                                       , 0)), str.push_back_atom("H");
  str.push_back(coord3d(rCC_sesqui*3*cos(M_PI/6)+(2*rCC_s+rCC_t)*cos(M_PI/6),       0.34                                                   , 0)), str.push_back_atom("H");
  Tile T(rCC_sesqui*4*cos(M_PI/6)+(4*rCC_s+2*rCC_t)*cos(M_PI/6), // x
         (2*rCC_sesqui+2*rCC_s+rCC_t)/2, // y
         0.5, // shift
         str,
         vector<size_t>({0, 1, 2, 3, 4, 5, 6, 7, 8}), // default
         vector<size_t>({0, 1, 2, 3, 4, 5,    7, 8, 9, 10}), // first
         vector<size_t>({0, 1, 2, 3, 4, 5,                 11, 12}) // last
         );
  return T;
};



  public:
    Tile get_tile(size_t c){
      switch(c){
      case 1: return tile_1();
      case 2: return tile_2();
      case 3: return tile_3();
      case 4: return tile_4();
      case 5: return tile_5();
      case 6: return tile_6();
      case 7: return tile_7();
      case 8: return tile_8();
      case 9: return tile_9();
      default:
        cerr << "invalid curve chosen, aborting ..." << endl;
        __builtin_unreachable();
        break;
      }
    }

};

#endif


