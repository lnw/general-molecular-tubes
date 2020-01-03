# Copyright (c) 2019, Lukas Wirz
# All rights reserved.

# This file is part of 'general-molecular-tubes' which is released under the
# BSD-2-clause license.  See file LICENSE in this project.

cmake_minimum_required(VERSION 3.2.0)
project(general-tube-gen)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

set(libs)

include_directories(AFTER SYSTEM .)

set(PROGRAMS1
  general-tubes
)

add_library(geom geometry3.cc tile.cc auxiliary.cc)

foreach(program ${PROGRAMS1})
  set(program_source ${program}.cc)
  add_executable(${program} ${program_source})
  target_link_libraries(${program} geom)
endforeach(program ${PROGRAMS})

