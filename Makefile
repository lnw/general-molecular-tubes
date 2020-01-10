# Copyright (c) 2020, Lukas Wirz
# All rights reserved.

# This file is part of 'general-molecular-tubes' which is released under the
# BSD-2-clause license.  See file LICENSE in this project.


#CXX=clang++
CXX=g++

# FLAGS=-O1 -g -std=c++17 -Wunused -Wshadow -Wall
FLAGS=-O3 -std=c++17 -g -DNDEBUG

HEADERS=geometry3.hh auxiliary.hh shape-gen.hh tile-gen.hh
OBJECTS= geometry3.o tile.o
OBJECTS_P=$(patsubst %.o, build/%.o, $(OBJECTS))

build/%.o: %.cc $(HEADERS) Makefile
	$(CXX) $(FLAGS) -c $< -o $@

all: general-tubes

general-tubes: Makefile general-tubes.cc $(OBJECTS_P) $(HEADERS)
	$(CXX) $(FLAGS) general-tubes.cc $(OBJECTS_P) -o $@

clean:
	rm -f output/*.xyz output/*.coord test debug final
	rm -f tube*.xyz tube-*coord

distclean: clean
	rm -f build/*o
	rm -f general-tubes

