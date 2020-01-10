# Generation of general molecular tubes

A program to generate structures of molecular tubes (open ended or closed),
consisting of (currently only rectangular) 3D tiles.  You can easily define
tiles of your own choice by adding an entry to `tile-gen.hh`.

## Licensing / Citing

This software is published under a very permissive license to make it as useful
as possible for a variety of usecases.  If it is useful to you, especially in
an academic context, you are encouraged to cite: [to be added later]


Some parts of `geometry3.{hh,cc}` and `auxiliary.hh` have been borrowed from
http://ctcp.massey.ac.nz/index.php?group=&page=fullerenes .


## Compiling

There are two ways to build this little program.  Should work on any Linux and
related systems as well as macOS.

### using the ordinary Makefile

```
make
```

### using CMake (probably more robust)

```
mkdir mybuild
cd mybuild
[CXX=clang++] cmake ..
make
```
## Usage

There are four input parameters:

### tile type

Currently there are 9 tile types defined (indices 1 through 9).  Edit tile-gen.hh to define
your own tiles.

### shape type

Currently there are 5 shapes defined (indices 1 through 5).  1: linear tube, 2:
torus, 3-5:  much less useful, a trefoil, a lemniscate, etc.

### number of tiles per layer

Too few tiles per layer lead to high distortions.

### number of layers

If choosing a shape other than 1, too few layers lead to high distortions.

