# Generation of general molecular tubes

a program to generate structures of molecular tubes (open ended or closed),
consisting of (currently only rectangular) 3D tiles.  You can easily define
tiles of your own choice by adding an entry to `tile-gen.hh`.

## Licensing / Citing

This software is published under a very permissive license to make it as useful
as possible for a variety of usecases.  If it is useful to you, especially in a
scientific context, you are encouraged to cite: ...


Some parts of ```geometry{2,3}.{hh,cc}``` and ```auxiliary.hh``` have been borrowed from
http://ctcp.massey.ac.nz/index.php?group=&page=fullerenes .


## Compiling

There are two ways to build these little programs.  Should work on any Linux
and related systems as well as macOS.

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

