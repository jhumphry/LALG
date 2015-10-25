# LALG

## Introduction

This is an Ada 2012 package that provides an interface to dense linear
algebra packages. Initially the focus is on
[BLAS](http://www.netlib.org/blas), which has been covered except for
routines working on banded matrices.

Other Ada interfaces to BLAS exist, but they tend to have limitations.
Any Ada compiler that supports `Interfaces.Fortran` can compile an
interface that calls the routines in BLAS, passing in Ada arrays with
`Convention => Fortran` set or, with a little more work, mutating the
operand order and the various `TRANS` parameters to allow Ada arrays
with the default `Convention` to be used. Such an interface can do
basic linear algebra operations, but it cannot operate on sub-slices or
sub-blocks of vectors or matrices without making temporary copies.

The aim of this project is to develop a set of vector and matrix types
in a tagged type hierarchy that encompasses both types that wrap
'concrete' in-memory arrays and types that provide a 'view' onto a
slice or sub-block of one of the 'concrete' types. If `A` is an object
of type `Concrete_Real_Matrix`, for example, then `A.Row(2)` will
return a `Real_Matrix_Vector` object which can be used by any routine
which requires a vector and will operate in-place on the second row of
the matrix `A`. This requires some assumptions about ABI compatibility
between Ada and Fortran on the platform in question, but in practice few
problems are expected.

The Ada subprograms in `LALG.Real_BLAS` have the same names as in BLAS,
but without the precision specifier `s` or `d` (which are unnecessary as
`LALG` will have been instantiated for one precision or the other).
Parameters relating to the size or stride of the inputs are not
required, as the vector and matrix types will provide the correct
values each time.

A set of unit tests has been provided. These do not test the BLAS
library for correctness or accuracy. Instead they run limited checks
against each routine in the Ada interface to make sure that the
interfacing with the Fortran code is functional.

Two `.gpr` project files have been provided for those using
GNAT/GPRBuild. There are two parameters - a `mode` parameters that can
be set to `debug` or `optimize` and a `blas` parameter that can be set
to `system`, `netlib` or `openblas` to select the BLAS library that you
wish to link to. It may be necessary to change the 'Linker' settings in
`lalg_tests.gpr` if your desired BLAS library is not in the expected
place.
