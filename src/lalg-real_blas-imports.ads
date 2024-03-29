-- LALG
-- An Ada 2012 binding to BLAS and other linear algebra routines

-- Copyright (c) 2015-2021, James Humphry - see LICENSE file for details
-- SPDX-License-Identifier: ISC

private generic
package LALG.Real_BLAS.Imports is

   subtype FI is IntFort.Fortran_Integer;
   subtype FP is FI range 1 .. FI'Last;
   subtype FC is IntFort.Character_Set;

   BLAS_Prefix : constant String := "";
   BLAS_Suffix : constant String := "_";

   -- *************
   -- *************
   -- ** Level 1 **
   -- *************
   -- *************

   procedure SROTG (SA, SB : in out Real; C, S : out Real) with
      Import        => True,
      Convention    => Fortran,
      External_Name => BLAS_Prefix & "srotg" & BLAS_Suffix;

   procedure DROTG (SA, SB : in out Real; C, S : out Real) with
      Import        => True,
      Convention    => Fortran,
      External_Name => BLAS_Prefix & "drotg" & BLAS_Suffix;

   procedure SROTMG
     (SD1, SD2 : in out Real;
      SX1      : in out Real;
      SY1      : in     Real;
      PARAMS   : in out Modified_Givens_Params) with
      Import        => True,
      Convention    => Fortran,
      External_Name => BLAS_Prefix & "srotmg" & BLAS_Suffix;

   procedure DROTMG
     (DD1, DD2 : in out Real;
      DX1      : in out Real;
      DY1      : in     Real;
      PARAMS   : in out Modified_Givens_Params) with
      Import        => True,
      Convention    => Fortran,
      External_Name => BLAS_Prefix & "drotmg" & BLAS_Suffix;

   procedure SROT
     (N    : FI;
      SX   : Real_Vector_Handle;
      INCX : FP;
      SY   : Real_Vector_Handle;
      INCY : FP;
      C    : Real;
      S    : Real) with
      Import        => True,
      Convention    => Fortran,
      External_Name => BLAS_Prefix & "srot" & BLAS_Suffix;

   procedure DROT
     (N    : FI;
      DX   : Real_Vector_Handle;
      INCX : FP;
      DY   : Real_Vector_Handle;
      INCY : FP;
      C    : Real;
      S    : Real) with
      Import        => True,
      Convention    => Fortran,
      External_Name => BLAS_Prefix & "drot" & BLAS_Suffix;

   procedure SROTM
     (N      :        FP;
      SX     :        Real_Vector_Handle;
      INCX   :        FP;
      SY     :        Real_Vector_Handle;
      INCY   :        FP;
      PARAMS : in out Modified_Givens_Params) with
      Import        => True,
      Convention    => Fortran,
      External_Name => BLAS_Prefix & "srotm" & BLAS_Suffix;

   procedure DROTM
     (N      :        FP;
      DX     :        Real_Vector_Handle;
      INCX   :        FP;
      DY     :        Real_Vector_Handle;
      INCY   :        FP;
      PARAMS : in out Modified_Givens_Params) with
      Import        => True,
      Convention    => Fortran,
      External_Name => BLAS_Prefix & "drotm" & BLAS_Suffix;

   procedure SSWAP
     (N    : FI;
      SX   : Real_Vector_Handle;
      INCX : FP;
      SY   : Real_Vector_Handle;
      INCY : FP) with
      Import        => True,
      Convention    => Fortran,
      External_Name => BLAS_Prefix & "sswap" & BLAS_Suffix;

   procedure DSWAP
     (N    : FI;
      DX   : Real_Vector_Handle;
      INCX : FP;
      DY   : Real_Vector_Handle;
      INCY : FP) with
      Import        => True,
      Convention    => Fortran,
      External_Name => BLAS_Prefix & "dswap" & BLAS_Suffix;

   procedure SSCAL (N : FP; SA : Real; SX : Real_Vector_Handle; INCX : FP) with
      Import        => True,
      Convention    => Fortran,
      External_Name => BLAS_Prefix & "sscal" & BLAS_Suffix;

   procedure DSCAL (N : FP; DA : Real; DX : Real_Vector_Handle; INCX : FP) with
      Import        => True,
      Convention    => Fortran,
      External_Name => BLAS_Prefix & "dscal" & BLAS_Suffix;

   procedure SCOPY
     (N    : FI;
      SX   : Real_Vector_Constant_Handle;
      INCX : FP;
      SY   : Real_Vector_Handle;
      INCY : FP) with
      Import        => True,
      Convention    => Fortran,
      External_Name => BLAS_Prefix & "scopy" & BLAS_Suffix;

   procedure DCOPY
     (N    : FI;
      DX   : Real_Vector_Constant_Handle;
      INCX : FP;
      DY   : Real_Vector_Handle;
      INCY : FP) with
      Import        => True,
      Convention    => Fortran,
      External_Name => BLAS_Prefix & "dcopy" & BLAS_Suffix;

   procedure SAXPY
     (N    : FI;
      SA   : Real;
      SX   : Real_Vector_Constant_Handle;
      INCX : FP;
      SY   : Real_Vector_Handle;
      INCY : FP) with
      Import        => True,
      Convention    => Fortran,
      External_Name => BLAS_Prefix & "saxpy" & BLAS_Suffix;

   procedure DAXPY
     (N    : FI;
      DA   : Real;
      DX   : Real_Vector_Constant_Handle;
      INCX : FP;
      DY   : Real_Vector_Handle;
      INCY : FP) with
      Import        => True,
      Convention    => Fortran,
      External_Name => BLAS_Prefix & "daxpy" & BLAS_Suffix;

   function SDOT
     (N    : FI;
      SX   : Real_Vector_Constant_Handle;
      INCX : FP;
      SY   : Real_Vector_Constant_Handle;
      INCY : FP) return Real with
      Import        => True,
      Convention    => Fortran,
      External_Name => BLAS_Prefix & "sdot" & BLAS_Suffix;

   function DDOT
     (N    : FI;
      DX   : Real_Vector_Constant_Handle;
      INCX : FP;
      DY   : Real_Vector_Constant_Handle;
      INCY : FP) return Real with
      Import        => True,
      Convention    => Fortran,
      External_Name => BLAS_Prefix & "ddot" & BLAS_Suffix;

   function SDSDOT
     (N    : FI;
      SB   : Real;
      SX   : Real_Vector_Constant_Handle;
      INCX : FP;
      SY   : Real_Vector_Constant_Handle;
      INCY : FP) return Real with
      Import        => True,
      Convention    => Fortran,
      External_Name => BLAS_Prefix & "sdsdot" & BLAS_Suffix;

   function SNRM2
     (N    : FI;
      X    : Real_Vector_Constant_Handle;
      INCX : FP) return Real with
      Import        => True,
      Convention    => Fortran,
      External_Name => BLAS_Prefix & "snrm2" & BLAS_Suffix;

   function DNRM2
     (N    : FI;
      X    : Real_Vector_Constant_Handle;
      INCX : FP) return Real with
      Import        => True,
      Convention    => Fortran,
      External_Name => BLAS_Prefix & "dnrm2" & BLAS_Suffix;

   function SASUM
     (N    : FP;
      SX   : Real_Vector_Constant_Handle;
      INCX : FP) return Real with
      Import        => True,
      Convention    => Fortran,
      External_Name => BLAS_Prefix & "sasum" & BLAS_Suffix;

   function DASUM
     (N    : FP;
      DX   : Real_Vector_Constant_Handle;
      INCX : FP) return Real with
      Import        => True,
      Convention    => Fortran,
      External_Name => BLAS_Prefix & "dasum" & BLAS_Suffix;

   function ISAMAX
     (N    : FI;
      SX   : Real_Vector_Constant_Handle;
      INCX : FP) return FI with
      Import        => True,
      Convention    => Fortran,
      External_Name => BLAS_Prefix & "isamax" & BLAS_Suffix;

   function IDAMAX
     (N    : FI;
      DX   : Real_Vector_Constant_Handle;
      INCX : FP) return FI with
      Import        => True,
      Convention    => Fortran,
      External_Name => BLAS_Prefix & "idamax" & BLAS_Suffix;
--
      -- *************
      -- *************
      -- ** Level 2 **
      -- *************
      -- *************

   procedure SGEMV
     (TRANS : FC;
      M     : FP;
      N     : FP;
      ALPHA : Real;
      A     : Real_Matrix_Constant_Handle;
      LDA   : FP;
      X     : Real_Vector_Constant_Handle;
      INCX  : FP;
      BETA  : Real;
      Y     : Real_Vector_Handle;
      INCY  : FP) with
      Import        => True,
      Convention    => Fortran,
      External_Name => BLAS_Prefix & "sgemv" & BLAS_Suffix;

   procedure DGEMV
     (TRANS : FC;
      M     : FP;
      N     : FP;
      ALPHA : Real;
      A     : Real_Matrix_Constant_Handle;
      LDA   : FP;
      X     : Real_Vector_Constant_Handle;
      INCX  : FP;
      BETA  : Real;
      Y     : Real_Vector_Handle;
      INCY  : FP) with
      Import        => True,
      Convention    => Fortran,
      External_Name => BLAS_Prefix & "dgemv" & BLAS_Suffix;

   procedure SSYMV
     (UPLO  : FC;
      N     : FP;
      ALPHA : Real;
      A     : Real_Matrix_Constant_Handle;
      LDA   : FP;
      X     : Real_Vector_Constant_Handle;
      INCX  : FP;
      BETA  : Real;
      Y     : Real_Vector_Handle;
      INCY  : FP) with
      Import        => True,
      Convention    => Fortran,
      External_Name => BLAS_Prefix & "ssymv" & BLAS_Suffix;

   procedure DSYMV
     (UPLO  : FC;
      N     : FP;
      ALPHA : Real;
      A     : Real_Matrix_Constant_Handle;
      LDA   : FP;
      X     : Real_Vector_Constant_Handle;
      INCX  : FP;
      BETA  : Real;
      Y     : Real_Vector_Handle;
      INCY  : FP) with
      Import        => True,
      Convention    => Fortran,
      External_Name => BLAS_Prefix & "dsymv" & BLAS_Suffix;

   procedure SSPMV
     (UPLO  : FC;
      N     : FP;
      ALPHA : Real;
      AP    : Packed_Real_Matrix_Constant_Handle;
      X     : Real_Vector_Constant_Handle;
      INCX  : FP;
      BETA  : Real;
      Y     : Real_Vector_Handle;
      INCY  : FP) with
     Import        => True,
     Convention    => Fortran,
     External_Name => BLAS_Prefix & "sspmv" & BLAS_Suffix;

   procedure DSPMV
     (UPLO  : FC;
      N     : FP;
      ALPHA : Real;
      AP    : Packed_Real_Matrix_Constant_Handle;
      X     : Real_Vector_Constant_Handle;
      INCX  : FP;
      BETA  : Real;
      Y     : Real_Vector_Handle;
      INCY  : FP) with
     Import        => True,
     Convention    => Fortran,
     External_Name => BLAS_Prefix & "dspmv" & BLAS_Suffix;

   procedure STRMV
     (UPLO  : FC;
      TRANS : FC;
      DIAG  : FC;
      N     : FP;
      A     : Real_Matrix_Constant_Handle;
      LDA   : FP;
      X     : Real_Vector_Handle;
      INCX  : FP) with
      Import        => True,
      Convention    => Fortran,
      External_Name => BLAS_Prefix & "strmv" & BLAS_Suffix;

   procedure DTRMV
     (UPLO  : FC;
      TRANS : FC;
      DIAG  : FC;
      N     : FP;
      A     : Real_Matrix_Constant_Handle;
      LDA   : FP;
      X     : Real_Vector_Handle;
      INCX  : FP) with
      Import        => True,
      Convention    => Fortran,
      External_Name => BLAS_Prefix & "dtrmv" & BLAS_Suffix;

   procedure STPMV
     (UPLO  : FC;
      TRANS : FC;
      DIAG  : FC;
      N     : FP;
      AP    : Packed_Real_Matrix_Constant_Handle;
      X     : Real_Vector_Handle;
      INCX  : FP) with
     Import        => True,
     Convention    => Fortran,
     External_Name => BLAS_Prefix & "stpmv" & BLAS_Suffix;

   procedure DTPMV
     (UPLO  : FC;
      TRANS : FC;
      DIAG  : FC;
      N     : FP;
      AP    : Packed_Real_Matrix_Constant_Handle;
      X     : Real_Vector_Handle;
      INCX  : FP) with
     Import        => True,
     Convention    => Fortran,
     External_Name => BLAS_Prefix & "dtpmv" & BLAS_Suffix;

   procedure STRSV
     (UPLO  : FC;
      TRANS : FC;
      DIAG  : FC;
      N     : FP;
      A     : Real_Matrix_Constant_Handle;
      LDA   : FP;
      X     : Real_Vector_Handle;
      INCX  : FP) with
      Import        => True,
      Convention    => Fortran,
      External_Name => BLAS_Prefix & "strsv" & BLAS_Suffix;

   procedure DTRSV
     (UPLO  : FC;
      TRANS : FC;
      DIAG  : FC;
      N     : FP;
      A     : Real_Matrix_Constant_Handle;
      LDA   : FP;
      X     : Real_Vector_Handle;
      INCX  : FP) with
      Import        => True,
      Convention    => Fortran,
     External_Name => BLAS_Prefix & "dtrsv" & BLAS_Suffix;

   procedure STPSV
     (UPLO  : FC;
      TRANS : FC;
      DIAG  : FC;
      N     : FP;
      AP    : Packed_Real_Matrix_Constant_Handle;
      X     : Real_Vector_Handle;
      INCX  : FP) with
      Import        => True,
      Convention    => Fortran,
      External_Name => BLAS_Prefix & "stpsv" & BLAS_Suffix;

   procedure DTPSV
     (UPLO  : FC;
      TRANS : FC;
      DIAG  : FC;
      N     : FP;
      AP    : Packed_Real_Matrix_Constant_Handle;
      X     : Real_Vector_Handle;
      INCX  : FP) with
     Import        => True,
     Convention    => Fortran,
     External_Name => BLAS_Prefix & "dtpsv" & BLAS_Suffix;

   procedure SGER
     (M     : FI;
      N     : FI;
      ALPHA : Real;
      X     : Real_Vector_Constant_Handle;
      INCX  : FP;
      Y     : Real_Vector_Constant_Handle;
      INCY  : FP;
      A     : Real_Matrix_Handle;
      LDA   : FP) with
      Import        => True,
      Convention    => Fortran,
      External_Name => BLAS_Prefix & "sger" & BLAS_Suffix;

   procedure DGER
     (M     : FI;
      N     : FI;
      ALPHA : Real;
      X     : Real_Vector_Constant_Handle;
      INCX  : FP;
      Y     : Real_Vector_Constant_Handle;
      INCY  : FP;
      A     : Real_Matrix_Handle;
      LDA   : FP) with
      Import        => True,
      Convention    => Fortran,
      External_Name => BLAS_Prefix & "dger" & BLAS_Suffix;

   procedure SSYR
     (UPLO  : FC;
      N     : FI;
      ALPHA : Real;
      X     : Real_Vector_Constant_Handle;
      INCX  : FP;
      A     : Real_Matrix_Handle;
      LDA   : FP) with
      Import        => True,
      Convention    => Fortran,
      External_Name => BLAS_Prefix & "ssyr" & BLAS_Suffix;

   procedure DSYR
     (UPLO  : FC;
      N     : FI;
      ALPHA : Real;
      X     : Real_Vector_Constant_Handle;
      INCX  : FP;
      A     : Real_Matrix_Handle;
      LDA   : FP) with
      Import        => True,
      Convention    => Fortran,
      External_Name => BLAS_Prefix & "dsyr" & BLAS_Suffix;

   procedure SSPR
     (UPLO  : FC;
      N     : FI;
      ALPHA : Real;
      X     : Real_Vector_Constant_Handle;
      INCX  : FP;
      AP    : Packed_Real_Matrix_Handle) with
      Import        => True,
      Convention    => Fortran,
      External_Name => BLAS_Prefix & "sspr" & BLAS_Suffix;

   procedure DSPR
     (UPLO  : FC;
      N     : FI;
      ALPHA : Real;
      X     : Real_Vector_Constant_Handle;
      INCX  : FP;
      AP    : Packed_Real_Matrix_Handle) with
      Import        => True,
      Convention    => Fortran,
      External_Name => BLAS_Prefix & "dspr" & BLAS_Suffix;

   procedure SSYR2
     (UPLO  : FC;
      N     : FI;
      ALPHA : Real;
      X     : Real_Vector_Constant_Handle;
      INCX  : FP;
      Y     : Real_Vector_Constant_Handle;
      INCY  : FP;
      A     : Real_Matrix_Handle;
      LDA   : FP) with
      Import        => True,
      Convention    => Fortran,
      External_Name => BLAS_Prefix & "ssyr2" & BLAS_Suffix;

   procedure DSYR2
     (UPLO  : FC;
      N     : FI;
      ALPHA : Real;
      X     : Real_Vector_Constant_Handle;
      INCX  : FP;
      Y     : Real_Vector_Constant_Handle;
      INCY  : FP;
      A     : Real_Matrix_Handle;
      LDA   : FP) with
      Import        => True,
      Convention    => Fortran,
      External_Name => BLAS_Prefix & "dsyr2" & BLAS_Suffix;

   procedure SSPR2
     (UPLO  : FC;
      N     : FI;
      ALPHA : Real;
      X     : Real_Vector_Constant_Handle;
      INCX  : FP;
      Y     : Real_Vector_Constant_Handle;
      INCY  : FP;
      AP    : Packed_Real_Matrix_Handle) with
      Import        => True,
      Convention    => Fortran,
      External_Name => BLAS_Prefix & "sspr2" & BLAS_Suffix;

   procedure DSPR2
     (UPLO  : FC;
      N     : FI;
      ALPHA : Real;
      X     : Real_Vector_Constant_Handle;
      INCX  : FP;
      Y     : Real_Vector_Constant_Handle;
      INCY  : FP;
      AP    : Packed_Real_Matrix_Handle) with
      Import        => True,
      Convention    => Fortran,
      External_Name => BLAS_Prefix & "dspr2" & BLAS_Suffix;

      -- *************
      -- *************
      -- ** Level 3 **
      -- *************
      -- *************

   procedure SGEMM
     (TRANA : FC;
      TRANB : FC;
      M     : FI;
      N     : FI;
      K     : FI;
      ALPHA : Real;
      A     : Real_Matrix_Constant_Handle;
      LDA   : FI;
      B     : Real_Matrix_Constant_Handle;
      LDB   : FI;
      BETA  : Real;
      C     : Real_Matrix_Handle;
      LDC   : FI) with
      Import        => True,
      Convention    => Fortran,
      External_Name => BLAS_Prefix & "sgemm" & BLAS_Suffix;

   procedure DGEMM
     (TRANA : FC;
      TRANB : FC;
      M     : FI;
      N     : FI;
      K     : FI;
      ALPHA : Real;
      A     : Real_Matrix_Constant_Handle;
      LDA   : FI;
      B     : Real_Matrix_Constant_Handle;
      LDB   : FI;
      BETA  : Real;
      C     : Real_Matrix_Handle;
      LDC   : FI) with
      Import        => True,
      Convention    => Fortran,
      External_Name => BLAS_Prefix & "dgemm" & BLAS_Suffix;

   procedure SSYMM
     (SIDE  : FC;
      UPLO  : FC;
      M     : FI;
      N     : FI;
      ALPHA : Real;
      A     : Real_Matrix_Constant_Handle;
      LDA   : FI;
      B     : Real_Matrix_Constant_Handle;
      LDB   : FI;
      BETA  : Real;
      C     : Real_Matrix_Handle;
      LDC   : FI) with
      Import        => True,
      Convention    => Fortran,
      External_Name => BLAS_Prefix & "ssymm" & BLAS_Suffix;

   procedure DSYMM
     (SIDE  : FC;
      UPLO  : FC;
      M     : FI;
      N     : FI;
      ALPHA : Real;
      A     : Real_Matrix_Constant_Handle;
      LDA   : FI;
      B     : Real_Matrix_Constant_Handle;
      LDB   : FI;
      BETA  : Real;
      C     : Real_Matrix_Handle;
      LDC   : FI) with
      Import        => True,
      Convention    => Fortran,
      External_Name => BLAS_Prefix & "dsymm" & BLAS_Suffix;

   procedure SSYRK
     (UPLO  : FC;
      TRANS : FC;
      N     : FI;
      K     : FI;
      ALPHA : Real;
      A     : Real_Matrix_Constant_Handle;
      LDA   : FI;
      BETA  : Real;
      C     : Real_Matrix_Handle;
      LDC   : FI) with
      Import        => True,
      Convention    => Fortran,
      External_Name => BLAS_Prefix & "ssyrk" & BLAS_Suffix;

   procedure DSYRK
     (UPLO  : FC;
      TRANS : FC;
      N     : FI;
      K     : FI;
      ALPHA : Real;
      A     : Real_Matrix_Constant_Handle;
      LDA   : FI;
      BETA  : Real;
      C     : Real_Matrix_Handle;
      LDC   : FI) with
      Import        => True,
      Convention    => Fortran,
      External_Name => BLAS_Prefix & "dsyrk" & BLAS_Suffix;

   procedure SSYR2K
     (UPLO  : FC;
      TRANS : FC;
      N     : FI;
      K     : FI;
      ALPHA : Real;
      A     : Real_Matrix_Constant_Handle;
      LDA   : FI;
      B     : Real_Matrix_Constant_Handle;
      LDB   : FI;
      BETA  : Real;
      C     : Real_Matrix_Handle;
      LDC   : FI) with
      Import        => True,
      Convention    => Fortran,
      External_Name => BLAS_Prefix & "ssyr2k" & BLAS_Suffix;

   procedure DSYR2K
     (UPLO  : FC;
      TRANS : FC;
      N     : FI;
      K     : FI;
      ALPHA : Real;
      A     : Real_Matrix_Constant_Handle;
      LDA   : FI;
      B     : Real_Matrix_Constant_Handle;
      LDB   : FI;
      BETA  : Real;
      C     : Real_Matrix_Handle;
      LDC   : FI) with
      Import        => True,
      Convention    => Fortran,
      External_Name => BLAS_Prefix & "dsyr2k" & BLAS_Suffix;

   procedure STRMM
     (SIDE   : FC;
      UPLO   : FC;
      TRANSA : FC;
      DIAG   : FC;
      M      : FI;
      N      : FI;
      ALPHA  : Real;
      A      : Real_Matrix_Constant_Handle;
      LDA    : FI;
      B      : Real_Matrix_Handle;
      LDB    : FI) with
      Import        => True,
      Convention    => Fortran,
      External_Name => BLAS_Prefix & "strmm" & BLAS_Suffix;

    procedure DTRMM
     (SIDE   : FC;
      UPLO   : FC;
      TRANSA : FC;
      DIAG   : FC;
      M      : FI;
      N      : FI;
      ALPHA  : Real;
      A      : Real_Matrix_Constant_Handle;
      LDA    : FI;
      B      : Real_Matrix_Handle;
      LDB    : FI) with
      Import        => True,
      Convention    => Fortran,
      External_Name => BLAS_Prefix & "dtrmm" & BLAS_Suffix;

   procedure STRSM
     (SIDE   : FC;
      UPLO   : FC;
      TRANSA : FC;
      DIAG   : FC;
      M      : FI;
      N      : FI;
      ALPHA  : Real;
      A      : Real_Matrix_Constant_Handle;
      LDA    : FI;
      B      : Real_Handle; -- Can take matrices or vectors
      LDB    : FI) with
      Import        => True,
      Convention    => Fortran,
      External_Name => BLAS_Prefix & "strsm" & BLAS_Suffix;

   procedure DTRSM
     (SIDE   : FC;
      UPLO   : FC;
      TRANSA : FC;
      DIAG   : FC;
      M      : FI;
      N      : FI;
      ALPHA  : Real;
      A      : Real_Matrix_Constant_Handle;
      LDA    : FI;
      B      : Real_Handle; -- Can take matrices or vectors
      LDB    : FI) with
      Import        => True,
      Convention    => Fortran,
      External_Name => BLAS_Prefix & "dtrsm" & BLAS_Suffix;

end LALG.Real_BLAS.Imports;
