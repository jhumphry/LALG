-- aBLAS
-- An Ada 2012 binding to BLAS

private generic package aBLAS.Real_BLAS.Imports is

   subtype FI is IntFort.Fortran_Integer;
   subtype FC is IntFort.Character_Set;

   BLAS_Prefix : constant String := "";
   BLAS_Suffix : constant String := "_";

   -- *************
   -- *************
   -- ** Level 1 **
   -- *************
   -- *************

   procedure SROTG(SA, SB : in out Real; C, S : out Real)
     with Import => True,
     Convention => Fortran,
     External_Name => BLAS_Prefix & "srotg" & BLAS_Suffix;

   procedure DROTG(SA, SB : in out Real; C, S : out Real)
     with Import => True,
     Convention => Fortran,
     External_Name => BLAS_Prefix & "drotg" & BLAS_Suffix;

   procedure SROT(N : FI;
                  SX : in out Real_Vector;
                  INCX : Increment;
                  SY : in out Real_Vector;
                  INCY : Increment;
                  C : Real;
                  S : Real)
   with Import => True,
   Convention => Fortran,
   External_Name => BLAS_Prefix & "srot" & BLAS_Suffix;

   procedure DROT(N : FI;
                  SX : in out Real_Vector;
                  INCX : Increment;
                  SY : in out Real_Vector;
                  INCY : Increment;
                  C : Real;
                  S : Real)
   with Import => True,
   Convention => Fortran,
   External_Name => BLAS_Prefix & "drot" & BLAS_Suffix;

   procedure SSWAP(N : FI;
                   SX : in out Real_Vector;
                   INCX : Increment;
                   SY : in out Real_Vector;
                   INCY : Increment)
     with Import => True,
     Convention => Fortran,
     External_Name => BLAS_Prefix & "sswap" & BLAS_Suffix;

   procedure DSWAP(N : FI;
                   DX : in out Real_Vector;
                   INCX : Increment;
                   DY : in out Real_Vector;
                   INCY : Increment)
     with Import => True,
     Convention => Fortran,
     External_Name => BLAS_Prefix & "dswap" & BLAS_Suffix;

   procedure SSCAL(N : FI;
                   SA : Real;
                   SX : in out Real_Vector;
                   INCX : Increment)
     with Import => True,
     Convention => Fortran,
     External_Name => BLAS_Prefix & "sscal" & BLAS_Suffix;

   procedure DSCAL(N : FI;
                   DA : Real;
                   DX : in out Real_Vector;
                   INCX : Increment)
     with Import => True,
     Convention => Fortran,
     External_Name => BLAS_Prefix & "dscal" & BLAS_Suffix;

   procedure SCOPY(N : FI;
                   SX : Real_Vector;
                   INCX : Increment;
                   SY : out Real_Vector;
                   INCY : Increment)
     with Import => True,
     Convention => Fortran,
     External_Name => BLAS_Prefix & "scopy" & BLAS_Suffix;

   procedure DCOPY(N : FI;
                   DX : Real_Vector;
                   INCX : Increment;
                   DY : out Real_Vector;
                   INCY : Increment)
     with Import => True,
     Convention => Fortran,
     External_Name => BLAS_Prefix & "dcopy" & BLAS_Suffix;

   procedure SAXPY(N : FI;
                   SA : Real;
                   SX : Real_Vector;
                   INCX : Increment;
                   SY : in out Real_Vector;
                   INCY : Increment)
     with Import => True,
     Convention => Fortran,
     External_Name => BLAS_Prefix & "saxpy" & BLAS_Suffix;

   procedure DAXPY(N : FI;
                   DA : Real;
                   DX : Real_Vector;
                   INCX : Increment;
                   DY : in out Real_Vector;
                   INCY : Increment)
     with Import => True,
     Convention => Fortran,
     External_Name => BLAS_Prefix & "daxpy" & BLAS_Suffix;

   function SDOT(N : FI;
                  SX : Real_Vector;
                  INCX : Increment;
                  SY : Real_Vector;
                  INCY : Increment) return Real
     with Import => True,
     Convention => Fortran,
     External_Name => BLAS_Prefix & "sdot" & BLAS_Suffix;

   function DDOT(N : FI;
                 DX : Real_Vector;
                 INCX : Increment;
                 DY : Real_Vector;
                 INCY : Increment) return Real
     with Import => True,
     Convention => Fortran,
     External_Name => BLAS_Prefix & "ddot" & BLAS_Suffix;

   function SDSDOT(N : FI;
                   SB : Real;
                   SX : Real_Vector;
                   INCX : Increment;
                   SY : Real_Vector;
                   INCY : Increment) return Real
     with Import => True,
     Convention => Fortran,
     External_Name => BLAS_Prefix & "sdsdot" & BLAS_Suffix;

   function SNRM2(N : FI;
                  X : Real_Vector;
                  INCX : Increment) return Real
     with Import => True,
     Convention => Fortran,
     External_Name => BLAS_Prefix & "snrm2" & BLAS_Suffix;

   function DNRM2(N : FI;
                  X : Real_Vector;
                  INCX : Increment) return Real
     with Import => True,
     Convention => Fortran,
     External_Name => BLAS_Prefix & "dnrm2" & BLAS_Suffix;

   function SASUM(N : FI;
                  SX : Real_Vector;
                  INCX : Increment) return Real
     with Import => True,
     Convention => Fortran,
     External_Name => BLAS_Prefix & "sasum" & BLAS_Suffix;

   function DASUM(N : FI;
                  DX : Real_Vector;
                  INCX : Increment) return Real
     with Import => True,
     Convention => Fortran,
     External_Name => BLAS_Prefix & "dasum" & BLAS_Suffix;

   function ISAMAX(N : FI;
                   SX : Real_Vector;
                   INCX : Increment) return FI
     with Import => True,
     Convention => Fortran,
     External_Name => BLAS_Prefix & "isamax" & BLAS_Suffix;

   function IDAMAX(N : FI;
                   DX : Real_Vector;
                   INCX : Increment) return FI
     with Import => True,
     Convention => Fortran,
     External_Name => BLAS_Prefix & "idamax" & BLAS_Suffix;

   -- *************
   -- *************
   -- ** Level 2 **
   -- *************
   -- *************

   procedure SGEMV(TRANS : FC;
                   M: FI;
                   N : FI;
                   ALPHA : Real;
                   A : Real_Matrix;
                   LDA : FI;
                   X : Real_Vector;
                   INCX : Increment;
                   BETA : Real;
                   Y : in out Real_Vector;
                   INCY : Increment
                  )
     with Import => True,
     Convention => Fortran,
     External_Name => BLAS_Prefix & "sgemv" & BLAS_Suffix;

   procedure DGEMV(TRANS : FC;
                   M: FI;
                   N : FI;
                   ALPHA : Real;
                   A : Real_Matrix;
                   LDA : FI;
                   X : Real_Vector;
                   INCX : Increment;
                   BETA : Real;
                   Y : in out Real_Vector;
                   INCY : Increment
                  )
     with Import => True,
     Convention => Fortran,
     External_Name => BLAS_Prefix & "dgemv" & BLAS_Suffix;

   -- *************
   -- *************
   -- ** Level 3 **
   -- *************
   -- *************

   procedure SGEMM(TRANA : FC;
                   TRANB : FC;
                   M: FI;
                   N : FI;
                   K : FI;
                   ALPHA : Real;
                   A : Real_Matrix;
                   LDA : FI;
                   B : Real_Matrix;
                   LDB : FI;
                   BETA : Real;
                   C : in out Real_Matrix;
                   LDC : FI
                  )
     with Import => True,
     Convention => Fortran,
     External_Name => BLAS_Prefix & "sgemm" & BLAS_Suffix;

   procedure DGEMM(TRANA : FC;
                   TRANB : FC;
                   M: FI;
                   N : FI;
                   K : FI;
                   ALPHA : Real;
                   A : Real_Matrix;
                   LDA : FI;
                   B : Real_Matrix;
                   LDB : FI;
                   BETA : Real;
                   C : in out Real_Matrix;
                   LDC : FI
                  )
     with Import => True,
     Convention => Fortran,
     External_Name => BLAS_Prefix & "dgemm" & BLAS_Suffix;

end aBLAS.Real_BLAS.Imports;
