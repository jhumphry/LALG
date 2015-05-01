-- aBLAS
-- An Ada 2012 binding to BLAS

private generic package aBLAS.Real_BLAS.Imports is

   subtype FI is IntFort.Fortran_Integer;
   subtype FC is IntFort.Character_Set;

   -- *************
   -- *************
   -- ** Level 1 **
   -- *************
   -- *************

   function SDOT(N : FI;
                  SX : Real_Vector;
                  INCX : Increment;
                  SY : Real_Vector;
                  INCY : Increment) return Real
     with Import => True,
     Convention => Fortran,
     External_Name => "sdot_";

   function DDOT(N : FI;
                  SX : Real_Vector;
                  INCX : Increment;
                  SY : Real_Vector;
                  INCY : Increment) return Real
     with Import => True,
     Convention => Fortran,
     External_Name => "ddot_";

   function SNRM2(N : FI;
                  X : Real_Vector;
                  INCX : Increment) return Real
     with Import => True,
     Convention => Fortran,
     External_Name => "snrm2_";

   function DNRM2(N : FI;
                  X : Real_Vector;
                  INCX : Increment) return Real
     with Import => True,
     Convention => Fortran,
     External_Name => "dnrm2_";

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
     External_Name => "sgemv_";

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
     External_Name => "dgemv_";

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
     External_Name => "sgemm_";

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
     External_Name => "dgemm_";

end aBLAS.Real_BLAS.Imports;
