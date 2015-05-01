-- aBLAS
-- An Ada 2012 binding to BLAS

private generic package aBLAS.Real_BLAS.Imports is

   subtype FI is IntFort.Fortran_Integer;
   subtype FC is IntFort.Character_Set;

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

end aBLAS.Real_BLAS.Imports;
