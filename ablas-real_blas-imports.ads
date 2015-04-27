-- aBLAS
-- An Ada 2012 binding to BLAS

private generic package aBLAS.Real_BLAS.Imports is

   subtype FI is IntFort.Fortran_Integer;

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

end aBLAS.Real_BLAS.Imports;
