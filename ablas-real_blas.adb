-- aBLAS
-- An Ada 2012 binding to BLAS

with aBLAS.Real_BLAS.Imports;

package body aBLAS.Real_BLAS is

   package Fortran_Imports is new aBLAS.Real_BLAS.Imports;
   use Fortran_Imports;

   ---------
   -- dot --
   ---------

   function dot
     (SX, SY: in Real_Vector;
      INCX : in Increment := 1;
      INCY : in Increment := 1;
      N : in Vector_Size := 0)
      return Real
   is
      NN : Vector_Size := (if N = 0 then SX'Length else N);
   begin
      if Precision = Single then
         return SDOT(NN,SX,INCX,SY,INCY);
      else
         return DDOT(NN,SX,INCX,SY,INCY);
      end if;
   end dot;

   ----------
   -- nrm2 --
   ----------

   function nrm2
     (SX : in Real_Vector;
      INCX : in Increment := 1;
      N : in Vector_Size := 0)
      return Real
   is
      NN : Vector_Size := (if N = 0 then SX'Length else N);
   begin
      if Precision = Single then
         return SNRM2(NN,SX,INCX);
      else
         return DNRM2(NN,SX,INCX);
      end if;
   end nrm2;


end aBLAS.Real_BLAS;
