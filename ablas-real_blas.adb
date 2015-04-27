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
      NN : Vector_Size;
   begin
      if N = 0 then
         NN := SX'Length;
      else
         NN := N;
      end if;

      if Precision = Single then
         return SDOT(NN,SX,INCX,SY,INCY);
      else
         return DDOT(NN,SX,INCX,SY,INCY);
      end if;
   end dot;

end aBLAS.Real_BLAS;
