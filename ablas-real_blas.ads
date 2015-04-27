-- aBLAS
-- An Ada 2012 binding to BLAS

with Interfaces.Fortran;

generic
   type Real is digits <>;
   type Real_Vector is array (Integer range <>) of Real'Base;
   type Real_Matrix is array (Integer range <>, Integer range <>) of Real'Base;
   Default_Matrix_Convention : Matrix_Convention := Ada_Convention;
package aBLAS.Real_BLAS is

   package IntFort renames Interfaces.Fortran;

   Precision : constant Precision_Specification := (if Real'Base'Digits = IntFort.Real'Base'Digits and
                                                      Real'Base'Size = IntFort.Real'Base'Size then Single
                                                    elsif Real'Base'Digits = IntFort.Double_Precision'Base'Digits and
                                                      Real'Base'Size = IntFort.Double_Precision'Base'Size then Double
                                                    else Unsupported);


   function dot(SX, SY: in Real_Vector;
                INCX : in Increment := 1;
                INCY : in Increment := 1;
                N : in Vector_Size := 0) return Real
     with Inline;

end aBLAS.Real_BLAS;
