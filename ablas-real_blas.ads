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

   Precision : constant Precision_Specification :=
     (if Real'Base'Digits = IntFort.Real'Base'Digits and
        Real'Base'Size = IntFort.Real'Base'Size then Single
      elsif Real'Base'Digits = IntFort.Double_Precision'Base'Digits and
        Real'Base'Size = IntFort.Double_Precision'Base'Size then Double
      else raise Program_Error with "Precision not supported for interfacing with Fortran code");

   -- *************
   -- *************
   -- ** Level 1 **
   -- *************
   -- *************

   -- Y <-> X
   procedure swap(X : in out Real_Vector;
                  Y : in out Real_Vector;
                  INCX : in Increment := 1;
                  INCY : in Increment := 1;
                  N : in Vector_Size := 0)
     with Inline;

   -- X <- aX
   procedure scal(X : in out Real_Vector;
                  A : in Real := 1.0;
                  INCX : in Increment := 1;
                  N : in Vector_Size := 0)
     with Inline;

   -- Y <- X
   procedure copy(X : in Real_Vector;
                  Y : out Real_Vector;
                  INCX : in Increment := 1;
                  INCY : in Increment := 1;
                  N : in Vector_Size := 0)
     with Inline;

   -- Y <- aX + Y
   procedure axpy(X : in Real_Vector;
                  Y : in out Real_Vector;
                  A : in Real := 1.0;
                  INCX : in Increment := 1;
                  INCY : in Increment := 1;
                  N : in Vector_Size := 0)
     with Inline;

   -- dot <- X^T . Y
   function dot(X, Y : in Real_Vector;
                INCX : in Increment := 1;
                INCY : in Increment := 1;
                N : in Vector_Size := 0) return Real
     with Inline;

   -- sdsdot <- SX^T . SY + SB with accumulation done in extended precision
   -- If Real is already double precision, this is the same as using the regular
   -- dot function and adding SB
   function sdsdot(SX, SY : in Real_Vector;
                   SB : in Real := 0.0;
                   INCX : in Increment := 1;
                   INCY : in Increment := 1;
                   N : in Vector_Size := 0) return Real
     with Inline;

   -- nrm2 <- sqrt(X^T . X)
   function nrm2(X : in Real_Vector;
                 INCX : in Increment := 1;
                 N : in Vector_Size := 0) return Real
     with Inline;

   --  asum <- |X|_1
   function asum(X : in Real_Vector;
                 INCX : in Increment := 1;
                 N : in Vector_Size := 0) return Real
     with Inline;

   --  iamax <- 1st k where X_k = MAX(abs(X_k))
   function iamax(X : in Real_Vector;
                  INCX : in Increment := 1;
                  N : in Vector_Size := 0) return Integer
     with Inline;

   -- *************
   -- *************
   -- ** Level 2 **
   -- *************
   -- *************

   -- y <- alpha*A*x + beta*y
   procedure gemv(A : in Real_Matrix;
                  X : in Real_Vector;
                  Y : in out Real_Vector;
                  ALPHA : in Real := 1.0;
                  BETA : in Real := 0.0;
                  TRANS : in Real_Trans_Op := No_Transpose;
                  INCX, INCY : in Increment := 1;
                  M, N : in Vector_Size := 0;
                  Convention : in Matrix_Convention := Default_Matrix_Convention)
     with Inline;

   -- gemv <- alpha*A*x
   function gemv(A : in Real_Matrix;
                 X : in Real_Vector;
                 ALPHA : in Real := 1.0;
                 TRANS : in Real_Trans_Op := No_Transpose;
                 INCX : in Increment := 1;
                 M, N : in Vector_Size := 0;
                 Convention : in Matrix_Convention := Default_Matrix_Convention)
                 return Real_Vector
     with Inline;

   -- *************
   -- *************
   -- ** Level 3 **
   -- *************
   -- *************

   -- C <- alpha*A*B + beta*C
   procedure gemm(A : in Real_Matrix;
                  B : in Real_Matrix;
                  C : in out Real_Matrix;
                  ALPHA : in Real := 1.0;
                  BETA : in Real := 0.0;
                  TRANA : in Real_Trans_Op := No_Transpose;
                  TRANB : in Real_Trans_Op := No_Transpose;
                  M, N, K : in Vector_Size := 0;
                  Convention : in Matrix_Convention := Default_Matrix_Convention)
     with Inline;

   -- gemm <- alpha*A*B
   function gemm(A : in Real_Matrix;
                 B : in Real_Matrix;
                 ALPHA : in Real := 1.0;
                 TRANA : in Real_Trans_Op := No_Transpose;
                 TRANB : in Real_Trans_Op := No_Transpose;
                 M, N, K : in Vector_Size := 0;
                 Convention : in Matrix_Convention := Default_Matrix_Convention)
                 return Real_Matrix
     with Inline;

end aBLAS.Real_BLAS;
