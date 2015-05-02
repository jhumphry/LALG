

with Ada.Text_IO, Ada.Float_Text_IO;
use Ada.Text_IO, Ada.Float_Text_IO;
with Ada.Numerics.Generic_Real_Arrays;

with aBLAS, aBLAS.Real_BLAS, aBLAS.Real_BLAS.Util;
use aBLAS;

procedure Simple_Example is
   package GRA is new Ada.Numerics.Generic_Real_Arrays(Real => Float);
   use GRA;

   package BLAS is new aBLAS.Real_BLAS(Real => Float,
                                       Real_Vector => GRA.Real_Vector,
                                       Real_Matrix => GRA.Real_Matrix);

   package BLAS_Util is new BLAS.Util(Default_Aft => 3, Default_Exp => 0);
   use BLAS_Util;

   SX : GRA.Real_Vector(1..3) := (1.0, 2.0, 3.0);
   SY : GRA.Real_Vector(1..3) := (6.0, 5.0, 4.0);
   SZ : GRA.Real_Vector(1..3) := SY;

   X  : GRA.Real_Vector(1..3) := (1.0, 2.0, 3.0);
   A  : GRA.Real_Matrix(1..2, 1..3) := ((4.0, 5.0, 6.0),
                                        (7.0, 8.0, 9.0));
   Y  : GRA.Real_Vector(1..2) := (100.0, 100.0);

   B  : GRA.Real_Matrix(1..3, 1..2) := ((11.0, 12.0),
                                        (13.0, 14.0),
                                        (15.0, 16.0));

   C  : GRA.Real_Matrix(1..2, 1..2) := (others => (others => 99.0));

   BV : GRA.Real_Vector(1..6) := (others => 3.5);
   BM : GRA.Real_Matrix(1..6, 1..6) := ((1.0, -1.0, others => 0.7),
                                        (0.5, 2.0, others => 0.8),
                                       others => (0.33, 0.67, others => 0.9));

begin

   Put_Line("*** Level 1 ***");
   Put("SX is: ");
   Put(SX); New_Line;
   Put("SY is: ");
   Put(SY); New_Line;
   Put("SZ is: ");
   Put(SZ); New_Line;

   Put("SX+SZ via Ada library is: ");
   Put(SX+SZ); New_Line;
   Put("SX+SZ via BLAS is: ");
   BLAS.axpy(SX, SZ);
   Put(SZ); New_Line;
   New_Line;

   Put("Dot Product SX*SY via Ada library is: ");
   Put(Float'(SX * SY)); New_Line;
   Put("Dot Product SX*SY via BLAS is: ");
   Put(BLAS.dot(SX, SY)); New_Line;
   Put("Dot Product SX*SY via BLAS with extended precision accumulation is: ");
   Put(BLAS.sdsdot(SX, SY)); New_Line;
   Put("Dot Product SX*(SY reversed) is: ");
   Put(BLAS.dot(SX, SY, INCY => Increment(-1)));
   New_Line;

   Put("|SX|_1 is via BLAS library is: ");
   Put(BLAS.asum(SX)); New_Line;
   Put("Euclidean norm of SX is via Ada library is: ");
   Put(Float'(abs(SX))); New_Line;
   Put("Euclidean norm of SX is via BLAS library is: ");
   Put(BLAS.nrm2(SX)); New_Line;

   New_Line;
   Put_Line("*** Level 2 ***");
   Put("Y is: ");
   Put(Y); New_Line;
   Put("X is: ");
   Put(X); New_Line;
   Put_Line("A is: ");
   Put(A); New_Line;
   Put("A*x via Ada library: ");
   Put(A*x); New_Line;
   Put("Perform Y <- A*x with BLAS. Y is now: ");
   BLAS.gemv(A, X, Y);
   Put(Y); New_Line;

   New_Line;
   Put_Line("*** Level 3 ***");
   Put("A is: ");
   Put(A); New_Line;
   Put("B is: ");
   Put(B); New_Line;
   Put("C is: ");
   Put(C); New_Line;
   Put("A*B via Ada library: ");
   Put(A*B); New_Line;
   Put("Perform C <- A*B with BLAS. C is now: ");
   BLAS.gemm(A, B, C);
   Put(C); New_Line;
   Put("Perform C <- Bt*At with BLAS. C is now: ");
   BLAS.gemm(B, A, C, TRANA => Transpose, TRANB => Transpose);
   Put(C); New_Line;

   New_Line;
   Put("BV is: ");
   Put(BV);  New_Line;
   Put_Line("BM is: ");
   Put(BM); New_Line;

end Simple_Example;
