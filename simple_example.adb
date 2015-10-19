
with Ada.Text_IO, Ada.Float_Text_IO;
use Ada.Text_IO, Ada.Float_Text_IO;

with aBLAS, aBLAS.Real_BLAS, aBLAS.Real_BLAS.Util;

procedure Simple_Example is

   package BLAS is new aBLAS(Real => Float);
   package Real_BLAS is new BLAS.Real_BLAS;
   use Real_BLAS;

   package BLAS_Util is new Real_BLAS.Util(Default_Aft => 3, Default_Exp => 0);
   use BLAS_Util;

   P : aliased BLAS.Concrete_Real_Vector := BLAS.Make((1.0, 2.0, 3.0, 4.0, 5.0));
   P_View1 : aliased BLAS.Real_Vector_View := BLAS.Make(V => P'Access,
                                                        Start => 2,
                                                        Stride => 2);
   P_View2 : aliased constant BLAS.Real_Vector_View := BLAS.Make(V => P'Access,
                                                                 Start => 2,
                                                                 Stride => 1,
                                                                 Length => 3);

   A : aliased BLAS.Concrete_Real_Matrix := BLAS.Make(((3.0, 4.0), (5.0, 6.0)));
   B : aliased BLAS.Concrete_Real_Matrix := BLAS.Make(((1.0, 2.0), (3.0, 4.0)));
   C : aliased BLAS.Concrete_Real_Matrix := BLAS.Zeros(2,2);

   X : aliased BLAS.Concrete_Real_Vector := BLAS.Make((1.0, 2.0));
   Y : aliased BLAS.Concrete_Real_Vector := BLAS.Make((5.0, 5.0));

   SP : aliased constant BLAS.Symmetric_Real_Matrix := BLAS.Make((
                                                        (1.0, 2.0),
                                                        (3.0, 4.0)
                                                       ));

begin

   Put_Line("Test output");
   Put("Zero vector length 15: ");
   Put(BLAS.Concrete_Real_Vector'(BLAS.Zeros(15))); New_Line;
   Put("Identity matrix 3x3:");
   Put(BLAS.Identity(3)); New_Line;
   Put("Identity matrix 10x10:");
   Put(BLAS.Identity(10)); New_Line;

   New_Line;
   Put_Line("*** Level 1 ***");

   Put("P => "); Put(P); New_Line;
   Put("P_View1 => "); Put(P_View1); New_Line;
   Put("P_View2 => "); Put(P_View2); New_Line;

   Put("|P|_1 is via BLAS library is: ");
   Put(asum(P)); New_Line;
   Put("|P_View1|_1 is via BLAS library is: ");
   Put(asum(P_View1)); New_Line;

   Put_Line("Scaling P_View1 by 2.0");
   scal(P_View1, 2.0);
   Put("P => "); Put(P); New_Line;
   Put("P_View1 => "); Put(P_View1); New_Line;
   Put("P_View2 => "); Put(P_View2); New_Line;

   New_Line;
   Put_Line("*** Level 2 ***");
   Put("A => "); Put(A); New_Line;
   Put("X => "); Put(X); New_Line;
   Put("Y => "); Put(Y); New_Line;
   Put_Line("Perform Y <- A*X + Y");
   gemv(A, X, Y, ALPHA => 1.0, BETA => 1.0);
   Put("Y => "); Put(Y); New_Line;

   New_Line;
   Put_Line("*** Slicing vectors from matrices ***");
   Put("A => "); Put(A); New_Line;
   Put("A(2,1) => "); Put(A(2,1)); New_Line;
   Put("|A(1,j)|_1 => "); Put(asum(A.Row(1))); New_Line;
   Put("|A(i,2)|_1 => "); Put(asum(A.Column(2))); New_Line;
   Put("diagonal(A) => "); Put(A.Diagonal); New_Line;
   Put("|diagonal(A)|_1 => "); Put(asum(A.Diagonal)); New_Line;
   Put_Line("Scale diagonal(A) by 2.5.");
   declare
      T : BLAS.Real_Matrix_Vector := A.Diagonal;
   begin
      scal(T, 2.5);
   end;
   Put("A => "); Put(A); New_Line;
   Put("|diagonal(A)|_1 => "); Put(asum(A.Diagonal)); New_Line;

   New_Line;
   Put_Line("*** Level 3 ***");
   Put("A => "); Put(A); New_Line;
   Put("B => "); Put(B); New_Line;
   Put("C => "); Put(C); New_Line;
   Put_Line("Perform C <- 2.0*A*B + C");
   gemm(A, B, C, ALPHA => 2.0, BETA => 1.0);
   Put("C => "); Put(C); New_Line;
   Put_Line("Perform C <- 2.0*A*B + C");
   gemm(A, B, C, ALPHA => 2.0, BETA => 1.0);
   Put("1.0*At*B => "); Put(gemm(A,B,1.0,Transpose,No_Transpose));
   New_Line;

   New_Line;
   Put_Line("Symmetrical Packed matrices");
   Put("SP => "); Put(SP); New_Line;

end Simple_Example;
